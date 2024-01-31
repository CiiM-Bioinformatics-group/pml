#!/usr/bin/env python3
# File: scpred.py
# Author:  Zhenhua Zhang
# E-mail:  zhenhua.zhang217@gmail.com
# Created: 2022 Mar 01
# Updated: Jan 28, 2024

import os, json, logging, argparse
from datetime import datetime as dt

import matplotlib.pyplot as plt
import joblib as jbl
import pickle as pkl
import pandas as pd
import numpy as np
import shap as sp

import torch
import torch.nn as tc_nn
import torch.nn.functional as tc_func
import torch.optim as tc_optim
from torch.utils.data.dataloader import DataLoader as tc_DataLoader
from torch.utils.data.dataloader import Dataset as tc_Dataset

import pyro
import pyro.distributions as pr_dist
import pyro.distributions.constraints as pr_const
import pyro.infer as pr_infer
import pyro.optim as pr_optim

from sklearn.base import BaseEstimator, ClassifierMixin
from sklearn.cluster import AffinityPropagation
from sklearn.decomposition import PCA, IncrementalPCA
from sklearn.ensemble import GradientBoostingClassifier, RandomForestClassifier
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.metrics import precision_recall_curve, precision_score, accuracy_score, roc_auc_score, recall_score, roc_curve
from sklearn.model_selection import RandomizedSearchCV, train_test_split
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


# Global random state for reproducibility. Ajustable by CLI options -s/--random-seed.
RANDOM_STATE = np.random.RandomState(42)


class LogManager(logging.Logger):
    def __init__(self, name, level=logging.INFO, logstream: bool = True, logfile: str = ""):
        super(LogManager, self).__init__(name)
        fmt = logging.Formatter("{levelname: >8} | {asctime} | {name: <9} | {message}", style="{", datefmt="%Y-%m-%d %H:%M:%S")
        if logstream: self._add_handler(logging.StreamHandler(), level, fmt)
        if logfile: self._add_handler(logging.FileHandler(logfile), level, fmt)

    def _add_handler(self, hdl, lvl, fmt):
        hdl.setLevel(lvl)
        hdl.setFormatter(fmt)
        self.addHandler(hdl)


class ExpDataSet(tc_Dataset):
    def __init__(self, X, y=None):
        super(ExpDataSet, self).__init__()
        self.x_mat, self.y_vec = X, y

    def __len__(self):
        return self.x_mat.shape[0]

    def __getitem__(self, idx: int):
        x_vec = torch.Tensor(self.x_mat[idx, :]).unsqueeze(0)

        if self.y_vec is None:
            y_vec = torch.tensor(torch.nan)
        else:
            y_vec = torch.tensor(self.y_vec.iloc[idx])

        return x_vec, y_vec

    def __next__(self):
        for x in range(len(self)):
            yield self[x]


class ScNN(tc_nn.Module):
    def __init__(self, fs: int = 32, cs: int = 2, insize: int = 32):
        super(ScNN, self).__init__()

        assert fs % 2 == 0, "The fully-connected input size should be even."

        self.ll_1 = tc_nn.Linear(insize, 128)
        self.ll_2 = tc_nn.Linear(128, 64)
        self.ll_3 = tc_nn.Linear(64, 32)
        self.ll_4 = tc_nn.Linear(32, 16)
        self.dp_1 = tc_nn.Dropout(0.25)
        self.ll_5 = tc_nn.Linear(16, cs)

    def forward(self, x):
        x = tc_func.relu(self.ll_1(x))
        x = tc_func.relu(self.ll_2(x))
        x = tc_func.relu(self.ll_3(x))
        x = tc_func.relu(self.ll_4(x))
        x = self.dp_1(x)
        x = tc_func.relu(self.ll_5(x))
        q_x = tc_func.softmax(x, dim=1)

        return q_x


class NNClassifier(BaseEstimator, ClassifierMixin):
    def __init__(self, lrate: float = 5e-6, neps: int = 500, log_pne: int = 25):
        super(NNClassifier, self).__init__()
        self.model: tc_nn.Module
        self._learning_rate = lrate
        self._neps = neps
        self._log_pne = log_pne

    @property
    def loss_func(self):
        return tc_nn.CrossEntropyLoss()

    @property
    def optimizer(self):
        return tc_optim.Adam(self.model.parameters(), lr=self._learning_rate)

    def fit(self, X, y, batch_size: int = 32, shuffle=True):
        _, in_size = X.shape
        self.model = ScNN(in_size)

        dtld = tc_DataLoader(ExpDataSet(X, y), batch_size=batch_size, shuffle=shuffle)

        for epoch in range(self._neps):
            loss_pep = torch.tensor(0)
            for _, (x_mat, y_true) in enumerate(dtld):
                y_pred = self.model(x_mat)
                y_true = y_true.type(torch.LongTensor)
                loss = self.loss_func(y_pred, y_true)

                self.optimizer.zero_grad()
                loss.backward()
                self.optimizer.step()

            if (epoch + 1) % self._log_pne == 0:
                print(f"Epoch {epoch+1:5d}, total loss {loss_pep:.3f}")

        return self

    def predict(self, X, pos_idx: int = 0, use_prob=False):
        dtld = tc_DataLoader(ExpDataSet(X), batch_size=1, shuffle=False)

        y_pprob, y_plabel = torch.Tensor(), torch.Tensor()

        if self.model is None:
            raise ValueError("Model wasn't ready yet, please use fit() first.")

        with torch.no_grad():
            for _, (x_test, _) in enumerate(dtld):
                y_pred = self.model(x_test)
                y_pprob = torch.concat((y_pprob, y_pred[:, pos_idx]))
                y_plabel = torch.concat((y_plabel, y_pred.argmax(1)))

            if use_prob:
                return y_pprob.data.item()
            else:
                return y_plabel.data.numpy().astype("i8")

    def predict_proba(self, X, pos_idx: int = 0):
        return self.predict(X, pos_idx, True)


class PredEnsembler:
    """A class to ensemble the predicted probabilities."""
    def __init__(self, data, n_iters=5000, logman: LogManager = LogManager("PredEnsembler")):
        self._logman = logman
        self._data = data
        self._n_iters = n_iters
        self._loss, self._alpha, self._beta, self._theta = [], [], [], []
        self._pp_val = 0
        
    @property
    def pp_val(self):
        return self._pp_val

    @property
    def theta(self):
        return self._theta[-1]

    @property
    def alpha(self):
        return self._alpha[-1]

    @property
    def beta(self):
        return self._beta[-1]

    # The model function
    @staticmethod
    def _model(data):
        if not isinstance(data, torch.Tensor): data = torch.Tensor(data)
        data = data * 0.99
        a_init, b_init = torch.tensor([5, 5])
        alpha = pyro.param("alpha", a_init, constraint=pr_const.positive)
        beta = pyro.param("beta", b_init, constraint=pr_const.positive)
        with pyro.plate("obs", len(data)):
            return pyro.sample("theta", pr_dist.Beta(alpha, beta), obs=data)

    @staticmethod
    def _guide(data):
        if not isinstance(data, torch.Tensor): data = torch.Tensor(data)
        a_init, b_init = torch.tensor([15, 15])
        alpha = pyro.param("alpha", a_init, constraint=pr_const.positive)
        beta = pyro.param("beta", b_init, constraint=pr_const.positive)

    def infer(self, adam_param={"lr": 1e-3}):
        """Infer theta from given"""
        adam = pr_optim.Adam(adam_param) # Optimizer
        elbo = pr_infer.Trace_ELBO() # Trace by ELBO (evidence lower bound)

        # Stochastic Variational Inference
        svi = pr_infer.SVI(self._model, self._guide, adam, elbo)

        # Sampling
        for _ in range(self._n_iters):
            try:
                per_loss = svi.step(self._data) # Collect loss
                per_beta, per_alpha = pyro.param("beta").item(), pyro.param("alpha").item() # Collect beta and alpha

                if per_loss and per_alpha and per_beta:
                    self._loss.append(per_loss)
                    self._alpha.append(per_alpha)
                    self._beta.append(per_beta)
            except ValueError as _:
                pass

        # The theta determined by the learned alpha and beta.
        self._theta = [a / (a + b) for a, b in zip(self._alpha, self._beta)]
    
    def cdf(self, x=0.5, n_points=10000): # Cumulative density/mass function
        theta = pr_dist.Beta(self.alpha, self.beta)
        x_space = torch.linspace(0, x, n_points)[1:-1]
        return torch.trapz(theta.log_prob(x_space).exp(), x_space)

    def report(self, pp=0.5, prefix="Report-", save_to="./"):
        # X-axis
        x_set = np.arange(self._n_iters)

        fig = plt.figure(figsize=(11, 12), constrained_layout=True)
        spec = fig.add_gridspec(3, 2)

        # Observation, theta, the most important parameters.
        ax_mean_lf = fig.add_subplot(spec[0, 0])
        ax_mean_lf.hist(self._data, alpha=0.25, color="green")
        ax_mean_lf.axvline(0.5, linestyle="--", color="gray", alpha=0.75)
        
        ## PP > pp, the pp is given by argument, e.g., 0.5
        self._pp_val = 100 * (1 - self.cdf().item())
        ax_mean_lf.set_title(f"$PP_{{{pp}}}$ = {self._pp_val:.2f}%")

        ## Sampling using the learned parameters
        pd_vi = pr_dist.Beta(self.alpha, self.beta)
        x_values = np.linspace(0, 1, num=1000)[1:-1]
        y_values = torch.exp(pd_vi.log_prob(torch.tensor(x_values)))
        ax_mean_lf_tx = ax_mean_lf.twinx()
        ax_mean_lf_tx.plot(x_values, y_values, color="blue")
        
        # Estimated theta per iteration
        ax_mean_rt = fig.add_subplot(spec[0, 1])
        ax_mean_rt.plot(x_set, self._theta)
        ax_mean_rt.axhline(0.5, linestyle="--", color="gray", alpha=0.75)
        ax_mean_rt.set_xlabel("Iterations")
        ax_mean_rt.set_ylabel("Estimated theta")
        ax_mean_rt.set_title("Estimated theta per iteration")

        # Parameter, alpha, beta
        ax_para_lf = fig.add_subplot(spec[1, 0])
        ax_para_lf.plot(x_set, self._alpha)
        ax_para_lf.plot(x_set, self._beta)
        ax_para_lf.set_xlabel("Iterations")
        ax_para_lf.set_ylabel("Learned alpha or beta value")
        ax_para_lf.set_title("Learned alpha and beta per step (diagnosis)")

        ax_para_rt = fig.add_subplot(spec[1, 1])
        ax_para_rt.hist(self._alpha, alpha=0.25)
        ax_para_rt.hist(self._beta, alpha=0.25)
        ax_para_rt.set_xlabel("Alpha or beta value")
        ax_para_rt.set_ylabel("Frequency")
        ax_para_rt.set_title("Learned alpha and beta per step bin (diagnosis)")
   
        # Loss
        ax_loss = fig.add_subplot(spec[2, :])
        ax_loss.plot(x_set, self._loss)
        ax_loss.set_xlabel("Number of iteration")
        ax_loss.set_ylabel("ELBO loss")
        ax_loss.set_title("ELBO loss (for diagnosis)")

        fig.savefig(f"{save_to}/{prefix}Diagnosis.png")


class CategoryEncoder(BaseEstimator):
    """A category variable transformer applied for whole DataFrame."""
    def __init__(self):
        self.categories_ = {}

    def _transformer(self, series: pd.Series, inverse=False):
        trans_dict = self.categories_.get(series.name, {})
        if not inverse:
            trans_dict = {v: k for k, v in trans_dict.items()}

        if trans_dict:
            return series.apply(lambda x: trans_dict.get(x, None))

        return series

    def fit(self, X: pd.DataFrame, y=None):
        _is_cat = X.dtypes == "object"
        self.categories_ = X.loc[:, _is_cat] .reset_index(drop=True) .apply(lambda x: x.drop_duplicates().reset_index(drop=True)) .to_dict()

        return self

    def transform(self, X: pd.DataFrame, copy=None):
        if copy:
            return X.copy().apply(self._transformer)
        return X.apply(self._transformer)

    def inverse_transform(self, X, copy=None):
        if copy:
            return X.copy().apply(self._transformer, inverse=True)
        return X.apply(self._transformer, inverse=True)

    def fit_transform(self, X, y=None):
        self.fit(X, y)
        return self.transform(X)


def setup(options):
    global RANDOM_STATE
    proj_dir = options.proj_dir
    seed = options.random_seed

    for sub_dir in ["Metacells", "Train", "Explain", "Predict"]:
        os.makedirs(f"{proj_dir}/{sub_dir}", exist_ok=True)

    torch.manual_seed(seed) # Pyro use Pytorch's random state.
    np.random.seed(seed)  # Scikit-learn uses Numpy's random state.
    if seed is not None: RANDOM_STATE = np.random.RandomState(seed)


def load_expression_matrix(fpath, as_train=True, min_pct=0.2, cell_types=None, keep_cell_type=True, test_ratio=None, features=None, ignore_cols=["ReSamplingIndex", "ClusterLabelPerGroup"], **kwarg):
    """Load expression matrix."""
    ct_col, id_col, lb_col = "CellType", "SampleID", "SampleLabel"
    res_cols = [lb_col, id_col, ct_col]

    exp_tab: pd.DataFrame = pd.read_csv(fpath, header=0, **kwarg).drop(labels=ignore_cols, axis=1, errors="ignore")

    # If there is a 'CellType' column, it allows the function to subset cell types for the downstream analysis.
    if ct_col in exp_tab.columns and cell_types:
        kept_recs = [x in cell_types for x in exp_tab.loc[:, ct_col]]
        exp_tab = exp_tab.loc[kept_recs, :]

    # If the min_pct argument is greater than 0 but less equal than 1, genes expressed less than min_pct will be discarded for downstream analyses.
    # FIXME Genes lowly expressed in samples A could be filtered out.
    # TODO Logics to filter category columns by counting missing values.
    if 0 <= min_pct <= 1 and as_train:
        n_cell, _ = exp_tab.shape
        cat_cols = exp_tab.columns[exp_tab.dtypes == "object"].to_list()
        cat_cols = [x for x in cat_cols if x not in res_cols]
        kept = exp_tab.drop(labels=res_cols+cat_cols, axis=1, errors="ignore").ne(0).sum(0).div(n_cell).ge(min_pct)
        kept_cols = res_cols + kept.index[kept].to_list() + cat_cols
        exp_tab = exp_tab.loc[:, kept_cols]

    # If there is a 'SampleID' column, it's a multi-sample dataset. However, the reality should also depend on the number of samples in the column.
    cts_map = None # cell to sample map.
    if id_col in exp_tab.columns:
        cts_map = exp_tab.loc[:, id_col]
        n_samples = cts_map.drop_duplicates().shape[0]
        if n_samples <= 1: # One or empty.
            cts_map = None

    # Select a subset of features. The logics also deal with missing features in unseen data, i.e., data to be predicted.
    if features and isinstance(features, list):
        exist_features = exp_tab.columns.to_list()

        # Assign 1 to features not existing in the exp_tab.
        mis_features = [x for x in features if x not in exist_features]
        exp_tab.loc[:, mis_features] = 1 # FIXME: a way to handle missing features

        # Keep features existing in the exp_tab.
        keep_features = res_cols + [x for x in features if x in exist_features] + mis_features

        exp_tab = exp_tab.loc[:, keep_features] # Ensure the order of input features

    if not keep_cell_type and exp_tab is not None:
        exp_tab = exp_tab.drop(labels=ct_col, axis=1, errors="ignore")

    # If there is a 'SampleLabel' column, it's a training dataset, otherwise, it's a dataset to be predicted. When 'SampleLabel' is in the dataset, `as_train=False` can be used to indicate the dataset is for training.
    xmat_tn, xmat_tt, yvec_tn, yvec_tt = [None] * 4
    if lb_col in exp_tab.columns:
        yvec = exp_tab.loc[:, lb_col]
        xmat = exp_tab.drop(labels=res_cols, axis=1, errors="ignore")

        if test_ratio is None or test_ratio <= 0 or test_ratio >= 1 or not as_train:
            xmat_tn, yvec_tn = xmat, yvec
        else:
            splits = train_test_split(xmat, yvec, test_size=test_ratio, stratify=yvec)
            xmat_tn, xmat_tt, yvec_tn, yvec_tt = splits
    else:
        xmat_tn = exp_tab

    return xmat_tn, xmat_tt, yvec_tn, yvec_tt, cts_map


def apply_reduction(mtx, non_feature_cols=["CellBarcodes"], reduc_mode="pca", n_pcs=50, logman: LogManager = LogManager("ApplyReduction")):
    if isinstance(non_feature_cols, str):
        non_feature_cols = [non_feature_cols]

    if reduc_mode == "ipca":
        model = IncrementalPCA(n_pcs)
    else:
        if reduc_mode != "pca":
            logman.info(f"Unknown reduction mode: {reduc_mode}. Use 'pca' instead.")
        model = PCA(n_pcs, random_state=RANDOM_STATE)

    feature_cols = [x for x in mtx.columns if x not in non_feature_cols]
    non_feature_mtx = mtx.loc[:, non_feature_cols].reset_index(drop=True)
    feature_mtx = mtx.loc[:, feature_cols].reset_index(drop=True)
    pca_mtx = pd.DataFrame(model.fit_transform(feature_mtx), columns=[f"PC_{x+1}" for x in range(n_pcs)])
    pca_mtx = pd.merge(non_feature_mtx, pca_mtx, left_index=True, right_index=True)

    return pca_mtx


def apply_clustering(reduc_mtx, non_reduc_cols, **kwargs):
    ap_model = AffinityPropagation(**kwargs).fit(reduc_mtx.drop(non_reduc_cols, axis=1, errors="ignore"))

    _cluster_labels = [f"cluster_{x:05}" for x in ap_model.labels_]

    if isinstance(non_reduc_cols, str): non_reduc_cols = [non_reduc_cols]
    _cluster_cols = non_reduc_cols + ["ClusterLabelPerGroup"]

    return reduc_mtx.assign(ClusterLabelPerGroup=_cluster_labels).loc[:, _cluster_cols]


def create_metacells(reduc_mtx, meta_mtx=None, group_by=[], barcode_col="CellBarcodes", max_iter=500, random_state=31415926, logman: LogManager = LogManager("CreateMetaCells")):
    if isinstance(random_state, int):
        random_state = np.random.RandomState(random_state)

    if meta_mtx is not None and group_by:
        if isinstance(group_by, str): group_by = [group_by]
        meta_mtx = meta_mtx.loc[:, [barcode_col] + group_by]
        merged_mtx = pd.merge(reduc_mtx, meta_mtx, on=barcode_col) if meta_mtx is not None else reduc_mtx
    else:
        logman.info("Create meta-cells without meta data.")
        group_by = ["FakeGroup"]
        merged_mtx = reduc_mtx.assign(FakeGroup="FG")

    non_pc_cols = [barcode_col] + group_by
    return merged_mtx.groupby(group_by).apply(apply_clustering, non_pc_cols, max_iter=max_iter, random_state=random_state).drop("FakeGroup", axis=1, errors="ignore")


def quantify_expression_by_resampling(expr_mtx, metacell_mtx, barcode_col="CellBarcodes", extra_group_by=[], n_samples=15, n_resamples=200, random_state=31415926):
    if isinstance(random_state, int):
        random_state = np.random.RandomState(random_state)
    seed_pool = random_state.randint(3141592654, size=n_resamples)

    # Resampling cells per metacells
    group_by_cols = ["ClusterLabelPerGroup"] + extra_group_by
    metacell_gp_mtx = metacell_mtx.groupby(group_by_cols, group_keys=False)
    resam_mtx = pd.concat([metacell_gp_mtx.apply(lambda x: x.sample(n=n_samples, random_state=np.random.RandomState(per_seed), replace=True)).assign(ReSamplingIndex=f"resampling_idx_{idx:06}") for idx, per_seed in enumerate(seed_pool)], axis=0, ignore_index=True)

    # Obtain average expression based on resampled cells per metacells
    merged_mtx = pd.merge(resam_mtx, expr_mtx, on=barcode_col)
    avexp_group_cols = group_by_cols + ["ReSamplingIndex"]
    avexp_mtx = merged_mtx.drop(barcode_col, axis=1, errors="ignore").groupby(avexp_group_cols).mean().reset_index()
    avexp_mtx = avexp_mtx.assign(MetacellBarcodes=[f"metacell_{x:06}" for x in range(avexp_mtx.shape[0])])
    resam_mtx = resam_mtx.loc[:, [barcode_col, "ReSamplingIndex"]]

    return avexp_mtx, resam_mtx


def load_model(fpath, fmt="pickle"):
    with open(fpath, "rb") as fhand:
        ext = os.path.splitext(fpath)
        if ext in ["pkl", "pickle"] or fmt in ["pkl", "pickle"]:
            return pkl.load(fhand)
        elif ext in ["jbl", "joblib"] or fmt in ["jbl", "joblib"]:
            return jbl.load(fpath)

    raise ValueError("Unknown model format.")


def parse_config(fpath):
    conf = None
    with open(fpath, "r") as fhandle:
        conf = json.load(fhandle)

    return conf["pos_label"], conf["param_space"]


def overwrite(force=False, *args):
    return not any([os.path.exists(x) for x in args]) or force


def ensemble_probs(pred_pc, save_to, **kwargs):
    pe_report = pd.DataFrame(columns=["SampleID", "ppval", "theta", "alpha", "beta"])

    for pid in pred_pc.loc[:, "SampleID"].unique():
        cur_probs = pred_pc.query(f"SampleID == '{pid}'").loc[:, "y_prob"]
        pe = PredEnsembler(cur_probs, **kwargs)
        pe.infer()
        pe.report(prefix=f"{pid}-", save_to=save_to)

        pe_report.loc[pid] = pd.Series(dict(SampleID=pid, ppval=pe.pp_val, theta=pe.theta, alpha=pe.alpha, beta=pe.beta))

    return pe_report.reset_index(drop=True)


def eval_model(xmat, y_true, model, pos_lab="Case", cts_map=None, save_to="./"):
    # Best meta-parameters by the RandomizedSearchCV()
    bst_par = model.best_params_
    with open(f"{save_to}/Best_params.json", "w") as bphandle:
        json.dump(bst_par, bphandle, indent=2)

    # Test the model
    y_pred = model.predict(xmat)
    pos_idx = model.classes_.tolist().index(pos_lab)
    y_prob = model.predict_proba(xmat)[:, pos_idx]

    # Save the selected features.
    xmat.columns.to_frame().to_csv(f"{save_to}/Selected_features.csv", index=False, header=["Feature"])

    # Save the predicted results.
    pred_pc = xmat.assign(SampleID=cts_map, y_true=y_true, y_pred=y_pred, y_prob=y_prob)
    pred_pc.to_csv(f"{save_to}/Prediction-percell.csv", index=False)

    # Infer expected probability per sample.
    pred_ps = ensemble_probs(pred_pc, save_to)
    pred_ps.to_csv(f"{save_to}/Prediction-persample.csv", index=False)

    # Evaluation matrices, plain text
    precision = precision_score(y_true, y_pred, pos_label=pos_lab)
    accuracy = accuracy_score(y_true, y_pred)
    roc_auc = roc_auc_score(y_true, y_prob)
    recall = recall_score(y_true, y_pred, pos_label=pos_lab)
    with open(f"{save_to}/Evaluation_metrics.txt", mode="w") as ofhand:
        ofhand.write(f"Precision: {precision:0.3f}\nAccuracy: {accuracy:0.3f}\nROC AUC: {roc_auc:0.3f}\nRecall: {recall:0.3f}\n")

    # ROC and PR curve
    fig, (roc_axe, prc_axe) = plt.subplots(ncols=2)

    ## ROC curve
    fpr, tpr, _ = roc_curve(y_true, y_prob, pos_label=pos_lab)
    roc_axe.plot(fpr, tpr)
    roc_axe.plot([0, 1], [0, 1], linestyle="dashed")
    roc_axe.text(0.75, 0.25, "{:.3f}".format(roc_auc))
    roc_axe.set_title("ROC curve (test)")
    roc_axe.set_xlabel("False positive rate")
    roc_axe.set_ylabel("True positive rate")
    roc_axe.set_xlim(-0.05, 1.05)
    roc_axe.set_ylim(-0.05, 1.05)

    ## PR curve
    pre, rec, _ = precision_recall_curve(y_true, y_prob, pos_label=pos_lab)
    prc_axe.plot(rec, pre)
    prc_axe.set_title("Precision-recall curve (test)")
    prc_axe.set_xlabel("Recall")
    prc_axe.set_ylabel("Precision")
    prc_axe.set_xlim(-0.05, 1.05)
    prc_axe.set_ylim(-0.05, 1.05)

    ## Save the figure
    fig.set_figwidth(12)
    fig.set_figheight(6)
    fig.set_tight_layout(True)
    fig.savefig(f"{save_to}/PR_and_ROC_curve.pdf")


def pipe_last_step(pipe, xmat):
    pbe = pipe.best_estimator_

    encoded = pbe["encode"].transform(xmat)
    scaled = pbe["scale"].transform(encoded)
    selected = pd.DataFrame(pbe["select"].transform(scaled), index=xmat.index, columns=xmat.columns[pbe["select"].get_support()])

    return pbe["classify"], selected


def metacell(options, logman: LogManager = LogManager("Metacell")):
    """Create meta cells and quantify expression by resampling."""
    logman.info("Create meta-cells.")

    # CLI options
    proj_dir = options.proj_dir
    expr_mtx_path = options.expr_mtx
    meta_mtx_path = options.meta_mtx
    reduc_mtx_path = options.reduc_mtx
    barcode_col = options.barcode_col
    n_rows = options.n_rows
    n_samples = options.n_samples
    n_resamples = options.n_resamples
    extra_group_by = options.extra_group_by
    out_prefix = options.out_prefix
    force = options.force

    # IO
    save_to = f"{proj_dir}/Metacells"

    metacell_mtx_saveto = f"{save_to}/{out_prefix}.metacells_percell.csv"
    avexp_mtx_saveto = f"{save_to}/{out_prefix}.expression_matrix_permetacell.csv"
    resam_mtx_saveto = f"{save_to}/{out_prefix}.resampling.csv"
    if not overwrite(force, metacell_mtx_saveto, avexp_mtx_saveto, resam_mtx_saveto):
        logman.error("Output files already exist, use -p/--out-prefix to distinguish them or -f/--force to overwrite them.")
        raise FileExistsError(f"Output files already exist: {metacell_mtx_saveto} or {avexp_mtx_saveto} or {resam_mtx_saveto}.")

    # Loading expression matrix and meta data (if any)
    expr_mtx = pd.read_csv(expr_mtx_path).head(n_rows) if n_rows else pd.read_csv(expr_mtx_path)
    meta_mtx = pd.read_csv(meta_mtx_path) if meta_mtx_path else None

    # Reduction matrix. If not given, create one from the expression matrix
    if reduc_mtx_path is None:
        reduc_mode = options.reduc_mode
        n_components = options.n_components
        reduc_mtx = apply_reduction(expr_mtx, barcode_col, reduc_mode, n_components)
    else:
        reduc_mtx = pd.read_csv(reduc_mtx_path)

    if n_rows is not None: reduc_mtx = reduc_mtx.head(n_rows)
    metacell_mtx = create_metacells(reduc_mtx, meta_mtx, group_by=extra_group_by, barcode_col=barcode_col, max_iter=500, random_state=RANDOM_STATE).reset_index(drop=True)
    if metacell_mtx is not None:
        metacell_mtx.to_csv(metacell_mtx_saveto, index=False)

    avexp_mtx, resam_mtx = quantify_expression_by_resampling(expr_mtx, metacell_mtx, barcode_col=barcode_col, extra_group_by=extra_group_by, n_samples=n_samples, n_resamples=n_resamples, random_state=RANDOM_STATE)
    avexp_mtx.to_csv(avexp_mtx_saveto, index=False)
    resam_mtx.to_csv(resam_mtx_saveto, index=False)

    logman.info(f"Done! Check {save_to} for the results.")


def train(options, logman: LogManager = LogManager("Train")):
    """Train a model"""
    logman.info("Train a model on given inputs.")

    # CLI Options
    proj_dir = options.proj_dir
    keep_cell_type = options.keep_cell_type
    model_arch = options.model_arch
    test_ratio = options.test_ratio
    cell_types = options.cell_types
    bcode_col = options.barcode_col
    cv_times = options.cv_times
    n_iters = options.n_iters
    n_jobs = options.n_jobs
    n_rows = options.n_rows

    # IO
    save_to = f"{proj_dir}/Train"
    in_file = f"{proj_dir}/Metacells/train.expression_matrix_permetacell.csv" if options.in_file is None else options.in_file

    # Parameters for the training
    pos_lab, param_space = parse_config(options.config)

    # Load expression matrix
    x_tn, x_tt, y_tn, y_tt, cts_map = load_expression_matrix(in_file, test_ratio=test_ratio, index_col=bcode_col, cell_types=cell_types, keep_cell_type=keep_cell_type, nrows=n_rows)

    # A pipeline for scaling data, selecting features, classifying samples.
    pipe_steps = [("encode", CategoryEncoder()), ("scale", StandardScaler()), ("select", SelectKBest(f_classif))]

    if model_arch == "nn":
        pipe_steps.append(("classify", NNClassifier()))
    elif model_arch == "rfc":
        pipe_steps.append(("classify", RandomForestClassifier(random_state=RANDOM_STATE)))
    else:
        if model_arch != "gbc":
            logman.warning(f"Unsupported {model_arch}, using gbc by default.")
        pipe_steps.append(("classify", GradientBoostingClassifier(random_state=RANDOM_STATE)))

    pipe = Pipeline(steps=pipe_steps)

    # Searching the best hyper-parameters randomly.
    rscv = RandomizedSearchCV(pipe, param_space, n_iter=n_iters, n_jobs=n_jobs, cv=cv_times)
    rscv.fit(x_tn, y_tn)

    # Save the model
    with open(f"{save_to}/Model.pickle", "bw") as mf_handle:
        pkl.dump(rscv, mf_handle)

    # Evaluate the model
    eval_model(x_tt, y_tt, rscv, pos_lab, cts_map, save_to)

    logman.info(f"Done! Check {save_to} for the results.")


def explain(options, logman: LogManager = LogManager("Explain")):
    logman.info("Explain the model by SHAP values.")

    # CLI options
    cell_types = options.cell_types
    barcode_col = options.barcode_col
    proj_dir = options.proj_dir
    n_rows = options.n_rows

    # IO
    save_to = f"{proj_dir}/Explain"
    tar_features = f"{proj_dir}/Train/Selected_features.csv"
    in_file = f"{proj_dir}/Metacells/train.expression_matrix_permetacell.csv" if options.in_file is None else options.in_file

    # Load model and expression matrix
    model = load_model(f"{proj_dir}/Train/Model.pickle")
    xmat, _, _, _, _ = load_expression_matrix(in_file, as_train=False, index_col=barcode_col, cell_types=cell_types, nrows=n_rows)
    features = pd.read_csv(tar_features, header=0).loc[:, "Feature"].to_list()
    xmat = xmat.loc[:, features]

    # Plot SHAP values. Including bar and dot (beeswarm) plot
    pbe, xmat = pipe_last_step(model, xmat)
    explainer = sp.TreeExplainer(pbe)
    shap_vals = explainer(xmat)

    # Save SHAP values
    with open(f"{save_to}/SHAP_values.pickle", "bw") as sf_handle:
        pkl.dump(shap_vals, sf_handle)

    # Plots
    for ptype in ["dot", "bar"]:
        plt.clf() # We clear the figure.

        sp.summary_plot(shap_vals, plot_type=ptype, show=False)
        fig = plt.gcf()
        if ptype == "dot":
            _, axe_cb = fig.get_axes()
            axe_cb.set_visible(False)
            axe_cb = plt.colorbar()

        fig.set_figwidth(7)
        fig.set_figheight(7)
        fig.set_tight_layout(True)
        fig.savefig(f"{save_to}/{ptype.title()}_plot.pdf")

    logman.info(f"Done! Check {save_to} for the results.")


def predict(options, logman: LogManager = LogManager("Predict")):
    logman.info("Predict unseen samples.")

    # CLI options
    in_file = options.in_file
    cell_types = options.cell_types
    barcode_col = options.barcode_col
    pos_label = options.pos_label
    proj_dir = options.proj_dir
    n_iters = options.n_iters
    n_rows = options.n_rows
    out_subdir = options.out_subdir

    # Create sub-dir to store the predicted results.
    time_stamp = dt.now().strftime("%Y%m%d%H%M%S")
    if out_subdir is None:
        out_subdir = f"{proj_dir}/Predict/{time_stamp}"
    os.makedirs(out_subdir, exist_ok=True)

    # Load features used in the model
    tar_features = f"{proj_dir}/Train/Selected_features.csv"
    features = pd.read_csv(tar_features, header=0).loc[:, "Feature"].to_list()

    # Load model and expression matrix
    model_path = f"{proj_dir}/Train/Model.pickle"
    model = load_model(model_path)
    xmat, _, _, _, cts_map = load_expression_matrix(in_file, as_train=False, cell_types=cell_types, min_pct=0, features=features, index_col=barcode_col, nrows=n_rows, keep_cell_type=False)

    # Prediction
    y_pred = model.predict(xmat)
    y_prob = model.predict_proba(xmat)[:, pos_label]

    # Save predicted results
    pred_pc = xmat.assign(SampleID=cts_map, y_prob=y_prob, y_pred=y_pred)
    pred_pc.to_csv(f"{out_subdir}/Prediction-percell.csv")

    # Infer expected probability per sample
    pred_ps = ensemble_probs(pred_pc, out_subdir, n_iters=n_iters)
    pred_ps.to_csv(f"{out_subdir}/Prediction-persample.csv", index=False)

    logman.info(f"Done! Check {out_subdir} for the results")


def get_cli_opts():
    """Get command line options"""
    par = argparse.ArgumentParser()
    par.add_argument("-s", "--random-seed", default=31415, type=int, help="Random seed. Default: %(default)s")
    par.add_argument("-P", "--proj-dir", default="Scpred", help="The project dir containing train, explain, and predict results. Default: %(default)s")

    subpar = par.add_subparsers(dest="subcmd", required=True)

    psc_par = subpar.add_parser("metacell", help="Create expression matrix of pseudo-cells by resampling.")
    psc_par.add_argument("-e", "--expr-mtx", required=True, help="Expression matrix. Required")
    psc_par.add_argument("-m", "--meta-mtx", default=None, help="Meta data matrix.")
    psc_par.add_argument("-r", "--reduc-mtx", default=None, help="Reduction matrix, e.g., PCA.")
    psc_par.add_argument("-b", "--barcode-col", default="CellBarcodes", help="The column used as the index column, i.e., cell barcodes. Default: %(default)s")
    psc_par.add_argument("-g", "--extra-group-by", default=None, nargs="*", help="Extra group by columns. Default: %(default)s")
    psc_par.add_argument("-f", "--force", action="store_true", help="Force to overwrite existing files.")
    psc_par.add_argument("-p", "--out-prefix", default="train", help="Prefix of the output files. Default: %(default)s")
    psc_par.add_argument("--n-resamples", default=200, type=int, help="Number of re-samples. Default: %(default)s")
    psc_par.add_argument("--n-samples", default=15, type=int, help="Number of samples. Default: %(default)s")
    psc_par.add_argument("--n-rows", default=None, type=int, help="Number of rows to read. If not specified all rows will be loaded, useful for testing. Default: all rows")

    trn_par = subpar.add_parser("train", help="Train a model from expression data.")
    trn_par.add_argument("-c", "--config", required=True, help="Configuration file. Required")
    trn_par.add_argument("-i", "--in-file", default=None, help="Input file for training. Default: %(default)s")
    trn_par.add_argument("-b", "--barcode-col", default="MetacellBarcodes", help="The column used as the index column, i.e., cell barcodes. Default: %(default)s")
    trn_par.add_argument("-p", "--test-ratio", default=0.3, type=float, help="Ratio of data used for test. Default: %(default)s")
    trn_par.add_argument("-I", "--n-iters", default=15, type=int, help="Number of iterations used for the RandomizedSearchCV. Default: %(default)s")
    trn_par.add_argument("-x", "--cv-times", default=10, type=int, help="Number of cross validation. Default: %(default)s")
    trn_par.add_argument("-t", "--cell-types", nargs="*", default=None, help="Cell types on which the model will be trained. Default: all")
    trn_par.add_argument("-k", "--keep-cell-type", action="store_true", default=False, help="Whether keep cell types as a predictor variable. Default: %(default)s")
    trn_par.add_argument("-m", "--model-arch", choices=["gbc", "rfc", "nn"], default="gbc", help="Which model architecture to be used, i.e., 'gbc', 'rfc', or 'nn'. Default: %(default)s")
    trn_par.add_argument("-@", "--n-jobs", default=1, type=int, help="Number of jobs to be run in parallel. Default: %(default)s")
    trn_par.add_argument("--force-features", default=None, nargs="*", help="Features should be included by force. Default: %(default)s")
    trn_par.add_argument("--n-rows", default=None, type=int, help="Number of rows to read. If not specified all rows will be loaded, useful for testing. Default: all rows")

    exp_par = subpar.add_parser("explain", help="Explain a model by SHAP values.")
    exp_par.add_argument("-i", "--in-file", default=None, help="The input data used for training. Default: %(default)s")
    exp_par.add_argument("-b", "--barcode-col", default="MetacellBarcodes", help="The column used as the index column, i.e., cell barcodes. Default: %(default)s")
    exp_par.add_argument("-t", "--cell-types", nargs="*", default=None, help="Cell types on which the model will be trained. Default: all")
    exp_par.add_argument("--n-rows", default=2000, type=int, help="Number of rows to read. If not specified all rows will be loaded. Default: all rows")

    prd_par = subpar.add_parser("predict", help="Predict unseen expression data.")
    prd_par.add_argument("-i", "--in-file", required=True, help="Input file to be predict. Required.")
    prd_par.add_argument("-p", "--pos-label", default=1, type=int, help="The positive label. Default: %(default)s")
    prd_par.add_argument("-b", "--barcode-col", default="MetacellBarcodes", help="The column used as the index column, i.e., cell barcodes. Default: %(default)s")
    prd_par.add_argument("-t", "--cell-types", nargs="*", default=None, help="Cell types on which the model will be trained. Default: all")
    prd_par.add_argument("-n", "--n-rows", default=None, type=int, help="Number of rows to be read from the training in-file. Default: all")
    prd_par.add_argument("-o", "--out-subdir", default=None, help="The file name prefix used to save prediction results.")
    prd_par.add_argument("--n-iters", default=5000, type=int, help="Number interations of resampling to estimate the sample-level prediction. Default: %(default)s")

    return par.parse_args()


def main(logman: LogManager = LogManager("Scpred")):
    # Get the CLI options
    options = get_cli_opts()

    # Set up project dir, random state, etc.
    setup(options)

    # Sub-command
    subcmd = options.subcmd
    if subcmd == "metacell":
        metacell(options)
    elif subcmd == "train":
        train(options)
    elif subcmd == "explain":
        explain(options)
    elif subcmd == "predict":
        predict(options)
    else:
        logman.error("Unknown sub-command")


if __name__ == "__main__":
    main()
