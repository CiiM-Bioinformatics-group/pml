#!/usr/bin/env python
import os

import pickle as pkl

import matplotlib as mpl
import seaborn as sns
import numpy as np
from scipy.spatial.distance import cdist
from scipy.stats import spearmanr, pearsonr, kendalltau, false_discovery_control

import shap as sp
import pandas as pd
import joblib as jbl
import matplotlib.pyplot as plt
from matplotlib import colormaps

from sklearn.base import BaseEstimator


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


def pipe_last_step(pipe, xmat):
    pbe = pipe.best_estimator_

    encoded = pbe["encode"].transform(xmat)
    scaled = pbe["scale"].transform(encoded)
    selected = pd.DataFrame(pbe["select"].transform(scaled), index=xmat.index, columns=xmat.columns[pbe["select"].get_support()])

    return pbe["classify"], selected


def load_model(fpath, fmt="pickle"):
    with open(fpath, "rb") as fhand:
        ext = os.path.splitext(fpath)
        if ext in ["pkl", "pickle"] or fmt in ["pkl", "pickle"]:
            return pkl.load(fhand)
        elif ext in ["jbl", "joblib"] or fmt in ["jbl", "joblib"]:
            return jbl.load(fpath)

    raise ValueError("Unknown model format.")


def plot_shap_vals(shap_vals, feature=None, dep_feature=None, instance_idx=None, save_to=None, fig_width=7.0, fig_height=6.0, **kwargs):
    draw_sumplot = feature is None and dep_feature is None and instance_idx is None
    draw_scatter = feature is not None and dep_feature is None and instance_idx is None
    draw_waterfall = feature is None and dep_feature is None and instance_idx is not None
    draw_interaction = feature is not None and dep_feature is not None and instance_idx is None

    plt.clf() # Clear the figure.
    save_token = "sumplot"
    if draw_sumplot:
        sp.summary_plot(shap_vals, plot_type="dot", show=False, **kwargs)
    if draw_scatter:
        save_token = f"{feature}.scatter"
        sp.plots.scatter(shap_vals[:, feature], color=shap_vals[:, feature], show=False, **kwargs)
    if draw_waterfall:
        save_token = f"sample_{instance_idx}.waterfall"
        sp.plots.waterfall(shap_vals[instance_idx], show=False, **kwargs)
    if draw_interaction:
        save_token = f"{feature}_by_{dep_feature}.interaction"
        sp.plots.scatter(shap_vals[:, dep_feature], color=shap_vals[:, feature], show=False, **kwargs)

    if save_to is None:
        save_to = f"SHAP_values.{save_token}.pdf"

    fig = plt.gcf()
    fig.set_figwidth(fig_width)
    fig.set_figheight(fig_height)
    fig.set_tight_layout(True)
    fig.savefig(save_to, dpi=300)
    plt.close(fig)


omics_type = "cite_seq"
work_dir = "/home/zzhang/Documents/projects/wp_pml/outputs/analysis/prediction"


# UMAP of metacells
group_cols = ["SampleID", "SampleLabel", "CellType", "ClusterLabelPerGroup"]
all_mc_umap_mtx, all_rc_umap_mtx = pd.DataFrame(), pd.DataFrame()
for pct in ["cd4_t", "cd8_t", "monocyte", "nk_cell", "b_cell"]:
    # features = pd.read_csv(f"{work_dir}/models/{pct}/Train/Selected_features.csv", header=0).loc[:, "Feature"].to_list()
    # rsc_mtx = pd.read_csv(f"{work_dir}/models/{pct}/Metacells/train.resampling.csv")
    mcc_mtx = pd.read_csv(f"{work_dir}/models/{pct}/Metacells/train.metacells_percell.csv")
    umap_mtx = pd.read_csv(f"{work_dir}/datasets/{omics_type}/pbmc.cite_seq.integrated.reductions_umap.csv")
    rc_umap_mtx = pd.merge(left=umap_mtx, right=mcc_mtx, on="CellBarcodes", how="inner")
    mc_umap_mtx = rc_umap_mtx.groupby(group_cols).apply(lambda x: x.loc[:, ["wnnUMAP_1", "wnnUMAP_2"]].mean()).reset_index()
    all_mc_umap_mtx = pd.concat([all_mc_umap_mtx, mc_umap_mtx], axis=0).reset_index(drop=True)
    all_rc_umap_mtx = pd.concat([all_rc_umap_mtx, rc_umap_mtx], axis=0).reset_index(drop=True)

response_order = ["Responder", "Non-responder"]
g = sns.JointGrid(data=all_rc_umap_mtx, x="wnnUMAP_1", y="wnnUMAP_2", ratio=5, space=0.2)
_ = g.plot_joint(sns.scatterplot, s=7, hue="CellType", alpha=0.5, edgecolor="none", data=all_rc_umap_mtx)
_ = g.ax_joint.scatter(x="wnnUMAP_1", y="wnnUMAP_2", c="black", s=12, edgecolors="none", alpha=0.75, data=all_mc_umap_mtx)
handles, labels = [], []
for per_handle, per_label in zip(*g.ax_joint.get_legend_handles_labels()):
    per_handle.set_alpha(1)
    if isinstance(per_handle, mpl.lines.Line2D): per_handle.set_markersize(5)
    if per_label.startswith("wnnUMAP"): per_label = "Metacell"
    handles.append(per_handle)
    labels.append(per_label)
_ = g.ax_joint.legend(handles, labels, loc="upper left", title="Cell Type", title_fontsize="small", fontsize="x-small")
_ = sns.kdeplot(x="wnnUMAP_1", hue="SampleLabel", multiple="stack", palette="Pastel1", data=all_rc_umap_mtx, ax=g.ax_marg_x)
_ = g.ax_marg_x.legend(loc="upper left", labels=response_order, title=None, fontsize="x-small")
_ = sns.kdeplot(y="wnnUMAP_2", hue="CellType", multiple="stack", legend=False, data=all_rc_umap_mtx, ax=g.ax_marg_y)
g.figure.set_size_inches(5.5, 5.5)
g.savefig(f"{work_dir}/plots/metacells.BL.all_cell_type.umap.pdf", bbox_inches="tight")


# Summary plot of top 30 important features
cell_type = "nk_cell"
in_file = f"{work_dir}/models/{cell_type}/Metacells/train.expression_matrix_permetacell.csv"
features = pd.read_csv(f"{work_dir}/models/{cell_type}/Train/Selected_features.csv", header=0).loc[:, "Feature"].to_list()
raw_xmat, _, _, _, cts_map = load_expression_matrix(in_file, as_train=False, index_col="MetacellBarcodes", cell_types=None)
raw_xmat = raw_xmat.loc[:, features]

# Load model and expression matrix
model = load_model(f"{work_dir}/models/{cell_type}/Train/Model.pickle")

# Calculate SHAP values
pbe, xmat = pipe_last_step(model, raw_xmat) # Obtain the last step of the sklearn pipeline and the corresponding input data.
explainer = sp.Explainer(pbe, xmat) # Build a SHAP explainer.
shap_vals = explainer(xmat) # Compute SHAP values.

# Potential interactions based on SHAP values
# inds = sp.utils.potential_interactions(shap_vals[:, "CD2"], shap_vals)
# dep_feature = shap_vals.feature_names[inds[1]]

# Plot SHAP values based on given model and input matrix
plot_shap_vals(shap_vals, max_display=30, save_to=f"{work_dir}/models/{cell_type}/Explain/SHAP_values.sumplot.pdf")

# Scatter plots of SHAP values
if cell_type == "cd8_t":
    plot_shap_vals(shap_vals, feature="CD2", save_to=f"{work_dir}/models/{cell_type}/Explain/SHAP_values.scatter.CD2.pdf")
    plot_shap_vals(shap_vals, feature="DDX3Y", save_to=f"{work_dir}/models/{cell_type}/Explain/SHAP_values.scatter.DDX3Y.pdf")
    plot_shap_vals(shap_vals, feature="ITGB1", save_to=f"{work_dir}/models/{cell_type}/Explain/SHAP_values.scatter.ITGB1.pdf")
    plot_shap_vals(shap_vals, feature="CD226", save_to=f"{work_dir}/models/{cell_type}/Explain/SHAP_values.scatter.CD226.pdf")
    plot_shap_vals(shap_vals, feature="S100A4", save_to=f"{work_dir}/models/{cell_type}/Explain/SHAP_values.scatter.S100A4.pdf")
elif cell_type == "cd4_t":
    plot_shap_vals(shap_vals, feature="CD2", save_to=f"{work_dir}/models/{cell_type}/Explain/SHAP_values.scatter.CD2.pdf")
    plot_shap_vals(shap_vals, feature="IL6R", save_to=f"{work_dir}/models/{cell_type}/Explain/SHAP_values.scatter.IL6R.pdf")
    plot_shap_vals(shap_vals, feature="FOS", save_to=f"{work_dir}/models/{cell_type}/Explain/SHAP_values.scatter.FOS.pdf")
    plot_shap_vals(shap_vals, feature="LPAR6", save_to=f"{work_dir}/models/{cell_type}/Explain/SHAP_values.scatter.LPAR6.pdf")
elif cell_type == "monocyte":
    plot_shap_vals(shap_vals, feature="HLA-DQB1", save_to=f"{work_dir}/models/{cell_type}/Explain/SHAP_values.scatter.HLA-DQB1.pdf")
    plot_shap_vals(shap_vals, feature="CLEC4E", save_to=f"{work_dir}/models/{cell_type}/Explain/SHAP_values.scatter.CLEC4E.pdf")
    plot_shap_vals(shap_vals, feature="TPT1", save_to=f"{work_dir}/models/{cell_type}/Explain/SHAP_values.scatter.TPT1.pdf")
    plot_shap_vals(shap_vals, feature="DSE", save_to=f"{work_dir}/models/{cell_type}/Explain/SHAP_values.scatter.DSE.pdf")
    plot_shap_vals(shap_vals, feature="SELL", save_to=f"{work_dir}/models/{cell_type}/Explain/SHAP_values.scatter.SELL.pdf")
elif cell_type == "nk_cell":
    plot_shap_vals(shap_vals, feature="IFITM2", save_to=f"{work_dir}/models/{cell_type}/Explain/SHAP_values.scatter.IFITM2.pdf")


# Interaction plots of SHAP values
if cell_type == "cd8_t":
    plot_shap_vals(shap_vals, feature="S100A4", dep_feature="CD2", save_to=f"{work_dir}/models/{cell_type}/Explain/SHAP_values.interaction.S100A4_by_CD2.pdf")
    plot_shap_vals(shap_vals, feature="ITGB1", dep_feature="CD2", save_to=f"{work_dir}/models/{cell_type}/Explain/SHAP_values.interaction.ITGB1_by_CD2.pdf")
    plot_shap_vals(shap_vals, feature="CD226", dep_feature="CD2", save_to=f"{work_dir}/models/{cell_type}/Explain/SHAP_values.interaction.CD226_by_CD2.pdf")
    plot_shap_vals(shap_vals, feature="CD226", dep_feature="NFKBIA", save_to=f"{work_dir}/models/{cell_type}/Explain/SHAP_values.interaction.CD226_by_NFKBIA.pdf")
    plot_shap_vals(shap_vals, feature="CD2", dep_feature="NFKBIA", save_to=f"{work_dir}/models/{cell_type}/Explain/SHAP_values.interaction.CD2_by_NFKBIA.pdf")
    plot_shap_vals(shap_vals, feature="CD2", dep_feature="THEMIS", save_to=f"{work_dir}/models/{cell_type}/Explain/SHAP_values.interaction.CD2_by_THEMIS.pdf")


spearman_r_mtx = pd.DataFrame()
for pct in ["monocyte", "cd4_t", "cd8_t", "b_cell", "nk_cell"]:
    in_file = f"{work_dir}/models/{pct}/Metacells/train.expression_matrix_permetacell.csv"
    features = pd.read_csv(f"{work_dir}/models/{pct}/Train/Selected_features.csv", header=0).loc[:, "Feature"].to_list()
    raw_xmat, _, _, _, cts_map = load_expression_matrix(in_file, as_train=False, index_col="MetacellBarcodes", cell_types=None)
    raw_xmat = raw_xmat.loc[:, features]

    # Load model and expression matrix
    model = load_model(f"{work_dir}/models/{pct}/Train/Model.pickle")

    # Calculate SHAP values
    pbe, xmat = pipe_last_step(model, raw_xmat) # Obtain the last step of the sklearn pipeline and the corresponding input data.
    explainer = sp.Explainer(pbe, xmat) # Build a SHAP explainer.
    shap_vals = explainer(xmat) # Compute SHAP values.

    # Plot SHAP values based on given model and input matrix
    plot_shap_vals(shap_vals, max_display=30, save_to=f"{work_dir}/models/{pct}/Explain/SHAP_values.sumplot.pdf", fig_width=4.5)

    # Features that are highly correlated with the SHAP values
    nonconst_feature_idx = [i for i in range(shap_vals.data.shape[1]) if np.std(shap_vals.values[:, i]) != 0]
    tmp_mtx = pd.DataFrame([
        [shap_vals.feature_names[pid], *spearmanr(shap_vals.data[:, pid], shap_vals.values[:, pid])]
        for pid in nonconst_feature_idx
    ], columns=["Feature", "SpearmanRho", "SpearmanPvalue"]).assign(CellType=pct).sort_values("SpearmanPvalue", ignore_index=True)
    tmp_mtx = tmp_mtx.assign(SpearmanFDR=false_discovery_control(tmp_mtx.loc[:, "SpearmanPvalue"].to_list()))
    spearman_r_mtx = pd.concat([spearman_r_mtx, tmp_mtx], ignore_index=True)

spearman_r_mtx.to_csv(f"{work_dir}/function_analysis/correlation_between_shap_and_expression.csv", index=False)
