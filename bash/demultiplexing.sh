#!/usr/bin/env bash
# File: dm_atac.sh
# Author: Zhenhua Zhang
# E-mail: zhenhua.zhang217@gmail.com
# Created: Oct 06, 2023
# Updated: Oct 06, 2023

# A pipe line to dmultiplex single-cell ATAC-seq data, adopted from Souporcell pipeline
projdir=~/Documents/projects/wp_pml

#
## Demultiplexing scATAC-seq by Vireo
#
sbatch --parsable --mem 4G --time 7:59:00 --array 1-4 --cpus-per-task 12 --job-name scATAC-seq.demultiplex \
  --output $projdir/outputs/demultiplex/logs/%A_%a-%u.%x.out <<'EOF'
#!/usr/bin/bash
# set -Eeuo pipefail

n_cpus=${SLURM_CPUS_PER_TASK:-1}
pool_id=${SLURM_ARRAY_TASK_ID:-1}

declare -A samples_in_pool=(
  ["1"]="PML0045,PML0061,PML0004"
  ["2"]="PML0055,PML0058,PML0066"
  ["3"]="PML0007,PML0033,PML0070"
  ["4"]="PML0013,PML0060,PML0063"
)
target_pool=ATACB2P${pool_id}
n_cluster=$(echo ${samples_in_pool[$pool_id]} | sed 's/,/\n/g' | wc -l)

projdir=~/Documents/projects/wp_pml
known_genotypes=${projdir}/outputs/genotypes/batch_2/postimputation/quality_control/all.GRCh38.imputed.snps.biallelic.geno_3e-1.maf_1e-2.annotated_b150.PML_samples.vcf.gz

in_dir=${projdir}/outputs/readcounts/ATAC_seq/${target_pool}
out_dir=${projdir}/outputs/demultiplex/ATAC_seq/${target_pool}

mkdir -p ${out_dir}

# Call genotypes per cell
echo "[I]: Calling genotypes per cell"
singularity exec -B ${projdir} ${projdir}/scripts/tools/dmppl.sif /opt/tools/bin/cellsnp-lite \
  -s ${in_dir}/outs/possorted_bam.bam \
  -b ${in_dir}/outs/filtered_peak_bc_matrix/barcodes.tsv \
  -R ${known_genotypes} \
  -O ${out_dir}/cellsnp_outdir \
  --UMItag None \
  --minMAF 0.05 \
  --minCOUNT 15 \
  -p ${n_cpus}

# Compress cellSNP.base.vcf
singularity exec -B ${projdir} ${projdir}/scripts/tools/dmppl.sif bcftools view \
  -O z \
  -o ${out_dir}/cellsnp_outdir/cellSNP.base.vcf.gz \
  ${out_dir}/cellsnp_outdir/cellSNP.base.vcf

# Demultiplexing supposing donors' genotypes are unknown
echo "[I]: Demultiplexing supposing donors' genotypes are unknown"
singularity exec -B ${projdir} ${projdir}/scripts/tools/dmppl.sif vireo \
  -N ${n_cluster} \
  -c ${out_dir}/cellsnp_outdir \
  -o ${out_dir}/vireo_outdir_denovo

# Demultiplexing using donors' genotypes
echo "[I]: Demultiplexing using donors' genotypes"
mkdir -p $out_dir/known_genotypes
singularity exec -B ${projdir} ${projdir}/scripts/tools/dmppl.sif bcftools view \
  -s ${samples_in_pool[$pool_id]} \
  -R ${out_dir}/cellsnp_outdir/cellSNP.base.vcf.gz \
  -O z \
  -o $out_dir/known_genotypes/known_genotypes.vcf.gz \
  $known_genotypes
singularity exec -B ${projdir} ${projdir}/scripts/tools/dmppl.sif bcftools index -t ${out_dir}/known_genotypes/known_genotypes.vcf.gz
singularity exec -B ${projdir} ${projdir}/scripts/tools/dmppl.sif vireo \
  -t GP \
  -d ${out_dir}/known_genotypes/known_genotypes.vcf.gz \
  -c ${out_dir}/cellsnp_outdir \
  -o ${out_dir}/vireo_outdir_ref
EOF


# DM scRNA-seq data
sbatch --parsable --mem 4G --time 7:59:00 --array 1-4 --cpus-per-task 12 --job-name scRNA-seq.demultiplex \
  --output $projdir/outputs/demultiplex/logs/%A_%a-%u.%x.out <<'EOF'
#!/usr/bin/bash
# set -Eeuo pipefail

n_cpus=${SLURM_CPUS_PER_TASK:-1}
pool_id=${SLURM_ARRAY_TASK_ID:-1}

declare -A samples_in_pool=(
  ["1"]="PML0045,PML0061,PML0004"
  ["2"]="PML0055,PML0058,PML0066"
  ["3"]="PML0007,PML0033,PML0070"
  ["4"]="PML0013,PML0060,PML0063"
)
target_pool=RNAB2P${pool_id}
n_cluster=$(echo ${samples_in_pool[$pool_id]} | sed 's/,/\n/g' | wc -l)

projdir=~/Documents/projects/wp_pml
known_genotypes=${projdir}/outputs/genotypes/batch_2/postimputation/quality_control/all.GRCh38.imputed.snps.biallelic.geno_3e-1.maf_1e-2.annotated_b150.PML_samples.vcf.gz

in_dir=${projdir}/outputs/readcounts/RNA_seq/${target_pool}

out_dir=${projdir}/outputs/demultiplex/RNA_seq/${target_pool}
mkdir -p ${out_dir}

# Call genotypes per cell
echo "[I]: Calling genotypes per cell"
singularity exec -B ${projdir} ${projdir}/scripts/tools/dmppl.sif /opt/tools/bin/cellsnp-lite \
  -s ${in_dir}/outs/possorted_genome_bam.bam \
  -b ${in_dir}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz \
  -R ${known_genotypes} \
  -O ${out_dir}/cellsnp_outdir \
  --minMAF 0.05 \
  --minCOUNT 15 \
  -p ${n_cpus}

# Compress cellSNP.base.vcf
singularity exec -B ${projdir} ${projdir}/scripts/tools/dmppl.sif bcftools view \
  -O z \
  -o ${out_dir}/cellsnp_outdir/cellSNP.base.vcf.gz \
  ${out_dir}/cellsnp_outdir/cellSNP.base.vcf

singularity exec -B ${projdir} ${projdir}/scripts/tools/dmppl.sif bcftools index \
  -ft ${out_dir}/cellsnp_outdir/cellSNP.base.vcf.gz

# Demultiplexing supposing donors' genotypes are unknown
echo "[I]: Demultiplexing supposing donors' genotypes are unknown"
singularity exec -B ${projdir} ${projdir}/scripts/tools/dmppl.sif vireo \
  -N ${n_cluster} \
  -c ${out_dir}/cellsnp_outdir \
  -o ${out_dir}/vireo_outdir_denovo

# Demultiplexing using donors' genotypes
echo "[I]: Demultiplexing using donors' genotypes"
mkdir -p $out_dir/known_genotypes
singularity exec -B ${projdir} ${projdir}/scripts/tools/dmppl.sif bcftools view \
  -s ${samples_in_pool[$pool_id]} \
  -R ${out_dir}/cellsnp_outdir/cellSNP.base.vcf.gz \
  -O z \
  -o $out_dir/known_genotypes/known_genotypes.vcf.gz \
  $known_genotypes

singularity exec -B ${projdir} ${projdir}/scripts/tools/dmppl.sif bcftools index \
  -ft ${out_dir}/known_genotypes/known_genotypes.vcf.gz

singularity exec -B ${projdir} ${projdir}/scripts/tools/dmppl.sif vireo \
  -t GP \
  -d ${out_dir}/known_genotypes/known_genotypes.vcf.gz \
  -c ${out_dir}/cellsnp_outdir \
  -o ${out_dir}/vireo_outdir_ref
EOF


# DM CITE-seq data
sbatch --parsable --mem 4G --time 7:59:00 --array 1-2 --cpus-per-task 12 --job-name scCITE-seq.demultiplex \
  --output $projdir/outputs/demultiplex/logs/%A_%a-%u.%x.out <<'EOF'
#!/usr/bin/bash
# set -Eeuo pipefail

n_cpus=${SLURM_CPUS_PER_TASK:-1}
pool_id=${SLURM_ARRAY_TASK_ID:-1}

declare -A samples_in_pool=(
  ["1"]="PML16,PML17,PML22"
  ["2"]="PML2,PML8,PML9,PML22"
  ["3"]="PML8,PML9,PML22,PML35"
  ["4"]="PML2,PML8,PML17,PML25"
)
target_pool=CITEpool${pool_id}
n_cluster=$(echo ${samples_in_pool[$pool_id]} | sed 's/,/\n/g' | wc -l)

projdir=~/Documents/projects/wp_pml
known_genotypes=${projdir}/outputs/genotypes/batch_2/postimputation/quality_control/all.GRCh38.imputed.snps.biallelic.geno_3e-1.maf_1e-2.annotated_b150.PML_samples.vcf.gz

in_dir=${projdir}/outputs/readcounts/CITE_seq/${target_pool}

out_dir=${projdir}/outputs/demultiplex/CITE_seq/${target_pool}
mkdir -p ${out_dir}

# Call genotypes per cell
echo "[I]: Calling genotypes per cell"
singularity exec -B ${projdir} ${projdir}/scripts/tools/dmppl.sif /opt/tools/bin/cellsnp-lite \
  -s ${in_dir}/outs/possorted_genome_bam.bam \
  -b ${in_dir}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz \
  -R ${known_genotypes} \
  -O ${out_dir}/cellsnp_outdir \
  --minMAF 0.05 \
  --minCOUNT 15 \
  -p ${n_cpus}

# Compress cellSNP.base.vcf
singularity exec -B ${projdir} ${projdir}/scripts/tools/dmppl.sif bcftools view \
  -O z \
  -o ${out_dir}/cellsnp_outdir/cellSNP.base.vcf.gz \
  ${out_dir}/cellsnp_outdir/cellSNP.base.vcf

# Demultiplexing supposing donors' genotypes are unknown
echo "[I]: Demultiplexing supposing donors' genotypes are unknown"
singularity exec -B ${projdir} ${projdir}/scripts/tools/dmppl.sif vireo \
  -N ${n_cluster} \
  -c ${out_dir}/cellsnp_outdir \
  -o ${out_dir}/vireo_outdir_denovo

# Demultiplexing using donors' genotypes
echo "[I]: Demultiplexing using donors' genotypes"
mkdir -p $out_dir/known_genotypes
singularity exec -B ${projdir} ${projdir}/scripts/tools/dmppl.sif bcftools index ${out_dir}/cellsnp_outdir/cellSNP.base.vcf.gz
singularity exec -B ${projdir} ${projdir}/scripts/tools/dmppl.sif bcftools view \
  -s ${samples_in_pool[$pool_id]} \
  -R ${out_dir}/cellsnp_outdir/cellSNP.base.vcf.gz \
  -O z \
  -o $out_dir/known_genotypes/known_genotypes.vcf.gz \
  $known_genotypes
singularity exec -B ${projdir} ${projdir}/scripts/tools/dmppl.sif bcftools index -t ${out_dir}/known_genotypes/known_genotypes.vcf.gz

singularity exec -B ${projdir} ${projdir}/scripts/tools/dmppl.sif vireo \
  -t GP \
  -d ${out_dir}/known_genotypes/known_genotypes.vcf.gz \
  -c ${out_dir}/cellsnp_outdir \
  -o ${out_dir}/vireo_outdir_ref
EOF


# DM multiome data
sbatch --parsable --mem 4G --time 7:59:00 --array 3-4 --cpus-per-task 12 --job-name Multiome-gex.demultiplex \
  --output $projdir/outputs/demultiplex/logs/%A_%a-%u.%x.out <<'EOF'
#!/usr/bin/bash
# set -Eeuo pipefail

n_cpus=${SLURM_CPUS_PER_TASK:-1}
pool_id=${SLURM_ARRAY_TASK_ID:-1}

declare -A samples_in_pool=(
  ["1"]="PML16,PML17,PML22"
  ["2"]="PML2,PML8,PML9,PML22"
  ["3"]="PML8,PML9,PML22,PML35"
  ["4"]="PML2,PML8,PML17,PML25"
)
target_pool=MOpool${pool_id}
n_cluster=$(echo ${samples_in_pool[$pool_id]} | sed 's/,/\n/g' | wc -l)

projdir=~/Documents/projects/wp_pml
known_genotypes=${projdir}/outputs/genotypes/batch_2/postimputation/quality_control/all.GRCh38.imputed.snps.biallelic.geno_3e-1.maf_1e-2.annotated_b150.PML_samples.vcf.gz

in_dir=${projdir}/outputs/readcounts/Multiome/${target_pool}

out_dir=${projdir}/outputs/demultiplex/Multiome/${target_pool}
mkdir -p ${out_dir}

# Call genotypes per cell
echo "[I]: Calling genotypes per cell"
singularity exec -B ${projdir} ${projdir}/scripts/tools/dmppl.sif /opt/tools/bin/cellsnp-lite \
  -s ${in_dir}/outs/gex_possorted_bam.bam \
  -b ${in_dir}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz \
  -R ${known_genotypes} \
  -O ${out_dir}/cellsnp_outdir \
  --minMAF 0.05 \
  --minCOUNT 15 \
  -p ${n_cpus}

# Compress cellSNP.base.vcf
singularity exec -B ${projdir} ${projdir}/scripts/tools/dmppl.sif bcftools view \
  -O z \
  -o ${out_dir}/cellsnp_outdir/cellSNP.base.vcf.gz \
  ${out_dir}/cellsnp_outdir/cellSNP.base.vcf

# Demultiplexing supposing donors' genotypes are unknown
echo "[I]: Demultiplexing supposing donors' genotypes are unknown"
singularity exec -B ${projdir} ${projdir}/scripts/tools/dmppl.sif vireo \
  -N ${n_cluster} \
  -c ${out_dir}/cellsnp_outdir \
  -o ${out_dir}/vireo_outdir_denovo

# Demultiplexing using donors' genotypes
echo "[I]: Demultiplexing using donors' genotypes"
mkdir -p $out_dir/known_genotypes
singularity exec -B ${projdir} ${projdir}/scripts/tools/dmppl.sif bcftools index ${out_dir}/cellsnp_outdir/cellSNP.base.vcf.gz
singularity exec -B ${projdir} ${projdir}/scripts/tools/dmppl.sif bcftools view \
  -s ${samples_in_pool[$pool_id]} \
  -R ${out_dir}/cellsnp_outdir/cellSNP.base.vcf.gz \
  -O z \
  -o $out_dir/known_genotypes/known_genotypes.vcf.gz \
  $known_genotypes
singularity exec -B ${projdir} ${projdir}/scripts/tools/dmppl.sif bcftools index -t ${out_dir}/known_genotypes/known_genotypes.vcf.gz

singularity exec -B ${projdir} ${projdir}/scripts/tools/dmppl.sif vireo \
  -t GP \
  -d ${out_dir}/known_genotypes/known_genotypes.vcf.gz \
  -c ${out_dir}/cellsnp_outdir \
  -o ${out_dir}/vireo_outdir_ref
EOF


#
## Demultiplexing by Souporcell
#
# Demultiplexing of scRNA-seq using souporcell pipeline (v2.0)
sbatch --parsable --mem 128G --time 4:59:00 --array 1-4 --cpus-per-task 8 --job-name scRNA-seq.demultiplex \
  --output $projdir/outputs/demultiplex/logs/%A_%a-%u.%x.out <<'EOF'
#!/usr/bin/bash
set -Eeuo pipefail

pool_id=$SLURM_ARRAY_TASK_ID
n_cpus=$SLURM_CPUS_PER_TASK

n_cluster=3
declare -A samples_in_pool=(
  ["1"]="PML0045 PML0061 PML0004"
  ["2"]="PML0055 PML0058 PML0066"
  ["3"]="PML0007 PML0033 PML0070"
  ["4"]="PML0013 PML0060 PML0063"
)

projdir=~/Documents/projects/wp_pml
omics_bam=$projdir/outputs/readcounts/RNA_seq/RNAB2P$pool_id/outs/possorted_genome_bam.bam
omics_barcode=$projdir/outputs/readcounts/RNA_seq/RNAB2P$pool_id/outs/filtered_feature_bc_matrix/barcodes.tsv.gz
known_genotypes=$projdir/outputs/genotypes/batch_2/postimputation/quality_control/all.GRCh38.imputed.snps.biallelic.geno_3e-1.maf_1e-2.annotated_b150.PML_samples.vcf

out_dir=$projdir/outputs/demultiplex/rna_seq/RNAB2P$pool_id

singularity exec -B $projdir $projdir/scripts/tools/souporcell_latest.sif souporcell_pipeline.py \
  --bam $omics_bam \
  --barcodes $omics_barcode \
  --fasta $projdir/inputs/references/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa \
  --threads $n_cpus \
  --no_umi True \
  --clusters $n_cluster \
  --known_genotypes $known_genotypes \
  --known_genotypes_sample_names ${samples_in_pool[$pool_id]} \
  --out_dir $out_dir

# Check concordance between souporcell genotypes and known genotypes.
ref_gt_file=$projdir/outputs/demultiplex/rna_seq/RNAB2P$pool_id/common_variants_covered.vcf
spc_gt_file=$projdir/outputs/demultiplex/rna_seq/RNAB2P$pool_id/cluster_genotypes.vcf

bcftools view $ref_gt_file -Oz -o $ref_gt_file.gz && bcftools index -t $ref_gt_file.gz
bcftools view $spc_gt_file -Oz -o $spc_gt_file.gz && bcftools index -t $spc_gt_file.gz

bcftools gtcheck -u GT -g $ref_gt_file.gz $spc_gt_file.gz 2>/dev/null \
  | awk 'BEGIN{print "cluster\tdonor\tdismat\t-log(Phwe)\tsites"} $1~/^DC/ && $3~/PML/ {print $2"\t"$3"\t"$4"\t"$5"\t"$6}' \
  >| $projdir/outputs/demultiplex/rna_seq/RNAB2P$pool_id/genotype_concordance.txt
EOF
