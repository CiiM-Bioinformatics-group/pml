#!/usr/bin/env bash
# Author: Zhenhua Zhang
# E-mail: zhenhua.zhang3@helmholtz-hzi.de/zhenhua.zhang217@gmail.com
# Created: 2023 Mar 23
# Updated:
# Team: BiiM

#
## Library construction
#
cat <<EOF
CITE-seq library reagents TotalSeq™-B Human Universal Cocktail, V1.0
EOF

# Raw data path
cat <<EOF
#
## PML project raw data.
#

# Batch 1
# Multiome sample path (ATAC): /vol/data/gmak/runs/raw/novaseq/230411_A00278_0544_AH2HHGDRX3
Order-ID  Nr    Sample            Sample_name  Index
23-0066   4_1   Multiome-ATAC_1   MOpool1      SI-NA-D6
23-0066   4_2   Multiome-ATAC_2   MOpool2      SI-NA-D7
23-0066   4_3   Multiome-ATAC_3   MOpool3      SI-NA-D8
23-0066   4_4   Multiome-ATAC_4   MOpool4      SI-NA-D9
 
# Multiome sample path (GEX): /vol/data/gmak/runs/raw/novaseq/230425_A00278_0549_BHJCNYDMXY
Order-ID  Nr    Sample            Sample_name  Index
23-0064   8_1   GEX-TT-G1-Pool1   MOpool1      SI-TT-G1
23-0064   8_2   GEX-TT-G2-Pool2   MOpool2      SI-TT-G2
23-0064   8_3   GEX-TT-G3-Pool3   MOpool3      SI-TT-G3
23-0064   8_4   GEX-TT-G4-Pool4   MOpool4      SI-TT-G4

# CITE-seq sample (CITE) path: /vol/data/gmak/runs/raw/novaseq/230425_A00278_0549_BHJCNYDMXY
Order-ID  Nr    Sample            Sample_name  Index
23-0065   4_1   CITE-NT-A9_Pool1  CCpool1      SI-NT-A9
23-0065   4_2   CITE-NT-A10_Pool2 CCpool2      SI-NT-A10
23-0065   4_3   CITE-NT-A11_Pool3 CCpool3      SI-NT-A11
23-0065   4_4   CITE-NT-A12_Pool4 CCpool4      SI-NT-A12

# CITE-seq sample (GEX) path: /vol/data/gmak/runs/raw/novaseq/230425_A00278_0549_BHJCNYDMXY
Order-ID  Nr    Sample            Sample_name  Index
23-0064   8_5   GEX-TT-F8-Pool1   CGpool1      SI-TT-F8
23-0064   8_6   GEX-TT-F9-Pool2   CGpool2      SI-TT-F9
23-0064   8_7   GEX-TT-F10-Pool3  CGpool3      SI-TT-F10
23-0064   8_8   GEX-TT-F11-Pool4  CGpool4      SI-TT-F11

# Batch 2
# ATAC-seq sample path: /vol/data/gmak/runs/raw/novaseq/230802_A00278_0577_BHGCMVDRX3
Order-ID  Nr    Sample            Sample_name  Index
23-0175   4_1   PMLT1             ATACB2P1      SI-NN-G1
23-0175   4_2   PMLT2             ATACB2P2      SI-NN-G2
23-0175   4_3   PMLT3             ATACB2P3      SI-NN-G3
23-0175   4_4   PMLT4             ATACB2P4      SI-NN-G4

# RNA-seq sample path: 
# RNA-seq sample path: /vol/data/gmak/runs/raw/novaseq/230803_A00278_0578_AH7J7LDRX3
Order-ID  Nr    Sample            Sample_name
23-0174   4_1   PMLR1             RNAB2P1
23-0174   4_2   PMLR2             RNAB2P2
23-0174   4_3   PMLR3             RNAB2P3
23-0174   4_4   PMLR4             RNAB2P4

For more information contact Zhenhua Zhang <Zhenhua.Zhang3@helmholtz-hzi.de> or HZI NGS-Service <NGS-Service@helmholtz-hzi.de>
EOF


projdir=~/Documents/projects/wp_pml

#
## Obtain CellRanger-ARC
#
cellranger_pkg=$projdir/scripts/tools/cellranger-7.1.0.tar.gz
if [[ ! -e $cellranger_pkg ]]; then
  curl -sSL# -o $cellranger_pkg "https://cf.10xgenomics.com/releases/cell-vdj/cellranger-7.1.0.tar.gz?Expires=1683073004&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cHM6Ly9jZi4xMHhnZW5vbWljcy5jb20vcmVsZWFzZXMvY2VsbC12ZGovY2VsbHJhbmdlci03LjEuMC50YXIuZ3oiLCJDb25kaXRpb24iOnsiRGF0ZUxlc3NUaGFuIjp7IkFXUzpFcG9jaFRpbWUiOjE2ODMwNzMwMDR9fX1dfQ__&Signature=Qm-EKpvk3t-HRSaRCO51Va5z2UKdVxkrAtmReTcqNz2faNmypzm5hnoo-SZlAA7BWtlqFhTAVARiOhNXMD7ire~yrFf5s8abq3PZ-qzGeqx5ielwCMdfi2mnI5zqW8l429mrctw-z93~46nmb5czRxX2bjiCI8nBSblJR5Ls5VD44o03DlmIrnx-vISDUg7GiprbohwP8G9j1HeXqYyJI4DIBUEJw5fbkK6QKFdJimRftKyVCZMYXPbstDbFBS02W2VIEQ~FT8GMY3WHQ13OJPwqoe17IZitgSV~ctF-LB6uJch~2vdNqjcHdTAgPpT34kHecjsufZx4F3X1LSb-ng__&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA"
  tar xf $cellranger_pkg -C $(dirname $cellranger_pkg)
else
  echo "cellranger-7.1.0.tar.gz"
fi
export PATH=$projdir/scripts/tools/cellranger-7.1.0/bin:$PATH

cellranger_arc_pkg=$projdir/scripts/tools/cellranger-arc-2.0.2.tar.gz
if [[ ! -e $cellranger_arc_pkg ]]; then
  curl -sSL# -o $cellranger_arc_pkg "https://cf.10xgenomics.com/releases/cell-arc/cellranger-arc-2.0.2.tar.gz?Expires=1682983360&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cHM6Ly9jZi4xMHhnZW5vbWljcy5jb20vcmVsZWFzZXMvY2VsbC1hcmMvY2VsbHJhbmdlci1hcmMtMi4wLjIudGFyLmd6IiwiQ29uZGl0aW9uIjp7IkRhdGVMZXNzVGhhbiI6eyJBV1M6RXBvY2hUaW1lIjoxNjgyOTgzMzYwfX19XX0_&Signature=XNv-N4ZB4tOk9ENJqp5cIXdx98BWlyPh8a1xjYpgL1fAwx8TZiq6BMZI92W7pgAS~thBZ4mnm1OC7UQwtEYRMM~R~SgvR-AgHxzk6mIklNePsBaG4-2y3pfXU7qjJBbUPzp6ARCFeBeveKY3nfmn2CuPp9GsrUqzixTFmu~s1pzAVw4u2XO4-nHXFFNC7MxvF7wGDB1XiXkjN47bnRJ~82uSSVR14m2I9YC4qq9ccyNpLk25LXW3KYtAK-mLfmvz5LzHEhyZ-7-z95xNA047m-adVljOtP3mQ7Rlnh3LcpcH30OR09NnGQmcQcJMWYiEa47VPypbIZlkw9Nf0JukBw__&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA"
  tar xf $cellranger_arc_pkg -C $(dirname $cellranger_arc_pkg)
else
  echo "cellranger-arc-2.0.2.tar.gz is available"
fi
export PATH=$projdir/scripts/tools/cellranger-arc-2.0.2/bin:$PATH

cellranger_atac_pkg=$projdir/scripts/tools/cellranger-atac-2.1.0.tar.gz
if [[ ! -e $cellranger_atac_pkg ]]; then
  curl -sSL# -o $cellranger_atac_pkg "https://cf.10xgenomics.com/releases/cell-atac/cellranger-atac-2.1.0.tar.gz?Expires=1683856095&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cHM6Ly9jZi4xMHhnZW5vbWljcy5jb20vcmVsZWFzZXMvY2VsbC1hdGFjL2NlbGxyYW5nZXItYXRhYy0yLjEuMC50YXIuZ3oiLCJDb25kaXRpb24iOnsiRGF0ZUxlc3NUaGFuIjp7IkFXUzpFcG9jaFRpbWUiOjE2ODM4NTYwOTV9fX1dfQ__&Signature=J8Oz-egbfplvdoMtvxVaj7NDEAJ5V0NV35I103uXOvo7K21q5stcL7EdDgW52QFypv1rZxijORakP3lb-vk0IRgTRlCW2TG5qmyz4ilZhB67WgbzIDsyR9Hxx3dG95DtrqWQTwj8Dy7tYEnvsJLYGTT-0oeO~l9uynSvoKzJukWO~usalwcpd4koREQlFf3z6XDWESI3f4IwhpvAxWH~~cn4~05osePfCpyzhgSIEFIegjvv1j7Q80CEL~vcsrDReyBpodSbb6OIQnwJlzqYERhYRC5f-5k6dAfBAPeggujb7UyQRFLg~CiYpvZt8I6Fr8g4VNfvpMV0I8Ecv2B19w__&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA"
  tar xf $cellranger_atac_pkg -C $(dirname $cellranger_atac_pkg)
else
  echo "$cellranger_atac_pkg is available"
fi
export PATH=$projdir/scripts/tools/cellranger-atac-2.1.0/bin:$PATH


# Software information
# - cellranger, version 7.1.0
# - cellranger-arc, version 2.0.2
# - cellranger-atac, version 2.1.0
# - souporcell, version 2.0
# - samtools, version 1.12
# - bcftools, version 1.12

# Reference information
# - reference genome, GRCh38
# - reference genomic features, GRCh38/xxx


#
## CITE-seq
#
workdir=$projdir/outputs
mkdir -p $workdir/{logs,FASTQs,demultiplex,analysis,readcounts,aggregate}

# Generate FASTQ file from BCL (raw base call) files.
#
## ATAC batch 2
#
cd $workdir/FASTQs
echo -e "Make FASTQs from BCLs of PML ATAC-seq"
cellranger-atac mkfastq \
  --id=ATAC_B2 \
  --jobmode=$projdir/scripts/cellranger_hpc_config/mkfastq/slurm.template \
  --run=$projdir/inputs/sequencing/BCLs/230802_A00278_0577_BHGCMVDRX3 \
  --csv=$projdir/inputs/sequencing/sample_sheet/sample_sheet.batch2.atac.csv

#
## Generating read counts from FASTQ files for ATAC batch 2
#
cd $workdir/readcounts/ATAC_seq
echo -e "Count reads from FASTQs of PML ATAC-seq"
for pp in {1..4}; do
  nohup cellranger-atac count \
    --id=ATACB2P$pp \
    --sample=ATACB2P$pp \
    --jobmode=$projdir/scripts/cellranger_hpc_config/count/slurm.template \
    --reference=$projdir/inputs/references/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
    --fastqs=$projdir/outputs/FASTQs/ATAC_B2/outs 2>&1 1>ATACB2P$pp.count.log &
done

#
## Aggregate read counts
#
# --jobmode=$projdir/scripts/cellranger_hpc_config/count/slurm.template \
cd $workdir/aggregate
echo -e "Aggregate read counts from PML ATAC-seq"
cellranger-atac aggr \
  --id=ATACseq_pool_3_4 \
  --csv=$projdir/inputs/sequencing/sample_sheet/sample_sheet.batch2.atac.aggr.pool_3_4.csv \
  --reference=$projdir/inputs/references/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
  --normalize=depth


#
## RNA-seq, batch 2
#
cd $workdir/FASTQs
echo -e "Make FASTQs from BCLs of PML RNA-seq"
  # --jobmode=$projdir/scripts/cellranger_hpc_config/mkfastq/slurm.template \
cellranger mkfastq \
  --id=RNA_seq \
  --run=$projdir/inputs/sequencing/BCLs/230803_A00278_0578_AH7J7LDRX3 \
  --csv=$projdir/inputs/sequencing/sample_sheet/sample_sheet.batch2.gex.csv

cd $workdir/readcounts/RNA_seq
echo -e "Count reads from FASTQs of RNA-seq"
for pp in {1..4}; do
  cellranger count \
    --id=RNAB2P$pp \
    --sample=RNAB2P$pp \
    --jobmode=$projdir/scripts/cellranger_hpc_config/count/slurm.template \
    --transcriptome=$projdir/inputs/references/refdata-gex-GRCh38-2020-A \
    --fastqs=$projdir/outputs/FASTQs/RNA_seq/outs
done


#
## RNA-seq, including CITE-seq (GEX/CITE) and multiome RNA-seq
#
echo -e "Make FASTQs from BCLs of PML multiome GEX, CITE-seq GEX and CITE"
cellranger mkfastq \
  --id=GEX \
  --jobmode=$projdir/scripts/cellranger_hpc_config/mkfastq/slurm.template \
  --run=$projdir/inputs/sequencing/BCLs/230425_A00278_0549_BHJCNYDMXY \
  --csv=$projdir/inputs/sequencing/sample_sheet/sample_sheet.gex.csv


#
## Generate read counts from FASTQ files for multiome
#
cd $workdir/Multiome
echo -e "Count reads from FASTQs of PML multiome ATAC-seq"
for pp in MOpool{1..4}; do
  cellranger-arc count \
    --id=$pp \
    --jobmode=$projdir/scripts/cellranger_hpc_config/count/slurm.template \
    --reference=$projdir/inputs/references/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
    --libraries=$projdir/inputs/sequencing/sample_sheet/libraries_multiome.$pp.csv
done


#
## Generate read counts from FASTQ files of CITE-seq
#
cd $workdir/CITE_seq
echo -e "Count reads from FASTQs of CITE-seq GEX and CITE"
for pp in pool{1..4}; do
  cellranger count \
    --id=CITE$pp \
    --jobmode=$projdir/scripts/cellranger_hpc_config/count/slurm.template \
    --transcriptome=$projdir/inputs/references/refdata-gex-GRCh38-2020-A \
    --libraries=$projdir/inputs/sequencing/sample_sheet/libraries_citeseq.CITE$pp.csv \
    --feature-ref=$projdir/inputs/references/BioLegend/TotalSeq_B_Human_Universal_Cocktail_V1_399904_Antibody_reference_UMI_counting_CellRanger.csv
done


#
## Demultiplexing, Please check ./demultiplexing.sh
#
# # Demultiplexing of multiome using souporcell pipeline (v2.0)
# sbatch --parsable --mem 128G --time 4:59:00 --array 1-4 --cpus-per-task 8 --job-name multiome_ATACseq.demultiplex \
#   --output $projdir/outputs/demultiplex/logs/%A_%a-%u.%x.out <<'EOF'
# #!/usr/bin/bash
# # set -Eeuo pipefail
# 
# pool_id=$SLURM_ARRAY_TASK_ID
# n_cpus=$SLURM_CPUS_PER_TASK
# 
# if [[ $pool_id == 1 ]]; then n_cluster=3; else n_cluster=4; fi
# declare -A samples_in_pool=(
#   ["1"]="PML16_PML16 PML17_PML17 PML22_PML22"
#   ["2"]="PML2_PML2 PML8_PML8 PML9_PML9 PML22_PML22"
#   ["3"]="PML8_PML8 PML9_PML9 PML22_PML22 PML35_PML35"
#   ["4"]="PML2_PML2 PML8_PML8 PML17_PML17 PML25_PML25"
# )
# 
# projdir=~/Documents/projects/wp_pml
# omics_bam=$projdir/outputs/Multiome/MOpool$pool_id/outs/atac_possorted_bam.bam
# out_dir=$projdir/outputs/demultiplex/multiome/MApool$pool_id
# 
# singularity exec -B $projdir $projdir/scripts/tools/souporcell_latest.sif souporcell_pipeline.py \
#   --bam $omics_bam \
#   --fasta $projdir/inputs/references/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa \
#   --barcodes $projdir/outputs/Multiome/MOpool$pool_id/outs/filtered_feature_bc_matrix/barcodes.tsv.gz \
#   --threads $n_cpus \
#   --no_umi True \
#   --clusters $n_cluster \
#   --known_genotypes $projdir/outputs/genotypes/postimputation/all.GRCh38.imputed.maf5e-2.add_chr.annot.vcf \
#   --known_genotypes_sample_names ${samples_in_pool[$pool_id]} \
#   --out_dir $out_dir
# 
# # Check concordance between souporcell genotypes and known genotypes.
# ref_gt_file=$projdir/outputs/demultiplex/multiome/MApool$pool_id/common_variants_covered.vcf
# spc_gt_file=$projdir/outputs/demultiplex/multiome/MApool$pool_id/cluster_genotypes.vcf
# 
# bcftools view $ref_gt_file -Oz -o $ref_gt_file.gz && bcftools index -t $ref_gt_file.gz
# bcftools view $spc_gt_file -Oz -o $spc_gt_file.gz && bcftools index -t $spc_gt_file.gz
# 
# bcftools gtcheck -u GT -g $ref_gt_file.gz $spc_gt_file.gz 2>/dev/null \
#   | awk 'BEGIN{print "cluster\tdonor\tdismat\t-log(Phwe)\tsites"} $1~/^DC/ && $3~/PML/ {print $2"\t"$3"\t"$4"\t"$5"\t"$6}' \
#   >| $projdir/outputs/demultiplex/multiome/MApool$pool_id/genotype_concordance.txt
# EOF
# 
# 
# # Demultiplexing CITE-seq results
# sbatch --parsable --mem 32G --time 4:59:00 --array 1-4 --cpus-per-task 8 --job-name cite_seq.demultiplex \
#   --output $projdir/outputs/demultiplex/logs/%A_%a-%u.%x.out <<'EOF'
# #!/usr/bin/bash
# # set -Eeuo pipefail
# 
# pool_id=$SLURM_ARRAY_TASK_ID
# n_cpus=$SLURM_CPUS_PER_TASK
# 
# if [[ $pool_id == 1 ]]; then n_cluster=3; else n_cluster=4; fi
# declare -A samples_in_pool=(
#   ["1"]="PML16_PML16 PML17_PML17 PML22_PML22"
#   ["2"]="PML2_PML2 PML8_PML8 PML9_PML9 PML22_PML22"
#   ["3"]="PML8_PML8 PML9_PML9 PML22_PML22 PML35_PML35"
#   ["4"]="PML2_PML2 PML8_PML8 PML17_PML17 PML25_PML25"
# )
# 
# projdir=~/Documents/projects/wp_pml
# 
# # Souporcell pipeline with known genotypes
# singularity exec -B $projdir $projdir/scripts/tools/souporcell_latest.sif souporcell_pipeline.py \
#   --bam $projdir/outputs/CITE_seq/CITEpool$pool_id/outs/possorted_genome_bam.bam \
#   --fasta $projdir/inputs/references/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa \
#   --barcodes $projdir/outputs/CITE_seq/CITEpool$pool_id/outs/filtered_feature_bc_matrix/barcodes.tsv.gz \
#   --threads $n_cpus \
#   --clusters $n_cluster \
#   --known_genotypes $projdir/outputs/genotypes/postimputation/all.GRCh38.imputed.maf5e-2.add_chr.annot.vcf \
#   --known_genotypes_sample_names ${samples_in_pool[$pool_id]} \
#   --out_dir $projdir/outputs/demultiplex/cite_seq/CITEpool$pool_id
# 
# # Check concordance between souporcell genotypes and known genotypes.
# ref_gt_file=$projdir/outputs/demultiplex/cite_seq/CITEpool$pool_id/common_variants_covered.vcf
# spc_gt_file=$projdir/outputs/demultiplex/cite_seq/CITEpool$pool_id/cluster_genotypes.vcf
# 
# bcftools view $ref_gt_file -Oz -o $ref_gt_file.gz && bcftools index -t $ref_gt_file.gz
# bcftools view $spc_gt_file -Oz -o $spc_gt_file.gz && bcftools index -t $spc_gt_file.gz
# 
# bcftools gtcheck -u GT -g $ref_gt_file.gz $spc_gt_file.gz 2>/dev/null \
#   | awk 'BEGIN{print "cluster\tdonor\tdismat\t-log(Phwe)\tsites"} $1~/^DC/ && $3~/PML/ {print $2"\t"$3"\t"$4"\t"$5"\t"$6}' \
#   >| $projdir/outputs/demultiplex/cite_seq/CITEpool$pool_id/genotype_concordance.txt
# EOF
