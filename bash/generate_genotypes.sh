#!/usr/bin/env bash
# generate_genotypes.sh
# Author: Zhenhua Zhang
# E-mail: zhenhua.zhang217@gmail.com
# Created: 2023 May 17
# Updated: 2023 May 17

# Tools
# - PLINK v1.90b6.26 64-bit (2 Apr 2022)
# - BCFtools 1.12
# - HTSlib 1.12
#

# Problem samples (Sep 15, 2023)
# | Sample-ID |                                                        Issue | Removed |
# | :-------: | :----------------------------------------------------------: | :-----: |
# |      C168 |                          failed IBD check (PHI_HAT > 0.1825) | NO      |
# |     IBD23 |                             potential sex swap (F = 0.02764) | NO      |
# |     IBD29 |                           genotype call rate (0.3846) < 0.90 | YES     |
# |     IBD76 |                          failed IBD check duplicated to C076 | YES     |
# |     IBD28 |                         failed IBD check; duplicated to C028 | YES     |
# |     IBD49 | uncertain sex (F=0.4591); genotype call rate (0.5993) < 0.90 | YES     |


projdir=~/Documents/projects/wp_pml
workdir=$projdir/outputs/genotypes/batch_2
bed_files=$workdir/preimputation/PML_IBD_LongCovid_merge
wk_bed_files=$workdir/preimputation/PML_IBD_LongCovid_merge.wk

#
## Merge multiple ped/map, bed/bim/fam files, including all PML samples, part of IBD samples, part of LongCOVID samples
#
genotype_path=/vol/projects/BIIM/PML/inputs/sequencing/genotypes
cat <<EOF > $workdir/batch_2/genotype_merge_list.txt
$genotype_path/GSA2023_1088_helmholtz_TAT_V3/PLINK_130923_0401/GSA2023_1088_helmholtz_TAT_V3.ped $genotype_path/GSA2023_1088_helmholtz_TAT_V3/PLINK_130923_0401/GSA2023_1088_helmholtz_TAT_V3.map
$genotype_path/GSa2023_1050_HelmHoltzTAT_V3/PLINK_150523_0112/GSa2023_1050_HelmHoltzTAT_V3.ped $genotype_path/GSa2023_1050_HelmHoltzTAT_V3/PLINK_150523_0112/GSa2023_1050_HelmHoltzTAT_V3.map
$genotype_path/GSA2022_953_025_V3_TAT/PLINK_030822_1240/GSA2022_953_025_V3_TAT.ped $genotype_path/GSA2022_953_025_V3_TAT/PLINK_030822_1240/GSA2022_953_025_V3_TAT.map
EOF
plink --merge-list $workdir/genotype_merge_list.txt --make-bed --out $workdir/preimputation/PML_IBD_LongCovid_merge


#
## Pre-imputation QC
#

# Update sample ids
cat <<EOF >| $workdir/new_ids.txt
1 4100786679 4100786679 4100786679
1 C013 C013 C013
1 PML17 PML17 PML17
2 C169 C169 C169
2 IBD89 IBD89 IBD89
2 PML2 PML2 PML2
3 2909366720 2909366720 2909366720
3 C110 C110 C110
3 PML16 PML16 PML16
4 IBD29 IBD29 IBD29
4 PML0004 PML0004 PML0004
4 PML22 PML22 PML22
5 4100934317 4100934317 4100934317
5 C054 C054 C054
5 PML9 PML9 PML9
6 703019410A 703019410A 703019410A
6 C205 C205 C205
6 IBD68 IBD68 IBD68
7 2302725222 2302725222 2302725222
7 C156 C156 C156
7 PML35 PML35 PML35
8 1207531220A 1207531220A 1207531220A
8 IBD49 IBD49 IBD49
8 PML0058 PML0058 PML0058
9 4100938336 4100938336 4100938336
9 C014 C014 C014
9 PML8 PML8 PML8
10 906612711A 906612711A 906612711A
10 C171 C171 C171
10 IBD57 IBD57 IBD57
11 C127 C127 C127
11 IBD12 IBD12 IBD12
11 PML25 PML25 PML25
12 IBD48 IBD48 IBD48
12 PML0007 PML0007 PML0007
13 4100938478 4100938478 4100938478
13 C061 C061 C061
14 C207 C207 C207
14 IBD76 IBD76 IBD76
15 C157 C157 C157
15 IBD65 IBD65 IBD65
16 IBD51 IBD51 IBD51
16 PML0060 PML0060 PML0060
17 4100951226 4100951226 4100951226
17 C028 C028 C028
18 C174 C174 C174
18 IBD47 IBD47 IBD47
19 C129 C129 C129
19 IBD32 IBD32 IBD32
20 0509757320 0509757320 0509757320
20 PML0013 PML0013 PML0013
21 C063 C063 C063
21 IBD23 IBD23 IBD23
22 C216 C216 C216
22 IBD28 IBD28 IBD28
23 C159 C159 C159
24 PML0061 PML0061 PML0061
25 C031 C031 C031
26 C178 C178 C178
27 C145 C145 C145
28 PML0033 PML0033 PML0033
29 C076 C076 C076
30 C218 C218 C218
31 C161 C161 C161
32 PML0063 PML0063 PML0063
33 C037 C037 C037
34 C184 C184 C184
35 C149 C149 C149
36 PML0045 PML0045 PML0045
37 C101 C101 C101
38 C219 C219 C219
39 C166 C166 C166
40 PML0066 PML0066 PML0066
41 C043 C043 C043
42 C192 C192 C192
43 C154 C154 C154
44 PML0055 PML0055 PML0055
45 C105 C105 C105
46 C220 C220 C220
47 C168 C168 C168
48 PML0070 PML0070 PML0070
EOF
plink --bfile $bed_files --update-ids $workdir/new_ids.txt --make-bed --out $bed_files.tmp.new_ids

# Update sample sex
cat <<EOF >| $workdir/new_sex.txt
4100934317 4100934317 M
C166 C166 F
C168 C168 M
C207 C207 M
IBD23 IBD23 F
IBD49 IBD49 M
0509757320 0509757320 F
1207531220A 1207531220A F
2302725222 2302725222 F
2909366720 2909366720 F
4100786679 4100786679 F
4100938336 4100938336 M
4100938478 4100938478 M
4100951226 4100951226 M
703019410A 703019410A M
906612711A 906612711A M
C013 C013 F
C014 C014 F
C028 C028 M
C031 C031 F
C037 C037 M
C043 C043 F
C054 C054 M
C061 C061 M
C063 C063 F
C076 C076 F
C101 C101 F
C105 C105 M
C110 C110 F
C127 C127 M
C129 C129 M
C145 C145 F
C149 C149 F
C154 C154 F
C156 C156 F
C157 C157 F
C159 C159 F
C161 C161 M
C169 C169 M
C171 C171 M
C174 C174 F
C178 C178 M
C184 C184 F
C192 C192 F
C205 C205 M
C216 C216 F
C218 C218 F
C219 C219 F
C220 C220 F
PML0004 PML0004 M
PML0007 PML0007 M
PML0013 PML0013 F
PML0033 PML0033 M
PML0045 PML0045 M
PML0055 PML0055 M
PML0058 PML0058 M
PML0060 PML0060 M
PML0061 PML0061 M
PML0063 PML0063 M
PML0066 PML0066 M
PML0070 PML0070 M
IBD12 IBD12 M
IBD28 IBD28 M
IBD29 IBD29 F
IBD32 IBD32 M
IBD47 IBD47 F
IBD48 IBD48 F
IBD51 IBD51 M
IBD57 IBD57 F
IBD65 IBD65 F
IBD68 IBD68 M
IBD76 IBD76 F
IBD89 IBD89 F
PML16 PML16 M
PML17 PML17 M
PML2 PML2 M
PML22 PML22 F
PML25 PML25 F
PML35 PML35 F
PML8 PML8 M
PML9 PML9 M
EOF
plink --bfile $bed_files.tmp.new_ids --update-sex $workdir/new_sex.txt --make-bed --out $wk_bed_files

# Quality control at sample level
# IBD23, potential sex swap, correct gender to PLINK's estimation
plink --bfile $wk_bed_files --check-sex 0.2 0.75 --out $workdir/preimputation/sample_level.sex_status
grep -i "problem" $workdir/preimputation/sample_level.sex_status.sexcheck >| $workdir/preimputation/fail-sexcheck.txt

# Check missning data rates, genotype failure rate >= 0.05, NO sample was removed
plink --bfile $wk_bed_files --missing --out $workdir/preimputation/sample_level.genotype_missing_rate
awk '$1 !~ /FID/ && $6 > 0.1 {print}' $workdir/preimputation/sample_level.genotype_missing_rate.imiss >| $workdir/preimputation/fail-genotype_missing_rate.txt 

# Check heterozygosity rate, heterozygosity rate ± 3 stdev, NO sample was removed
# The code was from https://github.com/genepi-freiburg/gwas/blob/master/cleaning-pipeline/cleaning/03-heterozygosity.R
plink --bfile $wk_bed_files --het --out $workdir/preimputation/sample_level.heterozygosity_rate
Rscript $projdir/scripts/r/analyze_heterozygosity.r $workdir/preimputation/sample_level.heterozygosity_rate $workdir/preimputation/sample_level.heterozygosity_rate.estimate.txt


# Check relatness across samples. NO samples were removed due to high relatness (IBD > 0.1875)
# High-LD regions from https://dougspeed.com/high-ld-regions/
cat <<EOF >| $workdir/high-LD-regions.txt
1 48000000 52000000 1
2 86000000 100500000 2
2 183000000 190000000 3
3 47500000 50000000 4
3 83500000 87000000 5
5 44500000 50500000 6
5 129000000 132000000 7
6 25500000 33500000 8
6 57000000 64000000 9
6 140000000 142500000 10
7 55000000 66000000 11
8 8000000 12000000 12
8 43000000 50000000 13
8 112000000 115000000 14
10 37000000 43000000 15
11 87500000 90500000 16
12 33000000 40000000 17
20 32000000 34500000 18
EOF
plink --bfile $wk_bed_files --exclude $workdir/high-LD-regions.txt --range --indep-pairwise 50 5 0.2 --out $workdir/preimputation/sample_level.relatness
plink --bfile $wk_bed_files --extract $workdir/preimputation/sample_level.relatness.prune.in --genome --out $workdir/preimputation/sample_level.ibd

perl $projdir/scripts/perls/run-IBD-QC.pl $workdir/preimputation/sample_level.genotype_missing_rate $workdir/preimputation/sample_level.ibd
Rscript $projdir/scripts/r/plot_IBD.r $workdir/preimputation/sample_level.ibd $workdir/preimputation $workdir/preimputation/PML_IBD_LongCovid_merge.wk.fam


# Quality control at SNP level
# MAF
plink --bfile $wk_bed_files --freq --out $workdir/preimputation/snp_level.freq

# Probe (SNP) level missingness (duplicated to sample level genotype missingness check)
plink --bfile $wk_bed_files --missing --out $workdir/preimputation/snp_level.genotype_missing_rate

# High differential missingness between case and control, skipped due to the phenotypes are not binary.
# plink --bfile $wk_bed_files --test-missing --out $workdir/preimputation/snp_level.diffenertial_missingness

# HWE test.
plink --bfile $wk_bed_files --hardy --out $workdir/preimputation/snp_level.hardy_weinberg_equilibrium


# The final clean genotypes ready for imputation or association analysis.
# 77 individuals, 485,104 variants, total genotyping rate is 0.996079
cat <<EOF >| $workdir/preimputation/fail-sample_id.txt
IBD28 IBD28
IBD29 IBD29
IBD49 IBD49
IBD76 IBD76
EOF
plink --bfile $wk_bed_files --remove $workdir/preimputation/fail-sample_id.txt --maf 0.01 --geno 0.05 --hwe 1e-5 --not-chr 0 --snps-only just-acgt --make-bed --out $workdir/preimputation/PML_IBD_LongCovid_merge.maf_1e-2.hwe_1e-5.geno_5e-2

# Convert it into VCF format
plink --bfile $workdir/preimputation/PML_IBD_LongCovid_merge.maf_1e-2.hwe_1e-5.geno_5e-2 --recode vcf-iid bgz --out $workdir/preimputation/PML_IBD_LongCovid_merge.maf_1e-2.hwe_1e-5.geno_5e-2
bcftools index -ft $workdir/preimputation/PML_IBD_LongCovid_merge.maf_1e-2.hwe_1e-5.geno_5e-2.vcf.gz
bcftools annotate --rename-chrs <(for x in {1..26}; do echo -e $x"\tchr"$x; done) -Oz -o $workdir/preimputation/PML_IBD_LongCovid_merge.maf_1e-2.hwe_1e-5.geno_5e-2.with_chr.vcf.gz $workdir/preimputation/PML_IBD_LongCovid_merge.maf_1e-2.hwe_1e-5.geno_5e-2.vcf.gz

bcftools index $workdir/preimputation/PML_IBD_LongCovid_merge.maf_1e-2.hwe_1e-5.geno_5e-2.with_chr.vcf.gz
seq 1 26 | xargs -n1 -I% -P 26 bcftools view -Oz -o $workdir/preimputation/per_chr/chr%.vcf.gz $workdir/preimputation/PML_IBD_LongCovid_merge.maf_1e-2.hwe_1e-5.geno_5e-2.vcf.gz %

# Clean up
rm -f $workdir/preimputation/*.{tmp,wk}.* && mkdir -p QC && mv -f *.pdf *.txt snp_level.* sample_level.* QC


#
## Imputation. Done by Michigen Imputation Server
#
# Parameters:
# - Reference panel: TOPMed r2
# - Array Build: GRCh37/hg19
# - rsq Filter: off
# - Phasing: Egle v2.4
# - Population: vs. TOPMed Panel
#
mkdir -p chr1_3
cd chr1_3
curl -sL https://imputation.biodatacatalyst.nhlbi.nih.gov/get/1049075/78db15c7b1357ea779fe821964cd5606ea9efc99855fdafc84e103a682edc09f | bash
curl -sL https://imputation.biodatacatalyst.nhlbi.nih.gov/get/1049079/9e4b23138f4ac266403a1341f79baac787b6134afaf87bc0b9ab923a92b48638 | bash
curl -sL https://imputation.biodatacatalyst.nhlbi.nih.gov/get/1049081/f1ff0c58f0e8e0fac4c242f83a927f1ffb2d6a1ebe7656eea1ea954100b985c2 | bash
curl -sL https://imputation.biodatacatalyst.nhlbi.nih.gov/get/1049082/b067528bf6cf0ecb1ff5c7aa4974f5dcdf6a11926222b2b41305127d65d6181f | bash
cd -

mkdir -p chr4_6
cd chr4_6
curl -sL https://imputation.biodatacatalyst.nhlbi.nih.gov/get/1049141/936b645a1459c641f1fcba11934e472d14a5b186f7fce68e385016d6b8c22da9 | bash
curl -sL https://imputation.biodatacatalyst.nhlbi.nih.gov/get/1049145/bfc86d6de7da2fce9cca56c9fede25c6c55d9ee87d48069320d05ad12ba036d6 | bash
curl -sL https://imputation.biodatacatalyst.nhlbi.nih.gov/get/1049147/aad5ba20052437dcca934c391016940a3ce6316b40858c49c917500115d8b8b1 | bash
curl -sL https://imputation.biodatacatalyst.nhlbi.nih.gov/get/1049148/5289a2417b897c485c31b7d48290a595c1966133af823f1a0a767fc3ec326fe9 | bash
cd -

mkdir -p chr7_11
cd chr7_11
curl -sL https://imputation.biodatacatalyst.nhlbi.nih.gov/get/1049163/b3794a1a8af17584c57820f3eabc50586fbe038a0a6fd9b5b169ddad68cec057 | bash
curl -sL https://imputation.biodatacatalyst.nhlbi.nih.gov/get/1049167/0797cadd2e963e0a2eb5916422f7c1c6b0728b5739bae247951defb68e0886e7 | bash
curl -sL https://imputation.biodatacatalyst.nhlbi.nih.gov/get/1049169/766102c597c9f17ce1380f889c08339c7d109039116ce53c4d09410c71c1d258 | bash
curl -sL https://imputation.biodatacatalyst.nhlbi.nih.gov/get/1049170/d787bea8c61255a6c6f2fa5889f0676c45e1930d76d4ee86b690052882fcc989 | bash
cd -

mkdir -p chr12_18
cd chr12_18
curl -sL https://imputation.biodatacatalyst.nhlbi.nih.gov/get/1049207/c4ba635cb2324a73d80438bf3a4edcfd3be94b6ad8bbc9f98ebf8f0493daa504 | bash
curl -sL https://imputation.biodatacatalyst.nhlbi.nih.gov/get/1049211/0be622aa8a747f46f26d441e1a23a0b12075e4fc5b8303e9d9596f732062987d | bash
curl -sL https://imputation.biodatacatalyst.nhlbi.nih.gov/get/1049213/684363b8c8a5bc06e18163bb46f095a09001fa95e46244125d8c995cb3c01986 | bash
curl -sL https://imputation.biodatacatalyst.nhlbi.nih.gov/get/1049214/80ddf3051c5697c21637ea6f010e7dd8607eda5e4b0cb555354fbdb102c99a8b | bash
cd -

mkdir -p chr19_23
cd chr19_23
curl -sL https://imputation.biodatacatalyst.nhlbi.nih.gov/get/1049295/41b4c1299b3e852474245b893425244ae9eee749844c9b6d1ad2b64507761d4e | bash
curl -sL https://imputation.biodatacatalyst.nhlbi.nih.gov/get/1049299/fd96a0568af0065d6621a209c235ee0f275e3328861f556e3204235971115a63 | bash
curl -sL https://imputation.biodatacatalyst.nhlbi.nih.gov/get/1049301/751c547358d3fc05a875045f022822a621dc304a65d657e195ed0bd1b4d3af0e | bash
curl -sL https://imputation.biodatacatalyst.nhlbi.nih.gov/get/1049302/90b5b61cb74d42dd1467fdcaf37ff8c06b7ef7f23eeacdd4f993931388825237 | bash
cd -


while read -r chr passwd; do
  unzip -P "$passwd" zipped/$chr.zip -d unzipped
done < zipped/password.txt



#
## Post-imputation QC
#
for x in $projdir/outputs/genotypes/batch_2/postimputation/unzipped/chr*.dose.vcf.gz; do
  echo Indexing $x ...
  bcftools index -f -t --threads 8 $x
done


# Concat all chroms
bcftools concat \
  -a \
  -O z \
  -o $projdir/outputs/genotypes/batch_2/postimputation/quality_control/all.GRCh38.imputed.raw.vcf.gz \
  --threads 8 \
  --no-version \
  $projdir/outputs/genotypes/batch_2/postimputation/unzipped/{chr{1..22}.dose.vcf.gz,chrX.dose.vcf.gz}
bcftools index -f -t --threads 8 $projdir/outputs/genotypes/batch_2/postimputation/quality_control/all.GRCh38.imputed.raw.vcf.gz


# Filter the variants
bcftools view \
  -v snps \
  -m 2 -M 2 \
  -q 0.01:minor \
  -i 'R2 > 0.3 || ER2 > 0.3' \
  -O z \
  -o $projdir/outputs/genotypes/batch_2/postimputation/quality_control/all.GRCh38.imputed.snps.biallelic.geno_3e-1.maf_1e-2.vcf.gz \
  --no-version \
  $projdir/outputs/genotypes/batch_2/postimputation/quality_control/all.GRCh38.imputed.raw.vcf.gz
bcftools index -f -t --threads 4 $projdir/outputs/genotypes/batch_2/postimputation/quality_control/all.GRCh38.imputed.snps.biallelic.geno_3e-1.maf_1e-2.vcf.gz


# Annotate the variants using rsID from dbSNP GRCh38 b150
bcftools annotate \
  -c ID \
  -a /vol/projects/BIIM/resources/snpDB/GRCh38_b150/00-All.with_chr.vcf.gz \
  -o $projdir/outputs/genotypes/batch_2/postimputation/quality_control/all.GRCh38.imputed.snps.biallelic.geno_3e-1.maf_1e-2.annotated_b150.vcf.gz \
  --threads 4 \
  --no-version \
  $projdir/outputs/genotypes/batch_2/postimputation/quality_control/all.GRCh38.imputed.snps.biallelic.geno_3e-1.maf_1e-2.vcf.gz
bcftools index -f -t --threads 4 $projdir/outputs/genotypes/batch_2/postimputation/quality_control/all.GRCh38.imputed.snps.biallelic.geno_3e-1.maf_1e-2.annotated_b150.vcf.gz

# Obtain PML samples
bcftools view \
  --no-version \
  -O z \
  -o $projdir/outputs/genotypes/batch_2/postimputation/quality_control/all.GRCh38.imputed.snps.biallelic.geno_3e-1.maf_1e-2.annotated_b150.PML_samples.vcf.gz \
  -S <(zgrep -m1 CHROM $projdir/outputs/genotypes/batch_2/postimputation/quality_control/all.GRCh38.imputed.snps.biallelic.geno_3e-1.maf_1e-2.annotated_b150.vcf.gz | tr "\t" "\n" | grep "PML") \
  $projdir/outputs/genotypes/batch_2/postimputation/quality_control/all.GRCh38.imputed.snps.biallelic.geno_3e-1.maf_1e-2.annotated_b150.vcf.gz

# Cleanup
rm -fr $projdir/outputs/genotypes/batch_2/postimputation/unzipped/*
# rm -fr $projdir/outputs/genotypes/postimputation/quality_control/*
# rm -fr $projdir/outputs/genotypes/postimputation/TOPMed/unzipped/*
