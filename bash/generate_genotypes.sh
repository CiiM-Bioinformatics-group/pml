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

projdir=
workdir=$projdir/outputs/genotypes
ped_files=$projdir/inputs/sequencing/genotypes/GSa2023_1050_HelmHoltzTAT_V3/PLINK_150523_0112/GSa2023_1050_HelmHoltzTAT_V3
bed_files=$projdir/outputs/genotypes/preimputation/GSa2023_1050_HelmHoltzTAT_V3
wk_bed_file=$projdir/outputs/genotypes/preimputation/GSa2023_1050_HelmHoltzTAT_V3.tmp.new_ids.new_sex


#
## Pre-imputation QC
#

# Convert the ped/map into bed/bim/fam
plink --file $ped_files --make-bed --out $bed_files

# Update sample ids
cat <<EOF >| $workdir/new_ids.txt
1 PML17 PML17 PML17
2 PML2 PML2 PML2
3 PML16 PML16 PML16
4 PML22 PML22 PML22
5 PML9 PML9 PML9
6 703019410A 703019410A 703019410A
7 PML35 PML35 PML35
8 1207531220A 1207531220A 1207531220A
9 PML8 PML8 PML8
10 906612711A 906612711A 906612711A
11 PML25 PML25 PML25
EOF
plink --bfile $bed_files --update-ids $workdir/new_ids.txt --make-bed --out $bed_files.tmp.new_ids

# Update sample sex
cat <<EOF >| $workdir/new_sex.txt
PML17 PML17 M
PML2 PML2 M
PML16 PML16 M
PML22 PML22 F
PML9 PML9 M
703019410A 703019410A M
PML35 PML35 F
1207531220A 1207531220A F
PML8 PML8 M
906612711A 906612711A M
PML25 PML25 F
EOF
plink --bfile $bed_files.tmp.new_ids --update-sex $workdir/new_sex.txt --make-bed --out $wk_bed_file

# Quality control at sample level
# Check sex, NO sample was removed
plink --bfile $wk_bed_file --check-sex --out $workdir/preimputation/sample_level.sex_status
grep -i "problem" $workdir/preimputation/sample_level.sex_status.sexcheck >| $workdir/preimputation/fail-sexcheck.txt

# Check missning data rates, genotype failure rate >= 0.05, NO sample was removed
plink --bfile $wk_bed_file --missing --out $workdir/preimputation/sample_level.genotype_missing_rate
awk '$1 !~ /FID/ && $6 > 0.05 {print}' $workdir/preimputation/sample_level.genotype_missing_rate.imiss >| $workdir/preimputation/fail-genotype_missing_rate.txt 

# Check heterozygosity rate, heterozygosity rate ± 3 stdev, NO sample was removed
# The code was from https://github.com/genepi-freiburg/gwas/blob/master/cleaning-pipeline/cleaning/03-heterozygosity.R
plink --bfile $wk_bed_file --het --out $workdir/preimputation/sample_level.heterozygosity_rate
Rscript $projdir/scripts/r/analyze_heterozygosity.R $workdir/preimputation/sample_level.heterozygosity_rate $workdir/preimputation/sample_level.heterozygosity_rate.estimate.txt


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
plink --bfile $wk_bed_file --exclude $workdir/high-LD-regions.txt --range --indep-pairwise 50 50 0.2 --out $workdir/preimputation/sample_level.relatness
plink --bfile $wk_bed_file --extract $workdir/preimputation/sample_level.relatness.prune.in --genome --out $workdir/preimputation/sample_level.ibd

perl $projdir/scripts/perls/run-IBD-QC.pl $workdir/preimputation/sample_level.genotype_missing_rate $workdir/preimputation/sample_level.ibd
Rscript $projdir/scripts/r/plot-IBD.R $workdir/preimputation/sample_level.ibd $workdir/preimputation $workdir/preimputation/GSa2023_1050_HelmHoltzTAT_V3.tmp.new_ids.new_sex.fam


# Quality control at SNP level
# MAF
plink --bfile $wk_bed_file --freq --out $workdir/preimputation/snp_level.freq

# Probe (SNP) level missingness (duplicated to sample level genotype missingness check)
plink --bfile $wk_bed_file --missing --out $workdir/preimputation/snp_level.genotype_missing_rate

# High differential missingness between case and control, skipped due to the phenotypes are not binary.
# plink --bfile $wk_bed_file --test-missing --out $workdir/preimputation/snp_level.diffenertial_missingness

# HWE test.
plink --bfile $wk_bed_file --hardy --out $workdir/preimputation/snp_level.hardy_weinberg_equilibrium

# The final clean genotypes ready for imputation or association analysis
plink --bfile $wk_bed_file --maf 0.01 --geno 0.05 --hwe 1e-5 --make-bed --out $workdir/preimputation/GSa2023_1050_HelmHoltzTAT_V3.clean

# Convert it into VCF format
plink --bfile $workdir/preimputation/GSa2023_1050_HelmHoltzTAT_V3.clean --recode vcf --out $workdir/preimputation/GSa2023_1050_HelmHoltzTAT_V3.clean

# Clean up
rm -f $workdir/*.tmp.* && mkdir -p QC && mv *.pdf *.txt snp_level.* sample_level.* QC


#
## Imputation. Done by Michigen Imputation Severe
#
# Parameters:
# - Reference panel: TOPMed r2
# - Array Build: GRCh37/hg19
# - rsq Filter: off
# - Phasing: Egle v2.4
# - Population: vs. TOPMed Panel

# ChrX
mkdir -p chrx
cd chrx
curl -sL https://imputation.biodatacatalyst.nhlbi.nih.gov/get/931639/9fda4974c98de92e4d435648b2c92e926d6d893c814b1f9b6527889272e92f5f | bash 
curl -sL https://imputation.biodatacatalyst.nhlbi.nih.gov/get/931643/fc9c32bdcb86418712e0bdac418a7e7d8290f645a7e7ebf0626885e7b9a4e4b1 | bash
curl -sL https://imputation.biodatacatalyst.nhlbi.nih.gov/get/931645/e88f1b689b449c90166f062793c0fbbb245803f5b174c40de35711ef4f57af3c | bash
curl -sL https://imputation.biodatacatalyst.nhlbi.nih.gov/get/931646/709bf0b65a8653fe5efe1a38ed8f5522f47b665fc13bef3be4328440c849f48d | bash
cd -

# Chr 17-22
mkdir -p chr17to22
cd chr17to22
curl -sL https://imputation.biodatacatalyst.nhlbi.nih.gov/get/931331/4bd6e2bd9d16874929bb8f9eb4dea6f8ae6657fddb19a943605718967b9e7efe | bash
curl -sL https://imputation.biodatacatalyst.nhlbi.nih.gov/get/931335/f2ebec13550d83c98dd7d59463a7e1e6a12a59205038024867f2e040780fe692 | bash
curl -sL https://imputation.biodatacatalyst.nhlbi.nih.gov/get/931337/002d460f9ca5c8adcf43d3b0364ef2a90bf1d67a2936ec49b26e0755bc08f6a6 | bash
curl -sL https://imputation.biodatacatalyst.nhlbi.nih.gov/get/931338/fa4c0f633831f9fc17b581d0f3244df0814a8bb4474c3817adc8cbcfb13ae630 | bash
cd -

# Chr 11-16
mkdir -p chr11to16
cd chr11to16
curl -sL https://imputation.biodatacatalyst.nhlbi.nih.gov/get/931243/1d5eb2ed0c4e6ac9f3347225ed0b94dcf3b4089466dad11f88ed474ecf76dd96 | bash
curl -sL https://imputation.biodatacatalyst.nhlbi.nih.gov/get/931247/6468535dd77ddc10c4eb74db66809673fcf922f87f9a81ee84f51f2ad2df2825 | bash
curl -sL https://imputation.biodatacatalyst.nhlbi.nih.gov/get/931249/47b573f7a439fdbd2f9c3bfd8d0f5e39353111abd8070a248f2aae96e68ed342 | bash
curl -sL https://imputation.biodatacatalyst.nhlbi.nih.gov/get/931250/2ea27887de2486b7c79d91c2626464c52cca6021cf576d8d10c75c78d8ee2f09 | bash
cd -

# Chr 8-10
mkdir -p chr8to10
cd chr8to10
curl -sL https://imputation.biodatacatalyst.nhlbi.nih.gov/get/930253/6ff252aef5f355a545b6fe757994522d2763f844311bdbf389b34e0c2330ca75 | bash
curl -sL https://imputation.biodatacatalyst.nhlbi.nih.gov/get/930257/a40b703dbaf57f961caf7a7c8d415aac1f46e2eaf8863967070d97d3eed011a1 | bash
curl -sL https://imputation.biodatacatalyst.nhlbi.nih.gov/get/930259/35784b123ca5ce3ea6091a5ac5d4a1e0c2a341dc47b18c2ea69fba0f43c9f14a | bash
curl -sL https://imputation.biodatacatalyst.nhlbi.nih.gov/get/930260/ec66b5fdaf9e0403faa48422f985c0a36918c1128332a6e0148fc10919028d3e | bash
cd -

# Chr 6,7
mkdir -p chr6to7
cd chr6to7
curl -sL https://imputation.biodatacatalyst.nhlbi.nih.gov/get/930231/dca00da3a336e042996bb1e84cd1c5f199c323d797efc1807bf7c7b1c9e39d10 | bash
curl -sL https://imputation.biodatacatalyst.nhlbi.nih.gov/get/930235/f4bd68f60bffe72798028b2dc22c2caf10203ca1fbee96dbe2daef16d9e11737 | bash
curl -sL https://imputation.biodatacatalyst.nhlbi.nih.gov/get/930237/bb685e7873d8a12a9fa5d4725501aa87bfad9257681aa49dec71bc85d5cc36ec | bash
curl -sL https://imputation.biodatacatalyst.nhlbi.nih.gov/get/930238/d6d90c6f533480c5e16566787f11e9f357d8c682d05d1dd2fe9315468bd3b3ab | bash
cd -

# Chr 4,5
mkdir -p chr4to5
cd chr4to5
curl -sL https://imputation.biodatacatalyst.nhlbi.nih.gov/get/930209/83bea84567da83898b8644706ec035ee0441bb277ae86f724f0ca86a70081bb3 | bash
curl -sL https://imputation.biodatacatalyst.nhlbi.nih.gov/get/930213/8c981ca7e1ff14c8fca65d62ddc3db65053289d99a11cd81db6f81ca2d76a054 | bash
curl -sL https://imputation.biodatacatalyst.nhlbi.nih.gov/get/930215/fdad000e7b667946ead306dd67140a9dfba121a9f62cf9289235f7924ee929a2 | bash
curl -sL https://imputation.biodatacatalyst.nhlbi.nih.gov/get/930216/d961f7ea349e0574ee1d0828a44ba0673669b24bb65dea313ce2ea693f2561d4 | bash
cd -

# Chr 3
mkdir -p chr3
cd chr3
curl -sL https://imputation.biodatacatalyst.nhlbi.nih.gov/get/928845/585545c945bb8780a7d5681f635788c9da24b0aa24d0ebfe32d5d8e6cf24a777 | bash
curl -sL https://imputation.biodatacatalyst.nhlbi.nih.gov/get/928849/8e6a6d5551aa0e30ae3005fa383664005a6be7e051f4673aab5e7ecb2ac4b28a | bash
curl -sL https://imputation.biodatacatalyst.nhlbi.nih.gov/get/928851/a18f53407c89054bb538a642e9f9baafd213820964e2dd63555da96e573e9873 | bash
curl -sL https://imputation.biodatacatalyst.nhlbi.nih.gov/get/928852/b94c52675619342952e6a525c5f5d75d9dd2690310ce12b178e1d5f6ff13ae00 | bash
cd -

# Chr 2
mkdir -p chr2
cd chr2
curl -sL https://imputation.biodatacatalyst.nhlbi.nih.gov/get/928823/bddcbb8107cf7ac902b626c4223b3d4791aa05775d054b1c2c92f1d7d22ca451 | bash
curl -sL https://imputation.biodatacatalyst.nhlbi.nih.gov/get/928827/f679c3439f7ace0bdc7cdda879141ef5294366b0e98cd528ebe4fa195992a1df | bash
curl -sL https://imputation.biodatacatalyst.nhlbi.nih.gov/get/928829/c55a0203770bdcd4047f5df361282c6f51a5bd190ada430dcf4f76aaf3ffd54a | bash
curl -sL https://imputation.biodatacatalyst.nhlbi.nih.gov/get/928830/bdb755452016b204ef2a8018b8d05c999df0f12833d3d3fc759a07615114bb73 | bash
cd -

# Chr 1
mkdir -p chr1
cd chr1
curl -sL https://imputation.biodatacatalyst.nhlbi.nih.gov/get/928757/1100dada37b9c8ce20619524791ae57ac621b7e21791cb9b0a60064981378302 | bash
curl -sL https://imputation.biodatacatalyst.nhlbi.nih.gov/get/928761/fda174b5b9f28e7e4de1d978164863ca8b77343fb0a39adbd38d29cdc7c31b1f | bash
curl -sL https://imputation.biodatacatalyst.nhlbi.nih.gov/get/928763/4dae10e75f79000a00304609b756e42ac95bad28156fe52579b63a1ebc278b67 | bash
curl -sL https://imputation.biodatacatalyst.nhlbi.nih.gov/get/928764/52954c5ca80c42fc3b4553bc230f6c577b866e2fdf5657db1450fbab25433597 | bash
cd -

mkdir -p zipped unzipped
while read -r chr passwd; do
  unzip -P $passwd zipped/$chr.zip -d unzipped
done < password.txt



#
## Post-imputation QC
#
for x in $projdir/outputs/genotypes/postimputation/TOPMed/unzipped/chr*.dose.vcf.gz; do
  echo Indexing $x ...
  bcftools index -f -t --threads 8 $x
done


# Concat all chroms
bcftools concat \
  -a \
  -O z \
  -o $projdir/outputs/genotypes/postimputation/quality_control/all.GRCh38.imputed.raw.vcf.gz \
  --threads 8 \
  $projdir/outputs/genotypes/postimputation/TOPMed/unzipped/{chr{1..22}.dose.vcf.gz,chrX.dose.vcf.gz}
bcftools index -f -t --threads 8 $projdir/outputs/genotypes/postimputation/quality_control/all.GRCh38.imputed.raw.vcf.gz


# Filter the variants
bcftools view \
  -v snps \
  -m 2 \
  -M 2 \
  -q 0.05 \
  -i 'R2>0.3' \
  -O z \
  -o $projdir/outputs/genotypes/postimputation/quality_control/all.GRCh38.imputed.maf5e-2.vcf.gz \
  $projdir/outputs/genotypes/postimputation/quality_control/all.GRCh38.imputed.raw.vcf.gz
bcftools index -f -t --threads 8 $projdir/outputs/genotypes/postimputation/quality_control/all.GRCh38.imputed.maf5e-2.vcf.gz


# Change the chromosome name to bare number
rm -f $projdir/outputs/genotypes/postimputation/quality_control/chr_map.txt
for x in {1..22} X; do
  echo -e "chr$x\t$x" >> $projdir/outputs/genotypes/postimputation/quality_control/chr_map.txt
done

bcftools annotate \
  -o $projdir/outputs/genotypes/postimputation/quality_control/all.GRCh38.imputed.maf5e-2.rm_chr.vcf.gz \
  --rename-chrs $projdir/outputs/genotypes/postimputation/quality_control/chr_map.txt \
  --threads 4 \
  $projdir/outputs/genotypes/postimputation/quality_control/all.GRCh38.imputed.maf5e-2.vcf.gz
bcftools index -f -t --threads 8 $projdir/outputs/genotypes/postimputation/quality_control/all.GRCh38.imputed.maf5e-2.rm_chr.vcf.gz


# Annotate the variants using rsID from dbsnp
awk '$0~/^#/{next} {print $1"\t"$2"\t"$3}' \
  /vol/projects/BIIM/resources/snpDB/GRCh38_b154/GRCh38_chrname.vcf \
  | bgzip >| $projdir/outputs/genotypes/postimputation/quality_control/GRCh38_b154.tsv.gz
tabix -b2 -e2 GRCh38_b154.tsv.gz

bcftools annotate \
  -c CHROM,POS,ID \
  -a $projdir/outputs/genotypes/postimputation/quality_control/GRCh38_b154.tsv.gz \
  -o $projdir/outputs/genotypes/postimputation/quality_control/all.GRCh38.imputed.maf5e-2.rm_chr.annot.vcf.gz \
  --threads 4 \
  $projdir/outputs/genotypes/postimputation/quality_control/all.GRCh38.imputed.maf5e-2.rm_chr.vcf.gz
bcftools index -f -t --threads 8 $projdir/outputs/genotypes/postimputation/quality_control/all.GRCh38.imputed.maf5e-2.rm_chr.annot.vcf.gz

# Cleanup
mv $projdir/outputs/genotypes/postimputation/quality_control/all.GRCh38.imputed.maf5e-2.rm_chr.annot.vcf.gz* $projdir/outputs/genotypes/postimputation/
# rm -fr $projdir/outputs/genotypes/postimputation/quality_control/*
# rm -fr $projdir/outputs/genotypes/postimputation/TOPMed/unzipped/*
