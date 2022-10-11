#!/bin/sh
# This script will perform the following steps:
#   1. download 1000 Genomes Project phase 3 genotype array data for each chromosome (build hg37)
#   2. convert downloaded vcf files to plink format
#   3. filter out duplicate SNPs and long indels
#   4. merge all chromosome files

# Softwares used: 
#   1. tabix (v.1.16): http://www.htslib.org/doc/tabix.html
#   2. bcftools (v.1.16): https://samtools.github.io/bcftools/howtos/index.html
#   3. plink (v.1.90beta): https://www.cog-genomics.org/plink/

# Output:
#   1. kgp.bed, kgp.bim, kgp.fam files 

### Create a directory
dir=/u/project/gandalm/shared/refGenomes/1000genomes_GeneticsMakie_chrpos
mkdir -p $dir/chrs

### Download hg37 fasta file
cd $dir
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.fai
gunzip human_g1k_v37.fasta.gz

### Download vcf files for each chromosome
cd $dir/chrs
for chr in {1..22} X; do 
    wget https://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/b37.vcf/chr$chr.1kg.phase3.v5a.vcf.gz
done

### Create tabix indices for each vcf file
# get path to the software
tabix=/u/project/gandalm/shared/apps/bcftools/bcftools-1.16/htslib-1.16/tabix
bcftools=/u/project/gandalm/shared/apps/bcftools/bcftools-1.16/bcftools
plink=/u/project/gandalm/shared/apps/plink1/v1.90beta-20220402/plink

# note: "--threads 16" can be omitted
for chr in {1..22} X; do
    echo "processing chr$chr"
    echo "tabix: creating tabix indices for vcf file" 
    $tabix -f chr$chr.1kg.phase3.v5a.vcf.gz 
    
    echo "bcftools & plink: converting vcf to plink format"
    $bcftools norm --threads 16 -Ou -m -any chr$chr.1kg.phase3.v5a.vcf.gz | \
    $bcftools norm --threads 16 -Ou -f ../human_g1k_v37.fasta | \
    $bcftools annotate --threads 16 -Ob -x ID -I +'%CHROM:%POS:%REF:%ALT' | \
    $plink --bcf /dev/stdin --keep-allele-order --vcf-idspace-to _ --const-fid \
           --allow-extra-chr 0 --split-x b37 no-fail --make-bed --out kgp.chr$chr
    echo " "
done
$plink --bfile kgp.chrX --keep-allele-order --impute-sex .95 .95 --make-bed --out kgp.chrX && /bin/rm kgp.chrX.{bed,bim,fam}~

### Remove duplicate and long indel SNPs
# Get duplicate SNP list
for chr in {{1..22},X}; do cut -f2 kgp.chr$chr.bim | sort | uniq -c | awk '$1>=2 {print $2}'; done > kgp.dups

# Get long indel SNP list
cut -f2 kgp.chr{{1..22},X}.bim | awk 'length($1)>=150' | sort | uniq > kgp.longindels

# Remove duplicate and long indel SNPs using the above lists
for chr in {1..22} X; do 
    cat kgp.{dups,longindels} | \
    $plink --bfile kgp.chr$chr --keep-allele-order --exclude /dev/stdin --make-bed --out kgp.clean.chr$chr
done

### merge all kgp.clean.chr$chr plink files
cat kgp.clean.chrX.fam > $dir/kgp.fam
cat kgp.clean.chr{{1..22},X}.bim > $dir/kgp.bim
(echo -en "\x6C\x1B\x01"; tail -qc +4 kgp.clean.chr{{1..22},X}.bed) > $dir/kgp.bed

