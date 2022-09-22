#!/bin/sh
# This script will process the following steps:
#   1. download LD reference panels from the 1000 genomes project website for each chromosome
#   2. convert the downloaded vcf files to plink format
#   3. filter out duplicate SNPs and SNPs have long indels
# The original script was written by Minsoo Kim.



### create a directory
dir=/u/project/gandalm/shared/refGenomes/1000genomes
mkdir -p $dir/chr


### get reference fasta file from the 1000 genomes project website
cd $dir
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.fai

gunzip human_g1k_v37.fasta.gz


### move to the directory and download a vcf file for each chromosome
cd $dir/chr
for chrs in {1..22} X; do 
    wget http://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/b37.vcf/chr$chrs.1kg.phase3.v5a.vcf.gz
done


### create index file for each chromosome vcf file
# tools used: samtools(tabix, bcftools), plink
for chr in {1..22} X; do
    # index each vcf file
    tabix -f chr$chr.1kg.phase3.v5a.vcf.gz 
    
    # convert vcf file to plink format
    bcftools norm -Ou -m -any chr$chr.1kg.phase3.v5a.vcf.gz | \
    bcftools norm -Ou -f ../human_g1k_v37.fasta | \
    bcftools annotate -Ob -x ID -I +'%CHROM:%POS:%REF:%ALT' | \
    plink --bcf /dev/stdin --keep-allele-order --vcf-idspace-to _ --const-fid \
          --allow-extra-chr 0 --split-x b37 no-fail --make-bed --out kgp.chr$chr
done
plink --bfile kgp.chrX --keep-allele-order --impute-sex .95 .95 --make-bed --out kgp.chrX && /bin/rm kgp.chrX.{bed,bim,fam}~



### remove duplicate SNPs and long indel SNPs
# get duplicate SNP list
for chr in {{1..22},X,Y,MT}; do cut -f2 kgp.chr$chr.bim | sort | uniq -c | awk '$1>=2 {print $2}'; done > kgp.dups

# get long indel SNP list
cut -f2 kgp.chr{{1..22},X,Y,MT}.bim | awk 'length($1)>=150' | sort | uniq > kgp.longindels

# remove duplicate SNPs and long indel SNPs using the created lists
for chr in {1..22} X; do 
    cat kgp.{dups,longindels} | \
    plink --bfile kgp.chr$chr --keep-allele-order --exclude /dev/stdin --make-bed --out kgp.clean.chr$chr
done
# the kgp.clean.chr$chr files will be used as the base files for the next step (2.filter_ld_panels_for_each_ancestry.ipynb)

