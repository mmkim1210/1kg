#!/bin/sh
# This script will pursue the following steps:
#   1. download 1000 Genomes Project phase 3 genotype array data for each chromosome (build hg37)
#   2. convert downloaded vcf files to plink format
#   3. filter out duplicate SNPs and long indels

# Softwares used: 
#   1. tabix: http://www.htslib.org/doc/tabix.html
#   2. bcftools: https://samtools.github.io/bcftools/howtos/index.html
#   3. plink: https://www.cog-genomics.org/plink/

# Output:
#   1. kgp.bed, kgp.bim, kgp.fam files with all chr1-22 and chrX files merged

### Create a directory
dir=/u/project/gandalm/shared/refGenomes/1000genomes
mkdir -p $dir/chrs

### Download hg37 fasta file
cd $dir
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.fai
gunzip human_g1k_v37.fasta.gz

### Downnload vcf files for each chromosome
cd $dir/chr
for chrs in {1..22} X; do 
    wget https://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/b37.vcf/chr$chr.1kg.phase3.v5a.vcf.gz
done

### Create tabix indices for each vcf file
# get path to the software
# bcftools & plink
. /u/local/Modules/default/init/modules.sh
module load bcftools/1.11 plink/1.90b624
# tabix
tabix=/u/project/gandalm/shared/apps/tabix-0.2.6/tabix

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

### Remove duplicate and long indel SNPs
# Get duplicate SNP list
for chr in {{1..22},X,Y,MT}; do cut -f2 kgp.chr$chr.bim | sort | uniq -c | awk '$1>=2 {print $2}'; done > kgp.dups

# Get long indel SNP list
cut -f2 kgp.chr{{1..22},X,Y,MT}.bim | awk 'length($1)>=150' | sort | uniq > kgp.longindels

# Remove duplicate and long indel SNPs using the above lists
for chr in {1..22} X; do 
    cat kgp.{dups,longindels} | \
    plink --bfile kgp.chr$chr --keep-allele-order --exclude /dev/stdin --make-bed --out kgp.clean.chr$chr
done

### merge all kgp.clean.chr$chr plink files
cat kgp.clean.chrX.fam > $dir/chrs/kgp.fam
cat kgp.clean.chr{{1..22},X}.bim > $dir/chrs/kgp.bim
(echo -en "\x6C\x1B\x01"; tail -qc +4 kgp.clean.chr{{1..22},X}.bed) > $dir/chrs/kgp.bed

