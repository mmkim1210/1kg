# 1000 genomes project LD reference panels by ancestry groups

This directory stores scripts that prepared LD reference panels by ancestry groups to use for the `GeneticsMakie.jl`. The LD panels were downloaded from the [1000 genomes project (Phase3)](https://www.internationalgenome.org/) website.  

The notation of the ancestry groups is as follows:  
- AFR: African-American
- AMR: Admixture-American
- EAS: East Asian
- EUR: European
- SAS: South Asian

### Short descriptions of the scripts and the logs directory

- `1.download_and_prep_1kg_ld_panels.sh`: downloads and removes unnecessary SNPs (e.g. duplicated SNPs, SNPs have 
long indels)
- `2.filter_ld_panels_for_each_ancestry.ipynb`: prepares LD reference panels for each ancestry group with maf > 
0.05 and mac > 1 using `Julia`
- `3.optional_sort_files.sh`: optionally, this script sorts the files by ancestry groups
- `logs` directory: for each chromosome, the log file shows the number of samples and SNPs before and after the process for each ancestry group
