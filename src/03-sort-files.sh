#!/bin/sh
# This script will craete ancestry directories and mv files to each directory

dir=/Users/Emma/GitHub/1kg/data
cd $dir
for ancestry in AFR AMR EAS EUR SAS; do
    echo "moving $ancestry files"
    mkdir -p $ancestry
    
    files=`ls | grep $ancestry`
    mv $files $ancestry
    num_files=`ls -1 $ancestry/* | wc -l`
    echo $num_files moved
    echo " "
done