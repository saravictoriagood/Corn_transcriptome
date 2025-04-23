#!/bin/bash
#SBATCH --time=30:30:00
#SBATCH --account=def-saragood
#SBATCH --mem=960000M
#SBATCH --mail-user=dmitrii.trubetskoy@gmail.com
#SBATCH --mail-type=All
#Created by Dimitri Trubetskoy (dmitrii.trubetskoy@gmail.com)
#Version 1.0 from 06.11.2019
#---------------------------

# VERY IMPORTANT! - First of all test this script with echo to be sure everything goes as expected
#!!! Attention !!! - test the list creation first and choose the right  number for the |head| part
#!!! Attention !!! - it is not yet tested, pay attention to  "paired/$line\_R2-paired" part

module load nixpkgs/16.09  intel/2018.3
module load hisat2/2.1.0
module load samtools/1.9

#Create a list of files
ls -ltr | awk '{print$9}' | sed 's/_R.*$//g' | uniq | head -19 | sed -n  '1!p' > list
#ls -ltr | awk '{print$9}' | sed 's/_R.*$//g' | uniq | sed -n '1!p' >  list

# Runs HISAT2 on trimmed samples
############# ~~~ TEST TEST TEST VERSION FOR L003 IMAM ~~~ #############
#while read -r line
#do
#echo $line\_R2-paired.fastq.gz
#echo "%%%%%%%%%%%%%%%%"
#echo "hisat2 -p 8 --dta -q --phred33 -x  /scratch/saragood/dimitri/L002_zipped/imam002_hisat_repo/imam002 \ # this is a repository with .ht2 files from hisat2 build index step
#       -1  /scratch/saragood/dimitri/L003_zipped/paired/${line}_R2-paired.fastq.gz  \ # this is where paired files are (R2)
#       -2  /scratch/saragood/dimitri/L003_zipped/paired/${line}_R1-paired.fastq.gz  \ # this is where paired files are (R1)
#       -S /scratch/saragood/dimitri/L003_zipped/sams/$line.sam" # this is where sams are going
#done < list
#rm list

############# ~~~ REAL VERSION FOR L003 IMAM ~~~ #############
# Runs HISAT2 on trimmed samples
while read -r line
do
#       echo $line\_R2-paired.fastq.gz
#       echo $line\_R1-paired.fastq.gz
#       echo "%%%%%%%%%%%%%%%%"


hisat2 -p 8 --dta -q --phred33 -x  /scratch/saragood/dimitri/corn/corn_hisat_repo/corn_hisat2_index \
        -1 /scratch/saragood/dimitri/corn/corno/paired/${line}_R2-paired.fastq.gz  \
        -2 /scratch/saragood/dimitri/corn/corno/paired/${line}_R1-paired.fastq.gz  \
        -S /scratch/saragood/dimitri/corn/sams/$line.sam

done < list
rm list
