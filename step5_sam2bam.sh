#!/bin/bash
#!!! Attention !!! - test the list creation first and choose the right  number for the |head| part

#SBATCH --time=10:00:00
#SBATCH --account=def-saragood
#SBATCH --mem=200000M
#SBATCH --mail-user=dmitrii.trubetskoy@gmail.com
#SBATCH --mail-type=All

module load nixpkgs/16.09  intel/2018.3
module load samtools/1.9

#ls -ltr | awk '{print$9}' | sed 's/.sam//g' | head -n -2 | sed -n '1!p'  > list
ls -ltr | awk '{print$9}' | sed 's/.sam//g' | head -n -1 > list

while read -r line
do
#       echo "samtools sort -o $line.bam $line.sam"
        samtools sort -o $line.bam $line.sam
        done < list
        rm list
