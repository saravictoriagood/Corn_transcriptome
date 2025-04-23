#!/bin/bash
#SBATCH --time=20:00:00
#SBATCH --account=def-saragood
#SBATCH --mem=124000M
#SBATCH --mail-user=dmitrii.trubetskoy@gmail.com
#SBATCH --mail-type=All

module load nixpkgs/16.09
module load fastqc/0.11.8

# UNZIPPING EVERYTHING that has .GZ extension!!!
#for file in $(find /home/saragood/scratch/Medaka_Fasting/Raw_reads/ -name '*.fastq.gz')
for file in $(find /scratch/saragood/dimitri/corn/test_site -name '*.fastq.gz')
do
        echo "~ ~ ~ ~ UNZIPPING.. $file"
        echo
        gunzip "$file"
done

echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ DONE UNZIPPING ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"

# RUNNING FASTQC on EVERYTHING that has .FASTQ extension!!!
#for file in $(find /home/saragood/scratch/Medaka_Fasting/Raw_reads/  -name '*.fastq')
#for file in $(find /home/saragood/scratch/saragood/dimitri/ -name  '*.fastq')
for file in $(find /scratch/saragood/dimitri/corn/test_site -name '*.fastq')
do
        echo "--------------FastQC analysis for $file started..."
        echo
        fastqc "$file" -o /home/saragood/scratch/dimitri/corn/FastQC_out_corn
done
echo
echo "---DONE---"
