#!/bin/bash


#SBATCH --time=25:00:00
#SBATCH --account=def-saragood
#SBATCH --mem=960000M
#SBATCH --mail-user=dmitrii.trubetskoy@gmail.com
#SBATCH --mail-type=ALL

module load nixpkgs/16.09
module load java/1.8.0_192
module load trimmomatic/0.36

for file in  *R1.fastq.gz;
do
BASENAME=${file%%_R1.fastq.gz*}
java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.36.jar PE -phred33 ${BASENAME}_R1.fastq.gz \
${BASENAME}_R2.fastq.gz \
/scratch/saragood/dimitri/corn/corno/paired/${BASENAME}_R1-paired.fastq.gz \
/scratch/saragood/dimitri/corn/corno/UNpaired/${BASENAME}_R1-unpaired.fastq.gz \
/scratch/saragood/dimitri/corn/corno/paired/${BASENAME}_R2-paired.fastq.gz \
/scratch/saragood/dimitri/corn/corno/UNpaired/${BASENAME}_R2-unpaired.fastq.gz \
ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10 LEADING:5 TRAILING:5 SLIDINGWINDOW:4:5 MINLEN:36

done





# -= STEP 2 script - trimmomatic, trimming sequences =-
# Created by Trubetskoy Dimitri, 15.11.2021 - dmitrii.trubetskoy@gmail.com

#/scratch/saragood/dimitri/corn/test/${BASENAME}_R1-paired.fastq.gz \
