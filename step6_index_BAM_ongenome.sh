#!/bin/bash
#SBATCH --time=10:00:00
#SBATCH --account=def-saragood
#SBATCH --mem=200000M
#SBATCH --mail-user=dmitrii.trubetskoy@gmail.com
#SBATCH --mail-type=All

module load nixpkgs/16.09  gcc/5.4.0
module load nixpkgs/16.09  gcc/7.3.0
module load nixpkgs/16.09  intel/2018.3
#module load stringtie/1.3.6
#module load stringtie/2.0
module load stringtie/1.3.4d

ls -ltr | awk '{print$9}' | sed 's/.bam//g' | sed -n '1!p' | head -n -2 > list

while read -r line

do
#       echo "stringtie -p 8 -G /scratch/saragood/dimitri/corn/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.gtf -o $line.gtf -l $line $line.bam"

#/home/saragood/scratch/dimitri/corn/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.gtf corn/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.gtf

#        stringtie -p 8 -G /scratch/saragood/dimitri/L002_zipped/hisat2020test/ncbi-genomes-2020-04-22/GCF_000001405.39_GRCh38.p13_genomic.gtf -o $line.gtf -l $line $line.bam

#stringtie -p 8 -G /home/saragood/scratch/dimitri/corn/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.gtf -o $line.gtf -l $line $line.bam
stringtie -p 8 -G GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.gff -o $line.gtf -l $line $line.bam

done < list

rm list
