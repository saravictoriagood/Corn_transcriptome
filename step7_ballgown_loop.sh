#!/bin/bash
#SBATCH --time=15:00:00
#SBATCH --account=def-saragood
#SBATCH --mem=200000M
#SBATCH --mail-user=dmitrii.trubetskoy@gmail.com
#SBATCH --mail-type=All


module load nixpkgs/16.09  intel/2018.3
module load stringtie/1.3.4d


ls -ltr | awk '{print$9}' | sed 's/.bam//g' | sed -n '1!p' | head -n -3 > list


while read -r line

do

stringtie -e -B -p 8 -G GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.gff -o /scratch/saragood/dimitri/corn/ballgown/$line/$line.gtf $line.bam

#stringtie -e -B -p 8 -G /scratch/saragood/dimitri/corn/bams/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.gff \
#                -o /scratch/saragood/dimitri/corn/ballgown/$line/$line.gtf \
#                   /scratch/saragood/dimitri/corn/bams/$line.bam

done < list

rm list


