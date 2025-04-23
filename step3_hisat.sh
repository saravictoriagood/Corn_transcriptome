#!/bin/bash
#SBATCH --time=05:30:00
#SBATCH --account=def-saragood
#SBATCH --mem=200000M
#SBATCH --mail-user=dmitrii.trubetskoy@gmail.com
#SBATCH --mail-type=All

module load nixpkgs/16.09  intel/2018.3
module load hisat2/2.1.0

hisat2-build --ss corn.ss --exon corn.exon GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.fna corn_hisat2_index
