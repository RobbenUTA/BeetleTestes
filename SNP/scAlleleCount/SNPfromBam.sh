#!/bin/bash
#SBATCH --nodes=1 --ntasks-per-node=12            ######1 Node with 64 cores
#SBATCH --mem=50000                             ########16 GB of memory
#SBATCH --time=148:00:00
#SBATCH --job-name=bam_snps
#SBATCH --error=myjob.%j.errors
#SBATCH --mail-type=END
#SBATCH --mail-user=Michael.robben@uta.edu
#####Type Code below#########
source activate scsnp #not sure why source works but conda doesn't but whatevs


BAM="/home/data/beetle/Demetra/10x9/outs/possorted_genome_bam.bam"
BARCODES="/home/robbenm/LuberLab/Beetle/CellRanger_out/feature/barcodes.tsv"
cd /home/robbenm/LuberLab/Beetle/SNP/scAlleleCount/

echo "check paths:"
echo "$PWD"
echo $BAM
echo $BARCODES

echo 'start'

./scAlleleCount/scAlleleCount.py -v --snps ../snps_fixed.txt --barcodes $BARCODES --bamfile $BAM --output-format mm --max-depth 9999999 --output-prefix out

echo 'done'