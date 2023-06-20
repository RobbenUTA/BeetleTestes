#!/bin/bash
#SBATCH --nodes=1 --ntasks-per-node=42            ######1 Node with 64 cores
#SBATCH --mem=50000                             ########16 GB of memory
#SBATCH --time=148:00:00
#SBATCH --job-name=call_snps
#SBATCH --error=myjob.%j.errors
#SBATCH --mail-type=END
#SBATCH --mail-user=Michael.robben@uta.edu
#####Type Code below#########
source activate Test #not sure why source works but conda doesn't but whatevs


BAM="/home/data/beetle/Demetra/10x9/outs/possorted_genome_bam.bam" #Bam file name
GENOME="/home/data/beetle/Demetra/tribolium_genome/fasta/genome.fa" #Genome file directory

cd /home/robbenm/LuberLab/Beetle/SNP/

echo 'start'
#samtools mpileup -g -f $GENOME $BAM > ./raw.bcf
#bcftools view -bvcg ./raw.bcf > ./var.bcf
#bcftools view ./var.bcf | vcfutils.pl varFilter - > ./var-final.vcf
#bcftools call --threads 42 -mv ./raw.bcf
bcftools mpileup --threads 42 -Q 30 -A -x -Ou -f $GENOME $BAM | bcftools call --threads 42 -mv > ./sample.vcf
echo 'done'