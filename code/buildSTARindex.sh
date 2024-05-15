#!/bin/bash
#$ -pe smp 48
#$ -l h_vmem=1G
#$ -l h_rt=1:0:0
#$ -cwd
#$ -j y

module load star/2.7.9a

STAR --runThreadN 48 --runMode genomeGenerate --genomeDir /data/scratch/btx915/data/genome/STAR/ --genomeFastaFiles /data/scratch/btx915/data/genome/Mus_musculus.GRCm39.dna.primary_assembly.fa --sjdbGTFfile /data/scratch/btx915/data/genome/Mus_musculus.GRCm39.111.gtf
