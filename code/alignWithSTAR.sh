#!/bin/bash
#$ -cwd # happens in the current wd
#$ -pe smp 1
#$ -l h_vmem=64G
#$ -j y # error and output text together
#$ -o logs/ #concatenate output texts in logs folder
#$ -l h_rt=1:0:0
#$ -t 1-50

module load star/2.7.9a

SAMPLES=$(sed -n "${SGE_TASK_ID}p" /data/scratch/btx915/data/sampleList50.txt)
INPUT_FILE_1=$(echo "/data/scratch/btx915/data/rawfastq/${SAMPLES}_pass_1.fastq.gz")
INPUT_FILE_1=$(echo "/data/scratch/btx915/data/rawfastq/${SAMPLES}_pass_2.fastq.gz")

STAR --genomeDir /data/scratch/btx915/data/genome/STAR --readFilesIn $INPUT_FILE_1 $INPUT_FILE_2 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --outFileNamePrefix /data/scratch/btx915/data/alignments/${SAMPLES}
