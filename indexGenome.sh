#!/bin/bash

#SBATCH -p bluemoon # partition name
#SBATCH -t 0-20:00 # hours:minutes runlimit after which job will be killed
#SBATCH -c 6 # number of cores requested -- this needs to be greater than or equal to the number of cores you plan to use to run your job
#SBATCH --mem 16G
#SBATCH --job-name STAR_index # Job name
#SBATCH --output=%x_%j.out.txt                                                          
#SBATCH --mail-user=rbarrant@uvm.edu
#SBATCH --mail-type=ALL

cd /users/r/b/rbarrant/projects/nfcorePipelines/eclipseq/genome

singularity exec --bind /users/r/b/rbarrant/projects/nfcorePipelines/eclipseq:/eclipseq docker://brianyee/star:2.7.6a STAR --runThreadN 6 \
--runMode genomeGenerate \
--genomeDir /eclipseq/genome \
--genomeFastaFiles /eclipseq/genome/hg38.fa \
--sjdbGTFfile /eclipseq/genome/hg38.knownGene.gtf \
--sjdbOverhang 99
