#!/bin/bash
#SBATCH --partition=bigmemwk
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=256G   
#SBATCH --time=150:00:00
#SBATCH --job-name=eclip
#SBATCH --output=%x_%j.out.txt                                                          
#SBATCH --mail-user=rbarrant@uvm.edu
#SBATCH --mail-type=ALL

export sampleSheet=$1
export outDir=$2
workdir=${sampleSheet##*/}
workdir=./deleteme_work_${workdir%.*}
module load singularity/3.7.1
source ~/.bashrc
mamba activate env_nf

echo $sampleSheet
echo $outDir
echo $workdir
export NXF_SINGULARITY_CACHEDIR=/gpfs1/mbsr_tools/NXF_SINGULARITY_CACHEDIR
nextflow run main.nf -profile singularity --input ${sampleSheet} --outdir ${outDir} -work-dir ${workdir}

#nextflow clean -f
mamba deactivate
