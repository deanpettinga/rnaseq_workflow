#!/bin/bash
#SBATCH -J rnaseq-workflow
#SBATCH --partition=savio
#SBATCH --account=co_rosalind
#SBATCH --qos=rosalind_savio_normal
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=72:00:00
#SBATCH --mem=8000
#SBATCH -o logs/runs/rnaseq-workflow.%j.o
#SBATCH -e logs/runs/rnaseq-workflow.%j.e
#SBATCH --mail-user=pettinga@berkeley.edu
#SBATCH --mail-type=All

source ~/.bashrc
conda activate snakemake

# # specify which conda installation to use
# conda_setup='/secondary/projects/bbc/tools/miniconda3/etc/profile.d/conda.sh'
#
# # this make 'conda' callable and allows conda environments to be created.
# source $conda_setup

# save DAG job file with time stamp
TIME=$(date "+%Y-%m-%d_%H.%M.%S")

# make logs dir if it does not exist already. Without this, logs/ is automatically generate only after the first run of the pipeline
logs_dir="logs/runs"
[[ -d $logs_dir ]] || mkdir -p $logs_dir

snakemake --use-conda -n > logs/runs/workflow_${TIME}.txt
snakemake --dag | dot -Tpng > logs/runs/workflow_${TIME}.png


snakemake \
--use-conda \
--jobs \
--cluster "sbatch \
-J star-index \
--partition=savio \
--account=co_rosalind \
--qos=rosalind_savio_normal \
--nodes={resources.nodes}
--cpus-per-task={threads} \
--mem={resources.mem_gb}000 \
-l walltime=72:00:00 \
-o logs/runs/ \
-e logs/runs/"
