#!/bin/bash
#SBATCH -D /global/scratch/users/pettinga/refs/nebraska/indexes/star/
#SBATCH -J star-index_nebraska
#SBATCH --partition=savio
#SBATCH --account=co_rosalind
#SBATCH --qos=rosalind_savio_normal
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --time=72:00:00
#SBATCH --mem=64000
#SBATCH -o /global/scratch/users/pettinga/logs/star-index_nebraska.%j
#SBATCH -e /global/scratch/users/pettinga/logs/star-index_nebraska.%j
#SBATCH --mail-user=pettinga@berkeley.edu
#SBATCH --mail-type=All

module load STAR/2.7.1a

STAR \
--runThreadN 8 \
--runMode genomeGenerate \
--genomeDir /global/scratch/users/pettinga/refs/nebraska/indexes/star \
--genomeFastaFiles /global/scratch/users/pettinga/refs/nebraska/sequence/HanXRQr2.0-SUNRISE-2.1_nebraska_consensus.fa \
--sjdbGTFfile /global/scratch/users/pettinga/refs/nebraska/annotation/HanXRQr2.0-SUNRISE-2.1_nebraska_consensus.gtf \
--sjdbOverhang 99
