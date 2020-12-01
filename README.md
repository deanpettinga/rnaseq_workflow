# Sunflower Domestication RNAseq

## Author

* Dean Pettinga [@deanpettinga](https://github.com/deanpettinga)

## Analysis Outline
1. prepare STAR indexes for Arikara and Nebraska
  - copies corrected HanXRQr2 genomes into `references/`
  - makes STAR indexes for each genome for alignments
2. raw read QC
  - Trim_Galore on paired-end reads (min Quality =20)
  - fastQC on each library
3. Alignment
  - STAR alignment of adapter and quality-trimmed paired end libraries to reference genomes (arikara used for all domesticated samples, nebraska used for all wild samples)
    - output are BAM sorted by coordinates, indexed with SAMtools.
4. MultiQC
  - aggregated results from Trim_Galore (Cutadapt,FastQC) and STAR alignments into one page.
5. Differential Expression
  - this was edited and run in RStudio, outside of this snakemake automation workflow
  - analysis script in `bin/diffExp.Rmd`
    +  makes extensive use of `vari-bbc/bbcRNA` package wrapper functions for differential expression with edgeR.
    + some paths within this doc refer to local machine directory structure, will need modification to reproduce analysis (in particular, take note when importing STAR counts and STAR alignment metrics)
  - the HTML report knitted from the .Rmd will contain the original code in collapsable boxes for reference.
