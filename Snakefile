import pandas as pd
import os
from shutil import which
from snakemake.utils import validate, min_version
##### set minimum snakemake version #####
min_version("5.14.0")

##### check if conda is available. If not, we don't run rules requiring conda envs. Applies to dryruns and actual runs.
conda_avail = (which('conda') is not None)

##### load config and sample sheets #####

configfile: "bin/config.yaml"
#validate(config, schema="schemas/config.schema.yaml")

units = pd.read_table(config["units"]).set_index("sample", drop=False)

var_calling_units = pd.read_table("bin/variant_calling_units.tsv").set_index("unit", drop=False)

contrasts = pd.read_table(config["contrasts"]).set_index("name", drop=False)
#validate(contrasts, schema="schemas/contrasts.schema.yaml")

if not os.path.exists('tmp'):
    os.mkdir('tmp')

fai_file = config["ref"]["fai"]
contigs_file = "bin/grouped_contigs.tsv"

# read in file with col1 as contig group name and col2 as a commma separated list of contigs. This was written for use with GATK for splitting the variant calling steps by chromosome/contig. We group the unplaced contigs together since those tend to be small.
# we read in this file even if not calling variants. Otherwise variant-calling rules relying on this file will cause an error.
contig_groups = pd.read_table(contigs_file)

# if config set to run variants component of pipeline,  check that chromosomes were parsed correctly in the grouped_Contigs.tsv file.
if config["call_variants"]:

    # check chromosomes/contigs parsed correctly.
    chroms = contig_groups['contigs'].values[0:(contig_groups.shape[0]-1)].tolist()
    unanchored_contigs = contig_groups['contigs'].values[contig_groups.shape[0]-1].split(",")
    contig_groups_flattened = chroms + unanchored_contigs

    fai_pd = pd.read_table(fai_file, header=None)

    assert contig_groups_flattened == fai_pd.iloc[:,0].values.tolist(), "Chromosomes in .fai do not match those in 'grouped_contigs.tsv'."


##### target rules #####

rule all:
    input:
        #lambda wildcards: expand("analysis/05_haplotypecaller/{units.sample}.Aligned.sortedByCoord.out.addRG.mrkdup.splitncigar.baserecal.vcf.gz", units=var_calling_units.itertuples()) if config["call_variants"] else [],
        #lambda wildcards: expand("analysis/07_jointgeno/all.{contig_group}.vcf.gz", contig_group=contig_groups.name) if config["call_variants"] else [],
        #lambda wildcards: "analysis/08_merge_and_filter/all.merged.filt.vcf.gz" if config["call_variants"] else [],
        lambda wildcards: "analysis/09a_variant_annot/all.merged.filt.PASS.snpeff.vcf.gz.tbi" if config["call_variants"] else [],
        "analysis/09b_snp_pca_and_dendro/report.html" if (config["call_variants"] and conda_avail) else [],

        # # copy_references
        # "references/HanXRQr2.0-SUNRISE-2.1.genome.fasta",
        # "references/HanXRQr2.0-SUNRISE-2.1.genome.gff3",
        # "references/HanXRQr2.0-SUNRISE-2.1_nebraska_consensus.fa",
        # "references/HanXRQr2.0-SUNRISE-2.1_nebraska_consensus.gtf",
        # "references/HanXRQr2.0-SUNRISE-2.1_arikara_consensus.fa",
        # "references/HanXRQr2.0-SUNRISE-2.1_arikara_consensus.gtf",

        # # star_index
        # expand("references/star_indexes/{ref}/Genome", ref=["arikara","nebraska"]),
        # expand("references/star_indexes/{ref}/Log.out", ref=["arikara","nebraska"]),
        # expand("references/star_indexes/{ref}/SA", ref=["arikara","nebraska"]),
        # expand("references/star_indexes/{ref}/SAindex", ref=["arikara","nebraska"]),
        # expand("references/star_indexes/{ref}/chrLength.txt", ref=["arikara","nebraska"]),
        # expand("references/star_indexes/{ref}/chrName.txt", ref=["arikara","nebraska"]),
        # expand("references/star_indexes/{ref}/chrNameLength.txt", ref=["arikara","nebraska"]),
        # expand("references/star_indexes/{ref}/chrStart.txt", ref=["arikara","nebraska"]),
        # expand("references/star_indexes/{ref}/exonGeTrInfo.tab", ref=["arikara","nebraska"]),
        # expand("references/star_indexes/{ref}/exonInfo.tab", ref=["arikara","nebraska"]),
        # expand("references/star_indexes/{ref}/geneInfo.tab", ref=["arikara","nebraska"]),
        # expand("references/star_indexes/{ref}/genomeParameters.txt", ref=["arikara","nebraska"]),
        # expand("references/star_indexes/{ref}/sjdbInfo.txt", ref=["arikara","nebraska"]),
        # expand("references/star_indexes/{ref}/sjdbList.fromGTF.out.tab", ref=["arikara","nebraska"]),
        # expand("references/star_indexes/{ref}/sjdbList.out.tab", ref=["arikara","nebraska"]),
        # expand("references/star_indexes/{ref}/transcriptInfo.tab", ref=["arikara","nebraska"]),

        # # fastq_screen_get_genomes
        # expand("analysis/fastq_screen/FastQ_Screen_Genomes/Arabidopsis/Arabidopsis_thaliana.TAIR10.{suffix}", suffix=["1.bt2","2.bt2","3.bt2","4.bt2","rev.1.bt2","rev.2.bt2"]),
        # expand("analysis/fastq_screen/FastQ_Screen_Genomes/Adapters/Contaminants.{suffix}", suffix=["1.bt2","2.bt2","3.bt2","4.bt2","rev.1.bt2","rev.2.bt2"]),
        # expand("analysis/fastq_screen/FastQ_Screen_Genomes/Drosophila/BDGP6.{suffix}", suffix=["1.bt2","2.bt2","3.bt2","4.bt2","rev.1.bt2","rev.2.bt2"]),
        # expand("analysis/fastq_screen/FastQ_Screen_Genomes/E_coli/Ecoli.{suffix}", suffix=["1.bt2","2.bt2","3.bt2","4.bt2","rev.1.bt2","rev.2.bt2"]),
        # expand("analysis/fastq_screen/FastQ_Screen_Genomes/Human/Homo_sapiens.GRCh38.{suffix}", suffix=["1.bt2","2.bt2","3.bt2","4.bt2","rev.1.bt2","rev.2.bt2"]),
        # expand("analysis/fastq_screen/FastQ_Screen_Genomes/Lambda/Lambda.{suffix}", suffix=["1.bt2","2.bt2","3.bt2","4.bt2","rev.1.bt2","rev.2.bt2"]),
        # expand("analysis/fastq_screen/FastQ_Screen_Genomes/Lambda/Vectors.{suffix}", suffix=["1.bt2","2.bt2","3.bt2","4.bt2","rev.1.bt2","rev.2.bt2"]),
        # expand("analysis/fastq_screen/FastQ_Screen_Genomes/Mitochondria/mitochondria.{suffix}", suffix=["1.bt2","2.bt2","3.bt2","4.bt2","rev.1.bt2","rev.2.bt2"]),
        # expand("analysis/fastq_screen/FastQ_Screen_Genomes/Mouse/Mus_musculus.GRCm38.{suffix}", suffix=["1.bt2","2.bt2","3.bt2","4.bt2","rev.1.bt2","rev.2.bt2"]),
        # expand("analysis/fastq_screen/FastQ_Screen_Genomes/PhiX/phi_plus_SNPs.{suffix}", suffix=["1.bt2","2.bt2","3.bt2","4.bt2","rev.1.bt2","rev.2.bt2"]),
        # expand("analysis/fastq_screen/FastQ_Screen_Genomes/Rat/Rnor_6.0.{suffix}", suffix=["1.bt2","2.bt2","3.bt2","4.bt2","rev.1.bt2","rev.2.bt2"]),
        # expand("analysis/fastq_screen/FastQ_Screen_Genomes/Vectors/Vectors.{suffix}", suffix=["1.bt2","2.bt2","3.bt2","4.bt2","rev.1.bt2","rev.2.bt2"]),
        # expand("analysis/fastq_screen/FastQ_Screen_Genomes/Worm/Caenorhabditis_elegans.WBcel235.{suffix}", suffix=["1.bt2","2.bt2","3.bt2","4.bt2","rev.1.bt2","rev.2.bt2"]),
        # expand("analysis/fastq_screen/FastQ_Screen_Genomes/Yeast/Saccharomyces_cerevisiae.R64-1-1.{suffix}", suffix=["1.bt2","2.bt2","3.bt2","4.bt2","rev.1.bt2","rev.2.bt2"]),
        # expand("analysis/fastq_screen/FastQ_Screen_Genomes/rRNA/GRCm38_rRNA.{suffix}", suffix=["1.bt2","2.bt2","3.bt2","4.bt2","rev.1.bt2","rev.2.bt2"]),
        # "analysis/fastq_screen/FastQ_Screen_Genomes/fastq_screen.conf",

        # # download_SILVA
        # expand("analysis/fastq_screen/FastQ_Screen_Genomes/SILVA/SILVA_138.1_{subunit}Ref_NR99_tax_silva.fasta", subunit=["SSU","LSU"]),

        # # fastq_screen_SILVA_db
        # expand("analysis/fastq_screen/FastQ_Screen_Genomes/SILVA/SILVA.{suffix}", suffix=["1.bt2","2.bt2","3.bt2","4.bt2","rev.1.bt2","rev.2.bt2"]),

        # # fastq_screen
        # expand("analysis/fastq_screen/{units.sample}-R1_screen.html", units=units.itertuples()),
        # expand("analysis/fastq_screen/{units.sample}-R1_screen.txt", units=units.itertuples()),
        # expand("analysis/fastq_screen/{units.sample}-R2_screen.html", units=units.itertuples()),
        # expand("analysis/fastq_screen/{units.sample}-R2_screen.txt", units=units.itertuples()),

        # mergeLanesAndRename
        # expand("raw_data/{units.sample}-R1.fastq.gz", units=units.itertuples()),
        # expand("raw_data/{units.sample}-R2.fastq.gz", units=units.itertuples()),\

        # Trim_Galore
        # expand("analysis/trimmed_data/{units.sample}_R{read}_val_{read}.fq.gz", read=[1,2], units=units.itertuples()),
        # expand("analysis/trimmed_data/{units.sample}-R{read}_val_{read}_fastqc.html", read=[1,2], units=units.itertuples()),
        # expand("analysis/trimmed_data/{units.sample}-R{read}_val_{read}_fastqc.zip", read=[1,2], units=units.itertuples()),
        # expand("analysis/trimmed_data/{units.sample}-R{read}.fastq.gz_trimming_report.txt", read=[1,2], units=units.itertuples()),

        # # STAR alignment
        # expand("analysis/star/{units.sample}.Aligned.sortedByCoord.out.bam", units=units.itertuples()),
        # expand("analysis/star/{units.sample}.Log.out", units=units.itertuples()),

        # multiQC
        "analysis/multiqc/multiqc_report.html",

        # # star_dedup
        # expand("analysis/star_dedup/{units.sample}.Aligned.sortedByCoord.out.dedup.bam", units=units.itertuples()),
        # expand("analysis/star_dedup/{units.sample}.Aligned.sortedByCoord.out.dedup.bam.bai", units=units.itertuples()),
        #
        # # dedup_2_fastq
        # expand("analysis/star_dedup/{units.sample}-unpaired.fastq.gz", units=units.itertuples()),
        # expand("analysis/star_dedup/{units.sample}-R1.fastq.gz", units=units.itertuples()),
        # expand("analysis/star_dedup/{units.sample}-R2.fastq.gz", units=units.itertuples()),

        #expand("analysis/02_splitncigar/{units.sample}.Aligned.sortedByCoord.out.addRG.mrkdup.splitncigar.bam", units=var_calling_units.itertuples())
        # edgeR
        #"bin/diffExp.html",

rule copy_references:
    input:
    output:
        HanXRQr2_fasta = "references/HanXRQr2.0-SUNRISE-2.1.genome.fasta",
        HanXRQr2_gff3 =  "references/HanXRQr2.0-SUNRISE-2.1.genome.gff3",
        nebraska_fasta = "references/HanXRQr2.0-SUNRISE-2.1_nebraska_consensus.fa",
        nebraska_gtf =   "references/HanXRQr2.0-SUNRISE-2.1_nebraska_consensus.gtf",
        arikara_fasta =  "references/HanXRQr2.0-SUNRISE-2.1_arikara_consensus.fa",
        arikara_gtf =    "references/HanXRQr2.0-SUNRISE-2.1_arikara_consensus.gtf",
    resources:
        threads = 1,
        nodes =   1,
        mem_gb =  16,
    shell:
        """
        # HanXRQr2
        cp /clusterfs/vector/instrumentData/blackmanlab/Helianthus/HanXRQr2/HanXRQr2.0-SUNRISE-2.1.genome.fasta {output.HanXRQr2_fasta}
        cp /clusterfs/vector/instrumentData/blackmanlab/Helianthus/HanXRQr2/HanXRQr2.0-SUNRISE-2.1.gff3 {output.HanXRQr2_gff3}
        # nebraska
        cp /clusterfs/vector/instrumentData/blackmanlab/Helianthus/nina/Reference_Indexes/NB_Consensus.fa {output.nebraska_fasta}
        cp /clusterfs/vector/instrumentData/blackmanlab/Helianthus/nina/Reference_Indexes/NB_CrossMap_LiftedOver.gtf {output.nebraska_gtf}
        # arikara
        cp /clusterfs/vector/instrumentData/blackmanlab/Helianthus/nina/Reference_Indexes/Arikara_Consensus.fa {output.arikara_fasta}
        cp /clusterfs/vector/instrumentData/blackmanlab/Helianthus/nina/Reference_Indexes/Arikara_CrossMap_LiftedOver.gtf {output.arikara_gtf}
        """

rule STAR_index:
    input:
        fasta = "references/HanXRQr2.0-SUNRISE-2.1_{ref}_consensus.fa",
        gtf =   "references/HanXRQr2.0-SUNRISE-2.1_{ref}_consensus.gtf",
    params:
        genome_dir = "references/star_indexes/{ref}/",
    output:
        "references/star_indexes/{ref}/Genome",
        "references/star_indexes/{ref}/Log.out",
        "references/star_indexes/{ref}/SA",
        "references/star_indexes/{ref}/SAindex",
        "references/star_indexes/{ref}/chrLength.txt",
        "references/star_indexes/{ref}/chrName.txt",
        "references/star_indexes/{ref}/chrNameLength.txt",
        "references/star_indexes/{ref}/chrStart.txt",
        "references/star_indexes/{ref}/exonGeTrInfo.tab",
        "references/star_indexes/{ref}/exonInfo.tab",
        "references/star_indexes/{ref}/geneInfo.tab",
        "references/star_indexes/{ref}/genomeParameters.txt",
        "references/star_indexes/{ref}/sjdbInfo.txt",
        "references/star_indexes/{ref}/sjdbList.fromGTF.out.tab",
        "references/star_indexes/{ref}/sjdbList.out.tab",
        "references/star_indexes/{ref}/transcriptInfo.tab",
    resources:
        threads = 8,
        nodes =   1,
        mem_gb =  64,
    log:
                 "logs/STAR_index/{ref}.log"
    conda:       "envs/star.yaml"
    shell:
        """
        STAR \
        --runThreadN {resources.threads} \
        --runMode genomeGenerate \
        --genomeDir {params.genome_dir} \
        --genomeFastaFiles {input.fasta} \
        --sjdbGTFfile {input.gtf} \
        --sjdbOverhang 99 \
        --outFileNamePrefix {params.genome_dir}
        """
rule mergeLanesAndRename_PE:
    input:
    output:      "raw_data/{sample}-R1.fastq.gz",
                 "raw_data/{sample}-R2.fastq.gz"
    log:         "logs/mergeLanesAndRename/mergeLanesAndRename_PE-{sample}.log"
    resources:
        threads = 1,
        nodes =   1,
        mem_gb =  16,
    conda:       "envs/R.yaml"
    script:      "bin/mergeLanesAndRename.R"


rule fastq_screen_get_genomes:
    input:
    output:      expand("analysis/fastq_screen/FastQ_Screen_Genomes/Arabidopsis/Arabidopsis_thaliana.TAIR10.{suffix}", suffix=["1.bt2","2.bt2","3.bt2","4.bt2","rev.1.bt2","rev.2.bt2"]),
                 expand("analysis/fastq_screen/FastQ_Screen_Genomes/Adapters/Contaminants.{suffix}", suffix=["1.bt2","2.bt2","3.bt2","4.bt2","rev.1.bt2","rev.2.bt2"]),
                 expand("analysis/fastq_screen/FastQ_Screen_Genomes/Drosophila/BDGP6.{suffix}", suffix=["1.bt2","2.bt2","3.bt2","4.bt2","rev.1.bt2","rev.2.bt2"]),
                 expand("analysis/fastq_screen/FastQ_Screen_Genomes/E_coli/Ecoli.{suffix}", suffix=["1.bt2","2.bt2","3.bt2","4.bt2","rev.1.bt2","rev.2.bt2"]),
                 expand("analysis/fastq_screen/FastQ_Screen_Genomes/Human/Homo_sapiens.GRCh38.{suffix}", suffix=["1.bt2","2.bt2","3.bt2","4.bt2","rev.1.bt2","rev.2.bt2"]),
                 expand("analysis/fastq_screen/FastQ_Screen_Genomes/Lambda/Lambda.{suffix}", suffix=["1.bt2","2.bt2","3.bt2","4.bt2","rev.1.bt2","rev.2.bt2"]),
                 expand("analysis/fastq_screen/FastQ_Screen_Genomes/Lambda/Vectors.{suffix}", suffix=["1.bt2","2.bt2","3.bt2","4.bt2","rev.1.bt2","rev.2.bt2"]),
                 expand("analysis/fastq_screen/FastQ_Screen_Genomes/Mitochondria/mitochondria.{suffix}", suffix=["1.bt2","2.bt2","3.bt2","4.bt2","rev.1.bt2","rev.2.bt2"]),
                 expand("analysis/fastq_screen/FastQ_Screen_Genomes/Mouse/Mus_musculus.GRCm38.{suffix}", suffix=["1.bt2","2.bt2","3.bt2","4.bt2","rev.1.bt2","rev.2.bt2"]),
                 expand("analysis/fastq_screen/FastQ_Screen_Genomes/PhiX/phi_plus_SNPs.{suffix}", suffix=["1.bt2","2.bt2","3.bt2","4.bt2","rev.1.bt2","rev.2.bt2"]),
                 expand("analysis/fastq_screen/FastQ_Screen_Genomes/Rat/Rnor_6.0.{suffix}", suffix=["1.bt2","2.bt2","3.bt2","4.bt2","rev.1.bt2","rev.2.bt2"]),
                 expand("analysis/fastq_screen/FastQ_Screen_Genomes/Vectors/Vectors.{suffix}", suffix=["1.bt2","2.bt2","3.bt2","4.bt2","rev.1.bt2","rev.2.bt2"]),
                 expand("analysis/fastq_screen/FastQ_Screen_Genomes/Worm/Caenorhabditis_elegans.WBcel235.{suffix}", suffix=["1.bt2","2.bt2","3.bt2","4.bt2","rev.1.bt2","rev.2.bt2"]),
                 expand("analysis/fastq_screen/FastQ_Screen_Genomes/Yeast/Saccharomyces_cerevisiae.R64-1-1.{suffix}", suffix=["1.bt2","2.bt2","3.bt2","4.bt2","rev.1.bt2","rev.2.bt2"]),
                 expand("analysis/fastq_screen/FastQ_Screen_Genomes/rRNA/GRCm38_rRNA.{suffix}", suffix=["1.bt2","2.bt2","3.bt2","4.bt2","rev.1.bt2","rev.2.bt2"]),
    params:
        outdir = "analysis/fastq_screen/",

    log:
                 "logs/fastq_screen_get_genomes.log"
    resources:
        threads = 8,
        nodes =   1,
        mem_gb =  64,
    conda:       "envs/fastq_screen.yaml"
    shell:
        """
        fastq_screen --get_genomes --outdir {params.outdir} 2> {log}
        rm analysis/fastq_screen/FastQ_Screen_Genomes/fastq_screen.conf
        """

rule download_SILVA:
    # require fastq_screen_get_genomes to run first, creating the FastQ_Screen_Genomes/ where silva is stored
    input:       expand("analysis/fastq_screen/FastQ_Screen_Genomes/Arabidopsis/Arabidopsis_thaliana.TAIR10.{suffix}", suffix=["1.bt2","2.bt2","3.bt2","4.bt2","rev.1.bt2","rev.2.bt2"]),
                 expand("analysis/fastq_screen/FastQ_Screen_Genomes/Adapters/Contaminants.{suffix}", suffix=["1.bt2","2.bt2","3.bt2","4.bt2","rev.1.bt2","rev.2.bt2"]),
                 expand("analysis/fastq_screen/FastQ_Screen_Genomes/Drosophila/BDGP6.{suffix}", suffix=["1.bt2","2.bt2","3.bt2","4.bt2","rev.1.bt2","rev.2.bt2"]),
                 expand("analysis/fastq_screen/FastQ_Screen_Genomes/E_coli/Ecoli.{suffix}", suffix=["1.bt2","2.bt2","3.bt2","4.bt2","rev.1.bt2","rev.2.bt2"]),
                 expand("analysis/fastq_screen/FastQ_Screen_Genomes/Human/Homo_sapiens.GRCh38.{suffix}", suffix=["1.bt2","2.bt2","3.bt2","4.bt2","rev.1.bt2","rev.2.bt2"]),
                 expand("analysis/fastq_screen/FastQ_Screen_Genomes/Lambda/Lambda.{suffix}", suffix=["1.bt2","2.bt2","3.bt2","4.bt2","rev.1.bt2","rev.2.bt2"]),
                 expand("analysis/fastq_screen/FastQ_Screen_Genomes/Lambda/Vectors.{suffix}", suffix=["1.bt2","2.bt2","3.bt2","4.bt2","rev.1.bt2","rev.2.bt2"]),
                 expand("analysis/fastq_screen/FastQ_Screen_Genomes/Mitochondria/mitochondria.{suffix}", suffix=["1.bt2","2.bt2","3.bt2","4.bt2","rev.1.bt2","rev.2.bt2"]),
                 expand("analysis/fastq_screen/FastQ_Screen_Genomes/Mouse/Mus_musculus.GRCm38.{suffix}", suffix=["1.bt2","2.bt2","3.bt2","4.bt2","rev.1.bt2","rev.2.bt2"]),
                 expand("analysis/fastq_screen/FastQ_Screen_Genomes/PhiX/phi_plus_SNPs.{suffix}", suffix=["1.bt2","2.bt2","3.bt2","4.bt2","rev.1.bt2","rev.2.bt2"]),
                 expand("analysis/fastq_screen/FastQ_Screen_Genomes/Rat/Rnor_6.0.{suffix}", suffix=["1.bt2","2.bt2","3.bt2","4.bt2","rev.1.bt2","rev.2.bt2"]),
                 expand("analysis/fastq_screen/FastQ_Screen_Genomes/Vectors/Vectors.{suffix}", suffix=["1.bt2","2.bt2","3.bt2","4.bt2","rev.1.bt2","rev.2.bt2"]),
                 expand("analysis/fastq_screen/FastQ_Screen_Genomes/Worm/Caenorhabditis_elegans.WBcel235.{suffix}", suffix=["1.bt2","2.bt2","3.bt2","4.bt2","rev.1.bt2","rev.2.bt2"]),
                 expand("analysis/fastq_screen/FastQ_Screen_Genomes/Yeast/Saccharomyces_cerevisiae.R64-1-1.{suffix}", suffix=["1.bt2","2.bt2","3.bt2","4.bt2","rev.1.bt2","rev.2.bt2"]),
                 expand("analysis/fastq_screen/FastQ_Screen_Genomes/rRNA/GRCm38_rRNA.{suffix}", suffix=["1.bt2","2.bt2","3.bt2","4.bt2","rev.1.bt2","rev.2.bt2"]),
    output:
        LSU = "analysis/fastq_screen/FastQ_Screen_Genomes/SILVA/SILVA_138.1_LSURef_NR99_tax_silva.fasta",
        SSU = "analysis/fastq_screen/FastQ_Screen_Genomes/SILVA/SILVA_138.1_SSURef_NR99_tax_silva.fasta",
    params:
        SILVA_release = "138.1",
        dir = "analysis/fastq_screen/FastQ_Screen_Genomes/SILVA",
    log:
            "logs/download_SILVA.log"
    resources:
        threads = 8,
        nodes =   1,
        mem_gb =  64,
    conda:       "envs/fastq_screen.yaml"
    shell:
        """
        wget -P {params.dir} https://www.arb-silva.de/fileadmin/silva_databases/release_{params.SILVA_release}/Exports/SILVA_138.1_LSURef_NR99_tax_silva.fasta.gz
        wget -P {params.dir} https://www.arb-silva.de/fileadmin/silva_databases/release_{params.SILVA_release}/Exports/SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz
        gunzip {output.LSU}.gz
        gunzip {output.SSU}.gz
        """

rule fastq_screen_SILVA_db:
    input:
        LSU = "analysis/fastq_screen/FastQ_Screen_Genomes/SILVA/SILVA_138.1_LSURef_NR99_tax_silva.fasta",
        SSU = "analysis/fastq_screen/FastQ_Screen_Genomes/SILVA/SILVA_138.1_SSURef_NR99_tax_silva.fasta",
    params:
        prefix = "analysis/fastq_screen/FastQ_Screen_Genomes/SILVA/SILVA",
    output:   expand("analysis/fastq_screen/FastQ_Screen_Genomes/SILVA/SILVA.{suffix}", suffix=["1.bt2","2.bt2","3.bt2","4.bt2","rev.1.bt2","rev.2.bt2"]),
    log:      "logs/fastq_screen_SILVA_db.log",
    resources:
        threads = 8,
        nodes =   1,
        mem_gb =  64,
    conda:       "envs/fastq_screen.yaml"
    shell:
        """
        bowtie2-build --threads {resources.threads} {input.LSU},{input.SSU} {params.prefix} 2> {log}
        """

rule fastq_screen_sunflower_db:
    input:
        fasta = "references/HanXRQr2.0-SUNRISE-2.1.genome.fasta"
    params:
        prefix = "analysis/fastq_screen/FastQ_Screen_Genomes/Sunflower/HanXRQr2.0-SUNRISE-2.1",
    output:      expand("analysis/fastq_screen/FastQ_Screen_Genomes/Sunflower/HanXRQr2.0-SUNRISE-2.1.{suffix}", suffix=["1.bt2","2.bt2","3.bt2","4.bt2","rev.1.bt2","rev.2.bt2"]),
    log:      "logs/fastq_screen_sunflower_db.log",
    resources:
        threads = 8,
        nodes =   1,
        mem_gb =  64,
    conda:       "envs/fastq_screen.yaml"
    shell:
        """
        bowtie2-build --threads {resources.threads} {input.fasta} {params.prefix} 2> {log}
        """

rule fastq_screen:
    input:
        sunflower =   expand("analysis/fastq_screen/FastQ_Screen_Genomes/Sunflower/HanXRQr2.0-SUNRISE-2.1.{suffix}", suffix=["1.bt2","2.bt2","3.bt2","4.bt2","rev.1.bt2","rev.2.bt2"]),
        silva =       expand("analysis/fastq_screen/FastQ_Screen_Genomes/SILVA/SILVA.{suffix}", suffix=["1.bt2","2.bt2","3.bt2","4.bt2","rev.1.bt2","rev.2.bt2"]),
        arabidopsis = expand("analysis/fastq_screen/FastQ_Screen_Genomes/Arabidopsis/Arabidopsis_thaliana.TAIR10.{suffix}", suffix=["1.bt2","2.bt2","3.bt2","4.bt2","rev.1.bt2","rev.2.bt2"]),
        adapters =    expand("analysis/fastq_screen/FastQ_Screen_Genomes/Adapters/Contaminants.{suffix}", suffix=["1.bt2","2.bt2","3.bt2","4.bt2","rev.1.bt2","rev.2.bt2"]),
        fly =         expand("analysis/fastq_screen/FastQ_Screen_Genomes/Drosophila/BDGP6.{suffix}", suffix=["1.bt2","2.bt2","3.bt2","4.bt2","rev.1.bt2","rev.2.bt2"]),
        e_coli =      expand("analysis/fastq_screen/FastQ_Screen_Genomes/E_coli/Ecoli.{suffix}", suffix=["1.bt2","2.bt2","3.bt2","4.bt2","rev.1.bt2","rev.2.bt2"]),
        human =       expand("analysis/fastq_screen/FastQ_Screen_Genomes/Human/Homo_sapiens.GRCh38.{suffix}", suffix=["1.bt2","2.bt2","3.bt2","4.bt2","rev.1.bt2","rev.2.bt2"]),
        lam =         expand("analysis/fastq_screen/FastQ_Screen_Genomes/Lambda/Lambda.{suffix}", suffix=["1.bt2","2.bt2","3.bt2","4.bt2","rev.1.bt2","rev.2.bt2"]),
        lam_vec =     expand("analysis/fastq_screen/FastQ_Screen_Genomes/Lambda/Vectors.{suffix}", suffix=["1.bt2","2.bt2","3.bt2","4.bt2","rev.1.bt2","rev.2.bt2"]),
        mito =        expand("analysis/fastq_screen/FastQ_Screen_Genomes/Mitochondria/mitochondria.{suffix}", suffix=["1.bt2","2.bt2","3.bt2","4.bt2","rev.1.bt2","rev.2.bt2"]),
        mouse =       expand("analysis/fastq_screen/FastQ_Screen_Genomes/Mouse/Mus_musculus.GRCm38.{suffix}", suffix=["1.bt2","2.bt2","3.bt2","4.bt2","rev.1.bt2","rev.2.bt2"]),
        phiX =        expand("analysis/fastq_screen/FastQ_Screen_Genomes/PhiX/phi_plus_SNPs.{suffix}", suffix=["1.bt2","2.bt2","3.bt2","4.bt2","rev.1.bt2","rev.2.bt2"]),
        rat =         expand("analysis/fastq_screen/FastQ_Screen_Genomes/Rat/Rnor_6.0.{suffix}", suffix=["1.bt2","2.bt2","3.bt2","4.bt2","rev.1.bt2","rev.2.bt2"]),
        vector =      expand("analysis/fastq_screen/FastQ_Screen_Genomes/Vectors/Vectors.{suffix}", suffix=["1.bt2","2.bt2","3.bt2","4.bt2","rev.1.bt2","rev.2.bt2"]),
        worm =        expand("analysis/fastq_screen/FastQ_Screen_Genomes/Worm/Caenorhabditis_elegans.WBcel235.{suffix}", suffix=["1.bt2","2.bt2","3.bt2","4.bt2","rev.1.bt2","rev.2.bt2"]),
        yeast =       expand("analysis/fastq_screen/FastQ_Screen_Genomes/Yeast/Saccharomyces_cerevisiae.R64-1-1.{suffix}", suffix=["1.bt2","2.bt2","3.bt2","4.bt2","rev.1.bt2","rev.2.bt2"]),
        rRNA =        expand("analysis/fastq_screen/FastQ_Screen_Genomes/rRNA/GRCm38_rRNA.{suffix}", suffix=["1.bt2","2.bt2","3.bt2","4.bt2","rev.1.bt2","rev.2.bt2"]),
        #conf =               "analysis/fastq_screen/FastQ_Screen_Genomes/fastq_screen.conf",
        # reads
        R1 =      "raw_data/{sample}-R1.fastq.gz",
        R2 =      "raw_data/{sample}-R2.fastq.gz",
    output:
        R1_html = "analysis/fastq_screen/{sample}-R1_screen.html",
        R1_txt =  "analysis/fastq_screen/{sample}-R1_screen.txt",
        R2_html = "analysis/fastq_screen/{sample}-R2_screen.html",
        R2_txt =  "analysis/fastq_screen/{sample}-R2_screen.txt",
    log:
        R1 =      "logs/fastq_screen/fastq_screen.{sample}-R1.log",
        R2 =      "logs/fastq_screen/fastq_screen.{sample}-R2.log",
    resources:
        threads = 8,
        nodes =   1,
        mem_gb =  64,
    conda:   "envs/fastq_screen.yaml"
    shell:
        """
        fastq_screen \
        --conf bin/fastq_screen.conf \
        --outdir analysis/fastq_screen/ \
        {input.R1} 2> {log.R1}

        fastq_screen \
        --conf bin/fastq_screen.conf \
        --outdir analysis/fastq_screen/ \
        {input.R2} 2> {log.R2}
        """

def trim_galore_input(wildcards):
    if config["PE_or_SE"] == "SE":
        reads = "raw_data/{sample}-SE.fastq.gz".format(**wildcards)
        return reads
    elif config["PE_or_SE"] == "PE":
        R1 = "raw_data/{sample}-R1.fastq.gz".format(**wildcards)
        R2 = "raw_data/{sample}-R2.fastq.gz".format(**wildcards)
        return [R1,R2]

rule trim_galore_PE:
    input:
        trim_galore_input
    output:
        "analysis/trimmed_data/{sample}-R1_val_1.fq.gz",
        "analysis/trimmed_data/{sample}-R1_val_1_fastqc.html",
        "analysis/trimmed_data/{sample}-R1_val_1_fastqc.zip",
        "analysis/trimmed_data/{sample}-R1.fastq.gz_trimming_report.txt",
        "analysis/trimmed_data/{sample}-R2_val_2.fq.gz",
        "analysis/trimmed_data/{sample}-R2_val_2_fastqc.html",
        "analysis/trimmed_data/{sample}-R2_val_2_fastqc.zip",
        "analysis/trimmed_data/{sample}-R2.fastq.gz_trimming_report.txt"
    params:
        outdir="analysis/trimmed_data/"
    log:
        stdout="logs/trim_galore/{sample}.o",
        stderr="logs/trim_galore/{sample}.e"
    benchmark:
        "benchmarks/trim_galore/{sample}.txt"
    conda:
        "envs/trim-galore.yaml"
    resources:
        threads = 4,
        nodes =   1,
        mem_gb =  64,
    shell:
        """
        trim_galore \
        --paired \
        {input} \
        --output_dir {params.outdir} \
        --cores {threads} \
        -q 20 \
        --fastqc \
        2> {log.stderr} 1> {log.stdout}
        """

def STAR_input(wildcards):
    if config["PE_or_SE"] == "SE":
        fq1="analysis/trimmed_data/{sample}-SE_trimmed.fq.gz".format(**wildcards)
        return fq1
    elif config["PE_or_SE"] == "PE":
        fq1 = "analysis/trimmed_data/{sample}-R1_val_1.fq.gz".format(**wildcards)
        fq2 = "analysis/trimmed_data/{sample}-R2_val_2.fq.gz".format(**wildcards)
        return [fq1,fq2]

def STAR_params(wildcards):
    index = config["ref"]["index"][units.loc["{sample}".format(**wildcards)]["accession"]]
    outprefix = "analysis/star/{sample}."
    return [index,outprefix]
    outprefix = "analysis/star/{sample}.".format(**wildcards)
    return [index, outprefix]

def STAR_index(wildcards):
    samp = "{sample}".format(**wildcards)
    return config["ref"]["index"][units.loc[samp]["accession"]]

def STAR_outprefix(wildcards):
    prefix = "analysis/star/{sample}.".format(**wildcards)
    return prefix

rule star:
    input:
        STAR_input
    output:
        # see STAR manual for additional output files
        bam =                 "analysis/star/{sample}.Aligned.sortedByCoord.out.bam",
        bai =                 "analysis/star/{sample}.Aligned.sortedByCoord.out.bam.bai",
        log_final =           "analysis/star/{sample}.Log.final.out",
        log =                 "analysis/star/{sample}.Log.out",
        rpg =                 "analysis/star/{sample}.ReadsPerGene.out.tab",
        sj =                  "analysis/star/{sample}.SJ.out.tab",
        g_dir =     directory("analysis/star/{sample}._STARgenome"),
        pass1_dir = directory("analysis/star/{sample}._STARpass1"),
    params:
        index = lambda wildcards: config["ref"]["index"][units.loc["{sample}".format(**wildcards)]["status"]],
        outprefix = "analysis/star/{sample}."
	# STAR_outprefix
    log:
        "logs/star/{sample}.log"
    benchmark:
        "benchmarks/star/{sample}.txt"
    conda:
        "envs/star.yaml"
    resources:
        threads = 8,
        nodes =   1,
        mem_gb =  64,
    shell:
        """
        STAR \
        --runThreadN {threads} \
        --genomeDir {params.index} \
        --readFilesIn {input} \
        --twopassMode Basic \
        --readFilesCommand zcat \
        --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix {params.outprefix} \
        --quantMode GeneCounts \
        --outStd Log 2> {log}

        samtools index -b -@ {resources.threads} {output.bam}
        """

rule star_dedup:
    input:
        "analysis/star/{sample}.Aligned.sortedByCoord.out.bam",
    output:
        bam = "analysis/star_dedup/{sample}.Aligned.sortedByCoord.out.dedup.bam",
        bai = "analysis/star_dedup/{sample}.Aligned.sortedByCoord.out.dedup.bam.bai",
    log:
        "logs/star_dedup/{sample}.log"
    conda:
        "envs/samblaster.yaml"
    resources:
        threads = 8,
        nodes =   1,
        mem_gb =  64
    shell:
        """
        samtools sort -@ {resources.threads} -n {input} | \
        samtools view -h | \
        samblaster --removeDups  2> {log} | \
        samtools view -Sb - | samtools sort -@ {resources.threads} - > {output.bam}
        samtools index -b -@ {resources.threads} {output.bam}
        """

rule dedup_2_fastq:
    input:         "analysis/star_dedup/{sample}.Aligned.sortedByCoord.out.dedup.bam",
    output:
        unpaired = "analysis/star_dedup/{sample}-unpaired.fastq.gz",
        R1 =       "analysis/star_dedup/{sample}-R1.fastq.gz",
        R2 =       "analysis/star_dedup/{sample}-R2.fastq.gz",
    log:
        "logs/dedup_2_fastq/{sample}.log"
    conda:
        "envs/bamutil.yaml"
    resources:
        threads = 8,
        nodes =   1,
        mem_gb =  64,
    shell:
        """
        bam bam2FastQ --in {input} --gzip \
        --unpairedOut {output.unpaired} \
        --firstOut {output.R1} \
        --secondOut {output.R2} \
        2> {log}
        """

multiqc_input = []
if config["PE_or_SE"] =="SE":
    # multiqc_input.append(expand("analysis/fastq_screen/{units.sample}-SE_screen.txt", units=units.itertuples()))
    multiqc_input.append(expand("analysis/trimmed_data/{units.sample}-SE_trimmed.fq.gz", units=units.itertuples()))
    multiqc_input.append(expand("analysis/trimmed_data/{units.sample}-SE_trimmed_fastqc.zip", units=units.itertuples()))
    multiqc_input.append(expand("analysis/trimmed_data/{units.sample}-SE_trimmed_fastqc.html", units=units.itertuples()))
    multiqc_input.append(expand("analysis/star/{units.sample}.Log.final.out", units=units.itertuples()))
elif config["PE_or_SE"] =="PE":
    # multiqc_input.append(expand("analysis/fastq_screen/{units.sample}-R{read}_screen.txt", units=units.itertuples(), read=["1","2"]))
    multiqc_input.append(expand("analysis/trimmed_data/{units.sample}-R{read}_val_{read}.fq.gz", units=units.itertuples(), read=["1","2"]))
    multiqc_input.append(expand("analysis/trimmed_data/{units.sample}-R{read}_val_{read}_fastqc.html", units=units.itertuples(), read=["1","2"]))
    multiqc_input.append(expand("analysis/trimmed_data/{units.sample}-R{read}_val_{read}_fastqc.zip", units=units.itertuples(), read=["1","2"]))
    multiqc_input.append(expand("analysis/trimmed_data/{units.sample}-R{read}.fastq.gz_trimming_report.txt", units=units.itertuples(), read=["1","2"]))
    multiqc_input.append(expand("analysis/star/{units.sample}.Log.final.out", units=units.itertuples()))

rule multiqc:
    input:
        multiqc_input
    params:
        "analysis/star/",
        "analysis/trimmed_data/",
        # "analysis/fastq_screen/",
    output:
        "analysis/multiqc/multiqc_report.html",
        "analysis/multiqc/multiqc_report_data/multiqc.log",
        "analysis/multiqc/multiqc_report_data/multiqc_cutadapt.txt",
        "analysis/multiqc/multiqc_report_data/multiqc_fastqc.txt",
        "analysis/multiqc/multiqc_report_data/multiqc_general_stats.txt",
        "analysis/multiqc/multiqc_report_data/multiqc_sources.txt",
        "analysis/multiqc/multiqc_report_data/multiqc_star.txt",
    log:
        "logs/multiqc.log"
    benchmark:
        "benchmarks/multiqc/multiqc.txt"
    resources:
        threads = 1,
        nodes = 1,
        mem_gb = 64,
    conda:
        "envs/multiqc.yaml"
    shell:
        """
        multiqc -f {params} \
        -o analysis/multiqc \
        --ignore '*._STARpass1/*' \
        -n multiqc_report.html \
        --cl-config 'max_table_rows: 999999' \
        2> {log}
        """

rule edgeR:
    input:
        expand("analysis/star/{units.sample}.Aligned.sortedByCoord.out.bam", units=units.itertuples()),
        "analysis/multiqc/multiqc_report.html" # require multiQC to be run before this analysis
    output:
        "bin/diffExp.html"
    log:
        "logs/edgeR.log"
    benchmark:
        "benchmarks/edgeR/edgeR.txt"
    conda:
        "envs/R.yaml"
    resources:
        threads = 1,
        nodes = 1,
        mem_gb = 16,
    script:
        "bin/diffExp.Rmd"


rule addRG:
    input:
        "analysis/star/{sample}.Aligned.sortedByCoord.out.bam"
    output:
        bam="analysis/00_addRG/{sample}.Aligned.sortedByCoord.out.addRG.bam",
    log:
        out="logs/00_addRG/{sample}.o",
        err="logs/00_addRG/{sample}.e"
    benchmark:
        "benchmarks/00_addRG/{sample}.txt"
    envmodules:
        "picard/2.9.0"
    params:
        unit=lambda wildcards: var_calling_units[var_calling_units["unit"]==wildcards.sample]['unit'].values[0],
        lib=lambda wildcards: var_calling_units[var_calling_units["unit"]==wildcards.sample]['library'].values[0],
        sample=lambda wildcards: var_calling_units[var_calling_units["unit"]==wildcards.sample]['sample'].values[0],
        platform_unit=lambda wildcards: var_calling_units[var_calling_units["unit"]==wildcards.sample]['platform_unit'].values[0]
    resources:
        threads = 4,
        nodes =   1,
        mem_gb =  64
    shell:
        """
        java -Xms16g -Xmx64g -Djava.io.tmpdir=./tmp -jar $PICARD AddOrReplaceReadGroups \
        --INPUT {input} \
        --OUTPUT {output.bam} \
        --RGID {params.unit} \
        --RGLB {params.lib} \
        --RGPL ILLUMINA \
        --RGPU {params.platform_unit} \
        --RGSM {params.sample} 1>{log.out} 2>{log.err}
        """

#
# # variant calling based on the GATK best practices as documented at https://github.com/gatk-workflows/gatk4-rnaseq-germline-snps-indels/blob/master/gatk4-rna-best-practices.wdl -- accessed Aug 11, 2020
#
# rule markdups:
#     input:
#         lambda wildcards: expand("analysis/00_addRG/{sample}.Aligned.sortedByCoord.out.addRG.bam", sample=var_calling_units[var_calling_units["sample"]==wildcards.sample].index.values)
#     output:
#         bam="analysis/01_markdup/{sample}.Aligned.sortedByCoord.out.addRG.mrkdup.bam",
#         metrics="analysis/01_markdup/{sample}.Aligned.sortedByCoord.out.addRG.mrkdup.metrics"
#     params:
#         input_params=lambda wildcards: expand("--INPUT analysis/00_addRG/{sample}.Aligned.sortedByCoord.out.addRG.bam", sample=var_calling_units[var_calling_units["sample"]==wildcards.sample].index.values)
#     log:
#         out="logs/01_markdup/{sample}.o",
#         err="logs/01_markdup/{sample}.e"
#     benchmark:
#         "benchmarks/01_markdup/{sample}.txt"
#     envmodules:
#         "bbc/picard/picard-2.23.3"
#     resources:
#         threads = 4,
#         nodes =   1,
#         mem_gb =  64,
#     shell:
#         """
#         java -Xms16g -Xmx{resources.mem_gb}g -Djava.io.tmpdir=./tmp -jar $PICARD MarkDuplicates \
#         {params.input_params} \
#         --OUTPUT {output.bam} \
#         --CREATE_INDEX true \
#         --VALIDATION_STRINGENCY SILENT \
#         --METRICS_FILE {output.metrics} 1>{log.out} 2>{log.err}
#         """
#
# rule splitncigar:
#     input:
#         "analysis/01_markdup/{sample}.Aligned.sortedByCoord.out.addRG.mrkdup.bam"
#     output:
#         "analysis/02_splitncigar/{sample}.Aligned.sortedByCoord.out.addRG.mrkdup.splitncigar.bam"
#     log:
#         out="logs/02_splitncigar/{sample}.o",
#         err="logs/02_splitncigar/{sample}.e"
#     benchmark:
#         "benchmarks/02_splitncigar/{sample}.txt"
#     params:
#         ref_fasta=config["ref"]["sequence"],
#     envmodules:
#         "bbc/gatk/gatk-4.1.8.1"
#     resources:
#         threads = 4,
#         nodes =   1,
#         mem_gb =  64,
#     shell:
#         """
#         gatk --java-options "-Xms8g -Xmx{resources.mem_gb}g -Djava.io.tmpdir=./tmp" \
#         SplitNCigarReads \
#         -R {params.ref_fasta} \
#         -I {input} \
#         -O {output} 1>{log.out} 2>{log.err}
#
#         """
#
# rule base_recalibrate:
#     input:
#         "analysis/02_splitncigar/{sample}.Aligned.sortedByCoord.out.addRG.mrkdup.splitncigar.bam"
#     output:
#         "analysis/03_base_recal/{sample}.Aligned.sortedByCoord.out.addRG.mrkdup.splitncigar.bam.recal_data.table"
#     log:
#         out="logs/03_base_recal/{sample}.o",
#         err="logs/03_base_recal/{sample}.e"
#     benchmark:
#         "benchmarks/03_base_recal/{sample}.txt"
#     params:
#         ref_fasta=config["ref"]["sequence"],
#         known_snps=config["ref"]["known_snps"],
#         known_indels=config["ref"]["known_indels"]
#     envmodules:
#         "bbc/gatk/gatk-4.1.8.1"
#     resources:
#         threads = 4,
#         nodes =   1,
#         mem_gb =  64,
#     shell:
#         """
#         gatk --java-options "-Xms8g -Xmx{resources.mem_gb}g -XX:+UseParallelGC -XX:ParallelGCThreads={threads} -Djava.io.tmpdir=./tmp" \
#             BaseRecalibrator \
#             -R {params.ref_fasta} \
#             -I {input} \
#             -O {output} \
#             -known-sites {params.known_snps} \
#             -known-sites {params.known_indels} 1>{log.out} 2>{log.err}
#
#         """
#
# rule applyBQSR:
#     input:
#         bam="analysis/02_splitncigar/{sample}.Aligned.sortedByCoord.out.addRG.mrkdup.splitncigar.bam",
#         recal_table="analysis/03_base_recal/{sample}.Aligned.sortedByCoord.out.addRG.mrkdup.splitncigar.bam.recal_data.table"
#     output:
#         "analysis/04_apply_base_recal/{sample}.Aligned.sortedByCoord.out.addRG.mrkdup.splitncigar.baserecal.bam"
#     log:
#         out="logs/04_apply_base_recal/{sample}.o",
#         err="logs/04_apply_base_recal/{sample}.e"
#     benchmark:
#         "benchmarks/04_apply_base_recal/{sample}.txt"
#     params:
#         ref_fasta=config["ref"]["sequence"],
#     envmodules:
#         "bbc/gatk/gatk-4.1.8.1"
#     resources:
#         threads = 4,
#         nodes =   1,
#         mem_gb =  64,
#     shell:
#         """
#         gatk --java-options "-Xms8g -Xmx{resources.mem_gb}g  -XX:+UseParallelGC -XX:ParallelGCThreads={threads} -Djava.io.tmpdir=./tmp" \
#             ApplyBQSR \
#             --add-output-sam-program-record \
#             -R {params.ref_fasta} \
#             -I {input.bam} \
#             -O {output} \
#             --bqsr-recal-file {input.recal_table} 1>{log.out} 2>{log.err}
#
#         """
#
# rule haplotypecaller:
#     input:
#         bam="analysis/04_apply_base_recal/{sample}.Aligned.sortedByCoord.out.addRG.mrkdup.splitncigar.baserecal.bam"
#     output:
#         "analysis/05_haplotypecaller/{sample}.{contig_group}.Aligned.sortedByCoord.out.addRG.mrkdup.splitncigar.baserecal.g.vcf.gz"
#     log:
#         out="logs/05_haplotypecaller/{sample}.{contig_group}.o",
#         err="logs/05_haplotypecaller/{sample}.{contig_group}.e"
#     benchmark:
#         "benchmarks/05_haplotypecaller/{sample}.{contig_group}.txt"
#     params:
#         dbsnp=config["ref"]["known_snps"],
#         ref_fasta=config["ref"]["sequence"],
#         contigs = lambda wildcards: "-L " + contig_groups[contig_groups.name == wildcards.contig_group]['contigs'].values[0].replace(",", " -L "),
#
#     envmodules:
#         "bbc/gatk/gatk-4.1.8.1"
#     resources:
#         threads = 4,
#         nodes =   1,
#         mem_gb =  64,
#     shell:
#         """
#         gatk --java-options "-Xms8g -Xmx{resources.mem_gb}g -Djava.io.tmpdir=./tmp" \
#         HaplotypeCaller \
#         -R {params.ref_fasta} \
#         -I {input.bam} \
#         -O {output} \
#         -ERC GVCF \
#         -dont-use-soft-clipped-bases \
#         --native-pair-hmm-threads {threads} \
#         --standard-min-confidence-threshold-for-calling 20 \
#         --dbsnp {params.dbsnp} \
#         {params.contigs} 1>{log.out} 2>{log.err}
#
#         """
#
# rule combinevar:
#     input:
#         lambda wildcards: expand("analysis/05_haplotypecaller/{sample}.{contig_group}.Aligned.sortedByCoord.out.addRG.mrkdup.splitncigar.baserecal.g.vcf.gz", sample=var_calling_units['sample'].unique(), contig_group=wildcards.contig_group)
#
#     output:
#         touch=touch("analysis/06_combinevar/{contig_group}.done"),
#         genomicsdb=directory("analysis/06_combinevar/{contig_group}.genomicsdb"),
#     log:
#         out="logs/06_combinevar/all.{contig_group}.o",
#         err="logs/06_combinevar/all.{contig_group}.e"
#     benchmark:
#         "benchmarks/06_combinevar/{contig_group}.txt"
#     params:
#         sample_gvcfs = lambda wildcards, input: list(map("-V {}".format, input)),
#         contigs = lambda wildcards: "-L " + contig_groups[contig_groups.name == wildcards.contig_group]['contigs'].values[0].replace(",", " -L "),
#     envmodules:
#         "bbc/gatk/gatk-4.1.8.1"
#     resources:
#         threads = 4,
#         nodes =   1,
#         mem_gb =  64,
#     shell:
#         """
#         gatk --java-options "-Xms8g -Xmx{resources.mem_gb}g -Djava.io.tmpdir=./tmp" \
#         GenomicsDBImport \
#         {params.sample_gvcfs} \
#         --genomicsdb-workspace-path {output.genomicsdb} \
#         {params.contigs} 1>{log.out} 2>{log.err}
#
#         """
#
# rule jointgeno:
#     input:
#         "analysis/06_combinevar/{contig_group}.done"
#         #"analysis/06_combinevar/{contig_group}.genomicsdb"
#     output:
#         vcf="analysis/07_jointgeno/all.{contig_group}.vcf.gz",
#     log:
#         out="logs/07_jointgeno/all.{contig_group}.o",
#         err="logs/07_jointgeno/all.{contig_group}.e"
#     benchmark:
#         "benchmarks/07_jointgeno/all.{contig_group}.txt"
#     params:
#         ref_fasta=config["ref"]["sequence"],
#         genomicsdb="analysis/06_combinevar/{contig_group}.genomicsdb"
#     envmodules:
#         "bbc/gatk/gatk-4.1.8.1"
#     resources:
#         threads = 4,
#         nodes =   1,
#         mem_gb =  64,
#     shell:
#         """
#         gatk --java-options "-Xms8g -Xmx{resources.mem_gb}g -Djava.io.tmpdir=./tmp" \
#         GenotypeGVCFs \
#         -R {params.ref_fasta} \
#         -V gendb://{params.genomicsdb} \
#         -O {output.vcf} 1>>{log.out} 2>>{log.err}
#         """
#
# rule merge_and_filter_vcf:
#     input:
#         expand("analysis/07_jointgeno/all.{contig_grp}.vcf.gz", contig_grp=contig_groups.name)
#     output:
#         raw="analysis/08_merge_and_filter/all.merged.vcf.gz",
#         filt="analysis/08_merge_and_filter/all.merged.filt.vcf.gz",
#         pass_only="analysis/08_merge_and_filter/all.merged.filt.PASS.vcf.gz",
#         vt_peek_raw="analysis/08_merge_and_filter/all.merged.vcf.gz.vt_peek.txt",
#         vt_peek_pass="analysis/08_merge_and_filter/all.merged.filt.PASS.vcf.gz.vt_peek.txt"
#     log:
#         out="logs/08_merge_and_filter/out.o",
#         err="logs/08_merge_and_filter/out.e"
#     benchmark:
#         "benchmarks/08_merge_and_filter/benchmark.txt"
#     params:
#         ref_fasta=config["ref"]["sequence"],
#         in_vcfs=expand("--INPUT analysis/07_jointgeno/all.{contig_grp}.vcf.gz", contig_grp=contig_groups.name)
#     envmodules:
#         "bbc/gatk/gatk-4.1.8.1",
#         "bbc/vt/vt-0.1.16"
#     resources:
#         threads = 4,
#         nodes =   1,
#         mem_gb =  64,
#     shell:
#         """
#         gatk --java-options "-Xms8g -Xmx{resources.mem_gb}g -Djava.io.tmpdir=./tmp" \
#         MergeVcfs \
#         {params.in_vcfs} \
#         --OUTPUT {output.raw} \
#         1>>{log.out} 2>>{log.err}
#
#         vt peek -r {params.ref_fasta} {output.raw} 2> {output.vt_peek_raw} 1>>{log.out}
#
#         echo "mergeVcfs done." >> {log.out}
#         echo "mergeVcfs done." >> {log.err}
#
#         gatk --java-options "-Xms8g -Xmx{resources.mem_gb}g -Djava.io.tmpdir=./tmp" \
#         VariantFiltration \
#         --R {params.ref_fasta} \
#         --V {output.raw} \
#         --window 35 \
#         --cluster 3 \
#         --filter-name "FS" \
#         --filter "FS > 30.0" \
#         --filter-name "QD" \
#         --filter "QD < 2.0" \
#         -O {output.filt} \
#         1>>{log.out} 2>>{log.err}
#
#         echo "VariantFiltration done." >> {log.out}
#         echo "VariantFiltration done." >> {log.err}
#
#         gatk --java-options "-Xms8g -Xmx{resources.mem_gb}g -Djava.io.tmpdir=./tmp" \
#         SelectVariants \
#         -R {params.ref_fasta} \
#         -V {output.filt} \
#         --exclude-filtered \
#         -O {output.pass_only} \
#         1>>{log.out} 2>>{log.err}
#
#         echo "SelectVariants done." >> {log.out}
#         echo "SelectVariants done." >> {log.err}
#
#         vt peek -r {params.ref_fasta} {output.pass_only} 2> {output.vt_peek_pass} 1>>{log.out}
#         """
#
# rule variant_annot:
#     input:
#         "analysis/08_merge_and_filter/all.merged.filt.PASS.vcf.gz"
#     output:
#         html="analysis/09a_variant_annot/all.merged.filt.PASS.snpeff.html",
#         vcf="analysis/09a_variant_annot/all.merged.filt.PASS.snpeff.vcf.gz",
#         tbi="analysis/09a_variant_annot/all.merged.filt.PASS.snpeff.vcf.gz.tbi",
#         html_canon="analysis/09a_variant_annot/all.merged.filt.PASS.snpeff_canonical.html",
#         vcf_canon="analysis/09a_variant_annot/all.merged.filt.PASS.snpeff_canonical.vcf.gz",
#         tbi_canon="analysis/09a_variant_annot/all.merged.filt.PASS.snpeff_canonical.vcf.gz.tbi",
#     log:
#         out="logs/09a_variant_annot/out.o",
#         err="logs/09a_variant_annot/out.e"
#     benchmark:
#         "benchmarks/09a_variant_annot/benchmark.txt"
#     params:
#         db_id=config["ref"]["snpeff_db_id"],
#     envmodules:
#         "bbc/SnpEff/SnpEff-4.3t",
#         "bbc/htslib/htslib-1.10.2"
#     resources:
#         threads = 4,
#         nodes =   1,
#         mem_gb =  64,
#     shell:
#         """
#         java -Xms8g -Xmx{resources.mem_gb}g -Djava.io.tmpdir=./tmp -jar $SNPEFF/snpEff.jar eff \
#         -v \
#         -canon \
#         -onlyProtein \
#         -stats {output.html_canon} \
#         {params.db_id} \
#         {input} \
#         2>>{log.err} | \
#         bgzip > {output.vcf_canon}
#
#         tabix {output.vcf_canon} 2>>{log.err} 1>>{log.out}
#
#
#         java -Xms8g -Xmx{resources.mem_gb}g -Djava.io.tmpdir=./tmp -jar $SNPEFF/snpEff.jar eff \
#         -v \
#         -onlyProtein \
#         -stats {output.html} \
#         {params.db_id} \
#         {input} \
#         2>>{log.err} | \
#         bgzip > {output.vcf}
#
#         tabix {output.vcf} 2>>{log.err} 1>>{log.out}
#         """
#
#
# rule snprelate:
#     input:
#         "analysis/08_merge_and_filter/all.merged.filt.PASS.vcf.gz"
#     output:
#         "analysis/09b_snp_pca_and_dendro/report.html"
#     params:
#         gds="analysis/09b_snp_pca_and_dendro/all.gds"
#     conda:
#         "envs/R.yaml"
#     #envmodules:
#     #    "bbc/cairo/cairo-1.16.0",
#     #    "bbc/R/R-3.6.0",
#     #    "bbc/pandoc/pandoc-2.7.3",
#     resources:
#         threads = 4,
#         nodes =   1,
#         mem_gb =  64,
#     script:
#         "bin/snprelate.Rmd"
