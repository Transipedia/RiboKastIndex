# ---------------------------------------------------------------
#               RiboKastIndex Pipeline

# Created by: Safa MADDOURI 
# <I2BC/SSFA>
# ---------------------------------------------------------------
# Imports packages
from optparse import OptionParser
import gffutils
import re
import os


# Sets the number of threads for multi-threading steps
multi_threads_nbr = 3
mem_mb_resources = 10000
utr_threshold = "0.25"

# ---------------------------------------------------------------
#               Paths from the config file
# ---------------------------------------------------------------
# Read the paths from the config file
local_path = config["paths"]["local_path"]
RibokastIndex_tools = config["paths"]["RibokastIndex_tools"]
results_path = config["paths"]["results_path"]
stats_path = config["paths"]["stats_path"]
logs_path = config["paths"]["logs_path"]
snakemake_log_path = config["paths"]["snakemake_log_path"]
fastq_path = config["paths"]["fastq_path"]

# ---------------------------------------------------------------
#               Sample Wildcards and Length Definitions
# ---------------------------------------------------------------
# Wildcards definition
SAMPLES, = SAMPLES, = glob_wildcards(fastq_path + "{sample}.fastq.gz")
SAMPLES.sort()
LENGTHS = list(map(str,range(int(config['readsLength_min']),int(config['readsLength_max'])+1)))
# Strings with minimum and maximum read lengths to be used in file names
frag_length_S = "." + LENGTHS[0]
frag_length_L = "." + LENGTHS[0] + "-" + LENGTHS[len(LENGTHS)-1]

# ---------------------------------------------------------------
#               Genome Isoform and GFF Attributes
# ---------------------------------------------------------------
# Mean number of isoforms in the genome (for transcriptome multi-alignment parameter)
isoform_nbr = "10"

# Check if attributes specified by the user are present in the gff file
htseq_additional = ""
htseq_header = 'ID\t'
fields = "1,2"

f = open(local_path + "database/" + config['gff'],"r")
for l in f:
    is_name = re.search("^([^\t]+\t){8}.*" + config['gff_name_attribut'],l)
    if is_name:
        name_in_gff = True

        features_of_interest = [config['gff_cds_feature']]
        if config['UTR'] == "yes":
            features_of_interest.extend([config['gff_5UTR_feature'],config['gff_3UTR_feature']])

        htseq_additional = "--additional-attr " + config['gff_name_attribut']
        htseq_header = 'ID\tName\t'
        fields = "1,3"
        break
    else:
        name_in_gff = False

f.close()
# ---------------------------------------------------------------
#               Onstart / Onsuccess / Onerror Instructions
# ---------------------------------------------------------------
# Create logs directory when RiboDoc starts
onstart:
    shell("mkdir -p " + logs_path)

# Actions to be taken upon successful execution
onsuccess:
    shell("cp " + local_path + "config.yaml " + results_path)  # Copy config file

# Actions in case of errors
onerror:
    shell("cp " + local_path + "config.yaml " + results_path)


def get_kmer_tsv_paths(wildcards):
    checkpoint_output = checkpoints.get_samples_for_kmer.get(**wildcards).output[0]
    with open(checkpoint_output) as f:
        samples = [line.strip() for line in f if line.strip()]
    return expand(results_path + "kmerCount/Kmer/{sample}.tsv", sample=samples)

# ---------------------------------------------------------------
#               Main Workflow
# ---------------------------------------------------------------
rule all:
    input:
        #results_path + "kmerCount/Kamrat/merged-res.tsv",
        expand(results_path + "fastqc/fastqc_before_trimming/{sample}_fastqc.html", sample=SAMPLES),
        expand(results_path + "fastqc/fastqc_after_trimming/{sample}.cutadapt" + frag_length_L + "_fastqc.html", sample=SAMPLES),
        expand(results_path + "fastqc/fastqc_after_trimming/{sample}.cutadapt" + frag_length_L + "_fastqc.zip", sample=SAMPLES),
        expand(results_path + "no-outRNA/{sample}" + frag_length_L + ".no-outRNA.fastq.gz",sample=SAMPLES),
        expand(results_path + "riboWaltz" + frag_length_L + "/{sample}/psite_offset.csv", sample=SAMPLES),
        expand(results_path + "riboWaltz" + frag_length_L + "/{sample}/frame_psite_length.csv", sample=SAMPLES),
        expand(results_path + "riboWaltz" + frag_length_L + "/{sample}/frame_psite.csv", sample=SAMPLES),
        expand(results_path + "riboWaltz" + frag_length_L + "/{sample}/best_offset.csv", sample=SAMPLES),
        results_path + "riboWaltz" + frag_length_L + "/frame_psite_length.csv",
        results_path + "riboWaltz" + frag_length_L + "/psite_offset.csv",
        results_path + "riboWaltz" + frag_length_L + "/frame_psite.csv",
        results_path + "riboWaltz" + frag_length_L + "/best_offset.csv",
        results_path + "riboWaltz" + frag_length_L + "/psite_table_forKmerCount.txt",
        results_path + "kmerCount/samples_kmer.txt",
        results_path + "kmerCount/Kamrat/merged-res.tsv"
        
########################################################################################################################
#                                               PIPELINE RULES
########################################################################################################################
# ---------------------------------------------------------------
#               Rule: Find Adapter Sequence
# ---------------------------------------------------------------
# Find the adapter sequence if not set in config file
rule find_adapter_sequence:
    input:
        fastq = fastq_path + "{sample}.fastq.gz"
    output:
        adapter = results_path + "adapter_lists/{sample}.txt"
    log:
        rscript = logs_path + "find_adapter_sequence/{sample}.rscript.log",
        sed = logs_path + "find_adapter_sequence/{sample}.sed.log",
        echo = logs_path + "find_adapter_sequence/{sample}.echo.log",
        touch = logs_path + "find_adapter_sequence/{sample}.touch.log"
    params:
        Result_path = results_path  # Assuming results_path is the desired local path
    shell:
        # Touch the output file to ensure it exists
        "touch {output.adapter} 2> {log.touch} ;"
        
        # If no adapter sequence is provided, run the R script to find the adapter
        "if [ -z " + config['adapt_sequence'] + " ]; then "
        "Rscript " + RibokastIndex_tools + "others/find_adapter_sequencearg.R {input.fastq} {params.Result_path} 2> {log.rscript} ;"
        
        # If an adapter sequence is already known, echo it into the output file
        "elif [ '" + config['already_trimmed'] + "' = 'no' ]; then echo " + config['adapt_sequence'] + " 1> {output.adapter} 2> {log.echo} ;"
        "fi"
# ---------------------------------------------------------------
#               Rule: Name CDS
# ---------------------------------------------------------------
# Adds transcript names and gene IDs to the CDS and exon lines if possible
rule name_CDS:
    input:
        gff = local_path + "database/" + os.path.basename(config['gff'])
    output:
        gff_namedCDS = results_path + "annex_database/NamedCDS_" + os.path.basename(config['gff'])
    run:
        gene_id_bool = True
        if name_in_gff == True:
            db = gffutils.create_db(input.gff, ':memory:', merge_strategy='create_unique', keep_order=True)
            with open(output.gff_namedCDS, 'w') as fout:
                for d in db.directives:
                    fout.write('##{0}\n'.format(d))
                for feature in db.all_features():
                    if feature.featuretype in features_of_interest or feature.featuretype == 'exon':
                        parent = list(db.parents(feature, level=1))
                        if len(parent) > 0:
                            parent = parent[0]
                            if parent.attributes.get(config['gff_name_attribut']) and not feature.attributes.get(config['gff_name_attribut']):
                                feature.attributes[config['gff_name_attribut']] = [i.replace("mRNA","cds") for i in parent.attributes.get(config['gff_name_attribut'])]
                                feature.attributes[config['gff_name_attribut']][0] + "_name"
                            if parent.attributes.get('ID') and not feature.attributes.get('ID'):
                                feature.attributes["ID"] = parent.attributes["ID"]
                                feature.attributes['ID'] = feature.attributes['ID'][0] + "_CDS"
                            if parent.attributes.get('ID') and not parent.attributes.get(config['gff_name_attribut']):
                                feature.attributes[config['gff_name_attribut']] = parent.attributes["ID"]
                                feature.attributes[config['gff_name_attribut']][0] + "_name"
                        if feature.attributes.get(config['gff_name_attribut']):
                            fout.write(str(feature) + ';\n')
                    else:
                        fout.write(str(feature) + '\n')
        else:
            shell("cp {input.gff} {output.gff_namedCDS}")
        if gene_id_bool:
            print("'gene_id' attributes are present.")
        else:
            print("Missing at least some 'gene_id' attribute in this gff.")
        shell("sed -i -E 's/\\s/\\t/8' {output.gff_namedCDS}")
        shell('sed -i -E "s/\\"//g" {output.gff_namedCDS}')

# ---------------------------------------------------------------
#               Rule: FastQC Before Trimming
# ---------------------------------------------------------------
# Quality control of data : build of the fastqc on raw data
rule make_fastqc_before_trimming:
    input:
        fastq_path + "{sample}.fastq.gz"
    output:
        results_path + "fastqc/fastqc_before_trimming/{sample}_fastqc.zip",
        results_path + "fastqc/fastqc_before_trimming/{sample}_fastqc.html"
    log:
        logs_path + "make_fastqc_before_trimming/{sample}.log"
    benchmark:
        local_path + "benchmarks/make_fastqc_before_trimming/{sample}.benchmark.txt"
    params:
       outdir = results_path + "fastqc/fastqc_before_trimming/"
    shell:
        "fastqc {input} --outdir {params.outdir} 2> {log}"

# ---------------------------------------------------------------
#               Rule: Adapter Trimming
# ---------------------------------------------------------------
# Removes/cuts potential adapters on the reads
rule adapt_trimming:
    input:
        fastq = fastq_path + "{sample}.fastq.gz",
        adapt_seq = results_path + "adapter_lists/{sample}.txt"
    output:
        cut_fastq = results_path + "cutadapt/{sample}.cutadapt" + frag_length_L + ".fastq.gz"
    log:
        trim_value = logs_path + "adapt_trimming/{sample}_trim_value.log",
        cutadapt = logs_path + "adapt_trimming/{sample}_cutadapt.log",
        cutadapt_out = stats_path + "{sample}_adapt_trimming.log"
    benchmark:
        local_path + "benchmarks/adapt_trimming/{sample}.benchmark.txt"
    resources:
        mem_mb = mem_mb_resources
    threads:
        multi_threads_nbr
    shell:
        # Read adapter sequence
        "adapter_sequence=`cat {input.adapt_seq}` ;"
        # Check if adapter_sequence is empty or if 'already_trimmed' is set to 'yes'
        "if [ -z \"$adapter_sequence\" ] || [ '" + config['already_trimmed'] + "' = 'yes' ]; then "
        "trim=''; "
        "else "
        "trim=\"-a ${{adapter_sequence}} --trimmed-only\"; "
        "fi 2> {log.trim_value} ;"
        # Run cutadapt with appropriate trimming options
        "cutadapt ${{trim}} -e 0.125 -j {threads} --max-n=1 -m " + config['readsLength_min'] + " -M " + config['readsLength_max'] + " -o {output.cut_fastq} {input.fastq} 1>> {log.cutadapt_out} 2> {log.cutadapt} ;"

# ---------------------------------------------------------------
#               Rule: FastQC After Trimming
# ---------------------------------------------------------------
# Quality control of data : build of the fastqc after adapter trimming
rule make_fastqc_after_trimming:
    input:
        results_path + "cutadapt/{sample}.cutadapt" + frag_length_L + ".fastq.gz"
    output:
        results_path + "fastqc/fastqc_after_trimming/{sample}.cutadapt" + frag_length_L + "_fastqc.zip",
        results_path + "fastqc/fastqc_after_trimming/{sample}.cutadapt" + frag_length_L + "_fastqc.html"
    log:
        logs_path + "make_fastqc_after_trimming/{sample}.log"
    benchmark:
        local_path + "benchmarks/make_fastqc_after_trimming/{sample}.benchmark.txt"
    params:
       outdir = results_path + "fastqc/fastqc_after_trimming/"
    shell:
        "fastqc {input} --outdir {params.outdir} 2> {log}"

# ---------------------------------------------------------------
#               Rule: Bowtie2 Build for OutRNA
# ---------------------------------------------------------------
# Builds the index of bowtie2 mapping on sequences for reads remove
rule bowtie2_build_outRNA:
    input:
        outRNA = local_path + "database/" + os.path.basename(config['fasta_outRNA'])
    output:
        results_path + "annex_database/index_files/outRNA_bowtie2.1.bt2"
    log:
        logs_path + "bowtie2_build_outRNA/bowtie2_build_outRNA.log"
    benchmark:
        local_path + "benchmarks/bowtie2_build_outRNA/bowtie2_build_outRNA.benchmark.txt"
    threads:
        multi_threads_nbr
    params:
        outNames = results_path + "annex_database/index_files/outRNA_bowtie2"
    shell:
        "bowtie2-build --threads {threads} {input.outRNA} {params.outNames} &> {log}"

# ---------------------------------------------------------------
#               Rule: Bowtie2 Run for OutRNA Depletion
# ---------------------------------------------------------------
# Mapping of non-coding RNA
rule bowtie2_run_outRNA:
    input:
        results_path + "annex_database/index_files/outRNA_bowtie2.1.bt2",
        fastq = results_path + "cutadapt/{sample}.cutadapt" + frag_length_L + ".fastq.gz"
    output:
        results_path + "no-outRNA/{sample}" + frag_length_L + ".no-outRNA.fastq.gz"  
    log:
        bt2 = stats_path + "{sample}_bowtie2_run_outRNA.log"
    benchmark:
        local_path + "benchmarks/bowtie2_run_outRNA/{sample}.benchmark.txt"
    resources:
        mem_mb = mem_mb_resources
    threads:
        multi_threads_nbr
    shell:
        "bowtie2 -x " + results_path + "annex_database/index_files/outRNA_bowtie2 --threads {threads} -U {input.fastq} --un-gz {output} > /dev/null 2>> {log.bt2}"

# ---------------------------------------------------------------
#               Rule: FastQC After OutRNA Depletion
# ---------------------------------------------------------------
# Quality control of data : build of the fastqc after depletin rRNA as ribosomal RNA can have an impact on the data's profiles (ATGC content)
rule make_fastqc_after_outRNA_depletion:
    input:
        results_path + "no-outRNA/{sample}" + frag_length_L + ".no-outRNA.fastq.gz"
    output:
        results_path + "fastqc/make_fastqc_after_outRNA_depletion/{sample}" + frag_length_L + ".no-outRNA_fastqc.zip",
        results_path + "fastqc/make_fastqc_after_outRNA_depletion/{sample}" + frag_length_L + ".no-outRNA_fastqc.html"
    log:
        logs_path + "make_fastqc_after_outRNA_depletion/{sample}.log"
    benchmark:
        local_path + "benchmarks/make_fastqc_after_outRNA_depletion/{sample}.benchmark.txt"
    params:
       outdir = results_path + "fastqc/make_fastqc_after_outRNA_depletion/"
    shell:
        "fastqc {input} --outdir {params.outdir} 2> {log}"

# ---------------------------------------------------------------
#               Rule: Bowtie2 Build for Transcriptome
# ---------------------------------------------------------------
# Builds the index of bowtie2 mapping for all RNA
rule bowtie2_build:
    input:
        fasta = local_path + "database/" + os.path.basename(config['fasta'])
    output:
        results_path + "annex_database/index_files/index_bowtie2.1.bt2"
    log:
        logs_path + "bowtie2_build/bowtie2_build.log"
    benchmark:
        local_path + "benchmarks/bowtie2_build/bowtie2_build.benchmark.txt"
    threads:
        multi_threads_nbr
    params:
        outNames = results_path + "annex_database/index_files/index_bowtie2"
    shell:
        "bowtie2-build --threads {threads} {input.fasta} {params.outNames} &> {log}"

# ---------------------------------------------------------------
#               Rule: Hisat2 Build for Transcriptome
# ---------------------------------------------------------------
# Builds the index of hisat2 mapping for all RNA
rule hisat2_build:
    input:
        fasta = local_path + "database/" + os.path.basename(config['fasta'])
    output:
        results_path + "annex_database/index_files/index_hisat2.1.ht2"
    log:
        logs_path + "hisat2_build/hisat2_build.log"
    benchmark:
        local_path + "benchmarks/hisat2_build/hisat2_build.benchmark.txt"
    threads:
        multi_threads_nbr
    params:
        outNames = results_path + "annex_database/index_files/index_hisat2"
    shell:
        "hisat2-build --threads {threads} {input.fasta} {params.outNames} &> {log}"


# ---------------------------------------------------------------
#               Rule: Filter 5-prime UTR
# ---------------------------------------------------------------
# Quality controls :
# Filter the transcripts to only keep those with a 5-prime UTR if there are enough of them
rule filter_by_5UTR:
    input:
        gff_namedCDS = results_path + "annex_database/NamedCDS_" + os.path.basename(config['gff'])
    output:
        gff_5UTR_filtered = temp(results_path + "annex_database/gff_5UTR_filtered_" + os.path.basename(config['gff'])),
        gff_features = results_path + "annex_database/gff_features_counts.txt"
    log:
        logs_path + "filter_by_5UTR/filtering.log"
    benchmark:
        local_path + "benchmarks/filter_by_5UTR/filtering.benchmark.txt"
    params:
        annex_path = results_path + "annex_database/"
    shell:
        "bash " + RibokastIndex_tools + "others/filter_5UTR.sh -g {input.gff_namedCDS} -o {output.gff_5UTR_filtered} -p {params.annex_path} -r " + config["gff_mRNA_feature"] + " -u " + config["gff_5UTR_feature"] + " -t " + utr_threshold + " 2> {log}"

# ---------------------------------------------------------------
#               Rule: ORFget
# ---------------------------------------------------------------
# In case there are no UTRs in the original GFF, call for ORFget functions
rule ORFget:
    input:
        fasta = local_path + "database/" + os.path.basename(config['fasta']),
        gff = results_path + "annex_database/gff_5UTR_filtered_" + os.path.basename(config['gff'])
    output:
        fasta = results_path + "annex_database/transcriptome_elongated.nfasta",
        gff = results_path + "annex_database/transcriptome_elongated.gff"
    log:
        orf_get = logs_path + "ORFget/orf_get.log"
    benchmark:
        local_path + "benchmarks/ORFget/ORFget.benchmark.txt"
    params:
        path = results_path + "annex_database/transcriptome"
    shell:
        "python3 " + RibokastIndex_tools + "orfget/ORFget.py -fna {input.fasta} -gff {input.gff} -features_include " + config['gff_cds_feature'] + " -name_attribute " + config['gff_name_attribut'] + " -o {params.path} -type nucl -elongate 50 -check 2> {log.orf_get}"
# ---------------------------------------------------------------
#               Rule: Create GTF file for riboWaltz
# ---------------------------------------------------------------
#Create GTF file for riboWaltz
rule transcriptome_construction_gtf:
    input:
        fasta = results_path + "annex_database/transcriptome_elongated.nfasta",
        gff = results_path + "annex_database/transcriptome_elongated.gff"
    output:
        fasta = results_path + "annex_database/transcriptome_elongated.exons_" + os.path.basename(config['fasta']),
        gtf = results_path + "annex_database/transcriptome_elongated.exons_" + os.path.basename(config['gff']) + ".gtf"
    log:
        samtools_index = logs_path + "transcriptome_construction_gtf/samtools_index.log",
        gffread_gtf = logs_path + "transcriptome_construction_gtf/gffread.log",
        sed = logs_path + "transcriptome_construction_gtf/sed.log",
        awk = logs_path + "transcriptome_construction_gtf/awk.log"
    benchmark:
        local_path + "benchmarks/transcriptome_construction_gtf/transcriptome_construction_gtf.benchmark.txt"
    params:
        tmp_gtf = results_path + "annex_database/transcriptome_elongated.exons_tmp.gtf"
    shell:
        "samtools faidx {input.fasta} 2> {log.samtools_index} ;"
        "gffread -F -T -w {output.fasta} -o {params.tmp_gtf} -g {input.fasta} {input.gff} 2> {log.gffread_gtf} ;"
        "sed -i 's/description[^\;]*\;//' {params.tmp_gtf} 2>> {log.sed} ;"
        "sed -i 's/\\t[A-Z]*[_]*gene_segment\\t/\\ttranscript\\t/' {params.tmp_gtf} 2>> {log.sed} ;"
        """awk -F '\\t' '{{if(NF<=9) {{print($0);}} else {{for(field=1;field<9;field++) {{printf("%s\\t",$field);}} for(field=9;field<=NF;field++) {{printf("%s ",$field);}} printf("\\n");}}}}' {params.tmp_gtf} > {output.gtf} 2>> {log.awk} ;"""
        "rm -f {params.tmp_gtf} ;"
        "sed -i 's/\\s$//' {output.gtf} 2>> {log.sed} ;"
# ---------------------------------------------------------------
#               Rule: ORFget
# ---------------------------------------------------------------
# Builds the index of bowtie2 mapping for all RNA
rule bowtie2_build_transcriptome:
    input:
        fasta = results_path + "annex_database/transcriptome_elongated.nfasta"
    output:
        results_path + "annex_database/index_files/transcriptome_elongated.transcriptome_index_bowtie2.1.bt2"
    log:
        logs_path + "bowtie2_build_transcriptome/bowtie2_build.log"
    benchmark:
        local_path + "benchmarks/bowtie2_build_transcriptome/transcriptome_index_bowtie2.benchmark.txt"
    params:
        index_names = results_path + "annex_database/index_files/transcriptome_elongated.transcriptome_index_bowtie2"
    threads:
        multi_threads_nbr
    shell:
        "bowtie2-build --threads {threads} {input.fasta} {params.index_names} &> {log}"
        
# ---------------------------------------------------------------
#               Rule: Bowtie2 Build for Transcriptome
# ---------------------------------------------------------------
# Builds the index of hisat2 mapping for all RNA
rule hisat2_build_transcriptome:
    input:
        fasta = results_path + "annex_database/transcriptome_elongated.nfasta"
    output:
        results_path + "annex_database/index_files/transcriptome_elongated.transcriptome_index_hisat2.1.ht2"
    log:
        logs_path + "hisat2_build_transcriptome/hisat2_build.log"
    benchmark:
        local_path + "benchmarks/hisat2_build_transcriptome/hisat2_build_transcriptome.benchmark.txt"
    params:
        index_names = results_path + "annex_database/index_files/transcriptome_elongated.transcriptome_index_hisat2"
    threads:
        multi_threads_nbr
    shell:
        "hisat2-build --threads {threads} {input.fasta} {params.index_names} &> {log}"

# ---------------------------------------------------------------
#               Rule: Hisat2 Build for Transcriptome
# ---------------------------------------------------------------
# Performs mapping on transcriptome
rule run_mapping_transcriptome:
    input:
        results_path + "annex_database/index_files/transcriptome_elongated.transcriptome_index_hisat2.1.ht2",
        results_path + "annex_database/index_files/transcriptome_elongated.transcriptome_index_bowtie2.1.bt2",
        fastq = results_path + "cutadapt/{sample}.cutadapt" + frag_length_L + ".fastq.gz" if config['fasta_outRNA']=="" else results_path + "no-outRNA/{sample}" + frag_length_L + ".no-outRNA.fastq.gz"
    output:
        sam_hisat2 = temp(results_path + "BAM_transcriptome" + frag_length_L + "/{sample}" + frag_length_L + ".hisat2.sam"),
        sam_bowtie2 = temp(results_path + "BAM_transcriptome" + frag_length_L + "/{sample}" + frag_length_L + ".bowtie2.sam"),
        fastq = temp(results_path + "no-outRNA/{sample}" + frag_length_L + ".no-outRNA.notAlign.transcriptome.fastq.gz")
    log:
        hisat2_out = logs_path + "run_mapping_transcriptome/{sample}_run_mapping_transcriptome_hisat2.log",
        bowtie2_out = logs_path + "run_mapping_transcriptome/{sample}_run_mapping_transcriptome_bowtie2.log"
    benchmark:
        local_path + "benchmarks/run_mapping_transcriptome/{sample}.benchmark.txt"
    resources:
        mem_mb = mem_mb_resources
    params:
        index_names_hisat2 = results_path + "annex_database/index_files/transcriptome_elongated.transcriptome_index_hisat2",
        index_names_bowtie2 = results_path + "annex_database/index_files/transcriptome_elongated.transcriptome_index_bowtie2"
    threads:
        multi_threads_nbr
    shell:
        "hisat2 -x {params.index_names_hisat2} --threads {threads} -k " + isoform_nbr + " --no-softclip --no-unal -U {input.fastq} --un-gz {output.fastq} -S {output.sam_hisat2} 2>> {log.hisat2_out} ;"
        "bowtie2 -x {params.index_names_bowtie2} --threads {threads} -k " + isoform_nbr + " --end-to-end --no-unal --no-hd -U {output.fastq} -S {output.sam_bowtie2} 2>> {log.bowtie2_out}"

# ---------------------------------------------------------------
#               Rule: Samtools Filter and Sort BAM
# ---------------------------------------------------------------
# Creates bam and sam files
rule samtools_filter_transcriptome:
    priority: 50
    input:
        sam_hisat2 = results_path + "BAM_transcriptome" + frag_length_L + "/{sample}" + frag_length_L + ".hisat2.sam",
        sam_bowtie2 = results_path + "BAM_transcriptome" + frag_length_L + "/{sample}" + frag_length_L + ".bowtie2.sam"
    output:
        bam = results_path + "BAM_transcriptome" + frag_length_L + "/transcriptome_elongated.{sample}" + frag_length_L + ".bam"
    log:
        view_header_hisat2 = logs_path + "samtools_filter_transcriptome/{sample}.grep_hisat2.log",
        uniq_header = logs_path + "samtools_filter_transcriptome/{sample}.uniq_header.log",
        view_core_hisat2 =  logs_path + "samtools_filter_transcriptome/{sample}.grep_core_hisat2.log",
        XM_filter_hisat2 = logs_path + "samtools_filter_transcriptome/{sample}.XM_filter_hisat2.log",
        XM_filter_bowtie2 = logs_path + "samtools_filter_transcriptome/{sample}.XM_filter_bowtie2.log",
        view_bam = logs_path + "samtools_filter_transcriptome/{sample}.view_bam.log",
        sort_bam = logs_path + "samtools_filter_transcriptome/{sample}.sort_bam.log"
    benchmark:
        local_path + "benchmarks/samtools_filter_transcriptome/{sample}.benchmark.txt"
    resources:
        mem_mb = round(mem_mb_resources / 3)
    params:
        sample = "{sample}",
        sam = results_path + "BAM_transcriptome" + frag_length_L + "/{sample}" + frag_length_L + ".sam"
    threads:
        multi_threads_nbr
    shell:
        "set +o pipefail ;"
        "samtools view -H {input.sam_hisat2}  2> {log.view_header_hisat2} | uniq 2> {log.uniq_header} 1> {params.sam} ;"
        "samtools view {input.sam_hisat2} 2> {log.view_core_hisat2} | egrep -vi 'XM:i:[2-9]' 2> {log.XM_filter_hisat2} 1>> {params.sam} ;"
        "egrep -vi 'XM:i:[2-9]' {input.sam_bowtie2} 2> {log.XM_filter_bowtie2} 1>> {params.sam} ;"
        "samtools view -@ {threads} -F 3588 -h -b {params.sam} 2> {log.view_bam} | samtools sort -@ {threads} -o {output.bam} 2> {log.sort_bam} ;"
        "rm -f {params.sam} ;"
        
# ---------------------------------------------------------------
#               Rule: BAM Indexing
# ---------------------------------------------------------------
# Index BAMs
rule index_bam_transcriptome:
    input:
        bam = results_path + "BAM_transcriptome" + frag_length_L + "/transcriptome_elongated.{sample}" + frag_length_L + ".bam"
    output:
        bai = results_path + "BAM_transcriptome" + frag_length_L + "/transcriptome_elongated.{sample}" + frag_length_L + ".bam.bai"
    log:
        index_bam = logs_path + "index_bam_transcriptome/{sample}.index_bam.log"
    benchmark:
        local_path + "benchmarks/index_bam_transcriptome/{sample}.benchmark.txt"
    shell:
        "samtools index {input.bam} 2> {log.index_bam}"

# ---------------------------------------------------------------
#               Rule: riboWaltz - Qualitative Analysis
# ---------------------------------------------------------------
rule riboWaltz_transcriptome:
    input:
        transcriptome_gtf = results_path + "annex_database/transcriptome_elongated.exons_" + os.path.basename(config['gff']) + ".gtf",
        transcriptome_bam = results_path + "BAM_transcriptome" + frag_length_L + "/transcriptome_elongated.{sample}" + frag_length_L + ".bam"
    output:
        psite_table = results_path + "riboWaltz" + frag_length_L + "/{sample}/psite_offset.csv",
        psite_length = results_path + "riboWaltz" + frag_length_L + "/{sample}/frame_psite_length.csv",
        frame_psite = results_path + "riboWaltz" + frag_length_L + "/{sample}/frame_psite.csv",
        best_offset = results_path + "riboWaltz" + frag_length_L + "/{sample}/best_offset.csv"
    log:
        periodicity = logs_path + "riboWaltz_transcriptome/{sample}.riboWaltz.log"
    resources:
        mem_mb = mem_mb_resources * 2
    benchmark:
        local_path + "benchmarks/riboWaltz_transcriptome/{sample}.riboWaltz.benchmark.txt"
    shell:
        """
        /usr/bin/time -v bash -c '
        Rscript {RibokastIndex_tools}ribowaltz/periodicity_riboWaltz_transcriptome.R \
            -w {local_path} \
            -f {RibokastIndex_tools}ribowaltz/ \
            -g {input.transcriptome_gtf} \
            --path_prefix {results_path} \
            -b {input.transcriptome_bam} \
            -n {wildcards.sample}
        ' 2> {log.periodicity}
        """

rule concat_riboWaltz_outputs:
    input:
        frame_psite_length = expand(results_path + "riboWaltz" + frag_length_L + "/{sample}/frame_psite_length.csv", sample=SAMPLES),
        psite_offset       = expand(results_path + "riboWaltz" + frag_length_L + "/{sample}/psite_offset.csv", sample=SAMPLES),
        frame_psite        = expand(results_path + "riboWaltz" + frag_length_L + "/{sample}/frame_psite.csv", sample=SAMPLES),
        best_offset        = expand(results_path + "riboWaltz" + frag_length_L + "/{sample}/best_offset.csv", sample=SAMPLES)
    output:
        frame_psite_length_merged = results_path + "riboWaltz" + frag_length_L + "/frame_psite_length.csv",
        psite_offset_merged       = results_path + "riboWaltz" + frag_length_L + "/psite_offset.csv",
        frame_psite_merged        = results_path + "riboWaltz" + frag_length_L + "/frame_psite.csv",
        best_offset_merged        = results_path + "riboWaltz" + frag_length_L + "/best_offset.csv"
    shell:
        # Frame P-site Length
        "(head -n 1 {input.frame_psite_length[0]} && tail -n +2 -q {input.frame_psite_length}) > {output.frame_psite_length_merged} ; " 
        # P-site Offset
        "(head -n 1 {input.psite_offset[0]} && tail -n +2 -q {input.psite_offset}) > {output.psite_offset_merged} ; " 
        # Frame P-site
        "(head -n 1 {input.frame_psite[0]} && tail -n +2 -q {input.frame_psite}) > {output.frame_psite_merged} ; " 
        # Best Offset
        "(head -n 1 {input.best_offset[0]} && tail -n +2 -q {input.best_offset}) > {output.best_offset_merged} ;"
       

#####################################################################################################################################################################################
#                                                         Kmer_Counts_Index
##################################################################################################################################################################################### 

checkpoint get_samples_for_kmer:
    input:
        psite_table = results_path + "riboWaltz" + frag_length_L + "/psite_table_forKmerCount.txt"
    output:
        samples_list = results_path + "kmerCount/samples_kmer.txt"
    run:
        samples = set()
        with open(input.psite_table, "r") as infile:
            for line in infile:
                fields = line.strip().split()
                if len(fields) < 1:
                    continue
                samples.add(fields[0])
        with open(output.samples_list, "w") as outfile:
            for sample in sorted(samples):
                outfile.write(sample + "\n")

rule creationKmer_trigger:
    input:
        get_kmer_tsv_paths
# ---------------------------------------------------------------
#               Rule: Create Phase Table
# ---------------------------------------------------------------      
# File for kmerCount

rule create_phase_table:
    input:
        psite_table  = results_path + "riboWaltz" + frag_length_L + "/psite_offset.csv",
        psite_length = results_path + "riboWaltz" + frag_length_L + "/frame_psite_length.csv"
    output:
        psite_table_offset      = results_path + "riboWaltz" + frag_length_L + "/psite_table_offset.csv",
        psite_table_for_KmerCount = results_path + "riboWaltz" + frag_length_L + "/psite_table_forKmerCount.txt"
    resources:
        mem_mb = mem_mb_resources * 100
    benchmark:
        local_path + "benchmarks/riboWaltz_transcriptome/riboWaltz.benchmark.txt"
    shell:
        """
        python3 {RibokastIndex_tools}ribowaltz/filter_merge_for_KmerCount.py \
        {input.psite_length} {input.psite_table} \
        {output.psite_table_offset} {output.psite_table_for_KmerCount} \
        {config[kmercount_pct_threshold]}
        """


# ---------------------------------------------------------------
#               Rule: K-mer Creation
# ---------------------------------------------------------------

rule creationKmer:
    input:
        fastq_norRNA = results_path + "no-outRNA/{sample}" + frag_length_L + ".no-outRNA.fastq.gz",
        psite_table_for_KmerCount = results_path + "riboWaltz" + frag_length_L + "/psite_table_forKmerCount.txt",
        sample_list = results_path + "kmerCount/samples_kmer.txt"
    output:
        vector = results_path + "kmerCount/Kmer/{sample}.tsv"
    log:
        creationKmer = logs_path + "creationKmer/{sample}_creationKmer.log"
    params:
        # forced_phase can be int, or null dans YAML
        forced_phase = config.get("forced_phase", None),
        mode = config["mode"]
    shell:
        r"""
        FORCE_ARG=""
        if [ "{params.mode}" = "phase" ] && [ "{params.forced_phase}" != "None" ] && [ -n "{params.forced_phase}" ]; then
            FORCE_ARG="--forced-phase {params.forced_phase}"
        fi

        /usr/bin/time -v python3 {RibokastIndex_tools}kmers_index/processCountingKmer_phase.py \
            -i {input.fastq_norRNA} \
            -o {results_path}kmerCount \
            -k {config[kmerSize]} \
            -n {wildcards.sample} \
            -m {config[mode]} \
            -t {input.psite_table_for_KmerCount} \
            --fastq-name {wildcards.sample} \
            $FORCE_ARG \
            2> {log.creationKmer}
        """

# ---------------------------------------------------------------
#               Rule: Join Counts (KaMRaT)
# ---------------------------------------------------------------

pathJoinCounts=config['pathJoinCounts']
rule joincounts:
    input:
        get_kmer_tsv_paths
    output:
        matrix = temp(results_path + "kmerCount/Matrix/matrix.tsv")
    log:
        joinCounts = logs_path + "joinCounts/joinCounts.log"
    shell:
        """
        export PATH={pathJoinCounts};
        joinCounts {input} > {output.matrix} 2> {log.joinCounts}
        """
       
# ---------------------------------------------------------------
#               Rule: Filter Matrix
# ---------------------------------------------------------------
rule suppressionLettre:
    input:
         matrix= results_path + "kmerCount/Matrix/matrix.tsv"
    output:
         matrix= temp(results_path + "kmerCount/Matrix/matrixDNA.tsv")
    shell:
        (
            "python3 " + RibokastIndex_tools + "kmers_index/removeNonDnaLetter.py "
            "-i {input} "
            "-o {output} "
        )

rule filter:
    input:
        matrix = results_path + "kmerCount/Matrix/matrixDNA.tsv"
    output:
        matrix_filtered = temp(results_path + "kmerCount/Matrix/matrixFiltered.tsv")
    log:
        sort_log = logs_path + "filter/sort.log",
        filter_log = logs_path + "filter/filter_py.log"
    shell:
        """
        sort -k1,1 {input.matrix} > {input.matrix}.sorted 2> {log.sort_log}
        python3 {RibokastIndex_tools}kmers_index/filter.py {input.matrix}.sorted {output.matrix_filtered} 2> {log.filter_log}
        rm -f {input.matrix}.sorted
        """

rule header:
    input:
        matrix_filtered=results_path + "kmerCount/Matrix/matrixFiltered.tsv"
    output:
        matrix_filtered_header=results_path + "kmerCount/Matrix/matrixFilteredHeader.tsv"

    shell:
        (
          "bash " + RibokastIndex_tools + "kmers_index/scriptAddHeader.sh "
          + results_path + "kmerCount/Kmer "
          + "{input.matrix_filtered} "
          + "{output.matrix_filtered_header}"
        )

# ---------------------------------------------------------------
#               Rule: KaMRaT - K-mer Index Merging
# ---------------------------------------------------------------
kamratImg = config['kamratImg']

rule kamrat:
    input:
        matrix_filtered_header = results_path + "kmerCount/Matrix/matrixFilteredHeader.tsv"
    output:
        results_path + "kmerCount/Kamrat/merged-res.tsv"
    log:
        kamrat_index_merge = logs_path + "joinCounts/kamrat_index_merge.log"

    shell:
        r"""
        bash {RibokastIndex_tools}kmers_index/scriptKamrat.sh \
          {kamratImg} \
          {input.matrix_filtered_header} \
          {results_path}kmerCount/Kamrat \
          {config[kmerSize]} \
          {config[kamrat_normalize]} \
          {config[kamrat_nfbase]} \
          > {log.kamrat_index_merge} 2>&1
        """
