#This pipeline will produce a number of test data sets and perfom analysis steps to determine the development of cancer over time

configfile: 'config_ctSimu.yaml'

patients = [str(n+1) for n in range(0, config['patients'])]
trends = config['trends']
timepoints = config['timepoints']
varcallers = config['varcallers']

#run all the rest
rule all:
    input:
        expand("../lineplots/{trend}_{patient}/{varcaller}.png", trend = trends, patient = patients, varcaller = varcallers),
        expand("../after_amplification/{trend}_{patient}/{trend}_{patient}_{timepoint}/ampli_bwa_sorted.qcml", trend = trends, patient = patients, timepoint = timepoints),
        expand("../lineplots/{trend}_{patient}/truth.png", trend = trends, patient = patients)

#create BED file for LoFreq and VarDict from input VCF
rule vcf_to_bed:
    input:
        config['template_vcf']
    output:
        config['template_bed']
    params:
        config['vcf_to_bed_extend']
    shell:
        """
        {params} -i {input} -o {output} -l 1 -d u
        """

rule simulate_vcfs:
    input:
        config['template_vcf']
    output:
        files = "../simulated_vcfs/{trend}_{patient}/{trend}_{patient}_{timepoint}_simu.vcf"
    params:
        sim = config['simulate_development'],
        var = config['variation']
    shell:
        """
        {params.sim} -i {input} -o {output} -t {wildcards.trend} -p {wildcards.timepoint} -v {params.var}
        """

#run NEAT for the simulation of healthy reads
rule neat_healthy:
    input:
        config['regions_bed']
    params:
        neat = config['neat'],
        coverage = config['neat_coverage'],
        reference = config['neat_reference_sim'],
        read_len = config['neat_read_length'],
        pe = config['neat_pe'],
        mut_rate = config['neat_mut_rate'],
        err_rate = config['neat_err_rate'],
        offtarget = config['neat_offtarget'],
        outpatient = "../neat_output/neat_healthy"
    output:
        f = "../neat_output/neat_healthy_read1.fq.gz",
        r = "../neat_output/neat_healthy_read2.fq.gz",
        sorted = "../neat_output/neat_healthy_golden_sorted.bam"
    shell:
        """
        {params.neat} -c {params.coverage} -r {params.reference} -R {params.read_len} \
        -o {params.outpatient} -tr {input} --pe {params.pe} --bam --vcf -M {params.mut_rate} -E {params.err_rate} \
        -to {params.offtarget}
        samtools sort {params.outpatient}_golden.bam > {output.sorted}
        samtools index {output.sorted}
        """
    
#run NEAT for the simulation of cancer reads
rule neat_cancer:
    input:
        bed = config['regions_bed'],
        vcf = config['template_vcf']
    params:
        neat = config['neat'],
        coverage = config['neat_coverage'],
        reference = config['neat_reference_sim'],
        read_len = config['neat_read_length'],
        pe = config['neat_pe'],
        mut_rate = config['neat_mut_rate'],
        err_rate = config['neat_err_rate'],
        offtarget = config['neat_offtarget'],
        outpatient = "../neat_output/neat_cancer"
    output:
        f = "../neat_output/neat_cancer_read1.fq.gz",
        r = "../neat_output/neat_cancer_read2.fq.gz",
        sorted = "../neat_output/neat_cancer_golden_sorted.bam"
    shell:
        """
        {params.neat} -c {params.coverage} -r {params.reference} -R {params.read_len} \
        -o {params.outpatient} -tr {input.bed} --pe {params.pe} -v {input.vcf} --bam --vcf -M {params.mut_rate}\
        -E {params.err_rate} -to {params.offtarget}
        samtools sort {params.outpatient}_golden.bam > {output.sorted}
        samtools index {output.sorted}
        """

#create sample FASTQs with variant frequencies based on the simulation
rule mix_sample:
    input:
        f = "../neat_output/neat_healthy_read1.fq.gz",
        r = "../neat_output/neat_healthy_read2.fq.gz",
        fc = "../neat_output/neat_cancer_read1.fq.gz",
        rc = "../neat_output/neat_cancer_read2.fq.gz",
        hg = "../neat_output/neat_healthy_golden_sorted.bam",
        cg =  "../neat_output/neat_cancer_golden_sorted.bam",
        umis = config['umis'],
        muts = "../simulated_vcfs/{trend}_{patient}/{trend}_{patient}_{timepoint}_simu.vcf"
    output:
        fo = "../before_amplification/{trend}_{patient}/{trend}_{patient}_{timepoint}/L001_R1_001.fastq.gz",
        ro = "../before_amplification/{trend}_{patient}/{trend}_{patient}_{timepoint}/L001_R2_001.fastq.gz",
        truth = "../call_vcfs/{trend}_{patient}/{trend}_{patient}_{timepoint}/truth.vcf"
    params:
        config['mix_cancer_reads']
    shell:
        """
        {params} -f {input.f} -r {input.r} -fc {input.fc} -rc {input.rc} -fo {output.fo} -ro {output.ro} \
        -hg {input.hg} -cg {input.cg} -u {input.umis} -m {input.muts} -t {output.truth}
        """

#run in silico amplification (and log mutations)
rule amplify:
    input:
        f = "../before_amplification/{trend}_{patient}/{trend}_{patient}_{timepoint}/L001_R1_001.fastq.gz",
        r = "../before_amplification/{trend}_{patient}/{trend}_{patient}_{timepoint}/L001_R2_001.fastq.gz"
    output:
        fo = "../after_amplification/{trend}_{patient}/{trend}_{patient}_{timepoint}/ampli_L001_R1_001.fastq.gz",
        ro = "../after_amplification/{trend}_{patient}/{trend}_{patient}_{timepoint}/ampli_L001_R2_001.fastq.gz",
        fl = "../mutation_logs/{trend}_{patient}/{trend}_{patient}_{timepoint}/L001_R1_001.txt",
        rl = "../mutation_logs/{trend}_{patient}/{trend}_{patient}_{timepoint}/L001_R2_001.txt"
    params:
        amplification = config['amplification'],
        multi = config['multi'],
        preerror = config['preerror'],
        error = config['error']
    shell:
        """
        {params.amplification} -f {input.f} -r {input.r} -fo {output.fo} -ro {output.ro} -fl {output.fl} -rl {output.rl} -m {params.multi} -p {params.preerror} -e {params.error}
        """

#map and sort amplified reads
rule mapsort_amplified:
    input:
        f = "../after_amplification/{trend}_{patient}/{trend}_{patient}_{timepoint}/ampli_L001_R1_001.fastq.gz",
        r = "../after_amplification/{trend}_{patient}/{trend}_{patient}_{timepoint}/ampli_L001_R2_001.fastq.gz"
    output:
        "../after_amplification/{trend}_{patient}/{trend}_{patient}_{timepoint}/ampli_bwa_sorted.bam"
    threads:
        config['threads']
    params:
        config['reference']
    shell:
        """
        bwa mem -t {threads} {params} {input.f} {input.r} | samtools sort -@{threads} -o {output} -
        """

#add UMIs (those are removed by bwa), index and dedup amplified reads
rule index_dedup_amplified:
    input:
        bam = "../after_amplification/{trend}_{patient}/{trend}_{patient}_{timepoint}/ampli_bwa_sorted.bam",
        fq1 = "../before_amplification/{trend}_{patient}/{trend}_{patient}_{timepoint}/L001_R1_001.fastq.gz",
        fq2 = "../before_amplification/{trend}_{patient}/{trend}_{patient}_{timepoint}/L001_R2_001.fastq.gz"
    output:
        bai = "../after_amplification/{trend}_{patient}/{trend}_{patient}_{timepoint}/ampli_bwa_sorted.bam.bai",
        deduped = "../after_amplification/{trend}_{patient}/{trend}_{patient}_{timepoint}/bwa_sorted_umis_dedup.bam",
        deduped_bai = "../after_amplification/{trend}_{patient}/{trend}_{patient}_{timepoint}/bwa_sorted_umis_dedup.bam.bai"
    params:
        config['dedup']
    shell:
        """
        samtools index {input.bam}
        {params} --infile {input.bam} --outfile {output.deduped}
        samtools index {output.deduped}
        """

#run mappingQC for a quality check
rule mappingQC:
    input:
        bam = "../after_amplification/{trend}_{patient}/{trend}_{patient}_{timepoint}/ampli_bwa_sorted.bam",
        bai = "../after_amplification/{trend}_{patient}/{trend}_{patient}_{timepoint}/ampli_bwa_sorted.bam.bai"
    output:
        qcml = "../after_amplification/{trend}_{patient}/{trend}_{patient}_{timepoint}/ampli_bwa_sorted.qcml"
    params:
        mappingQC = config['mappingQC'],
        roi = config['regions_bed'],
        reference = config['reference']
    shell:
        """
        {params.mappingQC} -in {input.bam} -out {output.qcml} -roi {params.roi} -ref {params.reference} -cfdna
        """

#map and sort unamplified reads for error finding
rule mapsort_unamplified:
    input:
        f = "../before_amplification/{trend}_{patient}/{trend}_{patient}_{timepoint}/L001_R1_001.fastq.gz",
        r = "../before_amplification/{trend}_{patient}/{trend}_{patient}_{timepoint}/L001_R2_001.fastq.gz"
    output:
        "../before_amplification/{trend}_{patient}/noampli_{timepoint}.bam"
    threads:
        config['threads']
    params:
        config['reference']
    shell:
        """
        bwa mem -t {threads} {params} {input.f} {input.r} | samtools sort -@{threads} -o {output} -
        """

############################################################################### Variant Calling!
#the following four rules each run the variant calling with a different program
rule lofreq_call:
    benchmark:      ########################################################### Remove????
        "../call_logs/{trend}_{patient}/{trend}_{patient}_{timepoint}/lofreq_log.tsv"
    input:
        bed = config['template_bed'],
        bam = "../after_amplification/{trend}_{patient}/{trend}_{patient}_{timepoint}/bwa_sorted_umis_dedup.bam",
        bai = "../after_amplification/{trend}_{patient}/{trend}_{patient}_{timepoint}/bwa_sorted_umis_dedup.bam.bai"
    output:
        "../call_vcfs/{trend}_{patient}/{trend}_{patient}_{timepoint}/lofreq.vcf"
    params:
        lofreq = config['lofreq'],
        reference = config['reference']
    shell:
        """
        {params.lofreq} -f {params.reference} -o {output} -l {input.bed} {input.bam}
        """

rule umivar_call:
    benchmark:      ########################################################### Remove????
        "../call_logs/{trend}_{patient}/{trend}_{patient}_{timepoint}/umivar_log.tsv"
    input:
        moni = "../simulated_vcfs/{trend}_{patient}/{trend}_{patient}_before_simu.vcf",
        bam = "../after_amplification/{trend}_{patient}/{trend}_{patient}_{timepoint}/bwa_sorted_umis_dedup.bam",
        bai = "../after_amplification/{trend}_{patient}/{trend}_{patient}_{timepoint}/bwa_sorted_umis_dedup.bam.bai"
    output:
        folder = directory("../call_vcfs/{trend}_{patient}/{trend}_{patient}_{timepoint}/umivar"),
        bed = "../call_vcfs/{trend}_{patient}/{trend}_{patient}_{timepoint}/extended_60.bed",
        vcf = "../call_vcfs/{trend}_{patient}/{trend}_{patient}_{timepoint}/umivar.vcf"
    params:
        vcf_to_bed_extend = config['vcf_to_bed_extend'],
        umivar = config['umivar'],
        reference = config['reference'],
        umi_file_name = config['umi_file_name']
    shell:
        """
        {params.vcf_to_bed_extend} -i {input.moni} -o {output.bed} -d b -l 60
        {params.umivar} -tbam {input.bam} -r {params.reference} -o {output.folder} -m {input.moni} -b {output.bed}
        cp {output.folder}/{params.umi_file_name} {output.vcf}
        """

rule vardict_call:
    benchmark:      ########################################################### Remove????
        "../call_logs/{trend}_{patient}/{trend}_{patient}_{timepoint}/vardict_log.tsv"
    input:
        bed = config['template_bed'],
        bam = "../after_amplification/{trend}_{patient}/{trend}_{patient}_{timepoint}/bwa_sorted_umis_dedup.bam",
        bai = "../after_amplification/{trend}_{patient}/{trend}_{patient}_{timepoint}/bwa_sorted_umis_dedup.bam.bai"
    output:
        "../call_vcfs/{trend}_{patient}/{trend}_{patient}_{timepoint}/vardict.vcf"
    params:
        vardict = config['vardict'],
        teststrandbias = config['teststrandbias'],
        var2vcf_valid = config['var2vcf_valid'],
        reference = config['reference']
    shell:
        """
        {params.vardict} -G {params.reference} -f 0.001 -N {wildcards.patient} -b {input.bam} -c 1 -S 2 -E 3 {input.bed} \
        | {params.teststrandbias} | {params.var2vcf_valid} -N {wildcards.patient} -E -f 0.001 > {output}
        """

#create lineplots of the allele frequencies of the variants across the decreasing concentrations (yf sets the y axis max to 0.07 for easier comparison of plots)
rule lineplotting:
    input:
        expand("../call_vcfs/{{trend}}_{{patient}}/{{trend}}_{{patient}}_{timepoint}/{{varcaller}}.vcf", timepoint = timepoints)
    output:
        "../lineplots/{trend}_{patient}/{varcaller}.png"
    params:
        plot_AF = config['plot_AF'],
        ymax = config['plot_ymax']
    shell:
        """
        python {params.plot_AF} -i {input} -o {output} -s before th1 th2 th3 after -yf {params.ymax} -l
        """