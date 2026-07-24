# ==============================================================================
# VARIANT CALLING PIPELINE (PRODUCTION)
# ==============================================================================
# This Snakemake pipeline performs alignment, duplicate marking, variant calling 
# (GATK & Freebayes), consensus generation, and annotation.
# ==============================================================================

# --- Configuration & Inputs ---
# In production, replace this hardcoded list with a dynamic glob.
# Example: SAMPLES, = glob_wildcards("data/{sample}_R1.fastq.gz")
SAMPLES = ["Sample1", "Sample2"] 

REF = "reference/genome.fa"
SNPEFF_DB = "GRCh38.105" # Replace with your actual target organism database

rule all:
    input:
        expand("results/{sample}_consensus/consensus.vcf", sample=SAMPLES),
        expand("results/{sample}.annotated.vcf", sample=SAMPLES),
        "results/multiqc_report.html"

# --- 0. Reference Indexing ---
rule index_reference:
    """
    Creates necessary indices for the reference genome.
    Required by BWA, GATK, and Samtools.
    """
    input:
        REF
    output:
        "reference/genome.dict",
        "reference/genome.fa.fai",
        "reference/genome.fa.bwt",
        "reference/genome.fa.amb",
        "reference/genome.fa.ann",
        "reference/genome.fa.pac",
        "reference/genome.fa.sa"
    conda:
        "envs/mapping.yaml"  
    shell:
        """
        samtools faidx {input}
        gatk CreateSequenceDictionary -R {input}
        bwa index {input}
        """

# --- 1. Alignment & Preprocessing ---
rule bwa_map:
    """
    Maps paired-end reads to the reference genome.
    Injects Read Group (@RG) headers which are strictly required by GATK.
    """
    input:
        ref = REF,
        r1 = "data/{sample}_R1.fastq.gz",
        r2 = "data/{sample}_R2.fastq.gz"
    output:
        temp("results/{sample}.sam")
    conda:
        "envs/mapping.yaml"
    threads: 8
    shell:
        "bwa mem -t {threads} -R '@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\\tPL:ILLUMINA' {input.ref} {input.r1} {input.r2} > {output}"

rule samtools_sort:
    """
    Converts SAM to BAM and sorts by genomic coordinate.
    """
    input:
        "results/{sample}.sam"
    output:
        "results/{sample}.sorted.bam"
    conda:
        "envs/mapping.yaml"
    threads: 4
    shell:
        "samtools sort -@ {threads} -o {output} {input}"

rule mark_duplicates:
    """
    Identifies and marks PCR duplicates to prevent amplification bias.
    """
    input:
        "results/{sample}.sorted.bam"
    output:
        bam = "results/{sample}.dedup.bam",
        bai = "results/{sample}.dedup.bam.bai",
        metrics = "results/{sample}.dup_metrics.txt"
    conda:
        "envs/mapping.yaml"
    shell:
        """
        picard MarkDuplicates \
            I={input} \
            O={output.bam} \
            M={output.metrics} \
            REMOVE_DUPLICATES=true \
            VALIDATION_STRINGENCY=LENIENT
        
        samtools index {output.bam} {output.bai}
        """

# --- 2. Variant Calling Paths ---
rule freebayes_call:
    """
    Calls variants using Freebayes (haplotype-based caller).
    """
    input:
        bam = "results/{sample}.dedup.bam",
        bai = "results/{sample}.dedup.bam.bai",
        ref = REF
    output:
        temp("results/{sample}.freebayes.vcf")
    conda:
        "envs/calling.yaml"
    shell:
        "freebayes -f {input.ref} {input.bam} > {output}"

rule gatk_haplotype_caller:
    """
    Calls variants using GATK HaplotypeCaller (local de-novo assembly).
    """
    input:
        bam = "results/{sample}.dedup.bam",
        bai = "results/{sample}.dedup.bam.bai",
        ref = REF,
        ref_dict = "reference/genome.dict",
        ref_idx = "reference/genome.fa.fai"
    output:
        temp("results/{sample}.gatk.vcf")
    conda:
        "envs/calling.yaml"
    shell:
        """
        gatk HaplotypeCaller \
            -R {input.ref} \
            -I {input.bam} \
            -O {output}
        """

# --- 3. Compression & Indexing ---
rule bgzip_and_index:
    """
    Compresses and indexes VCF files for efficient programmatic reading.
    """
    input:
        "results/{sample}.{caller}.vcf"
    output:
        vcf = temp("results/{sample}.{caller}.vcf.gz"),
        tbi = temp("results/{sample}.{caller}.vcf.gz.tbi")
    conda:
        "envs/calling.yaml"
    shell:
        """
        bgzip -c {input} > {output.vcf}
        tabix -p vcf {output.vcf}
        """

# --- 4. Consensus Calling (Python Implementation) ---
rule consensus_call:
    """
    Finds the intersection of variants called by both Freebayes and GATK.
    Implemented in pure Python to bypass strict VCF header compatibility issues
    often encountered when merging different callers with bcftools.
    """
    input:
        fb = "results/{sample}.freebayes.vcf.gz",
        gatk = "results/{sample}.gatk.vcf.gz"
    output:
        vcf = "results/{sample}_consensus/consensus.vcf"
    params:
        outdir = "results/{sample}_consensus"
    run:
        import gzip
        import os
        
        os.makedirs(params.outdir, exist_ok=True)
        
        # Parse Freebayes variants into a set of tuples: (CHROM, POS, REF, ALT)
        fb_variants = set()
        with gzip.open(input.fb, 'rt') as f:
            for line in f:
                if not line.startswith("#"):
                    parts = line.split()
                    fb_variants.add((parts[0], parts[1], parts[3], parts[4]))
        
        # Parse GATK variants
        gatk_variants = set()
        with gzip.open(input.gatk, 'rt') as f:
            for line in f:
                if not line.startswith("#"):
                    parts = line.split()
                    gatk_variants.add((parts[0], parts[1], parts[3], parts[4]))
        
        # Find pure intersection
        common = fb_variants.intersection(gatk_variants)
        
        # Write clean consensus VCF
        with open(output.vcf, "w") as f:
            f.write("##fileformat=VCFv4.2\n")
            f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            # Sort variants by chromosome, then position
            for var in sorted(list(common), key=lambda x: (x[0], int(x[1]))):
                f.write(f"{var[0]}\t{var[1]}\t.\t{var[2]}\t{var[3]}\t.\tPASS\t.\n")

# --- 5. Functional Annotation ---
rule annotate_variants:
    """
    Annotates the consensus variants with biological impact (e.g., missense, synonymous)
    using SnpEff.
    """
    input:
        vcf = "results/{sample}_consensus/consensus.vcf"
    output:
        vcf = "results/{sample}.annotated.vcf",
        stats = "results/{sample}.snpeff_summary.html"
    params:
        db = SNPEFF_DB
    conda:
        "envs/qc.yaml"
    shell:
        """
        # Ensure database is downloaded locally or available in your snpEff config
        snpEff -Xmx8g {params.db} {input.vcf} -s {output.stats} > {output.vcf}
        """

# --- 6. Quality Control ---
rule multiqc:
    """
    Aggregates all log files (Picard, SnpEff, etc.) into a single interactive HTML report.
    """
    input:
        dup_metrics = expand("results/{sample}.dup_metrics.txt", sample=SAMPLES),
        snpeff_stats = expand("results/{sample}.snpeff_summary.html", sample=SAMPLES)
    output:
        report = "results/multiqc_report.html"
    params:
        scan_dir = "results/",
        out_dir = "results/"
    conda:
        "envs/qc.yaml"
    shell:
        "multiqc {params.scan_dir} -o {params.out_dir} -n multiqc_report.html --force"