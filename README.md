# High-Confidence Variant Calling Pipeline

A robust, reproducible, and scalable Snakemake pipeline designed to identify high-confidence single nucleotide polymorphisms (SNPs) and short insertions/deletions (indels) from short-read sequencing data.

To drastically reduce false-positive rates, this pipeline utilizes an orthogonal calling strategy. It runs both local de-novo assembly (GATK HaplotypeCaller) and Bayesian haplotype-based calling (Freebayes), intersecting the results via a custom Python implementation to produce a highly accurate consensus variant set.

🧬 **Pipeline Architecture**

```mermaid
graph TD
    %% Define visual styles for different file types and processes
    classDef input fill:#e1f5fe,stroke:#0288d1,stroke-width:2px,color:#000
    classDef reference fill:#e8f5e9,stroke:#388e3c,stroke-width:2px,color:#000
    classDef bam fill:#fff3e0,stroke:#f57c00,stroke-width:2px,color:#000
    classDef vcf fill:#f3e5f5,stroke:#7b1fa2,stroke-width:2px,color:#000
    classDef python fill:#fffde7,stroke:#fbc02d,stroke-width:2px,color:#000
    classDef qc fill:#eceff1,stroke:#546e7a,stroke-width:2px,color:#000
    classDef process fill:#ffffff,stroke:#9e9e9e,stroke-width:1px,color:#000

    %% Inputs
    R1([Raw FASTQ R1]):::input
    R2([Raw FASTQ R2]):::input
    REF[(Reference Genome FASTA)]:::reference

    %% Step 0
    subgraph Step 0: Reference Indexing
        direction LR
        REF --> BWA_IDX(bwa index):::process
        REF --> FAI(samtools faidx):::process
        REF --> DICT(gatk CreateSequenceDictionary):::process
    end

    %% Step 1
    subgraph Step 1: Alignment & Preprocessing
        R1 --> BWA(bwa mem):::process
        R2 --> BWA
        BWA_IDX -.-> BWA
        BWA --> |SAM| SORT(samtools sort):::process
        SORT --> |BAM| DEDUP(picard MarkDuplicates):::process
        DEDUP --> |Dedup BAM| BAM_OUT[(Analysis-Ready BAM)]:::bam
    end

    %% Step 2 & 3
    subgraph Step 2 & 3: Orthogonal Variant Calling
        BAM_OUT --> FB(Freebayes):::process
        BAM_OUT --> GATK(GATK HaplotypeCaller):::process
        REF -.-> FB
        REF -.-> GATK
        FB --> |VCF| BGZ1(bgzip & tabix):::process
        GATK --> |VCF| BGZ2(bgzip & tabix):::process
        BGZ1 --> FB_VCF[(Freebayes VCF.gz)]:::vcf
        BGZ2 --> GATK_VCF[(GATK VCF.gz)]:::vcf
    end

    %% Step 4
    subgraph Step 4: Consensus Generation
        FB_VCF --> CONSENSUS{Pure Python Script<br/>Strict Set Intersection}:::python
        GATK_VCF --> CONSENSUS
        CONSENSUS --> |Consensus VCF| CONS_VCF[(High-Confidence VCF)]:::vcf
    end

    %% Step 5
    subgraph Step 5: Functional Annotation
        CONS_VCF --> SNPEFF(SnpEff):::process
        SNPEFF --> ANN_VCF[(Annotated VCF)]:::vcf
    end

    %% Step 6
    subgraph Step 6: Quality Control
        DEDUP -.-> |Duplicate Metrics| MULTIQC(MultiQC):::process
        SNPEFF -.-> |Annotation Stats| MULTIQC
        MULTIQC --> REPORT([Interactive HTML Report]):::qc
    end
```


**Features**

Fully Automated: Managed end-to-end by Snakemake.

Reproducible Environments: Tool dependencies are strictly isolated using Conda/Mamba environments (envs/).

Orthogonal Validation: Cross-references GATK and Freebayes to eliminate caller-specific artifacts.

Robust Merging: Uses a native Python set-intersection script to bypass common bcftools merge segmentation faults caused by differing VCF header types.

Comprehensive QC: Aggregates mapping, duplication, and annotation metrics into a single interactive MultiQC report.

**Prerequisites**

You will need Conda (preferably mamba for much faster environment solving) and Snakemake installed on your system.

**Install Mamba (if not already installed)**
```conda install -n base -c conda-forge mamba```

**Install Snakemake**
```mamba create -c conda-forge -c bioconda -n snakemake snakemake```


**Usage**

1. Clone the Repository & Setup Directory

```git clone https://github.com/YourUsername/VariantCallingPipeline.git```

```cd VariantCallingPipeline```


Ensure your directory structure looks like this:

```
VariantCallingPipeline/
├── Snakefile
├── envs/
│   ├── mapping.yaml
│   ├── calling.yaml
│   └── qc.yaml
├── data/               <-- Place your raw .fastq.gz files here
└── reference/          <-- Place your genome.fa file here
```

2. Configure the Pipeline

Open the Snakefile in a text editor and update the inputs at the top of the file:

**Update with your actual sample names**

```SAMPLES = ["Patient1", "Patient2", "Control1"]``` 

**Point to your reference genome**
```REF = "reference/hg38.fa"```

**Specify the SnpEff database for your organism**
```SNPEFF_DB = "GRCh38.105"``` 


Note: Make sure to download your SnpEff database before running the pipeline on real data:

```conda activate qc -> snpEff download GRCh38.105```

3. Execution

Activate your Snakemake environment and run the pipeline.

Dry Run (Highly Recommended):
This will build the DAG and show you exactly what will be executed without actually running anything.

``conda activate snakemake
snakemake -n``


**Full Production Run:**
Run the pipeline, allowing it to use all available cores. We recommend running this inside a tmux or screen session on your Linux server.

```snakemake --use-conda --conda-frontend mamba --cores all```


**Key Output**

Upon successful completion, the results/ directory will contain:

```results/{sample}_consensus/consensus.vcf```: The raw, unannotated high-confidence variants agreed upon by both callers.

```results/{sample}.annotated.vcf```: The consensus VCF enriched with biological impact predictions (missense, frameshift, etc.).

```results/multiqc_report.html```: A beautiful, interactive dashboard summarizing sequence quality, duplication rates, and variant impacts across all your samples.
