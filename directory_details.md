## Snakemake Project Directory Structure and Content

Following the official Snakemake best-practices and standard workflow recommendations, a genomic variant-calling pipeline (such as the alignment, orthogonal calling, and consensus workflow outlined previously) is organized into a clean, modular directory structure.

---

### **Project Directory Tree**

```text
my_snakemake_project/
├── Snakefile                 # Main entry point directing the workflow
├── config/
│   ├── config.yaml           # Global parameters, paths, and tool settings
│   └── samples.tsv           # Tabular metadata/sample sheet mapping R1/R2 reads
├── workflow/
│   ├── Snakefile             # Core workflow logic (if decoupled from root)
│   ├── rules/
│   │   ├── common.smk        # Helper functions, input selectors, and path builders
│   │   ├── indexing.smk      # Step 0: Reference indexing rules (bwa, samtools, gatk)
│   │   ├── alignment.smk     # Step 1: Alignment, sorting, and duplicate marking
│   │   ├── calling.smk       # Steps 2 & 3: Orthogonal variant callers (Freebayes, GATK)
│   │   ├── consensus.smk     # Step 4: Python-based set intersection consensus script
│   │   ├── annotation.smk    # Step 5: Functional annotation via SnpEff
│   │   └── qc.smk            # Step 6: MultiQC aggregation and reporting
│   ├── scripts/
│   │   └── consensus_caller.py # Pure Python script for strict VCF set intersection
│   ├── envs/
│   │   ├── alignment.yaml    # Conda environment for bwa, samtools, picard
│   │   ├── calling.yaml      # Conda environment for freebayes and gatk
│   │   └── annotation.yaml   # Conda environment for snpeff and multiqc
│   └── report/
│       └── quality.rst       # Text/RST caption templates for interactive HTML reports
├── resources/                # Static or raw input data directory
│   ├── ref/                  # Reference genome FASTA and index files
│   └── fastq/                # Raw paired-end sequencing inputs (R1, R2)
└── results/                  # Generated workflow output directory
    ├── bams/                 # Analysis-ready BAM files and indexes
    ├── vcfs/                 # Intermediate and high-confidence consensus VCFs
    ├── annotated/            # Final functionally annotated VCF outputs
    └── qc/                   # MultiQC HTML report and individual metrics

```

---

### **Folder and File Breakdown**

#### **1. Root Directory (`my_snakemake_project/`)**

* **`Snakefile`**: The central master file invoked by the Snakemake execution engine. It typically sets configuration paths and includes modular rule files from `workflow/rules/`.

#### **2. Configuration Directory (`config/`)**

* **`config.yaml`**: Contains global pipeline settings, resource thresholds, tool flags, and references to genome builds.
* **`samples.tsv`**: A tabular metadata file mapping sample identifiers to their respective raw FASTQ file locations (`R1`/`R2`).

#### **3. Workflow Core Directory (`workflow/`)**

* **`rules/`**: Houses modular rule files ending in `.smk`. Each file governs a logical block of the pipeline (indexing, alignment, variant calling, consensus generation, annotation, and QC).
* **`scripts/`**: Contains external execution scripts (e.g., custom Python scripts handling strict set intersections for consensus VCF generation).
* **`envs/`**: Contains isolated Conda YAML definition files specifying exact software dependencies and versions for individual rules to ensure reproducibility.
* **`report/`**: Contains text or restructured text (`.rst`) files providing captions, descriptions, and metadata for elements featured in Snakemake’s automatic interactive HTML reports.

#### **4. Data Directories (`resources/` & `results/`)**

* **`resources/`**: Houses upstream input data that are read-only to the pipeline, including raw FASTQ reads and reference genome FASTA files.
* **`results/`**: The designated output directory where Snakemake deposits all intermediate indices, aligned BAMs, VCF files, annotated outputs, and MultiQC HTML summary dashboards.