# ngs-variant-calling-mini-pipeline

**Clinical-Style Variant Discovery and Quality Control Mini-Pipeline**

## Problem
Align reads to a reference genome, and visualise the alignment 
Variant calling to identify indels and SNPs 

## Methods
 * Reference Indexing with BWA
 * Read Alignment with BWA-MEM
 * Convert to BAM Format with Samtools
 * Sort BAM File with Samtools
 * Index BAM File with Samtools
 * Variant Discovery with Freebayes
 * Visualisation of VCF and indexed BAM file with IGV 


## Results
*   **Output Table (VCF):** A Variant Call Format file containing all detected and filtered SNPs and Indels.
    *   `results/tables/variants_filtered.vcf`
*   **Figures (Visualization Prep):** The resulting alignment file and index, which are the inputs for IGV visualization.
    *   `results/figures/ind1_mapped_sorted.bam`
    *   `results/figures/ind1_mapped_sorted.bam.bai`




## Reproducibility
   *To reproduce this analysis, clone the repository and set up the necessary       tools using Conda:

      ```bash
      # Create and activate the Conda environment
      conda env create -f environment.yml
      conda activate ngs-pipeline

      # Example: Run the full pipeline script
      # bash scripts/run_pipeline.sh

      # Open IGV for visual validation
      # igv
      environment.yml
      name: ngs-pipeline
      channels:
        - bioconda
        - conda-forge
        - defaults
      dependencies:
        - bwa=0.7.17
        - samtools=1.16.1
        - freebayes=1.3.6
        - igv=2.17.0
        - python=3.10

## Skills demonstrated 
   1. NGS Pipeline Engineering
      joining reference genome indexing, paired-ends reads alignment,                   visualisation
      Proficiency in designing, implementing, and managing a robust                   bioinformatics pipeline
      using fundamental command-line tools (Bash/Linux).
      This includes sequential data processing steps, file format manipulation          (SAM to indexed BAM),
      and pipeline optimization using multi-threading (e.g., -t 2 and -@ 2 in          BWA-MEM and Samtools).

   2. Clinical-Grade Variant Discovery & Filtering
      Expertise in using industry-standard tools like BWA-MEM for read alignment
      and Freebayes for accurate SNP/Indel calling.
      This skill specifically includes applying stringent quality control             metrics,
      such as enforcing a minimum Mapping Quality (MQ) of 30,
      to ensure only high-confidence variants are reported.

   3. Data Visualization and Quality Control (QC)
      Ability to prepare and utilize specialized genomic visualization tools          (IGV) for critical manual quality          assurance. This demonstrates          the skill to interpret read-level data, validate computational results,         confirm Allelic Frequency, and effectively rule out technical artifacts          or sequencing errors.



