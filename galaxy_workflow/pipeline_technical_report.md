Technical Report: RNA-seq Read Processing & Quantification Pipeline
Overview 

This repository demonstrates an alignment of rna-seq reads to a reference genome through bioinformatic tools on the galaxy platform. 
To produce gene counts from NGS data. 


Raw sequencing reads (FASTQ format) were accessed from a Zenodo data repository, for sequence reads obtained from biological replicates of Basal and Luminal cells 
The raw data is processed into a read count matrix which can be passed downstream to differential expression analysis in the pipeline.



Input 
raw FASTQ reads organism, mouse mammary gland, cell type - basal vs luminal cells (e.g. MCL1-DL, MCL1-DH).
Output 
count matrix 



step 1. Data acquisition & Integrity 
Data Source: Zenodo repositories of sequence reads from Basal and Luminal cells (urls in sample.tsv)
Data management: Galaxy Dataset collections used to group and manage biological replicates and ensure parallel processing 


step 2. Quality Control 
Tools: FastQC and MultiQC 
Objective: Need to ensure the reads meet quality standards to ensure any biological interpretations from the output data is valid 
Quality criteria: Per base sequence quality (Phred scores), GC Content distribution, and adapter contamination levels. 
Outcome: Identified read degradation at 3'ends and potential adapter read through, to remove adapter sequences from reads, trimming step is performed 


step 3. Read Trimming & Artifact Removal 
Tool: Cutadapt 
Parameter Configuration:
	Quality Cutoff: 20 (Phred score) Bases must reach this read quality to ensure high confidence mapping of reads to reference genome and gene count 
	Illumina Universal Adapters were scanned and removed 
	Minimum length must be at least 20bp to limit multimapping ambiguity 

step 4. Splice Aware Genomic Alignment 
Tool: HISAT2 (Hierarchical Indexing for Spliced Alignment of Transcripts)
Rationale: HISAT2 will account for introns only present in reference genome but not reads, accurately mapping across exon-exon junctions to prevent overlapping genes which would occur if regular DNA aligners. e.g. Bowtie2 was used. 
• Reference Genome: Mus musculus (mm10/GRCm38).
5 Gene Quantification 
Tool - FeatureCounts 
Method - Mapped reads (BAM files) are overlapped with reference genome annotation (GTF) to count quantity of reads that mapped to exons. 
Configuration 
Strandedness - Unstranded (choice during library prep)
Summarisation: Multi-mapping reads excluded to ensure quantification accuracy

Output 
Raw count matrix (Genes X Samples)
Downstream Statistical Considerations
Raw count matrix is not normally distributed, and influenced by sequencing depth. 
Need to prepare data before differential expression analysis. 
1. Filtering - Genes with low levels of expression (CPM < 0.5) are removed to reduce noise 
2. Normalisation: Trimmed Mean of M-values (TMM). we account for the scaling factor, the library size difference and compositional bias that exists between samples, so comparisons done through Limma-voom or edgeR are reflective of biological differences and not technical or random variance between cells 
