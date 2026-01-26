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


Quality Control Configuration:
• FastQC: Executed using default parameters to assess per-base sequence quality and detect standard Illumina adapter contamination.
• MultiQC: Configured to aggregate "Raw Data" logs from the FastQC collection to generate a unified quality report.
	Type of FastQC output? Raw data
	FastQC output = FastQC on collection 1: RawData


	FastQC Configuration Rationale:
• Contaminant & Adapter Lists (Left Empty): Relied on FastQC’s internal standard library of common adapters (e.g., Illumina Universal, Nextera). Manually supplying a list was unnecessary as standard Illumina adapters were expected and automatically detected.
• Limits & Submodules (Default): Retained default processing limits to ensure an unbiased assessment of the raw data distribution without artificially filtering out low-quality reads at the visualization stage.
• Kmer Length (Default 7): Maintained the default Kmer size of 7. This length provides an optimal balance between sensitivity and specificity for detecting overrepresented repetitive motifs and biases in RNA-seq data without introducing excessive noise.
• Workflow Logic: Crucially, no filtering or trimming was applied within FastQC. This step was strictly reserved for diagnosis, ensuring that the subsequent +adapt step addressed identified issues (adapters/quality) based on the raw, unaltered reports.



step 3. Read Trimming & Artifact Removal 
Tool: Cutadapt 
Parameter Configuration:
	Quality Cutoff: 20 (Phred score) Bases must reach this read quality to ensure high confidence mapping of reads to reference genome and gene count 
	Illumina Universal Adapters were scanned and removed 
	Minimum length must be at least 20bp to limit multimapping ambiguity 

Trimming Strategy & Rationale:
• Adapter Selection (3' End): Configured Cutadapt to remove the Illumina Universal Adapter (AGATCGGAAGAG) from the 3' end of reads. This addresses "adapter read-through," a common phenomenon where the sequencing instrument reads beyond the biological insert into the 3' adapter sequence when the RNA fragment is shorter than the read length.
• Quality Filtering (Phred > 20): Applied a quality cutoff of 20 (99% base call accuracy). This algorithmically trims bases from the 3' end until the quality score rises above 20, preventing low-confidence base calls from interfering with the alignment algorithm (HISAT2).
• Length Filtering (> 20bp): Implemented a minimum length filter of 20bp. Reads reduced to <20bp after trimming were discarded, as extremely short sequences lack the complexity required for unique genomic mapping and would introduce noise into the count matrix.


Adapter Handling & Error Tolerance:
• Trimming Behavior: Configured the algorithm to Trim detected adapters. This removes the non-biological adapter sequence and any downstream bases, which is critical for accurate mapping; failing to remove these would lead to unmapped reads or mismatches at the 3' end.
• Error Rate (0.1): Utilized a maximum error rate of 10%. This parameter makes the pipeline robust against sequencing errors, which typically increase toward the 3' end of reads (Phred score decay). A stricter rate might miss adapters containing sequencing errors, leading to contamination.
• Minimum Overlap (3bp): Enforced a minimum overlap of 3 bases. This threshold balances sensitivity and specificity; it allows the detection of partial adapters (where the read ends inside the adapter) while preventing false-positive trimming caused by random sequence matches of 1 or 2 nucleotides.








• What to do if a match is found: Select Trim.
    - Will allow for the adapter sequence and everything following it to be removed from the read.
• Maximum error rate: Set to 0.1 (Default).
    - This allows for a 10% error rate (1 mismatch per 10 bases), which accounts for sequencing errors that are common at the ends of reads.
• Do not allow indels: Select No (Default).
    ◦ Why: Sequencing isn't perfect; small insertions or deletions might occur in the adapter sequence. You still want to identify and trim it.
• Match times: Set to 1 (Default).
    ◦ Why: Illumina adapters usually appear only once at the 3' end of the insert.
• Minimum overlap length: Set to 3 (Default).
    ◦ Why: Requiring at least 3 bases to match prevents the tool from randomly trimming the read just because the last base happens to match the first base of the adapter (false positives).
• Match wildcards in reads: No (Default).
• Match wildcards in adapters: Yes (Default).
• Look for adapters in the reverse complement: No (Default).
    ◦ Why: In standard Illumina libraries, the adapter orientation is fixed relative to the read.






Quality Trimming Configuration:
• Quality Cutoff (Phred 20): A Phred score threshold of 20 (corresponding to a 1% error probability) was applied to the 3' end of reads. This removes low-confidence base calls that typically accumulate at the end of sequencing reads, preventing sequencing errors from negatively impacting the alignment rate in HISAT2.
• Fixed Trimming (Disabled): "Bases to cut" was maintained at 0 to preserve maximum read length, as FastQC 	did not indicate significant positional bias or artifacts at the 5' end requiring unconditional removal.


	"Bases to cut from R1..." left at the default of 0.
	: This setting blindly removes bases from the start (5' end) of every read. Unless your FastQC report 		showed a specific problem at the very beginning of the reads (like severe nucleotide bias), you should not 	remove these bases. Keeping them maximizes the sequence length available for the aligner to map the read 	accurately.



	
• Platform Specifics: NextSeq trimming was disabled as the dataset was processed using standard Illumina quality encoding, not the 2-color chemistry specific to NextSeq platforms.






--------------------------------------------------------------------------------

Additional Trimming Parameters:
• Poly-A Trimming (Disabled): Poly-A tail trimming was disabled. While the library preparation involved Poly-A selection, aggressive algorithmic removal of Poly-A tails is unnecessary when standard adapter trimming is sufficient to handle read-through events.
• Fixed Length Trimming (Disabled): Reads were not shortened to a fixed length. Unlike some genomic DNA workflows that require uniform lengths, RNA-seq relies on maximizing the sequence overlap with reference exons and splice junctions. Hard-clipping reads to a specific length would unnecessarily discard high-quality bases, reducing the probability of unique mapping in HISAT2.









• Discard Trimmed Reads: No (Default).
• Discard Untrimmed Reads: No (Default).
• Minimum length (R1): 20.
    ◦ Note: This is a critical setting mentioned in your workshop notes.
• Minimum length (R2): (None).
    ◦ i.e. use the same 20bp threshold as R1.
• Maximum length (R1): (None /Default).

Read Filtering Strategy:
• Retention Policy: Both "Trimmed" and "Untrimmed" reads were retained. "Untrimmed" reads simply represent biological fragments longer than the sequencing read length (no adapter read-through), which are high-quality data points. Discarding them would result in significant data loss.
• Minimum Length Filter (20bp): A minimum length threshold of 20bp was enforced. Reads shorter than this after trimming lack sufficient sequence complexity to map uniquely to the reference genome (mm10). Short reads significantly increase the risk of multi-mapping (mapping to multiple genomic loci) and false-positive alignments, introducing noise into the quantification step.








 Maximum length (R2): (Leave Empty).
• Max N: (Leave Empty).
• Max expected errors: (Leave Empty).
• Max average expected errors: (Leave Empty) (Default).
• Discard CASAVA-filtered reads: No (Default).
• Pair filter: Any (Default).
    ◦ Note: This is critical. It ensures that if either the Forward (R1) or Reverse (R2) read fails a quality check (e.g., becomes <20bp), the entire pair is discarded. This keeps the Forward and Reverse files in sync, which is required for the HISAT2 aligner.





Paired-End Filtering Logic:
• Pair Synchronization (Filter 'Any'): The pipeline was configured to discard the entire read pair if either mate failed the quality or length thresholds (Pair filter: 'Any'). This prevents "broken pairs" (singleton reads) which would cause synchronization errors during the alignment step (HISAT2) and downstream counting.
• Ambiguity Handling: Explicit "N" filtering was disabled in favor of the Phred score quality cutoff. Since "N" bases typically have a quality score of 0, the previously applied quality cutoff of 20 effectively handles ambiguous bases without requiring a redundant filter.
• Length Preservation: No maximum length cutoff was applied. Retaining maximum sequence length is advantageous for splice-aware alignment, as longer reads provide more anchor points for mapping across exon-exon junctions.





 Strip suffix: (Leave Empty).
• Length tag: (Leave Empty).
• Rename reads: (Leave Empty).
• Change negative quality values to zero: No (Default).

Read Modification & Header Integrity:
• Header Preservation: Read modification options (suffix stripping, renaming, and length tagging) were disabled. Retaining the original FASTQ headers is essential for maintaining the traceability of biological replicates and ensuring that downstream tools (HISAT2) can accurately identify read pair information based on standard Illumina CASAVA headers.
• Quality Encoding: The "Change negative quality values to zero" option was disabled. The dataset utilizes standard Phred+33 encoding (Sanger/Illumina 1.8+), which does not generate negative quality scores. This option is reserved for resolving offsets in older legacy formats (e.g., Solexa) and is unnecessary for modern NGS data.






Trimming Metrics & Reporting:
• Report Generation: The "Cutadapt's per-adapter statistics" report was explicitly generated. This textual log captures essential metrics—including total reads processed, number of reads with adapters, and base pairs removed—which serves as the mandatory input for MultiQC.
• Pipeline Integration: Generating this report allows for a comparative "Pre-trim vs. Post-trim" quality assessment, verifying that the adapter contamination identified in the initial FastQC stage was successfully ameliorated without excessive data loss.

















	
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
