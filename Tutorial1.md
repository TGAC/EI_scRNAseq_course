#Single Cell RNA-seq using Galaxy

##Single Cell RNA-SeqSingle Cell mRNA sequencing (scRNA-seq) offers the ability to discover new genes and transcripts and measure transcript expression. The main difference between bulk RNA-seq and single cell RNA-seq is that each library originates from a single cell, instead of a population of cells. Because of this, the sequencing libraries can suffer from low read counts, non-uniform PCR amplification and a number of other confounding factor. Therefore, we pay significant attention to the quality control of each sequenced cell to remove low quality data from our analyses.

In this tutorial, you will start by walking through the individual steps of taking one scRNA-seq Illumina library, mapping the reads to the Arabidopsis genome and calculating expression for each annotated gene. 
After this, you'll use a ready-made Galaxy workflow to analyse multiple sequencing libraries and create an expression matrix for them. You'll go on to learn how to visualise your data to identify low-quality data and filter the data so that it's ready for analysis.  
You'll also be given a couple of datasets to 'troubleshoot', where you must decide what, if any, problems have occured with the generation of the data.
By the end of the tutorial you'll understand the complexities involved in scRNA-seq and how to examine and filter your data to get analysis-ready high-quality data.

Each step is numbered (e.g 'Step 1:, Step 2:, etc.). Once you reach the end of a step, stop working and place your green post-it note on the back of your laptop so we can see when everyone has finished. When you start the next step, take your post-it note down.  ##Software for scRNA-seqIn today's tutorial you will use a number of software applications, using Galaxy as the interface. Note: This tutorial is best executed using the Google Chrome web-browser.###STAR alignerMapping mRNA (which contains a string of concatenated exon sequences) to a genome (which contains exons and introns) means that reads from the RNA-seq experiment may span exons. Therefore, when RNA-seq reads are mapped to the genome, they may need to be split so that each end of the read maps to an exon with an intron in between them.   STAR is a fast splice junction mapper, which means that as well as aligning reads to a genome, it also accounts for intron-exon boundaries and splits reads apart in an attempt to map reads that bridge these intron-exon boundaries.
   

###Subread FeatureCounts
FeatureCounts is a software program developed for counting reads aligned to genomic features such as genes, exons, promoters, or any other pre-defined genomic locations.

###Galaxy Text Manipulation
Many of the tools in Galaxy's "Text Manipulation" tool group can be used to manipulate data in order to bring it into alignment with a standarised specification. You will use this tool to combine gene-read-counts for multiple samples and create a scRNA-seq expression matrix.##Exercise 1 - A single sample walkthroughYou will start by analysing just one sample of scRNA-seq data. Make sure you understand the process involved and the type of data used and generated. Upscaling this knowledge into a multiple-sample experiment will be quite straightforward once you understand how one sample is processed.##Getting started
In your history panel, click on the '+' sign to create a new history. A new, empty history will be created. Now, click on 'Unnamed history' and name it 'Single Cell 1' and hit return.##Step 1: Getting the data1. You will start by obtaining your data in Galaxy. Click on 'Shared Data' (in the panel at the top of the page) and then click on 'Data Libraries'1. In the list of Data Libraries you will see one called 'Galaxy courses'. Click on this and you will see a couple of sub folders. You want the 'Single Cell 2020' folder. Click on this and you'll see a number of folders. Select the 'Single Cell 1' folder. 
2. You should see 3 files: `A1_R1.fastq.gz` and `A1_R2.fastq.gz` which are Illumina paired-end sequences from a scRNA-seq experiment and `Arabidopsis_thaliana.TAIR10.37.gtf` which is an annotation.1. Select all the files and then click on 'Export to History' on the menu at the top of the page and then 'as Datasets' in the dropdown box. Make sure the 'Single Cell 1' history is selected in the 'Select history' box below and click on 'Import'. Click on 'Analyze Data' at the the top of the page to return to the main interface.1. The three items of data should now be in your history.

Take a look at one of the `fastq.gz` files, but clicking on the 'eye' icon to the right of the history item name. 

Each of the paired-end sequences are stored in different files, denoted by 'R1' or 'R2' in the file name. R1 are the forward-facing of the two paired-end reads and R2 are the reverse-facing of the pairs. The Illumina data consists of paired-end 101bp Illumina sequences and is in FASTQ format, a four-line format.  

 *  Line 1 is the read ID and always starts with an `@` sign. Whether the read is forward or reverse is denoted by `1:N:0` (forward) or `2:N:0` (revers) towards the end of the ID. Any letters after this are barcodes used in demultiplexing. 
*  Line 2 is the actual nucleotide sequence of the read. 
*  Line 3 is fairly redundant and often only contains a `+` sign. Otherwise, it would contain the same information as the read ID. 
*  Line 4 is the ASCII-formatted quality scores, as explained in the introductory presentation.

Stop here. Remember to let the instructor know you've finished this step.
##Step 2: Mapping with STAR###MappingThe term 'mapping' refers to assigning short stretches of DNA (like those found in Illumina reads) to identical (or near-identical) positions on a reference genome. Now you have your reads, you can go ahead and map your scRNA-seq reads to the Arabidopsis reference with the splice-aware-aligner STAR. You will start off by going through the steps of analysing one library of scRNA-seq reads. Of course, normally, you would have tens, or hundreds (or even thousands!) of libraries.1.  At the top of the Tools panel, type 'STAR' in the 'search tools' box. From the search results, select '**RNA STAR** Gapped-read mapper for RNA-seq data'. The interface for the tool should appear in the main window.
2. Next, click on 'Versions' at the top of the tool interface and select version 2.7.2b1. There are a few options that you need to set in STAR. 1. Starting at the top of the page, make sure that your 'Single-end or paired-end reads' is selected as `Paired-end (as individual datasets)`.
1. You need to supply the correct forward and reverse reads to the program, so under 'RNA-Seq FASTQ/FASTA file, forward reads' select the `A1_R1.fastq.gz` file and `A1_R2.fastq.gz` for the 'reverse reads' file.
1. Under 'Custom or built-in reference genome' make sure 'Use built-in index' is selected and two boxes down select 'Arabidopsis thaliana (TAIR10)' in the 'Select reference genome' box. For the box above entitled 'Reference genome with or without an annotation', make sure 'use genome reference **without** builtin gene-models'. Selecting these options means that the STAR aligner will align your reads to the reference without guidance from any annotation. Doing this will avoid reads being 'forced' to align to transcripts when they might better align to other regions in the genome.
1. The remaining parameters can be left as default, which is the general rule-of-thumb unless you have a specific reason to change any of them.Click 'Execute'.   Three new history items will appear in the history panel. They will be grey at first and then turn yellow to say that the job is executing. When the job has executed successfully the history items will turn green. If they turn red, seek help. This step should take a few minutes to complete.  
Once finished, you will see three new history items (output files), sequentially numbered, all starting with something like `4. RNA STAR on data 2 and data 1:`. There should be a `log` file, a `splice junctions.bed` file and a `mapped.bam` file.  
1. Click on the eye icon of the `log` file to view the contents. You should see the mapping stats for your alignment file (`mapped.bam` which you'll look at later). It's worth giving this file a cursory look over. Some important statistics to look at are 'Number of input reads', 'Uniquely mapped reads %', 'Mismatch rate per base %' and the 'Multi-mapping/Unmapped reads' sections. You can use this information to confirm that your read mapping experiment was a success. You can use this information to check the number of reads in the input FASTQ file, how many mapped to the genome, and how many mismatches there are between the reads and reference (if you're not mapping to the same species from which your reads originated, then you'd expect a relatively high mismatch rate).  
2. Click on the eye icon of the `splice junctions.bed` file. This BED-format file provides information about the splice junctions (intron-exon boundaries) in the reference sequence. The meaning of each column is as follows:  
 
*   column 1: chromosome
*   column 2: first base of the intron (1-based)
*   column 3: last base of the intron (1-based)
*   column 4: strand (0: undefined, 1: +, 2: -)
*   column 5: splice junction motif: 0: non-canonical; 1: GT/AG, 2: CT/AC, 3: GC/AG, 4: CT/GC, 5:AT/AC, 6: GT/AT
*   column 6: 0: unannotated, 1: annotated (only if splice junction database is used, so yours will be zeros)
*   column 7: number of uniquely mapping reads crossing the junction
*   column 8: number of multi-mapping reads crossing the junction
*   column 9: maximum spliced alignment overhang. For example, if you had a situation like this: 

```
Spliced read:	ACGTACGT----------ACGTATA
				Exon 1   intron   Exon 2
```
then the spliced alignment overhang will be 7, i.e. the length of `ACGTATA` on exon 2.

3.Finally, click on the eye icon on the `mapped.bam` file. This is the alignment file in SAM format, which represents the reads mapped to the reference.  

**SAM header line**  

The first few lines in the file, all starting with '@' represent information about the alignment. Each of these lines are then followed by key-pair values, separated by colons. These lines are optional, so might not be present in all SAM files. The first of these lines will be the header line, starting with '@HD' and then contain information about the SAM version number (e.g. 'VN:1.4' meaning SAM version 1.4) and sort order 'SO'. 'SO:coordinate' means that the SAM file is sorted by genome co-ordinates, with reads mapping to the start of the first scaffold (sorted alpha-numerically) being at the start of the file and reads mapping to the end of the last scaffold being placed at the end of the file. Other values are 'unknown', 'unsorted', and 'queryname' (sorted by readname).  
Next there are the '@SQ' lines, which are the sequence dictionary lines, providing information about the scaffolds in the SAM file. There is an '@SQ' line for each scaffold in the file. Each line has two mandatory key-pair values, 'SN' - the sequence name, and 'LN' the length of that sequence.  
Finally, you should see two lines, one starting with '@PG' and another starting with '@CO'. The information in these lines are quite similar. The '@CO' line is the command line used to run the program (when you hit 'Execute' in Galaxy, all the parameters are collected from the tool interface, turned into a command line and then run on the Earlham Institute's High-performance Computer. The '@CO' line is this command. The '@PG' line provides information about the 'program' being used. Here You can see that the program name ('PN') is 'STAR' and the version ('VN') is 2.7.5b. The 'CL' tag gives the command line, which is also given in the '@CO' line, described above.

**SAM alignment section** 

The remaining lines in a SAM file provides information about each read aligned to the reference sequence. It has a 11 mandatory fields.  

1. QNAME: The read name of the FASTQ sequence. 
2. Bitwise FLAG. The number in this field represents a cumulative total of other numbers. Each number represents a property of the read, such as which read pair it's from, whether it's mapped, whether it's mapped uniquely, is it a PCR duplicate, etc.
The values that make up this total are in hexadecimal format, rather than decimal.  
TASK: Go to http://www.samformat.info/sam-format-flag-single and input some of the bitwise flag values.  
3. RNAME: Reference sequence NAME of the alignment. If the @SQ header lines are present, RNAME must be present in one of the SQ-SN tag. An unmapped read has a `*' at this field.  
4. POS: 1-based leftmost mapping position of the first matching base. The first base in a reference sequence has coordinate 1. POS is set as 0 for an unmapped read without coordinate.  
5. MAPQ: mapping Quality. This provides information as to whether the read is mapped to one (unique), or more than one (not unique) position.
A value 60 indicates that the read is uniquely mapped and is the maximum value. For multi-mapping reads, the MAPQ is calculated by int(-10*log10(1-
1/Nmap)).
6. CIGAR string. The CIGAR string provides per-base mapping information between the read and the reference.
The CIGAR string consists of one or more continuous key-value pairs. For example 101M, represents 101 continuous aligned bases between the read and the reference, i.e. it maps perfectly for its full length.
55M1D46M denotes an alignment of 55 base pairs, followed by a 1 base pair deletion (D), followed by 46 base pair of alignment.
94M7S denotes an alignment of 94 base pairs, followed by 7 soft-clipped bases. The term 'soft-clipping' refers to bases at the 5' and 3' end of the reads, not included in the alignment, possible due to mismatches at the start or end of the read. Unique to RNA-seq data is the 'N' CIGAR character. This represents "a region of nucleotides not present in the read", usually referring to an intron, so this would most likely denote a spliced read.
7. RNEXT: In paired end sequencing, RNEXT is the scaffold name on which the corresponding paired read is mapped. If both reads map to the same scaffold, it is denoted by '='
8. PNEXT: In paired end sequencing, PNEXT is the position (see 4. POS, above) that the corresponding paired read is mapped. 
9. TLEN: observed template length.  If all segments are mapped to the same reference, the template length equals the number of bases between the first base in leftmost mapped read to the last base in the rightmost mapped read, with leftmost read having positive value and rightmost read having negative value.
10. SEQ: The sequence of the read.
11. QUAL: The ASCII base quality scores of the read.

SAM format also allows for one additional optional field, represented by one or more TAG:TYPE:VALUE entries, each entry separated by a space. The TAG represents a pre-defined property of the alignment. In an optional field, TYPE is a single letter which denotes the format of VALUE, such as a string, integer, byte, array, etc. For example, to store the information of a read's barcode, one would see the entry 'BC:Z:ACGATCT'. 'BC' represents a barcode TAG, 'Z' represents a string TYPE and the value 'ACGATCT' represents the actual barcode.


##Step 3: Calculating gene expression
You now want to calculate the expression of each gene in your SAM alignment file. The co-ordinates of each gene is provided by the `Arabidopsis_thaliana.TAIR10.37.gtf` file that you imported into your history earlier. Click on the eye icon to view it. After a few comment lines (starting with #), each 'feature' of the genome is described over 9 fields as follows:

1.	seqname - name of the chromosome or scaffold
2.	source - name of the program that generated this feature, or the data source.
3.	feature - feature type name, e.g. gene, transcript, exon, etc
4.	start - Start position of the feature on 'seqname'
5.	end - End position of the feature
6.	score - degree of confidence in the feature's existence and coordinates. It may be replaced with a dot when no score is generated.
7.	strand - defined as + (forward) or - (reverse).
8.	frame - One of '0', '1' or '2'. '0' indicates that the first base of the feature is the first base of a codon, '1' that the second base is the first base of a codon, and so on..
9.	attribute - A semicolon-separated list of key-value pairs, providing additional information about each feature, e.g gene\_id "AT1G01010"; transcript\_id "AT1G01010.1"; etc

You will use the Galaxy tool 'FeatureCounts', part of the SubRead package, to count the number of reads aligned to each gene annotated in the `Arabidopsis_thaliana.TAIR10.37.gtf` file. In the tools search box, search for 'feature' and select the '**featureCounts**' tool. 

1. You will be examining the reads mapped to your `RNA STAR....mapped.bam` file, so in the 'Alignment file' box, make sure that file is selected. 
2. The data is unstranded (i.e. it originates from both strands of the genome, so make sure 'unstranded' is selected for 'Specify strand information'.
3. Your 'Gene annotation file' is 'in your history', so change this to the correct option. 
4. The 'Gene annotation file' should be the 'Arabidopsis_thaliana.TAIR10.37.gtf' file.
5. The 'Output format' is 'Gene-ID read-count (MultiQC/DESeq2/edgeR/limma-voom compatible', so make sure that's selected.
6. You don't need to 'Create gene-length file' for this tutorial, so leave that as 'No'.
7. Hit Execute. This should finish in under a minute.

The featureCounts tool creates two output files, a `Summary` file and a `Counts` file. Take a look at the `Summary` file by clicking on the eye icon. You should see that the majority of reads fall under the 'Assigned' category, with smaller numbers (often up to about 20%) of reads being 'Unassigned_' for low mapping quality, ambiguity (e.g. overlapping features, such as overlapping genes on different strands), multi-mapping, or no features to map to.  
Once you're happy that the number reads assigned to your gene models is sufficient, then go ahead and look at the `Counts` file. This is a fairly simple to understand format. The first column in the file represents each gene IDs from the `Arabidopsis_thaliana.TAIR10.37.gtf` file and the second column represent the number of reads aligned to each gene. The combined total of the counts will equal the 'Assigned' value in the `Summary` file. Take a look at the read counts. You'll notice that there's a lot of variation between the counts for each gene. Some genes might have thousands of reads mapped to them, others zero. Obviously, if you were looking at gene expression analysis, this is a good thing.


##Exercise 2 - A multi-sample example
You've gone across the process of calculating gene expression for a single example, but obviously in real life You'll have many samples. Galaxy is an ideal tool for this as you can work on many samples at the same time. Now you'll go across the process of calculating gene expression of a number of samples. You will only be running 5 samples in this exercise as to run hundreds of samples would take too long. 

###Galaxy workflows
One of the many great features of Galaxy is that you can create pre-defined pipelines called 'workflows' that set out how data flows from one tool to another, what version of each tool is required and what type of data is allowed to enter and exit each tool. These workflows can then be shared or published, for example as part of the supplementary data in a paper. Researchers can then download the workflow and replicate what you've done.


You'll use a workflow to run the previous analysis across five samples.  

####Importing the data

1. Create a new history by clicking on the 'cog' icon at the top of the History panel and then click on 'Create New'. Then click on the words 'Unnamed history' and rename it 'Single Cell 2'. Hit return.
2. Navigate to Shared Data > Data Libraries > Galaxy Courses > Single Cell 2020 > Single Cell 2
3. Select all of the `fastq.gz` files and then click on 'To History' > as a Collection. 
4. Change the 'Collection type' to 'List of Pairs' and click 'Continue'. This will bring you to the pairing dialog. In the '0 unpaired forward' field, change '\_1' to 'R1'. In the '0 unpaired reverse' field, change '\_2' to 'R2'. Then click on 'Auto-pair'. The R1 and R2 pairs for each of the A1 - A5 datasets should now be paired (double-check that they are). At the bottom of the dialog enter a 'Name' for your data collection.  Call it '5 sample pairs' and click on 'Create list'
5. Next uncheck all of the datasets and select the `Arabidopsis_thaliana.TAIR10.37.gtf` file. Again 'To History' > As dataset > Import
6. Next, click on 'Analyze Data' to take you back to the home screen. 


####Running the workflow
You now want to run the same steps on the collection of fastq files that you ran in your single example earlier. You can automate this by using workflows. 

1. Navigate to Shared Data > Workflows
2. You will see a list of published workflows, one of which is called 'Single Cell reads to expression matrix'.
3. Click on the downward facing arrow to the right of the box and select 'Run'
4. The workflow interface will appear. There's only two options that need selecting here - the reads and the annotation.
5. Under '1. Input dataset', make sure the Arabidopsis GTF annotation file is selected.
6. Under 2. Input dataset collection, make sure your '5 sample pairs' is selected. 
7. If you click on 'Expand to full workflow form", you can see the tools that you used on the single sample in the rest of the workflow below, along with some text manipulation tools to create the expression matrix.
8. Click on 'Run workflow' at the top of the page.

In a few seconds new jobs will appear in the History and some will start running (turn yellow). These have been submitted by the workflow. As the input of some of the tools are dependent on the output of others, downstream tools will remain grey until their upstream tools have finished running. Wait for all history items to turn green (this should take about 20 minutes) and then take a look at the last file. This is an example of a small expression matrix, the input of all further scRNA-seq analyses. The rows of the expression matrix start with gene names (e.g. AT1G01010), and the columns start with sample names (e.g. A1, A2, etc.). The rest of the matrix represents read counts for each gene in each sample.


##SummaryYou have learned how to use Galaxy to take a single sequencing library, map the reads to a genome and calculate per-gene expression. Further, you have used a Galaxy workflow to do the same thing over multiple samples. For the remainder of the course, you will use pre-calculated expression matrices for larger sample sizes, but that have been produced using exactly the same methods. 