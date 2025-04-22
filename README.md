# ChIP-seq-analysis-and-visualization-using-Galaxy
This is a ChIP-seq analysis workflow using Galaxy using a [dataset on tissues taken from E18.5 mouse embryo](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA886671&o=acc_s%3Aa) from [Zhang et al.'s study](https://www.nature.com/articles/s41467-024-49684-1) on CTCF mutations. 

## Table of contents

- [Step 1: Importing data](#step-1-importing-data)
- [Step 2: Quality control using FastQC](#step-2-quality-control-using-fastqc)
- [Step 3: Trimming using Trimmomatic ](#step-3-trimming-using-trimmomatic)
- [Step 4: Mapping reads to mouse(mm10) genome using Bowtie2](#step-4-mapping-reads-to-mouse(mm10)-genome-using-bowtie2)
- [Step 5: Filter alignment based on quality using Samtools view ](#step-5-filter-alignment-based-on-quality-using-samtools-view)
- [Step 6: Find peaks using MACS2 callpeak](#step-6-find-peaks-using-macs2-callpeak)
- [Step 7: Mapping peaks to known genomic features using ChIPseeker](#step-7-mapping-peaks-to-known-genomic-features-using-chipseeker)
- [Step 8: Get the profile of the peaks](#step-8-get-profile-of-the-peaks)
- [Step 9: Visualizing the peaks using IGV](#step-9-visualizing-peaks-using-IGV)
- [Step 10: Motif analysis using memeChIP](#step-9-motif-analysis-using-memeChIP)
- [Step 11: Gene Ontology](#step-9-gene-ontology)

## workflow
### Step 1: Importing data 
For this analysis, we only used the lung tissue data. 
the corresponsindg SRR numbers are:

Homozygous Ctcf R567 mutant
Input → SRR21787371 
IP → SRR21787377 

Wildtype
Input → SRR21787372 
IP → SRR21787378 

From [this website](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA886671&o=acc_s%3Aa), click on the boxes next to the 4 SRR accession and then press the the galaxy button shown in the picture below. 

Go to **tools** --> **Get Data** --> **Download and Extract Reads in FASTQ format from NCBI SRA**

Use the following settings:
* select input type: list of SRA accession, one per line
* Under sra accession list, input your SRA collection
* select output format: gzip compressed fastqc


### Step 2: Quality control using FastQC 

Run FastQC twice: Once with the Paired-end data (fastq-dump) wt as the input and once with Paired-end data (fastq-dump) ctcf mutant as the input 


### Step 3: Trimming using Trimmomatic 

Run Trimmomatic twice: Once with the *Paired-end data (fastq-dump) wt* as the input and once with *Paired-end data (fastq-dump) ctcf mutant* as the input 

Use the following settings:
* Perform initial ILLUMINACLIP step? : Yes
* Adapter sequence: truSeq3(paired-end)
* Average quality required : 30
* Quality score encoding: phred 33

### Step 4: Mapping reads to mouse(mm10) genome using Bowtie2

Run Bowtie2 twice: Once with *trimmomatic on wt* as the input and once with *trimmomatic on ctcf mutant* as the input 

Use the following settings:
*  set paired-end options: yes
* --no-mixed
* -no-discordant
* Reference genome: Mouse (mus musculus) : mm10
* Select analysis mode
* Presets: Very sensitive end-to-end





### Step 5: Filter alignment based on quality using Samtools view 

Use the following settings:


### Step 6: Find peaks using MACS2 callpeak 

Use the following settings:

### Step 7: Mapping peaks to known genomic features using ChIPseeker 

Use the following settings:

### Step 8: Get the profile of the peaks 

Use the following settings:

### Step 9: Visualizing the peaks using IGV 

Use the following settings:

### Step 10: Motif analysis using memeChIP 

Use the following settings:

### Step 11: Gene Ontology  

Use the following settings:
