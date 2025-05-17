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
- [Step 8: Visualizing the peaks using IGV](#step-9-visualizing-peaks-using-IGV)
- [Step 9: Get the profile of the peaks](#step-8-get-profile-of-the-peaks)
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

Run **FastQC** twice: Once with the *Paired-end data (fastq-dump) wt* as the input and once with *Paired-end data (fastq-dump) ctcf mutant* as the input 


### Step 3: Trimming using Trimmomatic 

Run **Trimmomatic** twice: Once with the *Paired-end data (fastq-dump) wt* as the input and once with *Paired-end data (fastq-dump) ctcf mutant* as the input 

Use the following settings:
* Perform initial ILLUMINACLIP step? : Yes
* Adapter sequence: truSeq3(paired-end)
* Average quality required : 30
* Quality score encoding: phred 33

Name the outputs: *trimmomatic on wt* and *trimmomatic on ctcf mutant*

### Step 4: Mapping reads to mouse(mm10) genome using Bowtie2

Run **Bowtie2** twice: Once with *trimmomatic on wt* as the input and once with *trimmomatic on ctcf mutant* as the input 

Use the following settings:
*  set paired-end options: yes
* --no-mixed
* -no-discordant
* Reference genome: Mouse (mus musculus) : mm10
* Select analysis mode
* Presets: Very sensitive end-to-end

Name the outputs: *Bowtie2 on wt* and *Bowtie2 on ctcf mutant*


### Step 5: Filter alignment based on quality using Samtools view 

Run **Samtools view** twice: Once on  *Bowtie2 on wt* and once on *Bowtie2 on ctcf mutant* 

Use the following settings:
- What would you like to look at?:A filtered/subsampled section of reads
     - Configure filters
       - Filter by quality : 30 //Only uniquely mapped reads with MAPQ > 30 were retained

Name the outputs: *Samtools view on wt* and *Samtools view on ctcf mutant*

### Step 6: Find peaks using MACS2 callpeak 

Run **MACS2 callpeak** twice: 
- Once on WT
    - ChIP-Seq Treatment File : Result of  *Samtools view on  wt IP*
    - ChIP-Seq Control File : Result of  *Samtools view on  wt Input*
- Once on ctcf mutation
    - ChIP-Seq Treatment File : Result of * Samtools view on  ctcf mutant IP*
    - ChIP-Seq Control File : Result of  *Samtools view on  ctcf mutant Input*

Use the following settings:
- Format of Input Files : paired-end BAM
- Effective genome size : M.musculus (1.87e9)
  
Name the output: *MACS2 callpeak on wt* and once on *MACS2 callpeak on ctcf mutant*

### Step 7: Mapping peaks to known genomic features using ChIPseeker 

Download a gtf file of mouse basic gene annotation from GENCODE["https://www.gencodegenes.org/mouse/release_M10.html"]. 
- Content: Basic gene annotation
- Region: ALL
- Download: GTF

Upload this GTF file onto galaxy

Run **ChIPseeker** twice: once on *MACS2 callpeak on wt* and once on *MACS2 callpeak on ctcf mutant*.

Use the following settings:
- Annotation source : Use a GTF from history
    - M10(GRCm38.p4)_annotation.gtf //the GTF file from GENCODE
- Output Format : tabular
- Output PDF of plots?: yes

### Step 8: Visualizing the peaks using IGV 

Run **BamCoverage** twice: Once on *Samtools view on wt* and once on *Samtools view on ctcf mutant*

Use the following settings:
- Bin size: 10
- Scaling/Normalization method : Normalize to reads per kilobase per million
-  Show advanced options : yes
     - Scale factors
       - Wt = 1
       - Ctcf Mutant = 0.70

### Step 9: Get the profile of the peaks 

Create and upload list of genes for both wt and ctcf mutant
- Download the annotated peaks output of chIPseeker, delete the duplicate genes, create a txt file with just the list of genes.
- Upload this txt file to galaxy

Run **Filter GTF data by attribute values_list**
Use the following settings:
- Filter : M10(GRCm38.p4)_annotation.gtf
- Using attribute name: gene_Id
- attribute values : txt with gene ids 

Run *computeMatrix*
Use the following settings:
- Regions to plot : result of filter GTF data by wt_geneId
- Score file : 
    -  *bamCoverage on wt
    - *bamCoverage on ctcf mutant*
computeMatrix has two main output options : reference-point
The reference point for the plotting : beginning of region
--beforeRegionStartLength : 1000
--afterRegionStartLength: 1000
--binSize : 10


plotProfile
Input : result of computeMatrix
--plotHeight : 10
--plotWidth : 20
--plotType: lines
Make one plot per group of regions : Yes
  

### Step 10: Motif analysis using memeChIP 

Use the following settings:

### Step 11: Gene Ontology  

Use the following settings:
