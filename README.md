# From Raw To Coverage: An automatic pipeline for RNA-Seq Analysis.

## Purpose and General Info

This script will perform every task from trimming Illumina RNASeq raw data (fastq files) to write coverage files (bedgraph and IGV files). It also counts the start position of reads and end position if a paired-end experiment is supplied:

```
single-end
                   end position (@ threeprime directory; only meaningful in a few cases; use with caution)
R1                 |
------------------->
|
start position (@ fiveprime directory)

paired-end
                     RNA FRAGMENT (INSERT)
5prime                                                3prime
------------------------------------------------------------

R1                                                         end position (@ threeprime directory)
--------------->                                           |
|                                       <-------------------
start position (@ fiveprime directory)                    R2
```

**Warning**  

Be aware this feature doesn't work for dUTP library preparations and be aware that uniq mode is not supported for fiveprime, threeprime and TSSAR input results.  

Furthermore, it will filter bam files in order to keep only the R1 files, for they are useful to find TSS in dRNASeq experiments (e.g. TEXminus vs. TEXplus experiments) the filtered bams are placed within tssarinput directory.  

This program features a few modules. Each step usually relies on the previous one, but they will be executed only if the output directory is not created. This is a simple way to enable the rerun of the late steps without starting from the beginning.  

Those modules (or steps) are summarized below:  

1. Trimming files
2. Downloading reference genome and annotation from NCBI RefSeq; building HISAT2 index and aligning to ref. genome
3. Filtering uniquely aligned reads
4. Converting SAM files to BAM and sorting them by read name (performed by SAMtools)
5. Adjusting position of multi-mappers using MMR and removing pairs that do not match each other (custom Rscript)
6. Creating coverage files (bedgraph and igv format) using deepTools
7. Creating TSSAR input BAM files
8. Creating five prime profiling coverage (bedgraph format) using bedtools
9. Creating three prime profiling coverage (bedgraph format) using bedtools

## Usage

```{shell}
frtc.sh <threads> <maxfragsize> <read_size> <spp> <url>
```

e.g.:

```{shell}
frtc.sh 6 1000 150 Hsalinarum ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/805/GCF_000006805.1_ASM680v1/GCF_000006805.1_ASM680v1_genomic.fna.gz
```

There is a directory tree prerequisite to run frtc and you have to create it.  
An adequate directory tree to run frtc must look like:

```
yourDirectory
├── removeInconsistentPairs.R	provided
├── frtc.sh			provided
├── misc			notProvided
│   ├── adap.fa			
└── raw 			notProvided
    ├── S1_R1.fastq		
    ├── S1_R2.fastq
    ├── S2_R1.fastq
    ├── S2_R2.fastq
    ├── S3_R1.fastq
    ├── S3_R2.fastq
    ├── S4_R1.fastq
    └── S4_R2.fastq
```

Raw file names must be ended with R1 or R2 and must have fastq extension (i.e. uncompressed).  

If the libraries are paired-end:  

e.g. S1_R1.fastq S1_R2.fastq  

If the libraries are single-end:  

e.g. S1_R1.fastq  

**adap.fa** must be a fasta file containing what adapters would look like if they are sequenced. For more information, check https://support.illumina.com/bulletins/2016/12/what-sequences-do-i-use-for-adapter-trimming.html  

The file must look like:  

```
>indexedAdapter
TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC
>nonIndexedAdapter
GATCGTCGGACTGTAGAACTCTGAACGTGTAGA
```

Furthermore, you have to install manually the program requisites. Maybe newer versions are compatible, but these were the ones used to build and test this script:  

trimmomatic v0.36 (must be located @ /opt/Trimmomatic-0.36/trimmomatic-0.36.jar)  
hisat2 v2.1.0 (@ PATH)  
curl 7.47 (also tested with v7.37) (@ PATH)  
samtools v1.3.1 (@ PATH)  
mmr default version (@ PATH)  
deeptools 2.5.3 (@PATH)  
bedtools v2.26.0 (also tested with v2.21.0) (@ PATH)  

All the prerequisites will be checked before start the processing.
