# From Raw To Coverage: An automated pipeline for RNA-Seq Analysis.

## Purpose and General Information

This script will perform every task from trimming Illumina RNA-Seq raw data (fastq.gz files) to write coverage files (bedgraph and IGV files). It also counts the start position of reads and end position if a paired-end experiment is supplied:

```
single-end

                     RNA FRAGMENT (INSERT)
5prime                                                3prime
------------------------------------------------------------

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

---

Furthermore, it will filter bam files in order to keep only the R1 files, for they are useful to find TSS in dRNASeq experiments (e.g. TEXminus vs. TEXplus experiments) the filtered bams are placed within tssarinput directory.  

This program features a few modules. Each step usually relies on the previous one, but they will be executed only if the output directory is not created. This is a simple way to enable the rerun of the late steps without starting from the beginning.  

Those modules (or steps) are summarized below:  

1. Trimming files
2. Downloading reference genome and annotation from NCBI RefSeq; building HISAT2 index and aligning to ref. genome
3. Filtering uniquely aligned reads
4. Converting SAM files to BAM and sorting them by read name (performed by SAMtools)
5. Adjusting position of multi-mappers using MMR and [removing pairs that do not match each other](https://github.com/ratschlab/mmr/issues/5) (custom Rscript)
6. Creating coverage files (bedgraph and igv format) using deepTools
7. Creating TSSAR input BAM files
8. Creating five prime profiling coverage (bedgraph format) using bedtools
9. Creating three prime profiling coverage (bedgraph format) using bedtools

## Usage and Requisites

```
bash frtc.sh <threads> <maxfragsize> <read_size> <spp> <url>

threads [INT]:     number of threads to be passed to nested programs
                   maximum value: 99

maxfragsize [INT]: maximum insert size from the leftmost end to
                   the rightmost of an paired-end alignment.
                   This parameter is ignored if single-end libraries
                   are being processed.

read_size [INT]:   approximate size of reads. For an Illumina paired-end
                   experiment of 150x2 cycles, use 150.

spp [CHAR]:        prefix for genome and annotation files

url [CHAR]:        URL pointing to the genome file in NCBI RefSeq FTP Server.
                   One may find the link for RefSeq directory on the top right
                   corner of the NCBI Assembly page for a given genome.
```

e.g.:

```
bash frtc.sh 6 1000 150 Hsalinarum ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/805/GCF_000006805.1_ASM680v1/GCF_000006805.1_ASM680v1_genomic.fna.gz
```

This command line will run frtc.sh with six threads, with a maximum insert size of 1000, with read size of 150 nt, creating Hsalinarum.fa and Hsalinarum.gff files under misc directory to store genome and annotation, respectively.  

There is a directory tree prerequisite to run frtc and you have to create it.  
An adequate directory tree to run frtc must look like:

```
yourDirectory
├── removeInconsistentPairs.R	provided
├── frtc.sh			provided
├── misc			notProvided
│   └── adap.fa			
└── raw				notProvided
    ├── S1_R1.fastq.gz		
    ├── S1_R2.fastq.gz
    ├── S2_R1.fastq.gz
    ├── S2_R2.fastq.gz
    ├── S3_R1.fastq.gz
    ├── S3_R2.fastq.gz
    ├── S4_R1.fastq.gz
    └── S4_R2.fastq.gz
```

Raw file names must be ended with R1 or R2 and must have fastq.gz extension (i.e. gzip compressed).  

If the libraries are paired-end:  

e.g. S1_R1.fastq.gz S1_R2.fastq.gz  

If the libraries are single-end:  

e.g. S1_R1.fastq.gz  

Single-end or paired-end mode can be set in frtc.sh. This script was written using paired-end as default mode but you are able to change to single-end mode by setting **pairedend="n"** on **"CUSTOM VARIABLES"** section.

Also, if you are processing paired-end libraries, check if R1 and R2 corresponding reads have the same name. For example, if the first entry in R1 file has the header "@SRR9999999.1", the first entry in R2 file must have the same header "@SRR9999999.1". Everything should be the same before the first space character. NCBI SRA commonly provides files that present headers like "@SRR9999999.1.1" for R1 and "@SRR9999999.1.2" for R2, which may cause problems during the multi-mapper resolution step of this pipeline.

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
curl v7.58.0 (@ PATH)  
samtools v1.7 (@ PATH)  
mmr default version (@ PATH)  
deeptools v3.1.1 (@PATH)  
bedtools v2.26.0 (also tested with v2.21.0) (@ PATH)  

All the prerequisites will be checked before start the processing.

```
   __o
 _ \<_
(_)/(_)
```

