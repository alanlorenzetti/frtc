# From Raw To Coverage: An automatic pipeline for RNA-Seq Analysis.

## USAGE

```{bash}
frtc.sh <threads> <maxfragsize> <read_size> <spp> <url>
```

e.g.:
frtc.sh 6 1000 150 Hsalinarum ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/805/GCF_000006805.1_ASM680v1/GCF_000006805.1_ASM680v1_genomic.fna.gz

there is a directory tree prerequisite to run frtc.
an adequate directory tree to run frtc must look like

yourDirectory
├── removeInconsistentPairs.R			 provided
├── frtc.sh					 provided
├── misc					 notProvided
│   ├── adap.fa
└── raw 					 notProvided
    ├── S1_R1.fastq
    ├── S1_R2.fastq
    ├── S2_R1.fastq
    ├── S2_R2.fastq
    ├── S3_R1.fastq
    ├── S3_R2.fastq
    ├── S4_R1.fastq
    └── S4_R2.fastq

raw file names must be ended with R1 or R2 and
must have fastq extension (i.e. uncompressed)
if the libraries are paired-end
e.g. S1_R1.fastq
     S1_R2.fastq

if the libraries are single-end
e.g. S1_R1.fastq

adap.fa must be a fasta file containing
what adapters would look like if they are sequenced.
for more information, check https://support.illumina.com/bulletins/2016/12/what-sequences-do-i-use-for-adapter-trimming.html
the file must look like:

>indexedAdapter
TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC
>nonIndexedAdapter
GATCGTCGGACTGTAGAACTCTGAACGTGTAGA

but with no hashes or spaces

furthermore, you have to install manually the program requisites.
maybe newer versions are compatible, but these were the ones used to build and test
this script

trimmomatic v0.36 (must be located @ /opt/Trimmomatic-0.36/trimmomatic-0.36.jar)
hisat2 v2.1.0 (@ PATH)
curl 7.47 (also tested with v7.37) (@ PATH)
samtools v1.3.1 (@ PATH)
mmr default version (@ PATH)
deeptools 2.5.3 (@PATH)
bedtools v2.26.0 (also tested with v2.21.0) (@ PATH)

all the prerequisites will be checked before running


PURPOSE AND GENERAL INFO


this script will perform every task
from trimming Illumina RNASeq raw data (fastq files)
to write coverage files (bedgraph and IGV files)

it also counts the start position of reads
and end position if a paired-end experiment is supplied

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


\>>>>>>>>>>>>>>>>>>>>>>>>>>>>WARNING<<<<<<<<<<<<<<<<<<<<<<<<<<<<
be aware this feature doesn't work for dUTP library preparations;
be aware that uniq mode is not supported for
fiveprime, threeprime and TSSAR input results
---------------------------------------------------------------

furthermore, it will filter bam files in order to keep
only the R1 files, for they are useful to find TSS
in dRNASeq experiments (e.g. TEXminus vs. TEXplus experiments)
the filtered bams are placed within tssarinput directory

there are a few modules.
each step usually relies on the previous one
but they will be executed only if the output
directory is not created.
this is a simple way to enable the rerun of the
late steps without starting from the beginning

those modules (or steps) are summarized below:
1. trimming files
2. downloading reference genome and annotation from NCBI RefSeq; building HISAT2 index and aligning to ref. genome
3. filtering uniquely aligned reads
4. converting SAM files to BAM and sorting them by read name (performed by SAMtools)
5. adjusting position of multi-mappers using MMR and removing pairs that do not match each other (custom Rscript)
6. creating coverage files (bedgraph and igv format) using deepTools
7. creating TSSAR input BAM files
8. creating five prime profiling coverage (bedgraph format) using bedtools
9. creating three prime profiling coverage (bedgraph format) using bedtools


