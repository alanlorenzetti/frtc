#!/bin/bash

# alorenzetti
# built and tested on Ubuntu 16.04.3 LTS 64bit
# also tested on Debian GNU/Linux jessie/sid and Ubuntu 18.04.1

version=0.6.4
lastupdate=20191007

# please, check the README.md file before using this script
# there is also a version of the manual on the end of this file
# which can be accessed using ./frtc.sh --help

# starting an if statement to show the help
# which is presented on the end of the file
# this if statement only ends on the end of
# this script
if [ "$1" != "--help" ] ; then

# showing usage hints if no arguments are supplied
if [ $# -ne 5 ] ; then echo "
From Raw To Coverage (frtc):
A tool to process Illumina RNA-Seq data.
https://github.com/alanlorenzetti/frtc
Version: $version
Last update: $lastupdate

Usage:

bash frtc.sh <threads> <maxfragsize> <read_size> <spp> <url>

threads [INT]:     number of threads to be passed to nested programs.
                   maximum value: 99

maxfragsize [INT]: maximum insert size from the leftmost end to
                   the rightmost of an paired-end alignment.
                   This parameter is ignored if single-end libraries
                   are being processed.

read_size [INT]:   approximate size of reads. For an Illumina paired-end
                   experiment of 150x2 cycles, use 150.

spp [VARCHAR]:     prefix for genome and annotation files

url [VARCHAR]:     URL pointing to the genome file in NCBI RefSeq FTP Server.
                   One may find the link for RefSeq directory on the top right
                   corner of the NCBI Assembly page for a given genome.

e.g.:
bash frtc.sh 6 1000 150 Hsalinarum ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/805/GCF_000006805.1_ASM680v1/GCF_000006805.1_ASM680v1_genomic.fna.gz

this command line will run frtc.sh with six threads,
with a maximum insert size of 1000, with read size of 150 nt,
creating Hsalinarum.fa and Hsalinarum.gff files under misc
directory to store genome and annotation, respectively.

Check the extended documentation using: ./frtc.sh --help

" && exit 1

fi

####################################
# CUSTOM VARIABLES
####################################
# paired-end (y or n)
pairedend="y"

# should the program create uniquely alignment files? (y or n)
# warning: fiveprime, threeprime and TSSAR input will not be uniq
# this feature is available only for coverage
uniqaln="y"

# should the output be normalized by the amount of aligned reads
# of the library aligning the highest amount of reads? (y or n)
normalize="n"

####################################
# ARGUMENT VARIABLES
####################################
# number of threads to run the applications
threads=$1
if [ "$threads" -gt 99 ] ; then echo >&2 "Threads argument value can not be greater than 99" ; exit 1 ; fi

# max fragment size for hisat2 and bamCoverage
maxfragsize=$2

# read size for mmr
readsize=$3

# spp of interest and url to download the genome from NCBI RefSeq
spp=$4
url=$5
urlannot=$(echo $url | sed 's/fna.gz/gff.gz/')

####################################
# HARD CODED VARIABLES
####################################
prefixes=`ls raw/*.fastq.gz | sed 's/^raw\///;s/_R[12].*$//' | sort | uniq`
rawdir="raw"
miscdir="misc"
trimmeddir="trimmed"
samdir="sam"
samuniqdir="sam_uniq"
bamdir="bam"
bamuniqdir="bam_uniq"
mmrdir="bam_mmr"
coveragedir="coverage"
coverageuniqdir="coverage_uniq"
fiveprimedir="fiveprime"
threeprimedir="threeprime"
tssarinputdir="tssarinput"

# if the output must be normalized
# the pipeline will create coverage,
# coverage_uniq, fiveprime and threeprime
# directories with the suffix _norm

if [ "$normalize" == "y" ] ; then
        coveragedir="coverage_norm"
        coverageuniqdir="coverage_uniq_norm"
        fiveprimedir="fiveprime_norm"
        threeprimedir="threeprime_norm"
fi

####################################
# PROGRAM STAMP
####################################
echo "From Raw To Coverage (frtc):
A tool to process Illumina RNA-Seq data.
https://github.com/alanlorenzetti/frtc
Version: $version
Last update: $lastupdate"

####################################
# CALL AND DATE
####################################
dateAndTime=`date`
echo "$dateAndTime"
echo "Call: $0 $@"

####################################
# CHECKING DEPENDENCIES
####################################
echo "Checking dependencies..."

# checking files
if [ -z "$prefixes" ] ; then echo >&2 "FASTQ files not found. Aborting" ; exit 1 ; fi
if [ ! -e misc/adap.fa ] ; then echo >&2 "Adapter sequences not found. Aborting" ; exit 1 ; fi

# checking programs

for i in curl hisat2 samtools mmr deeptools R bedtools; do
        command -v $i > /dev/null >&1 || { echo >&2 "$i is not installed. Aborting" ; exit 1; }
done

R --slave -e 'if(!require("rbamtools", quietly=T)){quit(save="no", status=1)}else{quit(save="no", status=0)}'
if [ $? == 1 ] ; then echo >&2 "rbamtools package is not installed. Aborting" ; exit 1 ; fi
if [ ! -e /opt/Trimmomatic-0.39/trimmomatic-0.39.jar ] ; then echo >&2 "trimmomatic is not installed. Aborting" ; exit 1; fi

# checking scripts
if [ ! -e removeInconsistentPairs.R ] && [ "$pairedend" == "y" ] ; then echo >&2 "Missing removeInconsistentPairs.R script. Aborting" ; exit 1 ; fi

echo "Done!"

####################################
# PROCESSING STARTS HERE
####################################

######################
# trimming
######################
# if used in paired-end mode
# this step is going to output four files for each library
# e.g.
# S1_R1-paired.fastq.gz and its pair S1_R2-paired.fastq.gz
# S1_R1-unpaired.fastq.gz containing reads R1 that had its R2 eliminated by the trimming (unpaired or orphan R1)
# S1_R2-unpaired.fastq.gz containgin reads R2 that had its R1 eliminated by the trimming (unpaired or orphan R2)
#
# if used in single-end mode
# this step will provide just one file as output
# e.g.
# S1_R1-unpaired.fastq.gz

if [ ! -d $trimmeddir ] ; then
        mkdir $trimmeddir

        echo "Step 1: Starting Trimming"

        if [ "$pairedend" == "y" ] ; then
                for prefix in $prefixes ; do
                        echo "Trimming $prefix"
                        R1=$rawdir/$prefix"_R1.fastq.gz"
                        R2=$rawdir/$prefix"_R2.fastq.gz"
                        outpairedR1=$trimmeddir/$prefix"-paired_R1.fastq.gz"
                        outpairedR2=$trimmeddir/$prefix"-paired_R2.fastq.gz"
                        outunpairedR1=$trimmeddir/$prefix"-unpaired_R1.fastq.gz"
                        outunpairedR2=$trimmeddir/$prefix"-unpaired_R2.fastq.gz"
                        logfile=$trimmeddir/$prefix".log"

                        java -jar /opt/Trimmomatic-0.39/trimmomatic-0.39.jar PE \
                        -threads $threads \
                        $R1 $R2 \
                        $outpairedR1 $outunpairedR1 \
                        $outpairedR2 $outunpairedR2 \
                        ILLUMINACLIP:$miscdir/adap.fa:1:30:10 \
                        SLIDINGWINDOW:4:30 \
                        MINLEN:20 > $logfile 2>&1
                done
        else
                for prefix in $prefixes ; do
                        echo "Trimming $prefix"
                        R1=$rawdir/$prefix"_R1.fastq.gz"
                        outunpairedR1=$trimmeddir/$prefix"-unpaired_R1.fastq.gz"
                        logfile=$trimmeddir/$prefix".log"

                        java -jar /opt/Trimmomatic-0.39/trimmomatic-0.39.jar SE \
                        -threads $threads \
                        $R1 \
                        $outunpairedR1 \
                        ILLUMINACLIP:$miscdir/adap.fa:1:30:10 \
                        SLIDINGWINDOW:4:30 \
                        MINLEN:20 > $logfile 2>&1
                done
        fi

        echo "Done!"
fi

######################
# aligning
######################

if [ ! -d $samdir ] ; then
        mkdir $samdir

        echo "Step 2: Starting Alignment to ref. Genome"

        if [ ! -f $miscdir/$spp".1.ht2" ] ; then
                # downloading ref genome and annotation from NCBI RefSeq
                # it will remove region features from GFF, for they are annoying
                # to observe on genome browser
                echo "Downloading genome"
                curl -u anonymous: $url 2> /dev/null | zcat > $miscdir/$spp".fa" || { echo >&2 "Genome file download failed. Aborting" ; exit 1; }
                echo "Downloading annotation"
                curl -u anonymous: $urlannot 2> /dev/null | zcat | \
                awk -v OFS="\t" -v FS="\t" '{if(/^#/){print}else{if($3 != "region"){print}}}' > $miscdir/$spp".gff" || { echo >&2 "Annotation file download failed. Aborting" ; exit 1; }
                echo "Building HISAT2 index"
                hisat2-build $miscdir/$spp".fa" $miscdir/$spp > /dev/null 2>&1
                echo "Done!"
        fi

        if [ "$pairedend" == "y" ] ; then
                # aligning reads to ref genome
                # rdg and rfg are set to 1000,1 to make it impossible
                # to report reads with insertions or deletions considering the ref genome
                # no discordant, mixed, spliced or softclipped alignments are allowed
                # we allow the reporting of at most 1000 alignments per read
                # because we are gonna choose the best one using MMR downstream
                for prefix in $prefixes ; do
                        echo "Aligning $prefix (paired)"
                        # aligning paired reads
                        hisat2 \
                        --no-discordant --no-mixed \
                        --no-spliced-alignment \
                        --no-softclip \
                        --rdg 1000,1 \
                        --rfg 1000,1 \
                        --rna-strandness FR \
                        -k 1000 \
                        -X $maxfragsize \
                        -p $threads \
                        -x $miscdir/$spp \
                        -1 $trimmeddir/$prefix"-paired_R1.fastq.gz" \
                        -2 $trimmeddir/$prefix"-paired_R2.fastq.gz" \
                        --summary-file $samdir/$prefix"-paired.log" | grep "^@\|YT:Z:CP" > $samdir/$prefix"-paired.sam"

                        # aligning unpaired R1
                        echo "Aligning $prefix (unpaired R1)"
                        # unpaired R1 comes from the same strand of the RNA molecule in the sample
                        hisat2 \
                        --no-spliced-alignment \
                        --no-softclip \
                        --rdg 1000,1 \
                        --rfg 1000,1 \
                        --rna-strandness F \
                        -k 1000 \
                        -p $threads \
                        -x $miscdir/$spp \
                        -U $trimmeddir/$prefix"-unpaired_R1.fastq.gz" \
                        --summary-file $samdir/$prefix"-unpaired_R1.log" | grep "^@\|YT:Z:UU" > $samdir/$prefix"-unpaired_R1.sam"

                        # aligning unpaired R2
                        echo "Aligning $prefix (unpaired R2)"
                        # unpaired R2 comes from an artificial opposite strand of the RNA molecule in the sample (first cDNA strand synthesized)
                        hisat2 \
                        --no-spliced-alignment \
                        --no-softclip \
                        --rdg 1000,1 \
                        --rfg 1000,1 \
                        --rna-strandness R \
                        -k 1000 \
                        -p $threads \
                        -x $miscdir/$spp \
                        -U $trimmeddir/$prefix"-unpaired_R2.fastq.gz" \
                        --summary-file $samdir/$prefix"-unpaired_R2.log" | grep "^@\|YT:Z:UU" > $samdir/$prefix"-unpaired_R2.sam"

                        samtools view -@ $threads -h $samdir/$prefix"-unpaired_R2.sam" | \
                        awk -v OFS="\t" -v FS="\t" '{if(/^@/){print}\
                                                else{if($2 == 0 || $2 == 256){$1 = $1"_R2" ; $2 = $2+16; print}\
                                                else{if($2 == 16 || $2 == 272){$1 = $1"_R2" ; $2 = $2-16; print}\
                                                else{$1 = $1"_R2" ; print}}}}' > $samdir/$prefix"-unpaired_R2.tmp"

                        mv $samdir/$prefix"-unpaired_R2.tmp" $samdir/$prefix"-unpaired_R2.sam"
                done
        else
                for prefix in $prefixes ; do
                        # aligning unpaired R1
                        echo "Aligning $prefix (unpaired R1)"

                        hisat2 \
                        --no-spliced-alignment \
                        --no-softclip \
                        --rdg 1000,1 \
                        --rfg 1000,1 \
                        --rna-strandness F \
                        -k 1000 \
                        -p $threads \
                        -x $miscdir/$spp \
                        -U $trimmeddir/$prefix"-unpaired_R1.fastq.gz" \
                        --summary-file $samdir/$prefix"-unpaired_R1.log" > $samdir/$prefix"-unpaired_R1.sam"
                done
        fi

        echo "Done!"
fi

######################
# filtering unique alignments
######################

if [ "$uniqaln" == "y" ] ; then
        if [ ! -d $samuniqdir ] ; then
                mkdir $samuniqdir

                echo "Step 3: Filtering uniquely aligned reads"

                if [ "$pairedend" == "y" ] ; then
                        for prefix in $prefixes ; do
                                grep "^@\|NH:i:1$" $samdir/$prefix"-paired.sam" > $samuniqdir/$prefix"-paired.sam"
                                grep "^@\|NH:i:1$" $samdir/$prefix"-unpaired_R1.sam" > $samuniqdir/$prefix"-unpaired_R1.sam"
                                grep "^@\|NH:i:1$" $samdir/$prefix"-unpaired_R2.sam" > $samuniqdir/$prefix"-unpaired_R2.sam"
                        done

                        echo "Done!"
                else
                        for prefix in $prefixes ; do
                                grep "^@\|NH:i:1$" $samdir/$prefix"-unpaired_R1.sam" > $samuniqdir/$prefix"-unpaired_R1.sam"
                        done
                fi
        fi
fi

######################
# converting to bam and sorting (mult)
######################

if [ ! -d $bamdir ] ; then
        mkdir $bamdir

        echo "Step 4.1: Converting to BAM and sorting by read name"

        if [ "$pairedend" == "y" ] ; then
                for prefix in $prefixes ; do
                        samtools merge \
                        -@ $threads \
                        $bamdir/$prefix"-merged.bam" \
                        $samdir/$prefix"-paired.sam" \
                        $samdir/$prefix"-unpaired_R1.sam" \
                        $samdir/$prefix"-unpaired_R2.sam" 2> /dev/null

                        # q 1 grants a good quality of alignment
                        samtools view \
                        -@ $threads \
                        -h \
                        -b \
                        -q 1 $bamdir/$prefix"-merged.bam" | \
                        samtools sort \
                        -@ $threads \
                        -n \
                        -o $bamdir/$prefix"-sorted.bam" 2> /dev/null
                done
        else
                for prefix in $prefixes ; do
                        samtools view \
                        -@ $threads \
                        -h \
                        -b \
                        -q 1 $samdir/$prefix"-unpaired_R1.sam" | \
                        samtools sort \
                        -@ $threads \
                        -n \
                        -o $bamdir/$prefix"-sorted.bam" 2> /dev/null
                done
        fi

        echo "Done!"
fi

######################
# converting to bam and sorting (uniq)
######################

if [ "$uniqaln" == "y" ] ; then
        if [ ! -d $bamuniqdir ] ; then
                mkdir $bamuniqdir

                echo "Step 4.2: Converting uniquely aligned to BAM and sorting"

                if [ "$pairedend" == "y" ] ; then
                        for prefix in $prefixes ; do
                                samtools merge \
                                -@ $threads \
                                $bamuniqdir/$prefix"-merged.bam" \
                                $samuniqdir/$prefix"-paired.sam" \
                                $samuniqdir/$prefix"-unpaired_R1.sam" \
                                $samuniqdir/$prefix"-unpaired_R2.sam" 2> /dev/null

                                samtools view \
                                -@ $threads \
                                -h \
                                -b \
                                -q 1 $bamuniqdir/$prefix"-merged.bam" | \
                                samtools sort \
                                -@ $threads \
                                -o $bamuniqdir/$prefix"-sorted.bam" 2> /dev/null

                                samtools index \
                                -b \
                                $bamuniqdir/$prefix"-sorted.bam"
                        done
                else
                        for prefix in $prefixes ; do
                                samtools view \
                                -@ $threads \
                                -h \
                                -b \
                                -q 1 $samuniqdir/$prefix"-unpaired_R1.sam" | \
                                samtools sort \
                                -@ $threads \
                                -o $bamuniqdir/$prefix"-sorted.bam" 2> /dev/null

                                samtools index \
                                -b \
                                $bamuniqdir/$prefix"-sorted.bam"
                        done
                fi

                echo "Done!"
        fi
fi

######################
# MMR
######################
# MMR will remap the multi aligned reads
# it will avoid the inflation of reads and
# allow us to compute the read coverage of
# repetitive regions on the genome

if [ ! -d $mmrdir ] ; then
        mkdir $mmrdir

        echo "Step 5: Adjusting position of multi aligned reads"
        # we did not check if the files are sorted because they were sorted before
        # -A is max no of valid pairs before not using pair modus
        # we want to maximize the number of valid pairs used (i.e. all)
        if [ "$pairedend" == "y" ] ; then
                for prefix in $prefixes ; do
                        mmr \
                        -t $threads \
                        --no-sort-check \
                        -S \
                        -b \
                        -p \
                        -i $maxfragsize \
                        -A 1000000000 \
                        -o $mmrdir/$prefix".bam" \
                        -R $readsize \
                        $bamdir/$prefix"-sorted.bam"

                        samtools sort \
                        -@ $threads \
                        -o $mmrdir/$prefix"-sorted.bam" \
                        $mmrdir/$prefix".bam"

                        samtools index \
                        -b \
                        $mmrdir/$prefix"-sorted.bam"
                done
        else
                for prefix in $prefixes ; do
                        mmr \
                        -t $threads \
                        --no-sort-check \
                        -S \
                        -b \
                        -o $mmrdir/$prefix".bam" \
                        -R $readsize \
                        $bamdir/$prefix"-sorted.bam"

                        samtools sort \
                        -@ $threads \
                        -o $mmrdir/$prefix"-sorted.bam" \
                        $mmrdir/$prefix".bam"

                        samtools index \
                        -b \
                        $mmrdir/$prefix"-sorted.bam"
                done
        fi


        echo "Done!"
fi

######################
# Creating coverage files (mult)
######################
# coverage files are tricky
# deeptools does a good job extending
# counts in between two sequenced reads
# but sometimes presents strange behavior
# for example: it reports 2 counts for each infered insert
# in this way, we divide by 2 the counts for paired-end coverage
#
# this discourage us to use gapped or softclipped alignments
# because deeptools also attribute different values to regions of gap and softclip
# and after dividing by two, it will present non-integer values
#
# note that bedgraph output is 0-based and a few conversions have to be done
# in order to use it as 1-based
#
# if the normalize option is provided
# the software will multiply the values by correction factors
# which are derived from the total number of alignments
#
# briefly, it calculates a correction factor for every library
# dividing the total number of alignments of the *library with
# the highest number of alignments* by the total number of alignments
# of the *library that is subject to correction*.
#
# in this step we also execute a Rscript (removeInconsistentPairs.R)
# to read BAM files and remove entries of paired-end alignments that
# lack consistency after being remapped by MMR. If the original RNA
# in the sample is too small and align to direct repeats on the genome,
# MMR may not report them properly, and the downstream tools like
# deepTools might not work as expected for these particular cases.
# a brief explanation of this MMR issue can be checked at
# https://github.com/ratschlab/mmr/issues/5

if [ ! -d $coveragedir ] ; then
        mkdir $coveragedir

        echo "Step 6.1: Creating coverage files"

        if [ ! -f $miscdir/readCounts.txt ] ; then
                touch $miscdir/readCounts.txt
        else
                rm $miscdir/readCounts.txt
                touch $miscdir/readCounts.txt
        fi

        if [ "$pairedend" == "y" ] ; then
                for prefix in $prefixes ; do
                        samtools view -@ $threads -h $mmrdir/$prefix"-sorted.bam" | \
                        grep "^@\|YT:Z:CP" | \
                        samtools sort -n -@ $threads -O SAM -o $coveragedir/$prefix"-paired.sam"

                        samtools view -@ $threads -h $mmrdir/$prefix"-sorted.bam" | \
                        grep "^@\|YT:Z:UU" | \
                        samtools view -@ $threads -b > $coveragedir/$prefix"-unpaired.bam"

                        pairedReadsAligned=`cat $coveragedir/$prefix"-paired.sam" | wc -l`
                        perThreads=`echo $(($pairedReadsAligned/$threads))`

                        if [ $(($perThreads%2)) == 1 ] ; then
                                perThreads=`echo $(($perThreads+1))`
                        fi

                        split -d -l $perThreads <(grep -v "^@" $coveragedir/$prefix"-paired.sam") $coveragedir/$prefix"-paired_tmp" --additional-suffix "_noheader.sam"
                        
                        for tmpfile in `seq -w 00 $(($threads-1))` ; do
                                cat <(samtools view -H $coveragedir/$prefix"-paired.sam") \
                                $coveragedir/$prefix"-paired_tmp"$tmpfile"_noheader.sam" | \
                                samtools view -@ $threads -b > $coveragedir/$prefix"-paired_tmp"$tmpfile".bam"

                                R -q -f ./removeInconsistentPairs.R --args $coveragedir/$prefix"-paired_tmp"$tmpfile".bam" > /dev/null &
                        done
                        wait

                        find $coveragedir -name "*tmp*adjusted.bam" -exec \
                        samtools merge -@ $threads -n $coveragedir/$prefix"-paired-adjusted.bam" {} + 

                        find $coveragedir -name "*tmp*failPairs.bam" -exec \
                        samtools merge -@ $threads -n $coveragedir/$prefix"-paired-failPairs.bam" {} +

                        samtools sort -@ $threads $coveragedir/$prefix"-paired-adjusted.bam" > $coveragedir/$prefix"-paired.bam"
                        failreads=`samtools view -@ $threads $coveragedir/$prefix"-paired-failPairs.bam" | wc -l`
                        failpairs=`echo "$failreads/2" | bc`
                        echo "$prefix has $failpairs pair(s) containing reads that do not match each other. Removing."

                        rm $coveragedir/$prefix"-paired.sam"
                        rm $coveragedir/$prefix*tmp*
                        rm $coveragedir/$prefix"-paired-adjusted.bam" $coveragedir/$prefix"-paired-failPairs.bam"

                        samtools index -b $coveragedir/$prefix"-paired.bam"
                        samtools index -b $coveragedir/$prefix"-unpaired.bam"

                        pairedReadsAligned=`samtools view -@ $threads $coveragedir/$prefix"-paired.bam" | wc -l`
                        pairsAligned=`echo "$pairedReadsAligned / 2" | bc`
                        unpairedReadsAligned=`samtools view -@ $threads $coveragedir/$prefix"-unpaired.bam" | wc -l`
                        totalAligned=`echo "$pairsAligned + $unpairedReadsAligned" | bc`

                        printf "$prefix\t$totalAligned\n" >> $miscdir/readCounts.txt
                done

                maxTotalAligned=`awk -v FS="\t" -v OFS="\t" -v max=0 '{if($2 > max){max=$2}}END{print max}' $miscdir/readCounts.txt`

                while IFS="	" read prefix totalAligned ; do

                        if [ "$normalize" == "y" ]; then
                                correctionFactor=`echo "scale=3; $maxTotalAligned / $totalAligned" | bc`
                        else
                                correctionFactor=1
                        fi

                        bamCoverage \
                        --numberOfProcessors $threads \
                        -b $coveragedir/$prefix"-paired.bam" \
                        -o $coveragedir/$prefix"-paired-fwd.bedgraph" \
                        --binSize 1 \
                        --extendReads $maxfragsize \
                        --outFileFormat bedgraph \
                        --filterRNAstrand reverse 2> /dev/null

                        bamCoverage \
                        --numberOfProcessors $threads \
                        -b $coveragedir/$prefix"-paired.bam" \
                        -o $coveragedir/$prefix"-paired-rev.bedgraph" \
                        --binSize 1 \
                        --extendReads $maxfragsize \
                        --outFileFormat bedgraph \
                        --filterRNAstrand forward 2> /dev/null

                        bamCoverage \
                        --numberOfProcessors $threads \
                        -b $coveragedir/$prefix"-unpaired.bam" \
                        -o $coveragedir/$prefix"-unpaired-fwd.bedgraph" \
                        --binSize 1 \
                        --outFileFormat bedgraph \
                        --filterRNAstrand reverse 2> /dev/null

                        bamCoverage \
                        --numberOfProcessors $threads \
                        -b $coveragedir/$prefix"-unpaired.bam" \
                        -o $coveragedir/$prefix"-unpaired-rev.bedgraph" \
                        --binSize 1 \
                        --outFileFormat bedgraph \
                        --filterRNAstrand forward 2> /dev/null

                        bedtools unionbedg \
                        -i \
                        $coveragedir/$prefix"-paired-fwd.bedgraph" \
                        $coveragedir/$prefix"-unpaired-fwd.bedgraph" | \
                        awk -v FS="\t" -v OFS="\t" -v correctionFactor=$correctionFactor \
                        '{print $1,$2,$3,(($4/2)+$5)*correctionFactor}' > $coveragedir/$prefix"-fwd.bedgraph"

                        bedtools unionbedg \
                        -i \
                        $coveragedir/$prefix"-paired-rev.bedgraph" \
                        $coveragedir/$prefix"-unpaired-rev.bedgraph" | \
                        awk -v FS="\t" -v OFS="\t" -v correctionFactor=$correctionFactor \
                        '{print $1,$2,$3,(($4/2)+$5)*correctionFactor}' > $coveragedir/$prefix"-rev.bedgraph"

                        bedtools unionbedg \
                        -i \
                        $coveragedir/$prefix"-paired-fwd.bedgraph" \
                        $coveragedir/$prefix"-unpaired-fwd.bedgraph" | \
                        awk -v FS="\t" -v OFS="\t" -v correctionFactor=$correctionFactor -v prefix=$prefix \
                        '{print $1,$2+1,$3,prefix"-fwd",(($4/2)+$5)*correctionFactor}' > $coveragedir/$prefix"-fwd.igv"

                        bedtools unionbedg \
                        -i \
                        $coveragedir/$prefix"-paired-rev.bedgraph" \
                        $coveragedir/$prefix"-unpaired-rev.bedgraph" | \
                        awk -v FS="\t" -v OFS="\t" -v correctionFactor=$correctionFactor -v prefix=$prefix \
                        '{print $1,$2+1,$3,prefix"-rev",(($4/2)+$5)*correctionFactor}' > $coveragedir/$prefix"-rev.igv"

                done < $miscdir/readCounts.txt
        else
                for prefix in $prefixes ; do
                        samtools view -@ $threads -h $mmrdir/$prefix"-sorted.bam" | \
                        grep "^@\|YT:Z:UU" | \
                        samtools view -@ $threads -b > $coveragedir/$prefix"-unpaired.bam"

                        samtools index -b $coveragedir/$prefix"-unpaired.bam"

                        totalAligned=`samtools view -@ $threads $coveragedir/$prefix"-unpaired.bam" | wc -l`

                        printf "$prefix\t$totalAligned\n" >> $miscdir/readCounts.txt
                done

                maxTotalAligned=`awk -v FS="\t" -v OFS="\t" -v max=0 '{if($2 > max){max=$2}}END{print max}' $miscdir/readCounts.txt`

                while IFS="	" read prefix totalAligned ; do

                        if [ "$normalize" == "y" ]; then
                                correctionFactor=`echo "scale=3; $maxTotalAligned / $totalAligned" | bc`
                        else
                                correctionFactor=1
                        fi

                        bamCoverage \
                        --numberOfProcessors $threads \
                        -b $coveragedir/$prefix"-unpaired.bam" \
                        -o $coveragedir/$prefix"-unpaired-fwd.bedgraph.tmp" \
                        --binSize 1 \
                        --outFileFormat bedgraph \
                        --filterRNAstrand reverse 2> /dev/null

                        bamCoverage \
                        --numberOfProcessors $threads \
                        -b $coveragedir/$prefix"-unpaired.bam" \
                        -o $coveragedir/$prefix"-unpaired-rev.bedgraph.tmp" \
                        --binSize 1 \
                        --outFileFormat bedgraph \
                        --filterRNAstrand forward 2> /dev/null

                        awk -v FS="\t" -v OFS="\t" -v correctionFactor=$correctionFactor \
                        '{print $1,$2,$3,$4*correctionFactor}' $coveragedir/$prefix"-unpaired-fwd.bedgraph.tmp" > $coveragedir/$prefix"-fwd.bedgraph"

                        awk -v FS="\t" -v OFS="\t" -v correctionFactor=$correctionFactor \
                        '{print $1,$2,$3,$4*correctionFactor}' $coveragedir/$prefix"-unpaired-rev.bedgraph.tmp" > $coveragedir/$prefix"-rev.bedgraph"

                        awk -v FS="\t" -v OFS="\t" -v correctionFactor=$correctionFactor -v prefix=$prefix \
                        '{print $1,$2+1,$3,prefix"-fwd",$4*correctionFactor}' $coveragedir/$prefix"-unpaired-fwd.bedgraph.tmp" > $coveragedir/$prefix"-fwd.igv"

                        awk -v FS="\t" -v OFS="\t" -v correctionFactor=$correctionFactor -v prefix=$prefix \
                        '{print $1,$2+1,$3,prefix"-rev",$4*correctionFactor}' $coveragedir/$prefix"-unpaired-rev.bedgraph.tmp" > $coveragedir/$prefix"-rev.igv"

                        rm $coveragedir/*.tmp

                done < $miscdir/readCounts.txt
        fi

        echo "Done!"
fi

######################
# Creating coverage files (uniq)
######################
# similar behavior described on the previous chunk
# but in this step there is no need to use removeInconsistentPairs.R
#
# normalization is done using only the number of uniquely aligned reads
# in the same aforementioned fashion
if [ "$uniqaln" == "y" ] ; then
        if [ ! -d $coverageuniqdir ] ; then
                mkdir $coverageuniqdir

                if [ ! -f $miscdir/readCountsUniq.txt ] ; then
                        touch $miscdir/readCountsUniq.txt
                else
                        rm $miscdir/readCountsUniq.txt
                        touch $miscdir/readCountsUniq.txt
                fi

                echo "Step 6.2: Creating coverage files for uniquely aligned reads"

                if [ "$pairedend" == "y" ] ; then
                        for prefix in $prefixes ; do
                                samtools view -@ $threads -h $bamuniqdir/$prefix"-sorted.bam" | \
                                grep "^@\|YT:Z:CP" | \
                                samtools view -@ $threads -b > $coverageuniqdir/$prefix"-paired.bam"

                                samtools view -@ $threads -h $bamuniqdir/$prefix"-sorted.bam" | \
                                grep "^@\|YT:Z:UU" | \
                                samtools view -@ $threads -b > $coverageuniqdir/$prefix"-unpaired.bam"

                                samtools index -b $coverageuniqdir/$prefix"-paired.bam"
                                samtools index -b $coverageuniqdir/$prefix"-unpaired.bam"

                                pairedReadsAligned=`samtools view -@ $threads $coverageuniqdir/$prefix"-paired.bam" | wc -l`
                                pairsAligned=`echo "$pairedReadsAligned / 2" | bc`
                                unpairedReadsAligned=`samtools view -@ $threads $coverageuniqdir/$prefix"-unpaired.bam" | wc -l`
                                totalAligned=`echo "$pairsAligned + $unpairedReadsAligned" | bc`

                                printf "$prefix\t$totalAligned\n" >> $miscdir/readCountsUniq.txt
                        done

                        maxTotalAligned=`awk -v FS="\t" -v OFS="\t" -v max=0 '{if($2 > max){max=$2}}END{print max}' $miscdir/readCountsUniq.txt`

                        while IFS="	" read prefix totalAligned ; do

                                if [ "$normalize" == "y" ]; then
                                        correctionFactor=`echo "scale=3; $maxTotalAligned / $totalAligned" | bc`
                                else
                                        correctionFactor=1
                                fi

                                bamCoverage \
                                --numberOfProcessors $threads \
                                -b $coverageuniqdir/$prefix"-paired.bam" \
                                -o $coverageuniqdir/$prefix"-paired-fwd.bedgraph" \
                                --binSize 1 \
                                --extendReads $maxfragsize \
                                --outFileFormat bedgraph \
                                --filterRNAstrand reverse 2> /dev/null

                                bamCoverage \
                                --numberOfProcessors $threads \
                                -b $coverageuniqdir/$prefix"-paired.bam" \
                                -o $coverageuniqdir/$prefix"-paired-rev.bedgraph" \
                                --binSize 1 \
                                --extendReads $maxfragsize \
                                --outFileFormat bedgraph \
                                --filterRNAstrand forward 2> /dev/null

                                bamCoverage \
                                --numberOfProcessors $threads \
                                -b $coverageuniqdir/$prefix"-unpaired.bam" \
                                -o $coverageuniqdir/$prefix"-unpaired-fwd.bedgraph" \
                                --binSize 1 \
                                --outFileFormat bedgraph \
                                --filterRNAstrand reverse 2> /dev/null

                                bamCoverage \
                                --numberOfProcessors $threads \
                                -b $coverageuniqdir/$prefix"-unpaired.bam" \
                                -o $coverageuniqdir/$prefix"-unpaired-rev.bedgraph" \
                                --binSize 1 \
                                --outFileFormat bedgraph \
                                --filterRNAstrand forward 2> /dev/null

                                bedtools unionbedg \
                                -i \
                                $coverageuniqdir/$prefix"-paired-fwd.bedgraph" \
                                $coverageuniqdir/$prefix"-unpaired-fwd.bedgraph" | \
                                awk -v FS="\t" -v OFS="\t" -v correctionFactor=$correctionFactor \
                                '{print $1,$2,$3,(($4/2)+$5)*correctionFactor}' > $coverageuniqdir/$prefix"-fwd.bedgraph"

                                bedtools unionbedg \
                                -i \
                                $coverageuniqdir/$prefix"-paired-rev.bedgraph" \
                                $coverageuniqdir/$prefix"-unpaired-rev.bedgraph" | \
                                awk -v FS="\t" -v OFS="\t" -v correctionFactor=$correctionFactor \
                                '{print $1,$2,$3,(($4/2)+$5)*correctionFactor}' > $coverageuniqdir/$prefix"-rev.bedgraph"

                                bedtools unionbedg \
                                -i \
                                $coverageuniqdir/$prefix"-paired-fwd.bedgraph" \
                                $coverageuniqdir/$prefix"-unpaired-fwd.bedgraph" | \
                                awk -v FS="\t" -v OFS="\t" -v correctionFactor=$correctionFactor -v prefix=$prefix \
                                '{print $1,$2+1,$3,prefix"-fwd",(($4/2)+$5)*correctionFactor}' > $coverageuniqdir/$prefix"-fwd.igv"

                                bedtools unionbedg \
                                -i \
                                $coverageuniqdir/$prefix"-paired-rev.bedgraph" \
                                $coverageuniqdir/$prefix"-unpaired-rev.bedgraph" | \
                                awk -v FS="\t" -v OFS="\t" -v correctionFactor=$correctionFactor -v prefix=$prefix \
                                '{print $1,$2+1,$3,prefix"-rev",(($4/2)+$5)*correctionFactor}' > $coverageuniqdir/$prefix"-rev.igv"

                        done < $miscdir/readCountsUniq.txt
                else
                        for prefix in $prefixes ; do
                                samtools view -@ $threads -h $bamuniqdir/$prefix"-sorted.bam" | \
                                grep "^@\|YT:Z:UU" | \
                                samtools view -@ $threads -b > $coverageuniqdir/$prefix"-unpaired.bam"

                                samtools index -b $coverageuniqdir/$prefix"-unpaired.bam"

                                totalAligned=`samtools view -@ $threads $coverageuniqdir/$prefix"-unpaired.bam" | wc -l`

                                printf "$prefix\t$totalAligned\n" >> $miscdir/readCountsUniq.txt
                        done

                        maxTotalAligned=`awk -v FS="\t" -v OFS="\t" -v max=0 '{if($2 > max){max=$2}}END{print max}' $miscdir/readCountsUniq.txt`

                        while IFS="	" read prefix totalAligned ; do
                                if [ "$normalize" == "y" ]; then
                                        correctionFactor=`echo "scale=3; $maxTotalAligned / $totalAligned" | bc`
                                else
                                        correctionFactor=1
                                fi

                                bamCoverage \
                                --numberOfProcessors $threads \
                                -b $coverageuniqdir/$prefix"-unpaired.bam" \
                                -o $coverageuniqdir/$prefix"-unpaired-fwd.bedgraph.tmp" \
                                --binSize 1 \
                                --outFileFormat bedgraph \
                                --filterRNAstrand reverse 2> /dev/null

                                bamCoverage \
                                --numberOfProcessors $threads \
                                -b $coverageuniqdir/$prefix"-unpaired.bam" \
                                -o $coverageuniqdir/$prefix"-unpaired-rev.bedgraph.tmp" \
                                --binSize 1 \
                                --outFileFormat bedgraph \
                                --filterRNAstrand forward 2> /dev/null

                                awk -v FS="\t" -v OFS="\t" -v correctionFactor=$correctionFactor \
                                '{print $1,$2,$3,$4*correctionFactor}' $coverageuniqdir/$prefix"-unpaired-fwd.bedgraph.tmp" > $coverageuniqdir/$prefix"-fwd.bedgraph"

                                awk -v FS="\t" -v OFS="\t" -v correctionFactor=$correctionFactor \
                                '{print $1,$2,$3,$4*correctionFactor}' $coverageuniqdir/$prefix"-unpaired-rev.bedgraph.tmp" > $coverageuniqdir/$prefix"-rev.bedgraph"

                                awk -v FS="\t" -v OFS="\t" -v correctionFactor=$correctionFactor -v prefix=$prefix \
                                '{print $1,$2+1,$3,prefix"-fwd",$4*correctionFactor}' $coverageuniqdir/$prefix"-unpaired-fwd.bedgraph.tmp" > $coverageuniqdir/$prefix"-fwd.igv"

                                awk -v FS="\t" -v OFS="\t" -v correctionFactor=$correctionFactor -v prefix=$prefix \
                                '{print $1,$2+1,$3,prefix"-rev",$4*correctionFactor}' $coverageuniqdir/$prefix"-unpaired-rev.bedgraph.tmp" > $coverageuniqdir/$prefix"-rev.igv"

                        rm $coverageuniqdir/*.tmp

                        done < $miscdir/readCountsUniq.txt
                fi

                echo "Done!"
        fi
fi

######################
# Creating TSSAR input files
######################
# this step will create BAM files containing only aligned
# R1 from pairs and unpaired R1 (orphans)
# they are useful to compute enrichment between two
# enriched libraries. eg. TEX+ vs TEX- which
# allow to find primary transcripts
if [ ! -d $tssarinputdir ] ; then
        mkdir $tssarinputdir

        echo "Step 7: Creating TSSAR input files"

        replicons=`grep ">" $miscdir/$spp".fa" | sed 's/^>//;s/ .*$//'`

        if [ "$pairedend" == "y" ] ; then
                for prefix in $prefixes ; do
                        samtools view -@ $threads -h -b -f 0x40 $coveragedir/$prefix"-paired.bam" > $tssarinputdir/$prefix"-paired-R1.bam"

                        samtools view -@ $threads -h -F 0x10 $coveragedir/$prefix"-unpaired.bam" | \
                        awk -v OFS="\t" -v FS="\t" '{if(/^@/){print}else{if($1 !~ /_R2$/){print}}}' | grep "^@\|XS:A:+" | \
                        samtools view -@ $threads -h -b > $tssarinputdir/$prefix"-unpaired-R1-flag0.bam"

                        samtools view -@ $threads -h -f 0x10 $coveragedir/$prefix"-unpaired.bam" | \
                        awk -v OFS="\t" -v FS="\t" '{if(/^@/){print}else{if($1 !~ /_R2$/){print}}}' | grep "^@\|XS:A:-" | \
                        samtools view -@ $threads -h -b > $tssarinputdir/$prefix"-unpaired-R1-flag16.bam"

                        samtools merge \
                        -@ $threads \
                        $tssarinputdir/$prefix"-tssar-input.bam" \
                        $tssarinputdir/$prefix"-paired-R1.bam" \
                        $tssarinputdir/$prefix"-unpaired-R1-flag0.bam" \
                        $tssarinputdir/$prefix"-unpaired-R1-flag16.bam" 2> /dev/null

                        for replicon in $replicons ; do
                                samtools view -@ $threads -h $tssarinputdir/$prefix"-tssar-input.bam" | \
                                grep "^@HD	VN:1.0	SO:coordinate\|^@SQ	SN:$replicon\|^@PG	ID:hisat2	PN:hisat2\|	$replicon	" | \
                                samtools view -@ $threads -b > $tssarinputdir/$prefix"-tssar-input-"$replicon".bam"
                        done
                done
        else
                for prefix in $prefixes ; do
                        samtools view -@ $threads -h -F 0x10 $coveragedir/$prefix"-unpaired.bam" | \
                        awk -v OFS="\t" -v FS="\t" '{if(/^@/){print}else{if($1 !~ /_R2$/){print}}}' | grep "^@\|XS:A:+" | \
                        samtools view -@ $threads -h -b > $tssarinputdir/$prefix"-unpaired-R1-flag0.bam"

                        samtools view -@ $threads -h -f 0x10 $coveragedir/$prefix"-unpaired.bam" | \
                        awk -v OFS="\t" -v FS="\t" '{if(/^@/){print}else{if($1 !~ /_R2$/){print}}}' | grep "^@\|XS:A:-" | \
                        samtools view -@ $threads -h -b > $tssarinputdir/$prefix"-unpaired-R1-flag16.bam"

                        samtools merge \
                        -@ $threads \
                        $tssarinputdir/$prefix"-tssar-input.bam" \
                        $tssarinputdir/$prefix"-unpaired-R1-flag0.bam" \
                        $tssarinputdir/$prefix"-unpaired-R1-flag16.bam" 2> /dev/null

                        for replicon in $replicons ; do
                                samtools view -@ $threads -h $tssarinputdir/$prefix"-tssar-input.bam" | \
                                grep "^@HD	VN:1.0	SO:coordinate\|^@SQ	SN:$replicon\|^@PG	ID:hisat2	PN:hisat2\|	$replicon	" | \
                                samtools view -@ $threads -b > $tssarinputdir/$prefix"-tssar-input-"$replicon".bam"
                        done
                done
        fi

        echo "Done!"
fi

######################
# Creating 5 prime profile files
######################
# using only R1 files (from pairs and also unpaired)
# we can compute the start of RNA fragments in the sample
# this could help to infer transcription start sites or persistent
# transcripts (resistant to degradation) when combined to 3 prime profiles
if [ ! -d $fiveprimedir ] ; then
        mkdir $fiveprimedir

        echo "Step 8: Creating five prime profiling files"

        if [ "$pairedend" == "y" ] ; then
                for prefix in $prefixes ; do
                        bedtools genomecov -5 -strand + -bga -ibam $tssarinputdir/$prefix"-paired-R1.bam" > $fiveprimedir/$prefix"-5primeprofile-paired-fwd.bedgraph"

                        bedtools genomecov -5 -strand - -bga -ibam $tssarinputdir/$prefix"-paired-R1.bam" > $fiveprimedir/$prefix"-5primeprofile-paired-rev.bedgraph"

                        bedtools genomecov -5 -strand + -bga -ibam $tssarinputdir/$prefix"-unpaired-R1-flag0.bam" > $fiveprimedir/$prefix"-5primeprofile-unpaired-fwd.bedgraph"

                        bedtools genomecov -5 -strand - -bga -ibam $tssarinputdir/$prefix"-unpaired-R1-flag16.bam" > $fiveprimedir/$prefix"-5primeprofile-unpaired-rev.bedgraph"
                done

                maxTotalAligned=`awk -v FS="\t" -v OFS="\t" -v max=0 '{if($2 > max){max=$2}}END{print max}' $miscdir/readCounts.txt`

                while IFS="	" read prefix totalAligned ; do

                        if [ "$normalize" == "y" ]; then
                                correctionFactor=`echo "scale=3; $maxTotalAligned / $totalAligned" | bc`
                        else
                                correctionFactor=1
                        fi

                        bedtools unionbedg \
                        -i \
                        $fiveprimedir/$prefix"-5primeprofile-paired-fwd.bedgraph" \
                        $fiveprimedir/$prefix"-5primeprofile-unpaired-fwd.bedgraph" | \
                        awk -v FS="\t" -v OFS="\t" -v correctionFactor=$correctionFactor \
                        '{print $1,$2,$3,($4+$5)*correctionFactor}' > $fiveprimedir/$prefix"-5primeprofile-fwd.bedgraph"

                        bedtools unionbedg \
                        -i \
                        $fiveprimedir/$prefix"-5primeprofile-paired-rev.bedgraph" \
                        $fiveprimedir/$prefix"-5primeprofile-unpaired-rev.bedgraph" | \
                        awk -v FS="\t" -v OFS="\t" -v correctionFactor=$correctionFactor \
                        '{print $1,$2,$3,($4+$5)*correctionFactor}' > $fiveprimedir/$prefix"-5primeprofile-rev.bedgraph"

                done < $miscdir/readCounts.txt

        else
                maxTotalAligned=`awk -v FS="\t" -v OFS="\t" -v max=0 '{if($2 > max){max=$2}}END{print max}' $miscdir/readCounts.txt`

                while IFS="	" read prefix totalAligned ; do

                        if [ "$normalize" == "y" ]; then
                                correctionFactor=`echo "scale=3; $maxTotalAligned / $totalAligned" | bc`
                        else
                                correctionFactor=1
                        fi

                        bedtools genomecov -5 -strand + -bga -ibam $tssarinputdir/$prefix"-unpaired-R1-flag0.bam" | \
                        awk -v FS="\t" -v OFS="\t" -v correctionFactor=$correctionFactor \
                        '{print $1,$2,$3,$4*correctionFactor}' > $fiveprimedir/$prefix"-5primeprofile-fwd.bedgraph"

                        bedtools genomecov -5 -strand - -bga -ibam $tssarinputdir/$prefix"-unpaired-R1-flag16.bam" | \
                        awk -v FS="\t" -v OFS="\t" -v correctionFactor=$correctionFactor \
                        '{print $1,$2,$3,$4*correctionFactor}' > $fiveprimedir/$prefix"-5primeprofile-rev.bedgraph"

                done < $miscdir/readCounts.txt
        fi

        echo "Done!"
fi

######################
# Creating 3 prime profile files
######################
# using only R2 files (from pairs and also unpaired orphans)
# we can compute the end of RNA fragments in the sample
# this could help to infer the end of fragments.
#
# the end of a fragment is the five prime position of a R2
# read on the opposite strand (note that R2 comes from an artificial
# complementary DNA synthesized based on the real RNA molecule)
if [ ! -d $threeprimedir ] ; then
        mkdir $threeprimedir

        echo "Step 9: Creating three prime profiling files"

        if [ "$pairedend" == "y" ] ; then
                for prefix in $prefixes ; do
                        samtools view -@ $threads -h -b -f 0x80 $coveragedir/$prefix"-paired.bam" | \
                        bedtools genomecov -5 -strand - -bga -ibam stdin > $threeprimedir/$prefix"-3primeprofile-paired-fwd.bedgraph"

                        samtools view -@ $threads -h -b -f 0x80 $coveragedir/$prefix"-paired.bam" | \
                        bedtools genomecov -5 -strand + -bga -ibam stdin > $threeprimedir/$prefix"-3primeprofile-paired-rev.bedgraph"

                        samtools view -@ $threads -h -F 0x10 $coveragedir/$prefix"-unpaired.bam" | \
                        awk -v OFS="\t" -v FS="\t" '{if(/^@/){print}else{if($1 ~ /_R2$/){print}}}' | grep "^@\|XS:A:+" | \
                        samtools view -@ $threads -h -b | \
                        bedtools genomecov -3 -strand + -bga -ibam stdin > $threeprimedir/$prefix"-3primeprofile-unpaired-fwd.bedgraph"

                        samtools view -@ $threads -h -f 0x10 $coveragedir/$prefix"-unpaired.bam" | \
                        awk -v OFS="\t" -v FS="\t" '{if(/^@/){print}else{if($1 ~ /_R2$/){print}}}' | grep "^@\|XS:A:-" | \
                        samtools view -@ $threads -h -b | \
                        bedtools genomecov -3 -strand - -bga -ibam stdin > $threeprimedir/$prefix"-3primeprofile-unpaired-rev.bedgraph"
                done

                maxTotalAligned=`awk -v FS="\t" -v OFS="\t" -v max=0 '{if($2 > max){max=$2}}END{print max}' $miscdir/readCounts.txt`

                while IFS="	" read prefix totalAligned ; do

                        if [ "$normalize" == "y" ]; then
                                correctionFactor=`echo "scale=3; $maxTotalAligned / $totalAligned" | bc`
                        else
                                correctionFactor=1
                        fi

                        bedtools unionbedg \
                        -i \
                        $threeprimedir/$prefix"-3primeprofile-paired-fwd.bedgraph" \
                        $threeprimedir/$prefix"-3primeprofile-unpaired-fwd.bedgraph" | \
                        awk -v FS="\t" -v OFS="\t" -v correctionFactor=$correctionFactor \
                        '{print $1,$2,$3,($4+$5)*correctionFactor}' > $threeprimedir/$prefix"-3primeprofile-fwd.bedgraph"

                        bedtools unionbedg \
                        -i \
                        $threeprimedir/$prefix"-3primeprofile-paired-rev.bedgraph" \
                        $threeprimedir/$prefix"-3primeprofile-unpaired-rev.bedgraph" | \
                        awk -v FS="\t" -v OFS="\t" -v correctionFactor=$correctionFactor \
                        '{print $1,$2,$3,($4+$5)*correctionFactor}' > $threeprimedir/$prefix"-3primeprofile-rev.bedgraph"

                done < $miscdir/readCounts.txt
        else

                maxTotalAligned=`awk -v FS="\t" -v OFS="\t" -v max=0 '{if($2 > max){max=$2}}END{print max}' $miscdir/readCounts.txt`

                while IFS="	" read prefix totalAligned ; do

                        if [ "$normalize" == "y" ]; then
                                correctionFactor=`echo "scale=3; $maxTotalAligned / $totalAligned" | bc`
                        else
                                correctionFactor=1
                        fi

                samtools view -@ $threads -h -F 0x10 $coveragedir/$prefix"-unpaired.bam" | \
                grep "^@\|XS:A:+" | \
                samtools view -@ $threads -h -b | \
                bedtools genomecov -3 -strand + -bga -ibam stdin | \
                awk -v FS="\t" -v OFS="\t" -v correctionFactor=$correctionFactor \
                '{print $1,$2,$3,$4*correctionFactor}' > $threeprimedir/$prefix"-3primeprofile-unpaired-fwd.bedgraph"

                samtools view -@ $threads -h -f 0x10 $coveragedir/$prefix"-unpaired.bam" | \
                grep "^@\|XS:A:-" | \
                samtools view -@ $threads -h -b | \
                bedtools genomecov -3 -strand - -bga -ibam stdin | \
                awk -v FS="\t" -v OFS="\t" -v correctionFactor=$correctionFactor \
                '{print $1,$2,$3,$4*correctionFactor}' > $threeprimedir/$prefix"-3primeprofile-unpaired-rev.bedgraph"

                done < $miscdir/readCounts.txt
        fi

        echo "Done!"
fi

# continuing the help if-else statement
else
echo '
####################################
# USAGE MANUAL
####################################

bash frtc.sh <threads> <maxfragsize> <read_size> <spp> <url>

threads [INT]:     number of threads to be passed to nested programs
                   maximum value: 99

maxfragsize [INT]: maximum insert size from the leftmost end to
                   the rightmost of an paired-end alignment.
                   This parameter is ignored if single-end libraries
                   are being processed.

read_size [INT]:   approximate size of reads. For an Illumina paired-end
                   experiment of 150x2 cycles, use 150.

spp [VARCHAR]:     prefix for genome and annotation files

url [VARCHAR]:     URL pointing to the genome file in NCBI RefSeq FTP Server.
                   One may find the link for RefSeq directory on the top right
                   corner of the NCBI Assembly page for a given genome.

e.g.:
bash frtc.sh 6 1000 150 Hsalinarum ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/805/GCF_000006805.1_ASM680v1/GCF_000006805.1_ASM680v1_genomic.fna.gz

this command line will run frtc.sh with six threads,
with a maximum insert size of 1000, with read size of 150 nt,
creating Hsalinarum.fa and Hsalinarum.gff files under misc
directory to store genome and annotation, respectively.

there is a directory tree prerequisite to run frtc.
an adequate directory tree to run frtc must look like

yourDirectory
├── removeInconsistentPairs.R      # provided
├── frtc.sh                        # provided
├── misc                           # notProvided
│   └── adap.fa
└── raw                            # notProvided
    ├── S1_R1.fastq.gz
    ├── S1_R2.fastq.gz
    ├── S2_R1.fastq.gz
    ├── S2_R2.fastq.gz
    ├── S3_R1.fastq.gz
    ├── S3_R2.fastq.gz
    ├── S4_R1.fastq.gz
    └── S4_R2.fastq.gz

raw file names must be ended with R1 or R2 and
must have fastq.gz extension (i.e. gzip compressed)
if the libraries are paired-end
e.g. S1_R1.fastq.gz
     S1_R2.fastq.gz

if the libraries are single-end
e.g. S1_R1.fastq.gz

single-end or paired-end mode can be set in frtc.sh.
this script was written using paired-end as default mode
but you are able to change to single-end mode by setting
pairedend="n" on "CUSTOM VARIABLES" section.

Also, if you are processing paired-end libraries, check if
R1 and R2 corresponding reads have the same name. For example,
if the first entry in R1 file has the header "@SRR9999999.1",
the first entry in R2 file must have the same header "@SRR9999999.1".
Everything should be the same before the first space character.
NCBI SRA commonly provides files that present headers like "@SRR9999999.1.1"
for R1 and "@SRR9999999.1.2" for R2, which may cause problems during
the multi-mapper resolution step of this pipeline.

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
maybe newer versions are compatible, but these were the ones used
to build and test this script:

trimmomatic v0.39 (must be located @ /opt/Trimmomatic-0.39/trimmomatic-0.39.jar)
hisat2 v2.1.0 (@ PATH)
curl v7.58.0 (@ PATH)
samtools v1.7 (@ PATH)
mmr default version (@ PATH)
deeptools v3.1.1 (@PATH)
bedtools v2.26.0 (also tested with v2.21.0) (@ PATH)

all the prerequisites will be checked before running

**Warning**  

* rbamtools is a required R package.
This package lost its support recently, so the best
workaround now is to install it using a copy of the source
code stored at an unnofficial github repo. To do this,
you should run the following commands on your R console:

library("devtools") ; devtools::install_github("cran/rbamtools")

###################################
# PURPOSE AND GENERAL INFO
###################################

this script will perform every task from trimming
prokaryotic Illumina RNA-Seq raw data (fastq.gz files)
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


>>>>>>>>>>>>>>>>>>>>>>>>>>>>WARNING<<<<<<<<<<<<<<<<<<<<<<<<<<<<
be aware this feature doesnt work for dUTP library preparations;
be aware that uniq mode is not supported for
fiveprime, threeprime and TSSAR input results

This pipeline was conceived to process prokaryotic RNA-Seq data,
so by default spliced alignments are not allowed. Feel free to
adapt and test it using eukaryotic datasets and let me know if
everything is working well.
---------------------------------------------------------------

furthermore, it will filter bam files in order to keep
only the R1 files, for they are useful to find TSS
in dRNA-Seq experiments (e.g. TEXminus vs. TEXplus experiments)
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
5. adjusting position of multi-mappers using MMR and removing pairs that do not match each other (custom Rscript; see https://github.com/ratschlab/mmr/issues/5)
6. creating coverage files (bedgraph and igv format) using deepTools
7. creating TSSAR input BAM files
8. creating five prime profiling coverage (bedgraph format) using bedtools
9. creating three prime profiling coverage (bedgraph format) using bedtools

###################################
# CITATION
###################################

If you make use of frtc in your research, please cite:  

ten-Caten, F., Vencio, R. Z., Lorenzetti, A. P., Zaramela, L. S., Santana, A. C., & Koide, T. (2018). Internal RNAs overlapping coding sequences can drive the production of alternative proteins in archaea. RNA Biology.
http://dx.doi.org/10.1080/15476286.2018.1509661
'

fi
