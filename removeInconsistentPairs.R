# alorenzetti 20180319
# this script will remove paired-end reads that
# make no sense after being remapped by MMR
# we detected a bug in this application (MMR)
# some pairs doesnt have mateposition matching 
# the entry for its mate. also the insert sizes
# might be messed up

# loading libs
library("rbamtools", quietly = T)

# getting args
args = commandArgs(trailingOnly = T)

# reading bam file and getting info about alignment refs
bamnsortedfile = args[1]
bamnsorted = bamReader(filename = bamnsortedfile)
df = getRefData(bamnsorted)

# parsing file names for output purposes
prefix = sub(".bam", "", bamnsortedfile)

# creating output file names
out = paste0(prefix, "-adjusted.bam")
outfail = paste0(prefix, "-failPairs.bam")

# we need to create bamheader from scracth
# the function to read original header is not working and I dont know why
# setting header lines
headl = new("headerLine")
setVal(headl, "SO", "queryname")
setVal(headl, "VN", "1.0")

# setting reference seq dict header lines
dict = new("refSeqDict")
addseqtoheader = function(sn, ln){addSeq(dict, SN=sn, LN=ln)}
for(i in 1:dim(df)[1]){
  addSeq(dict, SN=df[i,2], LN=df[i,3])
}

# setting program header lines
prog <- new("headerProgram")
setVal(prog, "ID", "hisat2")
setVal(prog, "PN", "hisat2")

# assigning objects to header object
bamheader = bamHeaderText(head=headl, dict=dict, prog=prog)
bamheader = bamHeader(bamheader)

# opening file to write header and alignments
writer = bamWriter(bamheader, out)

# opening file to write header and alignments that present something strange
fail = bamWriter(bamheader, outfail)

# setting counter
ntotal = 0

# getting aligns
rewind(bamnsorted)
align1 = getNextAlign(bamnsorted)

while(!is.null(align1)){

  align2 = getNextAlign(bamnsorted)
  name1eqname2 = name(align1) == name(align2)
  
  if(name1eqname2){

	  pos1eqmpos2 = position(align1) == matePosition(align2)
	  mpos1eqpos2 = position(align2) == matePosition(align1)
	  
	  if(position(align1) < position(align2)){
	    insertlen1 = position(align2) + nchar(alignSeq(align2)) - position(align1)
	    insertlen2 = insertlen1 * -1
	  } else if(position(align1) > position(align2)){
	    insertlen1 = (position(align1) + nchar(alignSeq(align1)) - position(align2)) * -1
	    insertlen2 = insertlen1 * -1
	  } else {
	    if(nchar(alignSeq(align1)) < nchar(alignSeq(align2))){
	      insertlen1 = nchar(alignSeq(align2))
	      insertlen2 = insertlen1 * -1
	    } else if(nchar(alignSeq(align1)) > nchar(alignSeq(align2))){
	      insertlen1 = nchar(alignSeq(align1)) * -1
	      insertlen2 = insertlen1 * -1
	    } else {
	      insertlen1 = nchar(alignSeq(align1)) * -1
	      insertlen2 = insertlen1
	    }
	  }
	  
	  insertlengthok1 = insertlen1 == insertSize(align1)
	  insertlengthok2 = insertlen2 == insertSize(align2)

          # check if a big insert might be actually a small rna
          # we would rather remove these big unreal inserts 
          seq1eqseq2 = alignSeq(align1) == alignSeq(align2)
          seq1leinsertlength1 = nchar(alignSeq(align1)) < abs(insertSize(align1))
          shouldbekept = !(seq1eqseq2 & seq1leinsertlength1)
	  
	  if((name1eqname2 & pos1eqmpos2 & mpos1eqpos2 & insertlengthok1 & insertlengthok2 & shouldbekept) == T){
	    bamSave(object = writer, value = align1, refid = refID(align1))
	    bamSave(object = writer, value = align2, refid = refID(align2))
	  } else {
	    ntotal = ntotal + 1
	    bamSave(object = fail, value = align1, refid = refID(align1))
	    bamSave(object = fail, value = align2, refid = refID(align2))
	  }

    align1 = getNextAlign(bamnsorted)

  } else {
    align1 = align2
  }
  
}

# closing bamfile
bamClose(writer)

# closing bamfile fail
bamClose(fail)
