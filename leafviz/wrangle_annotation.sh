#!/bin/sh

if [[ $1 == "" || $1 == "--help" || $1 == "-h" ]];then
	echo "
wrangle_annotation.sh
-------------

finds introns, exons and splice sites in a GTF file

usage 
sh wrangle_annotation.sh <GTF> <outfile_prefix>

GTF can be gzipped
outfile_prefix can be file name or full path

WARNING: this has only been tested on GENCODE V26 mouse (mm10) and human (hg19/hg38)
"
exit 0 
fi

if [[ $1 =~ gz$ ]];then
	command=zless
else
	command=cat
fi

# $1 is the file
# $2 is the prefix code
CODE=$2

# This is the position of the gene_name in the last gtf column
SPECIES=3

#take a list of exons in a GTF and return the introns
echo finding introns

$command $1 | gawk -F '\t' \
		-v species=$SPECIES  '
	BEGIN{ OFS = "\t" } 
	$3 == "gene" {
		split( $9, a, ";" ) # split off the metadata
		split( a[species], b, " ") # split the gene name column
		GENE=b[2]
		split( a[1], b, " ")
		GENEID=b[2] 
		#print GENE > "/dev/stderr" # print to stderr
	}

	$3 == "transcript" {
		split( $9, a, ";" )
		split( a[2], b, " ")
		TRANSCRIPT=b[2] # this is the transcript name
		split( a[10], c, " ")
		TAG=c[2]
		if ( TAG == ""){
			TAG="NA"
		}
		split( a[3], d, " ")
		TRANSCRIPTTYPE=d[2]
		INTRONCOUNT=0
		INTRONSTART=0; INTRONEND=0 # reset the coordinates 
	}
	$3 == "exon" {
		if ( INTRONCOUNT != 0 ){ # cannot happen on the first exon
			if ( $7 == "+" ){ positive strand
				INTRONEND=$4
				print $1, INTRONSTART, INTRONEND, GENE, GENEID, $7, TRANSCRIPT, TRANSCRIPTTYPE, INTRONCOUNT,TAG
			}
			if ( $7 == "-" ){ # negative strand
				INTRONSTART=$5
				print $1, INTRONSTART, INTRONEND, GENE, GENEID, $7, TRANSCRIPT, INTRONCOUNT, TRANSCRIPTTYPE, TAG
			}
		}
		INTRONCOUNT+=1
		if ( $7 == "+" ){
			INTRONSTART=$5 # begins on the first exon
		}
		if ( $7 == "-" ){
			INTRONEND=$4 # begins on the first exon
		}
	}' | gzip > $CODE"_all_introns.bed.gz"


# Take a list of exons and return all the 5' splice sites to one file and all the 3' splice sites to another.
echo finding splice sites

$command $1 | gawk -F '\t' \
	-v CODE=$CODE \
	-v species=$SPECIES '
	BEGIN{ OFS = "\t" } 
	$3 == "gene" {
		split( $9, a, ";" )
		split( a[species], b, " ")
		GENE=b[2]
		split( a[1], b, " ")
		GENEID=b[2]  
	}

	$3 == "transcript" {
		split( $9, a, ";" )
		split( a[2], b, " ")
		TRANSCRIPT=b[2] # this is the transcript name
		split( a[10], c, " ")
		TAG=c[2] #
		if ( TAG == ""){
			TAG="NA"
		}
		split( a[3], d, " ")
		TRANSCRIPTTYPE=d[2]
		INTRONCOUNT=0
		INTRONSTART=0; INTRONEND=0 # reset the coordinates 
	}
	$3 == "exon" {
		if ( INTRONCOUNT != 0 ){ # cannot happen on the first exon
			if ( $7 == "+" ){ positive strand
				INTRONEND=$4				
			}
			if ( $7 == "-" ){ # negative strand
				INTRONSTART=$5
			}
				print $1, INTRONSTART, INTRONSTART + 1, GENE, GENEID, $7, TRANSCRIPT, INTRONCOUNT, TRANSCRIPTTYPE, TAG > CODE"_fiveprime.bed"
				print $1, INTRONEND, INTRONEND + 1, GENE, GENEID, $7, TRANSCRIPT, INTRONCOUNT, TRANSCRIPTTYPE, TAG > CODE"_threeprime.bed"
		}
		INTRONCOUNT+=1
		if ( $7 == "+" ){
			INTRONSTART=$5 # begins on the first exon
		}
		if ( $7 == "-" ){
			INTRONEND=$4 # begins on the first exon
		}
	}' 

gzip ${CODE}_fiveprime.bed
gzip ${CODE}_threeprime.bed

echo creating exon list
# TODO: also create the exon lists for the graphing step.
$command $1 | gawk -F '\t'   \
		-v species=$SPECIES ' # exons have the gene name in a different column to the gene entries. 
BEGIN{
	print "chr start end strand gene_name"
			} $3 == "exon" {
				split( $NF, a, "; " );
				split( a[species + 1], b, "\"" );
				print $1,$4,$5,$7,b[2]
			}' | gzip > ${CODE}_all_exons.txt.gz
