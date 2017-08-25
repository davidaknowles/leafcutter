#!/usr/bin/env perl

use strict;
use warnings;

use Data::Dumper;
use Carp;
use File::Basename;
use Getopt::Long;


my $prog = basename ($0);
my $separateBed = 0;
my $printUniqOnly = 0;
my $verbose = 0;
my $useRNAStrand = 0; # use the strand of the RNA instead of the read

GetOptions (
	"u|uniq"=>\$printUniqOnly,
	"r|use-RNA-strand"=>\$useRNAStrand,
#	"s|separate-bed"=>\$separateBed, 
	"v|verbose"=>\$verbose);

if (@ARGV != 2 && @ARGV != 3)
{
	print STDERR "Convert OLego SAM format to BED format (for both paired-end and single-end data)\n";
	print STDERR "Usage: $prog [options] <in.sam> <out1.bed> [out2.bed]\n";
	print STDERR " <in.sam> : gzip compressed input file with .gz extension is allowed\n";
	print STDERR " <out1.bed> [out2.bed]: specify both out1.bed and out2.bed to output results of PE data to separate BED files.\n";
	print STDERR " You can also use - to specify STDIN for input or STDOUT for output\n\n";
	print STDERR "options:\n";
	print STDERR "-u,--uniq:            print uniquely mapped reads only\n";
	print STDERR "-r,--use-RNA-strand:  force to use the strand of the RNA based on the XS tag \n";
	print STDERR "-v,--verbose:         verbose\n";
	exit (1);
}

my ($inSAMFile, $outBedFile) = @ARGV;
my $outBedFile2 = "";
if (@ARGV == 3)
{
	$outBedFile2 = $ARGV[2];
	die "Please specify different names for the seperate bed files.\n" if ($outBedFile eq $outBedFile2);
	$separateBed = 1;
}


my ($fin, $fout, $fout2);

if ( $inSAMFile eq "-")
{
    $fin = *STDIN;
}
else
{
	if ($inSAMFile =~/\.gz$/)
	{
		open ($fin, "gunzip -c $inSAMFile | ") || Carp::croak "cannot open file $inSAMFile to read\n";
	}
	else
	{
    	open ($fin, "<$inSAMFile") || Carp::croak "cannot open file $inSAMFile to read\n";
	}
}
if ( $outBedFile eq "-")
{
     $fout = *STDOUT;
}
else
{
    open ($fout, ">$outBedFile") || Carp::croak "cannot open file $outBedFile to write\n";
}

if ($separateBed)
{
    if ($outBedFile2 eq  "-")
    {
	$fout2 = *STDOUT;
    }
    else
    {
	open ($fout2, ">$outBedFile2") || Carp::croak "cannot open file $outBedFile2 to write\n";
    }
}


my $i = 0;
my $found = 0;

while (my $line = <$fin>)
{
	chomp $line;

	next if $line=~/^\s*$/;
	next if $line=~/^\@/;

	print STDERR "$i ...\n" if $verbose && $i % 50000 == 0;
	$i++;

	my $sam = lineToSam ($line);
	my $bed = samToBed ($sam, $useRNAStrand);
	next unless $bed; #no alignment

	my $flagInfo = $bed->{"flagInfo"};
	Carp::croak "inconsistency in specifying PE or SE data\n" if ($flagInfo->{'PE'}==0 &&  $separateBed == 1);

	my $read1_or_2 = $flagInfo->{'read_1_or_2'};

	my $uniq = 0;
	$uniq = 1 if $sam->{"TAGS"}=~/XT:A:U/;
	
	if ($printUniqOnly == 0 || $uniq == 1)
	{
		if ($separateBed && $read1_or_2 == 2)
		{
			print $fout2 bedToLine ($bed), "\n"  unless $flagInfo->{'query_nomap'};
		}
		else
		{
			print $fout bedToLine ($bed), "\n" unless $flagInfo->{'query_nomap'};
		}
	}
}

print STDERR "Done! Totally $i lines processed! \n" if $verbose;

close ($fin) if $inSAMFile ne '-';
close ($fout) if $outBedFile ne '-';
close ($fout2) if $separateBed && $outBedFile2 ne '-';




#subroutines in Align.pm
sub lineToSam
{
	my $line = $_[0];
	my ($QNAME, $FLAG, $RNAME, $POS, $MAPQ, $CIGAR, $MRNM, $MPOS, $ISIZE, $SEQ, $QUAL, $TAGS) = split (/\s+/, $line, 12);
	return {
	QNAME => $QNAME,
	FLAG=> $FLAG,
	RNAME=>$RNAME,
	POS=>$POS,
	MAPQ=>$MAPQ,
	CIGAR=>$CIGAR,
	MRNM=>$MRNM,
	MPOS=>$MPOS,
	ISIZE=>$ISIZE,
	SEQ=>$SEQ,
	QUAL=>$QUAL,
	TAGS=>$TAGS
	};
}


#return 0 if no alignment

sub samToBed
{
	my ($sam, $useRNAStrand) = @_;
	$useRNAStrand = 0 unless $useRNAStrand;

	return 0 if $sam->{"CIGAR"} eq '*'; #no alignment
	
	my $flagInfo = decodeSamFlag ($sam->{"FLAG"});
	
	my $strand = $flagInfo->{'query_strand'};
	
	my $TAGS = "";
	$TAGS = $sam->{"TAGS"} if $sam->{"TAGS"};
	if ($useRNAStrand)
	{
		if ($TAGS=~/XS\:\S*\:([\-\+\.])/)
		{
			$strand = $1;
			$strand = '+' if $strand eq '.';
		}
	}
	my $read1_or_2 = $flagInfo->{'read_1_or_2'};
	
	my $name = $sam->{"QNAME"};
	my $chrom = $sam->{"RNAME"};
	my $chromStart = $sam->{"POS"} - 1;

	my $score = 0;
	if ($TAGS=~/NM\:\S*\:(\d+)/)
	{
		$score = $1;
	}

	my $CIGAR = $sam->{"CIGAR"};
	my $QNAME = $sam->{"QNAME"};
	my $SEQ = $sam->{"SEQ"};

	#remove soft cliped nucleotides
	if ($CIGAR =~/^\d+S(.*?)$/)
	{
		$CIGAR = $1;
	}
	elsif ($CIGAR =~/^(.*?)\d+S$/)
	{
		$CIGAR = $1;
	}

	#deal with the rest
	if ($CIGAR=~/[^\d+|M|N|I|D]/g)
	{
		Carp::croak "unexpected CIGAR string: $CIGAR in $QNAME: $SEQ\n";
	}

	my (@blockSizes, @blockStarts);
	
	my $currLen = 0;
	my $extendBlock = 0;

	while ($CIGAR =~/(\d+)([M|N|I|D])/g)
	{
		my ($size, $type) = ($1, $2);
		if ($type eq 'I' || $type eq 'D')
		{
			#insertion in reads
			$extendBlock = 1;
			if ($type eq 'D')
			{
				my $n = @blockSizes;
				if ($n < 1)
				{
					$chromStart += $size;	
				}
				else
				{
					$blockSizes[$#blockSizes] += $size;
					$currLen += $size;
				}
			}
			next;
		}

		if ($type eq 'M')
		{
			if ($extendBlock && @blockSizes > 0)
			{
				#extend the previous block
				my $n = @blockSizes;
				$blockSizes[$n-1] += $size;
			}
			else
			{
				push @blockSizes, $size;
				push @blockStarts, $currLen;
			}
			$extendBlock = 0;
		}
		$currLen += $size;
	}
	
	my $blockCount = @blockSizes;
	my $chromEnd = $chromStart + $blockStarts[$blockCount-1] + $blockSizes[$blockCount-1] - 1;
	
	my $bed = {
		chrom=>$chrom,
		chromStart=>$chromStart,
		chromEnd=>$chromEnd,
		name=>$name,
		score=>$score,
		strand=>$strand,
		thickStart=>$chromStart,
		thickEnd=>$chromEnd,
		itemRgb=>0,
		blockCount=>$blockCount,
		blockSizes=>\@blockSizes,
		blockStarts=>\@blockStarts,
		flagInfo=>$flagInfo
	};
}	

sub decodeSamFlag
{
	my $flag = $_[0];
	$flag = sprintf ("%012b", $flag);
	my @flags = split (//, $flag);

	my $flagInfo = {
		PE=>$flags[11],					#1 means paired-end data
		PE_map=>$flags[10],				#1 means each end is properly aligned according to the aligner
		query_nomap=>$flags[9],				#1 means this read is unmapped
		mate_nomap=>$flags[8],				#1 means its mate is unmapped
		query_strand=>$flags[7] == 0 ? '+' : '-',	#1 means the strand of this read is on the negative strand
		mate_strand=>$flags[6] == 0 ? '+' : '-',	#1 means its mate is on the negative strand
		read_1_or_2=> $flags[5] == 1 ? 1 : 2 		#1 means this is read1
	};
	return $flagInfo;
}	



#subroutine in Bed.pm
sub bedToLine
{
	my $region = $_[0];
	
	my @colNames = qw (chrom chromStart chromEnd name score strand thickStart thickEnd itemRgb blockCount blockSizes blockStarts);
	
	my $colNum = 12; #keys %$region;
	
	my %rCopy = %$region;
	$rCopy{"chromEnd"} += 1;
	if (exists $rCopy{'thickEnd'})
	{
		$rCopy{'thickEnd'} += 1;
	}
	
	if (exists $rCopy{'blockCount'})
	{
		Carp::croak "no blockSizes\n" unless exists $rCopy {'blockSizes'};
		Carp::croak "no blockStarts\n" unless exists $rCopy {'blockStarts'};
			
		$rCopy{'blockSizes'} = join (",", @{$rCopy{'blockSizes'}});
		$rCopy{'blockStarts'} = join (",", @{$rCopy{'blockStarts'}});
	}
	
	my $ret = join ("\t", $rCopy{"chrom"}, $rCopy{"chromStart"}, $rCopy{"chromEnd"});
	for (my $i = 3; $i < $colNum; $i++)
	{
		my $col = $colNames[$i];
		if (exists $rCopy{$col})
		{
			$ret .= "\t" . $rCopy{$col};
		}
		else
		{
			last;
			#Carp::croak "col=$col is not defined\n"; 
		}
	}
	return $ret;
}



