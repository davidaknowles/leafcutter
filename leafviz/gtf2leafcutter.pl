#!/usr/bin/env perl

# MODULES
use strict;
use warnings;
use Getopt::Long;

# GET PARAMETERS
my $sHelp         = 0;
my $sOutputPrefix = "leafviz-annotations";
GetOptions("help!"   => \$sHelp, "output:s" => \$sOutputPrefix); 
# GLOBALS
my %hMultiCopyTags = ('tag'=>'','Dbxref'=>'','gbkey'=>'','Name'=>'','product'=>'','Note'=>'','tss_id'=>'', 'ont'=>'');


# PRINT HELP
$sHelp = 1 unless(@ARGV>0);
if ($sHelp) {
   my $sScriptName = ($0 =~ /^.*\/(.+$)/) ? $1 : $0;
   die <<HELP

   Usage: $sScriptName <gtf-file>
   
   Convert a gtf file to a set of annotation files to be used with leafviz.
   
   Arguments: 
    -o --output <string>
      Output file prefix. Default: $sOutputPrefix
    
    -help
      This help message
   
HELP
}


##########
## MAIN ##
##########

# Collect attribute combination counts
my %hOut;
my $nCountSkips = 0;
my $nCountParse = 0;

my $infile=$ARGV[0];

if ($infile =~ /.gz$/) {
    open(IN, "gunzip -c $infile |") || die "can't open pipe to $infile";
}
else {
    open(IN, $infile) || die "can’t open $infile";
}

while (<IN>){
   next if (/^\s*$/);
   next if (/^ *#/);
   s/[\n\r]+$//;
   my ($sChr, $sSource, $sFeature, $nStart, $nEnd, $nScore, $sStrand, $nFrame, $sGroup) = split /\t/, $_, -1;
   if ( lc($sFeature) eq 'exon'){
      my $rhAnnots = gtf_annots_to_hash($sGroup, \%hMultiCopyTags);
      
      # Get required annotation data
      my $sGeneName     = exists $rhAnnots->{'gene_name'}     ? $rhAnnots->{'gene_name'} : "";
      my $sGeneID       = exists $rhAnnots->{'gene_id'}       ? $rhAnnots->{'gene_id'} : "";
      my $sTranscriptID = exists $rhAnnots->{'transcript_id'} ? $rhAnnots->{'transcript_id'} : "";
      my $sTag          = exists $rhAnnots->{'tag'}           ? join("|", @{$rhAnnots->{'tag'}}) : "";
      
      # Fall back options if gene name is empty
      unless ($sGeneName){
         $sGeneName = $sGeneID ? $sGeneID : 'Unknown';
      }
      
      # Try to get the biotype; there are differences between gtf versions so let's traverse the various options in order of preference
      my $sBiotype = "";
      if ( exists $rhAnnots->{'transcript_type'} ){
         $sBiotype = $rhAnnots->{'transcript_type'};
      }
      elsif ( exists $rhAnnots->{'gene_type'} ){
         $sBiotype = $rhAnnots->{'gene_type'};
      }
      elsif ( exists $rhAnnots->{'gene_biotype'} ){
         $sBiotype = $rhAnnots->{'gene_biotype'};
      }
      else{
         $sBiotype = "Unknown";
      }
      
      # Write to hash with transcript_id as the key
      if ($sTranscriptID){
         if ( exists $hOut{$sTranscriptID} ){
            # Check for basic consistency between exon annotations
            if ( ($sChr eq $hOut{$sTranscriptID}{'chr'}) and ($sStrand eq $hOut{$sTranscriptID}{'strand'}) and ($sGeneName eq $hOut{$sTranscriptID}{'gene_name'}) and ($sGeneID eq $hOut{$sTranscriptID}{'gene_id'}) ){
               push @{ $hOut{$sTranscriptID}{'exons'} }, [sort {$a <=> $b} ($nStart, $nEnd) ]; # Push exon block as an array of arrays.
            }
            else{
               $nCountSkips++;
            }
         }
         else{
            $hOut{$sTranscriptID}{'chr'}           = $sChr;
            $hOut{$sTranscriptID}{'strand'}        = $sStrand;
            $hOut{$sTranscriptID}{'gene_name'}     = $sGeneName;
            $hOut{$sTranscriptID}{'gene_id'}       = $sGeneID;
            $hOut{$sTranscriptID}{'transcript_id'} = $sTranscriptID;
            $hOut{$sTranscriptID}{'biotype'}       = $sBiotype;
            $hOut{$sTranscriptID}{'tag'}           = $sTag;
            push @{ $hOut{$sTranscriptID}{'exons'} }, [sort {$a <=> $b} ($nStart, $nEnd) ]; # Push exon block as an array of arrays.
            $nCountParse++;
         }
      }
      else{
         $nCountSkips++;
      }
      
   }
}
close IN;

# Print some basic stats
warn ("Read $nCountParse transcripts\n");
warn ("Skipped $nCountSkips incorrectly formatted exon entries\n") if ($nCountSkips);


# At this point we gathered everything we need from the GTF file; time to write the outputs
if (keys %hOut){
   
   # Open output file handles
   open ALLEXONS,   "| gzip -c > ${sOutputPrefix}_all_exons.txt.gz"   or die "Error can't write to '${sOutputPrefix}_all_exons.txt.gz': $!\n";
   open ALLINTRONS, "| sort -s -t '\t' -k5,5 -k7,7 -k8,8n | gzip -c > ${sOutputPrefix}_all_introns.bed.gz" or die "Error can't write to '${sOutputPrefix}_all_introns.bed.gz': $!\n";
   open FIVEPRIME,  "| sort -s -t '\t' -k5,5 -k7,7 -k8,8n | gzip -c > ${sOutputPrefix}_fiveprime.bed.gz"   or die "Error can't write to '${sOutputPrefix}_all_fiveprime.bed.gz': $!\n";
   open THREEPRIME, "| sort -s -t '\t' -k5,5 -k7,7 -k8,8n | gzip -c > ${sOutputPrefix}_threeprime.bed.gz"  or die "Error can't write to '${sOutputPrefix}_all_threeprime.bed.gz': $!\n";
   
   # Write a header for the exon file
   print ALLEXONS "chr\tstart\tend\tstrand\tgene_name\n";
   
   # Start cycling through transcripts
   foreach my $sTranscriptID (keys %hOut){
      my %hT = %{$hOut{$sTranscriptID}}; # Transcript details
      my @aE = @{$hT{'exons'}};          # Array of arrays with exons
      
      # Sort exons by start coordinate, ascending
      @aE = sort { $a->[0] <=> $b->[0] } @aE;
      
      # Write to all_exon file
      foreach my $rExon (@aE){
         my ($nStart, $nEnd) = @$rExon;
         print ALLEXONS join ("\t", $hT{'chr'}, $nStart, $nEnd, $hT{'strand'}, $hT{'gene_name'}), "\n";
      }
      
      # Write to all_intron file and write 5' and 3' ends
      for ( my $i = 0 ; $i < $#aE ; $i++){
         my $nIstart = $aE[$i][1];
         my $nIend   = $aE[$i+1][0];
         if ( ($nIend - $nIstart) > 0 ){
            my $nIntronID = $hT{'strand'} eq '+' ? $i+1 : $#aE-$i;
            print ALLINTRONS join ("\t", $hT{'chr'}, $nIstart, $nIend, $hT{'gene_name'}, $hT{'gene_id'}, $hT{'strand'}, $hT{'transcript_id'}, $nIntronID, $hT{'biotype'}, $hT{'tag'} ), "\n";
            print FIVEPRIME  join ("\t", $hT{'chr'}, $nIstart, $nIstart + 1, $hT{'gene_name'}, $hT{'gene_id'}, $hT{'strand'}, $hT{'transcript_id'}, $nIntronID, $hT{'biotype'}, $hT{'tag'} ), "\n";
            print THREEPRIME join ("\t", $hT{'chr'}, $nIend,   $nIend + 1,   $hT{'gene_name'}, $hT{'gene_id'}, $hT{'strand'}, $hT{'transcript_id'}, $nIntronID, $hT{'biotype'}, $hT{'tag'} ), "\n";
            
            # Set 5' and 3' ends based on whether the feature is on the forward on reverse strand; seems to break leafviz visualization so disabled for now
            #if ( $hT{'strand'} eq '+' ){
               #print FIVEPRIME  join ("\t", $hT{'chr'}, $nIstart, $nIstart + 1, $hT{'gene_name'}, $hT{'gene_id'}, $hT{'strand'}, $hT{'transcript_id'}, $nIntronID, $hT{'biotype'}, $hT{'tag'} ), "\n";
               #print THREEPRIME join ("\t", $hT{'chr'}, $nIend,   $nIend + 1,   $hT{'gene_name'}, $hT{'gene_id'}, $hT{'strand'}, $hT{'transcript_id'}, $nIntronID, $hT{'biotype'}, $hT{'tag'} ), "\n";
            #}
            #else{
               #print FIVEPRIME  join ("\t", $hT{'chr'}, $nIend,   $nIend + 1,   $hT{'gene_name'}, $hT{'gene_id'}, $hT{'strand'}, $hT{'transcript_id'}, $nIntronID, $hT{'biotype'}, $hT{'tag'} ), "\n";
               #print THREEPRIME join ("\t", $hT{'chr'}, $nIstart, $nIstart + 1, $hT{'gene_name'}, $hT{'gene_id'}, $hT{'strand'}, $hT{'transcript_id'}, $nIntronID, $hT{'biotype'}, $hT{'tag'} ), "\n";
            #}
         }
      }
   }
   
   # Close filehandles
   close ALLEXONS;
   close ALLINTRONS;
   close FIVEPRIME;
   close THREEPRIME;
}
else{
   warn("No exon entries found in input file\n");
}



#################
## SUBROUTINES ##
#################

# gtf_annots_to_hash
#
# Parse gtf annotations and return key-value pairs
sub gtf_annots_to_hash {
   my ($sAnnots, $rMultiCopy) = @_;
   my %hReturn;
   
   my @asPairs = split / *; */, $sAnnots;
   foreach my $sPair (@asPairs){
      my ($sKey, $sVal) = split(/ /, $sPair, 2);
      $sVal =~ s/"//g;
      if ( exists $rMultiCopy->{$sKey} ){
         push @{$hReturn{$sKey}}, $sVal;
      }
      else{
         die "Error: Duplicate key entry '$sKey' found with values '$hReturn{$sKey}' and '$sVal'\n" if (exists $hReturn{$sKey});
         $hReturn{$sKey} = $sVal;
      }
   }
   return \%hReturn;
}

