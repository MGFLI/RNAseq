#!/usr/bin/perl -w

# count reads per sequence/transcript, optional normalization (RPM/RPKM)

# Author: Marco Groth
# date: 05.2010
# last update:	07.2012



use Data::Dumper 'Dumper';
use warnings;
use File::Basename;
use Getopt::Std;
use POSIX;

my $line;
my @linea;
my $mReads;

my($script) = basename($0);
my $usage = <<EOF;

USAGE: $script -e <export file> -m <modus of counting> [-a <transcript info file>] [-h]

DESCRIPTION
	 count reads per sequence/transcript, optional normalization (RPM/RPKM)

INPUT:	
	-e <export file>		mapping data in Illumina's export format
	
	-m <modus of counting>		how to count
					  -m rpt	count reads per transcript 
					  		(only transcripts with at least one hit will arise)
							= reads per transcript
							
					  -m rpm	equal to -m rpt, additional normalize for library size
					  		= reads per transcript per 1.000.000 mappable reads
							
					  -m rpkm	equal to -m rpm, additional normalize for transcript length
					  		transcripts without hits will arise also 
					  		additional information needed by -a
							= reads per 1kb transcript per 1.000.000 mappable reads
	
OPTIONS:
	-a <transcript info file>	additional information about transcript lengths if -m rpkm is choosen
					format: transcriptID as used in reference \t transcript length 
					e.g. DDB0349044|DDB_G0349042	1050
					
	-p 				set pseudo-count for each transcript (only possible with -m rpkm)

	-h				this help
	

	
EOF



if ($#ARGV < 3){die $usage;}


my %Options;
getopts('e:m:a:ph', \%Options);

if (not exists $Options{'e'}){die print("Please provide export file...\n");}

if ($Options{'h'}){die $usage;} 

# export file
my($fasta) = $Options{'e'};


# check whether transcript info file is needed
if (($Options{'m'} eq 'rpkm') and (not exists $Options{'a'}))
{
 die print("Please provide additional transcript info file using option -a or see usage (-h)...\n\n");
}

if (($Options{'m'} eq 'rpkm') and (exists $Options{'a'}))
{
 if (not exists $Options{'p'})
 {
  print("\nPseudo-count OFF\n\n");
 }
 else
 {
  print("\nPseudo-count ON\n\n");
 }
 
 my($dataset) = $Options{'a'};
 
 unless (open(DATASET,$dataset))
 {
  die print("Unable to open $dataset\n");
 }
}

unless (open(FASTA,$fasta))
{
 print("Unable to open $fasta\n");
 exit;
}


my %RPKM;

if ($Options{'m'} eq 'rpkm')
{
 print STDERR ("\n...get ID & length from dataset file...\n\n\n");


 # Trancripts getting a pseudocount if -p is set

 while($line = <DATASET>)
 {
  chomp($line);
  @linea = split(/\t/,$line);
  $RPKM{$linea[0]}{length} = $linea[1];

  if (exists $Options{'p'}){$RPKM{$linea[0]}{counts} = '1';}
  else {$RPKM{$linea[0]}{counts} = '0';}
 }
}


my %rest;
# count values for each reference in export file
$rest{mHits} = 0;
$rest{QC} = 0;
$rest{NM} = 0;

print STDERR ("\n...count hits from solexa export file...\n");
while($line = <FASTA>)
{
 chomp($line);
 @linea = split(/\t/,$line);
 
 if ($linea[10] eq 'NM'){$rest{NM}++;}
 elsif ($linea[10] eq 'QC'){$rest{QC}++;}
 elsif ($linea[10] =~ m/.*:.*/){$rest{mHits}++;}
 else
 {
  $RPKM{$linea[10]}{counts}++;
  $mReads++;
 }
}


# output of data

if ($Options{'m'} eq 'rpt')
{
 print("mappable reads\t$mReads\nnot mappable read\t$rest{NM}\nmultiple hits\t$rest{mHits}\nbad quality (QC)\t$rest{QC}\n\n\n");
 print("ID\tcounts\n");

 my $ID;

 foreach $ID (keys %RPKM)
 {
  print("$ID\t$RPKM{$ID}{counts}\n");
 }
}

if ($Options{'m'} eq 'rpm')
{
 print("mappable reads\t$mReads\nnot mappable read\t$rest{NM}\nmultiple hits\t$rest{mHits}\nbad quality (QC)\t$rest{QC}\n\n\n");
 print("ID\tcounts\tRPM\n");

 my $ID;

 foreach $ID (keys %RPKM)
 {
  print("$ID\t$RPKM{$ID}{counts}\t",$RPKM{$ID}{counts} / $mReads * 1000000,"\n");
 }
}


if ($Options{'m'} eq 'rpkm')
{
 print("mappable reads\t$mReads\nnot mappable read\t$rest{NM}\nmultiple hits\t$rest{mHits}\nbad quality (QC)\t$rest{QC}\n\n\n");
 print("ID\tcounts\tRPM\tlength of reference\tRPKM\n");

 my $ID;

 foreach $ID (keys %RPKM)
 {
  $RPKM{$ID}{RPK} = $RPKM{$ID}{counts} / $RPKM{$ID}{length} * 1000; 
  $RPKM{$ID}{RPKM} = $RPKM{$ID}{RPK} / $mReads * 1000000;
  print("$ID\t$RPKM{$ID}{counts}\t",$RPKM{$ID}{counts} / $mReads * 1000000,"\t$RPKM{$ID}{length}\t$RPKM{$ID}{RPKM}\n");
 }
}



close FASTA;
close DATASET;
exit;
