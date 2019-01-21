#! /usr/bin/perl

use strict;
use warnings;

##parser for HLAminer to write to standardised format
##parses file name for sampleID (splits on first "." keeping [0])
##format outputs:
#HLA-XYZ XYZ*01:001 for all precitions
##strategy: if previous line started ^HLA, \
#and this line contains "Prediction", \
#then until next line !~ "\*", \
#split on "," and return [0]=~s/\s+//;

##file naming
my $fileinp = $ARGV[0];
my @outname = split(/\./, $fileinp);
my $fileout = $outname[0] . "." . $outname[1] . ".tsv";

##iterate over file
open(IN, $fileinp);

##to write output
open(OUT, ">$fileout");
print OUT "SAMPLEID\tHLAminer_Call\n";

##flags to determine course of action
my $HLAline=0;
my $PRDline=0;
my $KEPline=0;

while(<IN>){
  if($_=~m/^HLA/){
    $HLAline=1;
    next;
  }
  if($_=~m/Prediction/){
    $PRDline=1;
    next;
  }
  if(($PRDline==1) && ($_ ne "")){
    if($_=~m/:/){
      my @sp=split(/\,/,$_);
      my $HLAcall=$sp[0];
      $HLAcall=~s/\s+//g;
      print OUT $outname[0] . "\t" . $HLAcall . "\n";
      next;
    }
  }
}

close IN;
close OUT;
