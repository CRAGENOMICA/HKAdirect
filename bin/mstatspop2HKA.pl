#!/usr/local/bin/perl
use Getopt::Long;
use File::Find;
use strict;

my (%options,@HKAtable,$nloci,$np,$miss);

#options
GetOptions(
     'in=s'      => \$options{'in'},
     'np=s'      => \$options{'np'},
     'f=s'       => \$options{'f'},
);

#Usage
if($options{'in'} eq undef) {
   print "\nUsage: \nperl mstats2HKA.pl -in (inputfile in mstatspop one-line format)  -np (def=0, number of population used) -f (def=1,factor_chromosome) \n\n";
   exit;
}
if($options{'np'} eq undef) {
   $options{'np'} = 0;
}
if($options{'f'} eq undef) {
   $options{'f'} = 1;
}
$np = $options{'np'};

#read the mstatspop and keep the necessary columns
$nloci = 0;
open(INPUT,"$options{'in'}") or die;
while(<INPUT>) {
   if($_=~/Ratio\_Missing:\t(\S+)\t.+nsam\[$np\]:\t(\S+)\t.+Eff\_length2\_pop\[$np\]:\t(\S+)\t.+Eff\_length1\_pop\_outg\[$np\]:\t(\S+)\t.+S\[$np\]:\t(\S+)\t.+Divergence\[$np\]:\t(\S+)\t.+/) {
      push @HKAtable, ($nloci."\t".$2."\t".$5."\t".$6."\t".$3."\t".$4."\t".$options{'f'}."\t".$1);
      $nloci = $nloci + 1;
   }
   elsif($_=~/Ratio\_Missing:\t(\S+)\t.+nsam\[$np\]:\t(\S+)\t.+Eff\_length1\_pop\_outg\[$np\]:\t(\S+)\tEff\_length2\_pop\_outg\[$np\]:\t(\S+)\t.+S\[$np\]:\t(\S+)\t.+Divergence\[$np\]:\t(\S+)\t.+/) {
      if($1 eq "NA") {
        $miss=0;
      } else {
        $miss = $1;
      }
      push @HKAtable, ($nloci."\t".$2."\t".$5."\t".$6."\t".$4."\t".$3."\t".$options{'f'}."\t".$miss);
      $nloci = $nloci + 1;
   }

}
close(INPUT);

#print HKA table to be used as input file in HKAdirect
print "input file for HKAdirect: Original reference from HKA test in Hudson et al. Genetics 1987. See also Ferretti et al. Genetics 2013.\n";
print "nloci ".$nloci."\n";
print "IDloci\tnsam\tS\tD\tLp\tLd\tfchrn\tpmissing\n";
foreach my $locus (@HKAtable){
   print $locus."\n";
}
print "\n";


