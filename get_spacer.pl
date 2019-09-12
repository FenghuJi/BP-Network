#!/usr/bin/perl -w
use strict;

## *********************************************
## ** GET opts **
## *********************************************
use Getopt::Long;
my %opts=();
GetOptions(\%opts,"i:s","o:s");

if (!$opts{i} ){
    print "--------------------------------------------------------------------------------------------------
    USAGE: perl $0
        -i input file
--------------------------------------------------------------------------------------------------\n";
    exit;
}

my $input_file = $opts{i};

my $infile = "01_spacers_result/$input_file";
my $out = "01_spacers/$input_file";

print "$opts{i} \n\ninputfile: $infile \noutput file: $out\n";

open (RAW,"$infile") or  die "can't open the file12";
my @file=<RAW>;
chomp @file;
close RAW;

open OF,">$out";
for (my $i=0; $i <@file; $i++){
	if($file[$i] =~ /^Array\s+(\d+)/){ 
		my $num = $1; 
		if( $file[$i+1] =~ />(\S+)length_(\d+)(\S+)/){
			my $name=$1."length_".$2.$3;
			my $test=0;
			for( my $k =$i+2;$k<@file;$k++){
				if($file[$k] =~ /^\s+(\d+)\s+(\d+)\s+(\S+)\s+(.+)\s+(\S+)\s+(\S+)\s+(\S+)/){
					my $spacer =$7;
					print OF ">".$name."__Array".$num."__".$test."\n".$spacer."\n";
					$test++;
				}elsif($file[$k] =~ /^Array/){
					$i = $k-1; last;
				}
			}
		}			
	}elsif($file[$i] =~ /^SUMMARY BY SIMILARITY/){

		last;
	}
}
close OF;


