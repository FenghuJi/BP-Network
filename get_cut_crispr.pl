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

my $infile = "/home/jfh/data/Virome/Gajewski/01_fasta/$input_file.fasta";
my $left = "03_cut/$input_file"."_left";
my $right = "03_cut/$input_file"."_right";
#print "$opts{i} \n\ninputfile: $infile \noutput file: $outfile\n";

open (RAW,"01_spacers_result/$input_file") or  die "can't open the file12";
my @file=<RAW>;
chomp @file;
close RAW;


open (FASTA, "$infile") or  die "can't open the blast!";
local $/ = ">";
		
my (%hash,$heads,$seqs);
while (<FASTA>){
	chomp; chomp;
	($heads,$seqs) = split(/\n/,$_,2);
        next unless($heads && $seqs);
        $seqs=~s/\s+//g;
       	$heads=~s/\s+$//;#
	$hash{$heads}=$seqs;
}#print keys %hash;exit;
close (FASTA);



my ($j,$ij);
for ($j =0; $j < @file; $j++){
	if($file[$j] =~ /^SUMMARY(.+)POSITION/){
		$ij = $j;
		#print $file[$j]; print $j;exit;
	}
}
open OF1 ,">$left";
open OF2 ,">$right";
for (my $i=$ij; $i <@file; $i++){
	if( $file[$i] =~ />(\S+)length_(\d+)(\S+)/){
		my $name=$1."length_".$2.$3;
		#my $seq = $hash{$name};
		my $length = $2;
		#my $tab =$2;
		my $head = $name;
		my $seq = $hash{$head};
		for( my $k =$i+1;$k<@file;$k++){
			if($file[$k] =~ /^\s+(\d+)\s+(\S+)\s+(\d+)\s+(\d+)/){
				my $num = $1; 
				my $start = $3;
				my $sequence = $4;
				my $stop = ($start + $sequence);
				if($start >= 300){
					my $pos = $start-300;
					my $left = substr($seq,($start-300),200);
					print OF1 ">".$name."__Array".$num."__".$length."__".$start."__".$stop."__".$pos."\n".$left."\n";
				}elsif(($start < 300) && ($start >= 150)){
					my $pos = 0;
					my $left = substr($seq,0,($start-50));
					print OF1 ">".$name."__Array".$num."__".$length."__".$start."__".$stop."__".$pos."\n".$left."\n";
				}elsif(($stop+300) < $length){
					my $pos = $stop+100;
					my $right = substr($seq,($stop+100),200);
					print OF2 ">".$name."__Array".$num."__".$length."__".$start."__".$stop."__".$pos."\n".$right."\n";
				}elsif((($stop+300) > $length) && ($length >= ($stop+150)) ){
					my $pos = $stop+50;
					my $right = substr($seq,($stop+50),($length-($stop+50)));
					print OF2 ">".$name."__Array".$num."__".$length."__".$start."__".$stop."__".$pos."\n".$right."\n";
				}
			}elsif($file[$k] =~ />/){
					$i = $k-1; last;
			}
		}		
	}
}
close OF1;
close OF2;



