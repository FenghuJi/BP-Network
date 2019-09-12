#!/usr/bin/perl  
open(INFILE, "results");
open(OUTFILE, ">cluster_results.txt")|| die "Cannot open the newfile: $!\n";
while (<INFILE>) {
        @a = split(" ");
              print OUTFILE "$a[0]\t $a[1]\t $a[2]\t $a[3]\t $a[4]\t $a[5]\t\n";
                 } 
close INFILE;
close OUTFILE;
