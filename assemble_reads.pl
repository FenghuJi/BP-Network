use strict;
use POSIX;
use Data::Dumper;
use Bio::SeqIO;

my @aInfiles = ("/home/jfh/data/Virome/Gajewski/03_cdhit/gajewski_phage.fna",
		"/home/jfh/data/Virome/Gajewski/03_cdhit/roux.fna",
		"/home/jfh/data/Virome/Gajewski/03_cdhit/imgvr_gut.fna",
		"/home/jfh/data/Virome/Gajewski/03_cdhit/gvd.fna",
                "/home/jfh/data/Virome/Gajewski/03_cdhit/viral.2.1.genomic.fna");

my $outstream = Bio::SeqIO->new(-file=>">prophage_virus.fa", -format=>'fasta');
open OUT, ">prophage_virus.names.new2old.txt" or die;
my $seqid_new = 'S000000';
foreach my $seqfile ( @aInfiles ){
    
    print STDERR "\t$seqfile ... \n";
    
    my $instream = Bio::SeqIO->new(-file=>$seqfile, -format=>'fasta');
    while( my $seq = $instream->next_seq() ){
        my $seq_id_old = $seq->id();
        my $seq_length = $seq->length();
        
        $seqid_new ++;
        
        ## -- make new seq and write to out file 
        my $new_seq_obj = Bio::Seq->new(-id=>$seqid_new, -seq=>$seq->seq);
        $outstream->write_seq( $new_seq_obj );
        
        ## -- save new name to old name;
        print OUT join("\t", $seqid_new, $seq_id_old, $seq_length), "\n";
    }
}
close OUT;
