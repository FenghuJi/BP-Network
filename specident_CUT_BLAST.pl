#!/usr/bin/perl -w
use strict;
use POSIX;
use Data::Dumper;
use Config::Auto;
use String::Util qw(trim);

## *********************************************
## ** version history **
## *********************************************
my $created_on      = "Aug 1, 2017";
my $ver             = '1.0';
my $last_modified   = 'Aug 1, 2017';

## *********************************************
## ** GET opts **
## *********************************************
use Getopt::Long;
my %opts=();
GetOptions(\%opts,"i:s","o:s", "n:s", "c:s", "debug", "ident:f", "hsplen:i", "cov:f");

if (!$opts{i} or !$opts{o}  or !$opts{c} ){
    print "--------------------------------------------------------------------------------------------------
    \t\tby Weihua Chen; created on $created_on
    \t\t\tlast modified : $last_modified; current version: $ver
--------------------------------------------------------------------------------------------------
    USAGE: perl $0
        -c input config file
        -i input a list of accession (file names) to process 
        -o output species name to taxonomy results
      [optional]
        -ident 0~1, default = 0.95; sequence identity
        -hsplen 200 or min-HSP coverage, whichever is shorter, **not used**
        -cov 0~1, default = 0.8; minimal % of the query should be covered by HSP
--------------------------------------------------------------------------------------------------\n";
    exit;
}

###########################################################################
###  -- global variables --
###########################################################################
## -- user-supplied parameters --
my $config_file         = $opts{c};
my $out_file            = $opts{o};
my $in_file             = $opts{i};

use DBI;
my $dbh         = DBI->connect("dbi:mysql:biosql:localhost;mysql_socket=/home/databases/mysql/mysql.sock","wchen","mylonelyboy",{PrintError=>0,RaiseError=>1}) or die "Can't connect to mysql database: $DBI::errstr\n";
my $sth_get_majorranks_by_ncbi_taxon_id = $dbh->prepare( "select species_taxon_id, genus_taxon_id, family_taxon_id, order_taxon_id, class_taxon_id, phylum_taxon_id, superkindom_taxon_id from taxon2majorranks where ncbi_taxon_id = ? limit 1" );
my $sth_get_sciname_by_taxon_id         = $dbh->prepare( "select name from taxon_name where taxon_id = ? and name_class = 'scientific name' limit 1" );
my $sth_get_sciname_by_ncbi_taxon_id    = $dbh->prepare( "select t1.name from taxon_name as t1, taxon as t2 where t1.taxon_id = t2.taxon_id and t2.ncbi_taxon_id = ? and t1.name_class = 'scientific name' limit 1" );
my $sth_get_ncbi_taxon_id_by_taxon_id   = $dbh->prepare(" select ncbi_taxon_id from taxon where taxon_id = ? limit 1");

my @aRanks = qw(species genus family order class phylum superkingdom);

my $debug_mode          = defined $opts{ debug };
if( $debug_mode ){
    print STDERR "\tyou are running at DEBUG mode, the program prints detailed debugging info on STDERR & STDOUT ... \n\n";
}

## -- arbituary cutoffs for the subsequent analyses --
my $CUTOFF_SEQIDENT_STRINGENT           = defined $opts{ident}      ? $opts{ident}      : 0.95; ## min seq identity --
my $CUTOFF_MIN_ALIGN_COV_PCT_STRINGENT  = defined $opts{ cov }      ? $opts{ cov }      : 0.8; ## min ALIGN coverage of the query sequence--

## ======================= all global variables will be declared here ========================
my %hQuery2Parts2Hits        = (); ## $hash{ $query_name }{ $parts } = $num_hits --

print STDERR "\tload and parse the config file ... \n";
my $hrConfigFile = Config::Auto::parse( $config_file );

## -- do some consistency check --
die "DIE: both 'cut_seq_dir' and 'cut_blast_host_dir' are mandatory fields in the configuration file: $config_file ... system quit now ... \n\n"
if( !( exists $$hrConfigFile{ cut_seq_dir } and exists $$hrConfigFile{ cut_blast_host_dir } ) );

###########################################################################
###  -- read and parse BLAST file --
###########################################################################
print STDERR "\tloading BLAST files and process ... \n";
open OUT, ">$out_file" or die;
print OUT join("\t", qw(prophage_seq_name host_taxon_id source contig_file)), "\n";
open LIST, $in_file or die;
while(<LIST>){
    chomp;
    next if( /^#/ or !/\S/ );
    
    my $file_name = trim($_);
    
    ## -- prepare files --
    my %aParts = ( $$hrConfigFile{ cut_seq_dir } . "/" . $file_name . "_left" =>  $$hrConfigFile{ cut_blast_host_dir } . "/" . $file_name . "_left" ,
                  $$hrConfigFile{ cut_seq_dir } . "/" . $file_name . "_right" => $$hrConfigFile{ cut_blast_host_dir } . "/" . $file_name . "_right");
    
    my %hQuery2TaxonIDs = (); ## $hash{ $query }{ $taxons } ++;
    while( my ( $seq_file, $blast_file ) = each %aParts ){
        $seq_file = trim($seq_file);
        $blast_file = trim($blast_file);
        print STDERR "\t\tinfile: ", -e $seq_file ? "Yes" : "No --> '",  $seq_file , "'\n" if( $debug_mode );
        print STDERR "\t\tinfile: ", -e $blast_file ? "Yes": "No --> '", $blast_file, "'\n" if( $debug_mode );
        
        if( -e $seq_file  and -e $blast_file ){
            print STDERR "\t\t\tboth input file exists ... \n";
            ## -- 读取序列文件，以取得序列长度 -- 
            ## -- read from seq file --
            my $hrSeqID2Len = &get_seq2length_from_fasta( $seq_file );
            
            ## -- 读取BLAST结果 --
            ## -- read from blast all file --
            open BLAST, $blast_file or die;
            my $count_blast_results     = 0;
            my $query                   = '';
            my $first_hit               = 1; ##
            my $first_hit_ident_98      = 0;
            while(<BLAST>){
                if( /^# BLASTN/ ){
                    ## reset some parameters --
                    $first_hit              = 1;
                    $first_hit_ident_98     = 0;
                    $count_blast_results    ++; ## 计数
                } elsif (!/^#/) {
                    my ( $query_id, $subject_id, $align_ident, $align_len ) = split(/\t/, $_);
                    my $query_length    = $$hrSeqID2Len{ $query_id };
                    my $align_coverage  = $align_len / $query_length * 100;
                    $align_ident        /= 100;
                    
                    ## -- 直接跳过不合格的 --
                    next if( $align_coverage < $CUTOFF_MIN_ALIGN_COV_PCT_STRINGENT or $align_ident < $CUTOFF_SEQIDENT_STRINGENT );
                    my ( $host_taxon_id ) = $subject_id =~ /^(\d+)/;
                    
                    ## -- 如果是第一条记录 --
                    if( $first_hit ){
                        $first_hit = 0;
                        if( $align_ident >= 0.98 ){
                            $first_hit_ident_98 = 1;
                        }
                    }
                    
                    ## -- 判断第一条记录是否有非常好的比对结果 --
                    ## -- check push into hash  --
                    if( $first_hit_ident_98 ){
                        $hQuery2TaxonIDs{ $query_id }{ $host_taxon_id } ++ if( $align_ident >= 0.98 );
                    } else {
                        $hQuery2TaxonIDs{ $query_id }{ $host_taxon_id } ++;
                    }
                }
            } ## end of BLAST --
            close BLAST;
            print STDERR "\t\t\t loaded ", scalar keys %hQuery2TaxonIDs,  " qualified quries from in total ", $count_blast_results , " BLAST queries loaded ... \n\n";

            ## -- 后期处理 --
            ## -- calculate final taxons --
            if( scalar keys %hQuery2TaxonIDs > 0 ){
                print STDERR "\t\t\tnow parsing the results for NCBI taxonomy ... \n";            
                while( my ( $query_id, $hrNCBITaxonIDs ) = each %hQuery2TaxonIDs ){
                    my @aNCBITaxonIDs = keys %{$hrNCBITaxonIDs};
                    my ( $lca_ncbi_taxon_id, $lca_rank, $lca_sci_name ) = &get_lca_by_a_list_of_ncbi_taxon_ids(@aNCBITaxonIDs);
                    print OUT join("\t", $query_id, $lca_ncbi_taxon_id, "prophage flanking regions", $seq_file, $lca_rank, $lca_sci_name), "\n" if( $lca_ncbi_taxon_id > 0 and ($lca_rank eq 'species' or $lca_rank eq 'genus') );
                }
            }
        }
        ## if both required files exists 
    } ## 循环每个文件 
}
close OUT;


###################################
### sub functions
## -- 输入： array of NCBI taxonomy IDs --
## -- 输出： lca of the input ID --
sub get_lca_by_a_list_of_ncbi_taxon_ids{
    my @aNCBITaxonIDs = @_;
    my $ncbi_taxon_id_count = scalar @aNCBITaxonIDs;
    
    ## -- 
    my ( $lca_ncbi_taxon_id, $lca_rank, $sci_name ) = (0, '', '');
    
    ## -- 
    if( $debug_mode ){
        print STDERR "\t\t\tinput: ", $ncbi_taxon_id_count, "--", join(",", @aNCBITaxonIDs ), "-------- \n";
    }
    
    ## -- 如果 array 只有一个成员 -- 
    if( scalar @aNCBITaxonIDs == 1 ){
        ## -- 1. --
        $lca_ncbi_taxon_id  = $aNCBITaxonIDs[0];
        ## -- 2  &  3 --
        $sth_get_sciname_by_ncbi_taxon_id->execute( $lca_ncbi_taxon_id );
        while( my ($this_sci_name) =  $sth_get_sciname_by_ncbi_taxon_id->fetchrow_array){
            $sci_name = $this_sci_name if( defined $this_sci_name );
            $lca_rank = &get_majorrank_by_ncbi_taxon_id( $lca_ncbi_taxon_id );
        }
    } else {
        ## -- if there are multiple --
        ## -- first get major ranks --
        print STDERR "\t\t\t--- major ranks ---\n" if( $debug_mode );
        my @aMajorRanks = ();
        foreach my $ncbi_taxon_id ( @aNCBITaxonIDs ){
            $sth_get_majorranks_by_ncbi_taxon_id->execute( $ncbi_taxon_id );
            while( my @arr = $sth_get_majorranks_by_ncbi_taxon_id->fetchrow_array() ){
                push @aMajorRanks, \@arr;
                print STDERR "\t\t\t", $ncbi_taxon_id, "==>", join(",", @arr), " ---\n" if( $debug_mode );
            }
        }
        
        ## -- get lca_taxon_id, lca_rank
        my $lca_taxon_id = 0;
        for( my $i = 0; $i < @aRanks; $i ++){
            my %hTaxIDs = ();
            foreach my $arrref ( @aMajorRanks ){
                $hTaxIDs{ $$arrref[ $i ] } ++;
                print STDERR "\t\t\t\t$i --> ", $$arrref[ $i ], "\n" if( $debug_mode ); 
            }
            
            if( scalar keys %hTaxIDs == 1 ){
                ($lca_taxon_id) = keys %hTaxIDs;
                $lca_rank = $aRanks[$i];
                
                $sth_get_sciname_by_taxon_id->execute( $lca_taxon_id );
                while( my ( $this_sci_name ) = $sth_get_sciname_by_taxon_id->fetchrow_array() ){
                    $sci_name = $this_sci_name if( defined $this_sci_name);
                }
                
                $sth_get_ncbi_taxon_id_by_taxon_id->execute( $lca_taxon_id );
                while( my ( $this_ncbi_taxon_id ) = $sth_get_ncbi_taxon_id_by_taxon_id->fetchrow_array() ){
                    $lca_ncbi_taxon_id = $this_ncbi_taxon_id if( defined $this_ncbi_taxon_id );
                }
                
                last;
            }
        }
    }
    
    ## -- debug --
    if( $debug_mode  ) {
        print STDERR "\t\t\tresult: ", join(",", $lca_ncbi_taxon_id, $lca_rank, $sci_name ), "-------- \n\n";
    }
    
    ## -- return --
    return ($lca_ncbi_taxon_id, $lca_rank, $sci_name);
}

sub get_majorrank_by_ncbi_taxon_id{
    my ( $ncbi_taxon_id ) = @_;
    
    my $major_rank = 'NA';
    $sth_get_majorranks_by_ncbi_taxon_id->execute( $ncbi_taxon_id );
    while( my @arr = $sth_get_majorranks_by_ncbi_taxon_id->fetchrow_array() ){
        for( my $i = 0; $i < @arr; $i ++ ){
            if( $arr[$i] > 0 ){
                $major_rank = $aRanks[ $i ];
                last;
            }
        }
    } # -- 
    
    ## --  -- 
    return $major_rank;
}


sub get_seq2length_from_fasta{
    my ($infile) = @_;
    
    my %hash = ();
    use Bio::SeqIO;
    my $instream = Bio::SeqIO->new(-file=>$infile, -format=>'fasta');
    while(my $seq = $instream->next_seq()){
        $hash{ $seq->id() } = $seq->length();
    }
    
    return \%hash;
}
