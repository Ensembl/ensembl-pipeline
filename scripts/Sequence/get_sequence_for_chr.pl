#!/usr/local/ensembl/bin/perl -w

# dump_seq_into_fastA.pl
# it reads a bit of sequence and dump it into a fasA file, to eb able to view it in Apollo

use strict;
use diagnostics;

use Bio::SeqIO;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Getopt::Long;

# get a contig with a piece-of/entire  chromosome

my $mask;
my $softmask;
my $dust;
my $outfile;
my $chr;
my $maskall;

my $dbhost;
my $dbuser = 'ensro';
my $dbname;

&GetOptions(
	    'chr:s'          => \$chr,
	    'dbname:s'       => \$dbname,
	    'dbhost:s'       => \$dbhost,
	    'mask'           => \$mask,
	    'dust'           => \$dust,
	    'softmask'       => \$softmask,
	    'maskall'        => \$maskall,
	    'outfile:s'      => \$outfile,
	    
	    );

unless( $dbhost && $dbname && $chr && $outfile ){
  print STDERR "Usage: $0 -chr -dbname -dbhost [ -mask -dust -softmask -maskall ] -outfile\n";
  exit(0);
}

# connect to the database
my $db= new Bio::EnsEMBL::DBSQL::DBAdaptor(
					   -host  => $dbhost,
					   -user  => $dbuser,
					   -dbname=> $dbname
					  );


open OUT, ">$outfile";
my $out = Bio::SeqIO->new(-format=>'Fasta',
			  -fh =>  \*OUT,
			 );

my $slice = $db->get_SliceAdaptor->fetch_by_chr_name( $chr );


my @logic_names;
my $soft = 0;
my $string  = '';
if ($mask){
    push( @logic_names, 'RepeatMask' );
}
if ( $dust ){
    push (@logic_names, 'Dust' );
}
if ($softmask){
    $soft = 1;
    $string = 'softmask';
}

my $seq;
if ($maskall){
    $string .= " maskall ";
    $seq = $slice->get_repeatmasked_seq([''],$soft);
}
elsif( @logic_names ){
 print STDERR "getting sequence soft - repeatmaksed and dusted\n"; 
 $seq = $slice->get_repeatmasked_seq(\@logic_names,$soft);
}
else{
    $seq = $slice;
}

$seq->display_id($slice->display_id." @logic_names ".$string);
$out->write_seq($seq);	       
close OUT;

############################################################
