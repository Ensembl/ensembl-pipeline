#!/usr/local/bin/perl -w

# dump_seq_into_fastA.pl
# it reads a bit of sequence and dump it into a fasA file, to eb able to view it in Apollo

use strict;
use diagnostics;

use Bio::SeqIO;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Getopt::Long;



# get a contig with a piece-of/entire  chromosome

my $global_start;
my $global_end; 
my $chr_name;   
my $input_id;
my $masked;
my $softmasked;
my $dust;
my $dbhost;
my $dbname;
my $maskall;


&GetOptions(
	    'input_id:s' => \$input_id,
	    'masked'     => \$masked,
	    'softmasked' => \$softmasked,
	    'maskall'    => \$maskall,
	    'dust'       => \$dust,
	    'dbname:s'   => \$dbname,
	    'dbhost:s'   => \$dbhost,
	   );

unless ( $input_id && $dbname && $dbhost ){
  print STDERR "Usage: $0 -input_id chr_name.chr_start-chr_end -dbname -dbhost ( -masked -dust -maskall -softmasked )\n";
  exit(0);
}


# connect to the database
my $db= new Bio::EnsEMBL::DBSQL::DBAdaptor(-host  => $dbhost,
					   -user  => 'ensro',
					   -dbname=> $dbname);


my $sgp = $db->get_SliceAdaptor;

my $outfile = "$input_id.fa";
$input_id =~/(\S+)\.(\d+)-(\d+)/;
$chr_name     = $1;
$global_start = $2;
$global_end   = $3;

my $slice = $sgp->fetch_by_chr_start_end($chr_name,$global_start,$global_end);

# slice is a Bio::EnsEMBL::PrimarySeq


my @logic_names;
my $soft = 0;
my $string  = '';
if ($masked){
  push( @logic_names, 'RepeatMask' );
}
if ( $dust ){
  push (@logic_names, 'Dust' );
}
if ($softmasked){
  $soft = 1;
  $string = 'softmask';
}

open OUT, ">$outfile";
# get a Bio::SeqIO object
my $out = Bio::SeqIO->new(-format=>'Fasta',
			  -fh =>  \*OUT,
			 );


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




