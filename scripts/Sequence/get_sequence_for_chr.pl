#!/usr/local/ensembl/bin/perl -w

# dump_seq_into_fastA.pl
# it reads a bit of sequence and dump it into a fasA file, to eb able to view it in Apollo

use strict;

use Bio::SeqIO;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Getopt::Long;

# get a contig with a piece-of/entire  chromosome

my $mask;
my $softmask;
my $dust;
my $outfile;
my $name;
my $maskall;
my $coord_system;
my $dbhost;
my $dbuser = 'ensro';
my $dbname;
my $dbpass;
my $dbport = 3306;

&GetOptions(
            'seq_region_name:s' => \$name,
            'dbname:s'       => \$dbname,
            'dbhost:s'       => \$dbhost,
            'dbpass:s' => \$dbpass,
            'dbuser:s' => \$dbuser,
            'dbport:s' => \$dbport,
            'mask'           => \$mask,
            'dust'           => \$dust,
            'softmask'       => \$softmask,
            'maskall'        => \$maskall,
            'outfile:s'      => \$outfile,
            'coord_system:s' => \$coord_system,
	    
	    );

unless( $dbhost && $dbname && $name && $outfile ){
  print STDERR ("Usage: -dbname $dbname -dbhost $dbhost -dbuser $dbuser ".
                "-seq_region_name $name -coord_system $coord_system ".
                "[ -mask -dust -softmask -maskall ] -outfile $outfile\n");
  exit(0);
}

# connect to the database
my $db= new Bio::EnsEMBL::DBSQL::DBAdaptor(
                                           -host  => $dbhost,
                                           -user  => $dbuser,
                                           -port => $dbport,
                                           -dbname=> $dbname,
                                           -dbpass => $dbpass,
                                          );


open OUT, ">$outfile";
my $out = Bio::SeqIO->new(-format=>'Fasta',
			  -fh =>  \*OUT,
			 );

my $slice = $db->get_SliceAdaptor->fetch_by_region($coord_system, $name );


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
