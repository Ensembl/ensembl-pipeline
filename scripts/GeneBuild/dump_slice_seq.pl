#!/usr/local/bin/perl -w

=head1 NAME

  dump_vc_seq.pl

=head1 SYNOPSIS
 
  chr_genedump.pl

=head1 DESCRIPTION

  Dumps the sequence of a virtual contig to file, masked or unmasked.
  db and dnadb can be the same or different databases - just make sure all the relevant arguments are passed in.
  The user for each db can (and probably should) be read-only.
  Specify chr_name, start and end.

=head1 OPTIONS

  -dbname
  -dbhost
  -dnadbname
  -dnadbhost
  -path
  -start
  -end
  -chr_name
  -out
  -path
  -masked

=cut

use strict;
use Getopt::Long;
use Bio::SeqIO;
use Bio::EnsEMBL::DBSQL::DBAdaptor;

my $dbname     = '';
my $dbuser     = 'ensro'; # always use read only
my $dbhost     = '';
my $dnadbname  = '';
my $dnadbuser  = 'ensro'; # always use ensro
my $dnadbhost  = '';
my $path       = '';
my $start;
my $end;
my $outfile    = 'out.fa';
my $chr_name    = '';
my $masked     = 0;

&GetOptions( 
	    'dbname:s'     => \$dbname,
	    'dbhost:s'     => \$dbhost,
	    'dnadbname:s'  => \$dnadbname,
	    'dnadbhost:s'  => \$dnadbhost,
	    'outfile:s'    => \$outfile,
	    'chrname:s'    => \$chr_name,
	    'path:s'       => \$path,
	    'start:i'      => \$start,
	    'end:i'        => \$end,
	    'masked'       => \$masked,
	   );

# usage
if(!defined $dbname    ||
   !defined $dnadbname ||
   !defined $dbhost    ||
   !defined $dnadbhost ||
   !defined $chr_name   ||
   !defined $start     ||
   !defined $end       ||
   !defined $path    
  ){
  print  "USAGE: dump_slice_seq.pl -dbname dbname -host host -chrname chr -start start_pos -end end_pos -path path\n optional:\n -masked, for repmasked sequence;\n -outfile to specify a filename other than out.fa\n";
  exit(1);
}

# global stuff
my $dnadb =  new Bio::EnsEMBL::DBSQL::DBAdaptor(
						-host   => $dnadbhost,
						-dbname => $dnadbname,
						-user   => $dnadbuser,
					       );

my $db =  new Bio::EnsEMBL::DBSQL::DBAdaptor(
					     -host   => $dbhost,
					     -dbname => $dbname,
					     -user   => $dbuser,
					     -dnadb  =>$dnadb
					    );

$db->assembly_type($path);

my $sgpa = $db->get_SliceAdaptor;

print STDERR "about to fetch $chr_name $start $end\n";

print "fetching virtual contig for ".$chr_name." ".$start." ".$end."\n";
my $vc = $sgpa->fetch_by_chr_start_end($chr_name,$start,$end);


my $seqout = Bio::SeqIO->new( '-format' => 'fasta',
                              '-file'   => ">>$outfile");
my $did = $chr_name . "." . $start . "-" . $end;
if($masked){
  print STDERR "about to get repeatmasked sequence\n";
  my $rmseq = $vc->get_repeatmasked_seq;
  $rmseq->display_id($did);
  print STDERR "about to write sequence to $outfile\n";

  $seqout->write_seq($rmseq);	       
}

else {
#  my $did = $chr_name . "." . $start . "-" . $end;
  $vc->display_id($did);
  print STDERR "about to write sequence to $outfile\n";
  $seqout->write_seq($vc);
}
