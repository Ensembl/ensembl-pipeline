#!/usr/local/ensembl/bin/perl -w

use strict;
use IO::File;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Getopt::Long;

my $fh = new IO::File;

# connect to the database
my $dbhost = 'ecs2f';
my $dbuser = 'ensro';
my $dbname = 'genewisedb_mouse';
my $outdir;


&GetOptions(
	    'dbname:s'       => \$dbname,
	    'dbhost:s'       => \$dbhost,
	    'outdir:s'       => \$outdir,
	    );

my $db= new Bio::EnsEMBL::DBSQL::DBAdaptor(
					   -host  => $dbhost,
					   -user  => $dbuser,
					   -dbname=> $dbname
					   );



my $size = 100000;
my $chr_start = &get_chrstart($db);
my $chr_end   = &get_chrend(  $db);

foreach my $chr ( keys %$chr_start ){
  my $start = $chr_start->{$chr};
  my $end   = $chr_end->{$chr};
  
  my $slice_start = $start;
  my $slice_end   = $start + ( $size - 1 );
  while ( $slice_start < $end ){
    
    my $slice = "$chr.".$slice_start."-".$slice_end;

    my $command = "bsub -q acari -C0 -f  \"/ecs2/work1/eae/GeneWiseHMM/ChrisHMMs.pfam > /tmp/HMMs.pfam\" -o /ecs2/work1/eae/GeneWiseHMM/mouse_jobs_dir/$slice.out -e /ecs2/work1/eae/GeneWiseHMM/mouse_jobs_dir/$slice.err /nfs/acari/eae/ensembl/ensembl-pipeline/scripts/HMMs/run_genewisedb.pl -input_id $slice  -hmm /tmp/HMMs.pfam  ";
    print $command."\n";
    
    #    my $command = "get_sequence_for_chr.pl -chr $chr -outfile $outdir/$chr.fa -mask -dust -softmask -dbname $dbname -dbhost $dbhost";
    #    $fh->open("| bsub -q acari -m ecs2_hosts -o $outdir/jobs/lsf-out-$chr");
    #    $fh->print($command);
    #    $fh->close;
    
    $slice_start += $size;
    $slice_end    = $slice_start + ( $size - 1 );
    if ( $slice_end > $end){
      $slice_end = $end;
    }
  }
}


sub get_chrstart{
  my $db = shift;

  my $q = "SELECT c.name, min(a.chr_start) 
           FROM   chromosome c, assembly a
           WHERE  c.chromosome_id = a.chromosome_id
           GROUP BY c.name";

  my $sth = $db->prepare($q) || $db->throw("can't prepare: $q");
  my $res = $sth->execute || $db->throw("can't execute: $q");
  
  my %chr_start;
  while( my ($chr, $start) = $sth->fetchrow_array) {
    $chr_start{$chr} = $start;
  }
  return \%chr_start;
}

sub get_chrend{
  my $db = shift;

  my $q = "SELECT c.name, max(a.chr_end) 
           FROM   chromosome c, assembly a
           WHERE  c.chromosome_id = a.chromosome_id
           GROUP BY c.name";

  my $sth = $db->prepare($q) || $db->throw("can't prepare: $q");
  my $res = $sth->execute || $db->throw("can't execute: $q");
  
  my %chr_end;
  while( my ($chr, $end) = $sth->fetchrow_array) {
    $chr_end{$chr} = $end;
  }
  return \%chr_end;
}
