#!/usr/local/ensembl/bin/perl -w

=head1 NAME


=head1 SYNOPSIS
 
  make_bsubs.pl
  Makes bsub entries for run_blat.pl, etc...
  bsubs can be submitted using submit.pl - they\'re not automatically 
  done from here as it\'s better to submit a few and check they come 
  back OK before sending a genome worth.

  Makes sure all the various scratch subdirectories needed are in place, 
  and makes them if necessary.

=head1 DESCRIPTION


=head1 OPTIONS

=cut

use strict;
use Getopt::Long;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Pipeline::Config::PseudoGenes::PseudoGenes;

my $gene_type;
&GetOptions(
	    'gene_type:s'       => \$gene_type,
	   );


if ( $gene_type ){
  print STDERR "restricting to genes of type $gene_type\n";
}
else{
  print STDERR "considering all gene types\n";
}

my %chrhash;

# declare these here so we can refer to them later
my $pseudogene_bsubdir  = "results/";

# get the gene ids from the database to be 'labelled'
my @gene_ids = &get_gene_ids();

# make output directories
&make_directories();

# create jobs file for Exonerate
&make_pseudogene_bsubs( \@gene_ids);

############################################################

sub get_gene_ids{
  
  my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
					       '-host'   => $LABEL_DBHOST,
					       '-user'   => 'ensro',
					       '-dbname' => $LABEL_DBNAME,
					      );
  if ( $gene_type ){
    print STDERR "Getting genes by ".$gene_type."\n";
    my @ids = &get_genes_by_type($db,$gene_type);
    return @ids;
  }
  else{
    my $gene_ids = $db->get_GeneAdaptor->list_geneIds;
    return @{$gene_ids};
  }
}

############################################################

sub make_directories {
  my $scratchdir =  $TMPDIR ;
  
  makedir($scratchdir);
  # bsub output directories
  my $bsubdir = $scratchdir . "/" . $pseudogene_bsubdir . "/";
  makedir($bsubdir);
  
}

############################################################

sub make_pseudogene_bsubs {
  my $gene_ids = shift;
  
  my $jobfile = $BSUBS_FILE;
  open (OUT, ">$jobfile") or die ("Can't open $jobfile for writing: $!");
  
  my $lsf_options   = $LSF_OPTIONS;
  my $check         = $LABEL_PRE_EXEC;
  my $pseudogene    = $LABEL_SCRIPT;
  my $bsubdir       = $TMPDIR . "/" . $pseudogene_bsubdir . "/";
  
  foreach my $id (@$gene_ids){
    my $num = int(rand(10));
    my $dir = $bsubdir."/".$num."/";
    if( ! -e $dir ) {
      system( "mkdir $dir" );
    }
    my $outfile   = $dir . $id. "_out";
    my $errfile   = $dir . $id. "_err";
    
    my $command = 
      "bsub $lsf_options -o $outfile -e $errfile -E \"$check \" $pseudogene -gene_id  $id";
    print OUT "$command\n";
  }
  
  close (OUT) or die (" Error closing $jobfile: $!");
}

############################################################

sub makedir{
  my ($dir) = @_;
  if(opendir(DIR, $dir)){ closedir(DIR); }
  else{ system("mkdir $dir") == 0 or die "error creating $dir\n"; }
}

############################################################

sub get_genes_by_type{
  my ($db, $type ) = @_;
  my $query = qq{
    SELECT gene_id
      FROM gene
     WHERE gene.type = "$type"
   };

  my $sth = $db->prepare( $query );
  $sth->execute();
  my @ids;
  while( my $id = $sth->fetchrow_array() ) {
    push( @ids, $id );
  }
  return @ids;
}


