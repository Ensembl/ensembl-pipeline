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


my %chrhash;

# declare these here so we can refer to them later
my $pseudogene_bsubdir  = "results/";

# get the gene ids from the database to be 'labelled'
my @gene_ids = &get_gene_ids();

# make output directories
&make_directories();

# create jobs file for Exonerate
&make_pseudogene_bsubs( @gene_ids);

############################################################

sub get_gene_ids{
    my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
						 '-host'   => $LABEL_DBHOST,
						 '-user'   => 'ensro',
						 '-dbname' => $LABEL_DBNAME,
						 );
    return $db->get_GeneAdaptor->list_geneIds;
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
    my @gene_ids = shift;
    
    my $jobfile = $BSUBS_FILE;
    open (OUT, ">$jobfile") or die ("Can't open $jobfile for writing: $!");
    
    my $lsf_options   = $LSF_OPTIONS;
    my $check         = $LABEL_PRE_EXEC;
    my $pseudogene    = $LABEL_SCRIPT;
    my $bsubdir       = $TMPDIR . "/" . $pseudogene_bsubdir . "/";
        
    foreach my $id (@gene_ids){
	my $outfile   = $bsubdir . $id. "_out";
	my $errfile   = $bsubdir . $id. "_err";
	
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
