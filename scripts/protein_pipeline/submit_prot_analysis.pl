#!/usr/local/bin/perl
# Author: Emmanuel Mongin
# Creation: 03.19.2001


=head1 Run_protein_RunnableDB

=head2 Description

This script will submit run_protein_RunnableDB

=cut

use strict;

BEGIN {
    require "Bio/EnsEMBL/Pipeline/Prot_analysis_Conf.pl";
}

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBLoader;

my %db_conf =  %::db_conf;
my %scripts_conf = %::scripts_conf;

my %runnables;

#Get the DB variables
my $dbname = $::db_conf{'dbname'};
my $dbhost = $::db_conf{'dbhost'};
my $dbuser = $::db_conf{'dbuser'};
my $dbpass = $::db_conf{'dbpass'};
my $dnadbhost = $::db_conf{'dnadbhost'};
my $dnadbname = $::db_conf{'dnadbname'};

#Get the DB handler

my $dnadb = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
        -user   => 'ensro',
        -dbname => $dnadbname,
        -host   => $dnadbhost,
        );

my $db =  new Bio::EnsEMBL::DBSQL::DBAdaptor(
					     -host             => $dbhost,
					     -user             => $dbuser,
					     -dbname           => $dbname,
					     -pass             => $dbpass,
					     -dnadb            => $dnadb,
					    );

#Get the location of the peptide file
my $pep_file = $::scripts_conf{'pep_file'};

#Get the location of the scratch directory
my $scratchdir =  $::scripts_conf{'tmpdir'};

print STDERR "SCRATCH: ".$scratchdir."\n";

#Give for each analysis its corresponding runnable

#prints, prosite, profile, pfam, scanprosite, tmhmm, coils, signalp, seg

$runnables{'Prints'} = "Bio::EnsEMBL::Pipeline::RunnableDB::Protein::Prints";
$runnables{'Prosite'} = "Bio::EnsEMBL::Pipeline::RunnableDB::Protein::ScanProsite";
$runnables{'Pfam'} = "Bio::EnsEMBL::Pipeline::RunnableDB::Protein::Hmmpfam   ";
$runnables{'Tmhmm'} = "Bio::EnsEMBL::Pipeline::RunnableDB::Protein::Tmhmm";
$runnables{'ncoils'} = "Bio::EnsEMBL::Pipeline::RunnableDB::Protein::Coil";
$runnables{'Signalp'} = "Bio::EnsEMBL::Pipeline::RunnableDB::Protein::Signalp";
$runnables{'Seg'} = "Bio::EnsEMBL::Pipeline::RunnableDB::Protein::Seg";
$runnables{'Profile'} = "Bio::EnsEMBL::Pipeline::RunnableDB::Protein::Profile";
$runnables{'Superfamily'} = "Bio::EnsEMBL::Pipeline::RunnableDB::Protein::Superfamily";

my @ids = &get_ids();
#my @ids = (28816);

&make_directories();
&chunk_pepfile();

&run_jobs(@ids);

sub get_ids {
    my @id;
    
#Get all of the protein ids for the given database
    
    my $q = "select translation_id from transcript";
    
    my $sth = $db->prepare($q) || $db->throw("can't prepare: $q");
    my $res = $sth->execute || $db->throw("can't execute: $q");
    
    while( my ($prot_id) = $sth->fetchrow_array) {
	push(@id,$prot_id);
    }
    return @id;
}


sub make_directories {
#Find out which analysis should be run
    my $run = $::scripts_conf{'2berun'};
    
    my @toberun = split /,/, $run;

    my $chunks = $scratchdir."/chunks";
    &makedir($chunks);

#First make the tmp directories for each analysis which is supposed to run
    foreach my $r(@toberun) {
	&makedir($scratchdir."/".$r);

	my $dir_err = $scratchdir."/$r/stderr";
	my $dir_out = $scratchdir."/$r/stdout";
	&makedir($dir_err);
	&makedir($dir_out);
    }
}

sub chunk_pepfile {
    open (PEPFILE, "$pep_file");
    my $count = 0;
    my $chunk = 1;

    my $size = $::scripts_conf{'chunk_size'};

    $/ = "\>";
    
    while(<PEPFILE>){
	
	if ($_ ne "\>") {
	    if ($count == 0) {
		open (CHUNK,">".$scratchdir."/chunks/chunk.$chunk");
	    }
	    
	    $_ =~ s/\>$//;  
	    
	    print CHUNK ">$_";
	    $count++;
	    if ($count == $size) {
		$count = 0;
		$chunk++;
	    }
	}
    }
    $/ = "\n";
}

sub run_jobs {
    my (@id) = @_;
    
    #Find out which analysis should be run
    my $run = $::scripts_conf{'2berun'};

    my @toberun = split /,/, $run;

    my $queue = $::scripts_conf{'queue'};
    my $runner = $::scripts_conf{'runner'};
    
    my $pep_file = $::scripts_conf{'pep_file'};
    
     foreach my $r(@toberun) {
	 #First get the analysisId of the module which is supposed to run

	 print STDERR "R: $r\n";
	 
	 my $q = "select analysisId from analysisprocess where module = '$r'";
    
	 my $sth = $db->prepare($q) || $db->throw("can't prepare: $q");
	 my $res = $sth->execute || $db->throw("can't execute: $q");
    
	 my ($analysis_id) = $sth->fetchrow_array;

	 unless ($analysis_id) {die "AnalysisId not defined\n"};

	 my $chunk_name = $r."_chunk";

	 my $chunk = $::scripts_conf{$chunk_name};

	 unless ($chunk) {die "No chunk option defined\n"};

	 my $runnabledb = $runnables{$r};

	 unless ($runnabledb) {die "No runnableDB defined\n"};

	 my $check = "-E \"".$runner." -check -runnable ".$runnabledb."\"";

	 if ($chunk == 1) {
	     foreach my $i(@id) {
		 my $command = "bsub -q ".$queue." -o ".$scratchdir."/".$r."/stdout/".$i.".out -e ".$scratchdir."/".$r."/stderr/".$i.".err ".$check." ".$runner." -runnable ".$runnabledb." -db ".$db." -input_id ".$i." -analysis ".$analysis_id;
				 
		 print STDERR "RUNNING: $command\n";
		 
		 system($command)==0 || die "$0\Error running '$command' : $!";

	     }
	 }

	 
	 if ($chunk == 3) {
	     my $command = "bsub -q ".$queue." -o ".$scratchdir."/".$r."/stdout/".$r.".out -e ".$scratchdir."/".$r."/stderr/".$r.".err ".$check." ".$runner." -runnable ".$runnabledb." -db ".$db." -input_id ".$pep_file." -analysis ".$analysis_id;
	    

	     print STDERR "RUNNING: $command\n";
	     
	     system($command)==0 || die "$0\Error running '$command' : $!";
	     
	    
	 }
	 
	 if ($chunk == 2) {
	     my $dir = $scratchdir."/chunks";     
	     opendir(DIR, $dir);
	     my @allfiles = readdir DIR;
	     closedir DIR;
	     
	     foreach my $f(@allfiles) {
		 if (($f ne ".") && ($f ne "..")) {
		      my $command = "bsub -q ".$queue." -o ".$scratchdir."/".$r."/stdout/".$f.".out -e ".$scratchdir."/".$r."/stderr/".$f.".err ".$check." ".$runner." -runnable ".$runnabledb." -db ".$db." -input_id ".$dir."/".$f." -analysis ".$analysis_id;
		  

		      print STDERR "RUNNING: $command\n";
	     
		      system($command)==0 || die "$0\Error running '$command' : $!";

		  }
	     }
	 }
		     
     }
}


sub makedir{
  my ($dir) = @_;
  if(opendir(DIR, $dir)){ closedir(DIR); }
  else{ system("mkdir $dir") == 0 or die "error creating $dir\n"; }
}
