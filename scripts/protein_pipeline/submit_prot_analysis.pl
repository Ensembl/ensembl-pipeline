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

#Get the location of the peptide file
my $pep_file = $::scripts_conf{'pep_file'};

#Get the location of the scratch directory
my $scratchdir =  $::scripts_conf{'tmpdir'};

#Give for each analysis its corresponding runnable

#prints, prosite, profile, pfam, scanprosite, tmhmm, coils, signalp, seg

$runnables{'prints'} = "Bio::EnsEMBL::Pipeline::RunnableDB::Protein::Prints";
$runnables{'prosite'} = "Bio::EnsEMBL::Pipeline::RunnableDB::Protein::ScanProsite";
$runnables{'scanprosite'} = "Bio::EnsEMBL::Pipeline::RunnableDB::Protein::ScanProsite";
$runnables{'pfam'} = "Bio::EnsEMBL::Pipeline::RunnableDB::Protein::ParacelHMM";
$runnables{'tmhmm'} = "Bio::EnsEMBL::Pipeline::RunnableDB::Protein::Tmhmm";
$runnables{'coils'} = "Bio::EnsEMBL::Pipeline::RunnableDB::Protein::Coil";
$runnables{'signalp'} = "Bio::EnsEMBL::Pipeline::RunnableDB::Protein::Signalp";
$runnables{'seg'} = "Bio::EnsEMBL::Pipeline::RunnableDB::Protein::Seg";

#my @ids = &get_ids();
&make_directories();
&chunk_pepfile();
&run_jobs();

sub get_ids {
    my @id;
    
#Get all of the protein ids for the given database
    my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host   => $::db_conf{'dbhost'},
						-user   => $::db_conf{'dbuser'},
						-dbname => $::db_conf{'dbname'},
					    );
    
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
	 my $chunk_name = $r."_chunk";

	 my $chunk = $::scripts_conf{$chunk_name};

	 my $runnabledb = $runnables{$r};

	 my $check = "-E \"".$runner." -check -runnable ".$runnabledb."\"";

	 if ($chunk == 1) {
	     foreach my $i(@id) {
		 my $command = "bsub -q ".$queue." -o ".$scratchdir."/".$r."/stdout/".$i.".out -o ".$scratchdir."/".$r."/stderr/".$i.".err ".$check." ".$runner." -runnable ".$runnabledb." -dbuser ".$dbuser." -pass ".$dbpass." -dbname ".$dbname." -host ".$dbhost." -input_id ".$i;
				 
	     }
	 }

	 
	 if ($chunk == 0) {
	     my $command = "bsub -q ".$queue." -o ".$scratchdir."/".$r."/stdout/".$r.".out -o ".$scratchdir."/".$r."/stderr/".$r.".err ".$check." ".$runner." -runnable ".$runnabledb." -dbuser ".$dbuser." -pass ".$dbpass." -dbname ".$dbname." -host ".$dbhost." -input_id ".$r;
	     
	     print STDERR "COMMAND: $command\n";
	 }
	 
	 if ($chunk == 2) {
	     my $dir = $scratchdir."/chunks";     
	     opendir(DIR, $dir);
	     my @allfiles = readdir DIR;
	     closedir DIR;
	     
	     foreach my $f(@allfiles) {
		 if (($f ne ".") && ($f ne "..")) {
		      my $command = "bsub -q ".$queue." -o ".$scratchdir."/".$r."/stdout/".$f.".out -o ".$scratchdir."/".$r."/stderr/".$f.".err ".$check." ".$runner." -runnable ".$runnabledb." -dbuser ".$dbuser." -pass ".$dbpass." -dbname ".$dbname." -host ".$dbhost." -input_id ".$f;
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
