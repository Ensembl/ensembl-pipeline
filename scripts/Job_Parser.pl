#!/usr/local/bin/perl

=head1 NAME

Job_Parser - Parses the output of submitted jobs

=head1 SYNOPSIS

    Job_Submitter

=head1 DESCRIPTION

Job_Parser does the following:

Looks in the analysis db for submitted jobs and for each job:

1)Checks status from status file

2)parses the STDOUT_file

3)Assigns the correct status (SUCCESSFUL, FAILED) in the Job_status table

4)Recreates the original Job module (e.g.Est2Genome)

5)Uses the Job module to parse the data output, 
i.e. thaws the frozen object containing the data

6)Uses the Job module to write the data to the database (e.g. features of contigs)

8)Finally, if all ok, sets status to DONE


=head1 OPTIONS

    -host      host name for database (gets put as host= in locator)

    -port      For RDBs, what port to connect to (port= in locator)

    -dbname    For RDBs, what name to connect to (dbname= in locator)

    -dbuser    For RDBs, what username to connect as (dbuser= in locator)

    -pass      For RDBs, what password to use (dbpass= in locator)

=cut


use strict;

use Getopt::Long;

use Bio::EnsEMBL::DBSQL::Obj;
use Bio::EnsEMBL::Pipeline::DBSQL::Obj;
use Bio::EnsEMBL::Pipeline::DBSQL::Job;
use Bio::EnsEMBL::Pipeline::SimpleJob;
use Bio::EnsEMBL::Pipeline::RunnableDBI;
use Bio::EnsEMBL::Pipeline::RunnableDB::Est2Genome;
use Bio::EnsEMBL::Pipeline::RunnableDB::Vert_Est2Genome;

my $anahost     = 'localhost';
my $anaport     = '410000';
my $anadbname   = 'pipeline';
my $anadbuser   = 'root';
my $anapass     =  undef;

&GetOptions(
	    'anahost:s'     => \$anahost,
	    'anaport:n'     => \$anaport,
	    'anadbname:s'   => \$anadbname,
	    'anadbuser:s'   => \$anadbuser,
	    'anapass:s'     => \$anapass,
	    );



my $anadb  = new Bio::EnsEMBL::Pipeline::DBSQL::Obj(-host   => $anahost,
						    -port   => $anaport,
						    -dbname => $anadbname,
						    -user   => $anadbuser,
						    -pass   => $anapass);


my @jobs = $anadb->get_JobsByCurrentStatus('SUBMITTED');

JOB:foreach my $job (@jobs) {
    
    $job->_dbobj($anadb);

    my $module   = $job->analysis->module;
    my $status   = $job->status_file;
    my $stdout   = $job->stdout_file();

    print(STDERR "Id is " . $job->id . "\t$module\t" . $job->input_id . "\n");


    next JOB if  is_processed($job);     
    next JOB if  -z $stdout;              # Need something in stdout to parse
    
    open(OUT,"<$stdout") || next;
    my $failed = 1;
    while(<OUT>){
	if(/Successfully/){
	    $failed=0;
	    $job->set_status('SUCCESSFUL');
	    print STDERR "Successful job!\n";  
	}
	if(/Read file <(\S+)> for stderr output/){
	    $job->stderr_file($1);
	}
    }
    close OUT;
    if ($failed == 1) {
	$job->set_status('FAILED');
	my $error=$job->stderr_file();
	print STDERR "Failed job! (STDERR file $error)\n";
	open(OUT,"<$error") || do {print STDERR "Could not open $error\n"; next;};
	while(<OUT>){
	    if(/Attempting to set the sequence to \[/){
		print STDERR "bad sequence setting\n";
		last;
	    }
	    if(/Error running \'est_genome/){
		print STDERR "Error runnning est2genome\n";
		last;
	    }
	    if(/MSG: Not enough disk space/){
		print STDERR "MSG: Not enough disk space\n";
		last;
	    }
	     if(/Can\'t use an undefined value as an ARRAY reference/){
		print STDERR "Empty array reference\n";
		last;
	    }
	}
	close OUT;
    }
    else {
	my $output= $job->output_file();
	print STDERR "Got output file $output\n\n";
	
	my $runnable = "$module"->new(-input_id => $job->input_id,
				      -dbobj    => $anadb);
	my @features=$runnable->fetch_output($output);
	print STDERR "Got features from frozen output (when it's fixed)\n";
	#$runnable->write_output(@features);
	#$job->set_status('DONE');    
    }
    
}


sub is_processed {
    my ($job) = @_;
    
    open(OUT,"<" . $job->status_file) || return 0;
    
    my $processed = 0;

    while(<OUT>){
	print $_;
	if(/PROCESSED/){
	    print STDERR "Found processed!\n";
	    $job->set_status('PROCESSED');
	    $processed = 1;
	}
    }

    close OUT;

    return $processed;
}
