# Object for storing details of an analysis job
#
# Cared for by Michele Clamp  <michele@sanger.ac.uk>
#
# Copyright Michele Clamp
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::DBSQL::Job

=head1 SYNOPSIS

=head1 DESCRIPTION

Stores run and status details of an analysis job

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::EnsEMBL::Pipeline::Job;
# use Data::Dumper;

BEGIN {
  require "Bio/EnsEMBL/Pipeline/pipeConf.pl";
}

use vars qw(@ISA);
use strict;

# use FreezeThaw qw(freeze thaw);

# Object preamble - inherits from Bio::Root::Object;

use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::Pipeline::Status;
use Bio::EnsEMBL::Pipeline::DBSQL::JobAdaptor;
use Bio::EnsEMBL::Pipeline::DB::JobI;


@ISA = qw(Bio::EnsEMBL::Pipeline::DB::JobI Bio::Root::RootI);

# following vars are static and not meaningful on remote side
# recreation of Job object. Not stored in db of course.
# hash with queue keys

my %batched_jobs;
my %batched_jobs_runtime;


=head2 new

 Arg [1]   : InputId (contig_id)  
 Function  : Creates a new Bio::EnsEMBL::Pipeline::Job
 Returntype: Bio::EnsEMBL::Pipeline::Job
 Exceptions: Throws if not passed an inputId or a Bio::EnsEMBL::Analysis object 
 Caller    : Bio::EnsEMBL::Pipeline::Job->new();
 Example   : my $job = Bio::EnsEMBL::Pipeline::Job->new( -input_id => $contig->id,
	       					         -analysis => $analysis);

=cut

#this method creates the Job object it requires an input id and an analysis object
#it also can take a jobadaptor, lsf id, class, job id, stderr, stdout, input_object_file and retry count

sub new {
    my ($class, @args) = @_;
    my $self = bless {},$class;

    my ($adaptor,$dbID,$lsfid,$input_id,$cls,$analysis,$stdout,$stderr,$input, $retry_count ) 
        = $self->_rearrange([qw(ADAPTOR
                                ID
                                LSF_ID
                                INPUT_ID
                                CLASS
                                ANALYSIS
                                STDOUT
                                STDERR
                                INPUT_OBJECT_FILE
                                RETRY_COUNT
                                )],@args);

    $dbID    = -1 unless defined($dbID);
    $lsfid = -1 unless defined($lsfid);
    $cls = 'contig' unless defined($cls);

    $input_id   || $self->throw("Can't create a job object without an input_id");
    $analysis   || $self->throw("Can't create a job object without an analysis object");

    $analysis->isa("Bio::EnsEMBL::Analysis") ||
        $self->throw("Analysis object [$analysis] is not a Bio::EnsEMBL::Analysis");

    $self->dbID         ($dbID);
    $self->adaptor  ($adaptor);
    $self->input_id   ($input_id);
    $self->class   ($cls);
    $self->analysis   ($analysis);
    $self->stdout_file($stdout);
    $self->stderr_file($stderr);
    $self->input_object_file($input);
    $self->retry_count( $retry_count );
    $self->LSF_id( $lsfid );

    return $self;
}



=head2 create_by_analysis_inputId

  Arg [1]   : Analysis object 
  Function  : creates an job given an input_id and analysis object
  Returntype: Bio::EnsEMBL::Pipeline::Job
  Exceptions: None
  Caller    : Bio::EnsEMBL::Pipeline::Job->create_by_analysis_inputId($analysis, $contig->id);
  Example   : 

=cut


#this method calls the new method creating a new job object with the given contig id and analysis object
#this is the recommended way of creating jobs without db connections

sub create_by_analysis_inputId {
  my $dummy = shift;
  my $analysis = shift;
  my $inputId = shift;
  my $class = shift;

  my $job = Bio::EnsEMBL::Pipeline::Job->new
    ( -input_id    => $inputId,
      -analysis    => $analysis,
      -retry_count => 0,
      -class       => $class
    );
  $job->make_filenames;
  return $job;
}





#################
#get/set methods#
#################


=head2 accessor methods

  Arg [1]   : what ever is being set 
  Function  : to set return appropriate variable
  Returntype: the appropriate value
  Exceptions: only the Analysis method throws if it hasn't been passed an analysis object'
  Caller    : Bio::EnsEMBL::Job
  Example   : $dbId = $self->dbID();

=cut


sub dbID {
    my ($self,$arg) = @_;

    if (defined($arg)) {
        $self->{'_dbID'} = $arg;
    }
    return $self->{'_dbID'};

}

sub adaptor {
    my ($self,$arg) = @_;

    if (defined($arg)) {
        $self->{'_adaptor'} = $arg;
    }
    return $self->{'_adaptor'};

}

sub class {
    my ($self,$arg) = @_;

    if (defined($arg)) {
        $self->{'_class'} = $arg;
    }
    return $self->{'_class'};
}


sub input_id {
    my ($self,$arg) = @_;

    if (defined($arg)) {
        $self->{'_input_id'} = $arg;
    }
    return $self->{'_input_id'};
}

sub analysis {
    my ($self,$arg) = @_;
    if (defined($arg)) {
        $self->throw("[$arg] is not a Bio::EnsEMBL::Analysis object" ) 
            unless $arg->isa("Bio::EnsEMBL::Analysis");

        $self->{'_analysis'} = $arg;
    }
    return $self->{'_analysis'};

}

sub stdout_file {
    my ($self,$arg) = @_;

    if (defined($arg)) {
        $self->{'_stdout_file'} = $arg;
    }
    return $self->{'_stdout_file'};
}

sub stderr_file {
    my ($self,$arg) = @_;

    if (defined($arg)) {
        $self->{'_stderr_file'} = $arg;
    }
    return $self->{'_stderr_file'};
}

sub input_object_file {
    my ($self,$arg) = @_;

    if (defined($arg)) {
        $self->{'_input_object_file'} = $arg;
    }
    return $self->{'_input_object_file'};
}

sub LSF_id {
  my ($self, $arg) = @_;
  (defined $arg) &&
    ( $self->{'_lsfid'} = $arg );
  $self->{'_lsfid'};
}

##used for redirecting lsf output to stderr and stdout files###

sub LSF_out {
  my ($self, $arg) = @_;
  (defined $arg) &&
    ( $self->{'_lsfout'} = $arg );
  $self->{'_lsfout'};
}

sub LSF_err {
  my ($self, $arg) = @_;
  (defined $arg) &&
    ( $self->{'_lsferr'} = $arg );
  $self->{'_lsferr'};
}

sub retry_count {
  my ($self, $arg) = @_;
  (defined $arg) &&
    ( $self->{'_retry_count'} = $arg );
  $self->{'_retry_count'};
}



#############
#RUN METHODS#
#############



=head2 flush_runs

  Arg [1]   : Job adaptor object
  Function  : When there are a defined number of jobs ready to run (set in pipeconf or RuleManger) this submits them to lsf using the adaptor to connect to the database. If there is more than one job being submitted the last job is used for stdout/stderr
  Returntype: none 
  Exceptions: throws if has no job adaptor or if no runner script is defined in pipeconf.pl
  Caller    : Bio::EnsEMBL::Pipeline::Job
  Example   : $self->flush_runs( $self->adaptor, $LSF_params );

=cut

#This checks if more than one queue should be used then if there is a runner script defined
#once this has been established for each queue jobs should be run on the last job id is got and the bsub line constructed and submitted 
#and the job status is set to submitted or failed depending on whether the job has an lsfId

sub flush_runs {
  my $self = shift;

  my $adaptor = shift;
  my $LSF_params = shift;
  my @queues;
  
  my $nodes   = $LSF_params->{'nodes'}   || undef;
  my $queue   = $LSF_params->{'queue'}   || undef;
  my $jobname = $LSF_params->{'jobname'} || undef;

  if( !defined $adaptor ) {
    $self->throw( "Cannot run remote without db connection" );
  }

  local *SUB;
  local *FILE;

#not sure why this does as %batched_jobs is  created using lsf_params in the first place
  if( ! defined $queue ) {
    @queues = keys %batched_jobs;
  } else {
    @queues = ( $queue );
  }
  
  my $db       = $adaptor->db;
  my $host     = $db->host;
  my $username = $db->username;
  my $dbname   = $db->dbname;
  my $pass     = $db->password;
  my $lsfid;

  # runner.pl: first look in pipeConf.pl,
  # then in same directory as Job.pm,
  # and fail if not found

  my $runner = $::pipeConf{'runner'} || undef;
##checking if runner is defined
  unless (-x $runner) {
    $runner = __FILE__;
    $runner =~ s:/[^/]*$:/runner.pl:;
    $self->throw("runner undefined - needs to be set in pipeConf.pl\n") unless defined $runner;
  }

  for my $queue ( @queues ) {
##for each queue jobs are to be submitted to get the la
    if (! (defined $batched_jobs{$queue}) || ! scalar (@{$batched_jobs{$queue}})) {
      next;
    }

    my $lastjob = $adaptor->fetch_by_dbID( $batched_jobs{$queue}->[$#{$batched_jobs{$queue}}] );

    if( ! defined $lastjob ) {
      $self->throw( "Last batch job not in db" );
    }
  
    my $cmd;
  ## -C 0 important as it prevents core dumps from filling up /tmp/
    $cmd = "bsub -C 0 -o ".$lastjob->stdout_file;
    if ($nodes) {
        # $nodes needs to be a space-delimited list
        $nodes =~ s/,/ /;
        $nodes =~ s/ +/ /;
        # undef $nodes unless $nodes =~ m{(\w+\ )*\w};
        $cmd .= " -m '$nodes' ";
    }
    $cmd .= " -q $queue " if defined $queue;
    $cmd .= " -J $jobname " if defined $jobname;
    $cmd .= " -r -e ".$lastjob->stderr_file." -E \"$runner -check\" ";

    # check if the password has been defined, and write the
    # "connect" command line accordingly (otherwise -pass gets the
    # first job id as password, instead of remaining undef)

    if ($pass) {
      $cmd .= $runner." -host $host -dbuser $username -dbname $dbname -pass $pass ".join( " ",@{$batched_jobs{$queue}} );
    }
    else {
      $cmd .= $runner." -host $host -dbuser $username -dbname $dbname ".join( " ",@{$batched_jobs{$queue}} );
    }
    
    print STDERR "$cmd\n";
    open (SUB,"$cmd 2>&1 |");
  
    while (<SUB>) {
      if (/Job <(\d+)>/) {
        $lsfid = $1;
      }
    }
    close(SUB);

    if( ! defined $lsfid ) {
      print STDERR ( "Couldnt submit ".join( " ",@{$batched_jobs{$queue}} )." to LSF" );
      foreach my $jobid ( @{$batched_jobs{$queue}} ) {
        my $job = $adaptor->fetch_by_dbID( $jobid );
        $job->set_status( "FAILED" );
      }
    } else {
    
      foreach my $jobid ( @{$batched_jobs{$queue}} ) {
        my $job = $adaptor->fetch_by_dbID( $jobid );
        if( $job->retry_count > 0 ) {
          for ( $job->stdout_file, $job->stderr_file ) {
            open( FILE, ">".$_ ); close( FILE );
          }
        }
        $job->LSF_id( $lsfid );
        # $job->create_lsflogfile;
        $job->retry_count( $job->retry_count + 1 );
        $job->set_status( "SUBMITTED" );
        $adaptor->update( $job );
      }
    }
    $batched_jobs{$queue} = [];
    $batched_jobs_runtime{$queue} = 0;
  }
}


=head2 batch_runRemote

  Arg [1]   : reference to hash of LSF_parameters
  Function  : create hash of jobids for a particular queue and submit these to flush_run
  Returntype: none
  Exceptions: none
  Caller    : Bio::EnsEMBL::Pipeline::Job
  Example   : $cj->batch_runRemote($LSF_params) (from RuleManger3.pl)

=cut



sub batch_runRemote {
  my $self = shift;
  my $LSF_params = shift;

  my $batchsize = $LSF_params->{'flushsize'} || 1;
  my $queue     = $LSF_params->{'queue'};
  
  # should check job->analysis->runtime
  # and add it to batched_jobs_runtime
  # but for now just
  push( @{$batched_jobs{$queue}}, $self->dbID );
  if ( scalar( @{$batched_jobs{$queue}} ) >= $batchsize ) {
    $self->flush_runs( $self->adaptor, $LSF_params );
  }
}


=head2 runLocally

  Arg [1]   : none 
  Function  : Runs a job locally rather than submitting to LSF
  Returntype: none
  Exceptions: none
  Caller    : Bio:EnsEMBL::Pipeline::Job
  Example   : $job->runLocally;

=cut

#if jobs are to be run locally it can be defined in rulemanager or pipeconf
sub runLocally {
  print STDERR "Running locally\n"; 
  my $self = shift;

  print STDERR "Running locally " . $self->stdout_file . " " . $self->stderr_file . "\n"; 

  local *STDOUT;
  local *STDERR;

  if( ! open ( STDOUT, ">".$self->stdout_file )) {
    $self->set_status( "FAILED" );
    return;
  }
        
  if( ! open ( STDERR, ">".$self->stderr_file )) {
    $self->set_status( "FAILED" );
    return;
  }
       print STDERR "Running inLSF\n"; 
  $self->runInLSF();
}



=head2 runRemote

  Arg [1]   : $queue, scalar of which queue should be used 
  Function  : submit single jobs to lsf
  Returntype: none
  Exceptions: thorws if has no job adaptor, no runner script or if useDB is set to 0
  Caller    : Bio::EnsEMBL::Pipeline::Job
  Example   : $self->runRemote($queue, $useDB);

=cut

#this does bacially the same things as flush_run but for one job with one queue but is probably never used

sub runRemote {
  my $self = shift;
  my $queue = shift;
  my $useDB = shift;

  local *SUB;
  local *FILE;

  if( !defined $self->adaptor ) {
    $self->throw( "Cannot run remote without db connection" );
  }

  my $db = $self->adaptor->db;
  my $host = $db->host;
  my $username = $db->username;
  my $dbname = $db->dbname;
  my $cmd;

  # runner.pl: first look in same directory as Job.pm
  # if it's not here use file defined in pipeConf.pl
  # otherwise fail

  my $runner = __FILE__;
  $runner =~ s:/[^/]*$:/runner.pl:;     

  unless (-x $runner) {
    $runner = $::pipeConf{'runner'} || undef;
    $self->throw("runner undefined - needs to be set in pipeConf.pl\n") unless defined $runner;
  }

  $cmd = "bsub -C0 -q ".$queue." -o ".$self->stdout_file.
#    " -q acarichunky " .
#    " -R osf1 ".
    " -e ".$self->stderr_file." -E \"$runner -check\" ";


  if( ! defined $useDB ) {
    $useDB = 1;
  }

  if( $self->retry_count > 0 ) {
    for ( $self->stdout_file, $self->stderr_file ) {
      open( FILE, ">".$_ ); close( FILE );
    }
  }

  if( $useDB ) {
    # find out db details from adaptor
    # generate the lsf call
    $cmd .= $runner." -host $host -dbuser $username -dbname $dbname ".$self->dbID;
    
  } else {
    # make the object
    # call its get input method..
    # freeze it
    $self->throw( "useDB=0 not implemented yet." );
  }
  print STDERR "$cmd\n";
  open (SUB,"$cmd 2>&1 |");
  
  while (<SUB>) {
    if (/Job <(\d+)>/) {
      $self->LSF_id($1);
      # print (STDERR $_);
    }
  }
  close(SUB);

  if( $self->LSF_id == -1 ) {
    print STDERR ( "Couldnt submit ".$self->dbID." to LSF" );
  } else {
    $self->retry_count( $self->retry_count + 1 );
    $self->set_status( "SUBMITTED" );
  }

  $self->adaptor->update( $self );
}


=head2 runInLSF

  Arg [1]   : name of runnabledb to be used 
  Function  : This creates and tries run a RunnableDB for the specified analysis
  Returntype: none
  Exceptions: throws if problem in any stage of running the RunnableDB or setting the status 
  Caller    : Bio::EnsEMBL::Pipeline::Job
  Example   : $self->runInLSF;

=cut



sub runInLSF {
  my $self = shift;
  my $module = $self->analysis->module;
  my $rdb;
  my $err;
  my $autoupdate = $::pipeConf{'autoupdate'};
  

  eval {
      if( $module =~ /::/ ) {
          $module =~ s/::/\//g;
          require "${module}.pm";
          $rdb = "${module}"->new
              ( -analysis => $self->analysis,
                -input_id => $self->input_id,
                -dbobj => $self->adaptor->db );
      } else {
          require "Bio/EnsEMBL/Pipeline/RunnableDB/${module}.pm";
          $module =~ s/\//::/g;
          $rdb = "Bio::EnsEMBL::Pipeline::RunnableDB::${module}"->new
              ( -analysis => $self->analysis,
                -input_id => $self->input_id,
                -dbobj => $self->adaptor->db );
      }
  };
  if ($err = $@) {
      print (STDERR "CREATE: Lost the will to live Error\n");
      $self->set_status( "FAILED" );
      $self->throw( "Problems creating runnable $module for " . $self->input_id . " [$err]\n");
  }
  eval {   
      $self->set_status( "READING" );
      $rdb->fetch_input;
  };
  if ($err = $@) {
      $self->set_status( "FAILED" );
      print (STDERR "READING: Lost the will to live Error\n");
      die "Problems with $module fetching input for " . $self->input_id . " [$err]\n";
  }
  if ($rdb->input_is_void) {
      $self->set_status( "VOID" );
      return;
  }
  eval {
      $self->set_status( "RUNNING" );
      $rdb->run;
  };
  if ($err = $@) {
      $self->set_status( "FAILED" );
      print (STDERR "RUNNING: Lost the will to live Error\n");
      die "Problems running $module for " . $self->input_id . " [$err]\n";
  }
  eval {
      $self->set_status( "WRITING" );
      $rdb->write_output;
      $self->set_status( "SUCCESSFUL" );
  }; 
  if ($err = $@) {
      $self->set_status( "FAILED" );
      print (STDERR "WRITING: Lost the will to live Error\n");
      die "Problems for $module writing output for " . $self->input_id . " [$err]" ;
  }
  if ($autoupdate) {
    eval {
      my $sic = $self->adaptor->db->get_StateInfoContainer;
      $sic->store_inputId_class_analysis(
        $self->input_id,
        $self->class,
        $self->analysis
      );
    };
    if ($err = $@) {
      print STDERR "Error updating successful job ".$self->dbID ."[$err]\n";
    }
    else {
      print STDERR "Updated successful job ".$self->dbID."\n";
      eval {
        $self->remove;
      };
      if ($err = $@) {
         print STDERR "Error deleting job ".$self->dbID." [$err]\n";
      }
      else {
         print STDERR "Deleted job ".$self->dbID."\n";
      }
    }
  }
}


################
#SETTING STATUS#
################



=head2 setting statuses

  Arg [1]   : scalar of status to be set 
  Function  : calls status setting methods in jobadaptor
  Returntype: a status object or an array of status objects  
  Exceptions: thows if no db connection or no status is passed
  Caller    : Bio::EnsEMBL::Pipeline::Job
  Example   : $self->set_status( "FAILED" );

=cut


sub set_status {
  my ($self,$arg) = @_;
  
  $self->throw("No status input" ) unless defined($arg);
  
  
  if (!(defined($self->adaptor))) {
    $self->warn("No database connection.  Can't set status to $arg");
    return;
  }
  
  return $self->adaptor->set_status( $self, $arg );
}

sub current_status {
  my ($self,$arg) = @_;
  
  if( ! defined( $self->adaptor )) {
    return undef;
  }
 
  return $self->adaptor->current_status( $self, $arg );
}

sub get_all_status {
  my ($self) = @_;
  
  if( $self->adaptor ) {
    return $self->adaptor->get_all_status( $self );
  } else {
    return undef;
  }
}


sub get_last_status {
  my ($self) = @_;
  
  if( $self->adaptor ) {
    return $self->adaptor->get_last_status( $self );
  } else {
    return undef;
  }
}

####################
#OUTPUTFILE METHODS#
####################

=head2 make_filenames

  Arg [1]   : none 
  Function  : creates the stdout and stderr file names
  Returntype: none
  Exceptions: none
  Caller    : $self;
  Example   : $self->make_filenames;

=cut



sub make_filenames {
  my ($self) = @_;
  
  my $num = int(rand(10));
  my $dir = $::pipeConf{'nfstmp.dir'} . "/$num/";
  if( ! -e $dir ) {
    system( "mkdir $dir" );
  }

# scp - one set of out files per job (even if batching together)
# this is a bit messy! added '.0' to $stub. This will be the master
# file containing LSF output. In runner.pl before each job is run
# replace 0 with the job ID to get one output file per $job.
# Change also Job::remove to do a glob on all these files. Yep it's
# nasty but it seems to work...


  my $stub = $self->input_id.".";
  $stub .= $self->analysis->logic_name.".";
  $stub .= time().".".int(rand(1000));
#   $stub .= time().".".int(rand(1000)) . '.0';

#  my @files = $self->get_files( $self->input_id,"obj","out","err" );
#  for (@files) {
#    open( FILE, ">".$_ ); close( FILE );
#  }
 
  $self->input_object_file($dir.$stub.".obj");
  $self->stdout_file($dir.$stub.".out");
  $self->stderr_file($dir.$stub.".err");
}




=head2 remove

  Arg [1]   : none
  Function  : deletes the output files
  Returntype: none
  Exceptions: none
  Caller    : Bio::EnsEMBL::Pipeline::Job
  Example   : $job->remove 

=cut


sub remove {
  my $self = shift;
  
  if( -e $self->stdout_file ) { unlink( $self->stdout_file ) };
  if( -e $self->stderr_file ) { unlink( $self->stderr_file ) };
  if( -e $self->input_object_file ) { unlink( $self->input_object_file ) };

   if( defined $self->adaptor ) {
   $self->adaptor->remove( $self );
   }
}

################
#UNUSED METHODS#
################


#this was used when complex data structure were to be stored in file similariy to the web ini file but it never actually worked 


sub write_object_file {
    my ($self,$arg) = @_;

    $self->throw("No input object file defined") unless defined($self->input_object_file);

    if (defined($arg)) {
        my $str = FreezeThaw::freeze($arg);
        open(OUT,">" . $self->input_object_file) || $self->throw("Couldn't open object file " . $self->input_object_file);
        print(OUT $str);
        close(OUT);
    }
}


####appears to be never used and obviously was never implemented 
##Simon not sure what it was for and why it was never implemented

sub resultToDb {
  my $self = shift;
  $self->throw( "Not implemented yet." );
}



#no  longer used so could potentially be removed as make_filenames seems to do the same thing

sub create_lsflogfile {
  my ($self) = @_;
  
  my $num = int(rand(10));
  my $dir = $::pipeConf{'nfstmp.dir'} . "/$num/";
  if( ! -e $dir ) {
    system( "mkdir $dir" );
  }

  my $stub = $self->LSF_id.".";
  $stub .= time().".".int(rand(1000));

  $self->LSF_out($dir.$stub.".out");
  $self->LSF_err($dir.$stub.".err");
}

# creates fileanames and does something else i don't understand be doesn't seem to be used at all

sub get_files {
  my ($self,$stub,@exts) = @_;
  my @result;
  my @files;
  my $dir;

  my $count = 0;

  while( 1 ) {
    my $num = int(rand(10));
    $dir = $::pipeConf{'nfstmp.dir'} . "/$num/";
    if( ! -e $dir ) {
      system( "mkdir $dir" );
    }
    opendir( DIR, $dir );
    @files = readdir( DIR );
    if( scalar( @files ) > 10000 ) {
      if( $count++ > 10 ) {
        $self->throw("10000 files in directory. Can't make a new file");
      } else {
        next;
      }
    } else {
      last;
    }
  }
  $count = 0;
  # Should check disk space here.
  
 OCCU: while( $count++ < 10 ) { 
    my $rand  = int(rand(100000));
    @result = ();
    foreach my $ext ( @exts ) {
      my $file  = $dir . $stub . "." . $rand . "." . $ext;
      if( -e $file ) {
        next OCCU;
      }
      push( @result, $file );
    }
    last;
  }

  if( $count > 9 ) {
    $self->throw( "File generation problem in $dir" );
  }
  return @result;
}


1;
