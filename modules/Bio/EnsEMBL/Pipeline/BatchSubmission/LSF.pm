=head1 LICENSE

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2020] EMBL-European Bioinformatics Institute
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#      http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

=head1 NAME

  Bio::EnsEMBL::Pipeline::BatchSubmission::LSF

=head1 SYNOPSIS

  my $batchjob = Bio::EnsEMBL::Pipeline::BatchSubmission::LSF->new(
              -STDOUT     => $stdout_file,
              -STDERR     => $stderr_file,
              -PARAMETERS => @args,
              -PRE_EXEC   => $pre_exec,
              -QUEUE      => $queue,
              -JOBNAME    => $jobname,
              -NODES      => $nodes,
              -RESOURCE   => $resource
              );

  $batch_job->construct_command_line('test.pl');
  $batch_job->open_command_line();

=head1 DESCRIPTION

  This module provides an interface to the Platform LSF load
  sharing software and its commands. It implements the method
  construct_command_line which is not defined in the base class
  and which enables the pipeline to submit jobs in a distributed
  environment using LSF.

  See base class Bio::EnsEMBL::Pipeline::BatchSubmission for more info.

=head1 APPENDIX

  The rest of the documentation details each of the object methods.
  Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Pipeline::BatchSubmission::LSF;

use warnings ;
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning info);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Pipeline::BatchSubmission;
use File::Copy;
use vars qw(@ISA);
use strict;

@ISA = qw(Bio::EnsEMBL::Pipeline::BatchSubmission);


sub new {
  my ( $class, @args ) = @_;
  my $self = $class->SUPER::new(@args);
  #print "CREATING ".$self." with ".join(" ", @args)."\n";

  return $self;

}

# Rather than the quite cumbersome
#   resource => 'select[mem>500] rusage[mem=500]',
#   sub_args => '-M 500000'
# we allow for
#   memory   => '500MB' # or '0.5GB' or '500000KB'
#
sub memstring_to_resource {
  my ( $self, $memstring ) = @_;

  my $resource_mem  = $self->memory_conversion( $memstring, 'MB' );
  # for farm 3
  my $softlimit_mem = $self->memory_conversion( $memstring, 'MB' );
  #my $softlimit_mem = $self->memory_conversion( $memstring, 'KB' );

  my $resource = $self->resource();

  my $new_resource_string = sprintf( "select[mem>%d] rusage[mem=%d]",
                                     $resource_mem, $resource_mem );

  if ( defined($resource) ) {
    if ( $resource =~ /^-R/ ) {
      # The resource string might already be prefixed with '-R', in this
      # case just add the memory resource with another '-R'.
      $self->resource(
             sprintf( "%s -R '%s'", $resource, $new_resource_string ) );
    }
    else {
      # ... otherwise just tag the resource onto the end of the string.
      $self->resource(
                  sprintf( "-R '%s' -R '%s'", $resource, $new_resource_string ) );
    }
  }
  else {
    # No previous resource spec.
    $self->resource($new_resource_string);
  }

  my $parameters = $self->parameters();

  my $new_parameter_string = sprintf( "-M %d", $softlimit_mem );

  if ( defined($parameters) ) {
    # Add the soft memory limit to the end of the parameter string.
    $self->parameters(
               sprintf( "%s %s", $parameters, $new_parameter_string ) );
  }
  else {
    # No previous parameter.
    $self->parameters($new_parameter_string);
  }
} ## end sub memstring_to_resource



######################
#command line methods#
######################

sub construct_command_line {
  my ( $self, $command, $stdout, $stderr ) = @_;
  # Command must be the first argument then if stdout or stderr aren't
  # definhed the objects own can be used

  if ( !$command ) {
    throw("cannot create bsub if nothing to submit to it : $!\n");
  }
  my $bsub_line;
  $self->command($command);
  if ($stdout) {
    $bsub_line = "bsub -o " . $stdout;
  } else {
    $bsub_line = "bsub -o " . $self->stdout_file;
  }
  if ( $self->nodes ) {
    my $nodes = $self->nodes;
    # $nodes needs to be a space-delimited list
    $nodes =~ s/,/ /;
    $nodes =~ s/ +/ /;
    # undef $nodes unless $nodes =~ m{(\w+\ )*\w};
    $bsub_line .= " -m '" . $nodes . "' ";
  }
  if ( my $res = $self->resource ) {
    $res = qq{-R '$res'} unless $res =~ /^-R/;
    $bsub_line .= " $res ";
  }
  $bsub_line .= " -q " . $self->queue         if $self->queue;
  $bsub_line .= " -J " . $self->jobname       if $self->jobname;
  $bsub_line .= " " . $self->parameters . " " if ( $self->parameters );
  if ($stderr) {
    $bsub_line .= " -e " . $stderr;
  } else {
    $bsub_line .= " -e " . $self->stderr_file;
  }
  $bsub_line .= " -E \"" . $self->pre_exec . "\"" if defined $self->pre_exec;
  ## must ensure the prexec is in quotes ##
  $bsub_line .= " " . $command;
  $self->command($bsub_line);

} ## end sub construct_command_line

sub open_command_line {
  my ( $self, $verbose ) = @_;

  my $lsf = '';

  if ( open( my $pipe, '-|' ) ) {
    while (<$pipe>) {
      next if /Mapping .* project definition to group/;
      if (/Job <(\d+)>/) {
        $lsf = $1;
      } else {
        warning("DEBUG: unexpected from bsub: '$_'");
      }
    }
    if ( close $pipe ) {
      if ( ( $? >> 8 ) == 0 ) {
        if ($lsf) {
          $self->id($lsf);
        } else {
          warning("Bsub worked but returned no job ID. Weird");
        }
      } else {
        throw( "Bsub failed : exit status " . $? >> 8 . "\n" );
      }
    } else {
      throw("Could not close bsub pipe : $!\n");
    }
  } else {
    # We want STDERR and STDOUT merged for the bsub process
    # open STDERR, '>&STDOUT';
    # probably better to do with shell redirection as above can fail
    exec( $self->command . ' 2>&1' ) || throw("Could not run bsub");
  }
} ## end sub open_command_line

sub get_pending_jobs {
  my ( $self, %args ) = @_;

  my ($user)    = $args{'-user'}    || $args{'-USER'}    || undef;
  my ($queue)   = $args{'-queue'}   || $args{'-QUEUE'}   || undef;
  my ($jobname) = $args{'-jobname'} || $args{'-JOBNAME'} || undef;

  my $cmd = "bjobs";
  $cmd .= " -q $queue"   if $queue;
  $cmd .= " -J $jobname" if $jobname;
  $cmd .= " -u $user"    if $user;
  $cmd .= " | grep -c PEND ";

  print STDERR "$cmd\n" if $args{'-debug'};

  my $pending_jobs = 0;
  if ( my $pid = open( my $fh, '-|' ) ) {
    eval {
      local $SIG{ALRM} = sub { kill 9, $pid; };
      alarm(60);
      while (<$fh>) {
        chomp;
        $pending_jobs = $_;
      }
      close $fh;
      alarm 0;
    };
  } else {
    exec($cmd );
    die q{Something went wrong here $!: } . $! . "\n";
  }
  print STDERR "FOUND $pending_jobs jobs pending\n" if $args{'-debug'};
  return $pending_jobs;
} ## end sub get_pending_jobs

sub get_job_time {
  my ($self) = @_;
  my $command = "bjobs -l";
  #print $command."\n";
  my %id_times;
  local *BJOB;
  open( BJOB, "$command |" ) or throw("couldn't open pipe to bjobs");
  my $job_id;
  while (<BJOB>) {
    chomp;
    if (/Job\s+\<(\d+)\>/) {
      $job_id = $1;
    } elsif (/The CPU time used/) {
      my ($time) = $_ =~ /The CPU time used is (\d+)/;
      $id_times{$job_id} = $time;
    }
  }
  close(BJOB);
  #or throw("couldn't close pipe to bjobs");
  return \%id_times;
}

sub check_existance {
  my ( $self, $id_hash, $verbose ) = @_;

  my %job_submission_ids = %$id_hash;
  my $command            = "bjobs";

  local *BJOB;
  open( BJOB, "$command 2>&1 |" )
    or throw("couldn't open pipe to bjobs");
  my %existing_ids;
LINE: while (<BJOB>) {
    print STDERR if ($verbose);
    chomp;
    if ( $_ =~ /No unfinished job found/ ) {
      last LINE;
    }
    my @values = split;
    if ( $values[0] =~ /\d+/ ) {
      if ( $values[2] eq 'UNKWN' ) {
        next LINE;
      }
      $existing_ids{ $values[0] } = 1;
    }
  }
  my @awol_jobs;
  foreach my $job_id ( keys(%job_submission_ids) ) {
    if ( !$existing_ids{$job_id} ) {
      push( @awol_jobs, @{ $job_submission_ids{$job_id} } );
    }
  }
  close(BJOB);
  #or throw("Can't close pipe to bjobs");
  return \@awol_jobs;
} ## end sub check_existance



sub job_stats {
  my ( $self, $verbose ) = @_;
  my $command = "bjobs";

  # Need to sleep to make sure lsf is displaying correct info
  sleep(20);

  local *BJOB;
  open( BJOB, "$command 2>&1 |" )
    or throw( "Couldn't open pipe to bjobs - "
            . "make sure you're using the correct BatchSubmision-module (LSF,GridEngine..)"
    );

  my %jobs;
LINE:
  while (<BJOB>) {
    chomp;
    if ( $_ =~ /No unfinished job found/ ) {
      last LINE;
    }
    my @values = split;
    $jobs{ $values[0] } = $values[2];
  }
  return \%jobs;
} ## end sub job_stats


#sub check_existance{
#  my ($self, $id, $verbose) = @_;
#  if(!$id){
#    die("Can't run without an LSF id");
#  }
#  my $command = "bjobs ".$id;
#  #print STDERR "Running ".$command."\n";
#  my $flag = 0; 
#  open(BJOB, "$command 2>&1 |") or throw("couldn't open pipe to bjobs");
#  while(<BJOB>){
#    print STDERR if($verbose);
#    chomp;
#    if ($_ =~ /No unfinished job found/) {
#      #print "Set flag\n";
#      $flag = 1;
#    } 
#    my @values = split;
#    if($values[0] =~ /\d+/){
#      return $values[0];
#    }
#  }
#  close(BJOB);
#  print STDERR "Have lost ".$id."\n" if($verbose);
#  return undef;
#}


sub kill_job {
  my ( $self, $job_id ) = @_;

  my $command = "bkill " . $job_id;
  system($command);
}

sub stdout_file {
  my ( $self, $arg ) = @_;

  if ($arg) {
    $self->{'stdout'} = $arg;
  }

  if ( !$self->{'stdout'} ) {
    $self->{'stdout'} = '/dev/null';
  }
  return $self->{'stdout'};
}

sub stderr_file {
  my ( $self, $arg ) = @_;

  if ($arg) {
    $self->{'stderr'} = $arg;
  }
  if ( !$self->{'stderr'} ) {
    $self->{'stderr'} = '/dev/null';
  }
  return $self->{'stderr'};
}

sub temp_filename {
  my ($self) = @_;

  $self->{'tmp_jobfilename'} = $ENV{'LSB_JOBFILENAME'};
  return $self->{'tmp_jobfilename'};
}

sub temp_outfile {
  my ($self) = @_;

  $self->{'_temp_outfile'} = $self->temp_filename . ".out";

  return $self->{'_temp_outfile'};
}

sub temp_errfile {
  my ($self) = @_;

  $self->{'_temp_errfile'} = $self->temp_filename . ".err";

  return $self->{'_temp_errfile'};
}

sub submission_host {
  my ($self) = @_;

  $self->{'_submission_host'} = $ENV{'LSB_SUB_HOST'};

  return $self->{'_submission_host'};
}

sub lsf_user {
  my ($self) = @_;

  $self->{'_lsf_user'} = $ENV{'LSFUSER'};

  return $self->{'_lsf_user'};
}

=head2 copy_output

  copy_output is used to copy the job's STDOUT and STDERR files using
  B<lsrcp>. This avoids using NFS'.

=cut

sub copy_output {
  my ( $self, $dest_err, $dest_out ) = @_;

  $dest_err ||= $self->stderr_file;
  $dest_out ||= $self->stdout_file;

  if ( !$self->temp_filename ) {
    my ( $p, $f, $l ) = caller;
    warning(   "The lsf environment variable LSB_JOBFILENAME is not defined"
             . " we can't copy the output files which don't exist $f:$l" );
    return;
  }

  # Unbuffer STDOUT so that data gets flushed to file
  # (It is OK to leave it unbuffered because this method
  # gets called after the job is finished.)
  my $old_fh = select(STDOUT);
  $| = 1;
  select($old_fh);

  my $temp_err = $self->temp_errfile;
  my $temp_out = $self->temp_outfile;

  my $command = $self->copy_command;
  my $remote  = $self->lsf_user . '@' . $self->submission_host;
  foreach my $set ( [ $temp_out, $dest_out ], [ $temp_err, $dest_err ] ) {
    my ( $temp, $dest ) = @$set;
    if ( -e $temp ) {
      if ( $command eq 'cp' || $dest =~ /^\/lustre/ ) {
        copy( $temp, $dest );
      } else {
        my $err_copy = "$command $temp $remote:$dest";
        unless ( system($err_copy) == 0 ) {
          warn "Error: copy '$err_copy' failed exit($?)";
        }
      }
    } else {
      warn "No such file '$temp' to copy\n";
    }
  }
} ## end sub copy_output

sub delete_output {
  my ($self) = @_;

  unlink $self->temp_errfile if ( -e $self->temp_errfile );
  unlink $self->temp_outfile if ( -e $self->temp_outfile );
}

sub copy_command {
  my ( $self, $arg ) = @_;

  if ($arg) {
    $self->{'_copy_command'} = $arg;
  }

  return $self->{'_copy_command'} || 'lsrcp ';
}

sub is_db_overloaded {
  my ( $self, $load_pending_cost ) = @_;

  my $resource = $self->resource;
  my ($select) = $resource =~ /select\[([^]]+)/;
  my ($rusage) = $resource =~ /rusage\[([^]]+)/;
  return 0 unless ( defined $select );
  my @a_dbs = $select =~ /my(\w+)\s*\W+\s*(\d+)/g;
  return 0 unless (@a_dbs);
  use Bio::EnsEMBL::Analysis::Config::Databases;
  for ( my $i = 0 ; $i < @a_dbs ; $i += 2 ) {
    my ( $host, $user, $passwd, $port );
    foreach my $db ( values %$DATABASES ) {
      next unless ( $db->{'-host'} eq $a_dbs[$i] );
      $host   = $db->{'-host'};
      $user   = $db->{'-user'};
      $passwd = $db->{'-pass'};
      $port   = $db->{'-port'};
    }
    my $dsn = "DBI:mysql:database=mysql;host=$host;port=$port";
    my $dbh = DBI->connect( $dsn, $user, $passwd )
      or die "Couldn't connect to database: " . DBI->errstr;

    my $sth = $dbh->prepare( 'SHOW STATUS WHERE Variable_name = "Threads_connected" OR Variable_name = "Queries"')
      or die "Couldn't prepare statement: " . $dbh->errstr;
    $sth->execute();
    my $t1 = time;
    sleep(1);
    my $t = $sth->fetchall_arrayref();
    $sth = $dbh->prepare("SHOW STATUS LIKE \'Queries\'")
      or die "Couldn't prepare statement: " . $dbh->errstr;
    $sth->execute();
    my @q           = $sth->fetchrow_array();
    my $time_diff   = time - $t1;
    my $queries_num = ( $q[1] - $t->[0][1] - 1 )/$time_diff;
    my $db_load =
      $t->[1][1] +
      ( $queries_num*10 ) +
      ( $self->get_pending_jobs( ( '-JOBNAME' => $self->jobname ) )*
        $load_pending_cost );
    return 1 if ( $a_dbs[ $i + 1 ] < $db_load );
  } ## end for ( my $i = 0 ; $i < ...
  return 0;
} ## end sub is_db_overloaded

1;
