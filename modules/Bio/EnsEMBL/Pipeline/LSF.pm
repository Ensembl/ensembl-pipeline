#
# Object for submitting jobs to and querying the LSF queue
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

Bio::EnsEMBL::Pipeline::LSF

=head1 SYNOPSIS

=head1 DESCRIPTION

Stores run and status details of an analysis job

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::EnsEMBL::Pipeline::LSF;

use vars qw(@ISA);
use strict;

use FreezeThaw qw(freeze thaw);

# Object preamble - inherits from Bio::Root::Object;

use Bio::EnsEMBL::Pipeline::LSFJob;

@ISA = qw(Bio::EnsEMBL::Pipeline::DB::JobI Bio::Root::Object);

sub _initialize {
    my ($self,@args) = @_;

    my $make = $self->SUPER::_initialize(@args);
    my ($user,$queue) = $self->_rearrange([qw(USER
					      QUEUE)],@args);

    $user  = "humpub"     unless defined($user);
    $queue = "blast_farm" unless defined($queue);

    $self->user ($user);
    $self->queue($queue);

    return $make; # success - we hope!
}

=head2 get_all_jobs

  Title   : get_all_jobs
  Usage   : my @jobs = get_all_jobs($user)
  Function: Gets all jobs owned by $user
  Returns : @Bio::EnsEMBL::Pipeline::LSFJob
  Args    : string

=cut


sub get_all_jobs {
    my ($self) = @_;

    my @jobs;

    my $cmd = "bjobs -w ";

    if (defined($self->user)) {
	$cmd .= " -u " . $self->user;
    }
    if (defined($self->queue)) {
	$cmd .= " -q " . $self->queue;
    }

    open(IN,$cmd . " |");

    
    while (<IN>) {
	if ($_ !~ /^JOBID/) {

	    my $job = $self->_parse_line($_);
	    push(@jobs,$job);
	    
	}
    }
    close(IN);

    return @jobs;
}

sub get_job {
    my ($self,$id) = @_;


    my $cmd = "bjobs -w ";

    if (defined($self->user)) {
	$cmd .= " -u " . $self->user;
    }

    $cmd .= " $id ";

    open(IN,$cmd . " |");
    
    while(<IN>) {
	if ($_ !~ /^JOBID/) {
	    my $job = $self->_parse_line($_);
	    return $job;
	}
    }
    close(IN);
}
	
sub _parse_line {
    my ($self,$line) = @_;

    my @f = split(' ',$line);
    
    my $id       = $f[0];
    my $user     = $f[1];
    my $stat     = $f[2];
    my $queue    = $f[3];
    my $from_host= $f[4];
    my $exec_host= $f[5];
    my $submit_time = $f[$#f-2] . " " .$f[$#f-1] . " " . $f[$#f];
    
    my $job_name = "";
    
    for (my $i = 6; $i < $#f-2; $i++) {
	$job_name = $job_name . $f[$i] . " ";
    }
    chomp($job_name);
    
    my $job = new Bio::EnsEMBL::Pipeline::LSFJob(-id              => $id,
						 -user            => $user,
						 -status          => $stat,
						 -queue           => $queue,
						 -from_host       => $from_host,
						 -exec_host       => $exec_host,
						 -job_name        => $job_name,
						 -submission_time => $submit_time);
    
    return $job;
}


sub user {
    my ($self,$arg) = @_;

    if (defined($arg)) {
	$self->{_user} = $arg;
    }

    return $self->{_user};
}

sub queue {
    my ($self,$arg) = @_;

    if (defined($arg)) {
	$self->{_queue} = $arg;
    }

    return $self->{_queue};

}

1;



