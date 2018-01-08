=head1 LICENSE

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2018] EMBL-European Bioinformatics Institute
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

=head1 AUTHORS

Michele Clamp  <michele@sanger.ac.uk>

=head1 NAME

Bio::EnsEMBL::Pipeline::Status - small object storing job status tags

=head1 SYNOPSIS

    my $obj    = new Bio::EnsEMBL::Pipeline::Status
    ('-jobid'              => $jobid,
     '-status'             => $status,
     '-created'            => $created,
     );

=head1 DESCRIPTION

Stores the status of a job at a certain time

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::EnsEMBL::Pipeline::Status;
use warnings ;
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning info);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );

use vars qw(@ISA);
use strict;


@ISA = qw();


sub new {
  my($class,@args) = @_;
  
   my $self = bless {},$class;

  my ($jobid,$status,$created)  =
      rearrange([qw(JOBID
			    STATUS
			    CREATED
			    )],@args);

  $jobid   || $self->throw("Can't create a status object with no jobid");
  $status  || $self->throw("Can't create a status object for job ".
			   $jobid." with no status string");
  $created || $self->throw("Can't create a status object for job ".
			   $jobid." with no created time");

  $self->jobid             ($jobid);
  $self->status            ($status);
  $self->created           ($created);

  return $self;
}


=head2 jobid

  Title   : jobid
  Usage   : $self->jobid
  Function: Get/set method for the jobid
  Returns : int
  Args    : int

=cut

sub jobid {
    my ($self,$arg) = @_;

    if (defined($arg)) {
	$self->{_jobid} = $arg;
    }

    return $self->{_jobid};
}

=head2 status

  Title   : status
  Usage   : $self->status
  Function: Get/set method for the status string
  Returns : string
  Args    : string

=cut

sub status {
    my ($self,$arg) = @_;

    if (defined($arg)) {
	$self->{_status} = $arg;
    }

    return $self->{_status};
}

=head2 created

  Title   : created
  Usage   : $self->created
  Function: Get/set method for the created time
  Returns : int
  Args    : int

=cut

sub created {
    my ($self,$arg) = @_;

    if (defined($arg)) {
	$self->{_created} = $arg;
    }

    return $self->{_created};
}

1;
