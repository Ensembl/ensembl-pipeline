# Mar 9, 2006 11:11:17 AM
#
# Created by Mustapha Larbaoui <ml6@sanger.ac.uk>

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Pipeline::DBSQL::Finished::JobAdaptor

=head1 SYNOPSIS

=head1 DESCRIPTION

Module which inherits from JobAdaptor and return a Finished Job object

=head1 FEEDBACK

=head1 AUTHOR - Mustapha Larbaoui

Mustapha Larbaoui E<lt>ml6@sanger.ac.ukE<gt>

=head1 CONTACT

Post general queries to B<anacode@sanger.ac.uk>

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _


=cut


# Let the code begin...

package Bio::EnsEMBL::Pipeline::DBSQL::Finished::JobAdaptor;

use vars qw(@ISA);
use strict;
use Bio::EnsEMBL::Pipeline::Finished::Job;
use Bio::EnsEMBL::Pipeline::DBSQL::JobAdaptor;

@ISA = qw(Bio::EnsEMBL::Pipeline::DBSQL::JobAdaptor);

=head2 _objFromHashref

  Title   : _objFromHashref
  Usage   : my $job = $self->objFromHashref( $queryResult )
  Function: Creates a Finished Job object from given hash reference.
            The hash contains column names and content of the column.
  Returns : the object or undef if that wasnt possible
  Args    : a hash reference

=cut

sub _objFromHashref {
  # create the appropriate job object

  my $self = shift;
  my $hashref = shift;
  my $job;
  my $analysis;
 
  
  $analysis =
    $self->db->get_AnalysisAdaptor->
      fetch_by_dbID( $hashref->{analysis_id} );
  $job = Bio::EnsEMBL::Pipeline::Finished::Job->new
    (
     '-dbobj'     => $self->db,
     '-adaptor'   => $self,
     '-id'        => $hashref->{'job_id'},
     '-submission_id' => $hashref->{'submission_id'},
     '-input_id'  => $hashref->{'input_id'},
     '-stdout'    => $hashref->{'stdout_file'},
     '-stderr'    => $hashref->{'stderr_file'},
     '-analysis'  => $analysis,
     '-retry_count' => $hashref->{'retry_count'},
     '-exec_host' => $hashref->{'exec_host'},
     '-temp_dir' => $hashref->{'temp_dir'},
  );
  return $job;
}

sub fetch_by_Status_not_like {
    my ($self, $status, $start, $end) = @_;

    throw("Require status for fetch_by_Status")
                            unless ($status);
    

    my $query = q{
	SELECT   j.job_id, j.input_id, j.analysis_id, j.submission_id,
	         j.stdout_file, j.stderr_file, j.retry_count, j.temp_dir, 
           j.exec_host
	FROM     job j, job_status js
        WHERE    j.job_id = js.job_id
        AND      js.status != ?
	AND      js.is_current = 'y'
        ORDER BY time desc
    };
    
    $query .= " LIMIT $start, $end" if ($start && $end);

    my $sth = $self->prepare($query);
    my $res = $sth->execute($status);

    my @jobs;

    while (my $row = $sth->fetchrow_hashref)
    {
	    my $job = $self->_objFromHashref($row);
	    push(@jobs,$job);
    }
    return @jobs;
}

1;
