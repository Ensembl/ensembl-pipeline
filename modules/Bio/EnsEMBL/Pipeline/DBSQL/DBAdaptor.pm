#
# Written by Simon Potter <scp@sanger.ac.uk>
# Based on Michele Clamp's Bio::EnsEMBL::Pipeline::DBSQL::Obj
#
# Copyright GRL/EBI
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code


=pod

=head1 NAME

Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor -
adapter class for EnsEMBL Pipeline DB

=head1 SYNOPSIS

 my $pipedb = Bio::EnsEMBL::Pipeline::DBAdaptor->new
    (-dbname => 'pipeline_db',
     -host   => 'ecs2a',
     -user   => 'ensadmin',
     -pass   => 'secret');

 my $job_adaptor = $pipedb->get_JobAdaptor();

=head1 DESCRIPTION

Maintains a persistant database connection and provides adaptors for retrieval
and storage of database objects relating to a running pipeline system

=head1 CONTACT

Post general queries to B<ensembl-dev@ebi.ac.uk>

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


package Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor;


use vars qw(@ISA);
use strict;
use Bio::EnsEMBL::DBSQL::DBConnection;

@ISA = qw(Bio::EnsEMBL::DBSQL::DBConnection);


sub get_JobAdaptor {
  my ($self) = @_;

  return $self->_get_adaptor('Bio::EnsEMBL::Pipeline::DBSQL::JobAdaptor');
}


sub pipeline_lock {
    my ($self, $string) = @_;

    my $sth;

    if ($string) {
	$sth = $self->prepare(qq{
	    INSERT into meta (meta_key, meta_value)
	    VALUES ('pipeline.lock', ?)
	});
	$sth->execute($string);
    }
    else {
	$sth = $self->prepare(qq{
	    SELECT meta_value
	    FROM   meta
	    WHERE  meta_key = 'pipeline.lock'
	});

        $sth->execute;
        my $row = $sth->fetchrow_arrayref;

        if ($row) {
	    return $row->[0];
        }
        else {
	    return undef;
        }
    }
    $sth->finish;
}


sub pipeline_unlock {
    my ($self) = @_;

    my $sth;

    $sth = $self->prepare(qq{
	DELETE
	FROM   meta
	WHERE  meta_key = 'pipeline.lock'
    });

    $sth->execute;
    $sth->finish;
}

1;
