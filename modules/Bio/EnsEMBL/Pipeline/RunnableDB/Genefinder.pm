#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::RunnableDB::Genefinder

=head1 SYNOPSIS

  my $db      = Bio::EnsEMBL::DBLoader->new($locator);
  my $genscan = Bio::EnsEMBL::Pipeline::RunnableDB::Genefinder->new ( 
                                                    -dbobj      => $db,
			                            -input_id   => $input_id
                                                    -analysis   => $analysis );
  $genscan->fetch_input();
  $genscan->run();
  $genscan->output();
  $genscan->write_output(); #writes to DB

=head1 DESCRIPTION

This object wraps Bio::EnsEMBL::Pipeline::Runnable::Genefinder to add
functionality to read and write to databases.
The appropriate Bio::EnsEMBL::Analysis object must be passed for
extraction of appropriate parameters. A Bio::EnsEMBL::Pipeline::DBSQL::Obj is
required for databse access.

=head1 CONTACT

Post general queries to B<ensembl-dev@ebi.ac.uk>

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Pipeline::RunnableDB::Genefinder;

use strict;

use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Pipeline::Runnable::Genefinder;
use Bio::EnsEMBL::Pipeline::Config::General;
use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDB);

=head2 fetch_input

=cut

sub fetch_input {
    my( $self) = @_;
    

    $self->throw("No input id") unless defined($self->input_id);
    
    my $contigid  = $self->input_id;
    my $contig    = $self->db->get_RawContigAdaptor->fetch_by_name($contigid);

    $self->query($contig->get_repeatmasked_seq($PIPELINE_REPEAT_MASKING));

    my $runnable = new Bio::EnsEMBL::Pipeline::Runnable::Genefinder(
	      -query   => $self->query,
              -genscan => $self->analysis->program_file,
              -matrix  => $self->analysis->db_file,
	      -args    => $self->arguments
    );

    $self->runnable($runnable);

    return 1;
}


1;
