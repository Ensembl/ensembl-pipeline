#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code
# author: mongin@ebi.ac.uk

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::RunnableDB::Snap

=head1 SYNOPSIS

  my $db      = Bio::EnsEMBL::DBLoader->new($locator);
  my $Snap = Bio::EnsEMBL::Pipeline::RunnableDB::Snap->new ( 
                                                    -dbobj      => $db,
			                            -input_id   => $input_id
                                                    -analysis   => $analysis );
  $snap->fetch_input();
  $snap->run();
  $snap->output();
  $snap->write_output(); #writes to DB

=head1 DESCRIPTION

This object wraps Bio::EnsEMBL::Pipeline::Runnable::Snap to add
functionality to read and write to databases.
The appropriate Bio::EnsEMBL::Analysis object must be passed for
extraction of appropriate parameters. A Bio::EnsEMBL::Pipeline::DBSQL::Obj is
required for databse access.

snap is a gene predictor written by Ian Korf (ik1@sanger.ac.uk) part the Zoe software library.

=head1 CONTACT

For general queries please contact <ensembl-dev@ebi.ac.uk>

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Pipeline::RunnableDB::Snap;

use strict;

use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Pipeline::Runnable::Snap;
use Bio::EnsEMBL::Pipeline::Config::General;
use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDB);

=head2 fetch_input

=cut

sub fetch_input {
    my( $self) = @_;
    

    $self->throw("No input id") unless defined($self->input_id);
   
    $self->fetch_sequence($SNAP_MASKING);

    
    my $runnable = new Bio::EnsEMBL::Pipeline::Runnable::Snap(
	      -query   => $self->query,
              -snap => $self->analysis->program_file,
              -hmmfile  => $self->analysis->db_file,
	      -args    => $self->arguments
    );

    $self->runnable($runnable);

    return 1;
}


1;
