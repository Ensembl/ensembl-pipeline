# Author: Marc Sohrmann (ms2@sanger.ac.uk)
# Copyright (c) Marc Sohrmann, 2001
# You may distribute this code under the same terms as perl itself
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

  Bio::EnsEMBL::Pipeline::RunnableDB::Protein::Signalp

=head1 SYNOPSIS

  my $signalp = Bio::EnsEMBL::Pipeline::RunnableDB::Protein::Signalp->new ( -db      => $db,
    	  	                                                            -input_id   => $input_id,
                                                                            -analysis   => $analysis,
                                                                          );
  $signalp->fetch_input;  # gets sequence from DB
  $signalp->run;
  $signalp->output;
  $signalp->write_output; # writes features to to DB

 NB: The input_id can either be a peptide id or the location for a protein file. 

=head1 DESCRIPTION

  This object wraps Bio::EnsEMBL::Pipeline::Runnable::Protein::Signalp
  to add functionality to read and write to databases.
  A Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor is required for database access (db).
  The query sequence is provided through the input_id.
  The appropriate Bio::EnsEMBL::Analysis object
  must be passed for extraction of parameters.

=head1 CONTACT

  Marc Sohrmann: ms2@sanger.ac.uk

=head1 APPENDIX

  The rest of the documentation details each of the object methods. 
  Internal methods are usually preceded with a _.

=cut

package Bio::EnsEMBL::Pipeline::RunnableDB::Protein::Signalp;

use strict;
use vars qw(@ISA);

use Bio::EnsEMBL::Pipeline::RunnableDB::Protein::Protein_Annotation;
use Bio::EnsEMBL::Pipeline::Runnable::Protein::Signalp;
use Bio::EnsEMBL::DBSQL::ProteinFeatureAdaptor;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDB::Protein::Protein_Annotation);




sub runnable {
    my ($self) = @_;
    
    if (!($self->{'_runnable'})) {
	print STDERR "QUERY: ".$self->query."\n";
	my $run = Bio::EnsEMBL::Pipeline::Runnable::Protein::Signalp->new(-query     => $self->query,
									-analysis  => $self->analysis	);
 
           
      $self->{'_runnable'} = $run;
    }
    
    return $self->{'_runnable'};
}


1;
