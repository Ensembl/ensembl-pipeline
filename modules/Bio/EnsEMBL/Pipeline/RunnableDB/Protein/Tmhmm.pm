# Author: Marc Sohrmann (ms2@sanger.ac.uk)
# Copyright (c) Marc Sohrmann, 2001
# You may distribute this code under the same terms as perl itself
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

  Bio::EnsEMBL::Pipeline::RunnableDB::Protein::Tmhmm

=head1 SYNOPSIS

  my $tmhmm = Bio::EnsEMBL::Pipeline::RunnableDB::Protein::Tmhmm->new ( -db      => $db,
	  	                                                        -input_id   => $input_id,
                                                                        -analysis   => $analysis,
                                                                      );
  $tmhmm->fetch_input;  # gets sequence from DB
  $tmhmm->run;
  $tmhmm->output;
  $tmhmm->write_output; # writes features to to DB

 NB: The input_id can either be a peptide id or the location for a protein file. 

=head1 DESCRIPTION

  This object wraps Bio::EnsEMBL::Pipeline::Runnable::Protein::Tmhmm
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

package Bio::EnsEMBL::Pipeline::RunnableDB::Protein::Tmhmm;

use strict;
use vars qw(@ISA);


use Bio::EnsEMBL::Pipeline::Runnable::Protein::Tmhmm;
use Bio::EnsEMBL::Pipeline::RunnableDB::Protein::Protein_Annotation;


@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDB::Protein::Protein_Annotation);




# runnable method

=head2 runnable

 Title    :  runnable
 Usage    :  $self->runnable($arg)
 Function :  sets a runnable for this RunnableDB
 Example  :
 Returns  :  Bio::EnsEMBL::Pipeline::RunnableI
 Args     :  Bio::EnsEMBL::Pipeline::RunnableI
 Throws   :

=cut

sub runnable {
  my ($self) = @_;
    
    if (!defined($self->{'_runnable'})) {
	
	my $run = Bio::EnsEMBL::Pipeline::Runnable::Protein::Tmhmm->new(-query     => $self->query,
									-analysis  => $self->analysis	);
	
	
	$self->{'_runnable'} = $run;
    }
  
  return $self->{'_runnable'};
}



1;
