#
# BioPerl module for Prints.pm
#
# Cared for by Emmanuel Mongin <mongin@ebi.ac.uk>
#
# Copyright Emmanuel Mongin
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Prints.pm - DESCRIPTION of Object

=head1 SYNOPSIS

 $self->new(-DB       => $db
                           -INPUT_ID    => $id
                           -ANALYSIS    => $analysis);

Where the analysis id can be either a translation internal id or the location of a file. That\'s feature is used only for the protein annotation. Some of the protein analysis are extremely fast to run (eg: seg) in that cases a full protein dataset will be given to the RunnableDB (running a protein at a time would be to expensive).

=head1 DESCRIPTION

 This object wraps Bio::EnsEMBL::Pipeline::Runnable::Protein::Seg
  to add functionality to read and write to databases.
  A Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor is required for database access (db).
  The query sequence is provided through the input_id.
  The appropriate Bio::EnsEMBL::Analysis object
  must be passed for extraction of parameters.

=head1 CONTACT

mongin@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::Pipeline::RunnableDB::Protein::Prints;
use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Pipeline::Runnable::Protein::Prints;
use Bio::EnsEMBL::DBSQL::ProteinAdaptor;
use Bio::EnsEMBL::DBSQL::ProteinFeatureAdaptor;
use Bio::EnsEMBL::Pipeline::Config::Protein_Annotation::General;
use Bio::EnsEMBL::Pipeline::RunnableDB::Protein::Protein_Annotation;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDB::Protein::Protein_Annotation);



#get/set for runnable and args
sub runnable {
  my ($self) = @_;
    
    if (!($self->{'_runnable'})) {
	
	my $run = Bio::EnsEMBL::Pipeline::Runnable::Protein::Prints->new(-query     => $self->query,
									  -analysis  => $self->analysis	);
	
	
	$self->{'_runnable'} = $run;
    }
    
    return $self->{'_runnable'};
}




1;










