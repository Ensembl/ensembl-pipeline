# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::RunnableDB::Pseudotest_v2_DB.pm

=head1 SYNOPSIS



=head1 DESCRIPTION

This object wraps Bio::EnsEMBL::Pipeline::Runnable::Pseudotest_v2.pm to add
functionality to read from databases (so far).

=head1 CONTACT

Describe contact details here

=head1 APPENDIX



=cut

package Bio::EnsEMBL::Pipeline::RunnableDB::Pseudotest_v2_DB;

use strict;

use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Pipeline::Runnable::Pseudotest_v2;
use Bio::EnsEMBL::Pipeline::Config::PseudoGenes::Pseudotest_v2_config;

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDB);

=head2 fetch_input

    Title   :   fetch_input
    Usage   :   $self->fetch_input
    Function:   Fetches input data for Pseudotest_v2.pm from the database
    Returns :   none
    Args    :   none

=cut



sub fetch_input {
    my( $self) = @_;
    
    $self->throw("No input id") unless defined($self->input_id);

    $self->fetch_sequence;

    my %parameters = $self->parameter_hash;
    my %repeats;
    $parameters{'-query'} = $self->query;



    my $runname = "Bio::EnsEMBL::Pipeline::Runnable::Pseudotest_v2";

  my $rep_db = new Bio::EnsEMBL::DBSQL::DBAdaptor
    (
       '-host'   => $REPEAT_DBHOST,
       '-user'   => $REPEAT_DBUSER,
       '-dbname' => $REPEAT_DBNAME,
       '-pass'   => $REPEAT_DBPASS,
       '-port'   => $REPEAT_DBPORT,
       '-dnadb'  => $self->db,
      );

    my $genes = $self->query->get_all_Genes;
    foreach my $gene(@{$genes}){
      my $rsa = $rep_db->get_SliceAdaptor;
      my $rep_gene_slice = $rsa->fetch_by_region(
						    'toplevel',
						    $self->query->chr_name,
						    $gene->start,
						    $gene->end,
						  );

      $repeats{$gene} = $rep_gene_slice->get_all_RepeatFeatures;
      $gene->transfer($rep_gene_slice);
    }
    my $runnable = $runname->new
      ( 
       '-max_intron_length' => $MAX_INTRON_LENGTH,
       '-max_intron_coverage' => $MAX_INTRON_COVERAGE,
       '-max_exon_coverage' => $MAX_EXON_COVERAGE,
       '-genes' => $genes,
       '-repeat_features' => \%repeats
      );

    $self->runnable($runnable);

    return 1;

}

#sub write_output {
#    my($self,@genes) = @_;
#    
#    # write genes out to a different database from the one we read genewise genes from.
#    my $dbname = $PSEUDO_DBNAME;
#   my $dbhost = $PSEUDO_DBHOST;
#   my $dbuser = $PSEUDO_DBUSER;
#   my $dbpass = $PSEUDO_DBPASS;
#   my $dbport = $PSEUDO_DBPORT;
#
#    my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
#						'-host'   => $dbhost,
#						'-user'   => $dbuser,
#						'-dbname' => $dbname,
#						'-pass'   => $dbpass,
#						'-port'   => $dbport,
#						'-dnadb'  => $self->db,
#					       );
#    # sort out analysis
#    my $analysis = $self->analysis;
#    unless ($analysis){
#	$self->throw("an analysis logic name must be defined in the command line");
#    }
#
#    my $gene_adaptor = $db->get_GeneAdaptor;
#	
#	foreach my $gene (@genes) { 
#	    $gene->analysis($analysis);
#	    $gene->type($genetype);
#
#	    # store
#	    eval {
#        $gene_adaptor->store($gene);
#        print STDERR "wrote gene " . $gene->dbID . " to database ".
#          $gene->adaptor->db->dbname."\n";
#	    }; 
#	    if( $@ ) {
#        $self->warn("UNABLE TO WRITE GENE:\n$@");
#      }
#	  }
#
#  }
#
#}


1;
