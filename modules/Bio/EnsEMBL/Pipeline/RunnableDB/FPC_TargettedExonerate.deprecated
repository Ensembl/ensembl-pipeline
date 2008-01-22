
# Ensembl module for Bio::EnsEMBL::Pipeline::RunnableDB::FPC_TargettedExonerate.pm
#
# Cared for by EnsEMBL  <ensembl-dev@ebi.ac.uk>
#
# Copyright GRL and EBI
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Pipeline::RunnableDB::FPC_TargettedExonerate

=head1 SYNOPSIS

=head1 DESCRIPTION

Runs all the TargettedExonerate jobs needed for a chunk of genomic sequence

=head1 CONTACT

Ensembl - ensembl-dev@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::Pipeline::RunnableDB::FPC_TargettedExonerate;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Pipeline::Config::GeneBuild::Databases qw (
                                                             GB_GW_DBNAME
                                                             GB_GW_DBHOST
                                                             GB_GW_DBUSER
                                                             GB_GW_DBPASS
                                                             GB_GW_DBPORT
                                                            );

use Bio::EnsEMBL::Pipeline::Config::GeneBuild::Sequences
  qw (
      GB_PROTEIN_INDEX
      GB_PROTEIN_SEQFETCHER
     );

use Bio::EnsEMBL::Pipeline::Config::GeneBuild::General
  qw (
      GB_INPUTID_REGEX
     );

use Bio::EnsEMBL::Pipeline::Config::GeneBuild::Pmatch
  qw (
      GB_FINAL_PMATCH_LOGICNAME
     );

use Bio::EnsEMBL::Pipeline::Config::GeneBuild::TargettedExonerate qw
  (
   GB_TARGETTEDXRATE_MASKING
   GB_TARGETTEDXRATE_SOFTMASK
  );

use Bio::EnsEMBL::Pipeline::RunnableDB::TargettedExonerate;
use Bio::EnsEMBL::Pipeline::RunnableDB::FPC_TargettedBase;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDB::FPC_TargettedBase);


=head2 make_targetted_runnables

 Title   : make_targetted_runnables
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub make_targetted_runnables {
  my ($self) = @_;

  # set up seqfetchers
  my $protein_fetcher = $self->make_seqfetcher($GB_PROTEIN_INDEX, $GB_PROTEIN_SEQFETCHER);

  # we need to find all the proteins that pmatch into this region
  # take a note of those that fall across the ends of the vc? and do what, precisely?
  # extend the VC? that will completely screw up the final genebuild. Hmmm.
  # do it, track it & see how many are affected.
  #input_id cb25.fpc4118.1-298757 has invalid format - expecting chr_name.start-end
  
  my $pmfa = $self->db->get_PmatchFeatureAdaptor();

  $self->fetch_sequence($GB_TARGETTEDXRATE_MASKING);
  my @features = $pmfa->get_PmatchFeatures_by_chr_start_end
    ($self->query->seq_region_name, $self->query->start, 
     $self->query->end, $GB_FINAL_PMATCH_LOGICNAME );
 
  my $targetted_db = new Bio::EnsEMBL::DBSQL::DBAdaptor
    (
     '-host'   => $GB_GW_DBHOST,
     '-user'   => $GB_GW_DBUSER,
     '-pass'   => $GB_GW_DBPASS,
     '-dbname' => $GB_GW_DBNAME,
     '-port' => $GB_GW_DBPORT,
    );
  
  
  $targetted_db->dnadb($self->db);
  $self->output_db($targetted_db);
  my %kill_list = %{$self->populate_kill_list};

  foreach my $feat(@features){
    #get the protein_id without a version
    my $feat_protein_id = $feat->protein_id;
    if (!defined $feat_protein_id) {
      throw("No feature protein_id in FPC_TargettedExonerate");
    } else {
      print "Feature has Protein_id $feat_protein_id\n";
    }
    $feat_protein_id =~ s/\.\d+//;
    #reject any proteins that are in the kill list
    if(defined $kill_list{$feat_protein_id}){
      #print STDERR "skipping " . $feat->protein_id . "\n";
      next;
    }
                                   
    my ($start, $end);
    
    $start = $feat->start;
    $end = $feat->end;

    my $input = ($self->query->coord_system->name.":".
                 $self->query->coord_system->version.":".
                 $self->query->seq_region_name.":".
                 $start.":".$end.":".
                 $self->query->strand.",".
                 $feat->protein_id);

    #print STDERR "\n\n***TGE input: $input***\n";

    my $tgr = new Bio::EnsEMBL::Pipeline::RunnableDB::TargettedExonerate
      (
       -db => $self->db,
       -input_id => $input,
       -seqfetcher => $protein_fetcher,
       -analysis => $self->analysis,
       -output_db => $self->output_db,
      );
    $self->targetted_runnable($tgr);
  }

}

1;
