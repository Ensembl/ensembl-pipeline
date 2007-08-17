#
# Cared for by EnsEMBL  <ensembl-dev@ebi.ac.uk>
#
# written by Michele Clamp
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::RunnableDB::Gene_Builder

=head1 SYNOPSIS

    my $obj = Bio::EnsEMBL::Pipeline::RunnableDB::Gene_Builder->new(
								    -db        => $db,
								    -input_id  => $id,
								    );
    $obj->fetch_input
    $obj->run

    my @newfeatures = $obj->output;


=head1 DESCRIPTION

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::EnsEMBL::Pipeline::RunnableDB::Havana_Adder;

use vars qw(@ISA);
use strict;

# Object preamble
use Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils;
use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Pipeline::HavanaAdder;

use Bio::EnsEMBL::Pipeline::Config::GeneBuild::General     qw (
							       GB_INPUTID_REGEX
							      );
use Bio::EnsEMBL::Pipeline::Config::GeneBuild::GeneBuilder qw (
							       GB_VCONTIG
							      );
use Bio::EnsEMBL::Pipeline::Config::HavanaAdder            qw (
                                                               GB_GENE_OUTPUT_BIOTYPE
                                                              );
use Bio::EnsEMBL::Pipeline::Config::GeneBuild::Databases   qw (
                                                               GB_FINALDBHOST
                                                               GB_FINALDBNAME
                                                               GB_FINALDBUSER
                                                               GB_FINALDBPASS
                                                               GB_FINALDBPORT
                                                               PSEUDO_DBHOST
                                                               PSEUDO_DBNAME
                                                               PSEUDO_DBUSER
                                                               PSEUDO_DBPASS
                                                               PSEUDO_DBPORT
                                                               GB_HAVANA_DBHOST
                                                               GB_HAVANA_DBNAME
                                                               GB_HAVANA_DBUSER
                                                               GB_HAVANA_DBPASS
                                                               GB_HAVANA_DBPORT
                                                              );

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDB);

############################################################

=head2 new

    Usage   :   $self->new(-DBOBJ       => $db,
                           -INPUT_ID    => $id,
                           -SEQFETCHER  => $sf,
                           -ANALYSIS    => $analysis,
                           -VCONTIG     => 1,
                          );

                           
    Function:   creates a Bio::EnsEMBL::Pipeline::RunnableDB::Gene_Builder object
    Returns :   A Bio::EnsEMBL::Pipeline::RunnableDB::Gene_Builder object
    Args    :   -dbobj:      A Bio::EnsEMBL::DBSQL::DBAdaptor (required), 
                -input_id:   Contig input id (required), 
                -seqfetcher: A Sequence Fetcher Object,
                -analysis:   A Bio::EnsEMBL::Analysis (optional) 
                -vcontig:    determines whether it is running on virtual contigs
                             or RawContigs
                -extend:     determines the extension of the virtual contig
                             note: not implemented yet!
                -golden_path: determines the name of the golden path to use
=cut

sub new {
    my ($class,@args) = @_;
    my $self = $class->SUPER::new(@args);    
           
    my( $use_vcontig) = $self->_rearrange([qw(VCONTIG)], @args);
       
    if (! defined $use_vcontig) {
      $use_vcontig = $GB_VCONTIG;
    }  
    
    $self->use_vcontig($use_vcontig);

    return $self;
}

############################################################

sub input_id {
    my ($self,$arg) = @_;
    
    if (defined($arg)) {
	$self->{_input_id} = $arg;
    }
    
    return $self->{_input_id};
}

############################################################

=head2 write_output

    Title   :   write_output
    Usage   :   $self->write_output
    Function:   Writes output data to db
    Returns :   array of exons (with start and end)
    Args    :   none

=cut
    
    
sub write_output {
  my($self,@genes) = @_;
  
  # write genes out to a different database from the one we read genewise genes from.
  my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
                                              '-host'   => $GB_FINALDBHOST,
                                              '-user'   => $GB_FINALDBUSER,
                                              '-dbname' => $GB_FINALDBNAME,
                                              '-pass'   => $GB_FINALDBPASS,
                                              '-port'   => $GB_FINALDBPORT,
                                              '-dnadb'  => $self->db,
                                              );
  # sort out analysis
  
  my $analysis = $self->analysis;
  unless ($analysis){
    $self->throw("an analysis logic name must be defined in the command line");
  }
  
  my %contighash;
  my $gene_adaptor = $db->get_GeneAdaptor;
  
  # this now assummes that we are building on a single VC.
  my $genebuilders = $self->get_genebuilders;
    
  foreach my $contig ( keys %$genebuilders ){
    my $vc = $genebuilders->{$contig}->query;
      
    @genes = $genebuilders->{$contig}->final_genes;
    
    return unless ($#genes >= 0);
    my @newgenes;
    
    foreach my $gene (@genes) { 
      $gene->analysis($analysis);
      $gene->type($GB_GENE_OUTPUT_BIOTYPE);
      # poke the caches
      my %s_pfhash;
      foreach my $tran (@{$gene->get_all_Transcripts}) {
        $tran->stable_id(undef);
        my @tsf = @{$tran->get_all_supporting_features};
        
#this may be a miostake - will proftein features get transferred?
#          if (defined($tran->translation) && defined($s_pfa)) {
#              $s_pfhash{get_transcript_id($tran)} = $s_pfa->fetch_by_translation_id($tran->translation->dbID);
#              print "Got " . scalar(@{$s_pfhash{get_transcript_id($tran)}}) . " pfs for translation " . get_translation_id($tran->translation) . "\n";
#          }
        
        my @exons= @{$tran->get_all_Exons};
        my $tln = $tran->translation;
        $tln->{'stable_id'} = undef;
        
        foreach my $exon (@exons) {
          $exon->{'stable_id'} = undef;
        }
      }  
      # store
      eval {
        $gene_adaptor->store($gene);
        #print STDERR "wrote gene " . $gene->dbID . " to database ".
        #   $gene->adaptor->db->dbname."\n";
      }; 
      if( $@ ) {
        $self->warn("UNABLE TO WRITE GENE:\n$@");
      }
    }   
  }    
}

############################################################

=head2 fetch_input

    Function:   It fetches the slice or contig according to the input_id, 
                and it defines the database where the
                previous annotations are stored and create a Bio::EnsEMBL::Pipeline::GeneBuilder
                object for that genomic, input_id and db
    Returns :   nothing
    Args    :   none

=cut

sub fetch_input {
    my( $self) = @_;
    
    $self->throw("No input id") unless defined($self->input_id);
    
    $self->fetch_sequence();
    # database where the genebuild produced genes are
    my $ensembl_db = new Bio::EnsEMBL::DBSQL::DBAdaptor
      (
       '-host'   => $PSEUDO_DBHOST,
       '-user'   => $PSEUDO_DBUSER,
       '-dbname' => $PSEUDO_DBNAME,
       '-pass'   => $PSEUDO_DBPASS,
       '-port'   => $PSEUDO_DBPORT,
       '-dnadb'  => $self->db,
      );

    my $havana_db = new Bio::EnsEMBL::DBSQL::DBAdaptor
      (
       '-host'   => $GB_HAVANA_DBHOST,
       '-user'   => $GB_HAVANA_DBUSER,
       '-dbname' => $GB_HAVANA_DBNAME,
       '-pass'   => $GB_HAVANA_DBPASS,
       '-port'   => $GB_HAVANA_DBPORT,
       '-dnadb'  => $self->db,
      );
    
    #print STDERR "reading genewise and combined genes from $GB_COMB_DBNAME : $GB_COMB_DBHOST\n";
    
    my $genebuilder = new Bio::EnsEMBL::Pipeline::HavanaAdder
      (
       '-slice'   => $self->query,
       '-input_id' => $self->input_id,
      );
    $genebuilder->ensembl_db($ensembl_db);
    $genebuilder->havana_db($havana_db);
    
    # store the object and the piece of genomic where it will run
    $self->addgenebuilder($genebuilder,$self->query);
    
}

############################################################

sub use_vcontig {
    my ($self,$arg) = @_;
    
    if (defined($arg)) {
	$self->{_vcontig} = $arg;
    }

    return $self->{_vcontig};
}

############################################################

sub addgenebuilder {
    my ($self,$arg,$contig) = @_;
    
    if (defined($arg) && defined($contig)) {
	$self->{_genebuilder}{$contig->id} = $arg;
    } 
    else {
	$self->throw("Wrong number of inputs [$arg,$contig]\n");
    }
}

############################################################

sub get_genebuilders {
    my ($self) = @_;
    
    return $self->{_genebuilder};
}

############################################################
	
sub run {
    my ($self) = @_;
    
    # get a hash, with keys = contig/slice and value = genebuilder object
    my $genebuilders = $self->get_genebuilders;
    
    my @genes;
    foreach my $contig (keys %{ $genebuilders } ) {
      my $query = $genebuilders->{$contig}->query;
      
      #print(STDERR "GeneBuilding for $contig\n");
      
      $genebuilders->{$contig}->build_Genes;
      
      @genes = $genebuilders->{$contig}->final_genes;
    }
    
    $self->output( @genes );
}

############################################################

# override the evil RunnableDB output method:

sub output{
    my ($self, @genes ) = @_;
    unless ( $self->{_output} ){
	$self->{_output} = [];
    }
    if (@genes){
	push( @{$self->{_output}}, @genes );
    }
    return @{$self->{_output}};
}

############################################################



1;
