#
#
# POD documentation - main docs before the code

=pod

=head1 NAME

Bio::EnsEMBL::Pipeline::RunnableDB::Syntenylite

=head1 SYNOPSIS

  my $synteny= Bio::EnsEMBL::Pipeline::Runnable::Syntenylite->new ( -prot_hits => $prot_hits );
  $synteny->run();
  $synteny->output();

=head1 DESCRIPTION

This modules takes a protein id eg. IPI:0000122 and the contig on which it was found on and finds a corrsponding chr if it exists.. it returns a AlignBlock(?) with the scaf, chr and the same protein as the corr. link.

=head1 CONTACT

Describe contact details here

Tania Oh (tania@fugu-sg.org) 

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Pipeline::Runnable::Syntenylite;


use strict;
use Bio::EnsEMBL::Pipeline::RunnableI;
use Bio::EnsEMBL::Compara::Synteny_Cluster;
use Data::Dumper;

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI);

=head2 new

    Title   :   new
    Usage   :   $self->new( -prot_match => $prot_match,
                            -contigid =>$contigid,
                          );
                           


    Function:   creates a Bio::EnsEMBL::Pipeline::Runnable::Syntenylite object
    Returns :   A Bio::EnsEMBL::Pipeline::RunnableDB::Syntenylite object
    Args    :   -prot_match: Bio::EnsEMBL::Compara::GenomicAlign containing align_id(protein id from the genomic align
                             block table, and the adaptor); 
                -parameters : ref. to a hash of parameters

=cut

sub new {
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);
 

  $self->{'_prot_hits'} = [];
  $self->{'_syntenyhits'} =undef;


  # Read the input parameters and set them in the object

  my ( $prot_hits, $adaptor) = $self->_rearrange ([qw(PROT_HITS ADAPTOR)], @args);
  if ($prot_hits) {
  my @prots = @$prot_hits;
    $self->protein_matches(@prots);
  } else {
    $self->throw("No protein_hits");
  }
  
  if (defined $adaptor){
      $self->adaptor($adaptor);
  }

 
  return $self;
}


=head2 run

    Title   :   run
    Usage   :   $self->run;
    Function:   Computes the clusters and returns them 
    Returns :   none
    Args    :   none

=cut

sub run {
  my ($self) = @_;
  
    #creates a cluster with the proteins
    #testing:

    
    my @prots = $self->protein_matches; 
    my $prot = \@prots;
    my $adaptor = $self->adaptor;
    my $cluster = Bio::EnsEMBL::Compara::Synteny_Cluster->new(-syn_hits => $prot, 
                                                              -syn_desc => "some kind of protein",
                                                              );



    $self->syntenyhits($cluster);
}


=head2 output

    Title   :   output
    Usage   :   $self->output();
    Function:   Runs Bio::EnsEMBL::Pipeline::Runnable::Synteny->output()
    Returns :   An array of Bio::EnsEMBL::Repeat objects (FeaturePairs)
    Args    :   none

=cut

sub output {
    my ($self) = @_;

    return $self->syntenyhits;
}




=head2 contigid 

    Title   :   contigid 
    Usage   :   $self->contigid('121.1');
    Function:   Get/set for the program to know which contig we are talking about 
    Returns :   String
    Args    :   String

=cut

sub contigid{
    my($self,$arg) = @_;

    if (defined($arg)) {
      $self->{'_contigid'} = $arg;
    }

    return $self->{'_contigid'};
}



=head2 syntenyhits 

    Title   :  syntenyhits 
    Usage   :   @out = $obj->alignedblocks
    Function:   Returns the output alignedblocks after they\'ve
                been converted from featurepairs 
    Returns :   @Bio::EnsEMBL::Compara::AlignBlocks
    Args    :   Bio::EnsEMBL::AlignedBlocks?

=cut



sub syntenyhits{
    my ($self, @sh) = @_;
   
  if (@sh){
    foreach my $sh (@sh)
    {
        $self->throw("$sh isn't a Bio::EnsEMBL::Compara::Synteny_Cluster")
                unless $sh->isa("Bio::EnsEMBL::Compara::Synteny_Cluster");
        push (@{$self->{'_syntenyhits'}}, $sh);
    }
   }
    return @{$self->{'_syntenyhits'}};
}

=head2 protein_matches


    Title   :  protein_matches
    Usage   :  $self->protein_matches(@protein_id);
    Function:
    Returns :  A simple  array of protein ids
    Args    :  none

=cut

sub protein_matches{
    my ($self, @proteins) = @_;

    if (@proteins)
    {
       foreach my $prot(@proteins)
          {
      		    $prot->isa("Bio::EnsEMBL::Compara::Synteny_Hit") or
              $self->throw("there are no proteins available-in protein_matches, Synteny!!");
          }
       push (@{$self->{'_protein_matches'}}, @proteins);
     }
     return @{$self->{'_protein_matches'}};
}

=head2 adaptor

 Title   : adaptor
 Usage   : $obj->adaptor($newval)
 Function: Getset for adaptor value
 Returns : value of adaptor
 Args    : newvalue (optional)


=cut

sub  adaptor{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'adaptor'} = $value;
    }
    return $obj->{'adaptor'};

}





1;





