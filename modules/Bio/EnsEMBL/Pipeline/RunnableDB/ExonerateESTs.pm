#
#
# Cared for by EnsEMBL  <ensembl-dev@ebi.ac.uk>
#
# Copyright GRL & EBI
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::RunnableDB::FPC_BlastMiniEst2Genome

=head1 SYNOPSIS

    my $obj = Bio::EnsEMBL::Pipeline::RunnableDB::FPC_BlastMiniEst2Genome->new(
					     -dbobj     => $db,
					     -input_id  => $id,
					     -estfile   => $estfile
                                             );
    $obj->fetch_input
    $obj->run

    my @genes = $obj->output;


=head1 DESCRIPTION
Just runs Exonerate over a  chunk of dbEST and spits the output to file to be assessed later.
File1 => all FeaturePairs produced by Exonerate

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::EnsEMBL::Pipeline::RunnableDB::ExonerateESTs;

use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::RootI;
use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Pipeline::Runnable::ExonerateESTs;

@ISA = qw( Bio::EnsEMBL::Pipeline::RunnableDB );


=head2 new

    Title   :   new
    Usage   :   $self->new(-DBOBJ       => $db
                           -INPUT_ID    => $id
                           -ANALYSIS    => $analysis);
                           
    Function:   creates a 
                Bio::EnsEMBL::Pipeline::RunnableDB::ExonerateESTs
                object
    Returns :   A Bio::EnsEMBL::Pipeline::RunnableDB::ExonerateESTs
                object
    Args    :   -dbobj:      A Bio::EnsEMBL::DB::Obj (required), 
                -input_id:   Contig input id (required), 
                -analysis:   A Bio::EnsEMBL::Pipeline::Analysis (optional) 
=cut

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);
           
    # dbobj, input_id, seqfetcher, and analysis objects are all set in
    # in superclass constructor (RunnableDB.pm)

    my( $estfile, $allresults ) = $self->_rearrange([qw(ESTFILE ALLRESULTS)], @args);

    $self->throw("No est file specified") unless defined($estfile);
    $self->estfile($estfile);

    $self->throw("No allresults file specified") unless defined($allresults);
    $self->allresults($allresults);

    return $self;
}

=head2 write_output

    Title   :   write_output
    Usage   :   $self->write_output
    Function:   Writes output to file
    Returns :   
    Args    :   none

=cut

sub write_output {

  my ($self) = @_;
  foreach my $feat(@{$self->output}) {

    print STDERR "$feat\n";
  }
}

=head2 fetch_input

    Title   :   fetch_input
    Usage   :   $self->fetch_input
    Function:   Fetches input data for ExonerateESTs and makes runnable
    Returns :   nothing
    Args    :   none

=cut

sub fetch_input {
  my ($self) = @_;
  
#  print STDERR "Fetching input \n";
  $self->throw("No input id") unless defined($self->input_id);
  
  my $chrid  = $self->input_id;
      $chrid =~ s/\.(.*)-(.*)//;

  my $chrstart = $1;
  my $chrend   = $2;
  
#  print STDERR "Chromosome id = $chrid , range $chrstart $chrend\n";
  
  $self->dbobj->static_golden_path_type('UCSC');
  
  my $stadaptor = $self->dbobj->get_StaticGoldenPathAdaptor();
  my $contig    = $stadaptor->fetch_VirtualContig_by_chr_start_end($chrid,$chrstart,$chrend);
  
  $contig->_chr_name($chrid);
  
  my $genseq    = $contig->get_repeatmasked_seq;
  
  my $runnable  = new Bio::EnsEMBL::Pipeline::Runnable::ExonerateESTs('-genomic'     => $genseq, 
								      '-ests'        => $self->estfile,
								      '-resfile'     => $self->allresults);
  
  $self->runnable($runnable);
  # at present, we'll only ever have one ...
  $self->vc($contig);
}

=head2 run

    Title   :   run
    Usage   :   $self->run
    Function:   Calls run method of each runnable, & converts output into remapped genes
    Returns :   Nothing
    Args    :   None

=cut

sub run {
  my ($self) = @_;

  $self->throw("Can't run - no runnable objects") unless defined($self->{_runnables});
  
  foreach my $runnable($self->runnable) {
    $runnable->run;

    my $out = $runnable->output;
    $self->output($out);
  }

  # something about filtering here?

}

=head2 output

    Title   :   output
    Usage   :   $self->output
    Function:   Returns output from this RunnableDB
    Returns :   Array of Bio::EnsEMBL::Gene
    Args    :   None

=cut

sub output {
   my ($self,$feat) = @_;

   if (!defined($self->{'_output'})) {
     $self->{'_output'} = [];
   }
    
   if(defined $feat){
     push(@{$self->{'_output'}},@{$feat});
   }

   return $self->{'_output'};# ref to an array
}

=head2 vc

 Title   : vc
 Usage   : $obj->vc($newval)
 Function: 
 Returns : value of vc
 Args    : newvalue (optional)

=head2 estfile

 Title   : estfile
 Usage   : $obj->estfile($newval)
 Function: 
 Returns : value of estfile
 Args    : newvalue (optional)


=cut

sub estfile {
   my ($self, $value) = @_;
   if( defined $value ) {
      $self->{'_estfile'} = $value;
    }
    return $self->{'_estfile'};

}

=head2 allresults

 Title   : allresults
 Usage   : $obj->allresults($newval)
 Function: get/set filename to which all exonerate results should be written
 Returns : 
 Args    : 


=cut

sub allresults {
   my ($self, $value) = @_;
   if( defined $value ) {
      $self->{'_allresults'} = $value;
    }
    return $self->{'_allresults'};

}

=head2 filtered

 Title   : filtered
 Usage   : $obj->filtered($newval)
 Function: gets/sets the name of a file to which the coverage filtered exonerate results should be written
 Returns : 
 Args    : 


=cut

sub filtered {
   my ($self, $value) = @_;
   if( defined $value ) {
      $self->{'_filtered'} = $value;
    }
    return $self->{'_filtered'};

}


1;


