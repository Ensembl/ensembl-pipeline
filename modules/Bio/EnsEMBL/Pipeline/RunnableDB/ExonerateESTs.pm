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

Bio::EnsEMBL::Pipeline::RunnableDB::ExonerateESTs

=head1 SYNOPSIS

    my $obj = Bio::EnsEMBL::Pipeline::RunnableDB::ExonerateESTs->new(
								     -dbobj     => $db,
								     -input_id  => $id,
								     -estfile   => $estfile,
								     -exonerate => $exonerate,
								     -exonerate_args => $exargs
								    );
    $obj->fetch_input
    $obj->run

    my @genes = $obj->output;


=head1 DESCRIPTION
Just runs Exonerate over a  chunk of dbEST and spits the output to STDOUT to be assessed later.
File1 => all FeaturePairs produced by Exonerate

Can run either over virtual contigs (chr_name.start-end) or on a file containing genomic sequence(s)

The running of this process is controlled by the scripts in ensembl-pipeline/scripts/EST

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

# Object preamble
use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Pipeline::Runnable::ExonerateESTs;

@ISA = qw( Bio::EnsEMBL::Pipeline::RunnableDB );


=head2 new

    Title   :   new
    Usage   :   $self->new(-DBOBJ       => $db,
                           -INPUT_ID    => $id,
                           -ANALYSIS    => $analysis,
                           -ESTFILE     => $estfile,
			   -EXONERATE   => $exonerate,
			   -EXONERATE_ARGS => $exargs,
			  );
                           
    Function:   creates a 
                Bio::EnsEMBL::Pipeline::RunnableDB::ExonerateESTs
                object
    Returns :   A Bio::EnsEMBL::Pipeline::RunnableDB::ExonerateESTs
                object
    Args    :   -db:           A Bio::EnsEMBL::DBSQL::DBAdaptor (required), 
                -input_id:        Contig input id (required), or filename
                -analysis:        A Bio::EnsEMBL::Analysis (optional)
                -estfile:         filename
                -exonerate:       path to exonerate executable (optional)
                -exonerate_args : string (optional) 
=cut

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);
#    print "@args\n";       
    # db, input_id, seqfetcher, and analysis objects are all set in
    # in superclass constructor (RunnableDB.pm)

    my( $estfile, $exonerate, $exargs ) = $self->_rearrange([qw(ESTFILE
								EXONERATE
								EXONERATE_ARGS)], @args);

    $self->throw("No est file specified") unless defined($estfile);
    $self->estfile($estfile);

    $self->exonerate($exonerate) if defined $exonerate;

    # ought not to hard code
    $exargs = " -w 14 -t 65 -H 100 -D 15 -m 500 " unless defined $exargs;
    $self->exonerate_args($exargs) if defined $exargs;

    return $self;
}

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

=head2 exonerate

 Title   : exonerate
 Usage   : $obj->exonerate($exonerate)
 Function: get/set for exonerate
 Returns : path to exonerate
 Args    : exonerate (optional)


=cut

sub exonerate {
   my ($self, $exonerate) = @_;

   if (!defined $self->{'_exonerate'}){
      $self->{'_exonerate'} = "";
   }

   if( defined $exonerate ) {
      $self->{'_exonerate'} = $exonerate;
    }
    return $self->{'_exonerate'};

}


=head2 exonerate_args

 Title   : exonerate_args
 Usage   : $obj->exonerate_args($exargs)
 Function: get/set for arguments to exonerate
 Returns : value of exonerate_args
 Args    : exargs (optional)


=cut

sub exonerate_args {
   my ($self, $exargs) = @_;

   if (!defined $self->{'_exonerate_args'}){
      $self->{'_exonerate_args'} = "";
   }

   if( defined $exargs ) {
      $self->{'_exonerate_args'} = $exargs;
    }
    return $self->{'_exonerate_args'};

}

=head2 write_output

    Title   :   write_output
    Usage   :   $self->write_output
    Function:   Does nothing. Output is written to STDOUT by Exonerate runnable
    Returns :   
    Args    :   none

=cut

sub write_output {

  my ($self) = @_;

}

=head2 fetch_input

    Title   :   fetch_input
    Usage   :   $self->fetch_input
    Function:   Fetches input data for ExonerateESTs and makes runnable. Input genomic data 
                can either be a chromosomal range specified as chr_name.start-end, or can be 
                the name of a file containing genomic sequence(s)
    Returns :   nothing
    Args    :   none

=cut

sub fetch_input {
  my ($self) = @_;
  
  $self->throw("No input id") unless defined($self->input_id);
  
  my $runnable;
  my $chrid  = $self->input_id;

  if($chrid =~ s/\.(.*)-(.*)//){
    # we're using virtual contigs
    my $chrstart = $1;
    my $chrend   = $2;
    
    my $stadaptor = $self->db->get_StaticGoldenPathAdaptor();
    my $contig    = $stadaptor->fetch_VirtualContig_by_chr_start_end($chrid,$chrstart,$chrend);
    
    $contig->_chr_name($chrid);
    
    my $genseq    = $contig->get_repeatmasked_seq;
  
    $runnable  = new Bio::EnsEMBL::Pipeline::Runnable::ExonerateESTs('-genomic'        => $genseq, 
								     '-ests'           => $self->estfile,
								     '-exonerate'      => $self->exonerate,
								     '-exonerate_args' => $self->exonerate_args);
    $self->vc($contig);
  }

  else {
    # assume it's a file of genomic sequences - does it exist?
    $self->throw("$chrid cannot be parsed and is not a file\n") unless -e $chrid;
    $runnable  = new Bio::EnsEMBL::Pipeline::Runnable::ExonerateESTs('-genomic'        => $chrid, 
								     '-ests'           => $self->estfile,
								     '-exonerate'      => $self->exonerate,
								     '-exonerate_args' => $self->exonerate_args);
  }

  $self->runnable($runnable);
}

=head2 run

    Title   :   run
    Usage   :   $self->run
    Function:   Calls run method of each runnable
    Returns :   Nothing
    Args    :   None

=cut

sub run {
  my ($self) = @_;

  $self->throw("Can't run - no runnable objects") unless defined($self->{_runnables});
  
  foreach my $runnable($self->runnable) {
    $runnable->run;

    # all the output goes to STDOUT
  }

}

=head2 output

    Title   :   output
    Usage   :   $self->output
    Function:   Returns output from this RunnableDB
    Returns :   Array of Bio::EnsEMBL::FeaturePair
    Args    :   None

=cut

sub output {
   my ($self, $feat) = @_;

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

=cut

1;
