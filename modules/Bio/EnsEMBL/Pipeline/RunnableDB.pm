#
#
# Cared for by Michele Clamp  <michele@sanger.ac.uk>
#
# Copyright Michele Clamp
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::RunnableDB

=head1 SYNOPSIS

# get a Bio::EnsEMBL::Pipeline::RunnableDB pbject somehow

  $runnabledb->fetch_input();
  $runnabledb->run();
  $runnabledb->output();
  $runnabledb->write_output(); #writes to DB

=head1 DESCRIPTION

This is the base implementation of
Bio::EnsEMBL::Pipeline::RunnableDBI.  This object encapsulates the
basic main methods of a RunnableDB which a subclass may override.

parameters to new
-db:        A Bio::EnsEMBL::DB::Obj (required), 
-input_id:   Contig input id (required), 
-analysis:  A Bio::EnsEMBL::Analysis (optional) 

This object wraps Bio::EnsEMBL::Pipeline::Runnable to add
functionality for reading and writing to databases.  The appropriate
Bio::EnsEMBL::Analysis object must be passed for extraction
of appropriate parameters. A Bio::EnsEMBL::Pipeline::DBSQL::Obj is
required for databse access.

=head1 CONTACT

ensembl-dev@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Pipeline::RunnableDB;

use strict;
use Bio::EnsEMBL::Pipeline::RunnableDBI;
use Bio::EnsEMBL::Pipeline::SeqFetcher;
use Bio::EnsEMBL::Pipeline::SeqFetcher::Pfetch;
use Bio::EnsEMBL::DBSQL::FeatureAdaptor;
use Bio::DB::RandomAccessI;

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDBI);



=head2 new

    Title   :   new
    Usage   :   $self->new(-DB          => $db,
                           -INPUT_ID    => $id,
                           -SEQFETCHER  => $sf,
			   -ANALYSIS    => $analysis);
                           
    Function:   creates a Bio::EnsEMBL::Pipeline::RunnableDB object
    Returns :   A Bio::EnsEMBL::Pipeline::RunnableDB object
    Args    :   -db:         A Bio::EnsEMBL::DB::Obj (required), 
                -input_id:   Contig input id (required), 
                -seqfetcher: A Bio::DB::RandomAccessI Object (required),
                -analysis:   A Bio::EnsEMBL::Analysis (optional) 
=cut

sub new {
    my ($class, @args) = @_;
    
    my $self = {};
    bless $self, $class;
    #print STDERR "@args\n";
    my ($db,$input_id, $seqfetcher, 
	$analysis) = $self->_rearrange([qw(DB
					   INPUT_ID
					   SEQFETCHER 
					   ANALYSIS )], 
				       @args);
    
    $self->throw("No database handle input") unless defined($db);
    $self->db($db);

    $self->throw("No input id input") unless defined($input_id);
    $self->input_id($input_id);
    
#    if(!defined $seqfetcher) {
#      $seqfetcher = new Bio::EnsEMBL::Pipeline::SeqFetcher::Pfetch;
#    }

    # we can't just default this to pfetch
    $seqfetcher && $self->seqfetcher($seqfetcher);

    # this is an optional field, (I think)
    $analysis && $self->analysis($analysis);
    return $self;
}

=head2 analysis

    Title   :   analysis
    Usage   :   $self->analysis($analysis);
    Function:   Gets or sets the stored Analusis object
    Returns :   Bio::EnsEMBL::Analysis object
    Args    :   Bio::EnsEMBL::Analysis object

=cut

sub analysis {
    my ($self, $analysis) = @_;
    
    if ($analysis)
    {
        $self->throw("Not a Bio::EnsEMBL::Analysis object")
            unless ($analysis->isa("Bio::EnsEMBL::Analysis"));
        $self->{'_analysis'} = $analysis;
        $self->parameters($analysis->parameters);
    }
    return $self->{'_analysis'};
}

=head2 parameters

    Title   :   parameters
    Usage   :   $self->parameters($param);
    Function:   Gets or sets the value of parameters
    Returns :   A string containing parameters for Bio::EnsEMBL::Runnable run
    Args    :   A string containing parameters for Bio::EnsEMBL::Runnable run

=cut

sub parameters {
    my ($self, $parameters) = @_;
    $self->analysis->parameters($parameters) if ($parameters);
    return $self->analysis->parameters();
}

=head2 db

    Title   :   db
    Usage   :   $self->db($obj);
    Function:   Gets or sets the value of db
    Returns :   A Bio::EnsEMBL::Pipeline::DB::ObjI compliant object
                (which extends Bio::EnsEMBL::DB::ObjI)
    Args    :   A Bio::EnsEMBL::Pipeline::DB::ObjI compliant object

=cut

sub db {
    my( $self, $value ) = @_;
    
    if ($value) 
    {
        $value->isa("Bio::EnsEMBL::DBSQL::DBAdaptor")
            || $self->throw("Input [$value] isn't a Bio::EnsEMBL::DBSQ::DBAdpator");
        $self->{'_db'} = $value;
    }
    return $self->{'_db'};
}

=head2 input_id

    Title   :   input_id
    Usage   :   $self->input_id($input_id);
    Function:   Gets or sets the value of input_id
    Returns :   valid input id for this analysis (if set) 
    Args    :   input id for this analysis 

=cut

sub input_id {
    my ($self, $input) = @_;
    if ($input)
    {
        $self->{'_input_id'} = $input;
    }
    return $self->{'_input_id'};
}

=head2 genseq

    Title   :   genseq
    Usage   :   $self->genseq($genseq);
    Function:   Get/set genseq
    Returns :   
    Args    :   

=cut

sub genseq {
    my ($self, $genseq) = @_;

    if (defined($genseq)){ 
	$self->{'_genseq'} = $genseq; 
    }
    return $self->{'_genseq'}
}


=head2 output

    Title   :   output
    Usage   :   $self->output()
    Function:   
    Returns :   Array of Bio::EnsEMBL::FeaturePair
    Args    :   None

=cut

sub output {
    my ($self) = @_;
   
    $self->{'_output'} = [];
    
    my @r = $self->runnable;

    if(@r && scalar(@r)){
      foreach my $r ($self->runnable){
	push(@{$self->{'_output'}}, $r->output);
      }
    }
    return @{$self->{'_output'}};
}

=head2 run

    Title   :   run
    Usage   :   $self->run();
    Function:   Runs Bio::EnsEMBL::Pipeline::Runnable::xxxx->run()
    Returns :   none
    Args    :   none

=cut

sub run {
    my ($self) = @_;

    foreach my $runnable ($self->runnable) {
      $self->throw("Runnable module not set") unless ($runnable);
      $self->throw("Input not fetched") unless ($self->genseq());
      $runnable->query($self->genseq());
      $runnable->run();
    }
}

=head2 runnable

    Title   :   runnable
    Usage   :   $self->runnable($arg)
    Function:   Sets a runnable for this RunnableDB
    Returns :   Bio::EnsEMBL::Pipeline::RunnableI
    Args    :   Bio::EnsEMBL::Pipeline::RunnableI

=cut

sub runnable {
  my ($self,$arg) = @_;

  if (!defined($self->{'_runnables'})) {
      $self->{'_runnables'} = [];
  }
  
  if (defined($arg)) {

      if ($arg->isa("Bio::EnsEMBL::Pipeline::RunnableI")) {
	  push(@{$self->{'_runnables'}},$arg);
      } else {
	  $self->throw("[$arg] is not a Bio::EnsEMBL::Pipeline::RunnableI");
      }
  }
  
  return @{$self->{'_runnables'}};  
}

=head2 vcontig

 Title   : vcontig
 Usage   : $obj->vcontig($newval)
 Function: 
 Returns : value of vcontig
 Args    : newvalue (optional)


=cut

sub vcontig {
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'_vcontig'} = $value;
    }
    return $obj->{'_vcontig'};

}


=head2 slice

 Title   : slice
 Usage   : $obj->slice($newval)
 Function: 
 Returns : value of slice
 Args    : newvalue (optional)


=cut

sub slice {
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'_slice'} = $value;
    }
    return $obj->{'_slice'};

}


=head2 write_output

    Title   :   write_output
    Usage   :   $self->write_output
    Function:   Writes output data to db
    Returns :   array of repeats (with start and end)
    Args    :   none

=cut

sub write_output {
    my($self) = @_;

    my $db=$self->db();
    my @features = $self->output();
    my $contig;

    eval {
      $contig = $db->get_RawContigAdaptor->fetch_by_name($self->input_id);
    };

    if ($@) {
	print STDERR "Contig not found, skipping writing output to db: $@\n";
	return 1;
    }
  
    foreach my $f (@features) {
	$f->analysis($self->analysis);
	$f->attach_seq($contig);
    }

    if (@features) {
	print STDERR "Writing features to database\n";

        my $feat_adp=Bio::EnsEMBL::DBSQL::FeatureAdaptor->new($db);
	print STDERR $feat_adp." is a feature adpator\n";
	$feat_adp->store(@features);
    }
    return 1;
}

=head2 seqfetcher

    Title   :   seqfetcher
    Usage   :   $self->seqfetcher($seqfetcher)
    Function:   Get/set method for SeqFetcher
    Returns :   Bio::DB::RandomAccessI object
    Args    :   Bio::DB::RandomAccessI object

=cut

sub seqfetcher {
  my( $self, $value ) = @_;    
  if (defined($value)) {
    #need to check if passed sequence is Bio::DB::RandomAccessI object
    #$value->isa("Bio::DB::RandomAccessI") || 
    #  $self->throw("Input isn't a Bio::DB::RandomAccessI");
    $self->{'_seqfetcher'} = $value;
  }
    return $self->{'_seqfetcher'};
}

=head2 input_is_void

    Title   :   input_is_void
    Usage   :   $self->input_is_void(1)
    Function:   Get/set flag for sanity of input sequence
                e.g. reject seqs with only two base pairs
    Returns :   Boolean
    Args    :   Boolean

=cut

sub input_is_void {
    my ($self, $value) = @_;

    if ($value) {
	$self->{'_input_is_void'} = $value;
    }
    return $self->{'_input_is_void'};

}

1;
