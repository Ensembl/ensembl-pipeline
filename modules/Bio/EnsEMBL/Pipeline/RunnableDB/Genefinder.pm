
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

Bio::EnsEMBL::Pipeline::RunnableDB::Genefinder

=head1 SYNOPSIS

my $db          = Bio::EnsEMBL::DBLoader->new($locator);
my $genefinder     = Bio::EnsEMBL::Pipeline::RunnableDB::Genefinder->new ( 
                                                    -dbobj      => $db,
			                                        -input_id   => $input_id
                                                    -analysis   => $analysis );
$genefinder->fetch_input();
$genefinder->run();
$genefinder->output();
$genefinder->write_output(); #writes to DB

=head1 DESCRIPTION

This object wraps Bio::EnsEMBL::Pipeline::Runnable::Genefinder to add
functionality for reading and writing to databases. The appropriate
Bio::EnsEMBL::Analysis object must be passed for extraction
of appropriate parameters. A Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor
is required for database access.

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Pipeline::RunnableDB::Genefinder;

use strict;
use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Pipeline::Runnable::Genefinder;
use Data::Dumper;

use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDB);

=head2 new

    Title   :   new
    Usage   :   $self->new(-DBOBJ       => $db
                           -INPUT_ID    => $id
                           -ANALYSIS    => $analysis);      
                          
    Function:   creates a Bio::EnsEMBL::Pipeline::RunnableDB::Genefinder object
    Returns :   A Bio::EnsEMBL::Pipeline::RunnableDB::Genefinder object
    Args    :   -db:     A Bio::EnsEMBL::DBSQL::DBAdaptor, 
                input_id:   Contig input id , 
                -analysis:  A Bio::EnsEMBL::Analysis 

=cut

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);
    
    $self->{'_fplist'}      = [];
    $self->{'_genseq'}      = undef;
    $self->{'_runnable'}    = undef;
    $self->{'_parameters'}  = undef;
    
    # db input_id mandatory and read in by BlastableDB

    # anlaysis not mandatory for BlastableDB, so we check here 
    $self->throw("Analysis object required") unless ($self->analysis);
    $self->init('Bio::EnsEMBL::Pipeline::Runnable::Genefinder');
    
    return $self;
}


=head2 fetch_input

    Title   :   fetch_input
    Usage   :   $self->fetch_input
    Function:   Fetches input data for repeatmasker from the database
    Returns :   none
    Args    :   none

=cut

sub fetch_input {
    my($self) = @_;
    
    $self->throw("No input id") unless defined($self->input_id);

    my $contigid  = $self->input_id;
    my $contig    = $self->db->get_RawContigAdaptor->fetch_by_name($contigid);
    my $genseq    = $contig->get_repeatmasked_seq() or $self->throw("Unable to fetch contig");
    $self->{'contig'} = $contig;
    $self->genseq($genseq);

}

#get/set for runnable and args
sub init {
    my ($self, $runnable) = @_;
    my %parameters;
    if ($runnable) {
      #extract parameters into a hash
    
      my ($parameter_string) = $self->parameters();
      if ($parameter_string) {
	my @pairs = split (/,/, $parameter_string);
	foreach my $pair (@pairs) {
	  
	  my ($key, $value) = split (/=>/, $pair);
	  $key =~ s/\s+//g;
	  $parameters{$key} = $value;
	}
	
      }
   
      #creates empty Bio::EnsEMBL::Runnable::Genefinder object
      my $runnable = $runnable->new(%parameters);
      $self->runnable($runnable);
    }
}


=head2 result_quality_tag

    Title   :   result_quality_tag
    Usage   :   $self->result_quality_tag
    Function:   Returns an indication of whether the data is suitable for 
                further analyses. Allows distinction between failed run and 
                no hits on a sequence.
    Returns :   'VALID' or 'INVALID'
    Args    :   none

=cut
#a method of writing back result quality
sub result_quality_tag {
    my ($self) = @_;
    
    if ($self->output)
    {
        return 'VALID';
    }
    else
    {
        return 'INVALID';
    }
}


sub write_output {
   my $self = shift;
  
   my $genefinder_runnable = ($self->runnable())[0];
   my @transcripts = $genefinder_runnable->each_Transcript();
   if( ! @transcripts ) { return; }

   my $ptransAdaptor = $self->db()->get_PredictionTranscriptAdaptor();
   print "there are ".@transcripts." transcripts\n";
   for my $trans ( @transcripts ) {
     my $ptrans = Bio::EnsEMBL::PredictionTranscript->new();
     my @exons = $trans->get_all_Exons();

     if ($exons[0]->strand == 1) {
       @exons = sort {$a->start <=> $b->start } @exons;
     } else {
       @exons = sort {$b->start <=> $a->start } @exons;
     }
     #print "ANALYSIS: ",$self->analysis()->dbID,"\n";

     $ptrans->analysis( $self->analysis() );
     for my $exon ( @exons ) {
       print STDERR "adding contig ".$self->{'contig'}." to exon\n";
       $exon->contig( $self->{'contig'} );
       $ptrans->add_Exon( $exon );
     }
     $ptransAdaptor->store( $ptrans );
   }
  
}

1;






