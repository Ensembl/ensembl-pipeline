
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
                                                    -db      => $db,
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
use Bio::EnsEMBL::PredictionTranscript;
use Data::Dumper;

use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDB);


=head2 new

  Arg [1]   : RunnableDB standard args
  Function  : creates a Genefinder Runnable DB
  Returntype: self
  Exceptions: is no analysis object is passed in
  Caller    : 
  Example   : my $runnabledb = Bio::EnsEMBL::Pipeline::RunnableDB::GeneFinder->new(-input_id => $seq,
                                                                                   -analysis => $analysis
                                                                                   -db => $db);
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
    #$self->init('Bio::EnsEMBL::Pipeline::Runnable::Genefinder');
    
    return $self;
}


=head2 fetch_input

  Arg [1]   : none
  Function  : get data for runnable from database
  Returntype: non
  Exceptions: throws if no input id is set
  Caller    : 
  Example   : 

=cut




sub fetch_input {
    my($self) = @_;
    #print STDERR "in genefinder\n";
    #print STDERR "using db ".$self->db->dbname."\n";
    #print STDERR "input id ".$self->input_id."\n";
    $self->throw("No input id") unless ($self->input_id);
    #my $contigid  = $self->input_id;
    my $contig    = $self->db->get_RawContigAdaptor->fetch_by_name($self->input_id);
    #print STDERR "have ".$contig." contig\n";
    my $genseq    = $contig->get_repeatmasked_seq() or $self->throw("Unable to fetch contig");
    
#    print STDERR "have ".$genseq." repeatmasked seq\n";
    $self->{'contig'} = $contig;
    $self->query($genseq);

    my $runnable = Bio::EnsEMBL::Pipeline::Runnable::Genefinder->new('-query'=> $self->query);
    $self->runnable($runnable);

}

#get/set for runnable and args



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


=head2 create_PredictionTranscripts

  Arg [1]   : none
  Function  : makes Prediction transcripts from results from runnable
  Returntype: the valid transcripts
  Exceptions: none
  Caller    : $self
  Example   : $self->create_PredictionTranscripts

=cut


sub create_PredictionTranscripts{
  my $self = shift;

  my @ptranscripts;
  my @checked_transcripts;
  my $genefinder_runnable = ($self->runnable())[0];
  my @transcripts = $genefinder_runnable->each_Transcript();
  if( ! @transcripts ) { return; }
  
  for my $trans ( @transcripts ) {
    my $ptrans = Bio::EnsEMBL::PredictionTranscript->new();
    my @exons = @{$trans->get_all_Exons()};
    
    if ($exons[0]->strand == 1) {
      @exons = sort {$a->start <=> $b->start } @exons;
    } else {
      @exons = sort {$b->start <=> $a->start } @exons;
    }
    #print "ANALYSIS: ",$self->analysis()->dbID,"\n";
    
    $ptrans->analysis( $self->analysis() );
    for my $exon ( @exons ) {
      #print STDERR "adding contig ".$self->{'contig'}." to exon\n";
      $exon->contig( $self->{'contig'} );
      $ptrans->add_Exon( $exon );
    }
    push(@ptranscripts, $ptrans);
  }
  
  foreach my $pt(@ptranscripts){
    my $checked = $self->check_translation($pt);
    if($checked){
      push(@checked_transcripts, $checked);
    }
  }
  return @checked_transcripts;
}


=head2 write_output

  Arg [1]   : none
  Function  : writes output to db
  Returntype: none
  Exceptions: none
  Caller    : 
  Example   : 

=cut


sub write_output {
   my $self = shift;
  
  
   my @transcripts = $self->create_PredictionTranscripts;
   if( ! @transcripts ) { return; }

   my $ptransAdaptor = $self->db()->get_PredictionTranscriptAdaptor();
   for my $trans ( @transcripts ) {
     $ptransAdaptor->store( $trans );
   }
  
}


=head2 check_translation

  Arg [1]   : Bio::EnsEMBL::PredictionTranscript
  Function  : check if transcript translates
  Returntype: valid transcript or undef
  Exceptions: none
  Caller    : 
  Example   : 

=cut


sub check_translation{
  my ($self, $transcript) = @_; 
  my @exons = @{$transcript->get_all_Exons};
  if($exons[0]->strand == 1){
    @exons = sort{$a->start <=> $b->start} @exons;
  }else{
    @exons = sort{$b->start <=> $a->start} @exons;
  }
  print STDERR "checking translation\n";
  #foreach my $e(@exons){
  #  print STDERR $e->seqname."\t ".$e->start."\t ".$e->end."\t ".$e->strand."\t ".$e->phase."\t ".$e->end_phase."\n";
  #}

  my $pep = $transcript->translate->seq;
 
  my $translate = 1;
  if($pep =~ /\*/){
    print STDERR "translation ".$pep." doesn't translate\n";
    return undef;
  }
  
  if(!$translate){
    $self->throw("there isn't a translation without stops the code shouldn't have got here $!");
  }else{
    return $transcript;
  }
}




1;






