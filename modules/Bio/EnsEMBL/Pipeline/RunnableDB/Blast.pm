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

Bio::EnsEMBL::Pipeline::RunnableDB::Blast

=head1 SYNOPSIS

my $db      = Bio::EnsEMBL::DBLoader->new($locator);
my $blast   = Bio::EnsEMBL::Pipeline::RunnableDB::Blast->new ( 
                                                    -dbobj      => $db,
			                                        -input_id   => $input_id
                                                    -analysis   => $analysis );
$blast->fetch_input();
$blast->run();
$blast->output();
$blast->write_output(); #writes to DB

=head1 DESCRIPTION

This object wraps Bio::EnsEMBL::Pipeline::Runnable::Blast to add
functionality for reading and writing to databases.
The appropriate Bio::EnsEMBL::Analysis object must be passed for
extraction of appropriate parameters. A Bio::EnsEMBL::Pipeline::DBSQL::Obj is
required for databse access.

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Pipeline::RunnableDB::Blast;

use strict;
use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Pipeline::Runnable::Blast;

use vars qw(@ISA);

BEGIN {
    require "Bio/EnsEMBL/Pipeline/pipeConf.pl";
    require "Bio/EnsEMBL/Pipeline/Blast_conf.pl";
}

@ISA = qw (Bio::EnsEMBL::Pipeline::RunnableDB);

=head2 new

    Title   :   new
    Usage   :   $self->new(-DBOBJ       => $db
                           -INPUT_ID    => $id
                           -ANALYSIS    => $analysis);
                           
    Function:   creates a Bio::EnsEMBL::Pipeline::RunnableDB::Blast object
    Returns :   A Bio::EnsEMBL::Pipeline::RunnableDB::Blast object
    Args    :   -dbobj:     A Bio::EnsEMBL::DBSQL::DBAdaptor, 
                -input_id:   Contig input id , 
                -analysis:  A Bio::EnsEMBL::Analysis 

=cut

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);
    
    $self->{'_fplist'}      = [];
    $self->{'_genseq'}      = undef;
    $self->{'_runnable'}    = undef;            
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

    
    $self->genseq($genseq);
  
# input sequence needs to contain at least 3 consecutive nucleotides
    my $seq = $self->genseq->seq;
    if ($seq =~ /[CATG]{3}/) {
        $self->input_is_void(0);
    }
    else {
        $self->input_is_void(1);
        $self->warn("Need at least 3 nucleotides");
    }
   
}

#get/set for runnable and args
sub runnable {
    my ($self) = @_;
    
    if (!defined($self->{'_runnable'})) {
      my $ungapped;
      if($::pipeConf{'ungapped'}){
	$ungapped = 1;
      }else{
	$ungapped = undef;
      }
      #print STDERR "ungapped = ".$ungapped."\n";
      my $run = Bio::EnsEMBL::Pipeline::Runnable::Blast->new(-query     => $self->genseq,
							     -database  => $self->analysis->db_file,
							     -program   => $self->analysis->program,
							     -options => $self->analysis->parameters,
							     -threshold_type => 'PVALUE',
							     -threshold => 1,
							     -ungapped => $ungapped,
							    );

      $self->{'_runnable'} = $run;
    }
    
    return $self->{'_runnable'};
}

=head2 run

    Title   :   run
    Usage   :   $self->run();
    Function:   Runs Bio::EnsEMBL::Pipeline::Runnable::Blast->run()
    Returns :   none
    Args    :   none

=cut

sub run {
    my ($self,$dir) = @_;
   
    $self->throw("Runnable module not set") unless ($self->runnable());
    $self->throw("Input not fetched")       unless ($self->genseq());

    $self->runnable->run($dir);
  
}


=head2 output

    Title   :   output
    Usage   :   $self->output();
    Function:   Runs Bio::EnsEMBL::Pipeline::Runnable->output()
    Returns :   An array of Bio::EnsEMBL::FeaturePair objects
    Args    :   none

=cut

sub output {
    my ($self) = @_;

    my $runnable = $self->runnable;
    $runnable || $self->throw("Can't return output - no runnable object");

    return $runnable->output;
}


sub write_output{
  my ($self) = @_;

  my @features = $self->output();
  my $dna_f_a = $self->db->get_DnaAlignFeatureAdaptor();
  my $pep_f_a = $self->db->get_ProteinAlignFeatureAdaptor();
  my $contig;
  eval 
    {
      $contig = $self->db->get_RawContigAdaptor->fetch_by_name($self->input_id);
    };

  if ($@) {
      print STDERR "Contig not found, skipping writing output to db: $@\n";
      return 1;
  }
  foreach my $f(@features){
    $f->analysis($self->analysis);
    $f->attach_seq($contig);
    if($f->isa('Bio::EnsEMBL::DnaDnaAlignFeature')){
      $dna_f_a->store($f);
    }elsif($f->isa('Bio::EnsEMBL::DnaPepAlignFeature')){
      $pep_f_a->store($f);
    }else{
      $self->throw("don't know how to store $f\n");
    }
  }


}

1;
