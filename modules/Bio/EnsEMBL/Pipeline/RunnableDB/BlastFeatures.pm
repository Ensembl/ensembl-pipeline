#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::RunnableDB::BlastFeatures

=head1 SYNOPSIS

my $db      = Bio::EnsEMBL::DBAdaptor->new(@args);
my $blast   = Bio::EnsEMBL::Pipeline::RunnableDB::BlastFeatures->new ( 
    -db           => $db,
    -input_id     => $input_id
    -analysis     => $analysis
    -source => $source_logic_name,
    -seqtype  => 'dna',
    -seqfetcher => $seqfetcher
    @other_blast_options );

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

package Bio::EnsEMBL::Pipeline::RunnableDB::BlastFeatures;

use strict;
use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Pipeline::Runnable::Blast;
use Bio::EnsEMBL::Pipeline::Config::Blast;
use Bio::EnsEMBL::Pipeline::Config::General;
use vars qw(@ISA);



@ISA = qw (Bio::EnsEMBL::Pipeline::RunnableDB);

my %UNGAPPED;
my %UNMASKED;
my %SEQFETCHER;
my %INDEX;

foreach my $db (@$DB_CONFIG) {
  my ($name, $ungapped, $unmasked, $seqfetcher, $index) = ($db->{'name'}, 
							   $db->{'ungapped'}, 
							   $db->{'min_unmasked'},
							   $db->{'seqfetcher'},
							   $db->{'index'});
  
  if($db && $name){
      $UNGAPPED{$name} = $ungapped;
      $UNMASKED{$name} = $unmasked;
      $SEQFETCHER{$name} = $seqfetcher;
      $INDEX{$name} = $index;
  }else{
      my($p, $f, $l) = caller;
      warn("either db ".$db." or name ".$name." isn't defined so can't work $f:$l\n");
  }
}

=head2 new

    Title   :   new
    Usage   :   
    Function:   
    Returns :   
    Args    :   

=cut

sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);
  
  my ($logic, $logic_type, $seqfetcher, $seqfetchermod, $indexfile);

  ($logic, $logic_type, $seqfetcher ) =  $self->_rearrange([qw(SOURCE
							       SEQTYPE
							       SEQFETCHER])], @args);
  
  $self->source_logic( $logic );
  $self->seq_type( $logic_type );
  $self->seqfetcher( $seqfetcher );
  
  # analysis table 'parameters' field is looked to to obtain missing parameters
  ($logic, $logic_type, $seqfetchermod, $indexfile ) =  $self->_rearrange([qw(SOURCE
									      SEQTYPE
									      SEQFETCHER
									      INDEX)], $self->parameter_hash);
  							   
  $self->source_logic( $logic ) if not $self->source_logic;
  $self->seq_type( $logic_type ) if not $self->seq_type;
  
  $self->throw("You must give a source logic_name") if not $self->source_logic;
  $self->seq_type('PROTEIN') if not $logic_type;

  if ($seqfetchermod and $indexfile) {
      $self->seqfetcher( $self->make_seqfetcher($indexfile, $seqfetchermod) );
  }

  # if we stil have no seqfetcher, look to the Blast Config

  if (not $self->seqfetcher and 
      $SEQFETCHER{$self->analysis->db_file} and
      $INDEX{$self->analysis->db_file}) {

      $self->seqfetcher( $self->make_seqfetcher($INDEX{$self->analysis->db_file},
						$SEQFETCHER{$self->analysis->db_file}) );
  }

  $self->throw("You must give a seqfetcher (options: object creation, parameters in analysis table, or Blast config)")
      if not $self->seqfetcher;
	    
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
    
    my $contig    = $self->db->get_RawContigAdaptor->fetch_by_name($self->input_id);
    my $genseq;
    if(@$PIPELINE_REPEAT_MASKING){
	my $genseq    = $contig->get_repeatmasked_seq($PIPELINE_REPEAT_MASKING) or $self->throw("Unable to fetch contig");
	$self->query($genseq);
    }else{
	$self->query($contig);
    }
    
    my $seq = $self->query->seq;
    my $unmasked;
    if($UNMASKED{$self->analysis->db_file}){
	$unmasked = $UNMASKED{$self->analysis->db_file};
    } else {
	$unmasked = 3;
    }
    if ($seq =~ /[CATG]{$unmasked}/) {
        $self->input_is_void(0);
    } else {
        $self->input_is_void(1);
        $self->warn("Need at least $UNMASKED{$self->analysis->db_file} nucleotides");
    }
    
    my $ungapped;
    
    if($UNGAPPED{$self->analysis->db_file}){
	$ungapped = 1;
    } else {
	$ungapped = undef;
    }
    
    my (%id_list, @flist);
    if ($self->seq_type eq 'PROTEIN') {
	@flist = @{$contig->get_all_ProteinAlignFeatures( $self->source_logic )};
    }
    else {
	@flist = @{$contig->get_all_DnaAlignFeatures( $self->source_logic )} 
    }
    
    map { $id_list{$_->hseqname} = 1 } @flist;
	  
    my $run = Bio::EnsEMBL::Pipeline::Runnable::Blast->new(-query          => $self->query,
							   -database       => $self->analysis->db_file,
							   -ids            => [keys %id_list],
							   -seqfetcher     => $self->seqfetcher,
							   -program        => $self->analysis->program,
							   -threshold_type => 'PVALUE',
							   -threshold      => 1,
							   -ungapped       => $ungapped,
							   $self->parameter_hash
							   # allows parameters to be passed 
							   # from analysis 'parameters' field
							   );

    $self->runnable($run);

    return 1;
}


sub source_logic {
    my ($self, $logic) = @_;

    if ($logic) {
	$self->{'_blast_restricted_logic'} = $logic;
    }    
    return $self->{'_blast_restricted_logic'};
}



sub seq_type {
    my ($self, $type) = @_;

    if ($type) {
	$self->{'_blast_restricted_seqtype'} = $type;
    }    
    return $self->{'_blast_restricted_seqtype'};
}



sub make_seqfetcher {
  my ( $self, $index, $seqfetcher_class ) = @_;
  my $seqfetcher;

  (my $class = $seqfetcher_class) =~ s/::/\//g;
  require "$class.pm";

  if(defined $index && $index ne ''){
    my @db = ( $index );
    
    # make sure that your class is compatible with the index type
    $seqfetcher = "$seqfetcher_class"->new('-db' => \@db, );
    
  }
  else{
      $self->throw("Can't make seqfetcher\n");
  }

  return $seqfetcher;
}



1;


