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

Bio::EnsEMBL::Pipeline::RunnableDB::Exonerate2Array

=head1 SYNOPSIS

    my $obj = Bio::EnsEMBL::Pipeline::RunnableDB::Exonerate2Array->new(
					     -dbobj     => $db,
					     -input_id  => $id,
					     -analysis   => $analysis 		 
                                             );
    $obj->fetch_input();
    $obj->run();

    my @newfeatures = $obj->output();
    
    $obj->write_output();

=head1 DESCRIPTION

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::EnsEMBL::Pipeline::RunnableDB::Exonerate2Array;

use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::EnsEMBL::Root;
use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Pipeline::Runnable::ExonerateArray;
use Bio::SeqIO;
use Bio::EnsEMBL::Root;
use Data::Dumper;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDB);

sub new {
  my ($class, @args) = @_;
  my $self = bless {}, $class;

  my ($db,$query_file,$target_file,$exonerate,$options) = $self->_rearrange([qw(
										DB
										QUERY_FILE
										TARGET_FILE
										EXONERATE
										OPTIONS
									       )],@args);
  
  $self->db($db);
  $self->query_file($query_file) if defined $query_file;
  $self->target_file($target_file) if defined $target_file;
  $self->exonerate($exonerate) if defined $exonerate;
  $self->options($options) if defined $options;

  return $self; 
}

=head2 fetch_input

    Title   :   fetch_input
    Usage   :   $self->fetch_input
    Function:   Fetches input data for exonerate from the query_file
    Returns :   nothing
    Args    :   none

=cut

sub fetch_input {
  my( $self) = @_;
  print STDERR "Fetching input \n";
  
  $self->throw("No exonerate executable") unless defined($self->exonerate);
  
  my $query_file = $self->query_file();
  my $target_file = $self->target_file();
  my $options = $self->options();

  ###Runnable::ExonerateArray take a array of query_seq_obj, so it's need to be generated here###

  my @query_seqs;

  my $in = Bio::SeqIO->newFh(
			     -FILE => $query_file,
			     -FORMAT => 'Fasta',
			    );

  while (my $seq = <$in>) {
    push (@query_seqs, $seq);
  }

  # prepare runnable
  
  $self->throw("Can't run Exonerate without both query and target sequences") 
    unless (defined($query_file) && defined($target_file));
  
  my $executable =  $self->exonerate();
  my $exonerate = new Bio::EnsEMBL::Pipeline::Runnable::ExonerateArray(
								       '-db'           => $self->db,
								       '-database'     => $target_file,
								       '-query_seqs'   => \@query_seqs,
								       '-query_type'   => 'dna',
								       '-target_type'  => 'dna',
								       '-exonerate'    => $executable,
								       '-options'      => $options,
								       '-analysis'     => $self->analysis,
								      );
  $self->runnable($exonerate);
  
}


=head2 run

    Title   :   run
    Usage   :   $self->run()
    Function:   Runs the exonerate analysis, producing Bio::EnsEMBL::MiscFeature
    Returns :   Nothing, but $self{_output} contains the features
    Args    :   None

=cut

sub run {
  my ($self) = @_;
  
  $self->throw("Can't run - no runnable objects") unless defined($self->runnable);
  
  $self->runnable->run();

}

=head2 output

    Title   :   output
    Usage   :   $self->output()
    Function:   Returns the contents of $self->{_output}, which holds predicted genes.
    Returns :   Array of Bio::EnsEMBL::DnaDnaAlignFeature
    Args    :   None

=cut

=head2 write_output

    Title   :   write_output
    Usage   :   $self->write_output()
    Function:   Writes contents of $self->{_output} into $self->dbobj
    Returns :   1
    Args    :   None

=cut

sub write_output {

  my($self) = @_;
  
  my @features = $self->output();
  my $db = $self->db();
  my $feat_adp = $db->get_DnaAlignFeatureAdaptor;
  
  my ($count,$only_25,$only_24,$both,%count,%done);

  #$feat_adp->store(@features);
  foreach my $feat (@features) {
    #print "seqname is ",$feat->seqname," seq_start is ",$feat->start," seq_end is ",$feat->end, " strand is ",$feat->strand," hseqname is ",$feat->hseqname," percent_id is ",$feat->percent_id,"\n";
    $count++;
    $done{$feat->hseqname}=1;
    if ($feat->percent_id == 100 and ($feat->end - $feat->start +1) ==25) {
      $count{$feat->hseqname}{'25f'}++;
    }
    elsif ($feat->percent_id >=96 and ($feat->end - $feat->start +1) ==25) {
      $count{$feat->hseqname}{'25p'}++;
    }
    elsif (($feat->end - $feat->start +1) ==24) {
      $count{$feat->hseqname}{'24'}++;
    }
  }

  foreach my $q_id (keys %done) {
    if ($count{$q_id}{'25f'} and !$count{$q_id}{'25p'} and !$count{$q_id}{'24'}) {
      $only_25++;
    }
    elsif ($count{$q_id}{'25f'} and $count{$q_id}{'25p'} or $count{$q_id}{'25f'} and $count{$q_id}{'24'}) {
      $both++;
    }
    elsif ($count{$q_id}{'24'} or $count{$q_id}{'25p'} or ($count{$q_id}{'25p'} and $count{$q_id}{'24'})) {
      $only_24++;
    }
  }

  my $tot_pass_ids = keys %done;
  my $ratio_25 = $only_25/$tot_pass_ids;
  my $ratio_24 = $only_24/$tot_pass_ids;
  my $ratio_both = $both/$tot_pass_ids;

  printf "total pass ids is $tot_pass_ids, 25 exact match is $only_25 (%.2f), 24 exact match is $only_24 (%.2f) and both is $both (%.2f)\n", $ratio_25,$ratio_24,$ratio_both;

  return 1;
}

############################################################
#
# get/set methods
#
############################################################

sub query_file {
  
  my ($self,$file) = @_;
  if( defined $file) {
  $self->{'_query_file'} = $file;
}
  return $self->{'_query_file'};
}

############################################################

sub target_file {
  
  my ($self,$file) = @_;
  if( defined $file) {
  $self->{'_target_file'} = $file;
}
  return $self->{'_target_file'};
}


############################################################

sub exonerate {
  my ($self, $location) = @_;
  if ($location) {
    $self->throw("Exonerate not found at $location: $!\n") unless (-e $location);
    $self->{_exonerate} = $location ;
  }
  return $self->{_exonerate};
}

############################################################

sub options {
  my ($self, $options) = @_;
  if ($options) {
    $self->{_options} = $options ;
  }
  return $self->{_options};
}

############################################################

sub db {
  my ($self, $db) = @_;
  if ($db) {
    $self->{_db} = $db ;
  }
  return $self->{_db};
}

############################################################

sub analysis {

  my ($self, $analysis) = @_;

  if ($analysis) {
    $self->{'_analysis'} = $analysis;
  }
  else {
    my ($program, $version) = split /\-/, $self->exonerate;
    
    my $analysis = new Bio::EnsEMBL::Analysis (
					       -db              => 'Affymetrix_array',
					       -db_version      => '',
					       -program         => $program,
					       -program_version => $version,
					       -gff_source      => $program,
					       -gff_feature     => 'mapping',
					       -logic_name      => 'mapping',
					       -module          => 'ExonerateArray',
					      );
    
    $self->{'_analysis'} = $analysis;
  }
  return $self->{'_analysis'};
}

############################################################
    
sub runnable{
   my ($self, $runnable) = @_;
   if($runnable  ) {
     $self->{'_runnable'} = $runnable;
    }
    return $self->{'_runnable'};
}


1;
