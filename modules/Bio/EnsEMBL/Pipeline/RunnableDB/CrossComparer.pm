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

Bio::EnsEMBL::Pipeline::RunnableDB::CrossComparer

=head1 SYNOPSIS

    my $obj = new Bio::EnsEMBL::Pipeline::RunnableDB::CrossComparer(-db => $db
								    -input_id => $id,
								    -alnprog => 'bl2seq',
								    -alntype => 'blastn',
								    -min_score => 40,
								    -masked => 1,
								    -filter => "Bio::EnsEMBL::Compara::Filter::Greedy",
								    -options => $analysis->parameters);

    $obj->fetch_input();
    $obj->run();

    my @newfeatures = $obj->output();
    
    $obj->write_output();

=head1 DESCRIPTION

    This runnabledb creates an alignment runnable needed to 
    compare two DNA fragments (dnafrag could be a RawContig or a VirtualContig as Slice objects)
    from databases of different organisms

    This runnabledb is more complex than usual in that it needs
    to be connected to several dbs:

    alnprog: string defining the alignment runnable used. Can be 'crossmatch', 
             'bl2seq', 'blastz' or 'exonerate' (the latter not implemented at the moment).
             'bl2seq' is the default.

    -db: this is where the output is finally written, 
    storing the hits between two Ensembl databases. 
    It is actually a Bio::EnsEMBL::Compara::DBSQL::DBAdaptor 
    database, not a normal pipeline EnsEMBL database


    The input id has the following format: 

    sequence_source1:species1:dnafrag_type1:dnafrag_name1::sequence_source2:species2:dnafrag_type2:dnafrag_name2

    species1 and species2: these are the two species for which the 
    input dna has to be compared

    e.g.
    ENSEMBL:Homo_sapiens:VirtualContig:1.1.250000::ENSEMBL:Mus_musculus:VirtualContig:4.151500001.151730910

    dnafrag_name1 being considered as the reference contig to define the 
    reference AlignBlockSet

    -min_score is an common argument for all alnprog, and throw any alignment below min_score
    -masked 0:unmasked 1:masked (repeat replaced by Ns) 2:soft masked (repeats in lower case)
    -filter specify a perl module which takes an array of DnaDnaAlignFeature objects and 
            returns an array of DnaDnaAlignFeature objects after filtering features
    -options takes a string corresponding to options/parameters to be added in the executable 
             command line.

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::EnsEMBL::Pipeline::RunnableDB::CrossComparer;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Pipeline::Runnable::CrossMatch;
use Bio::EnsEMBL::Pipeline::Runnable::Exonerate;
use Bio::EnsEMBL::Pipeline::Runnable::Bl2seq;
use Bio::EnsEMBL::Pipeline::Runnable::Blastz;
use Bio::PrimarySeq;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDB);

sub new {
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);
  
  $self->{'_filter'} = 0;
  $self->{'_masked'} = 0;
  $self->{'_alnprog'} = 'bl2seq';
  $self->{'_alntype'} = 'blastn';
  $self->{'_min_score'} = 100;

  my ($min_score, $alnprog, $alntype, $masked, $filter, $options) = 
    $self->_rearrange([qw(MIN_SCORE
			  ALNPROG
			  ALNTYPE
			  MASKED
			  FILTER
			  OPTIONS)],@args);
  
  bless $self, $class;

  $self->min_score($min_score) if (defined $min_score);
  $self->alnprog($alnprog) if (defined $alnprog);
  $self->alntype($alntype) if (defined $alntype);
  $self->masked($masked) if (defined $masked);
  $self->options($options) if (defined $options);

  if (defined $filter) {
    $self->filter($filter);
    eval "require $filter";
    if($@) {
      $self->warn("$filter cannot be found.\nException $@\n");
      return undef;
    }
  }

  return $self; 
}

=head2 fetch_input

    Title   :   fetch_input
    Usage   :   $self->fetch_input
    Function:   Fetches input dna for crossmatch from the databases
    Returns :   nothing
    Args    :   none

=cut

sub fetch_input {
    my($self) = @_; 

    $self->throw("No input id") unless defined($self->input_id);
    $self->throw("alnprog should be 'crossmatch' or 'bl2seq'or 'exonerate'")
      unless (defined $self->alnprog &&
	      ($self->alnprog eq 'crossmatch' ||
	       $self->alnprog eq 'bl2seq' ||
	       $self->alnprog eq 'exonerate' ||
	       $self->alnprog eq 'blastz'));

    my $input_id  = $self->input_id;

    my ($sequence_source1,$species1,$dnafrag_type1,$dnafrag_name1,$sequence_source2,$species2,$dnafrag_type2,$dnafrag_name2);
    if ($input_id =~ /^(\S+):(\S+):(\S+):(\S+)::(\S+):(\S+):(\S+):(\S+)$/ ) {
      $sequence_source1 = $1;
      $species1 = $2;
      $dnafrag_type1 = $3;
      $dnafrag_name1 = $4;
      $sequence_source2 = $5;
      $species2 = $6;
      $dnafrag_type2 = $7;
      $dnafrag_name2 = $8;
    }
    else {
	$self->throw("Input id not in correct format: got $input_id, should be parsable by 
\/^\\S+:\\S+:\\S+:\\S+::\\S+:\\S+:\\S+\:\\S+\$\/");
    }

    $self->_c1_id($dnafrag_name1);
    $self->_c2_id($dnafrag_name2);

    my $gadp = $self->db->get_GenomeDBAdaptor();
    
    my ($contig1,$contig2);
    my ($seq1,$seq2);

    if ($sequence_source1 eq "ENSEMBL") {

      my $core_species1adp = $gadp->fetch_by_species_tag($species1)->db_adaptor;

      if ($dnafrag_type1 eq "RawContig") {

	$contig1 = $core_species1adp->get_RawContigAdaptor->fetch_by_name($dnafrag_name1);

      } elsif ($dnafrag_type1 eq "VirtualContig") {

	my ($chr,$start,$end) = split /\./, $dnafrag_name1;
	$contig1 = $core_species1adp->get_SliceAdaptor->fetch_by_chr_start_end($chr,$start,$end);

      }

      if ($self->masked == 1) {
	$seq1 = $contig1->get_repeatmasked_seq;
	$seq1->display_id("seq1");
      } elsif ($self->masked == 2) {
	$seq1 = $contig1->get_repeatmasked_seq('RepeatMask',1);
	$seq1->display_id("seq1");
      } else {
	$seq1 = Bio::PrimarySeq->new( -display_id => 'seq1', -seq => $contig1->seq);
      }

    } elsif ($sequence_source1 eq "FASTA") {

      $seq1 = "$dnafrag_name1";

    }

    if ($sequence_source2 eq "ENSEMBL") {

      my $core_species2adp = $gadp->fetch_by_species_tag($species2)->db_adaptor;

      if ($dnafrag_type2 eq "RawContig") {

	$contig2 = $core_species2adp->get_RawContigAdaptor->fetch_by_name($dnafrag_name2);

      } elsif ($dnafrag_type2 eq "VirtualContig") {

	my ($chr,$start,$end) = split /\./, $dnafrag_name2;
	$contig2 = $core_species2adp->get_SliceAdaptor->fetch_by_chr_start_end($chr,$start,$end);

      } 

      if ($self->masked == 1) {
	$seq2 = $contig2->get_repeatmasked_seq;
	$seq2->display_id("seq2");
      } elsif ($self->masked == 2) {
	$seq2 = $contig2->get_repeatmasked_seq('RepeatMask',1);
	$seq2->display_id("seq2");
      } else {
	$seq2 = Bio::PrimarySeq->new( -display_id => 'seq2', -seq => $contig2->seq);
      }

    } elsif ($sequence_source2 eq "FASTA") {

      $seq2 = "$dnafrag_name2";

    }
    
    
    my $alnrunnable;

    if ($self->alnprog eq 'crossmatch') {
      $alnrunnable = Bio::EnsEMBL::Pipeline::Runnable::CrossMatch->new(-nocopy => 1,
								       -seq1 => $seq1,
								       -seq2 => $seq2,
								       -score => $self->min_score,
								       -minmatch => 14,
								       -masklevel => 101);
    } elsif ($self->alnprog eq 'bl2seq') {

# bl2seq runnable runs by default a 'blastn' alignment type.
      
      $alnrunnable = Bio::EnsEMBL::Pipeline::Runnable::Bl2seq->new(-seq1 => $seq1,
								   -seq2 => $seq2,
								   -alntype => $self->alntype,
								   -min_score => $self->min_score,
								   -options => $self->analysis->parameters
								  );
    } elsif ($self->alnprog eq 'exonerate') {
      $self->throw("exonerate aligment runnable not implemented yet");
#      $alnrunnable = Bio::EnsEMBL::Pipeline::Runnable::Exonerate->new(-exonerate => "/usr/local/ensembl/bin/exonerate",
#								      -genomic => $seq2,
#								      -est => [$seq1],
#								      -args => "-a 400 -p no -n -w 14 -t 65 -H 150 -D 5 -m 768");
    } elsif ($self->alnprog eq 'blastz') {

      $alnrunnable = Bio::EnsEMBL::Pipeline::Runnable::Blastz->new(-query => $seq1,
								   -database => $seq2,
								   -options => $self->analysis->parameters
								  );
    }
    
    $self->runnable($alnrunnable);
}

=head2 run

    Title   :   run
    Usage   :   $self->run()
    Function:   Runs the crossmatch analysis, producing feature pairs
    Returns :   Nothing, but $self->output is filled 
                contains the crossmatch feature pairs.
    Args    :   None

=cut

sub run {
    my ($self) = @_;
    $self->throw("Can't run - no runnable objects") unless $self->runnable;
    $self->runnable->run();
    $self->output($self->runnable->output);
}

=head2 output

    Title   :   output
    Usage   :   $self->output() or $self->output(@)
    Function:   Push an array of Bio::EnsEMBL::FeaturePairs if Args is given and
                Returns all output feature pairs
    Returns :   Array of Bio::EnsEMBL::FeaturePairs
    Args    :   Array of Bio::EnsEMBL::FeaturePairs (optional)

=cut

sub output {
   my ($self,@args) = @_;
   if (@args) {
     push @{$self->{'output'}}, @args;
   }

   if (defined $self->{'output'}) {
     return $self->{'output'};
   } else {
     return ();
   }
}

=head2 write_output

    Title   :   write_output
    Usage   :   $self->write_output()
    Function:   Writes contents of $self->output into $self->db
    Returns :   1
    Args    :   None

=cut

sub write_output {
  my ($self) = @_;

  if (! scalar @{$self->output}) {
      return 1;
  } 

  my $DnaDnaAlignFeatures;
  if ($self->filter) {
    print STDERR "Features before filtering : ",scalar @{$self->output},"\n";
    $DnaDnaAlignFeatures = $self->filter->filter($self->output);
    print STDERR "Features after filtering  : ",scalar @{$DnaDnaAlignFeatures},"\n";
  } else {
    $DnaDnaAlignFeatures = $self->output;
  }
  my $db = $self->db();
  my $gadb = $db->get_GenomeDBAdaptor();
  $db->get_DnaFragAdaptor();

  my $input_id = $self->input_id();
  my ($sequence_source1,$species_tag1,$dnafrag_type1,$dnafrag_name1,$sequence_source2,$species_tag2,$dnafrag_type2,$dnafrag_name2);

  if ($input_id =~ /^(\S+):(\S+):(\S+):(\S+)::(\S+):(\S+):(\S+):(\S+)$/) {
    ($sequence_source1,$species_tag1,$dnafrag_type1,$dnafrag_name1,$sequence_source2,$species_tag2,$dnafrag_type2,$dnafrag_name2) = ($1,$2,$3,$4,$5,$6,$7,$8);
  } else {
    die "\$input_id should be sequence_source1:species1:dnafrag_type1:dnafrag_name1::sequence_source2:species2:dnafrag_type2:dnafrag_name2 in CrossComparer.pm\n";
  }
  
  # Using $dnafrag_name1 as consensus sequence for the reference align_id

  my $galn = $db->get_GenomicAlignAdaptor();
  my $align_id = $galn->fetch_align_id_by_align_name($dnafrag_name1);

  # Defining the align_row_id

  my $current_align_row_id = 1;

  # Defining an alignement

  my $aln = Bio::EnsEMBL::Compara::GenomicAlign->new();
  
  # Using $dnafrag_name2 as query sequence

  my $dnafrag;

  if ($sequence_source2 eq "ENSEMBL") {
    my $gdb = $gadb->fetch_by_species_tag($species_tag2);
    $dnafrag = Bio::EnsEMBL::Compara::DnaFrag->new();
    $dnafrag->name($dnafrag_name2);
    $dnafrag->genomedb($gdb);
    $dnafrag->type($dnafrag_type2);
  }
  
  my $abs = Bio::EnsEMBL::Compara::AlignBlockSet->new();
  my $hseqname;

  foreach my $f (sort {$a->hseqname cmp $b->hseqname} @{$DnaDnaAlignFeatures}) {
    $hseqname = $f->hseqname unless (defined $hseqname);
    if ($hseqname ne $f->hseqname) {
      $aln->add_AlignBlockSet($current_align_row_id,$abs);
      $abs = Bio::EnsEMBL::Compara::AlignBlockSet->new();
      $hseqname = $f->hseqname;
      $current_align_row_id++;
    }
    if ($sequence_source2 eq "FASTA") {
      my $gdb = $gadb->fetch_by_species_tag($species_tag2);
      $dnafrag = Bio::EnsEMBL::Compara::DnaFrag->new();
      $dnafrag->name($f->hseqname);
      $dnafrag->genomedb($gdb);
      $dnafrag->type($dnafrag_type2);
    }
    my $ab = Bio::EnsEMBL::Compara::AlignBlock->new();

    $ab->align_start($f->start);
    $ab->align_end($f->end);
    $ab->start($f->hstart);
    $ab->end($f->hend);
    if ($f->strand == 1) {
      $ab->strand($f->hstrand);
    } elsif ($f->strand == -1) {
      $ab->strand(- $f->hstrand);
    }
    $ab->score($f->score);
    $ab->perc_id($f->percent_id);
# think here to revert cigar_string if strand==-1 !!
    $ab->cigar_string($f->cigar_string);
    $ab->dnafrag($dnafrag);
    
    $abs->add_AlignBlock($ab);
  }

  # Adding the reference AlignBlockSet to the alignment data

  $aln->add_AlignBlockSet($current_align_row_id,$abs);

  # Storing alignment in the corresponding database with the relevant align_id

  my $galnad = $db->get_GenomicAlignAdaptor();
  $galnad->store($aln,$align_id);

  return 1;
}

=head2 runnable

 Title   : runnable
 Usage   : $obj->runnable($runnable)
 Function: get/set for runnable
 Returns : path to runnable
 Args    : runnable (optional)


=cut

sub runnable {
   my ($self, $runnable) = @_;

   if( defined $runnable ) {
       $self->{'_runnable'} = $runnable;
   }
   return $self->{'_runnable'};
}


=head2 min_score

 Title   : min_score
 Usage   : $obj->min_score($newval)
 Function: Getset for min_score value
 Returns : value of min_score
 Args    : newvalue (optional)


=cut

sub min_score{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'min_score'} = $value;
    }
    return $obj->{'min_score'};

}

=head2 masked

 Title   : masked
 Usage   : $obj->masked($boolean)
 Function: Get/set the  value
 Returns : value of masked
 Args    : boolean, 0 or 1


=cut

sub masked {
   my ($self,@args) = @_;
   if (@args) {
      my ($value) = @args;
      $self->{'masked'} = $value;
    }
    return $self->{'masked'};
}

=head2 filter

 Title   : filter
 Usage   : $obj->filter($boolean)
 Function: Get/set the  value
 Returns : value of filter
 Args    : boolean, 0 or 1


=cut

sub filter {
   my ($self,@args) = @_;
   if (@args) {
      my ($value) = @args;
      $self->{'filter'} = $value;
    }
    return $self->{'filter'};
}

=head2 alnprog

 Title   : alnprog
 Usage   : $obj->alnprog($string)
 Function: Get/set the alnprog value
 Returns : value of alnprog
 Args    : string, could be 'crossmatch', 'exonerate' or 'bl2seq' (optional)


=cut

sub alnprog {
   my ($self,@args) = @_;
   if (@args) {
      my ($value) = @args;
      $self->{'alnprog'} = $value;
    }
    return $self->{'alnprog'};
}

=head2 alntype

 Title   : alntype
 Usage   : $obj->alntype($string)
 Function: Get/set the alntype value
 Returns : value of alntype
 Args    : string, some program could use different algorithm to align sequences.
           In the case of bl2seq, it could be 'blastp', 'blastn', 'blastx', 'tblastn' 
           or 'tblastx' (optional)

=cut

sub alntype {
   my ($self,@args) = @_;
   if (@args) {
      my ($value) = @args;
      $self->{'alntype'} = $value;
    }
    return $self->{'alntype'};
}

=head2 _c1_id

 Title   : _c1_id
 Usage   : $obj->_c1_id($newval)
 Function: Getset for _c1_id value
 Returns : value of _c1_id
 Args    : newvalue (optional)


=cut

sub _c1_id{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'_c1_id'} = $value;
    }
    return $obj->{'_c1_id'};

}

=head2 _c2_id

 Title   : _c2_id
 Usage   : $obj->_c2_id($newval)
 Function: Getset for _c2_id value
 Returns : value of _c2_id
 Args    : newvalue (optional)


=cut

sub _c2_id{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'_c2_id'} = $value;
    }
    return $obj->{'_c2_id'};

}

1;
