=pod

=head1 NAME
  
Bio::EnsEMBL::Pipeline::Runnable::GeneWise - aligns protein hits to genomic sequence

=head1 SYNOPSIS

    # Initialise the genewise module  
    my $genewise = new  Bio::EnsEMBL::Pipeline::Runnable::Genewise  ('genomic'  => $genomic,
																																		 'protein'  => $protein,
																																		 'memory'   => $memory);


   $genewise->run;
    
   my @genes = $genewise->output

=head1 DESCRIPTION


=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Pipeline::Runnable::Genewise;

use strict;  
use vars   qw(@ISA);

# Object preamble

use Bio::EnsEMBL::Root;
use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::Pipeline::RunnableI;
use Bio::SeqIO;
use Bio::EnsEMBL::Pipeline::Config::GeneBuild::Genewise;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::DnaPepAlignFeature

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI);

sub new {
  my($class,@args) = @_;
  my $self = $class->SUPER::new(@args);
  
  my ($genomic, $protein, $slice, $memory,$reverse,$endbias,$genewise,$gap, 
      $ext, $subs, $options) = $self->_rearrange([qw(GENOMIC 
						     PROTEIN 
						     SLICE
						     MEMORY 
						     REVERSE 
						     ENDBIAS 
						     GENEWISE 
						     GAP 
						     EXTENSION 
						     SUBS 
						     OPTIONS)], @args);
  
  $genewise = $GB_GENEWISE_EXE unless defined($genewise);

  $self->genewise($self->find_executable($genewise));

  $self->genomic($genomic) || $self->throw("No genomic sequence entered for blastwise");
  $self->protein($protein) || $self->throw("No protein sequence entered for blastwise");
  $self->slice($slice)     || $self->throw("No slice sequence entered for blastwise");

  $self->reverse($reverse)   if ($reverse);             
  $self->endbias($endbias)   if ($endbias);
  $memory = $GB_GENEWISE_MEMORY if(!$memory);
  $self->memory ($memory);
  $gap = $GB_GENEWISE_GAP unless($gap);
  $self->gap($gap);
  $ext = $GB_GENEWISE_EXTENSION unless($ext);
  $self->extension($ext);
  $subs = $GB_GENEWISE_SUBS unless($subs);
  $self->subs($subs);
  $options = $GB_GENEWISE_OPTIONS unless($options);
  $self->options($options);
  #print STDERR "Have genomic of length ".$self->genomic->length."\n";
  return $self;
}



sub run {
    my ($self) = @_;

    my($p, $f, $l) = caller;
    $self->check_environment;
    eval{
      $self->align_protein;
    };
    if($@){
      $self->warn("Genewise run failed caller = $f:$l $@");
    }
    return 1;
}



=pod

=head2 align_protein

  Title   : align_protein
  Usage   : $bw->align_protein
  Function: aligns a protein to a genomic sequence between the genomic coords $start,$end
  Returns : nothing
  Args    : Bio:Seq,int,int

=cut


sub align_protein {
  my ($self) = @_;

  my $genewise = $self->genewise;
  my $memory   = $self->memory;
  my $gap = $self->gap;
  my $ext = $self->extension;
  my $subs = $self->subs;
  my $options = $self->options;
  my $results_file  = $self->get_tmp_file('/tmp/',$self->protein->id,
                                          ".gen.fa");

  my $genfile = $self->write_sequence_to_file($self->genomic);
  my $pepfile = $self->write_sequence_to_file($self->protein);

  my $command = "$genewise $pepfile $genfile -genesf -kbyte $memory -ext $ext -gap $gap -subs $subs $options";


  if ($self->endbias == 1) {
    $command =~ s/-init endbias//;
    $command =~ s/-splice flat//;

    $command .= " -init endbias -splice flat ";

  }

  if (($self->reverse) && $self->reverse == 1) {
    $command .= " -trev ";
  }

  print "Running genewise with $command\n";

  open(GW, "$command | ") or $self->throw("error piping to genewise: $!\n");

  my @genesf_exons;
  my $curr_exon;
  my $curr_gene;


 GENESF: 
  while (<GW>) {
    chomp;
    my @f = split;
    next unless (defined $f[0] && ($f[0] eq 'Gene' || $f[0] eq 'Exon' || $f[0] eq 'Supporting'));
    $self->parse_genewise_output_line(\@f, \$curr_gene, \$curr_exon, \@genesf_exons);

  }

  close(GW) or $self->throw("Error running genewise with command line ".$command."\n $!");
  unlink $genfile;
  unlink $pepfile;

  # now make gene and transcript from array of exons
  my $transcript = $self->make_transcript(\@genesf_exons);

 if(defined $transcript){
    my $gene = new Bio::EnsEMBL::Gene;
    $gene->type("genewise");
    $gene->add_Transcript($transcript);

    $self->addGene($gene);

  }
  else{
   print "No valid transcript made, cannot make gene\n"; 
  }
}

=head2 parse_genewise_output_line

  Arg [1]   : $columns_ref, reference to array of strings
  Arg [2]   : $curr_gene_ref, reference to scalar
  Arg [3]   : $curr_exon_ref, reference to Bio::EnsEMBL::Feature
  Arg [4]   : $exons_ref, reference to array of Bio::EnsEMBL::Feature
  Function  : Parses genewise output lines and adds exons/supporting features to exons_ref
  Returntype: 
  Exceptions: 
  Caller    :
  Example   :

=cut

sub parse_genewise_output_line{
  my ($self, $columns_ref, $curr_gene_ref, $curr_exon_ref, $exons_ref) = @_;

  if($columns_ref->[0] eq 'Gene'){

    # flag a frameshift - ultimately we will do something clever here but for now ...
    # frameshift is here defined as two or more genes produced by Genewise from a single protein

    if (/^(Gene\s+\d+)$/){
      if(!defined($$curr_gene_ref)){
	$$curr_gene_ref = $1;
      }	elsif ($1 ne $$curr_gene_ref) {
	$self->warn("frameshift!\n");
      }
    }
  }

  elsif($columns_ref->[0] eq 'Exon'){
    my $start  = $columns_ref->[1];
    my $end    = $columns_ref->[2];
    my $phase  = $columns_ref->[4];
    my $strand = 1;

    if($columns_ref->[1] > $columns_ref->[2]){
      $strand = -1;
      $start  = $columns_ref->[2];
      $end    = $columns_ref->[1];
    }

    my $exon_length = $end - $start + 1;

    # end phase is the number of bases at the end of the exon which do not
    # fall in a codon and it coincides with the phase of the following exon.

    my $end_phase   = ( $exon_length + $phase ) %3;

    $$curr_exon_ref = new Bio::EnsEMBL::Exon;
    $$curr_exon_ref->start    ($start);
    $$curr_exon_ref->end      ($end);
    $$curr_exon_ref->strand   ($strand);
    $$curr_exon_ref->phase    ($phase);
    $$curr_exon_ref->end_phase($end_phase);

    $$curr_exon_ref->slice($self->slice);

    push(@$exons_ref, $$curr_exon_ref);
  }

  elsif($columns_ref->[0] eq 'Supporting') {

    my $gstart = $columns_ref->[1];
    my $gend   = $columns_ref->[2];
    my $pstart = $columns_ref->[3];
    my $pend   = $columns_ref->[4];

    my $strand = 1;

    if ($gstart > $gend){
      $gstart = $columns_ref->[2];
      $gend   = $columns_ref->[1];
      $strand = -1;
    }

    if($strand != $$curr_exon_ref->strand){
      $self->warn("incompatible strands between exon and supporting feature - cannot add suppfeat\n");
      return;
    }

    if($pstart > $pend){
      $self->warn("Protein start greater than end! Skipping this suppfeat\n");
      return;
    }

    my $fp = new Bio::EnsEMBL::FeaturePair();
    my @features;
    $fp->start   ($gstart);
    $fp->end     ($gend);
    $fp->strand  ($strand);
    $fp->seqname ($self->genomic->id);
    $fp->hseqname($self->protein->id);
    $fp->hstart  ($pstart);
    $fp->hend    ($pend);
    $fp->hstrand (1);
    push(@features, $fp);

    my $align = new Bio::EnsEMBL::DnaPepAlignFeature(-features => \@features);
    $align->seqname($self->slice->seq_region_name);
    $align->slice($self->slice);
    $align->score(100);

    $$curr_exon_ref->add_supporting_features($align);
  }
}

=head2 make_transcript

  Arg [1]   : $exons_ref, reference to array of Bio::EnsEMBL::Feature
  Function  : Turns array of exons into transcript & validates it.
  Returntype: Bio::EnsEMBL::Transcript
  Exceptions: 
  Caller    :
  Example   :

=cut

sub make_transcript{
  my ($self, $exonsref) = @_;
  my @exons = @$exonsref;

  my $transcript   = Bio::EnsEMBL::Transcript->new;
  my $translation  = Bio::EnsEMBL::Translation->new;
  $transcript->translation($translation);

  if ($#exons < 0) {
    print STDERR "Odd.  No exons found\n";
    return undef;
  }

  else {

    if ($exons[0]->strand == -1) {
      @exons = sort {$b->start <=> $a->start} @exons;
    } else {
      @exons = sort {$a->start <=> $b->start} @exons;
    }

    $translation->start_Exon($exons[0]);
    $translation->end_Exon  ($exons[$#exons]);

    # phase is relative to the 5' end of the transcript (start translation)
    if ($exons[0]->phase == 0) {
      $translation->start(1);
    } elsif ($exons[0]->phase == 1) {
      $translation->start(3);
    } elsif ($exons[0]->phase == 2) {
      $translation->start(2);
    }

    $translation->end  ($exons[$#exons]->end - $exons[$#exons]->start + 1);
    foreach my $exon(@exons){
      $transcript->add_Exon($exon);
    }

  }

  $transcript->slice($self->slice);
  return $transcript;
}

# These all set/get or initializing methods

sub reverse {
    my ($self,$arg) = @_;

    if (defined($arg)) {
	$self->{'_reverse'} = $arg;
    }
    return $self->{'_reverse'};
}

sub endbias {
    my ($self,$arg) = @_;

    if (defined($arg)) {
	$self->{'_endbias'} = $arg;
    }

    if (!defined($self->{'_endbias'})) {
      $self->{'_endbias'} = 0;
    }

    return $self->{'_endbias'};
}

sub memory {
    my ($self,$arg) = @_;

    if (defined($arg)) {
			$self->{'_memory'} = $arg;
    }

    return $self->{'_memory'} || 100000;
}

sub gap {
    my ($self,$arg) = @_;

    if(!$self->{'_gap'}){
      $self->{'_gap'} = undef;
    }
    if (defined($arg)) {
      $self->{'_gap'} = $arg;
    }

    return $self->{'_gap'} || 12;
}

sub extension {
    my ($self,$arg) = @_;

    if(!$self->{'_ext'}){
      $self->{'_ext'} = undef;
    }
    if (defined($arg)) {
      $self->{'_ext'} = $arg;
    }

    return $self->{'_ext'} || 2;
}

sub subs{
    my ($self,$arg) = @_;

    if(!$self->{'_subs'}){
      $self->{'_subs'} = undef;
    }
    if (defined($arg)) {
      $self->{'_subs'} = $arg;
    }

    return $self->{'_subs'} || 0.0000001;
}





sub genomic {
    my ($self,$arg) = @_;

    if (defined($arg)) {
			$self->throw("Genomic sequence input is not a Bio::SeqI or Bio::Seq or Bio::PrimarySeqI") unless
				($arg->isa("Bio::SeqI") || 
				 $arg->isa("Bio::Seq")  || 
				 $arg->isa("Bio::PrimarySeqI"));
			
			$self->{'_genomic'} = $arg;
    }
    return $self->{'_genomic'};
}

sub protein {
    my ($self,$arg) = @_;

    if (defined($arg)) {
			$self->throw("[$arg] is not a Bio::SeqI or Bio::Seq or Bio::PrimarySeqI") unless
				($arg->isa("Bio::SeqI") || 
				 $arg->isa("Bio::Seq")  || 
				 $arg->isa("Bio::PrimarySeqI"));
			
			$self->{'_protein'} = $arg;
    }
    return $self->{'_protein'};
}

sub slice {
    my ($self,$arg) = @_;

    if (defined($arg)) {
			$self->throw("slice sequence input is not a Bio::EnsEMBL::Slice") unless
				($arg->isa("Bio::EnsEMBL::Slice"));
			
			$self->{'_slice'} = $arg;
    }
    return $self->{'_slice'};
}

sub addGene {
  my ($self,$arg) = @_;

  if (!defined($self->{'_output'})) {
    $self->{'_output'} = [];
  }
  #print STDERR "Adding ".$arg." ".$arg->start." ".$arg->end."\n";
  push(@{$self->{'_output'}},$arg);

}

sub genewise {
  my ($self,$arg) = @_;

  if (defined($arg)) {
    $self->{_genewise} = $arg;
  }
  return $self->{_genewise};
}

sub check_environment {
  my($self,$genefile) = @_;

  if (! -d $ENV{WISECONFIGDIR}) {
    $self->throw("No WISECONFIGDIR ["  . $ENV{WISECONFIGDIR} . "]");
  }

}

1;
