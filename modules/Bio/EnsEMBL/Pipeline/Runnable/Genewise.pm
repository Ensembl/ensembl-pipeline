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
use Bio::EnsEMBL::SeqFeature;
use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::Pipeline::RunnableI;
use Bio::SeqIO;
use Bio::EnsEMBL::Pipeline::Config::GeneBuild::Genewise;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI);

sub new {
  my($class,@args) = @_;
  my $self = $class->SUPER::new(@args);
  
  my ($genomic, $protein, $memory,$reverse,$endbias,$genewise,$gap, 
      $ext, $subs, $options) = $self->_rearrange([qw(GENOMIC 
						     PROTEIN 
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
  #my $genfile  = $self->get_tmp_file('/tmp/',$self->protein->id,".gen.fa");
  #my $pepfile  = $self->get_tmp_file('/tmp/',$self->protein->id,".pro.fa");

  my $genfile = $self->write_sequence_to_file($self->genomic);
  my $pepfile = $self->write_sequence_to_file($self->protein);

  my $command = "$genewise $pepfile $genfile -genesf -kbyte $memory -ext $ext -gap $gap -subs $subs $options";

  #print STDERR $command."\n";
  if ($self->endbias == 1) {
    $command .= " -init endbias -splice flat ";
  }

  if (defined($self->reverse) && $self->reverse == 1) {
    $command .= " -trev ";
  }
  
  #print STDERR "GENEWISE running on ".$self->protein->id." with command ".$command."\n";
  my @time = times;
  #print STDERR "BEFORE GENEWISE @time\n"; 
  open(GW, "$command |") or $self->throw("error piping to genewise: $!\n");

  my @genesf_exons;
  my $curr_exon;
  my $curr_gene;

  # making assumption of only 1 gene prediction ... this will change once we start using hmms

 GENESF: 
  while (<GW>) {
   
    chomp;
    my @f = split;
   
    next unless (defined $f[0] && ($f[0] eq 'Gene' || $f[0] eq 'Exon' || $f[0] eq 'Supporting'));
    
    if($f[0] eq 'Gene'){

      # flag a frameshift - ultimately we will do something clever here but for now ...
      # frameshift is here defined as two or more genes produced by Genewise from a single protein

      if (/^(Gene\s+\d+)$/){
				if(!defined($curr_gene)){
					$curr_gene = $1;
				}	elsif ($1 ne $curr_gene) {
					$self->warn("frameshift!\n");
				}
      }
    }
    
    elsif($f[0] eq 'Exon'){
      my $start  = $f[1];
      my $end    = $f[2];
      my $phase  = $f[4];
      my $strand = 1;

      if($f[1] > $f[2]){
				$strand = -1;
				$start  = $f[2];
				$end    = $f[1];
      }

      my $exon_length = $end - $start + 1;

      # end phase is the number of bases at the end of the exon which do not 
      # fall in a codon and it coincides with the phase of the following exon.

      my $end_phase   = ( $exon_length + $phase ) %3;

      $curr_exon = new Bio::EnsEMBL::SeqFeature;
      $curr_exon->seqname  ($self->genomic->id);
      $curr_exon->id       ($self->protein->id);
      $curr_exon->start    ($start);
      $curr_exon->end      ($end);
      $curr_exon->strand   ($strand);
      $curr_exon->phase    ($phase);
      $curr_exon->end_phase($end_phase);

      $self->addExon($curr_exon);

      push(@genesf_exons, $curr_exon);
    
    }  elsif($f[0] eq 'Supporting') {

      my $gstart = $f[1];
      my $gend   = $f[2];
      my $pstart = $f[3];
      my $pend   = $f[4];

      my $strand = 1;

      if ($gstart > $gend){
				$gstart = $f[2];
				$gend   = $f[1];
				$strand = -1;
      }
      
      if($strand != $curr_exon->strand){
				$self->warn("incompatible strands between exon and supporting feature - cannot add suppfeat\n");
				next GENESF;
      }

      if($pstart > $pend){
				$self->warn("Protein start greater than end! Skipping this suppfeat\n");
				next GENESF;
      }
      
      my $fp = new Bio::EnsEMBL::FeaturePair();

      $fp->start   ($gstart);
      $fp->end     ($gend);
      $fp->strand  ($strand);
      $fp->seqname ($self->genomic->id);
      $fp->hseqname($self->protein->id);
      $fp->hstart  ($pstart);
      $fp->hend    ($pend);
      $fp->hstrand (1);

      $curr_exon->add_sub_SeqFeature($fp,'');
    }
  } 
  
  close(GW) or $self->throw("Error running genewise with command line ".$command."\n $!");
  @time = times;
  #print STDERR "AFTER GENEWISE @time\n"; 
  unlink $genfile;
  unlink $pepfile;
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

sub addExon {
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
