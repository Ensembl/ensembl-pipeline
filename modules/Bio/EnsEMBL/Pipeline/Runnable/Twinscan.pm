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
  
Bio::EnsEMBL::Pipeline::Runnable::Twinscan - uses twinscan to make gene predictions

=head1 SYNOPSIS

    # Initialise the Twinscan module  
    my $twinscan = new  Bio::EnsEMBL::Pipeline::Runnable::Twinscan  ('genomic' => $genomic,
								     'conseq'  => $conseq);


   $twinscan->run;
    
   my @genes = $twinscan->output

=head1 DESCRIPTION


=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Pipeline::Runnable::Twinscan;

use strict;  
use vars   qw(@ISA);

# Object preamble - inherits from Bio::Root::RootI

use Bio::Root::RootI;
use Bio::EnsEMBL::Analysis::Programs qw(genewise); 
use Bio::EnsEMBL::SeqFeature;
use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::Pipeline::RunnableI;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI);

sub new {
  my($class,@args) = @_;
  my $self = $class->SUPER::new(@args);
  
  my ($genomic, $conseq) = 
      $self->_rearrange([qw(GENOMIC CONSEQ)], @args);
  
  $self->genomic($genomic) || $self->throw("No genomic sequence entered for blastwise");
  $self->conservation_sequence($conseq) || $self->throw("No protein sequence entered for blastwise");

  return $self;
}

# RunnableI methods

sub run {
    my ($self) = @_;
    $self->set_environment;
    # assumes twinscan is in PATH
    my $sps_parf = "$ENV{TWINSCAN}/misc/sps3sg2x.smat"; # 2x shotgun
    my $genscan_parf = "$ENV{TWINSCAN}/misc/genscan.smat";
    my $seqfile = "/tmp/twin.$$.seq";
#    my $seqfile = "/scratch3/ensembl/vac/twinscan/AC008898/AC008898.00005.seq";
    
    # have to do it this way; it doesn't handle bioperl seqs properly - wants the sequence all on one line I *think*
    open (OUTSEQ, ">$seqfile");
    print OUTSEQ ">" . $self->genomic->id . "\n";
    print OUTSEQ $self->genomic->seq . "\n";
    close OUTSEQ;


    my $command = "twinscan -g $genscan_parf -m SPS -c $sps_parf $seqfile " . $self->conservation_sequence;
    open(TWIN, "$command |") or $self->throw("error running twinscan: $command\n");

    while(<TWIN>){
      print $_;
    }

    close TWIN;

}

sub output {
    my ($self) = @_;

    return $self->eachGene;
#    return $self->eachExon;
}

# Now the twinscan stuff

sub set_environment {
    my ($self) = @_;

#    $ENV{TWINSCAN} = '/usr/local/ensembl/data/twinscan/';
    $ENV{TWINSCAN} = '/scratch3/ensembl/vac/twinscan/';

#    if (! -d '/usr/local/ensembl/data/twinscan/') {
    if (! -d '/scratch3/ensembl/vac/twinscan/') {
#	$self->throw("No TWINSCAN /usr/local/ensembl/data/twinscan/");
	$self->throw("No TWINSCAN /scratch3/ensembl/vac/twinscan/");
    }
}

# These all set/get or initializing methods

sub reverse {
    my ($self,$arg) = @_;

    if (defined($arg)) {
	$self->{'_reverse'} = $arg;
    }
    return $self->{'_reverse'};
}

sub genomic {
    my ($self,$arg) = @_;

#    print ("Gen $arg\n");
    if (defined($arg)) {
	$self->throw("Genomic sequence input is not a Bio::SeqI or Bio::Seq or Bio::PrimarySeqI") unless
	    ($arg->isa("Bio::SeqI") || 
	     $arg->isa("Bio::Seq")  || 
	     $arg->isa("Bio::PrimarySeqI"));
		
	
	$self->{'_genomic'} = $arg;
    }
    return $self->{'_genomic'};
}

sub conservation_sequence {
  my ($self,$conseq) = @_;
  
  if (defined($conseq)) {
    $self->throw("$conseq does not exist") unless -e $conseq;
    $self->{'_conseq'} = $conseq;
  }
#  $self->{'_conseq'} = '/scratch3/ensembl/vac/twinscan/AC008898/conseq.file';
  return $self->{'_conseq'};
}

sub addGene {

  my ($self,$arg) = @_;
  if (!defined($self->{'_output'})) {
    $self->{'_output'} = [];
  }
  push(@{$self->{'_output'}},$arg);
  
}

sub eachGene {
    my ($self) = @_;

    if (!defined($self->{'_output'})) {
	$self->{'_output'} = [];
    }

    return @{$self->{'_output'}};
}

sub addExon {
    my ($self,$arg) = @_;

    if (!defined($self->{'_output'})) {
	$self->{'_output'} = [];
    }
    push(@{$self->{'_output'}},$arg);

}

sub eachExon {
    my ($self) = @_;

    if (!defined($self->{'_output'})) {
	$self->{'_output'} = [];
    }

    return @{$self->{'_output'}};
}

1;
