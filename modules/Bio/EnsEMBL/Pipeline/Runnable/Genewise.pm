package Bio::EnsEMBL::Pipeline::Runnable::Genewise;

=pod

=head1 NAME
  
Bio::EnsEMBL::Pipeline::Runnable::GeneWise - aligns protein hits to genomic sequence

=head1 SYNOPSIS

    # Initialise the genewise module  
    my $genewise = new  Bio::EnsEMBL::Pipeline::Runnable::Genewise  ('genomic'  => $genomic,
								       'protein'  => $protein,
								       'memory'   => $memory);


   $genewise->run;
    
   my @genes = $blastwise->output

=head1 DESCRIPTION


=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut

use strict;  
use vars   qw(@ISA);

# Object preamble - inherits from Bio::Root::Object

use Bio::Root::Object;
use Bio::EnsEMBL::Analysis::Programs qw(genewise); 
use Bio::EnsEMBL::SeqFeature;
use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::Pipeline::RunnableI;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI Bio::Root::Object);

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,@args) = @_;

  my $make = {};
  bless $make, ref($self);

  my ($genomic, $protein, $memory,$reverse,$endbias) = 
      $self->_rearrange(['GENOMIC','PROTEIN','MEMORY','REVERSE','ENDBIAS'], @args);
  
  print $genomic . "\n";

  $self->genomic($genomic) || $self->throw("No genomic sequence entered for blastwise");
  $self->protein($protein) || $self->throw("No protein sequence entered for blastwise");

  $self->reverse($reverse)   if (defined($reverse));             
  $self->endbias($endbias)   if (defined($endbias));             
  $self->memory ($memory)    if (defined($memory));

  return $make; # success - we hope!
}

# RunnableI methods

sub run {
    my ($self) = @_;

    $self->_wise_init;
    $self->_align_protein;
}

sub output {
    my ($self) = @_;

    return $self->eachGene;
}

# Now the blastwise stuff

sub set_environment {
    my ($self) = @_;

    $ENV{WISECONFIGDIR} = '/usr/local/ensembl/data/wisecfg/';

    if (! -d '/usr/local/ensembl/data/wisecfg/') {
	$self->throw("No WISECONFIGDIR /usr/local/ensembl/data/wisecfg/");
    }
}

sub _wise_init {
  my($self,$genefile) = @_;

  print("Initializing genewise\n");

  $self->set_environment;

}


=pod

=head2 _align_protein

  Title   : _align_protein
  Usage   : $bw->align_protein($pro,$start,$end)
  Function: aligns a protein to a genomic sequence between the genomic coords $start,$end
  Returns : nothing
  Args    : Bio:Seq,int,int

=cut


sub _align_protein {
  my ($self) = @_;
  my $memory   = $self->memory;
  
  my $gwfile   = "/tmp/gw." . $$ . ".out";
  my $genfile  = "/tmp/gw." . $$ . "." . $self->protein->id . ".gen.fa";
  my $protfile = "/tmp/gw." . $$ . "." . $self->protein->id . ".pro.fa";
  
  my $genio  = new Bio::SeqIO(-file   => ">$genfile",
			      '-format' => 'fasta');
  
  $genio->write_seq($self->genomic);
  $genio = undef;
  
  my $proio  = new Bio::SeqIO(-file   => ">$protfile",
			      '-format' => 'fasta');
  
  $proio->write_seq($self->protein);
  $proio = undef;

# my $command = "genewise $protfile $genfile -gff -alb -kbyte $memory -ext 1 -gap 8 -subs 0.0000001 -quiet";
#  my $command = "genewise $protfile $genfile -gff -alb -kbyte $memory -ext 2 -gap 12 -subs 0.0000001 -quiet -init endbias -splice flat";
  my $command = "genewise $protfile $genfile -gff -alb -kbyte $memory -ext 2 -gap 12 -subs 0.0000001 -quiet";
    
  if ($self->endbias == 1) {
    $command .= " -init endbias -splice flat ";
  }
  
  if ($self->reverse == 1) {
    $command .= " -trev ";
  }
  
  print STDERR "Command is $command\n";
  
  open(GW, "$command |");
  
  my @gff;
  my $sf;
  my $count = 1;
  
  my @gw_exons;
  my $curr_exon;
  
  my $genename;
  # need to parse both gff and alb output
  # making assumption of only 1 gene prediction ...
  # gff o/p always comes before alb o/p
 GW: while (<GW>) {
    chomp;
    my @f = split;
    
    #      print STDERR "Line is " . $_ . "\n";
    # may be foolish, but I *think* we only care about exon intron boundaries ... at least at this stage.
    next GW unless (defined($f[2]) &&($f[2] eq '"MATCH_STATE"' || $f[2] eq '"INTRON_STATE"' || $f[2] eq "cds"));
    
    if($f[2] eq 'cds'){
      # this is gff
      
      if ($f[8] ne $genename) {
	print STDERR "Found new gene\n";
	
	$genename = $f[8];
	$sf = new Bio::EnsEMBL::SeqFeature;
	$sf->seqname   ($f[8]);
	$sf->id        ($self->protein->id);
	$sf->source_tag('genewise');
	$sf->primary_tag('similarity');
	
	push(@gff,$sf);
	
	$count = 1;
      }
      
      my $subsf = new Bio::EnsEMBL::SeqFeature;
      
      my $strand = $f[6];
      
      if ($strand eq "+") {
	$subsf->strand(1);
      } else {
	$subsf->strand(-1);
      }
      
      $subsf->seqname($f[8]);
      
      if ($f[3] < $f[4]) {
	$subsf->start  ($f[3]);
	$subsf->end    ($f[4]);
      } else {
	$subsf->start  ($f[4]);
	$subsf->end    ($f[3]);
      }
      print STDERR "Found exon " . $sf->start . " " . $sf->end . "\n";
      
      my $phase = $f[7];
      
      if ($phase == 2) { 
	$phase = 1;
      }	elsif ($phase == 1) { 
	$phase = 2;
      }
      $subsf->{_phase} = $phase;
      $subsf->id ($self->protein->id . ".$count");
      
      print STDERR "SF " . $subsf->start . " " . $subsf->end . " " . $subsf->strand . " " . $subsf->{_phase} . "\n";
      $count++;
      
      $sf->add_sub_SeqFeature($subsf,'EXPAND');
      
    }
    else {
      # this is alb
      $f[1] =~ /\[\d+:(\d+)/;
      my $coord = $1;
      $coord++;
      
      if($f[2] ne '"MATCH_STATE"'){
	
	if($f[4] eq '"CENTRAL_INTRON"'){
	  $coord++;
	  # start a new exon
	  my $ef = new Bio::EnsEMBL::SeqFeature;
	  $ef->seqname   ($self->protein->id);
	  $ef->id        ('exon_pred');
	  $ef->source_tag('genewise');
	  $ef->primary_tag('similarity');
	  $ef->start($coord);
	  $ef->end($coord);
	  
	  $curr_exon = $ef;
	  push(@gw_exons,$curr_exon);
	}
	
      }
      else {
	if(!defined($curr_exon)){
	  # start a new exon
	  my $ef = new Bio::EnsEMBL::SeqFeature;
	  $ef->seqname   ($self->protein->id);
	  $ef->id        ('exon_pred');
	  $ef->source_tag('genewise');
	  $ef->primary_tag('similarity');
	  $curr_exon = $ef;
	  $curr_exon->start(0);
	  push(@gw_exons,$curr_exon);
	}
	
	# extend existing exon boundary
	$curr_exon->end($coord);
	# set start if we need to - this should only happen for the very first exon
	if(!$curr_exon->start) {
	  print STDERR "rigging start\n";
	  $curr_exon->start($coord);
	  
	}
	
      }
    }
    
  } 
  
  close(GW);
  
  # Now convert into feature pairs
  # makes hstart & hend for the feature pair in cDNA coordinates?
  if(scalar(@gff) > 1){
    $self->throw("rats; more than one gene predicted");
  }
  
  foreach my $gene (@gff) {
    my @sub = $gene->sub_SeqFeature;
    
    if (scalar(@sub) != scalar(@gw_exons)) {
      die "unequal number of features from gff and alb sections\n"; 
    }
    
    my $index=0;
    my @fp;
    for ($index=0; $index <= $#sub; $index++){
      my $feat_pair = new Bio::EnsEMBL::FeaturePair(-feature1 => $sub[$index],
						    -feature2 => $gw_exons[$index]);
      push(@fp,$feat_pair);
    }     
    
    # now sort them
    if ($#sub >=0) {
      if ($sub[0]->strand == 1) {
	@fp = sort {$a->start <=> $b->start} @fp;
      } else {
	@fp = sort {$b->start <=> $a->start} @fp;
      }
    }
    
    foreach my $pair(@fp){
      $self->addGene($pair);
    }
  }
  
      unlink $genfile;
      unlink $protfile;
      unlink $gwfile;
}


# These all set/get or initializing methods

sub reverse {
    my ($self,$arg) = @_;

    if (defined($arg)) {
	$self->{_reverse} = $arg;
    }
    return $self->{_reverse};
}

sub endbias {
    my ($self,$arg) = @_;

    if (defined($arg)) {
	$self->{_endbias} = $arg;
    }

    if (!defined($self->{_endbias})) {
      $self->{_endbias} = 0;
    }

    return $self->{_endbias};
}

sub memory {
    my ($self,$arg) = @_;

    if (defined($arg)) {
	$self->{_memory} = $arg;
    }

    return $self->{_memory} || 100000;
}

sub genomic {
    my ($self,$arg) = @_;

#    print ("Gen $arg\n");
    if (defined($arg)) {
	$self->throw("Genomic sequence input is not a Bio::SeqI or Bio::Seq or Bio::PrimarySeqI") unless
	    ($arg->isa("Bio::SeqI") || 
	     $arg->isa("Bio::Seq")  || 
	     $arg->isa("Bio::PrimarySeqI"));
		
	
	$self->{_genomic} = $arg;
    }
    return $self->{_genomic};
}

sub protein {
    my ($self,$arg) = @_;

    if (defined($arg)) {
	$self->throw("[$arg] is not a Bio::SeqI or Bio::Seq or Bio::PrimarySeqI") unless
	    ($arg->isa("Bio::SeqI") || 
	     $arg->isa("Bio::Seq")  || 
	     $arg->isa("Bio::PrimarySeqI"));

	$self->{_protein} = $arg;


    }
    return $self->{_protein};
}

sub addGene {
    my ($self,$arg) = @_;

    if (!defined($self->{_output})) {
	$self->{_output} = [];
    }
    push(@{$self->{_output}},$arg);

}

sub eachGene {
    my ($self) = @_;

    if (!defined($self->{_output})) {
	$self->{_output} = [];
    }

    return @{$self->{_output}};
}


1;







