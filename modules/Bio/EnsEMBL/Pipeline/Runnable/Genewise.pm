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

#    return $self->eachGene;
    return $self->eachExon;
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
  
  my @gff_exons;
  my $curr_gff;
  my $count = 1;
  
  my @alb_exons;
  my $curr_exon;
  my $prev_state = "";  
  my $albphase = "";
  
  # need to parse both gff and alb output
  # making assumption of only 1 gene prediction ...
  # gff o/p always comes before alb o/p
 GW: while (<GW>) {
    chomp;
    my @f = split;
    
    #      print STDERR "Line is " . $_ . "\n";
    # at this stage, just worry about exon boundaries
    next GW unless (defined($f[2]) &&( $f[2] eq '"MATCH_STATE"'  || $f[2] eq '"INTRON_STATE"' ||
				       $f[2] eq '"INSERT_STATE"' || $f[2] eq '"DELETE_STATE"' ||
				       $f[2] eq "cds"));
    
    if($f[2] eq 'cds'){
      # this is gff
      # make a new "exon"
      $curr_gff = new Bio::EnsEMBL::SeqFeature;
      $curr_gff->seqname   ($f[8]);
      $curr_gff->id        ($self->protein->id);
      $curr_gff->source_tag('genewise');
      $curr_gff->primary_tag('similarity');
      
      push(@gff_exons,$curr_gff);
      
      $count = 1;
      
      
      my $strand = $f[6];
      
      if ($strand eq "+") {
	$curr_gff->strand(1);
      } else {
	$curr_gff->strand(-1);
      }
      
      $curr_gff->seqname($f[8]);
      
      if ($f[3] < $f[4]) {
	$curr_gff->start  ($f[3]);
	$curr_gff->end    ($f[4]);
      } else {
	$curr_gff->start  ($f[4]);
	$curr_gff->end    ($f[3]);
      }
      
      my $phase = $f[7];
      
      if ($phase == 2) { 
	$phase = 1;
      }	elsif ($phase == 1) { 
	$phase = 2;
      }
      
      $curr_gff->{_phase} = $phase;
      $curr_gff->id ($self->protein->id . ".$count"); 
      
      $count++;
    }
    else {
      # this is alb
      $f[1] =~ /\[\d+:(\d+)/;
      my $coord = $1;
      $coord++;

      if ($f[4] eq '"SEQUENCE_INSERTION"' || $f[4] eq '"SEQUENCE_DELETION"') {

	# and make sure we don't increment it
	$albphase = '"3SS_PHASE_0"';
	# start a new exon
	$curr_exon = undef;
      }

      
      if ($f[2] eq '"MATCH_STATE"') {
	if ($f[4] =~ /PHASE/) {
	  $albphase = $f[4];
	}
	
	elsif ( $f[4] eq '"CODON"' || $f[4] eq '"3SS_PHASE_2"' ) {
	  # first, make sure there's a $curr_exon
	  if (!defined($curr_exon)) {
	    # start a new exon
	    my $ef = new Bio::EnsEMBL::SeqFeature;
	    $ef->seqname   ($self->protein->id);
	    $ef->id        ('exon_pred');
	    $ef->source_tag('genewise');
	    $ef->primary_tag('similarity');
	    $ef->start($coord);
	    $ef->end($coord);
	    
	    $curr_exon = $ef;
	    push(@alb_exons,$curr_exon);
	  }
	  
	  $curr_exon->end($coord);
	  
	}
	$prev_state = 'MATCH';
	
      }
      
      elsif ( $f[2] eq '"INTRON_STATE"' && $f[4] eq '"CENTRAL_INTRON"' ) {
	$prev_state = 'INTRON';
	$curr_exon = undef;
      }
      
    }
  } 
  
  close(GW);
  
  # Now convert into feature pairs
  # makes hstart & hend for the feature pair in cDNA coordinates?
  # Now convert into feature pairs
  if(scalar(@gff_exons) != scalar(@alb_exons) ){
    print STDERR "gff: " .scalar(@gff_exons) . "\talb: " . scalar(@alb_exons) .  "\n";
    unlink $genfile;
    unlink $protfile;
    unlink $gwfile;
    $self->throw("rats; exon numbers don't match");
  }
  
  my $index=0;
  my @fp;
  for ($index=0; $index <= $#gff_exons; $index++){
    my $feat_pair = new Bio::EnsEMBL::FeaturePair(-feature1 => $gff_exons[$index],
						  -feature2 => $alb_exons[$index]);
    push(@fp,$feat_pair);
  }     
  
  # now sort them
  if ($#gff_exons >=0) {
    if ($gff_exons[0]->strand == 1) {
      @fp = sort {$a->start <=> $b->start} @fp;
    } else {
      @fp = sort {$b->start <=> $a->start} @fp;
    }
  }
  
  foreach my $pair(@fp){
    $self->addExon($pair);
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

#sub addGene {
#    my ($self,$arg) = @_;

#    if (!defined($self->{_output})) {
#	$self->{_output} = [];
#    }
#    push(@{$self->{_output}},$arg);

#}

#sub eachGene {
#    my ($self) = @_;

#    if (!defined($self->{_output})) {
#	$self->{_output} = [];
#    }

#    return @{$self->{_output}};
#}
sub addExon {
    my ($self,$arg) = @_;

    if (!defined($self->{_output})) {
	$self->{_output} = [];
    }
    push(@{$self->{_output}},$arg);

}

sub eachExon {
    my ($self) = @_;

    if (!defined($self->{_output})) {
	$self->{_output} = [];
    }

    return @{$self->{_output}};
}



1;







