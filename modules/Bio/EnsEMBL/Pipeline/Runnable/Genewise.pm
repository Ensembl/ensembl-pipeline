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

  my ($genomic, $protein, $memory,$reverse) = 
      $self->_rearrange(['GENOMIC','PROTEIN','MEMORY','REVERSE'], @args);
  
  print $genomic . "\n";

  $self->genomic($genomic) || $self->throw("No genomic sequence entered for blastwise");
  $self->protein($protein) || $self->throw("No protein sequence entered for blastwise");

  $self->reverse($reverse)   if (defined($reverse));             
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

    $ENV{WISECONFIGDIR} = '/nfs/disk100/pubseq/wise/wisecfg/';

    if (! -d '/nfs/disk100/pubseq/wise/wisecfg/') {
	$self->throw("No WISECONFIGDIR /nfs/disk100/pubseq/wise/wisecfg/");
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
    my $genfile  = "/tmp/gw." . $$ . ".gen.fa";
    my $protfile = "/tmp/gw." . $$ . ".pro.fa";

    my $genio  = new Bio::SeqIO(-file   => ">$genfile",
				'-format' => 'fasta');

    $genio->write_seq($self->genomic);
    $genio = undef;

    my $proio  = new Bio::SeqIO(-file   => ">$protfile",
				'-format' => 'fasta');

    $proio->write_seq($self->protein);
    $proio = undef;

    my $command = "genewise $protfile $genfile -gff -kbyte $memory -ext 1 -gap 8 -subs 0.0000001 -quiet";

    if ($self->reverse == 1) {
	$command .= " -trev ";
    }

    print "Command is $command\n";
    
    open(GW, "$command |");
    
    my @out;
    my $sf;
    my $count = 1;

    my $genename;
    
    GW: while (<GW>) {
	chomp;
	my @f = split(/\t/,$_);
	
	print STDERR "Line is " . $_ . "\n";
	next GW if ($f[2] ne "cds");

	if ($f[8] ne $genename) {
	    print STDERR "Found new gene\n";
	    
	    $genename = $f[8];
	    $sf = new Bio::EnsEMBL::SeqFeature;
	    $sf->seqname   ($f[8]);
	    $sf->id        ($self->protein->id);
	    $sf->source_tag('genewise');
	    $sf->primary_tag('similarity');

	    push(@out,$sf);

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

    close(GW);

    # Now convert into feature pairs

    foreach my $gene (@out) {
	my @sub = $gene->sub_SeqFeature;

	my $strand = $sub[0]->strand;
	if ($strand == 1) {
	    @sub = sort {$a->start <=> $b->start} @sub;
	} else {
	    @sub = sort {$b->start <=> $a->start} @sub;
	}

	if ($#sub >=0) {
	    if ($sub[0]->strand == 1) {
		@sub = sort {$a->start <=> $b->start} @sub;
	    } else {
		@sub = sort {$b->start <=> $a->start} @sub;
	    }

	    my $count = 0;

	    foreach my $sf (@sub) {
		my $start = $count+1;
		my $end   = $start + ($sf->end - $sf->start);
#		print STDERR "Start/end $start $end :\n";
		$count = $end;

		my $gf = new Bio::EnsEMBL::SeqFeature(-start      => $start,
						      -end        => $end,
						      -strand     => $sf->strand,
						      -seqname    => $self->genomic->id,
						      -source_tag => 'genewise',
						      -primary_tag => 'similarity');

		my $fp = new Bio::EnsEMBL::FeaturePair(-feature1 => $sf,
						       -feature2 => $gf);
		
		$self->addGene($fp);
	    }
	}
    }
    #unlink $genfile;
    #unlink $protfile;

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







