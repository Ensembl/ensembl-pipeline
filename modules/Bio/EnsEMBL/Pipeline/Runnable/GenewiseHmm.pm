=pod

=head1 NAME
  
Bio::EnsEMBL::Pipeline::Runnable::GenewiseHmm - Hmms to genomic sequence

=head1 SYNOPSIS

# Initialise the GenewiseHmm module  
my $genewise = new  Bio::EnsEMBL::Pipeline::Runnable::GenewiseHmm  ('genomic'  => $genomic,
								    'hmmfile'  => $hmmfile,
								    'memory'   => $memory);


$genewise->run;

my @genes = $blastwise->output

=head1 DESCRIPTION


=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Pipeline::Runnable::GenewiseHmm;

use strict;  
use vars   qw(@ISA);

# Object preamble - inherits from Bio::EnsEMBL::Root

use Bio::EnsEMBL::Root;
use Bio::EnsEMBL::Analysis::Programs qw(genewise); 
use Bio::EnsEMBL::SeqFeature;
use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::Pipeline::RunnableI;
use Bio::SeqIO;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI);

sub new {
  my($class,@args) = @_;
  my $self = $class->SUPER::new(@args);
  
  my ($genomic, $memory,$reverse,$endbias,$genewise,$hmmfile) = 
    $self->_rearrange([qw(
			  GENOMIC 
			  MEMORY 
			  REVERSE 
			  ENDBIAS 
			  GENEWISE 
			  HMMFILE
			 )], @args);
  

  $genewise = 'genewise' unless defined($genewise);

  $self->genewise($self->find_executable($genewise));

  $self->genomic($genomic) || $self->throw("No genomic sequence entered for blastwise");
  $self->hmmfile($hmmfile) || $self->throw("No Hmm file entered for Hmmgenewise");

  $self->reverse($reverse)   if (defined($reverse));             
  $self->endbias($endbias)   if (defined($endbias));             
  $self->memory ($memory)    if (defined($memory));

  return $self;
}

sub genewise {
  my ($self,$arg) = @_;

  if (defined($arg)) {
    $self->{_genewise} = $arg;
  }
  return $self->{_genewise};
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
  
  if (! -d $ENV{WISECONFIGDIR}) {
    $self->throw("No WISECONFIGDIR ["  . $ENV{WISECONFIGDIR} . "]");
  }
}

sub _wise_init {
  my($self,$genefile) = @_;

#  print("Initializing genewise\n");

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
  
  my $genio  = Bio::SeqIO->new (-file   => ">$genfile",
			      '-format' => 'fasta');
  
  $genio->write_seq($self->genomic);
  $genio = undef;

  
  #print STDERR "hmmfile: ". $self->hmmfile."\n";
  my ($hmmname) = $self->hmmfile =~ /\w+.hmm$/g;

  #print STDERR "hmmname: ". $hmmname."\n";

  my $genewise = $self->genewise;

  my $hmm = $self->hmmfile;
  #my $command = "genewise -hmmer $hmm $genfile -alg 623L -genesf -kbyte $memory -ext 2 -gap 12 -subs 0.0000001 -quiet";
  my $command = "genewise -hmmer $hmm $genfile -alg 623L -genesf -both -quiet";
  #my $command = "genewise -hname OR -hmmer $hmm $genfile -alg 623L -genesf -quiet";
    
  if ($self->endbias == 1) {
    $command .= " -init endbias -splice flat ";
  }  
  if ($self->reverse) {
    $command .= " -trev ";
  }
  
  print STDERR "Command is $command\n";
  
  open(GW, "$command |") or $self->throw("error piping to genewise: $!\n");

  my $curr_exon;
  my $curr_gene;
  my @scores;
  my $score;
  my $genecount = 0;
  my $scorecount = 0;
 GENESF: 
  while (<GW>) {
    #print STDERR $_;
    
    chomp;
    my @f = split;
    for(my$i=0; $i<$#f;$i++){
      if ( $f[$i] eq 'Score' ){
	$scores[$scorecount] = $f[$i+1];
	$scorecount++;
	#print STDERR "%%%% score: $score\n";
	last;
      }
    }
    
    next unless (defined $f[0] && ($f[0] eq 'Gene' || $f[0] eq 'Exon' || $f[0] eq 'Supporting'));
    
    if($f[0] eq 'Gene' && $f[1] && $f[2] && !($f[1] eq 'Paras:')){
      #print STDERR "creating new gene\n";
      $curr_gene = Bio::EnsEMBL::SeqFeature->new();
      $genecount++;
      $self->addGene($curr_gene);
      if (/^(Gene\s+\d+)$/){
	#print STDERR "1: $1\n";
      }
    }
    
    elsif($f[0] eq 'Exon'){
      #      make a new "exon"
      #print STDERR "Making new exon\n";
      $curr_exon = new Bio::EnsEMBL::SeqFeature;
      $curr_exon->seqname  ($self->genomic->id);
      $curr_exon->id       ($hmmname);
      #$curr_exon->source_tag('genewise');
      #$curr_exon->primary_tag('similarity');
      
      $curr_exon->phase($f[4]);
      
      my $start  = $f[1];
      my $end    = $f[2];
      my $strand = 1;
      if($f[1] > $f[2]){
	$strand = -1;
	$start  = $f[2];
	$end    = $f[1];
      }
      
      $curr_exon->start($start);
      $curr_exon->end($end);
      $curr_exon->strand($strand);
      
      #$self->addExon($curr_exon);
      #push(@genesf_exons, $curr_exon);
      $curr_gene->add_sub_SeqFeature($curr_exon,'EXPAND');
      #print STDERR "gene contains ".scalar( $curr_gene->sub_SeqFeature )." sub features\n";
    }
    
    elsif($f[0] eq 'Supporting'){
      my $gstart = $f[1];
      my $gend   = $f[2];
      my $strand = 1;
      if ($gstart > $gend){
	$gstart = $f[2];
	$gend = $f[1];
	$strand = -1;
      }
      
      # check strand sanity
      if($strand != $curr_exon->strand){
	$self->warn("incompatible strands between exon and supporting feature - cannot add suppfeat\n");
	next GENESF;
      }
      
      my $pstart = $f[3];
      my $pend   = $f[4];
      if($pstart > $pend){
	$self->warn("Protein start greater than end! Skipping this suppfeat\n");
	next GENESF;
      }

      # start a new "alignment"
      #print STDERR "### adding score ".$scores[$genecount-1]." ###\n";
      my $pf = new Bio::EnsEMBL::SeqFeature( -start   => $pstart,
					     -end     => $pend,
					     -seqname => $hmmname,
					     -strand  => 1,
					     -score   => $scores[$genecount-1],
					   ); 
      my $gf  = new Bio::EnsEMBL::SeqFeature( -start   => $gstart,
					      -end     => $gend,
					      -seqname => 'genomic',
					      -strand  => $strand,
					      -score   => $scores[$genecount-1],
					    );
      my $fp = new Bio::EnsEMBL::FeaturePair( -feature2 => $pf,
					      -feature1 => $gf);
      
      $curr_exon->add_sub_SeqFeature($fp,'');
      #print STDERR "exon contains ".scalar( $curr_exon->sub_SeqFeature )." sub features\n";
    }
    
  } 
  
  close(GW) or $self->throw("Error running genewise:$!\n");
  
  unlink $genfile;
  unlink $gwfile;
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

sub genomic {
    my ($self,$arg) = @_;

    if (defined($arg)) {
      $self->throw("Genomic sequence input is not a Bio::SeqI or Bio::Seq or Bio::PrimarySeq") unless
	($arg->isa("Bio::SeqI") || 
	 $arg->isa("Bio::Seq")  || 
	 $arg->isa("Bio::PrimarySeq"));
		
	
	$self->{'_genomic'} = $arg;
    }
    return $self->{'_genomic'};
}

=head2 hmmfile

 Title   : hmmfile
 Usage   : $obj->hmmfile($newval)
 Function: 
 Returns : value of hmmfile
 Args    : newvalue (optional)


=cut

sub hmmfile{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'hmmfile'} = $value;
    }
    return $obj->{'hmmfile'};

}

#sub addExon {
#    my ($self,$arg) = @_;
#
#    if (!defined($self->{'_output'})) {
#	$self->{'_output'} = [];
#    }
#    push(@{$self->{'_output'}},$arg);
#}

sub addGene {
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
