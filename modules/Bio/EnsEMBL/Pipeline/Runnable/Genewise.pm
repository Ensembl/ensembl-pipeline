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

# Object preamble - inherits from Bio::Root::RootI

use Bio::Root::RootI;
use Bio::EnsEMBL::Analysis::Programs qw(genewise); 
use Bio::EnsEMBL::SeqFeature;
use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::Pipeline::RunnableI;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI Bio::Root::RootI);

# _initialize is where the heavy stuff will happen when new is called

sub new {
  my ($class,@args) = @_;
  my $self = bless {}, $class;

  my ($genomic, $protein, $memory,$reverse,$endbias) = 
      $self->_rearrange(['GENOMIC','PROTEIN','MEMORY','REVERSE','ENDBIAS'], @args);
  
  print $genomic . "\n";

  $self->genomic($genomic) || $self->throw("No genomic sequence entered for blastwise");
  $self->protein($protein) || $self->throw("No protein sequence entered for blastwise");

  $self->reverse($reverse)   if (defined($reverse));             
  $self->endbias($endbias)   if (defined($endbias));             
  $self->memory ($memory)    if (defined($memory));

  return $self; # success - we hope!
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

  # for gff parsing
  my @gff_exons;
  my $curr_gff;
  my $count = 1;

  # for alb parsing
  my @alb_exons;
  my $curr_exon;
  my $exon_pos=0; # current position in exon, relative to start=1
  my $curr_aln = undef; # feature pair
  my $prev_state = "";  
  my $albphase = 0;
  
  # need to parse both gff and alb output
  # making assumption of only 1 gene prediction ...
  # gff o/p always comes before alb o/p
 GW: while (<GW>) {
    chomp;
    my @f = split;
    
    #      print STDERR "Line is " . $_ . "\n";
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
      my $prot_pos = $1 + 1;
      
      if ($f[4] eq '"SEQUENCE_INSERTION"' || $f[4] eq '"SEQUENCE_DELETION"') {
	# start a new "exon"
	$curr_exon = undef;
	$exon_pos = 0;
	$albphase = 0;
	$prev_state = 'SEQINDEL';
      } 
      
      if ($f[4] =~ /PHASE/) {
	# adjustment to start position of sub-aligned block
	$albphase = 0;
	$albphase = 1 if ($f[4] eq '"3SS_PHASE_2"' || $f[4] eq '"5SS_PHASE_2"' );
	$albphase = 2 if ($f[4] eq '"3SS_PHASE_1"'|| $f[4] eq '"5SS_PHASE_1"');
      }
      
      if ($f[2] eq '"MATCH_STATE"' && $f[4] eq '"CODON"') {
	if (!defined($curr_exon)) {
	  # start a new "exon"
	  my $ef = new Bio::EnsEMBL::SeqFeature;
	  $ef->seqname   ($self->protein->id);
	  $ef->id        ('exon_pred');
	  $ef->source_tag('genewise');
	  $ef->primary_tag('similarity');
	  $ef->start($prot_pos);
	  $ef->end($prot_pos);
	  
	  $curr_exon = $ef;
	  push(@alb_exons,$curr_exon);
	}
	
	if ($prev_state eq 'INSERT' || $prev_state eq 'DELETE' || $exon_pos == 0) {
	  if($exon_pos == 0) {
	    $exon_pos += $albphase;
	  }
	  
	  # start a new "alignment"
	  my $pf = new Bio::EnsEMBL::SeqFeature( -start   => $prot_pos,
						 -end     => $prot_pos,
						 -seqname => $self->protein->id,
					       ); 
	  my $gf  = new Bio::EnsEMBL::SeqFeature( -start   => $exon_pos +1,
						  -end     => $exon_pos +1,
						  -seqname => 'genomic',
						);
	  my $fp = new Bio::EnsEMBL::FeaturePair( -feature1 => $pf,
						  -feature2 => $gf);
	  $curr_aln = $fp;
	  $curr_exon->add_sub_SeqFeature($fp,'EXPAND');
	}
	
	$exon_pos += 3;
	$curr_exon->end($prot_pos);
	$curr_aln->end($prot_pos);      
	$curr_aln->hend($exon_pos);      
	$prev_state = 'MATCH';
      }
      
      elsif ( $f[2] eq '"INSERT_STATE"' ) {
	$exon_pos += 3 unless ($prev_state eq 'INTRON' || $prev_state eq 'SEQINDEL');
	$exon_pos += $albphase if ($prev_state eq 'INTRON');

	$prev_state = 'INSERT';
      }
      
      elsif ( $f[2] eq '"DELETE_STATE"' ) {
	$prev_state = 'DELETE';
      }
      
      elsif ( $f[4] eq '"CENTRAL_INTRON"' ) {
	$exon_pos = 0;
	$prev_state = 'INTRON';
	$curr_exon = undef;
      }   
    }
  } 
  
  close(GW);
  
  # Now convert into feature pairs
  # makes hstart & hend for the feature pair in cDNA coordinates?
  # Now convert into feature pairs
  my $dummy_ex;
  if(scalar(@gff_exons) != scalar(@alb_exons) ){
    print STDERR "gff: " .scalar(@gff_exons) . "\talb: " . scalar(@alb_exons) .  "\n";
    # still write the exons, but the supporting feature data is rubbish now.
    $self->warn("rats; exon numbers don't match - no sensible supporting feature data");
    $dummy_ex = 1;
  }
  
  my $index=0;
  my @fp;
  for ($index=0; $index <= $#gff_exons; $index++){
    my $feat_pair = new Bio::EnsEMBL::FeaturePair(-feature1 => $gff_exons[$index]);
    if($dummy_ex){
      my $f = new Bio::EnsEMBL::SeqFeature;
      $f->seqname   ($self->protein->id);
      $f->id        ('exon_pred');
      $f->source_tag('genewise');
      $f->primary_tag('similarity');
      $f->start(0);
      $f->end(0);
      $feat_pair->feature2($f);
    }
    else{
      $feat_pair->feature2($alb_exons[$index]);
    }
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
  
  # and sort out the sub_seqfeature coordinates
  # and transfer them to the genomic feature.
  foreach my $pair(@fp){
    my $strand = $pair->feature1->strand;
    foreach my $aln($pair->feature2->sub_SeqFeature) {
      $aln->strand(1);
      $aln->hstrand($strand); # genomic is feature 2. hmmm
      if ($strand == 1) {
	my $ex_start = $pair->feature1->start;
	my $ex_end = $pair->feature1->end;
	my $sp = $aln->hstart - 1;
	my $ep = $aln->hend - 1;
	$aln->hstart($sp + $ex_start);
	$aln->hend($ep + $ex_start);
      }
      else {
	my $ex_start = $pair->feature1->end;
	my $ex_end = $pair->feature1->start;
	my $sp = $aln->hstart - 1;
	my $ep = $aln->hend - 1;
	$aln->hend($ex_start - $sp);
	$aln->hstart($ex_start - $ep);
      }
      
      my $nfp = new Bio::EnsEMBL::FeaturePair( -feature1 => $aln->feature2,
					       -feature2 => $aln->feature1);

      # final screening out of off by one errors
      if($nfp->feature1->start  >= $pair->feature1->start &&
	 $nfp->feature1->end    <= $pair->feature1->end    &&
	 $nfp->feature1->strand == $pair->feature1->strand){

	$pair->feature1->add_sub_SeqFeature($nfp,'');
      }
      
      else{
	print STDERR "\n***double rats: coordinate mismatch; can't write supporting_features\n"; 
	print STDERR "\tpair->feature1\n";
	print STDERR "\t" . $pair->feature1->start . " " . $pair->feature1->end . " " . $pair->feature1->strand. "\n";
	print STDERR "\tnfp\n";
	print STDERR "\t" . $nfp->feature1->start . " " . $nfp->feature1->end . " " . $nfp->feature1->strand. "\n";
      }
    }
    
    # and flush feature2
    $pair->feature2->flush_sub_SeqFeature;
    $self->addExon($pair);
  }
  
  # check them
  foreach my $pair(@fp){
#    print STDERR $pair->seqname . "\t" . $pair->start . "\t" . $pair->end . "\t" . 
      $pair->score . "\t" . $pair->strand . "\t" . $pair->hseqname . "\t" . 
	$pair->hstart . "\t" . $pair->hend. "\n";#
    
    foreach my $aln($pair->feature2->sub_SeqFeature) {
#      print "\t" . $aln->seqname  . " " . $aln->start  . "-" . $aln->end  . " " . $aln->strand ."\t" . 
	$aln->hseqname . " " . $aln->hstart . "-" . $aln->hend . " " . $aln->hstrand ."\n";  
    }
    
    foreach my $aln($pair->feature1->sub_SeqFeature) {
#      print "\t" . $aln->seqname  . " " . $aln->start  . "-" . $aln->end  . " " . $aln->strand ."\t" . 
	$aln->hseqname . " " . $aln->hstart . "-" . $aln->hend . " " . $aln->hstrand ."\n";  
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







