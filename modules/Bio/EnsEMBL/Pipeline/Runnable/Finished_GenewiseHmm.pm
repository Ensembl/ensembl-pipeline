

=pod
=head1 NAME
  
Bio::EnsEMBL::Pipeline::Runnable::Finsihed_GenewiseHmm - Hmms to genomic sequence

=head1 SYNOPSIS

    # Initialise the GenewiseHmm module  
    my $genewise = new  Bio::EnsEMBL::Pipeline::Runnable::Finsiehd_GenewiseHmm  ('query'  => $genomic,
								       'hmmfile'  => $hmmfile,
								       'memory'   => $memory);


   $genewise->run;
    
   my @genes = $blastwise->output

=head1 DESCRIPTION

This modeule is created and run by halfwisehmm.pm, it takes a pfam hmm and some genomic sequence and runns genewise using these then parses the output to produce feature pairs representing the genes with exons as subseq features and supporting evidence as subseqfeatures of the exons

=head1 CONTACT

ensembl-dev <ensembl-dev@ebi.ac.uk>

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Pipeline::Runnable::Finished_GenewiseHmm;

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


=head2 new

  Args      : Bio::Seq, int, boolean, boolean, path, filename   
  Function  : create a Finished_GenewiseHmm object
  Returntype: object created
  Exceptions: throws if give no query sequence or no hmmfile
  Caller    : 
  Example   : 

=cut


sub new {
  my($class,@args) = @_;
  my $self = $class->SUPER::new(@args);
  
  $self->{'_output'} = [];
  $self->{'_reverse'} = undef;
  my ($query, $memory,$reverse,$endbias,$genewise,$hmmfile) = 
      $self->_rearrange([qw(QUERY MEMORY REVERSE ENDBIAS GENEWISE HMMFILE)], @args);
  
  #  print $query . "\n";
 
  $genewise = 'genewise' unless defined($genewise);

  $self->genewise($self->find_executable($genewise));

  $self->query($query) || $self->throw("No query sequence entered for blastwise");
  $self->hmmfile($hmmfile) || $self->throw("No Hmm file entered for Hmmgenewise");
  #print $reverse."\n";
  $self->is_reverse($reverse)   if (defined($reverse));             
  $self->endbias($endbias)   if (defined($endbias));             
  $self->memory ($memory)    if (defined($memory));

  return $self;
}


##################
#accessor methods#
##################


=head2 accessor_methods

  Arg [1]   : varies
  Function  : set a varible to a particular value
  Returntype: the variable
  Exceptions: some throw if not passed the correct arguement
  Caller    : 
  Example   : 

=cut


sub genewise {
  my ($self,$arg) = @_;

  if (defined($arg)) {
    $self->{_genewise} = $arg;
  }
  return $self->{_genewise};
}



# These all set/get or initializing methods

sub is_reverse {
  my ($self,$arg) = @_;
  #print "running is reverse\n";
  if ($arg==1) {
    $self->{'_reverse'} = $arg;
  }else{
    $self->{'_reverse'} = 0;
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

sub query {
    my ($self,$arg) = @_;

    if (defined($arg)) {
	$self->throw("Genomic sequence input is not a Bio::SeqI or Bio::Seq or Bio::PrimarySeqI") unless
	    ($arg->isa("Bio::SeqI") || 
	     $arg->isa("Bio::Seq")  || 
	     $arg->isa("Bio::PrimarySeqI"));
		
	
	$self->{'_query'} = $arg;
    }
    return $self->{'_query'};
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

######################
#running the analysis#
######################



=head2 set_environment

  Arg [1]   : none
  Function  : check that the WISECONFIGDIR is set
  Returntype: none
  Exceptions: thows if variable isn't a directory'
  Caller    : 
  Example   : 

=cut


sub set_environment {
  my ($self) = @_;
  
  if (! -d $ENV{WISECONFIGDIR}) {
    $self->throw("No WISECONFIGDIR ["  . $ENV{WISECONFIGDIR} . "]");
  }
}


=head2 run

  Arg [1]   : path to directory outputfiles are to be written too 
  Function  : calls methods which run genewise
  Returntype: none
  Exceptions: none
  Caller    : 
  Example   : 

=cut



sub run {
    my ($self, $dir) = @_;
    
    $self->set_environment;
    $self->_align_protein($dir);
}


=head2 _align_protein

  Arg [1]   : path to directory output files are to be written too  
  Function  : runs genewise and parse output producing arrays of genes, exons and supporting features
  Returntype: none
  Exceptions: throws if problems opening or closing pipe to genewise 
  Caller    : 
  Example   : 

=cut



sub _align_protein {
 
  my ($self, $dir) = @_;
  my $memory   = $self->memory;
  $self->workdir('/tmp') unless ($self->workdir($dir));
  $self->checkdir();
  my $gwfile   = "gw." . $$ . ".out";
  my $genfile  = "gw." . $$ . ".gen.fa";
  
  my $genio  = Bio::SeqIO->new (-file   => ">$genfile",
			      '-format' => 'fasta');
  my $gene_count = 0;
  #print "running genewise\n";
  $genio->write_seq($self->query);
  $genio = undef;
  #print"running genewiseHMM\n";
  my ($hmmname) = $self->hmmfile =~ /PF\d+/g;

  my $genewise = $self->genewise;

  my $hmm = $self->hmmfile;
  my $command = "$genewise -hmmer $hmm $genfile -genesf -kbyte $memory -ext 2 -gap 12 -subs 0.0000001 -quiet";
    
  if ($self->endbias == 1) {
    $command .= " -init endbias -splice flat ";
  }
  
  if ($self->is_reverse == 1) {
    $command .= " -trev ";
  }
  my $outputfile = $self->hmmfile.".output";
  #print STDERR "Command is $command\n";
  #system("pwd");
  #open(GW, "$command | tee -a $outputfile |") or $self->throw("error piping to genewise: $!\n");
  open(GW, "$command |") or $self->throw("error piping to genewise: $!\n");
  # for genesf parsing
  my @genesf_exons;
  my $curr_exon;
  my $curr_gene;
  
 #print "parseing output\n"; # making assumption of only 1 gene prediction ... this will change once we start using hmms
 GENESF: while (<GW>) {
    chomp;
    my @f = split;

    next unless (defined $f[0] && ($f[0] eq 'Gene' || $f[0] eq 'Exon' || $f[0] eq 'Supporting') && ($f[1] ne 'Paras:'));
    
    if($f[0] eq 'Gene'){
      # flag a frameshift - ultimately we will do something clever here but for now ...
      if($f[0] eq 'Gene'){
	if (/^(Gene\s+\d+)$/){
	  $curr_gene  = new Bio::EnsEMBL::SeqFeature;
	  $curr_gene->seqname($1);
	  #print STDERR "CURRENT GENE: ".$curr_gene->seqname."\n";
	  $curr_gene->start(1);
	  $curr_gene->end($self->query->length());
	  #print STDERR "adding gene to output\n";
	  $self->addGene($curr_gene);
	  
	  $gene_count++;
	}
      }

    }elsif($f[0] eq 'Exon'){
      #      make a new "exon"
      $curr_exon = new Bio::EnsEMBL::SeqFeature;
      $curr_exon->seqname  ($self->query->id);
      $curr_exon->id        ($hmmname);
      

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

      $curr_gene->add_sub_SeqFeature($curr_exon, '');
      push(@genesf_exons, $curr_exon);
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
	  my $pf = new Bio::EnsEMBL::SeqFeature( -start   => $pstart,
						 -end     => $pend,
						 -seqname => $hmmname,
						 -strand  => 1
					       ); 
	  my $gf  = new Bio::EnsEMBL::SeqFeature( -start   => $gstart,
						  -end     => $gend,
						  -seqname => 'genomic',
						  -strand  => $strand,
						);
	  my $fp = new Bio::EnsEMBL::FeaturePair( -feature2 => $pf,
						  -feature1 => $gf);
	  $curr_exon->add_sub_SeqFeature($fp,'');
    }

  }
  
  close(GW) or $self->throw("Error running genewise:$!\n");
  #print "there are ",$gene_count." genes predicted using ".$hmmname." hmm\n";
  unlink $genfile;
  unlink $gwfile;
  #unlink $outputfile;
}


=head2 addGene

  Arg [1]   : Bio:EnsEMBL::SeqFeature
  Function  : adds the seqfeature to the array
  Returntype: none
  Exceptions: throws if no passed a seqfeature
  Caller    : 
  Example   : 

=cut



sub addGene {
    my ($self,$arg) = @_;
    #print STDERR "adding ".$arg->seqname." to ouput\n";
    if($arg->isa("Bio::EnsEMBL::SeqFeature")){
      push(@{$self->{'_output'}},$arg);
    }else{
      $self->throw("this, $arg, should be a seqfeature\n");
    }
    #foreach my $result(@{$self->{'_output'}}){
    #  print STDERR $result->seqname."\n";
    #}
    #print STDERR "there are ".scalar(@{$self->{'_output'}})." genes\n";
}


=head2 eachGene

  Arg [1]   : none 
  Function  : return the array of seqfeatures
  Returntype: 
  Exceptions: none
  Caller    : 
  Example   : 

=cut


sub eachGene {
    my ($self) = @_;
   
    #print STDERR "outputting all data got\nthere are ".scalar(@{$self->{'_output'}})." genes\n";
    foreach my $result(@{$self->{'_output'}}){
      #print STDERR $result->seqname."\n";
    }
    return @{$self->{'_output'}};
}


################
#output methods#
################


=head2 output

  Arg [1]   : none
  Function  : returns the output array
  Returntype: 
  Exceptions: none
  Caller    : 
  Example   : 

=cut


sub output {
    my ($self) = @_;
    #print "output from genewiseHmm\n";
    my @out = $self->eachGene;
       
    return @out;
}







1;
