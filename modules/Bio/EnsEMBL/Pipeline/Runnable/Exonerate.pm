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

Bio::EnsEMBL::Pipeline::Runnable::Exonerate

=head1 SYNOPSIS

    my $obj = Bio::EnsEMBL::Pipeline::Runnable::Exonerate->new(
                                             -genomic => $genseq,
                                             -est     => $estseq,
					     -outfile => $outfile
                                             );
    or
    
    my $obj = Bio::EnsEMBL::Pipeline::Runnable::Exonerate->new()

=head1 DESCRIPTION

Exonerate is a fast EST:genomic alignment program written by Guy Slater.
This object runs exonerate over input EST and genomic sequences, and stores the 
exonerate matches as an array of Bio::EnsEMBL::FeaturePair

The passed in $genseq and $estseq can be filenames or references to arrays of Bio::Seq; exonerate 
runs faster if given multiple query sequences at once. 

=head2 Methods:

 new,
 genomic_sequence,
 est_sequence,
 run,
 output.

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::EnsEMBL::Pipeline::Runnable::Exonerate;

use vars qw(@ISA);
use strict;
# Object preamble - inherits from Bio::Root::RootI;

use Bio::EnsEMBL::Pipeline::RunnableI;
use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::SeqFeature;
use Bio::EnsEMBL::Analysis;
#compile time check for executable - won't work till centrally installed ...
use Bio::EnsEMBL::Analysis::Programs qw(exonerate); 
use Bio::PrimarySeq;
use Bio::SeqIO;
use Bio::Root::RootI;
use FileHandle;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI Bio::Root::RootI );

sub new {
  my ($class,@args) = @_;
  my $self = bless {}, $class;
  
  $self->{'_fplist'}    = []; #create key to an array of feature pairs
  $self->{'_clone'}     = undef;        #location of Bio::Seq object
  $self->{'_exonerate'} = undef;     #location of exonerate
  $self->{'_workdir'}   = undef;     #location of temp directory
  $self->{'_genfilename'}  =undef;      #file to store Bio::Seq object
  $self->{'_estfilename'} = undef;   #file to store EST Bio::Seq object
  $self->{'_results'}     = undef;      #file to store results of analysis
  $self->{'_protected'}   = [];         #a list of files protected from deletion
  $self->{'_arguments'}   = undef;      #arguments for exonerate
  
  my( $genomic, $est, $exonerate, $arguments) = 
    $self->_rearrange([qw(GENOMIC EST EXONERATE ARGS)], @args);
  
  $self->throw("no genomic sequence given\n") unless defined $genomic;
  $self->genomic_sequence($genomic) if $genomic; #create & fill key to Bio::Seq

  # allow it to take a filename or an array ref
  $self->throw("no est sequence given\n") unless defined $est;
  $self->est_sequence($est) if defined $est; 

  if ($exonerate) {   
    $self->exonerate($exonerate) ;
  }
  else {   
    eval 
      { $self->exonerate($self->locate_executable('exonerate')); };
    if ($@)
      # need a central installation ...
      { $self->throw("Can't find exonerate!"); }
  }
  if ($arguments) 
    {   $self->arguments($arguments) ;}
  
  return $self;
}

#################
# get/set methods 
#################
=head2 genomic_sequence
  
    Title   :   genomic_sequence
    Usage   :   $self->genomic_sequence($seq)
    Function:   Get/set method for genomic sequences
    Returns :   Bio::PrimarySeqI, or filename
    Args    :   Bio::PrimarySeqI, or filename

=cut

sub genomic_sequence {
  my( $self, $value ) = @_;    


  
  if (defined $value) {
    if (!ref($value)){
      # assume it's a filename - check the file exists
      $self->throw("[$value] : file does not exist\n") unless -e $value;
      $self->genfilename($value);
      $self->{'_genomic_sequence'} = $value;
    }
    elsif( $value->isa("Bio::PrimarySeqI") ){
      $self->{'_genomic_sequence'} = $value;
      $self->genfilename("/tmp/genfile_.$$.fn");
    }
    else {
      $self->throw("$value is neither a Bio::Seq  nor a filename\n");
    }
  }
  
  return $self->{'_genomic_sequence'};
}

=head2 est_sequence

    Title   :   est_sequence
    Usage   :   $self->est_sequence($seq)
    Function:   Get/set method for est sequences
    Returns :   reference to an array of Bio::Seq
    Args    :   Either a filename or a reference to an array of Bio::Seq

=cut

sub est_sequence {
  my( $self, $value ) = @_;

  if ($value) {
    if (ref($value) eq 'ARRAY') {
      foreach my $est(@$value) {
	$est->isa("Bio::PrimarySeqI") || $self->throw("Input isn't a Bio::PrimarySeqI");
      }
      $self->{'_est_sequences'} = $value;
      $self->estfilename("/tmp/estfile_.$$.fn");
    }
    else {
      # it's a filename - check the file exists
      $self->throw("[$value] : file does not exist\n") unless -e $value;
      $self->estfilename($value);
      $self->{'_est_sequences'} = $value;
    }
  }
  
  #NB ref to an array of Bio::Seq
  return $self->{'_est_sequences'};
}

=head2 exonerate

    Title   :   exonerate
    Usage   :   $self->exonerate('/path/to/executable')
    Function:   Get/set method for exonerate executable path
    Returns :   
    Args    :   

=cut

sub exonerate {
  my ($self,$arg) = @_;
  
  if (defined($arg)) {
    $self->{'_exonerate'} = $arg;
  }
  return $self->{'_exonerate'};
}

=head2 arguments

    Title   :   arguments
    Usage   :   $self->est_sequence($args)
    Function:   Get/set method for exonerate arguments
    Returns :   
    Args    :   

=cut

sub arguments {
  my ($self, $args) = @_;
  if ($args)
    {
      $self->{'_arguments'} = $args ;
    }
  return $self->{'_arguments'};
}


=head2 estfilename

    Title   :   estfilename
    Usage   :   $self->estfilename($filename)
    Function:   Get/set method for estfilename
    Returns :   
    Args    :   

=cut

sub estfilename {
  my ($self, $estfilename) = @_;
  $self->{'_estfilename'} = $estfilename if ($estfilename);
  return $self->{'_estfilename'};
}


=head2 genfilename

    Title   :   genfilename
    Usage   :   $self->genfilename($filename)
    Function:   Get/set method for genfilename
    Returns :   
    Args    :   

=cut

sub genfilename {
  my ($self, $genfilename) = @_;
  $self->{'_genfilename'} = $genfilename if ($genfilename);
  return $self->{'_genfilename'};
}

=head2 run

  Title   : run
  Usage   : $self->run()
  Function: Runs exonerate and stores results as FeaturePairs
  Returns : TRUE on success, FALSE on failure.
  Args    : 

=cut

sub run {
  my ($self, @args) = @_;

  #check inputs
  my $genomicseq = $self->genomic_sequence ||
    $self->throw("Genomic sequences not provided");
  my $estseq = $self->est_sequence ||
    $self->throw("EST sequences not provided");
  
  #extract filenames from args and check/create files and directory
  my $genfile = $self->genfilename;
  my $estfile = $self->estfilename;
  
  # do we need to write out the genomic sequence?
  if($genomicseq ne $genfile){
    my $genOutput = Bio::SeqIO->new(-file => ">$genfile" , '-format' => 'Fasta')
      or $self->throw("Can't create new Bio::SeqIO from $genfile '$' : $!");
    
    $genOutput->write_seq($genomicseq);
  }
    
  # do we need to write out the est sequences?
  if(ref($estseq) eq 'ARRAY'){
    my $estOutput = Bio::SeqIO->new(-file => ">$estfile" , '-format' => 'Fasta')
      or $self->throw("Can't create new Bio::SeqIO from $estfile '$' : $!");
    
    foreach my $eseq(@$estseq) {
      $estOutput->write_seq($eseq);
    }
  }
  
  # finally we run the beast.
  my $command = $self->exonerate() . " -w 14 -t 65 -H 100 -D 15 -m 500 -n yes -A false -G yes --cdna $estfile --genomic $genfile";
  
STDOUT->autoflush(1);

  eval {
#    print (STDERR "Running command $command\n");
    open(EXONERATE, "$command |") or $self->throw("Error forking exonerate pipe on ".$self->genfilename."\n"); 
    # some constant strings
    my $source_tag  = "exonerate";
    
    #read output
    my $estname   = undef;
    my $genname   = "";
    my %gff;
    my $prevstate = "cigar";

    # exonerate output is split up by est.
    # gff comes first, then all the cigar lines
    # so we can process all the gff, then all the cigar lines.
    # if we keep track of when we flip from cigar to gff, we know we're moving on to the next est.
    while(<EXONERATE>){
      # gff
      if(/\S+\s+exonerate\s+exon\s+\d+/){

	# is this the first gff? ie have we moved onto a new est?
	if($prevstate eq 'cigar') {
	  $prevstate = 'gff';
	  $estname = undef;
	  foreach my $entry( keys %gff) {
	    delete $gff{$entry};
	  }
	}

	# M93650  exonerate       exon    14193   14292   .       +       .       gene_id 8       ;       
	# exon_id 1       ;       insertions      0       ;       deletions 0       ;       
	# percent_id      100.000
	my @cols = split /;/;
	my $pid = $cols[$#cols];
	$pid =~ /(\d+)\./;
	$pid = $1;
	
	my @coords = split /\s+/, $cols[0];
	$estname = $coords[0] unless defined $estname;
	$self->throw("$coords[0] does not match $estname") unless $coords[0] eq $estname;

	# key by start-end-strand
	my $start = $coords[3]; 
	$start--;# presently off by one wrt cigar
	my $idstring = $start . "-" . $coords[4] . "-" . $coords[6];
	$gff{$idstring} = $pid;
      }

      #cigar
      elsif ($_ =~ /cigar/) {
	next if($_ =~ /^Message/);
	
	$prevstate = 'cigar';

	#split on whitespace
	my @elements = split;
	
	next unless $elements[0] eq 'cigar:';
	
	my $primary_tag    = 'exonerate';
	my $source_tag     = 'exonerate';
	my $genomic_start  = $elements[6];
	my $genomic_end    = $elements[7];
	my $genomic_id     = $elements[5];
	my $est_id         = $elements[1];
	my $est_start      = $elements[2];
	my $est_end        = $elements[3];
	
	# can we continue?
	$self->throw("$estname(gff) and $est_id(cigar) do not match!") unless $estname eq $est_id;

	# find percent id from %gff
	my $querystr = $genomic_start . "-" . $genomic_end . "-" . $elements[8];
	my $pid = $gff{$querystr};

	# not every cigar prediction has a corresponding gff exon
	if(!defined($pid)) { $pid = -1; }

	# est seqname
	my $genomic_score  = $elements[9];
	if ($genomic_score eq ".") { $genomic_score = 0; } 
	my $est_score = $genomic_score;
	
	my $genomic_source = $source_tag;
	my $est_source = $source_tag;
	
	my $genomic_strand = 1;
	if ($elements[8] eq '-') {
	  $genomic_strand = -1;
	}
	my $est_strand = 1;
	if ($elements[4] eq '-') {
	  $est_strand = -1;
	}

	my $genomic_primary = $primary_tag;
	my $est_primary = $primary_tag;

	my $pair = $self->_createfeatures ($genomic_score, $pid, $genomic_start, $genomic_end, $genomic_id, 
					 $est_start, $est_end, $est_id, 
					 $genomic_source, $est_source, 
					 $genomic_strand, $est_strand, 
					 $genomic_primary, $est_primary);

# needed for overall run, otherwise not
	print $pair->seqname  . "\t" . $pair->start        . "\t" . $pair->end    . "\t" . 
	          $pair->score    . "\t" . $pair->percent_id   . "\t" . $pair->strand . "\t" . 
		  $pair->hseqname . "\t" . $pair->hstart       . "\t" . $pair->hend   . "\t" . 
                  $pair->hstrand  . "\n";

      }    
    }
    
    
    close (EXONERATE) or $self->throw("Error running exonerate: $command");
    
  };  

#  close (OUT);

  #clean up temp files
  if(ref($estseq) eq 'ARRAY' ) {
    $self->_deletefiles($estfile);
  }
  
  if($genomicseq ne $genfile) {
    $self->_deletefiles($genfile);
  }
  
  if ($@) {
    $self->throw("Error running exonerate :$@ \n");
  } 
  else {
    return (1);
  }
}

=head2 convert_output
  
    Title   :   convert_output
    Usage   :   $obj->convert_output
    Function:   Parses exonerate output to make "gene" features with exons as subfeatures
                At present pretty hopeless supporting data for exons themselves ...
    Returns :   none
    Args    :   

=cut

sub convert_output {
  my($self) = @_;

  my @tmpf = @{$self->{'_fplist'}};

  #ought to be done here but will probably break Exonerate RunnableDB ...
  my $curr_gene;
  my @genes;
  foreach my $fp(@tmpf){
    # if it's a gene line, make a gene
    if($fp->primary_tag eq 'gene'){
      $curr_gene = $fp->feature1;
      push (@genes,$curr_gene)
    }
    # if it's an exon line, make an exon and add it to the current gene
    elsif($fp->primary_tag eq 'exon'){
      $self->throw("Can't cope with an exon without having a gene first!") unless defined ($curr_gene);
      $curr_gene->add_sub_SeqFeature($fp->feature1,'');
      $fp->feature1->add_sub_SeqFeature($fp,'');
    }
    # otherwise skip it
  }
  
  if (!defined($self->{'_output'})) {
      $self->{'_output'} = [];
  }
  print "genes are @genes\n";
  push(@{$self->{'_output'}},@genes);

}

=head2 output

  Title   : output
  Usage   : $self->output
  Function: Returns results of exonerate as array of FeaturePair, one per gene, 
            with exons as sub_seqfeatures
  Returns : An array of Bio::EnsEMBL::FeaturePair
  Args    : none

=cut

sub output {
  my ($self) = @_;

  #  could return a list of featurepairs, one per gene, as with Genewise.

# TEMPORARY FIXME
#  return @{$self->{'_output'}};
#  return @{$self->{'_fplist'}};
  return $self->{'_fplist'}; # ref to an array
}

=head2 _create_features

  Title   : _create_features
  Usage   : $self->_create_features($f1score, $f1pid$f1start, $f1end, $f1id, 
				    $f2start, $f2end, $f2id, $f1source, 
				    $f2source, $f1strand, $f2strand, 
				    $f1primary, $f2primary)
  Function: Returns results of exonerate as array of FeaturePair
  Returns : Nothing, but $self->{'_fplist'} contains a new FeaturePair
  Args    : 

=cut

sub _createfeatures {
  my ($self, $f1score, $f1pid, $f1start, $f1end, $f1id, $f2start, $f2end, $f2id,
      $f1source, $f2source, $f1strand, $f2strand, $f1primary, $f2primary) = @_;

  #create analysis object
  my $analysis_obj    = new Bio::EnsEMBL::Analysis
    (-db              => "none",
     -db_version      => "none",
     -program         => "exonerate",
     -program_version => "1",
     -gff_source      => $f1source,
     -gff_feature     => $f1primary,);
  
  
  #create features
  my $feat1 = new Bio::EnsEMBL::SeqFeature  (-start      =>   $f1start,
					     -end         =>   $f1end,
					     -seqname     =>   $f1id,
					     -strand      =>   $f1strand,
					     -score       =>   $f1score,
					     -percent_id  =>   $f1pid, 
					     -source_tag  =>   $f1source,
					     -primary_tag =>   $f1primary,
					     -analysis    =>   $analysis_obj );
  
  my $feat2 = new Bio::EnsEMBL::SeqFeature  (-start       =>   $f2start,
					     -end         =>   $f2end,
					     -seqname     =>   $f2id,
					     -strand      =>   $f2strand,
					     -score       =>   $f1score,
					     -percent_id  =>   $f1pid, 
					     -source_tag  =>   $f2source,
					     -primary_tag =>   $f2primary,
					     -analysis    =>   $analysis_obj );
  #create featurepair
  my $fp = new Bio::EnsEMBL::FeaturePair  (-feature1 => $feat1,
					   -feature2 => $feat2) ;
  
  if ($fp) {
    $self->throw("Can't validate") unless $fp->validate();
    push(@{$self->{'_fplist'}}, $fp);
  }
  return $fp;
}

#####################################
# creating and clearing up temp files
#####################################

# a lot of this is shared with eg Vert_Est2Genome. Need a common parent ...
sub _createfiles {
  my ($self, $genfile, $estfile, $dirname)= @_;
  
  #check for diskspace
  my $spacelimit = 0.01; # 0.01Gb or about 10 MB
  my $dir ="./";
  unless ($self->_diskspace($dir, $spacelimit)) 
    {
      $self->throw("Not enough disk space ($spacelimit Gb required)");
    }
  
  #if names not provided create unique names based on process ID    
  $genfile = $self->_getname("genfile") unless ($genfile);
  $estfile = $self->_getname("estfile") unless ($estfile);    
  
  # Should check we can write to this directory 
  $self->throw("No directory $dirname") unless -e $dirname;
  
  return ($genfile, $estfile);
}

sub _getname {
  my ($self, $typename) = @_;
  return  $typename."_".$$.".fn"; 
}

sub _diskspace {
  my ($self, $dir, $limit) =@_;
  my $block_size; #could be used where block size != 512 ?
  my $Gb = 1024 ** 3;
  
  open DF, "df $dir |" or $self->throw ("Can't open 'du' pipe");
  while (<DF>) 
    {
      if ($block_size) 
        {
	  my @L = split;
	  my $space_in_Gb = $L[3] * 512 / $Gb;
	  return 0 if ($space_in_Gb < $limit);
	  return 1;
        } 
      else 
        {
	  ($block_size) = /(\d+).+blocks/i
	    || $self->throw ("Can't determine block size from:\n$_");
        }
    }
  close DF || $self->throw("Error from 'df' : $!");
}

sub _deletefiles {
  my ($self, @files) = @_;
  
  my $unlinked = unlink(@files);
  
  if ($unlinked == @files) {
    return 1;
  } else {
    my @fails = grep -e, @files;
    $self->throw("Failed to remove @fails : $!\n");
  }
}

1;
