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
# Object preamble

use Bio::EnsEMBL::Pipeline::RunnableI;
use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::SeqFeature;
use Bio::EnsEMBL::Analysis;
use Bio::PrimarySeq;
use Bio::SeqIO;
use Bio::EnsEMBL::Root;
use Bio::EnsEMBL::DnaDnaAlignFeature;
use FileHandle;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI Bio::EnsEMBL::Root );

sub new {
  my ($class,@args) = @_;
  my $self = bless {}, $class;
  
  my( $genomic, $est, $exonerate, $arguments, $print) = 
    $self->_rearrange([qw(GENOMIC EST EXONERATE ARGS PRINT)], @args);

  $self->throw("no genomic sequence given\n") unless defined $genomic;
  $self->throw("no est sequence given\n")     unless defined $est;

  $self->genomic_sequence($genomic) if $genomic; #create & fill key to Bio::Seq
  $self->est_sequence($est) if defined $est; 

  if(defined $exonerate){
    $self->exonerate($exonerate);
  } else {
    $self->exonerate($self->find_executable('exonerate')); 
  }

  if (defined $arguments) {   
    $self->arguments($arguments);
  }
  
  if (defined $print) {
    $self->print_results($print);
  }

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
      my $filename = "/tmp/genfile_$$.fn";
      $self->genfilename($filename);
      my $genOutput = Bio::SeqIO->new(-file => ">$filename" , '-format' => "Fasta")
      or $self->throw("Can't create new Bio::SeqIO from $filename '$' : $!");
    
      $self->throw ("problem writing genomic seqeunce to $filename\n" ) unless $genOutput->write_seq($value);
      
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
      
      my $filename = "/tmp/estfile_$$.fn";
      $self->estfilename($filename);
      my $estOutput = Bio::SeqIO->new(-file => ">$filename" , '-format' => 'Fasta')
	or $self->throw("Can't create new Bio::SeqIO from $filename '$' : $!");
      
      foreach my $eseq(@$value) {
	$self->throw ("problem writing est seqeunce to $filename\n" ) unless $estOutput->write_seq($eseq);
      }
      
    }
    else {
      print STDERR $value." must be a filename checking\n";
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
  if ($args) {
      $self->{'_arguments'} = $args ;
    }
  return $self->{'_arguments'};
}

=head2 print_results

    Title   :   print_results
    Usage   :   $self->print_results(1)
    Function:   Get/set method for a flag determining whether exonerate output should be dumped to STDOUT
    Returns :   0 or 1
    Args    :   1 or nothing

=cut

sub print_results {
  my ($self, $print) = @_;

  if (!defined $self->{'_print'}) {
     $self->{'_print'} = 0;
  }

  if ($print) {
      $self->{'_print'} = 1 ;
    }

  return $self->{'_print'};
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
  Usage   : $self->run
  Function: Runs exonerate; the optional $ungapped flag determines whether to run
            ungapped - default is gapped.
  Returns : nothing
  Args    : $ungapped flag - optional

=cut

sub run {  
  my ($self, $ungapped) = @_;

  if(defined $ungapped){
    $self->run_ungapped;
  }
  else{
    $self->run_gapped;
  }
}

=head2 run_gapped

  Title   : run_gapped
  Usage   : $self->run_gapped()
  Function: Runs exonerate in gapped mode and stores results as FeaturePairs
  Returns : TRUE on success, FALSE on failure.
  Args    : none 

=cut

sub run_gapped {
  my ($self) = @_;

  # this method needs serious tidying
  #check inputs
  my $genomicseq = $self->genomic_sequence ||
    $self->throw("Genomic sequences not provided");
  my $estseq = $self->est_sequence ||
    $self->throw("EST sequences not provided");
  
  #extract filenames from args and check/create files and directory
  my $genfile = $self->genfilename;
  my $estfile = $self->estfilename;
  
  # finally we run the beast.
  my $resfile = "/tmp/exonerate_resfile_" . $$;
  my $exonerate_command = $self->exonerate() . " -n true -A false --cdna $estfile --genomic $genfile >" . $resfile;
  
  eval {
    # better to do this as a pipe?
    print (STDERR "Running command $exonerate_command\n");
    $self->throw("Error running exonerate on ".$self->filename."\n") 
      if (system ($exonerate_command)); 
    $self->parse_results($resfile);
    $self->convert_output();
  };  
  
  #clean up temp files
  if(ref($estseq) eq 'ARRAY'){
    $self->_deletefiles($genfile, $estfile, $resfile);
  }
  else {
    $self->_deletefiles($genfile, $resfile);
  }
  if ($@) {
    $self->throw("Error running exonerate :$@ \n");
  } 
  else {
    return (1);
  }
  # store raw output FeaturePairs in $self->output
  # store "genes" in $self->add_genes via convert_output
  # $self->convert_output;
}

=head2 parse_results
  
    Title   :   parse_results
    Usage   :   $obj->parse_results($filename)
    Function:   Parses exonerate output to give a set of features
                parsefile can accept filenames, filehandles or pipes (\*STDIN)
    Returns :   none
    Args    :   optional filename

=cut

sub parse_results {
  my ($self, $resfile) = @_;

  # some constant strings
  my $source_tag  = "exonerate";
  #  my $primary_tag = "similarity";
  if (-e $resfile) {
    open (EXONERATE, "<$resfile") or $self->throw("Error opening ", $resfile, " \n");#
  }
  else {
    $self->throw("Can't open $resfile :$!\n");
  }
  
  #read output
  my $queryname = "";
  while (<EXONERATE>) {
    
    # VAC temporary-ish changes - we're not using exonerate fully at the moment.
    # output parsing needs to be rahashed once exonerate is stable
    
    #    if ($_ =~ /exonerate/) {
    if ($_ =~ /cigar/) {
      next if($_ =~ /^Message/);
      
      #split on whitespace
      my @elements = split;

      #      if( $elements[1] ne 'exonerate' ) { next; }
      next unless $elements[0] eq 'cigar:';
      
      if($_ =~ /query\s+\"(\w+)\"/) {
	$queryname = $1;
      }
      
      # cigar parsing needs a LOT of looking at
      # cigar: gi|550092|dbj|D29023.1|D29023 0 311 + static0 1996 2307 + 1555.00 M 311
      #extract values from output line [0] - [7]
      
      #      my $primary_tag = $elements[2];
      my $primary_tag = 'Gene';
      
      my $genomic_start  = $elements[6];
      my $genomic_end    = $elements[7];
      my $genomic_id     = $elements[5];
      # start & end on EST sequence are not currently given by exonerate output ...
      my $est_id     = $elements[1];
      my $est_start  = $elements[2];
      my $est_end    = $elements[3];
      
      # est seqname
      my $genomic_score  = $elements[9];
      if ($genomic_score eq ".") { $genomic_score = 0; } 
      my $est_score = $genomic_score;
      
      my $genomic_source = $source_tag;
      my $est_source = $source_tag;
      
      my $genomic_strand = 1;
      if ($elements[4] eq '-') {
	$genomic_strand = -1;
      }
      my $est_strand = 1;
      if ($elements[8] eq '-') {
	$est_strand = -1;
      }
      
      # currently doesn't deal well with - strand ... genes
      my $genomic_primary = $primary_tag;
      my $est_primary = $primary_tag;
      
      my $pair = $self->_create_featurepair ($genomic_score, $genomic_id, $genomic_start, $genomic_end, 
					     $est_id, $est_start, $est_end,
					     $genomic_source, $est_source, 
					     $genomic_strand, $est_strand, 
					     $genomic_primary, $est_primary);

      $self->output($pair);
    }    
  }
  close(EXONERATE);
}

=head2 run_ungapped

  Title   : run_ungapped
  Usage   : $self->run_ungapped()
  Function: Runs exonerate in ungapped mode (-n yes) and stores results as FeaturePairs
            Here we are only interested in having one FeaturePair per exon
  Returns : TRUE on success, FALSE on failure.
  Args    : none 

=cut

sub run_ungapped {
  my ($self) = @_;

  #check inputs
  my $genomicseq = $self->genomic_sequence ||    $self->throw("Genomic sequences not provided");
  my $estseq     = $self->est_sequence     ||    $self->throw("EST sequences not provided");
  
  #extract filenames from args and check/create files and directory
  my $genfile = $self->genfilename;
  my $estfile = $self->estfilename;
  
  # output parsing requires that we have both gff and cigar outputs. We don't want to do intron
  # prediction (ungapped) and we don't want to see alignments. Other options (wordsize, memory etc) 
  # are got from $self->arguments.

  my $command = $self->exonerate() . " " .  $self->arguments . " -n yes -A false -G yes --cdna $estfile --genomic $genfile";

  print STDERR "command is $command\n";

  # disable buffering of STDOUT
  STDOUT->autoflush(1);

  eval {
    open(EXONERATE, "$command |") or $self->throw("Error forking exonerate pipe on ".$self->genfilename."\n"); 

    my $source_tag  = "exonerate";
    my $prevstate   = "cigar";
    my $estname     = undef;
    my $genname     = "";

    my %gff;


    # exonerate output is split up by est.
    # gff comes first, then all the cigar lines, so we can process all the gff, then all the cigar lines.
    # if we keep track of when we flip from cigar to gff, we know we're moving on to the next est.

    while(<EXONERATE>){
      # gff
      if(/\S+\s+exonerate\s+exon\s+\d+/){

	# is this the first gff? ie have we moved onto a new est?
	if($prevstate eq 'cigar') {

	  $prevstate = 'gff';
	  $estname   = undef;

	  foreach my $entry( keys %gff) {
	    delete $gff{$entry};
	  }
	}

	my @cols = split /;/;
	my $pid  = $cols[$#cols];

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


	my ($cig, $hit, $query, $score, $cigar) = $_ =~ /^(cigar):\s+(\S+\s+\d+\s+\d+\s+[\+|-])\s+(\S+\s+\d+\s+\d+\s+[\+|-])\s+(\d+\.\d+)\s+(.+)/;
	next unless $cig eq 'cigar';

	my($qid, $qstart, $qend, $qstrand) = split / /, $query;
	my($hid, $hstart, $hend, $hstrand) = split / /, $hit;
	my $querystr = $qstart . "-" . $qend . "-" . $qstrand;
        my $pid = $gff{$querystr};
	if($qstrand eq "-"){
	  $qstrand = -1;
	}else{
	  $qstrand = 1;
	}
	if($hstrand eq "-"){
	  $hstrand = -1;
	}else{
	  $hstrand = 1;
	}
	my @cigars = split / /, $cigar;
	my $length = @cigars;
	if(!$length%2 == 1){
	  die("there are ".@cigars." elements in the cigar array very odd\n");
	}
	my $db_cigar;
	while(@cigars){
	  my $state = shift @cigars;
	  my $number = shift @cigars;
	  $db_cigar .= $number.$state;
	}

	$self->throw("$estname(gff) and $hid(cigar) do not match!") unless $estname eq $hid;
	my $pair = $self->_create_alignfeature($score, $qid, $qstart, $qend, $hid, $hstart, $hend, $qstrand, $hstrand, $pid, $db_cigar, 'exonerate', 'exonerate'); 
	if($self->print_results) {
	  print $pair->seqname  . "\t" . $pair->start        . "\t" . $pair->end    . "\t" . 
	        $pair->score    . "\t" . $pair->percent_id   . "\t" . $pair->strand . "\t" . 
	        $pair->hseqname . "\t" . $pair->hstart       . "\t" . $pair->hend   . "\t" . 
                $pair->hstrand  . "\t". $pair->cigar_string."\n";
	}
	$self->output($pair);
      }
     
    }
    
    
    close (EXONERATE) or $self->throw("Error running exonerate: $command");
    
  };  

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
    Returns :   none, but genes can be retrieved using each_gene
    Args    :   none

=cut

sub convert_output {
  my($self) = @_;

  my $curr_gene;
  my @genes;
  foreach my $fp($self->output){
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

  $self->add_genes(@genes);
}

=head2 add_genes

  Title   : add_genes
  Usage   : $self->add_genes
  Function: Adds an array of FeaturePairs representing genes, with exons as subFeatures, 
            into $self->{_'genelist'}. They can be retrieved using each_gene.
            NB this method is add_genes not add_Genes as these are not Gene objects ...
  Returns : nothing
  Args    : @genes: an array of Bio::EnsEMBL:FeaturePair

=cut

sub add_genes {

  my ($self, @genes) = @_;
    if (!defined($self->{'_genelist'})) {
      $self->{'_genelist'} = [];
  }

   if (@genes) {
    push(@{$self->{'_genelist'}},@genes);
  }
  
}

=head2 each_gene

  Title   : each_genes
  Usage   : $self->each_genes
  Function: Returns results of exonerate as an array of FeaturePair, one per gene. Exons
            are added as subSeqfeatures of genes, and supporting evidence as subSeqFeatures 
            of Exons.
            If you want just the raw Exonerate results as a set of FeaturePairs, use output 
            instead of each_gene
            NB this method is each_gene not each_Gene as these are not Gene objects ...
  Returns : An array of Bio::EnsEMBL::FeaturePair
  Args    : none

=cut

sub each_gene{
  my ($self) = @_;
      
  if (!defined($self->{'_genelist'})) {
    $self->{'_genelist'} = [];
  }
  
  return @{$self->{'_genelist'}};

}

=head2 output

  Title   : output
  Usage   : $self->output
  Function: Get/set method - can optionally add a FeaturePair to $self->{'_output'}
            Returns results of exonerate as array of FeaturePair, one per alignment
            If you want the results returned organised with one FeaturePair per "gene", 
            with associated exons and supporting features, use each_gene
  Returns : An array of Bio::EnsEMBL::FeaturePair
  Args    : either none, or $featpair: a Bio:EnsEMBL::FeaturePair

=cut

sub output {
  my ($self, $featpair) = @_;

  if (!defined($self->{'_output'})) {
    $self->{'_output'} = [];
  }

  if (defined $featpair) {
    push(@{$self->{'_output'}}, $featpair);
  }
  
  return @{$self->{'_output'}};

}

=head2 _create_featurepair

  Title   : _create_featurepair
  Usage   : $self->_create_featurepair($f1score, $f1pid, $f1start, $f1end, $f1id, 
				    $f2start, $f2end, $f2id, $f1source, 
				    $f2source, $f1strand, $f2strand, 
				    $f1primary, $f2primary)
  Function: Makes and returns a FeaturePair
  Returns : Bio::EnsEMBL::FeaturePair
  Args    : $f1score, $f1pid, $f1start, $f1end, $f1id, $f2start, $f2end, $f2id,
            $f1source, $f2source, $f1strand, $f2strand, $f1primary, $f2primary:
            strings and ints representing the various fileds to be filled in
            when creating a Bio::EnsEMBL::FeaturePair

=cut

sub _create_featurepair {
  my ($self, $f1score,  $f1id, $f1start, $f1end, $f2id, $f2start, $f2end,
      $f1source, $f2source, $f1strand, $f2strand, $f1primary, $f2primary, $f1pid) = @_;

  #create analysis object
  my $analysis_obj    = new Bio::EnsEMBL::Analysis
    (-db              => "none",
     -db_version      => "none",
     -program         => "exonerate",
     -program_version => "1",
     -gff_source      => $f1source,
     -gff_feature     => $f1primary,);
  
  $f1pid = 0 unless defined $f1pid;

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
#    push(@{$self->{'_fplist'}}, $fp);
  }
  return $fp;
}

sub _create_alignfeature {
  my ($self, $f1score,  $f1id, $f1start, $f1end, $f2id, $f2start, $f2end, $f1strand, $f2strand, $f1pid, $cigar, $gff_source, $gff_feature) = @_;
  
  #print "feature 1 ".$f1id." ".$f1start." ".$f1end." ".$f1strand." ".$f1pid." ".$f1score." \n";
  #print "feature 2 ".$f2id." ".$f2start." ".$f2end." ".$f2strand." ".$f1pid." ".$f1score." \n";
  #print "others ".$cigar." ".$gff_source." ".$gff_feature."\n";
  #create analysis object
  my $analysis_obj    = new Bio::EnsEMBL::Analysis
    (-db              => "none",
     -db_version      => "none",
     -program         => "exonerate",
     -program_version => "1",
     -gff_source      => $gff_source,
     -gff_feature     => $gff_feature);
  
  $f1pid = 0 unless defined $f1pid;

  #create features
  my $feat1 = new Bio::EnsEMBL::SeqFeature  (-start      =>   $f1start,
					     -end         =>   $f1end,
					     -seqname     =>   $f1id,
					     -strand      =>   $f1strand,
					     -score       =>   $f1score,
					     -percent_id  =>   $f1pid, 
					     -analysis    =>   $analysis_obj );
  
  my $feat2 = new Bio::EnsEMBL::SeqFeature  (-start       =>   $f2start,
					     -end         =>   $f2end,
					     -seqname     =>   $f2id,
					     -strand      =>   $f2strand,
					     -score       =>   $f1score,
					     -percent_id  =>   $f1pid, 
					     -analysis    =>   $analysis_obj );
  #create featurepair
  my $fp = new Bio::EnsEMBL::DnaDnaAlignFeature  (-feature1 => $feat1,
						  -feature2 => $feat2,
						  -cigar_string => $cigar,
						 );
  
  if ($fp) {
    $self->throw("Can't validate") unless $fp->validate();
#    push(@{$self->{'_fplist'}}, $fp);
  }
  return $fp;
}

#####################################
# clearing up temp files
#####################################


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
