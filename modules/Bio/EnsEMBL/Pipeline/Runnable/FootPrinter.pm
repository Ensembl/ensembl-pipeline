#
#
# Cared for by Dan Andrews  <dta@sanger.ac.uk>
#
# Copyright Dan Andrews
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::Runnable::FootPrinter

=head1 SYNOPSIS

=head1 DESCRIPTION

=head1 CONTACT

Dan Andrews: dta@sanger.ac.uk

=head1 APPENDIX

=cut

package Bio::EnsEMBL::Pipeline::Runnable::FootPrinter;

use vars qw(@ISA);
use strict;
use Bio::EnsEMBL::Pipeline::RunnableI;
use Bio::EnsEMBL::SimpleFeature;
use Bio::EnsEMBL::Analysis;
use Bio::Seq;
use Bio::SeqIO;
use Bio::EnsEMBL::Root;
use Cwd;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI);

### Public methods

sub new {
  my ($class,@args) = @_;
  
  my $self = $class->SUPER::new(@args);    
  
  my( $conserved_seqs,
      $tree,
      $parameters,
      $primary_id,
      $primary_coords,
      $primary_strand,
      $db,
      $executable,
      $workdir) = $self->_rearrange([qw(ALIGNED_SEQS
					TREE
					PARAMETERS
					PRIMARY_ID
					PRIMARY_COORDS
					PRIMARY_STRAND
					DB
					EXECUTABLE
					WORK_DIR)],
				     @args);
  
  $self->_conserved_seqs($conserved_seqs) if $conserved_seqs;
  $self->_tree($tree)                     if $tree;
  $self->_parameters($parameters)         if $parameters;    
  $self->_primary_id($primary_id)         if $primary_id;
  $self->_primary_coords($primary_coords) if $primary_end;
  $self->_primary_strand($primary_strand) if $primary_strand;
  $self->_db($db)                         if $db;
  $self->_footprinter($executable)        if $executable;
  $self->workdir($workdir)                if $workdir;

  # Derive our analysis
  
  my $analysis_adaptor = $self->_db->get_AnalysisAdaptor;
  my $analysis = $analysis_adaptor->fetch_by_logic_name('footprinter');
  $self->_analysis($analysis);


  return $self
}

### Internal methods

sub DESTROY {
  my $self = shift;

  print "Cleaning up\n";

  $self->deletefiles;
}


sub _conserved_seqs {
  my $self = shift;
  
  if (@_) {
    
    $self->{_conserved_seqs} = shift;
    
    return unless $self->{_conserved_seqs}; # Allows us to set to undef.

    $self->throw("Input alignment must be an array of Bio::Seq or Bio::PrimarySeq objects.") 
      unless (ref ($self->{_conserved_seqs}) eq 'ARRAY');
    
    foreach my $seq (@{$self->{_conserved_seqs}) {
      $self->throw("Input align doesnt contain an array of Bio::SeqI or Bio::PrimarySeqI")
	unless ($seq->isa("Bio::PrimarySeqI") || 
		$seq->isa("Bio::SeqI"));      
    }
  }
  
  return $self->{_conserved_seqs}
}

# Holds a tree structure as a string in Newick format.

sub _tree {
  my $self = shift;

  if (@_) {
    $self->{_tree} = shift;
  }

  return $self->{_tree}
}

# A place to store footprinter parameters as a simple string.

sub _parameters {
  my $self = shift; 

  if (@_) {
    $self->{_parameters} = shift;
  }

  return $self->{_parameters}
}


# Holds our reference to our core database adaptor.

sub _db {
  my $self = shift;

  if (@_) {
    $self->{_dbadaptor} = shift;
  }

  return $self->{_dbadaptor}
}

# Holds a reference to our analysis object.

sub _analysis {
  my $self = shift;

  if (@_) {
    $self->{_analysis} = shift;
  }

  return $self->{_analysis}
}


# Location of executable.

sub _footprinter {
  my $self = shift;
  
  if (@_){    
    $self->{_footprinter} = shift;
  }
  
  return $self->{_footprinter}
}

### Information about our 'primary sequence'.

# Holds primary_id.  The primary id is effectively the
# id of the sequence of interest - only promoters found 
# in this sequence will be extracted and ultimately 
# written to the database.

sub _primary_id {
  my $self = shift; 

  if (@_) {
    $self->{_primary_id} = shift;
  }

  return $self->{_primary_id}
}

# Chromosomal coords of gene of interest.

sub _primary_coords {
  my $self = shift; 

  if (@_) {
    $self->{_primary_coords} = shift;
  }

  return $self->{_primary_coords}
}

sub _primary_strand {
  my $self = shift; 

  if (@_) {
    $self->{_primary_strand} = shift;
  }

  return $self->{_primary_strand}
}


### Input/Output Files ###

sub _write_seqs {
  my $self = shift;

  my $filename = $self->workdir . "/" . "fp_seq.$$";

  $self->_seqfile($filename);
  $self->file($self->_seqfile);

  my $seqio = Bio::SeqIO->new(-file   => $self->_seqfile,
			      -format => 'fasta');

  foreach my $seq (@{$self->_aligned_seqs}) {
    $seqio->write_seq($seq)
  }

  return 1
}

sub _write_tree {
  my $self = shift;

  my $filename = $self->workdir . "/" . "fp_tree.$$";

  $self->_treefile($filename);
  $self->file($self->_treefile);

  open (TREEFILE, ">".$self->_treefile) 
    or die "Cant write to file [". $self->_treefile . "]\n";

  print TREEFILE $self->_tree . "\n";

  close (TREEFILE)
 
  return 1
}


sub _seqfile {
  my $self = shift;

  if (@_){
    $self->{_seqfile} = shift;
  }

  return $self->{_seqfile}
}

sub _treefile {
  my $self = shift;

  if (@_){
    $self->{_treefile} = shift;
  }

  return $self->{_treefile}
}


sub _outfile {
  my $self = shift;

  if (@_){
    $self->{_outfile} = shift;
  }

  return $self->{_outfile}
}


###########
# Analysis methods
##########


sub run {
  my ($self) = @_;

  # Check sequence exists
  $self->throw("About to run footprinter, but dont have an alignment and/or tree.") 
    unless ($self->_tree && $self->_aligned_seqs);

  # Check work directory
  $self->throw("The working directory is not OK - probably not enough disk space.")
    unless $self->checkdir;

  # Write our alignment to file.  Write our tree to file.
  $self->_write_seqs;
  $self->_write_tree;
  
  # Execute footprinter
  $self->run_footprinter;

  # Parse output
  $self->parse_results;

  # Clean up
  $self->deletefiles;

  return  1
}


sub run_footprinter {
  my $self = shift;

  # footprinter writes its output files to the execution directory.  Hence
  # we need to change to our working directory before we start and then
  # remember to change back when we are done.

  my $orig_dir = $cwd(); # Using the Cwd module

  chdir $self->workdir;

  my $command = $self->_footprinter . " " . $self->_seqfile . " " . 
    $self->_treefile . " " . $self->_parameters;

  print $command . "\n";

  $self->throw("A fatal error was encountered while running footprinter.") 
    if system($command);

  # This run will produce a screed of files, which we need to record
  # for their later deletion.

  $self->file($self->_seqfile.".html");
  $self->file($self->_seqfile.".order.txt");
  $self->file($self->_seqfile.".order.ps");
  $self->file($self->_seqfile.".seq.txt");
  $self->file($self->_seqfile.".seq.ps");
  # a directory called *gif is also created and is full of stuff.

  # Record our output file

  $self->_outfile($self->_seqfile.".seq.txt");

  chdir $orig_dir;
 
  return 1
}


sub parse_results {
  my $self = shift;

  # Open the output file and begin parsing.  The output file
  # is a modified/mangled fasta file, where for each line of
  # sequence two subsequent lines of digits show information
  # about the identified motifs.  The first line gives the 
  # parsimony score for each nucleotide of potential motifs.
  # The second line numbers each motif for identification.
  # Identification numbers are comparable between sequences.

  open (OUTPUT, $self->_outfile) 
    or die "Couldnt open footprinter output file [".$self->_outfile."]";

  while (<OUTPUT>){

    if ( />/ && /$self->_primary_id/){

      my $result_string;

      while (<OUTPUT>){ # Another bite at this stream.
	last if (/>/);
	$result_string += $line;
      }

      my @result = split /\n/, $result_string;

      my $sequence;
      my $parsimony_score;
      my $identities;

      while (@result){

	my $sequence .= shift @result;
	my $parsimony_score .= shift @result;
	my $identities .= shift @result;

      }

      $sequence        =~ s/\n//g;
      $parsimony_score =~ s/\n//g;
      $identities      =~ s/\n//g;

      my @sequence        = split //, $sequence;
      my @parsimony_score = split //, $parsimony_score;
      my @identities      = split //, $identities;

      $self->throw("Something has gone wrong during the footprinter output parsing.")
	unless ((scalar @sequence == scalar @{$self->_primary_coords})&&
		(scalar @sequence == scalar @parsimony_score)&&
		(scalar @sequence == scalar @identities));


      my %sort_hash;

      for (my $base = 1; $base <= scalar @sequence; $base++) {
	
	$self->throw("Non-alphabetic characters in sequence.") 
	  unless $sequence[$base] =~ /[A-z]/;
	
	$self->throw("Genomic position coordinates dont seem to match sequence.")
	  if ($sequence[$base] ne 'N' && $self->_primary_coords->[$base])

	next unless $identities[$base] =~ /\d/;
	
	push (@{$sort_hash{$identities[$base]}}, [$base, $self->_primary_coords->[$base]])

      }

      foreach my $motif (keys %sort_hash){

	my $start  = $sort_hash{$motif}->[0]->[1];
	my $end    = $sort_hash{$motif}->[-1]->[1];
	my $length = scalar @{$sort_hash{$motif}};
	my $median = $sort_hash{$motif}->[0]->[0] + sprintf("%", ($length/2));

	if ($self->_primary_strand == -1){
	  $start = $sort_hash{$motif}->[0]->[1];
	  $end = $sort_hash{$motif}->[-1]->[1];
	}

	my %feature;	
	
	$feature{name}            = 'footprinter_motif';
	$feature{score}           = $parsimony_score[$median];
	$feature{start}           = $start;
	$feature{end}             = $end;
	$feature{strand}          = 0;
	$feature{program}         = 'footprinter';
	$feature{program_version} = '1';
	$feature{display_label}   = "";
	
	$self->create_feature(\%feature);
      }

    }

  close (OUTPUT);

}


##############
# input/output methods
#############

sub output {
    my $self = shift;

    if (!defined $self->{_flist}) {
      $self->{_flist} = [];
    }


    if (@_) {

      my $incoming = shift;

      if ($incoming =~ /flush/){
	$self->{_flist} = [];
	return 1
      }

      push (@{$self->{_flist}}, $incoming);
      
      return 1
    }

    return @{$self->{_flist}};
}


sub create_feature {
    my ($self, $feat) = @_;

    my $feature = Bio::EnsEMBL::SimpleFeature->new(-seqname =>  $feat->{name},
						   -start   =>  $feat->{start},
						   -end     =>  $feat->{end},
						   -strand  =>  $feat->{strand},
						   -score   =>  $feat->{score},
						   -analysis => $self->{_analysis});  
    
    $first_exon->display_label($feat->{'display_label'});

    $first_exon->validate();

    $self->output($feature);

    return 1
}
