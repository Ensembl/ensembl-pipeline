#
 Trivial edit#
# Cared for by Dan Andrews  <dta@sanger.ac.uk>
#
# Copyright Dan Andrews
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::Runnable::Motifwise

=head1 SYNOPSIS

=head1 DESCRIPTION

=head1 CONTACT

Dan Andrews: dta@sanger.ac.uk

=head1 APPENDIX

=cut

package Bio::EnsEMBL::Pipeline::Runnable::Motifwise;

use vars qw(@ISA);
use strict;
use Bio::EnsEMBL::Pipeline::RunnableI;
use Bio::EnsEMBL::SimpleFeature;
use Bio::EnsEMBL::Analysis;
use Bio::Seq;
use Bio::SeqIO;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI);

### Public methods

sub new {
  my ($class,@args) = @_;
  
  my $self = $class->SUPER::new(@args);    
  
  my( $query_seq,
      $executable,
      $parameters,
      $motif_file,
      $workdir) = $self->_rearrange([qw(QUERY_SEQ
					EXECUTABLE
					PARAMETERS
					MOTIF_FILE
					WORK_DIR)],
				    @args);

  $self->query_seq($query_seq)   if defined($query_seq);
  $self->workdir($workdir)       if defined($workdir);
  $self->executable($executable) if defined($executable);
  $self->parameters($parameters) if defined($parameters);
  $self->motif_file($motif_file) if defined($motif_file);

  return $self
}

sub query_seq {
  my $self = shift;

  if (@_) {
   
    $self->{_query_seq} = shift;

    return unless $self->{_query_seq}; # Allows us to set to undef.

    $self->throw("Input isnt a Bio::Seq or Bio::PrimarySeqI")
      unless ($self->{_query_seq} && 
	      ($self->{_query_seq}->isa("Bio::Seq") ||
	       $self->{_query_seq}->isa("Bio::PrimarySeqI")));
  }

  $self->throw("There is currently no defined query sequence.")
    unless $self->{_query_seq};

  return $self->{_query_seq}
}

sub executable {
  my $self = shift;

  if (@_) {
    $self->{_executable} = shift;
  }

  $self->{_executable} = 'motifwise'
    unless $self->{_executable}; # Hoping motifwise is in our path.

  return $self->{_executable}
}

sub parameters {
  my $self = shift;

  if (@_) {
    $self->{_parameters} = shift;
  }

  return $self->{_parameters}
}


sub motif_file {
  my $self = shift;

  if (@_) {
    $self->{_motif_file} = shift;

    $self->throw("Motif file does not exist.")
      unless ($self->{_motif_file} && -e $self->{_motif_file});
  }

  $self->throw("Motif file has not been set.")
    unless $self->{_motif_file};

  return $self->{_motif_file}
}


sub run {
  my ($self) = @_;

  # Clean up files from any previous runs.
  $self->_flush;

  # Check work directory
  $self->throw("The working directory is not OK - probably not enough disk space.")
    unless $self->checkdir;

  # Put our sequence in the right places and generate input/output files.
  $self->_write_seq;

  # Execute motifwise
  $self->_run_motifwise;

  # Parse output
  $self->parse_results;

  return  1
}


sub parse_results {
  my $self = shift;

  open (OUTPUT, $self->_outfile) 
    or die "Couldnt open motifwise results file [".$self->_outfile."]";

  while (<OUTPUT>) {

    # E.G.:

    # Region	Contig:AC008066.4.1.174347	100326	100403	0.51	21.71
    # motif
    # Motif	Contig:AC008066.4.1.174347	100326	100336	1	motif_4	11.60	CTGCGCGGGCC
    # Motif	Contig:AC008066.4.1.174347	100326	100336	-1	motif_52	12.06	CTGCGCGGGCC
    # Motif	Contig:AC008066.4.1.174347	100334	100341	1	motif_5	11.91	GCCCCACC
    # Motif	Contig:AC008066.4.1.174347	100375	100383	1	motif_35	14.05	GGCTCCGCC
    # Motif	Contig:AC008066.4.1.174347	100377	100387	-1	motif_41	16.84	CTCCGCCCCCG
    # Motif	Contig:AC008066.4.1.174347	100377	100386	1	motif_3	11.50	CTCCGCCCCC
    # Motif	Contig:AC008066.4.1.174347	100393	100403	1	motif_4	12.06	GTGCGCAAGCG
    # Motif	Contig:AC008066.4.1.174347	100395	100404	1	motif_26	12.81	GCGCAAGCGC
    # Motif	Contig:AC008066.4.1.174347	100395	100404	-1	motif_26	12.81	GCGCAAGCGC
    # end motif
    # >Contig:AC008066.4.1.174347
    # CTGCGCGGGCCCCACCTACGGGCTTCGAGCTTCCACGTGCAGGGCCTCTGGCTCCGCCCC
    # CGGCGCAGTGCGCAAGCG
    # end region    

    if (/^Region\s+/){
      /Region\t[^\t]+\t(\d+)\t(\d+)\t[\d\.]+\t([\d\.]+)/;

      $self->throw("Cant properly parse motifwise (region) output.")
	unless ($1 && $2 && $3);

      my ($start, $end) = sort {$a <=> $b} ($1 * 1 , $2 * 1);

      my %tss_simplefeat;

      $tss_simplefeat{name}            = 'motifwise_tss';
      $tss_simplefeat{score}           = $3;
      $tss_simplefeat{start}           = $start;
      $tss_simplefeat{end}             = $end;
      $tss_simplefeat{strand}          = 0;
      $tss_simplefeat{display_label}   = "";

      $self->_create_feature(\%tss_simplefeat);

      next;
    }

    if (/^Motif\s+/){
      /Motif\t[^\t]+\t(\d+)\t(\d+)\t([1\-]+)\t([\w\_]+)\t([\d\.]+)/;

     $self->throw("Cant properly parse motifwise (motif) output.")
	unless ($1 && $2 && $3 && $4 && $5);

      my ($start, $end) = sort {$a <=> $b} ($1 * 1 , $2 * 1);

      my %motif_simplefeat;

      $motif_simplefeat{name}            = 'motifwise_motif';
      $motif_simplefeat{score}           = $5;
      $motif_simplefeat{start}           = $start;
      $motif_simplefeat{end}             = $end;
      $motif_simplefeat{strand}          = $3;
      $motif_simplefeat{display_label}   = $4;

      $self->_create_feature(\%motif_simplefeat);

      next;
    }

  }

  close (OUTPUT);
}

sub parse_file {
  my ($self, $filename) = shift;

  $self->throw("Can\'t open file for parsing [file:$filename]\n")
    unless (-e $filename);

  $self->_outfile($filename);

  $self->parse_results;
}


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




### Internal methods

sub DESTROY {
  my $self = shift;

  print "Cleaning up\n";

  $self->deletefiles;

}


sub _write_seq {
  my $self = shift;

  # Check that a sequence exists.

  $self->throw("No sequence has yet been passed to Motifwise runnable.")
    unless $self->query_seq;

  # Create/open a sequence file.

  my $seqfile = $self->workdir . '/motifwise_' . $$ . '.fa';

  my $seqio_out = Bio::SeqIO->new(-file   => ">$seqfile",
				  -format => 'fasta');

  # Write query sequence to file.

  $seqio_out->write_seq($self->query_seq);

  # Store filenames for housekeeping purposes.

  $self->_seqfile($seqfile);

  $self->file($seqfile); # For final cleanup

  return 1
}


sub _run_motifwise {
  my $self = shift;

  my $command = $self->executable . ' -lr ' . $self->motif_file . ' ' . 
    $self->parameters . ' ' . $self->_seqfile . ' > ' . $self->_outfile;

  print $command . "\n";

  $self->throw("A fatal error was encountered while running motifwise.") 
    if system($command);

  return 1
}


sub _create_feature {
    my ($self, $feat) = @_;

    my $simple_feat = Bio::EnsEMBL::SimpleFeature->new(-seqname => $feat->{name},
						       -start   => $feat->{start},
						       -end     => $feat->{end},
						       -strand  => $feat->{strand},
						       -score   => $feat->{score});  
    
    $simple_feat->display_label($feat->{'display_label'});

    $simple_feat->validate();

    $self->output($simple_feat);

    return 1
}


sub _flush {
  my $self = shift;

  unlink $self->{_outfile} if $self->{_outfile};
  unlink $self->{_seqfile} if $self->{_seqfile};

  $self->{_outfile} = '';
  $self->{_seqfile} = '';
  $self->output('flush');

  return 1
}


### Storage/retrieval

sub _seqfile {
  my $self = shift;

  if (@_){
    $self->{_seqfile} = shift;
    $self->file($self->{_seqfile});
  }

  $self->throw("Must not run until a sequence file has been set.")
    unless $self->{_seqfile};

  return $self->{_seqfile}
}

sub _outfile {
  my $self = shift;

  unless ($self->{_outfile}){
    $self->{_outfile} = $self->workdir . '/motifwise_' . $$ . '.out';
    $self->file($self->{_outfile});
  }

  return $self->{_outfile}
}

return 1;

