package Bio::EnsEMBL::Pipeline::Runnable::MinimalBlast;

# Blasts a sequence against a pre-formatted database.

# This is a low overhead blast module.  It assumes you know what kind of blast
# you want and expects you know the command line parameters that you want
# set.  If this is the case, this module should allow you to run blast jobs with
# a minimum of pre-configuration and will allow you total control over how your job
# is handled.

# Inputs are : query sequence(s) and a blastDB object.
# Output is : a BPLite parser object.


use vars qw(@ISA);
use strict;
use Bio::EnsEMBL::Root;
use Bio::Tools::BPlite;
use Bio::SeqIO;

@ISA = qw(Bio::EnsEMBL::Root);


sub new {
  my ($class, @args) = @_;

  my $self = bless {},$class;

  my ($program,
      $path,
      $blastdb,
      $queryseq,
      $options,
      $workdir,
      $identity_cutoff) = 
	$self->_rearrange([qw(PROGRAM
			      PATH
			      BLASTDB
			      QUERYSEQ
			      OPTIONS
			      WORKDIR
			      IDENTITY_CUTOFF
			     )],@args);

  $self->path($path)         if defined($path);
  $self->program($program)   if defined($program);
  $self->blastdb($blastdb)   if defined($blastdb);
  $self->queryseq($queryseq) if defined($queryseq);

  $self->workdir($workdir) if defined($workdir);
  $self->identity($identity_cutoff) if defined($identity_cutoff);

  return $self
}


sub DESTROY {
  my $self = shift;

  unlink $self->output_file  if $self->output_file;
  rmdir  $self->workdir      if $self->workdir;
}

sub run {
  my ($self) = @_;

  my $tag = time;

  my $input_file  = $self->workdir . 'blast_input.'  . $tag;
  my $output_file = $self->workdir . 'blast_output.' . $tag;

  my $seqio = Bio::SeqIO->new(-file   => ">$input_file",
			      -format => 'fasta');

  $seqio->write_seq($self->queryseq);

  my $command;

  if ($self->blastdb->index_type =~ /wu/) {

    $command = $self->program . ' ' . $self->blastdb->dbfile . ' ' 
      . $input_file . ' ' . $self->options . ' > ' . $output_file;

  } elsif ($self->blastdb->index_type eq 'ncbi') {

    $command = 'blastall -p ' . $self->program . ' -d ' 
      . $self->blastdb->dbfile . ' -i ' . $input_file . ' ' 
	. $self->options . ' > ' . $output_file;
    
  }
  print $command . "\n";

  my $exit_status = system($command);

  if ($exit_status) {
    $self->warn("Blast exited with a non-zero status.\n".
		"Execution command was [$command].");
    return 0
  }

  $self->output_file($output_file);
  unlink $input_file;

  my $blast_parser = Bio::Tools::BPlite->new(-file => $self->output_file);

  return $blast_parser
}

### Getters/Setters ###

sub output_file {
  my $self = shift;

  if (@_) {
    $self->{_output_file} = shift;
  }

  return $self->{_output_file}
}

sub path {
  my $self = shift;

  if (@_) {
    $self->{_path} = shift;

    $self->throw("Directory specified in the path does not\n" .
		 "exist [".$self->{_path}."].")
      unless -d $self->{_path};

    # Append a backslash if the path string doesnt already
    # terminate with one.
    $self->{_path} .= '/' 
      unless $self->{_path} =~ /\/$/;
  }

  return $self->{_program}
}

sub program {
  my $self = shift;

  if (@_) {
    $self->{_program} = shift;
  }

  if (! defined $self->{_program}){
    $self->throw("Blast program type not specified " .
		 "(eg. blastn, wublastp).");
  }

  return $self->{_program}
}

sub blastdb {
  my $self = shift;

  if (@_) {
    $self->{_blastdb} = shift;
  }

  return $self->{_blastdb}
}

sub queryseq {
  my $self = shift;

  if (@_) {
    $self->{_queryseq} = shift;
  }

  return $self->{_queryseq}
}

sub options {
  my $self = shift;

  if (@_) {
    $self->{_options} = shift;
  }

  $self->{_options} = ' ' unless $self->{_options};

  return $self->{_options}
}

sub identity {
  my $self = shift;

  if (@_) {
    $self->{_identity} = shift;
    $self->{_identity} /= 100 if $self->{_identity} > 1;
  }

  return $self->{_identity}
}

sub workdir {
  my $self = shift;

  if (@_ && $self->{_workdir}&& -d $self->{_workdir}) {
    $self->throw("Blast working directory already set")
  }

  if(!@_ && !$self->{_workdir}){
    $self->{_workdir} = '/tmp/';
  }

  if (@_) {
    $self->{_workdir} = shift;
    $self->{_workdir} .= '/';
  }

  if ($self->{_workdir} !~ /tempblast/) {
    my $tag = time;
    my $dir = $self->{_workdir} . "\/tempblast." . $tag . '/';

    mkdir $self->{_workdir};
    mkdir $dir;

    $self->throw("Failed to create a temporary working directory.") 
      unless (-d $dir);
    
    $self->{_workdir} = $dir;

  }

  return $self->{_workdir}
}

return 1;
