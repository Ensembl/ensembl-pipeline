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

Bio::EnsEMBL::Pipeline::Runnable::FirstEF

=head1 SYNOPSIS

=head1 DESCRIPTION

=head1 CONTACT

Dan Andrews: dta@sanger.ac.uk

=head1 APPENDIX

=cut

package Bio::EnsEMBL::Pipeline::Runnable::FirstEF;

use vars qw(@ISA);
use strict;
use Bio::EnsEMBL::Pipeline::RunnableI;
use Bio::EnsEMBL::SimpleFeature;
use Bio::EnsEMBL::Analysis;
use Bio::Seq;
use Bio::SeqIO;
use Bio::EnsEMBL::Root;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI);

### Public methods

sub new {
  my ($class,@args) = @_;
  
  my $self = $class->SUPER::new(@args);    
  
  my( $query_seq,
      $repeatmask,
      $firstef_bin,
      $param_dir,
      $parse_script,
      $workdir) = $self->_rearrange([qw(QUERY
					REPEATMASKED
					FIRSTEF_BIN
					PARAM_DIR
					PARSE_SCRIPT
					WORK_DIR)],
				     @args);

  $self->_query_seq($query_seq)       if $query_seq;
  $self->_firstef_bin($firstef_bin)   if $firstef_bin;
  $self->_param_dir($param_dir)       if $param_dir;
  $self->_parse_script($parse_script) if $parse_script;
  $self->workdir($workdir)            if $workdir;

  if (defined $repeatmask){
    $self->_repeatmask($repeatmask);
  } else {
    $self->_repeatmask(1);
  }

  return $self
}

### Internal methods

sub DESTROY {
  my $self = shift;

  print "Cleaning up\n";

  $self->deletefiles;

}


sub _query_seq {
  my $self = shift;

  if (@_) {

    $self->{_query_seq} = shift;

    return unless $self->{_query_seq}; # Allows us to set to undef.

    $self->throw("Input isnt a Bio::SeqI or Bio::PrimarySeqI")
      unless ($self->{_query_seq}->isa("Bio::PrimarySeqI") || 
	      $self->{_query_seq}->isa("Bio::SeqI"));

    $self->filename($self->{_query_seq}->id . ".$$.seq");

    $self->results($self->filename . ".out");
  }

  return $self->{_query_seq}
}

sub _repeatmask {
  my $self = shift;

  if (@_) {
    $self->{_repeatmask} = shift;

    return unless $self->{_repeatmask}; # Allows us to set to undef.
  }

  return $self->{_repeatmask}
}

sub _write_seqs_for_firstef {
  my $self = shift;

  # Write our sequence to file - assuming only one sequence per 
  # firstef run.

  my $seqfile = $self->workdir . '/firstef_' . $$ . '.fa';

  my $seqio_out = Bio::SeqIO->new(-file   => ">$seqfile",
				  -format => 'fasta');

     # NOTE, this little backwater is where the repeatmasked 
     # sequence is plonked in, or not.
  if ($self->_repeatmask) {
print "USING REPEATMASKED SEQUENCE\n";
    $seqio_out->write_seq($self->_query_seq->get_repeatmasked_seq);
  } else {
print "NOT USING REPEATMASKED SEQUENCE\n";
    $seqio_out->write_seq($self->_query_seq)
  }

  $self->_seqfile($seqfile);

  $self->file($seqfile); # For final cleanup

  # Write the listfile that firstef likes to use.

  my $listfile = $self->workdir . '/firstef_listfile_' . $$;

  open(LISTFILE, ">$listfile") or die "Cant write to file [$listfile]";

  print LISTFILE "$seqfile   -1500\n";

  close(LISTFILE);

  $self->_listfile($listfile);

  $self->file($listfile); # For final cleanup
  $self->file($listfile ."_domain");
  $self->file($listfile ."_domain_comp");

  # While we are at it, generate the outfile name.
  # This filename is not something that the user can control.
  # It is the name of the listfile with _out appended.

  my $outfile = $listfile . '_out';

  $self->_outfile($outfile);

  $self->file($outfile); # For final cleanup

  return 1
}


sub _listfile {
  my $self = shift;

  if (@_){
    $self->{_listfile} = shift;
  }

  return $self->{_listfile}
}


sub _seqfile {
  my $self = shift;

  if (@_){
    $self->{_seqfile} = shift;
  }

  return $self->{_seqfile}
}


sub _outfile {
  my $self = shift;

  if (@_){
    $self->{_outfile} = shift;
  }

  return $self->{_outfile}
}

sub _parsed_outfile {
  my $self = shift;

  if (@_){
    $self->{_parsed_outfile} = shift;
    $self->file($self->{_parsed_outfile});
  }

  return $self->{_parsed_outfile}
}


sub _firstef_bin {
  my $self = shift;

  if (@_){
    $self->{_firstef_bin} = shift;

    $self->throw("firstef binary not executable [trying to execute " .
		 $self->{_firstef_bin} . "]")
      unless (-x $self->{_firstef_bin});
  }
print "Executable : " . $self->{_firstef_bin} . "\n";
  $self->throw("FirstEF executable has not been set.") 
    unless defined($self->{_firstef_bin});

  return $self->{_firstef_bin}
}

sub _parse_script {
  my $self = shift;

  if (@_){
    $self->{_parse_script} = shift;

    $self->throw("FirstEF parser script not executable [trying to execute " .
		 $self->{_parse_script} . "]")
      unless (-e $self->{_parse_script});
  }

  $self->throw("FirstEF parser script location has not been set.") 
    unless defined($self->{_parse_script});

  return $self->{_parse_script}
}

sub _param_dir {
  my $self = shift;

  if (@_) {

    $self->{_param_dir} = shift;

    $self->{_param_dir} .= '/' unless $self->{_param_dir} =~ /\/$/;

    return unless $self->{_param_dir}; # Allows us to set to undef.

    my @known_param_files = ('donor.3mer_wts_GChighDown',
                             'donor.3mer_wts_GClowDown',
                             'donor.6mer_wts_GChighDown',
                             'donor.6mer_wts_GChighUp',
                             'donor.6mer_wts_GClowDown',
                             'donor.6mer_wts_GClowUp',
                             'donor.decisiontree',
                             'donor.decisiontree.orig',
                             'donor.qdamodel.GChigh',
                             'donor.qdamodel.GClow',
                             'exon.qdamodel.CpGpoor_GChigh',
                             'exon.qdamodel.CpGpoor_GClow',
                             'exon.qdamodel.CpGrich_GChigh',
                             'exon.qdamodel.CpGrich_GClow',
                             'promoter.5mer_wts_CpGpoor_430.510',
                             'promoter.5mer_wts_CpGpoor_490.570',
                             'promoter.5mer_wts_CpGrich_430.510',
                             'promoter.5mer_wts_CpGrich_490.570',
                             'promoter.6mer_wts_CpGpoor_1.250',
                             'promoter.6mer_wts_CpGpoor_1.450',
                             'promoter.6mer_wts_CpGpoor_200.450',
                             'promoter.6mer_wts_CpGrich_1.250',
                             'promoter.6mer_wts_CpGrich_1.450',
                             'promoter.6mer_wts_CpGrich_200.450',
                             'promoter.qdamodel.CpGpoor',
                             'promoter.qdamodel.CpGrich');

    my @missing_files;

    foreach my $param_file (@known_param_files) {

      unless (-e $self->{_param_dir} . "/$param_file"){
	push (@missing_files, $param_file);

	print "Missing parameter file : " . $self->{_param_dir} . "/$param_file\n";
      }
    }

    $self->throw("The above parameter files are missing.") if @missing_files;

  }

  return $self->{_param_dir}
}



###########
# Analysis methods
##########


sub run {
  my ($self) = @_;

  # Check sequence exists
  $self->throw("About to run firstef, but dont yet have a sequence.") 
    unless $self->_query_seq;

  # Check work directory
  $self->throw("The working directory is not OK - probably not enough disk space.")
    unless $self->checkdir;

  # Put our sequence in the right places and generate input/output files.
  $self->_write_seqs_for_firstef;

  # Execute firstef
  $self->run_firstef;

  # Parse output
  $self->parse_results;

  # Clean up
  $self->deletefiles;

  return  1
}


sub run_firstef {
  my $self = shift;

  # firstef insists on writing its output file to the same dir as the
  # listfile - it seems.

  my $command = $self->_firstef_bin . ' 1500 ' . $self->_listfile . ' ' . 
    $self->_param_dir . ' 0 0.4 0.4 0.5';

  print $command . "\n";

  $self->throw("A fatal error was encountered while running firstef.") 
    if system($command);

  return 1
}


sub parse_results {
  my $self = shift;

  # The output of FirstEF first has to be parsed by a script that accompanies
  # the firstef binary.  This parser determines the best predicted exons.
  # Hence, we have to run this parser then parse the output from the file
  # it generates.  This is messy, but it turned out to be preferable to 
  # importing (and no doubt maintaining) the external parser.

  # Run the external parser.
  my $parsed_file = $self->_first_parse;

  # Parse the output of the external parser.

  open (OUTPUT, $parsed_file) 
    or die "Couldnt open file produced by external firstef parser [$parsed_file]";

  my $strand;

  while (<OUTPUT>) {

    $strand = 1 if (/direct strand/);
    $strand = -1 if (/complementary strand/);


    # E.G.:
    #No.     Promoter     P(promoter)     Exon          P(exon) P(donor)  CpG Window        Rank
    #1  00023223..00023792  1.0000  00023723..00023784  1.0000  0.9940  00023779..00023980  	1
    #1  00023425..00023994  1.0000  00023925..00024627  1.0000  0.9863  00023779..00023980  	2


    if (/\d+\s+\S+\s+\S+([^\.]+)\.\.(\d+)\s+(\S+)\s+\S+\s+\S+\s+(\d+)/) {

      my %feature;

      my ($start, $end) = sort {$a <=> $b} ($1 * 1 , $2 * 1);

      $feature{name}            = 'firstef_exon';
      $feature{score}           = $3;
      $feature{start}           = $start;
      $feature{end}             = $end;
      $feature{strand}          = $strand;
      $feature{program}         = 'firstef';
      $feature{program_version} = '1';
      $feature{display_label}   = "rank = $4";

      $self->create_feature(\%feature);
    }

  }

  close (OUTPUT);

}


sub _first_parse {
  my ($self) = @_;

  # An essential part of the firstef run is a somewhat messy post-
  # parsing of the output.  This is done in a separate external
  # script called FirstEF_parser.pl.

  my $parsed_output = $self->workdir . '/first_parse.' . $$;
  $self->_parsed_outfile($parsed_output);

  my $command = $self->_parse_script . ' ' . $self->_outfile . ' ' . 
    $self->_parsed_outfile;

  $self->throw("Problem parsing output with external parser [" . 
	       $self->_parse_script . "]")
    if (system($command));

  return $self->_parsed_outfile
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

  my $first_exon = Bio::EnsEMBL::SimpleFeature->new(-seqname =>  $feat->{name},
						    -start   =>  $feat->{start},
						    -end     =>  $feat->{end},
						    -strand  =>  $feat->{strand},
						    -score   =>  $feat->{score});

  $first_exon->display_label($feat->{'display_label'});

  $first_exon->validate();

  $self->output($first_exon);

  return 1
}
