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
      $db,
      $firstef_dir,
      $param_dir,
      $workdir) = $self->_rearrange([qw(QUERY
					 DB
					 FIRSTEF_DIR
					 PARAM_DIR
					 WORK_DIR)],
				     @args);
  
  $self->_query_seq($query_seq)     if $query_seq;
  $self->_db($db)                   if $db;
  $self->_firstef_dir($firstef_dir) if $firstef_dir;
  $self->_param_dir($param_dir)     if $param_dir;
  $self->workdir($workdir)          if $workdir;

  # Derive our analysis
  
  my $analysis_adaptor = $self->_db->get_AnalysisAdaptor;
  my $analysis = $analysis_adaptor->fetch_by_logic_name('firstef');
  $self->_analysis($analysis);


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

sub _db {
  my $self = shift;

  if (@_) {
    $self->{_dbadaptor} = shift;
  }

  return $self->{_dbadaptor}
}

sub _analysis {
  my $self = shift;

  if (@_) {
    $self->{_analysis} = shift;
  }

  return $self->{_analysis}
}


sub _write_seqs_for_firstef {
  my $self = shift;

  # Write our sequence to file - assuming only one sequence per 
  # firstef run.

  my $seqfile = $self->workdir . '/firstef_' . $$ . '.fa';

  my $seqio_out = Bio::SeqIO->new(-file   => ">$seqfile",
				  -format => 'fasta');
  
  $seqio_out->write_seq($self->_query_seq);

  $self->_seqfile($seqfile);

  $self->file($seqfile); # For final cleanup

  # Write the listfile that firstef likes to use.

  my $listfile = $self->workdir . '/firstef_listfile_' . $$;

  open(LISTFILE, ">$listfile") or die "Cant write to file [$listfile]";

  print LISTFILE "$seqfile   -1500\n";

  close(LISTFILE);

  $self->_listfile($listfile);

  $self->file($listfile); # For final cleanup

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
  }

  return $self->{_parsed_outfile}
}


sub _firstef_dir {
  my $self = shift;
  
  if (@_){    
    $self->{_firstef_dir} = shift;
  }
  
  return $self->{_firstef_dir}
}


sub _firstef {
  my $self = shift;
  
  unless ($self->{_firstef}){
    
    ### Fiddling around to get the correct executable for the operating system

    open(HOSTTYPE, "echo \$HOSTTYPE |") 
      or die "Problem determining the runtime hosttype";

    my $hosttype = <HOSTTYPE>;
    chomp $hosttype;
    
    close (HOSTTYPE);

    $self->{_firstef} = $self->_firstef_dir . "/firstef.$hosttype";
    
    $self->throw("firstef not found at " . $self->{_firstef} . ": $!\n" .
		 "Please provide just the _directory_ where the executables are\n" .
		 "to be found.  The actual executable is determined according\n" .
		 "to the host operating system at runtime") 
      unless (-e $self->{_firstef});
    
    $self->throw("firstef binary not executable [trying to execute " .  
		 $self->{_firstef} . "]") 
      unless (-x $self->{_firstef});
  }
  
  return $self->{_firstef}
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

  my $command = $self->_firstef . ' 1500 ' . $self->_listfile . ' ' . 
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

  my $command = $self->_firstef_dir . '/FirstEF_parser.pl ' . 
    $self->_outfile . ' ' . $self->_parsed_outfile;

  $self->throw("Problem parsing output with external parser [" . 
	       $self->_firstef_dir . "/FirstEF_parser.pl]")
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
						      -score   =>  $feat->{score},
						      -analysis => $self->{_analysis});  
    
    $first_exon->display_label($feat->{'display_label'});

    $first_exon->validate();

    $self->output($first_exon);

    return 1
}
