
# To do:
# ------
# -Currently only implements codeml wrapper.  Should at least
#   include baseml too.
# -Some perfectly acceptable PAML options set to 0 are being 
#   treated by Perl as 'unset'.  This probably is harmless, 
#   but is a bug nonetheless.
# -Set/Get functions should check dependencies such that config 
#   errors are thrown by the perl layer rather than the paml 
#   application.  Tough.
# -Set functions should check that the set values are allowed.

package Bio::EnsEMBL::Pipeline::GeneDuplication::PAML;

use vars qw(@ISA);
use strict;
use Bio::EnsEMBL::Pipeline::RunnableI;
use Bio::Tools::Phylo::PAML;
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Exception qw(warning throw);
use Cwd;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI);

sub new {
  my ($class, @args) = @_;

  my $self = bless {},$class;

  my ($executable,
      $work_dir,
      $aligned_seqs,
      $seqfile,
      $treefile,
      $outfile,
      $noisy,
      $verbose,
      $runmode,
      $seqtype,
      $codonfreq,
      $aadist,
      $aaratefile,
      $model,
      $nssites,
      $icode,
      $mgene,
      $fix_kappa,
      $kappa,
      $fix_omega,
      $omega,
      $fix_alpha,
      $alpha,
      $malpha,
      $ncatg,
      $clock,
      $getse,
      $rateancestor,
      $small_diff) = rearrange([qw(EXECUTABLE
				   WORK_DIR
				   ALIGNED_SEQS
				   SEQFILE
				   TREEFILE
				   OUTFILE
				   NOISY
				   VERBOSE
				   RUNMODE
				   SEQTYPE
				   CODONFREQ
				   AADIST
				   AARATEFILE
				   MODEL
				   NSSITES
				   ICODE
				   MGENE
				   FIX_KAPPA
				   KAPPA
				   FIX_OMEGA
				   OMEGA
				   FIX_ALPHA
				   ALPHA
				   MALPHA
				   NCATG
				   CLOCK
				   GETSE
				   RATEANCESTOR
				   SMALL_DIFF)],@args);

  throw ("Have sequence input and input file - which do I use?") 
    if (defined $aligned_seqs && defined $seqfile);
  
  $self->_aligned_seqs($aligned_seqs) if $aligned_seqs;

  $self->_codeml_executable($executable) if $executable;
  $self->work_dir($work_dir)             if $work_dir;
  $self->_seqfile($seqfile)              if $seqfile;
  $self->_treefile($treefile)            if $treefile;

  $self->_pre_copying; # Set up our work dir the way we like it.

  $self->_outfile($outfile)              if $outfile;
  $self->_noisy($noisy)                  if $noisy;
  $self->_verbose($verbose)              if $verbose;
  $self->_runmode('set', $runmode);
  $self->_seqtype($seqtype)              if $seqtype;
  $self->_codonfreq($codonfreq)          if $codonfreq;
  $self->_aadist($aadist)                if $aadist;
  $self->_aaratefile($aaratefile)        if $aaratefile;
  $self->_model('set', $model);
  $self->_nssites($nssites);
  $self->_icode($icode);
  $self->_treefile($mgene)               if $mgene;
  $self->_treefile($fix_kappa)           if $fix_kappa;
  $self->_treefile($kappa)               if $kappa;
  $self->_treefile($fix_omega)           if $fix_omega;
  $self->_treefile($omega)               if $omega;
  $self->_treefile($fix_alpha)           if $fix_alpha;
  $self->_treefile($alpha)               if $alpha;
  $self->_treefile($malpha)              if $malpha;
  $self->_treefile($ncatg)               if $ncatg;
  $self->_treefile($clock)               if $clock;
  $self->_treefile($getse)               if $getse;
  $self->_treefile($rateancestor)        if $rateancestor;
  $self->_treefile($small_diff)          if $small_diff;
  
  return $self
}

sub DESTROY {
  my $self = shift;

  unlink $self->_config_file if $self->_config_file; 
  unlink $self->_seqfile if $self->_seqfile;

  unlink ($self->work_dir . "/rub",    $self->work_dir . "/rst",
	  $self->work_dir . "/rst1",   $self->work_dir . "/lnf",
	  $self->work_dir . "/2ML.dN", $self->work_dir . "/2ML.dS",
	  $self->work_dir . "/2ML.t",  $self->work_dir . "/2NG.dN", 
	  $self->work_dir . "/2NG.dS", $self->work_dir . "/2NG.t") 
    if $self->_outfile;

  unlink $self->_outfile if $self->_outfile;

  warning ("Paml execution resulted in a core file being dumped.")
    if (unlink $self->work_dir . "/core");

  warning ("Could not remove working directory [". $self->work_dir ."].\n$!") 
    if (! rmdir $self->work_dir);
}


sub run_codeml {
  my $self = shift;

  # Directory changing shenanigans to make PAML behave!

  my $user_dir = cwd();

  # Move to our working dir
  warning ($!) if (! chdir $self->work_dir);

  $self->_write_seqs($self->_seqfile);
  $self->_write_config_file;

  my $command = $self->_codeml_executable . " " . $self->_config_file;
  print STDERR $command . "\n";

  eval {
    system($command)
  };

  if ($@ or -e $self->work_dir . "/core") {
    throw ("Something went wrong when codeml was executed.\n" . $@);
  }

  # Move back to original dir
  print STDERR $! if (! chdir $user_dir);

  my $parser = Bio::Tools::Phylo::PAML->new('-file' => $self->_outfile, 
					    '-dir'  => $self->work_dir);

  return $parser;
}


# Make a temporary playground for PAML and some of its anti-social ways.
# Execute application from here such that all unwanted temp files can be 
# neatly deleted.  Also, can copy assorted external files, like grantham.dat
# here before execution.  Give temp dir a unique name to avoid conflicts
# from multiple simultaneous jobs.

sub _pre_copying {
  my $self = shift;

  #system("cp " . $self->_executable . " " . $self->work_dir);

  # Ultimately all the datafiles need to be copied into here.

  return 1;
}


sub _write_config_file {
  my ($self) = @_;

  my $config_string;

  $config_string .= '     seqfile = ' . $self->_seqfile      . "\n";
  $config_string .= '    treefile = ' . $self->_treefile     . "\n";
  $config_string .= '     outfile = ' . $self->_outfile      . "\n";
  $config_string .= '       noisy = ' . $self->_noisy        . "\n";
  $config_string .= '     verbose = ' . $self->_verbose      . "\n";
  $config_string .= '     runmode = ' . $self->_runmode      . "\n";
  $config_string .= '     seqtype = ' . $self->_seqtype      . "\n";
  $config_string .= '   CodonFreq = ' . $self->_codonfreq    . "\n";
  $config_string .= '      aaDist = ' . $self->_aadist       . "\n";
  $config_string .= '  aaRatefile = ' . $self->_aaratefile   . "\n";
  $config_string .= '       model = ' . $self->_model        . "\n";
  $config_string .= '     NSsites = ' . $self->_nssites      . "\n";
  $config_string .= '       icode = ' . $self->_icode        . "\n";
  $config_string .= '       Mgene = ' . $self->_mgene        . "\n";
  $config_string .= '   fix_kappa = ' . $self->_fix_kappa    . "\n";
  $config_string .= '       kappa = ' . $self->_kappa        . "\n";
  $config_string .= '   fix_omega = ' . $self->_fix_omega    . "\n";
  $config_string .= '       omega = ' . $self->_omega        . "\n";
  $config_string .= '   fix_alpha = ' . $self->_fix_alpha    . "\n";
  $config_string .= '       alpha = ' . $self->_alpha        . "\n";
  $config_string .= '      Malpha = ' . $self->_malpha       . "\n";
  $config_string .= '       ncatG = ' . $self->_ncatg        . "\n";
  $config_string .= '       clock = ' . $self->_clock        . "\n";
  $config_string .= '       getSE = ' . $self->_getse        . "\n";
  $config_string .= 'RateAncestor = ' . $self->_rateancestor . "\n";
  $config_string .= '  Small_Diff = ' . $self->_small_diff   . "\n";

  my $file_with_path = $self->work_dir . "paml_codeml.ctl";

  open(CONFIG, ">$file_with_path") or die "Cant open config file for writing.";

  print CONFIG $config_string;

  close(CONFIG);

  $self->_config_file($file_with_path);

  return 1
}


### Files and executables ###


sub _codeml_executable {
  my $self = shift;

  if (@_) {
    $self->{_codeml_executable} = shift;
  }

  unless ($self->{_codeml_executable}) {
    print STDERR "Explicit location of codeml executable not set.  Assuming\n" .
      "that codeml is in your path.";

    $self->{_codeml_executable} = 'codeml';
  }

  return $self->{_codeml_executable}
}


sub work_dir {
  my ($self) = shift;

  unless ($self->{_work_dir}) {

    my $work_dir;

    if (@_) {
      my $arg_dir = shift;
      $work_dir = $arg_dir . '/' . 'paml_temp_' . time . '/';
    } else {
      $work_dir = '/tmp/paml_temp_' . time . '/';
    }

    unless (-d $work_dir){
      mkdir $work_dir;
    }

    $self->{_work_dir} = $work_dir;
  }

  return $self->{_work_dir}
}


sub _write_seqs {
  my ($self, $filename) =  @_;

  throw ("No sequences to write.")
    unless $self->_aligned_seqs;

  system("rm -f $filename");

  open(OUT, ">$filename") or die "Cant write to file [$filename]\n";

  print OUT scalar @{$self->_aligned_seqs} . " " . $self->_aligned_seqs->[0]->length . "\n";

  foreach my $aligned_seq (@{$self->_aligned_seqs}){
    print OUT $aligned_seq->display_id . "\n" . $aligned_seq->seq . "\n";
  }

  close(OUT);

  return 1
}


sub _aligned_seqs {
  my $self = shift;

  if (@_){
    $self->{_aligned_seqs} = shift;

    throw ("No sequences to set.")
      unless $self->{_aligned_seqs} and 
	$self->{_aligned_seqs}->[0]->isa("Bio::Seq");
  }

  return $self->{_aligned_seqs}
}


sub _config_file {
  my ($self) = shift;


  if (@_){
    $self->{_config_file} = shift;
  }

#print "Config file for PAML : " . $self->{_config_file} . "\n";
  return $self->{_config_file}
}


### Config options ###

sub _seqfile {
  my ($self) = shift;

  if (@_) {
    $self->{_seqfile} = shift;
  }

  unless ($self->{_seqfile}) {
    $self->{_seqfile} = $self->work_dir . 'paml_in_' . time;
  }

  return $self->{_seqfile}
}



sub _treefile {
  my ($self) = shift;

  if (@_) {
    $self->{_treefile} = shift;
  }

  if (!$self->{_treefile} && $self->_runmode == 0) {
    throw ("Input tree file has not been defined.");
  }

  if (!$self->{_treefile}) {
    $self->{_treefile} = '';
  }

  return $self->{_treefile}
}



sub _outfile {
  my ($self) = shift;

  if (@_) {
    $self->{_outfile} = shift;
  }

  unless ($self->{_outfile}) {
    $self->{_outfile} = $self->work_dir . 'paml_out_' . time;
  }

  return $self->{_outfile}
}


sub _noisy {
  my ($self) = shift;

  if (@_) {
    $self->{_noisy} = shift;
  }

  if (!$self->{_noisy}){
    return 0
  } else {
    return $self->{_noisy}
  }
}

sub _verbose {
  my ($self) = shift;

  if (@_) {
    $self->{_verbose} = shift;
  }

  if (!$self->{_verbose}){
    return 0
  } else {
    return $self->{_verbose}
  }
}

sub _runmode {
  my ($self, $set, $value) = @_;

  if (defined $set) {
    $self->{_runmode} = $value + 1;  # Add one to cope with 0 as a valid option.
  }

  unless (defined $self->{_runmode}){
    throw ("Runmode has not been set.");
  } 

  return $self->{_runmode} - 1 # Remember to remove one to return correct option.
}



sub _seqtype {
  my ($self) = shift;

  if (@_) {
    $self->{_seqtype} = shift;
  }

  unless (defined $self->{_seqtype}){
    throw ("Seqtype has not been set.");
  } 

  return $self->{_seqtype}
}


sub _codonfreq {
  my ($self) = shift;

  if (@_) {
    $self->{_codonfreq} = shift;
  }

  if ($self->_seqtype == 2 && !defined $self->{_codonfreq}){
#    warning ("Codonfreq not set.  Defaulting to option 2: F3X4 frequencies.");
    return 2;
  } 

  if (! $self->{_codonfreq}){
    return 2;
  }

  return $self->{_codonfreq}
}


sub _aadist {
  my ($self) = shift;

  if (@_) {
    $self->{_aadist} = shift;
  }

  unless (defined $self->{_aadist}){
#    warning ("Aadist not set.  Defaulting to option 1: equal distances.");
    return 0;
  } 

  return $self->{_aadist}
}


sub _aaratefile {
  my ($self) = shift;

  if (@_) {
    $self->{_aaratefile} = shift;
  }

  if (!defined $self->{_aaratefile}){
    if ($self->_seqtype == 2 && $self->_model > 1) {
#      warning ("Aaratefile not set.  Defaulting to file wag.dat.");
    }
    return 'wag.dat';
  } 

  return $self->{_aaratefile}
}


sub _model {
  my ($self, $set, $value) = @_;

  if ($set) {
    $self->{_model} = $value + 1; # Add one to cope with 0 as a valid option.
  }

  unless (defined $self->{_model}){
    throw ("Model has not been set.");
  } 

  return $self->{_model} - 1 # Remove one to obtain correct value.
}

sub _nssites {
  my ($self) = shift;

  if (@_) {
    $self->{_nssites} = shift;
  }

  unless (defined $self->{_nssites}){
#    warning ("NSsites has not been set.  Defaulting to 0.");
  } 

  return $self->{_nssites}
}


sub _icode {
  my ($self) = shift;

  if (@_) {
    $self->{_icode} = shift;
  }

  unless (defined $self->{_icode}){
    warning ("PAML: Defaulting to Universal Genetic Code.\n");
    return 0 # Universal code
  } 

  return $self->{_icode}
}


sub _mgene {
  my ($self) = shift;

  if (@_) {
    $self->{_mgene} = shift;
  }

  unless (defined $self->{_mgene}){
#    print STDERR "Defaulting to mgene = 0\n";
    return 0 # default to separate
  } 

  return $self->{_mgene}
}


sub _fix_kappa {
  my ($self) = shift;

  if (@_) {
    $self->{_fix_kappa} = shift;
  }

  unless (defined $self->{_fix_kappa}){
#    print STDERR "Defaulting to unfixed kappa - will estimate kappa.\n";
    return 0 # default to unfixed kappa
  } 

  return $self->{_fix_kappa}
}


sub _kappa {
  my ($self) = shift;

  if (@_) {
    $self->{_kappa} = shift;
  }

  if ($self->_fix_kappa != 0 && !defined $self->{_kappa}){
    throw ("Youve asked for kappa to be fixed, but\n" .
		 "have not specified a value for it to be\n" .
		 "fixed to.  Try setting '-kappa' => x or\n" . 
		 "whatever");
  } 

  if (! $self->{_kappa}){
    return 2;
  }

  return $self->{_kappa}
}

sub _fix_omega {
  my ($self) = shift;

  if (@_) {
    $self->{_fix_omega} = shift;
  }

  unless (defined $self->{_fix_omega}){
#    print STDERR "Defaulting to unfixed omega - will estimate omega.\n";
    return 0 # default to unfixed omega
  } 

  return $self->{_fix_omega}
}


sub _omega {
  my ($self) = shift;

  if (@_) {
    $self->{_omega} = shift;
  }

  if ($self->_fix_omega != 0 && !defined $self->{_omega}){
    throw ("Youve asked for omega to be fixed, but\n" .
		 "have not specified a value for it to be\n" .
		 "fixed to.  Try setting '-omega' => x or\n" . 
		 "whatever");
  } 

  if (! defined $self->{_omega}){
    return 0;
  }

  return $self->{_omega}
}


sub _fix_alpha {
  my ($self) = shift;

  if (@_) {
    $self->{_fix_alpha} = shift;
  }

  unless (defined $self->{_fix_alpha}){
#    print STDERR "Defaulting to unfixed alpha - will estimate alpha.\n";
    return 0 # default to unfixed alpha
  } 

  return $self->{_fix_alpha}
}


sub _alpha {
  my ($self) = shift;

  if (@_) {
    $self->{_alpha} = shift;
  }

  if ($self->_fix_alpha != 0 && !defined $self->{_alpha}){
    throw ("Youve asked for alpha to be fixed, but\n" .
		 "have not specified a value for it to be\n" .
		 "fixed to.  Try setting '-alpha' => x or\n" . 
		 "whatever");
  } 

  if (! defined $self->{_alpha}){
    return 0;
  }

  return $self->{_alpha}
}

sub _malpha {
  my ($self) = shift;

  if (@_) {
    $self->{_malpha} = shift;
  }

  if (!defined $self->{_malpha}){
    # silently default to 0.  Ive never used this - what 
    # other parameters does it depend on?
    return 0 
  } 

  return $self->{_malpha}
}

sub _ncatg {
  my ($self) = shift;

  if (@_) {
    $self->{_ncatg} = shift;
  }

  if (!defined $self->{_ncatg}){
    # silently default to 10.  Ive never used this - what 
    # other parameters does it depend on?
    return 10 
  } 

  return $self->{_ncatg}
}


sub _clock {
  my ($self) = shift;

  if (@_) {
    $self->{_clock} = shift;
  }

  if (!defined $self->{_clock}){
    # What other parameters are set with clock - should check these before throwing a warning.
#    warning ("Clock parameter not set.  Defaulting to option 0: no clock.");
    return 0
  } 

  return $self->{_clock}

}


sub _getse {
  my ($self) = shift;

  if (@_) {
    $self->{_getse} = shift;
  }

  if (!defined $self->{_getse}){
#    print STDERR "By default, not estimating SE.\n";
    return 0
  } 

  return $self->{_getse}
}


sub _rateancestor {
  my ($self) = shift;

  if (@_) {
    $self->{_rateancestor} = shift;
  }

  if (!defined $self->{_rateancestor}){
#    print STDERR "By default, RateAncestor option set to 0 - not" .
#      "\nestimating ancestral states.\n";
    return 0
  } 

  return $self->{_rateancestor}
}


# What is this value for?
sub _small_diff {
  my ($self) = shift;

  if (@_) {
    $self->{_small_diff} = shift;
  }

  if (!defined $self->{_small_diff}){
    # Dont know what this does - default set directly from sample file.
    return '0.5e-6'
  } 

  return $self->{_small_diff}

}
