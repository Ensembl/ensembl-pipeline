#
#
# Cared for by Thomas Down <td2@sanger.ac.uk>
#
# Based on CPG.pm by Val Curwen
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::Runnable::EponineTSS

=head1 SYNOPSIS

#create and fill Bio::Seq object

my $seqfile = '/nfs/disk65/mq2/temp/bA151E14.seq';

my $seq = Bio::Seq->new();

my $seqstream = Bio::SeqIO->new(-file => $seqfile, -fmt => 'Fasta');

$seq = $seqstream->next_seq();

#create Bio::EnsEMBL::Pipeline::Runnable::EponineTSS object

my $eponine = Bio::EnsEMBL::Pipeline::Runnable::EponineTSS->new (-QUERY => $seq);

$eponine->workdir($workdir);

$eponine->run();

my @results = $eponine->output();

=head1 DESCRIPTION

EponineTSS takes a Bio::Seq (or Bio::PrimarySeq) object and runs the
standalone JAX version of the Eponine Transcription Start Site
predictor.

=head1 CONTACT

Thomas Down: td2@sanger.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Pipeline::Runnable::EponineTSS;

use vars qw(@ISA);
use strict;
# Object preamble - inherits from Bio::Root::RootI;

use Bio::EnsEMBL::Pipeline::RunnableI;
use Bio::EnsEMBL::SeqFeature;
use Bio::EnsEMBL::Analysis;
use Bio::Seq;
use Bio::SeqIO;
use Bio::Root::RootI;

BEGIN {
    require "Bio/EnsEMBL/Pipeline/pipeConf.pl";
}

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI);

=head2 new

    Title   :   new
    Usage   :   my obj =  Bio::EnsEMBL::Pipeline::Runnable::EponineTSS->new (-QUERY => $seq);
    Function:   Initialises an EponineTSS object
    Returns :   an EponineTSS Object
    Args    :   A Bio::Seq object (-QUERY), any arguments (-LENGTH, -GC, -OE) 

=cut

sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);    
  
  $self->{'_flist'}     = [];    # an array of Bio::SeqFeatures
  $self->{'_sequence'}  = undef; # location of Bio::Seq object
  $self->{'_java'}      = undef; # location of java vm
  $self->{'_epojar'}    = undef; # location of eponine-scan.jar executable JAR file.
  $self->{'_workdir'}   = undef; # location of temp directory
  $self->{'_filename'}  = undef; # file to store Bio::Seq object
  $self->{'_results'}   = undef; # file to store results of eponine-scan
  $self->{'_protected'} = [];    # a list of files protected from deletion ???
  $self->{'_threshold'} = 0.999; # minimum posterior for filtering predictions 
  
  print STDERR "args: ", @args, "\n";

  my( $sequence, $java, $epojar, $threshold) = $self->_rearrange([qw(QUERY 
								     JAVA
								     EPOJAR
								     THRESHOLD)],
								 @args);
  
  $self->query($sequence) if ($sequence);       

  my $bindir = $::pipeConf{'bindir'} || undef;

  if (-x $java) {
    # passed from RunnableDB (full path assumed)
    $self->java($java);
  }
  elsif ($::pipeConf{'bin_JAVA'} && -x ($java = "$::pipeConf{'bin_JAVA'}")) {
    $self->java($java);
  }
  elsif ($bindir && -x ($java = "/usr/opt/java130/bin/java")) {
    $self->java($java);
    }
  else {
    # search shell $PATH
    eval {
      $self->java($self->locate_executable('java'));
    };
    if ($@) {
      $self->throw("Can't find executable java");
    }
  }
  
  if (-e $epojar) {
      $self->epojar($epojar);
  } elsif ($::pipeConf{'bin_EPONINE'} && -e ($epojar = $::pipeConf{'bin_EPONINE'})) {
      $self->epojar($epojar);
  } else {
      $self->throw("Can't find eponine-scan.jar");
  }
  
  if (defined $threshold && $threshold >=0 ){
      $self->threshold($threshold);
  } else {
      $self->threshold(50); 
  }
  
  return $self; # success - we hope!
}

#################
# get/set methods 
#################
# really ough to be renamed "sequence" but this involves rewriting RunnableI::writefile and also any other modules that inherit from it.
# to do!
sub query {
  my ($self, $seq) = @_;
  if ($seq)
  {
      unless ($seq->isa("Bio::PrimarySeqI") || $seq->isa("Bio::SeqI")) 
      {
	  $self->throw("Input isn't a Bio::SeqI or Bio::PrimarySeqI");
      }
      $self->{'_sequence'} = $seq ;

      $self->filename($self->query->id.".$$.seq");
      $self->results($self->filename.".out");
  }
  return $self->{'_sequence'};
}

=head2 java

    Title   :   java
    Usage   :   $obj->java('/usr/opt/java130/bin/java');
    Function:   Get/set method for the location of java VM
    Args    :   File path (optional)

=cut

sub java {
    my ($self, $location) = @_;
    if ($location)
    {
      $self->throw("java not found at $location: $!\n") 
	unless (-e $location);
        $self->{'_java'} = $location ;
      $self->throw("java not executable: $!\n") 
      	unless (-x $location);
    }
    return $self->{'_java'};
}

=head2 epojar

    Title   :   epojar
    Usage   :   $obj->epojar('/some/path/to/eponine-scan.jar');
    Function:   Get/set method for the location of the eponine-scan executable JAR
    Args    :   Path (optional)

=cut

sub epojar {
    my ($self, $location) = @_;
    if ($location)
    {
      $self->throw("eponine-scan.jar not found at $location: $!\n") 
	  unless (-e $location);
      $self->{'_epojar'} = $location ;
    }
    return $self->{'_epojar'};
}

=head2 threshold

    Title   :   threshold
    Usage   :   $obj->threshold('50');
    Function:   Get/set method for posterior required before a TSS is reported
    Args    :   threshold (optional)

=cut

sub threshold {
    my ($self, $threshold) = @_;
    if (defined $threshold) {
      $self->{'_threshold'} = $threshold ;
    }
    return $self->{'_threshold'};
}

###########
# Analysis methods
##########

=head2 run

    Title   :  run
    Usage   :   $obj->run($workdir, $args)
    Function:   Runs eponine-scan script and creates array of features
    Returns :   none
    Args    :   optional $workdir and $args (e.g. '-ace' for ace file output)

=cut

sub run {
    my ($self, $dir, $args) = @_;
    #set arguments for epo
    #check query
    my $seq = $self->query() || $self->throw("Seq required for EponineTSS\n");
    #set directory if provided
    $self->workdir('/tmp') unless ($self->workdir($dir));
    $self->checkdir();

    # reset filename and results as necessary
    $self->filename($self->workdir()."/".$self->filename());
    $self->results($self->workdir()."/".$self->results());

    #write sequence to file
    $self->writefile();        
    $self->run_eponine();
    #parse output of eponine
    $self->parse_results();
    $self->deletefiles();
}


=head2 run_eponine

    Title   :  run_eponine
    Usage   :   $obj->run_eponine($filename)
    Function:   execs the Java VM to run eponine
    Returns :   none
    Args    :   none

=cut
sub run_eponine {
    my ($self) = @_;
    #run eponine
    print "Running eponine-scan\n";
    $self->throw("Error running eponine-scan on ".$self->filename."\n") 
        if (system ($self->java.' -fast -jar '.$self->epojar.' -seq '.$self->filename.' -threshold '.$self->threshold." > ".$self->results)); 
}

=head2 parse_results

    Title   :  parse_results
    Usage   :   $obj->parse_results($filename)
    Function:   Parses eponine output to give a set of features
                parsefile can accept filenames, filehandles or pipes (\*STDIN)
    Returns :   none
    Args    :   optional filename

=cut
sub parse_results {
    my ($self, $filehandle) = @_;
#   my $resfile = $self->workdir()."/".$self->results();
   my $resfile = $self->results();

    if (-e $resfile) {
      if (-z $self->results) {  
	print STDERR "eponine didn't find any TSSs\n";
	return; }       
      else {
        open (EPONINEOUT, "<$resfile") or $self->throw("Error opening ", $resfile, " \n");#
        $filehandle = \*EPONINEOUT;
      }
    }
    else { #it'a a filehandle
      $filehandle = $resfile;
    }
    
    #extract values
    while (<$filehandle>)
    {  
        if (! /^\#/) #ignore introductory lines
        {
            my @element = split;
            my (%feature);
            $feature {name} = $element[0];
            $feature {score} = $element[5];
            $feature {start} = $element[3];
            $feature {end} = $element[4];
            $feature {strand} = $element[6];
            $feature {source}= 'Eponine';
            $feature {primary}= 'TSS';
	    $feature {program} = 'eponine-scan';
	    $feature {program_version} = '2';
          
	    $self->create_feature(\%feature);

	  
        }
    }
    close $filehandle;   
}


##############
# input/output methods
#############

=head2 output

    Title   :   output
    Usage   :   obj->output()
    Function:   Returns an array of features
    Returns :   Returns an array of features
    Args    :   none

=cut

sub output {
    my ($self) = @_;
    return @{$self->{'_flist'}};
}

=head2 create_feature

    Title   :   create_feature
    Usage   :   obj->create_feature($feature)
    Function:   Returns an array of features
    Returns :   Returns an array of features
    Args    :   none

=cut
sub create_feature {
    my ($self, $feat) = @_;

    #create analysis object
    my $analysis_obj = Bio::EnsEMBL::Analysis->new
                        (   -db              => undef,
                            -db_version      => undef,
                            -program         => $feat->{'program'},
                            -program_version => $feat->{'program_version'},
                            -gff_source      => $feat->{'source'},
                            -gff_feature     => $feat->{'primary'});

    #create and fill Bio::EnsEMBL::Seqfeature object
    my $tss = Bio::EnsEMBL::SeqFeature->new
                        (   -seqname => $feat->{'name'},
                            -start   => $feat->{'start'},
                            -end     => $feat->{'end'},
                            -strand  => $feat->{'strand'},
                            -score   => $feat->{'score'},
                            -source_tag  => $feat->{'source'},
                            -primary_tag => $feat->{'primary'},
                            -analysis => $analysis_obj);  

    if ($tss)
      {
	$tss->validate();

	# add to _flist
	push(@{$self->{'_flist'}}, $tss);
      }
}
