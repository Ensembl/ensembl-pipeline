#
#
# Cared for by Val Curwen  <vac@sanger.ac.uk>
#
# Copyright Val Curwen
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::Runnable::CPG

=head1 SYNOPSIS

#create and fill Bio::Seq object

my $clonefile = '/nfs/disk65/mq2/temp/bA151E14.seq';

my $seq = Bio::Seq->new();

my $seqstream = Bio::SeqIO->new(-file => $clonefile, -fmt => 'Fasta');

$seq = $seqstream->next_seq();

#create Bio::EnsEMBL::Pipeline::Runnable::CPG object

my $cpg = Bio::EnsEMBL::Pipeline::Runnable::CPG->new (-QUERY => $seq);

$cpg->workdir($workdir);

$cpg->run();

my @results = $cpg->output();

=head1 DESCRIPTION

CPG takes a Bio::Seq (or Bio::PrimarySeq) object and runs Gos Micklems cpg 
on it. The resulting output file is parsed to produce a set of features.

=head1 CONTACT

Val Curwen: vac@sanger.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Pipeline::Runnable::CPG;

use vars qw(@ISA);
use strict;
# Object preamble - inherits from Bio::Root::RootI;

use Bio::EnsEMBL::Pipeline::RunnableI;
use Bio::EnsEMBL::SimpleFeature;
use Bio::EnsEMBL::Analysis;
use Bio::Seq;
use Bio::SeqIO;
use Bio::Root::RootI;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI);

=head2 new

    Title   :   new
    Usage   :   my obj =  Bio::EnsEMBL::Pipeline::Runnable::CPG->new (-QUERY => $seq);
    Function:   Initialises CPG object
    Returns :   a CPG Object
    Args    :   A Bio::Seq object (-QUERY), any arguments (-LENGTH, -GC, -OE) 

=cut

sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);    
  
  $self->{'_flist'}     = [];    # an array of Bio::EnsEMBL::SimpleFeatures
  $self->{'_sequence'}  = undef; # location of Bio::Seq object
  $self->{'_cpg'}       = undef; # location of cpg executable
  $self->{'_workdir'}   = undef; # location of temp directory
  $self->{'_filename'}  = undef; # file to store Bio::Seq object
  $self->{'_results'}   = undef; # file to store results of cpg
  $self->{'_protected'} = [];    # a list of files protected from deletion ???
  $self->{'_min_length'}= 0;     # minimum length required for a
                                 # predicted cpg island to be reported
  $self->{'_min_gc'}    = 0;     # minimum %GC required for a predicted 
                                 # cpg island to be reported
  $self->{'_min_oe'}    = 0;     # minimum observed/expected ratio 
                                 # required for a predicted cpg island to
                                 # be reported
  
  print STDERR "args: ", @args, "\n";

  my( $sequence, $len, $gc, $oe, $cpg) = $self->_rearrange([qw(QUERY 
							       LENGTH 
							       GC 
							       OE 
							       CPG)],
							   @args);


  $cpg = 'cpg' unless ($cpg);

  $self->query($sequence) if ($sequence);       
  $self->cpg($self->find_executable($cpg));
  
  if (defined $len && $len >=0 ) { 
    $self->min_length($len); }
  else { 
    $self->min_length(1000); } # for parsing output
  
  if (defined $gc && $gc >=0 ){
    $self->min_gc($gc);}
  else {
    $self->min_gc(50); }
  
  if (defined $oe && $oe >=0 ){
    $self->min_oe($oe); }
  else {
    $self->min_oe(0.6); }
  
  return $self; # success - we hope!
}

#################
# get/set methods 
#################

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

=head2 cpg

    Title   :   cpg
    Usage   :   $obj->cpg('/usr/local/pubseq/bin/cpg');
    Function:   Get/set method for the location of cpg script
    Args    :   File path (optional)

=cut

sub cpg {
    my ($self, $location) = @_;
    if ($location)
    {
      $self->throw("cpg not found at $location: $!\n") 
	unless (-e $location);
        $self->{'_cpg'} = $location ;
      $self->throw("cpg not executable: $!\n") 
      	unless (-x $location);
    }
    return $self->{'_cpg'};
}

=head2 min_length

    Title   :   min_length
    Usage   :   $obj->min_length(1000);
    Function:   Get/set method for minimum length required before a cpg island is reported
    Args    :   length(optional)

=cut

sub min_length {
    my ($self, $len) = @_;
    if (defined $len) {
      $self->{'_min_length'} = $len ;
    }
    return $self->{'_min_length'};
}

=head2 min_gc

    Title   :   min_gc
    Usage   :   $obj->min_gc('50');
    Function:   Get/set method for minimum %GC required before a cpg island is reported
    Args    :   gc(optional)

=cut

sub min_gc {
    my ($self, $gc) = @_;
    if (defined $gc) {
      $self->{'_min_gc'} = $gc ;
    }
    return $self->{'_min_gc'};
}

=head2 min_oe

    Title   :   min_oe
    Usage   :   $obj->min_oe(0.6);
    Function:   Get/set method for minimum observed/expected ratio required before a cpg island is reported
    Args    :   oe(optional)

=cut

sub min_oe {
    my ($self, $oe) = @_;
    if (defined $oe) {
      $self->{'_min_oe'} = $oe ; 
    }
    return $self->{'_min_oe'};
}

###########
# Analysis methods
##########

=head2 run

    Title   :  run
    Usage   :   $obj->run($workdir, $args)
    Function:   Runs cpg script and creates array of features
    Returns :   none
    Args    :   optional $workdir and $args (e.g. '-ace' for ace file output)

=cut

sub run {
    my ($self) = @_;
    #set arguments for cpg
    #check clone
    my $seq = $self->query() || $self->throw("Clone required for cpg\n");
    #set directory if provided
    $self->workdir('/tmp') unless $self->workdir();
    $self->checkdir();

    # reset filename and results as necessary
    $self->filename($self->workdir()."/".$self->filename());
    $self->results($self->workdir()."/".$self->results());

    #write sequence to file
    $self->writefile();        
    $self->run_cpg();
    #parse output of cpg
    $self->parse_results();
    $self->deletefiles();
}


=head2 run_cpg

    Title   :  run_cpg
    Usage   :   $obj->run_cpg($filename)
    Function:   makes the system call to cpg
    Returns :   none
    Args    :   optional filename

=cut
sub run_cpg {
    my ($self) = @_;
    #run cpg
    print "Running cpg\n";
    $self->throw("Error running cpg on ".$self->filename."\n") 
        if (system ($self->cpg.' '.$self->filename." > ".$self->results)); 
}

=head2 parse_results

    Title   :  parse_results
    Usage   :   $obj->parse_results($filename)
    Function:   Parses cpg output to give a set of features
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
	print STDERR "cpg didn't find any cpg islands\n";
	return; }       
      else {
        open (CPGOUT, "<$resfile") or $self->throw("Error opening ", $resfile, " \n");#
        $filehandle = \*CPGOUT;
      }
    }
    else { #it'a a filehandle
      $filehandle = $resfile;
    }
    
    #extract values
    while (<$filehandle>)
    {  
        if (/\d+/) #ignore introductory lines
        {
            my @element = split;
	    my $oe = $element[7];
	    if($oe eq "-") { $oe = 0; }
	    next unless (($element[2] - $element[1] + 1 >= $self->min_length)
			 && ($element[6] >= $self->min_gc) 
			 && $oe >= $self->min_oe);
	  
            my (%feature);
            $feature {name} = $element[0];
            $feature {score} = $element[3];
            $feature {start} = $element[1];
            $feature {end} = $element[2];
            $feature {strand} = 0;
            $feature {program} = 'cpg';
	    $feature {program_version} = '1';
	    $feature {display_label} = "oe = $oe";
          
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
    my $cpg = Bio::EnsEMBL::SimpleFeature->new
                        (   -seqname => $feat->{'name'},
                            -start   => $feat->{'start'},
                            -end     => $feat->{'end'},
                            -strand  => $feat->{'strand'},
                            -score   => $feat->{'score'},
                            -analysis => $analysis_obj);  

    $cpg->display_label($feat->{'display_label'});
    if ($cpg)
      {
	$cpg->validate();

	# add to _flist
	push(@{$self->{'_flist'}}, $cpg);
      }
}
