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

Bio::EnsEMBL::Pipeline::Runnable::tRNAscan_SE

=head1 SYNOPSIS

#create and fill Bio::Seq object
my $seqfile = '/nfs/disk65/mq2/temp/bA151E14.seq';
my $seq = Bio::Seq->new();
my $seqstream = Bio::SeqIO->new(-file => $seqfile, -fmt => 'Fasta');
$seq = $seqstream->next_seq();
#create Bio::EnsEMBL::Pipeline::Runnable::tRNAscan_SE object
my $trna = Bio::EnsEMBL::Pipeline::Runnable::tRNAscan_SE->new (-QUERY => $seq);
$trna->workdir($workdir);
$trna->run();
my @results = $trna->output();

=head1 DESCRIPTION

tRNAscan_SE takes a Bio::Seq (or Bio::PrimarySeq) object and runs tRNAscan-SE on it. 
The resulting output file is parsed to produce a set of features.

=head1 CONTACT

Val Curwen: vac@sanger.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Pipeline::Runnable::tRNAscan_SE;

use vars qw(@ISA);
use strict;
# Object preamble - inherits from Bio::Root::RootI;

use Bio::EnsEMBL::Pipeline::RunnableI;
use Bio::EnsEMBL::SeqFeature;
use Bio::EnsEMBL::Analysis;
use Bio::Seq;
use Bio::SeqIO;
use Bio::Root::RootI;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI Bio::Root::RootI );

=head2 new

    Title   :   new
    Usage   :   my obj =  Bio::EnsEMBL::Pipeline::Runnable::tRNAscan_SE->new (-QUERY => $seq);
    Function:   Initialises tRNAscan_SE object
    Returns :   a tRNAscan_SE Object
    Args    :   A Bio::Seq object (-QUERY) 

=cut

sub new {
  my ($class,@args) = @_;

  my $self = $class->SUPER::new(@_);    
  
  $self->{'_flist'} = [];           #an array of Bio::SeqFeatures
  $self->{'_sequence'}  = undef;    #location of Bio::Seq object
  $self->{'_tRNAscan_SE'} = undef;  #location of tRNAscan_SE executable
  $self->{'_workdir'}   = undef;    #location of temp directory
  $self->{'_filename'}  = undef;    #file to store Bio::Seq object
  $self->{'_results'}   = undef;    #file to store results of tRNAscan_SE
  $self->{'_protected'} = [];       #a list of files protected from deletion
  
  my( $sequence, $tRNAscan_SE) = $self->_rearrange([qw(QUERY 
						       TRNASCAN_SE)], 
						   @args);


  $tRNAscan_SE = 'tRNAscan-SE' unless $tRNAscan_SE;

  $self->query($sequence) if ($sequence);       
  $self->tRNAscan_SE($self->find_executable($tRNAscan_SE));

  return $self;
}

#################
# get/set methods 
#################
# really ought to be renamed "sequence" but this involves rewriting RunnableI::writefile and also any other modules that inherit from it.
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

	$self->queryname($self->query->id);
	$self->filename($self->query->id.".$$.seq");
	$self->results($self->filename.".out");
    }
    return $self->{'_sequence'};
}

=head2 tRNAscan_SE

    Title   :   tRNAscan_SE
    Usage   :   $obj->tRNAscan_SE('/usr/local/pubseq/bin/tRNAscan-SE');
    Function:   Get/set method for the location of tRNAscan-SE script
    Args    :   File path (optional)

=cut

sub tRNAscan_SE {
    my ($self, $location) = @_;
    if ($location)
    {
        $self->throw("tRNAscan-SE not found at $location: $!\n") 
	    unless (-e $location);
        $self->{'_tRNAscan_SE'} = $location ;
    }
    return $self->{'_tRNAscan_SE'};
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
    my ($self, $dir, $args) = @_;
    #set arguments for tRNAscan_SE
    #check seq
    my $seq = $self->query() || $self->throw("Seq required for tRNAscan_SE\n");
    #set directory if provided
    $self->workdir('/tmp') unless ($self->workdir($dir));
    $self->checkdir();

    # reset filename and results as necessary
    my $tmp = $self->workdir();
    my $input = $tmp."/".$self->filename();
    $self->filename($input);
    $tmp .= "/".$self->results();
    $self->results($tmp);

    #write sequence to file
    $self->writefile();        
    $self->run_tRNAscan_SE();
    #parse output of tRNAscan_SE
    $self->parse_results();
    $self->deletefiles();
}

=head2 run_tRNAscan_SE

    Title   :  run_tRNAscan_SE
    Usage   :   $obj->run_tRNAscan_SE($filename)
    Function:   makes the system call to tRNAscan_SE
    Returns :   none
    Args    :   optional filename

=cut

sub run_tRNAscan_SE {
    my ($self) = @_;
    #run tRNAscan_SE
    print "Running tRNAscan_SE\n";
    $self->throw("Error running tRNAscan_SE on ".$self->filename."\n") 
        if (system ($self->tRNAscan_SE.' -q '.$self->filename." > ".$self->results)); 
}

=head2 parse_results

    Title   :  parse_results
    Usage   :   $obj->parse_results($filename)
    Function:   Parses tRNAscan_SE output to give a set of features
                parsefile can accept filenames, filehandles or pipes (\*STDIN)
    Returns :   none
    Args    :   optional filename

=cut

sub parse_results {
    my ($self) = @_;
    my $filehandle;
    #   my $resfile = $self->workdir()."/".$self->results();
    my $resfile = $self->results();
    
    if (-e $resfile) {
      if (-z $self->results) {  
	print STDERR "tRNAscan-SE didn't find any tRNAs\n";
	return; }       
      else {
        open (TRNAOUT, "<$resfile") or $self->throw("Error opening ", $resfile, " \n");#
        $filehandle = \*TRNAOUT;
      }
    }
    else { #it'a a filehandle
      $filehandle = $resfile;
    }
    
    #extract values
    while (<$filehandle>)
      {  
        next if (/^Sequence/ || /^Name/ || /^---/); #ignore introductory lines
	my @element = split (/\s+/, $_); 
	# sort out start, end & strand
	my ($start,$end,$strand);
	$start = $element[2];
	$end = $element[3];
	$strand = 1;
	if ($start > $end)
	  {
	    # is this the right way?
	    $start = $element[3];
	    $end = $element[2];
	    $strand = -1;
	  }

	my (%feature);
	$feature {name} = $element[0];
	$feature {score} = $element[8];
	$feature {start} = $start;
	$feature {end} = $end;
	$feature {strand} = $strand;
	$feature {source}= 'tRNAscan-SE';
	$feature {primary}= 'tRNA';
	$feature {program} = 'tRNAscan-SE';
	$feature {program_version} = '1.11';
	
# need the tRNA type?

	$self->create_feature(\%feature);
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
    my @list = @{$self->{'_flist'}};
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
    my $tRNA = Bio::EnsEMBL::SeqFeature->new
                        (   -seqname => $feat->{'name'},
                            -start   => $feat->{'start'},
                            -end     => $feat->{'end'},
                            -strand  => $feat->{'strand'},
                            -score   => $feat->{'score'},
                            -source_tag  => $feat->{'source'},
                            -primary_tag => $feat->{'primary'},
                            -analysis => $analysis_obj);  

    if ($tRNA)
      {
	$tRNA->validate();

	# add to _flist
	push(@{$self->{'_flist'}}, $tRNA);
      }
}

1;
