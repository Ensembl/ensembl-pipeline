# Cared for by Abel Ureta-Vidal <abel@ebi.ac.uk>
#
# Copyright GRL & EBI
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::Runnable::Blastz

=head1 SYNOPSIS

  # To run a blastz job from scratch do the following.

  my $query = new Bio::Seq(-file   => 'somefile.fa',
                           -format => 'fasta');

  my $database = 'multifastafile.fa';

  my $blastz =  Bio::EnsEMBL::Pipeline::Runnable::Blastz->new 
    ('-query'     => $query,
     '-database'  => $database,
     '-options'   => 'T=2');

  $blastz->run();

  @featurepairs = $blast->output();

  foreach my $fp (@featurepairs) {
      print $fp->gffstring . "\n";
  }

  # Additionally if you have blast runs lying around that need parsing
  # you can use the EnsEMBL blastz parser module 
  # perldoc Bio::EnsEMBL::Pipeline::Runnable::Parser::Blastz


=head1 DESCRIPTION

Blastz takes a Bio::Seq (or Bio::PrimarySeq) object and runs blastz with against 
the specified multi-FASTA file database. Tthe output is parsed by 
Bio::EnsEMBL::Pipeline::Runnable::Parser::Blastz and stored as Bio::EnsEMBL::DnaDnaAlignFeature 

Other options can be passed to the blast program using the -options method

=head1 CONTACT

Describe contact details here

=head1 APPENDIX


=cut

package Bio::EnsEMBL::Pipeline::Runnable::Blastz;


use vars qw(@ISA);
use strict;

# Object preamble

use Bio::EnsEMBL::Pipeline::RunnableI;
use Bio::EnsEMBL::DnaDnaAlignFeature;
use Bio::EnsEMBL::Pipeline::Tools::Blastz;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI);

=head2 new

    Title   :   new
    Usage   :   my obj =  Bio::EnsEMBL::Pipeline::Runnable::Blast->new 
    (-query    => $seq,
     -database => $database,
     -options   => 'C=2 K=3000 H=2200');

    Function:   Initialises Blastz object
    Returns :   a Blastz Object
    Args    :   A Bio::Seq object (-query)
                A database file (-database)
                The blastz executable (-program)
                Options (-options)

=cut

sub new {
    my ($class,@args) = @_;

    my $self = $class->SUPER::new(@args);    

    $self->{'_query'}     = undef;     # location of Bio::Seq object
    $self->{'_filename'}  = undef;     # file to store Bio::Seq object
    $self->{'_program'}   = "blastz";     # location of Blast
    $self->{'_database'}  = undef;     # name of database filename
    $self->{'_options'}   = "";     # options for blastz
    $self->{'_fplist'}    = [];        # an array of feature pairs (the output)
    $self->{'_workdir'}   = "/tmp";     # location of temp directory
    $self->{'_results'}   = $self->{'_workdir'}."/results.".$$; # location of result file

    # Now parse the input options and store them in the object
    my($program,$query,$database,$options) = $self->_rearrange([qw(PROGRAM
								   QUERY 
								   DATABASE 
								   OPTIONS)], 
							       @args);

    if ($query) {
      $self->query($query);
    } else {
      $self->throw("No query sequence input.");
    }

    $self->program($self->find_executable($program)) if ($program);

    if ($database) {
      $self->database($database);
    } else {
      $self->throw("No database input");
    }
    
    if ($options) {
       $self->options($options);
    } 
    
    return $self; # success - we hope!
}

=head2 run

    Title   :  run
    Usage   :   $obj->run()
    Function:   Runs blast and BPLite and creates array of feature pairs
    Returns :   none
    Args    :   none

=cut

sub run {
  my ($self) = @_;
  
  $self->checkdir();
  
  #write sequence to file
  $self->writefile();

  $self->run_analysis();
  
  #parse output and create features
  $self->parse_results;  
  $self->deletefiles();
}


sub run_analysis {
  my ($self) = @_;
  print STDERR ("Running blastz...\n" . $self->program . " " .
		                        $self->filename . " " .
		                        $self->database . " " .
		                        $self->options . " > " .
		                        $self->results . "\n");

  $self->throw("Failed during blastz run, $!\n") unless (system ($self->program  . " " .
								 $self->filename . " " .
								 $self->database . " " .
								 $self->options . " > " .
								 $self->results) == 0);
}  
  
sub parse_results { 
  my ($self) = @_;
  
  my $BlastzParser = Bio::EnsEMBL::Pipeline::Tools::Blastz->new('-file' => $self->results);
  
  while (defined (my $alignment = $BlastzParser->nextAlignment)) { # nextHSP-like
    $self->_add_fp($alignment);
  }
}


#################
# get/set methods 
#################

sub query {
  my ($self, $seq) = @_;
  if ($seq) {
    unless ($seq->isa("Bio::PrimarySeqI") || $seq->isa("Bio::Seq")) {
      $self->throw("Input isn't a Bio::Seq or Bio::PrimarySeq");
    }
    
    $self->{'_query'} = $seq ;
    $self->filename($seq->id.".$$.seq");
  }
  return $self->{'_query'};
}

=head2 database
  
    Title   :   database
    Usage   :   $self->database($seq)
    Function:   Get/set method for database
    Returns :   Bio::PrimarySeqI, or filename
    Args    :   Bio::PrimarySeqI, or filename

=cut

sub database {
  my($self, $value) = @_;    
  
  if (defined $value) {

     if (! ref($value)) {
      # assume it's a filename - check the file exists
      $self->throw("[$value] : file does not exist\n") unless -e $value;
      $self->genfilename($value);
      $self->{'_database'} = $value;
    }
    elsif ($value->isa("Bio::PrimarySeqI")) {
      my $filename = "/tmp/genfile_$$.fn";
      $self->genfilename($filename);
      my $genOutput = Bio::SeqIO->new(-file => ">$filename" , '-format' => "Fasta")
	or $self->throw("Can't create new Bio::SeqIO from $filename '$' : $!");
    
      $self->throw ("problem writing genomic sequence to $filename\n" ) unless $genOutput->write_seq($value);
      $self->{'_database'} = $self->genfilename;
    }
    else {
      $self->throw("$value is neither a Bio::Seq  nor a filename\n");
    }
  }
  
  return $self->{'_database'};
}

=head2 options

    Title   :   options
    Usage   :   $obj->options(' -I ');
    Function:   Get/set method for blast arguments
    Args    :   File path (optional)

=cut

sub options {
  my ($self, $args) = @_;
  
  if (defined($args)) {
    $self->{'_options'} = $args ;
  }
  return $self->{'_options'};
}

=head2 output

    Title   :   output
    Usage   :   $self->output()
    Function:   Returns all output feature pairs
    Returns :   Array of Bio::EnsEMBL::FeaturePairs
    Args    :   None

=cut

sub output {
  my ($self) = @_;
  return @{$self->{'_fplist'}};
}

sub _add_fp {
  my ($self,@args) = @_;
  if (@args) {
    push(@{$self->{'_fplist'}},@args);
  } else {
    warn "WARN: Bio::EnsEMBL::Pipeline::Runnable::Blastz->_add_fp should have an argument\n";
  }
}

=head2 workdir

 Title   : workdir
 Usage   : $obj->workdir($newval)
 Function: 
 Example : 
 Returns : value of workdir
 Args    : newvalue (optional)


=cut

sub workdir{
   my ($self,$value) = @_;
   if( defined $value) {
       $self->{'_workdir'} = $value;
   }
   return $self->{'_workdir'};
}

=head2 program

    Title   :   program
    Usage   :   $obj->program('/usr/local/ensembl/bin/bl2seq');
    Function:   Get/set method for the location of the bl2seq executable
    Returns :   string
    Args    :   string

=cut

sub program {
  my ($self, $location) = @_;

  if ($location) {
    $self->throw("executable not found at $location: $!\n") unless (-e $location && -x $location);
    $self->{'_program'} = $location ;
  }
  return $self->{'_program'};
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
