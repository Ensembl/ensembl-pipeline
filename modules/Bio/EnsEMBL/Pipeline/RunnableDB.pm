#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::RunnableDB

=head1 SYNOPSIS

# get a Bio::EnsEMBL::Pipeline::RunnableDB pbject somehow

  $runnabledb->fetch_input();
  $runnabledb->run();
  $runnabledb->output();
  $runnabledb->write_output(); #writes to DB

=head1 DESCRIPTION

This is the base implementation of
This object encapsulates the basic main methods of a RunnableDB
which a subclass may override.

parameters to new
-db:        A Bio::EnsEMBL::DBSQL::DBAdaptor (required), 
-input_id:   Contig input id (required), 
-analysis:  A Bio::EnsEMBL::Analysis (optional) 

This object wraps Bio::EnsEMBL::Pipeline::Runnable to add
functionality for reading and writing to databases.  The appropriate
Bio::EnsEMBL::Analysis object must be passed for extraction
of appropriate parameters. A Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor is
required for databse access.

=head1 CONTACT

ensembl-dev@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Pipeline::RunnableDB;

use strict;
use Bio::EnsEMBL::Pipeline::SeqFetcher;
use Bio::EnsEMBL::Pipeline::SeqFetcher::Pfetch;
use Bio::EnsEMBL::Pipeline::RunnableI;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::DB::RandomAccessI;

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI);

=head2 new

    Title   :   new
    Usage   :   $self->new(-DB          => $db,
                           -INPUT_ID    => $id,
                           -SEQFETCHER  => $sf,
			   -ANALYSIS    => $analysis);

    Function:   creates a Bio::EnsEMBL::Pipeline::RunnableDB object
    Returns :   A Bio::EnsEMBL::Pipeline::RunnableDB object
    Args    :   -db:         A Bio::EnsEMBL::DBSQL::DBAdaptor (required), 
                -input_id:   Contig input id (required), 
                -seqfetcher: A Bio::DB::RandomAccessI Object (required),
                -analysis:   A Bio::EnsEMBL::Analysis (optional) 
=cut

sub new {
    my ($class, @args) = @_;

    my $self = {};
    bless $self, $class;

    my ($db,$input_id, $seqfetcher, $analysis, $parameters) = 
    $self->_rearrange([qw(DB INPUT_ID	SEQFETCHER ANALYSIS PARAMETERS)],@args);

    $self->{'_genseq'}      = undef;
    $self->{'_runnable'}    = undef;
    $self->{'_parameters'}  = undef;
    $self->{'_analysis'}    = undef;
    $self->{job_parameters} = undef;

    if($parameters && (!$db || !$analysis)){
      $self->job_parameters($parameters);
      my ($host, $port, $user, $pass, $dbname, $logic_name) = 
        split /:/, $parameters;
      $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host => $host,
					       -user => $user,
					       -pass => $pass,
					       -dbname => $dbname,
                 -port => $port);
      my $ana_adp = $db->get_AnalysisAdaptor;
      $analysis = $ana_adp->fetch_by_logic_name($logic_name);
    }
    $self->throw("No database handle input") unless defined($db);
    $self->db($db);

    $self->throw("No input id input")        unless defined($input_id);
    $self->input_id($input_id);

    # we can't just default this to pfetch
    $seqfetcher && $self->seqfetcher($seqfetcher);

    $self->throw("No analysis object input") unless defined($analysis);
    $self->analysis($analysis);

    return $self;
}

=head2 analysis

    Title   :   analysis
    Usage   :   $self->analysis($analysis);
    Function:   Gets or sets the stored Analusis object
    Returns :   Bio::EnsEMBL::Analysis object
    Args    :   Bio::EnsEMBL::Analysis object

=cut

sub analysis {
    my ($self, $analysis) = @_;
    
    if ($analysis) {
        $self->throw("Not a Bio::EnsEMBL::Analysis object")
            unless ($analysis->isa("Bio::EnsEMBL::Analysis"));
        $self->{'_analysis'} = $analysis;
        $self->parameters($analysis->parameters);
    }
    return $self->{'_analysis'};
}

=head2 parameters

    Title   :   parameters
    Usage   :   $self->parameters($param);
    Function:   Gets or sets the value of parameters
    Returns :   A string containing parameters for Bio::EnsEMBL::Runnable run
    Args    :   A string containing parameters for Bio::EnsEMBL::Runnable run

=cut

sub parameters {
    my ($self, $parameters) = @_;

    $self->analysis->parameters($parameters) if ($parameters);


    return $self->analysis->parameters();
}

sub arguments {
  my ($self) = @_;

  my %parameters = $self->parameter_hash;

  my $options = "";

  foreach my $key (keys %parameters) {
    if ($parameters{$key} ne "__NONE__") {
      $options .= " " . $key . " " . $parameters{$key};
    } else {
      $options .= " " . $key;
    }
  }
  return $options;
}

=head2 parameter_hash

 Title   : parameter_hash
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub parameter_hash{
   my ($self,@args) = @_;

    my ($parameter_string) = $self->analysis->parameters() ;

    my %parameters;

    if ($parameter_string) {

      my @pairs = split (/,/, $parameter_string);
      foreach my $pair (@pairs) {
	
	my ($key, $value) = split (/=>/, $pair);

	if ($key && $value) {
	  $key   =~ s/^\s+//g;
	  $key   =~ s/\s+$//g;
	  $value =~ s/^\s+//g;
	  $value =~ s/\s+$//g;
	  
	  $parameters{$key} = $value;
	} else {
          $parameters{$key} = "__NONE__";
	}
      }
    }
    return %parameters;
}

=head2 db

    Title   :   db
    Usage   :   $self->db($obj);
    Function:   Gets or sets the value of db
    Returns :   A Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor compliant object
                (which extends Bio::EnsEMBL::DBSQL::DBAdaptor)
    Args    :   A Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor compliant object

=cut

sub db {
    my( $self, $value ) = @_;

    if ($value) {
       $value->isa("Bio::EnsEMBL::DBSQL::DBConnection")
         || $self->throw("Input [$value] isn't a Bio::EnsEMBL::DBSQL::DBConnection");

       $self->{'_db'} = $value;
    }
    return $self->{'_db'};
}

=head2 input_id

    Title   :   input_id
    Usage   :   $self->input_id($input_id);
    Function:   Gets or sets the value of input_id
    Returns :   valid input id for this analysis (if set) 
    Args    :   input id for this analysis 

=cut

sub input_id {
    my ($self, $input) = @_;

    if ($input) {
        $self->{'_input_id'} = $input;
    }

    return $self->{'_input_id'};
}

=head2 query

    Title   :   query
    Usage   :   $self->query($query);
    Function:   Get/set query
    Returns :   
    Args    :   

=cut

sub query {
    my ($self, $query) = @_;

    if (defined($query)){ 
	$self->{'_query'} = $query; 
    }

    return $self->{'_query'}
}

=head2 output

    Title   :   output
    Usage   :   $self->output()
    Function:   
    Returns :   Array of Bio::EnsEMBL::FeaturePair
    Args    :   None

=cut

sub output {
    my ($self) = @_;
   
    $self->{'_output'} = [];
    
    my @r = $self->runnable;

    if(@r && scalar(@r)){
      foreach my $r ($self->runnable){
	push(@{$self->{'_output'}}, $r->output);
      }
    }
    return @{$self->{'_output'}};
}

=head2 run

    Title   :   run
    Usage   :   $self->run();
    Function:   Runs Bio::EnsEMBL::Pipeline::Runnable::xxxx->run()
    Returns :   none
    Args    :   none

=cut

sub run {
    my ($self) = @_;

    foreach my $runnable ($self->runnable) {

      $self->throw("Runnable module not set") unless ($runnable);

      # Not sure about this
      $self->throw("Input not fetched")       unless ($self->query);

      $runnable->run();
    }
    return 1;
}

=head2 runnable

    Title   :   runnable
    Usage   :   $self->runnable($arg)
    Function:   Sets a runnable for this RunnableDB
    Returns :   Bio::EnsEMBL::Pipeline::RunnableI
    Args    :   Bio::EnsEMBL::Pipeline::RunnableI

=cut

sub runnable {
  my ($self,$arg) = @_;

  if (!defined($self->{'_runnables'})) {
      $self->{'_runnables'} = [];
  }
  
  if (defined($arg)) {

      if ($arg->isa("Bio::EnsEMBL::Pipeline::RunnableI")) {
	  push(@{$self->{'_runnables'}},$arg);
      } else {
	  $self->throw("[$arg] is not a Bio::EnsEMBL::Pipeline::RunnableI");
      }
  }
  
  return @{$self->{'_runnables'}};  
}


=head2 write_output

    Title   :   write_output
    Usage   :   $self->write_output
    Function:   Writes output data to db
    Returns :   array of repeats (with start and end)
    Args    :   none

=cut

sub write_output {
    my ($self) = @_;

    my $db       = $self->db();
    my @features = $self->output();
    my $contig;

    my $sim_adp  = $self->db->get_SimpleFeatureAdaptor;
    my $mf_adp   = $self->db->get_MarkerFeatureAdaptor;
    my $pred_adp = $self->db->get_PredictionTranscriptAdaptor;
    my $dna_adp  = $self->db->get_DnaAlignFeatureAdaptor;
    my $rept_adp = $self->db->get_RepeatFeatureAdaptor;
    my $pep_adp  = $self->db->get_ProteinAlignFeatureAdaptor;
    my $gene_adp = $self->db->get_GeneAdaptor;


    $self->warn("shouldn't be using the write_output method of runnabledb");

    eval {
      $contig = $db->get_RawContigAdaptor->fetch_by_name($self->input_id);
    };

    if ($@) {
      $self->throw("Can't find contig " . $self->input_id . " . Can't write output");
    }
  
    my %features;

    foreach my $f (@features) {

      $f->analysis($self->analysis);

      unless ($f->isa("Bio::EnsEMBL::PredictionTranscript")) {
        print $f->gffstring . "\n";
        $f->attach_seq($contig);
      }

      if ($f->isa("Bio::EnsEMBL::PredictionTranscript")) {
	foreach my $exon (@{$f->get_all_Exons}) {
	  $exon->contig($contig);
	}

	if (!defined($features{prediction})) {
	  $features{prediction} = [];
	}

	push(@{$features{prediction}},$f);

	print "F " . $features{prediction} . "\n";

      } elsif ($f->isa("Bio::EnsEMBL::SimpleFeature")) {

	if (!defined($features{simple})) {
	  $features{simple} = [];
	}

	push(@{$features{simple}},$f);

      } elsif ($f->isa("Bio::EnsEMBL::Map::MarkerFeature")) {

	if (!defined($features{marker})) {
	  $features{marker} = [];
	}

	push(@{$features{marker}},$f);

      } elsif ($f->isa("Bio::EnsEMBL::DnaPepAlignFeature")) {

	if (!defined($features{dnapep})) {
	  $features{dnapep} = [];
	}

	push(@{$features{dnapep}},$f);

      } elsif ($f->isa("Bio::EnsEMBL::DnaDnaAlignFeature")) {

	if (!defined($features{dnadna})) {
	  $features{dnadna} = [];
	}

	push(@{$features{dnadna}},$f);

      } elsif ($f->isa("Bio::EnsEMBL::RepeatFeature")) {

	if (!defined($features{repeat})) {
	  $features{repeat} = [];
	}

	push(@{$features{repeat}},$f);

      } elsif ($f->isa("Bio::EnsEMBL:Gene")) {

	  foreach my $exon (@{$f->get_all_Exons}) {
	    $exon->contig($contig);
	  }
	  if (!defined($features{gene})) {
	    $features{gene} = [];
	  }
	  
	push(@{$features{gene}},$f);
	  
      }
    }

    if ($features{prediction}) {
      $pred_adp->store(@{$features{prediction}});
      print "Storing " . @{$features{prediction}} . "\n";
    }
    if ($features{simple}) {
      $sim_adp->store(@{$features{simple}});
    }
    if ($features{marker}) {
      $mf_adp->store(@{$features{marker}});
    }
    if ($features{dnadna}) {
      $dna_adp->store(@{$features{dnadna}});
    }
    if ($features{dnapep}) {
      $pep_adp->store(@{$features{dnapep}});
    }
    if ($features{repeat}) {
      $rept_adp->store(@{$features{repeat}});
    }
    if ($features{gene}) {
      $gene_adp->store(@{$features{gene}});
    }

    return 1;
}

=head2 seqfetcher

    Title   :   seqfetcher
    Usage   :   $self->seqfetcher($seqfetcher)
    Function:   Get/set method for SeqFetcher
    Returns :   Bio::DB::RandomAccessI object
    Args    :   Bio::DB::RandomAccessI object

=cut

sub seqfetcher {
  my( $self, $value ) = @_;    

  if (defined($value)) {

    #need to check if passed sequence is Bio::DB::RandomAccessI object
    #$value->isa("Bio::DB::RandomAccessI") || 
    #  $self->throw("Input isn't a Bio::DB::RandomAccessI");

    $self->{'_seqfetcher'} = $value;
  }
    return $self->{'_seqfetcher'};
}

=head2 input_is_void

    Title   :   input_is_void
    Usage   :   $self->input_is_void(1)
    Function:   Get/set flag for sanity of input sequence
                e.g. reject seqs with only two base pairs
    Returns :   Boolean
    Args    :   Boolean

=cut

sub input_is_void {
    my ($self, $value) = @_;

    if ($value) {
	$self->{'_input_is_void'} = $value;
    }
    return $self->{'_input_is_void'};

}

sub job_parameters{
  my $self = shift;

  if(@_){
    $self->{job_parameters} = shift;
  }

  return $self->{job_parameters};
  
}






1;
