#
#
# Cared for by Abel Ureta-Vidal <abel@ebi.ac.uk>
#
# Copyright Abel Ureta-Vidal
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::Runnable::Promoterwise

=head1 SYNOPSIS

  my $query_seqfile = '/path/to/query.fasta';
  my $target_seqfile = '/path/to/target.fasta';

  #create Bio::EnsEMBL::Pipeline::Runnable::Promoterwise object
  my $pw = Bio::EnsEMBL::Pipeline::Runnable::Promoterwise->new (-QUERY => $query,
                                                                -TARGET => $target
                                                                -OPTIONS => '-lhreject both');
  $pw->workdir($workdir);
  $pw->run();
  my @results = $pw->output();

=head1 DESCRIPTION

Promoterwise takes 2 fasta files paths, query and target respectively,
and an option string. The output is parsed to procudes a set of DnaDnaAlignFeature pairs.
Options can be passed to Promoterwise through the options() method.

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Pipeline::Runnable::Promoterwise;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Pipeline::Tools::Promoterwise;
use Bio::EnsEMBL::Pipeline::RunnableI;
use Bio::EnsEMBL::Analysis; 
use Bio::Seq;
use Bio::SeqIO;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI);


=head2 new

  Arg [1]    : Bio::PrimarySeq $query
  Arg [2]    : Bio::PrimarySeq $target
  Arg [3]    : string $options (optional)
               e.g. "-lhreject both", default is ""
  Arg [4]    : string $program (optional)
               e.g. "/usr/local/ensembl/bin/promoterwise"
  Arg [5]    : string $workdir (optional)
               e.g. "/my/dir", default is /tmp
  Example    : my $pw = new Bio::EnsEMBL::Pipeline::Runnable::Promoterwise(-QUERY => $query,
                                                                           -TARGET => $target,
                                                                           -OPTIONS => $options);
  Description: Initialises Promoterwise object
  Returntype : Bio::EnsEMBL::Pipeline::Runnable::Promoterwise
  Exceptions : 
  Caller     : 

=cut

sub new {
  my ($class,@args) = @_;
  
  my $self = $class->SUPER::new(@args);    
  
  $self->{'_fplist'} = [];      #an array of feature pairs
  $self->{'_query'} = undef;    #location of query sequence file
  $self->{'_target'} = undef;   #location of target sequence file
  $self->{'_program'} = "promoterwise";  #location of promoterwise executable
  $self->{'_workdir'} = "/tmp"; #location of temp directory
  $self->{'_results'} = $self->{'_workdir'}."/pwresults.".$$;  # location of result file
  $self->{'_options'} = "";  #options for Promoterwise
  
  my ($query, $target, $options, $program, $workdir) = $self->_rearrange([qw(QUERY
                                                                             TARGET
                                                                             OPTIONS
                                                                             PROGRAM
                                                                             WORKDIR)], 
                                                               @args);
  
  $self->program($self->find_executable($program)) if ($program);

  # check that query and target exists

  $self->query($query) if ($query);
  $self->target($target) if ($target);
  $self->options($options) if ($options);
  if ($workdir) {
    $self->workdir($workdir);
    $self->results($self->workdir."/pwresults.".$$);
  }

  return $self;
}

=head2 query

  Arg [1]    : 
  Example    : 
  Description: 
  Returntype : 
  Exceptions : 
  Caller     : 

=cut

sub query {
  my $self = shift;

  if (@_) {
    my $ref = shift;
    if (! ref $ref || ! $ref->isa('Bio::Seq')) {
      $self->throw("$ref is not a PrimarySeq object.");
    }
    if ($ref->length == 0) {
      $self->throw("attempting to promoterwise seemingly 0 length sequence!");
    }
    
    $self->{'_query'} = $ref;
  }
  return $self->{'_query'};
}

=head2 target

  Arg [1]    : 
  Example    : 
  Description: 
  Returntype : 
  Exceptions : 
  Caller     : 

=cut

sub target {
  my $self = shift;

  if (@_) {
    my $ref = shift;
    if (! ref $ref || ! $ref->isa('Bio::Seq')) {
      $self->throw("$ref is not a PrimarySeq object.");
    }
    if ($ref->length == 0) {
      $self->throw("attempting to promoterwise seemingly 0 length sequence!");
    }
    
    $self->{'_target'} = $ref;
  }

  return $self->{'_target'};
}

=head2 results

  Arg [1]    : 
  Example    : 
  Description: 
  Returntype : 
  Exceptions : 
  Caller     : 

=cut

sub results {
  my $self = shift;
  $self->{'_results'} = shift if(@_);
  return $self->{'_results'};
}


=head2 program

  Arg [1]    : 
  Example    : 
  Description: 
  Returntype : 
  Exceptions : 
  Caller     : 

=cut


sub program {
  my ($self, $location) = @_;

  if ($location) {
    $self->throw("executable not found at $location: $!\n") unless (-e $location && -x $location);
    $self->{'_program'} = $location ;
  }
  return $self->{'_program'};
}

=head2 options

  Arg [1]    : 
  Example    : 
  Description: 
  Returntype : 
  Exceptions : 
  Caller     : 

=cut

sub options {
  my $self = shift;
  $self->{'_options'} = shift if(@_);
  return $self->{'_options'};
}

=head2 run

  Arg [1]    : 
  Example    : 
  Description: 
  Returntype : 
  Exceptions : 
  Caller     : 

=cut

sub run {
  my ($self) = @_;
  $self->checkdir;
  $self->run_analysis;
  $self->parse_results;
  $self->deletefiles;

  return 1;
}

=head2 run_analysis

  Arg [1]    : 
  Example    : 
  Description: 
  Returntype : 
  Exceptions : 
  Caller     : 

=cut

sub run_analysis {
    my ($self) = @_;

    my $query_file = $self->write_sequence_to_file($self->query);
    $self->file($query_file);
    my $target_file = $self->write_sequence_to_file($self->target);
    $self->file($target_file);
 
    print STDERR "Running promoterwise...\n" . $self->program . " $query_file $target_file " . $self->options . " > " . $self->results . "\n";

    $self->throw("Failed during promoterwise run, $!\n") unless (system ($self->program .
                                                                         " $query_file $target_file " .
                                                                         $self->options .
                                                                         "> " .
                                                                         $self->results) == 0);
    $self->file($self->results);
}

=head2 parse_results

  Arg [1]    : 
  Example    : 
  Description: 
  Returntype : 
  Exceptions : 
  Caller     : 

=cut

sub parse_results {
  my ($self) = @_;
  
  open PW, $self->results || 
    $self->throw("Coudn't open file ".$self->results.", $!\n");
  my $filehandle = \*PW;
  
  my $pw_parser = Bio::EnsEMBL::Pipeline::Tools::Promoterwise->new('-fh' => $filehandle,
                                                                   '-query_id' => $self->query->id,
                                                                   '-target_id' => $self->target->id);
  
  while (my $DnaDnaAlignFeature = $pw_parser->nextAlignment) {
    push @{$self->{'_fplist'}}, $DnaDnaAlignFeature;
  }
}

=head2 output

  Arg [1]    : None
  Example    : $pw->output
  Description: 
  Returntype : a array reference containing Bio::EnsEMBL::DnaDnaAlginFeature objects
  Exceptions : 
  Caller     : 

=cut

sub output {
    my ($self) = @_;
    return $self->{'_fplist'};
}

1;
