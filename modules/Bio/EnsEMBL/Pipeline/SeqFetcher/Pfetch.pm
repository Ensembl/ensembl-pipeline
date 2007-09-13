#
# Cared for by EnsEMBL  <ensembl-dev@ebi.ac.uk>
#
# Copyright GRL & EBI
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

  Bio::EnsEMBL::Pipeline::SeqFetcher::Pfetch

=head1 SYNOPSIS

    my $obj = Bio::EnsEMBL::Pipeline::SeqFetcher::Pfetch->new(
							      -executable => $exe
							     );
    my $seq = $obj->get_Seq_by_acc($acc);

=head1 DESCRIPTION

  Object to retrieve sequences as Bio::Seq, using pfetch.

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

Method Bio::EnsEMBL::Root::_rearrange is deprecated.
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
rearrange(order, list); #instead

=cut

# Let the code begin...
package Bio::EnsEMBL::Pipeline::SeqFetcher::Pfetch;

use strict;
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::DB::RandomAccessI;
use Bio::Seq;

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Root Bio::DB::RandomAccessI);

sub new {
  my ($class, @args) = @_;
  my $self = bless {}, $class;

  my ($exe, $options) = rearrange(['EXECUTABLE', 'OPTIONS'], @args);
  
  if (!defined $exe) {
    $exe = 'pfetch';
  }
  $self->executable($exe);

  if (defined $options) {
    $self->options($options);
  } 
  # caching of sequences 
  $self->{_seq_cache}={};    

  return $self; # success - we hope!
}

=head2 executable

  Title   : executable
  Usage   : $self->executable('/path/to/executable');
  Function: Get/set for the path to the executable being used by the module. If not set, the executable is looked for in $PATH.
  Returns : string
  Args    : string

=cut

sub executable {
  my ($self, $exe) = @_;
  if ($exe)
    {
      $self->{'_exe'} = $exe;
    }
  return $self->{'_exe'};  
}

=head2 options

  Title   : options
  Usage   : $self->options('tc');
  Function: Get/set for options to pfetch
  Returns : string
  Args    : string

=cut

sub options {

  my ($self, $options) = @_;
  if ($options)
    {
      $self->{'_options'} = $options;
    }
  return $self->{'_options'};  

}

=head2 get_Seq_by_acc

  Title   : get_Seq_by_acc
  Usage   : $self->get_eq_by_acc($accession);
  Function: Does the sequence retrieval via pfetch
  Returns : Bio::Seq
  Args    : 

=cut

sub  get_Seq_by_acc {
  my ($self, $acc) = @_;
  
  if (!defined($acc)) {
    $self->throw("No accession input");
  }  
 
  if ( defined ( $self->{_seq_cache}{$acc})) {  
     return $self->{_seq_cache}{$acc};
  } 
  # seqence not cached, needs to be fetched 
  my $seqstr;
  my $seq;
  my $pfetch = $self->executable;
  my $options = $self->options;
  if (defined($options)) { $options = '-' . $options  unless $options =~ /^-/; }

  my $command = "$pfetch -q ";
  if (defined $options){
    $command .= "$options ";
  }

  $command .= $acc;

  #print STDERR "$command\n";

  open(IN,"$command |") or $self->throw("Error opening pipe to pfetch for accession [$acc]: $pfetch");
  $seqstr = <IN>;
  close IN or $self->throw("Error running pfetch for accession [$acc]: $pfetch");
  
  chomp($seqstr);
  eval{
    if(defined $seqstr && $seqstr ne "no match") {
      $seq = new Bio::Seq('-seq'               => $seqstr,
			  '-accession_number'  => $acc,
			  '-display_id'        => $acc);
    }
  };

  if($@){
    print STDERR "$@\n";
  }
  
  $self->throw("Could not pfetch sequence for [$acc]\n") unless defined $seq;

  $self->{_seq_cache}{$acc}=$seq;
  return $seq;
}


=head2 batch_fetch

  Title   : batch_fetch
  Usage   : $self->batch_retrieval(@accession_list);
  Function: Retrieves multiple sequences via pfetch
  Returns : reference to a list of Bio::Seq objects
  Args    : array of accession strings

=cut

sub  batch_fetch {
  my $self = shift @_;

  my @sequence_list;
  
  unless (scalar @_) {
    $self->throw("No accession input");
  }  

  my $accession_concatenation;
  my @accession_list;

  while (my $acc = shift @_) {
    push (@accession_list, $acc);
    $accession_concatenation .= $acc . ' ';
  }
  
  my $pfetch = $self->executable;
  my $options = $self->options;
  if (defined($options)) { $options = '-' . $options  unless $options =~ /^-/; }

  my $command = "$pfetch -q ";
  if (defined $options){
    $command .= "$options ";
  }

  $command .= $accession_concatenation;

  #print STDERR "$command\n";

  open(IN,"$command |") or $self->throw("Error opening pipe to pfetch : $pfetch");

  my $seq_placemarker = -1;

 SEQ:
  while (my $seqstr = <IN>) {
    
    $seq_placemarker++;

    chomp($seqstr);
    
    my $seq;
    
    unless (defined $seqstr && $seqstr eq 'no match') {

      eval{
	if(defined $seqstr && $seqstr ne "no match") {
	  $seq = new Bio::Seq('-seq'               => $seqstr,
			      '-accession_number'  => $accession_list[$seq_placemarker],
			      '-display_id'        => $accession_list[$seq_placemarker]);
	}
      };
      
      if($@){
	print STDERR "$@\n";
      }
    }

    unless (defined $seq){
      $self->warn("PFetch Error : Could not pfetch sequence for " . 
		  $accession_list[$seq_placemarker] . "\n");
      next SEQ;
    }

    push (@sequence_list, $seq);
  }
  
  close IN or $self->throw("Error running pfetch : $pfetch");
  
  return \@sequence_list;
}



1;
