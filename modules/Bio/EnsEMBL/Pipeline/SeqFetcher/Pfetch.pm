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

=cut

# Let the code begin...
package Bio::EnsEMBL::Pipeline::SeqFetcher::Pfetch;

use strict;
use Bio::Root::RootI;
use Bio::DB::RandomAccessI;
use Bio::Seq;

use vars qw(@ISA);

@ISA = qw(Bio::Root::RootI Bio::DB::RandomAccessI);

sub new {
  my ($class, @args) = @_;
  my $self = bless {}, $class;

  my ($exe, $options) = $self->_rearrange(['EXECUTABLE', 'OPTIONS'], @args);

  if (!defined $exe) {
    $exe = 'pfetch';
  }
  $self->executable($exe);

  if (defined $options) {
    $self->options($options);
  }
  
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

#  print STDERR "$command\n";

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

  return $seq;
}

=head2 get_Seqs_by_accs

  Title   : get_Seqs_by_accs
  Usage   : $self->get_Seqs_by_accs(@accession);
  Function: Does the sequence retrieval for an array of accesions in one pfetch call
  Returns : Array of Bio::Seq
  Args    : 

=cut

sub  get_Seqs_by_accs {
  my ($self, @acc) = @_;
  
  if ((! @acc) || scalar(@acc < 1)) {
    $self->throw("No accession input");
  }  
  
  my @seq;
  my $newseq;
  my $tracker = 0;
  my $pfetch = $self->executable;
  my $options = $self->options;

  if (defined($options)) { $options = '-' . $options  unless $options =~ /^-/; }

  my $command = "$pfetch -q ";
  if (defined $options){
    $command .= "$options ";
  }
  
  $command .= join " ", @acc;
  
#  print STDERR "$command\n";
  
  open(IN,"$command |") or $self->throw("Error opening pipe to pfetch for $pfetch");
  while(<IN>){

    chomp;
    eval{
      if(defined $_ && $_ ne "no match") {
	$newseq = new Bio::Seq('-seq'               => $_,
			       '-accession_number'  => $acc[$tracker],
			       '-display_id'        => $acc[$tracker]);
      }
    };
    
    if($@){
      print STDERR "$@\n";
    }
    
    if (defined $newseq){
      push (@seq, $newseq);
    }
    else{
      $self->warn("Could not even pfetch sequence for [" . $acc[$tracker] . "]\n");
    }
    $tracker++;
  }

  close IN or $self->throw("Error running pfetch for $pfetch");
  return @seq;
}

1;
