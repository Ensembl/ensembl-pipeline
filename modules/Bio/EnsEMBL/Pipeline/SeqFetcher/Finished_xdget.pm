=pod

=head1 NAME

Bio::EnsEMBL::Pipeline::SeqFetcher::Finished_xdget

=head1 SYNOPSIS

  my $obj = Bio::EnsEMBL::Pipeline::SeqFetcher::xdget->new(
      -executable => '/blah/xdget',
      -db         => '/data/db'
  );
  my $seq = $obj->get_Seq_by_accs(\@acc);

=head1 DESCRIPTION

Object to retrieve sequences using xdget (Wash U).
Database must be formatted with xdformat. Sequence type (protein
or nucleotide) is guessed, based on file extensions of the database
files. Returns hashref with accession/Bio::Seq as key/value.
(Undef value if accession not found)

Additional options for xdget can be specifed though no checking
is performed for compatibility.

Note that, at the time of writing, xdget is case-insensitive:
retrieved sequence is in upper case, irrespective of what was
in the original unformatted fasta file.

=head1 CONTACT

B<anacode@sanger.ac.uk>

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


package Bio::EnsEMBL::Pipeline::SeqFetcher::Finished_xdget;

use strict;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Root;
use Bio::DB::RandomAccessI;
use Bio::Seq;

use base 'Bio::EnsEMBL::Pipeline::SeqFetcher::xdget';

my $debug = 0;

=head2 get_Seq_by_accs

  Title   : get_Seq_by_accs
  Usage   : $self->get_eq_by_acc($accessions);
  Function: retrieves sequence via xdget
  Returns : hashref of Bio::Seq
  Args    : array ref of Sequence identifiers

=cut

sub get_Seq_by_accs {
  my ($self, $accs) = @_;

  throw("Accessions input is not an ARRAY") unless ref($accs) eq 'ARRAY';
  throw("No database defined") unless $self->db;


  my $xdget   = $self->executable;
  my $db      = $self->db;
  local       *FH;
  my %acc_hash = map{ $_ => undef } @$accs;


  # maybe should have some checking here to see if -n/-p have
  # already been specified in options

  DB: foreach my $db (@{$self->db}) {
  	my $acc;
  	my $seq;
	my $seqstr;
	my $desc;
	my $command;
    my $options = $self->options;

    if ($self->_moltype($db) eq 'n') {
      $options .= " -n";
    }
    else {
      $options .= " -p";
    }

    $command = "$xdget $options $db ".join(" ",@$accs);
    print STDOUT $command."\n" if $debug;
    undef $/;           # enable "slurp" mode
    open FH, "$command 2> /dev/null |" or throw("Error retrieving ".scalar(@$accs)." accessions from $db with $xdget");
	my $out = <FH>;     # whole output now here

    BLOCK:foreach(split(/\n>/,$out)) {
      my @rows = split(/\n/,$_);
      $desc = shift @rows;
      $desc =~ s/^>//;
      # EMBL header
      # >AF001541.1 Homo sapiens clone alpha_S628 mRNA sequence.
      ($acc) = $desc =~ /^(\w+\.\w+)/;
      # RefSeq header
      # >gi|8923664|ref|NM_017949.1| Homo sapiens CUE domain containing 1 (CUEDC1), mRNA
      ($acc) = $desc =~ /\|(\w+_\w+\.\w+)\|/ unless $acc;
	  $seqstr = join('',@rows);
	  $seqstr =~ s/\s//g;
	  print STDOUT "Making Bio::Seq accession[$acc] description[$desc] seq[".length($seqstr)."]\n" if $debug;
	  $seq = Bio::Seq->new(
			    -seq              => $seqstr,
			    -display_id       => $acc,
			    -accession_number => $acc,
			    -desc             => $desc
	  );
	  $acc_hash{$acc} = $seq;

    }
    close FH;
    $/ = "\n";

	$accs = [];
	foreach(keys %acc_hash) {
		push @$accs, $_ unless $acc_hash{$_};
	}
	last DB unless @$accs;
  }

  return \%acc_hash;
}

1;
