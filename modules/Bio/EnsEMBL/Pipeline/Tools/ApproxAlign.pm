#
# Written by Jan-Hinnerk Vogel
#
# jhv [at] sanger.ac.uk
#
# Copyright GRL/EBI 2004
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod

=head1 NAME

Bio::EnsEMBL::Pipeline::Runnable::Tools::ApproxAlign

=head1 SYNOPSIS
$slice1       = the first Bio::Seq object,
$slice2       = the second Bio::Seq object,


my $obj = Bio::EnsEMBL::Pipeline::Runnable::Tools::ApproxAlign->new(
                                                      -slice1        => $slice1,
                                                      -slice2        => $slice2,
                                                     );

=head1 DESCRIPTION

ApproxAlign reads the output of an approximate alignment. It can be used to manipulate the alignment as well as 
to write it to a file. ApproxAlign inherits it's methods from a module called aat.pm (Copyright by L.Pachter, 
see http://baboon.math.berkeley.edu/~syntenic/slam.html or http://math.berkeley.edu/~lpachter/).


=head1 CONTACT

ensembl-dev@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Pipeline::Tools::ApproxAlign;

use vars qw(@ISA);
use strict;
use aat;

 @ISA = qw (aat Bio::EnsEMBL::Pipeline::RunnableI);




### TODO:
### Constructing new ApproxAlign-Object (inherits aat.pm)
### option for workdir can be added, Super the constructor !!
### option for filename can be added, Super the constructor !!
### needs in -ApproxAlign: new constructor
###          -run-method which reads the supplied files and the alignment
###          -DESTROY method which destroys written file
###          -filename which returns name of written file for handling over to slam


# Reads in approximate alignment from supplied filename

sub read {
  my($this, $file) = @_;

  my($i,@fields,$aatLenY,$aatLenZ);
  $i = 0;
  $this->{'aat'} = [];

  open(IN,"<$file") || die "Could not open file: $file\n";
  while(<IN>) {
    @fields = split;
    die "aat->read(): bad format in line ".($i+1).", should have 3 fields\n" if(scalar(@fields) != 3);
    die "aat->read(): expected first entry of line ".($i+1)." to be $i.\n" if($fields[0] != $i);
    push(@{$this->{'aat'}},[$fields[1],$fields[2]]);
    $this->{'weight'} += ($fields[2] - $fields[1] + 1);
    $i++;
  }
  close(IN);
  $aatLenY = scalar(@{$this->{'aat'}});
  $aatLenZ = 1 + $this->{'aat'}->[$aatLenY-1][1];
  if(defined($this->{seqY}) && ($aatLenY != scalar(@{$this->{seqY}}))) {
    printf STDERR ("Unequal lengths for aatSeqY (%d) and seqY (%d).\n",$aatLenY,scalar(@{$this->{'seqY'}}));
  }
  if(defined($this->{seqZ}) && ($aatLenZ != scalar(@{$this->{seqZ}}))) {
    printf STDERR ("Unequal lengths for aatSeqZ (%d) and seqZ (%d).\n",$aatLenZ,scalar(@{$this->{'seqZ'}}));
  }
  return;
}


# Writes sequence to a file (filename will be created by Bio::EnsEMBL::Pipeline::RunnableI->get_tmp_filename
# Bio::EnsEMBL::Pipeline::RunnableI must be in your path
# and return the filename of the written file

sub write {
  my($this, $out) = @_;
  my($n,$i);

  $n = scalar(@{$this->{'aat'}});
  open(OUT,">$out") || die "Could not write to file $out\n";

  for($i=0; $i < $n; $i++) {
    print OUT $i . " " . ${$this->{'aat'}}[$i][0] . " " . ${$this->{'aat'}}[$i][1] . "\n";
  }
 close(OUT);
  return;
}


1;
