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

ApproxAlign reads the output of an approximate alignment. It can be used to manipulate the
alignment as well as to write it to a file. ApproxAlign inherits from a module called aat.pm
(Copyright by L.Pachter,see http://baboon.math.berkeley.edu/~syntenic/slam.html or
http://math.berkeley.edu/~lpachter/). To use this module, you have to include the path to
the aat.pm-module in your PERL5LIB-variable.


=head1 CONTACT

ensembl-dev@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Pipeline::Tools::ApproxAlign;

use aat;
use strict;
use vars qw(@ISA);
use Bio::EnsEMBL::Pipeline::RunnableI;

@ISA = qw (aat Bio::EnsEMBL::Pipeline::RunnableI);

### TODO:
### Constructing new ApproxAlign-Object (inherits aat.pm)
### option for workdir can be added, Super the constructor !!
### option for filename can be added, Super the constructor !!
### needs in -ApproxAlign: new constructor
###          -run-method which reads the supplied files and the alignment
###          -DESTROY method which destroys written file
###          -filename which returns name of written file for handling over to slam


# super the constructor for aat-object and set default-values for expanding the approximate Alignment

sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);

  my (
      $fatten_value,
      $fwdonly,
      $exonbounds,
      $workdir
     ) = $self->_rearrange([qw(
                               FATTEN_VALUE
                               FWDONLY
                               EXONBOUNDS
                               WORKDIR
                              )
                           ], @args
                          );

  # setting defaults
  $self->fwdonly('0');
  $self->fatten_value(undef);

  # overriding default if optional value is passed
  $self->fwdonly($fwdonly);
  $self->fatten_value($fatten_value);

  $self->exonbounds($exonbounds);
  $self->workdir($workdir);
  $self->checkdir();

  return $self;
}

sub exonbounds {
  my ($self,$exonbounds) = @_;

  if (!defined $exonbounds && !defined $self->{_exonbounds}) {

    # default parameters for expanding the exon bounds 
    # of approximate Aligement (aatFill.pl)

      my $acc = 20;             # acc-slam-default: 20
      my $don = 10;             # don-slam-default: 10
      my $stp = 5;              # stp-slam-default: 5

      # fwd-exon-boundaries
      my $fwdStp1 = ['TAA', 0, 'zOnly', $stp];
      my $fwdStp2 = ['TAG', 0, 'zOnly', $stp];
      my $fwdStp3 = ['TGA', 0, 'zOnly', $stp];

      # rev-exon-boundaries
      my $revStp1 = ['TTA', 2, 'zOnly', $stp];
      my $revStp2 = ['CTA', 2, 'zOnly', $stp];
      my $revStp3 = ['TCA', 2, 'zOnly', $stp];

      #my $fwdDon = ['GT',  2, 'points',  $don];
      #my $revDon = ['AC', -1, 'points', -$don, -($don+1) ],
      #my $fwdAcc = ['AG', -1, 'points', -$acc, -($acc+1) ];
      #my $revAcc = ['CT',  2, 'points',  $acc];

      my $fwdDon = ['GT',  2, 'line',   $don]; # 10
      my $revDon = ['AC', -1, 'line', -($don + 1)]; # $don=10--> -11

      my $fwdAcc = ['AG', -1, 'line', -($acc + 1)]; # $acc=20--> -21
      my $revAcc = ['CT',  2, 'line',   $acc]; # 20


      if (defined $self->fwdonly && $self->fwdonly=~m/(1|true|t)/) {
        $exonbounds = [$fwdStp1, $fwdStp2, $fwdStp3, $fwdDon, $fwdAcc ];
      } else {
        $exonbounds = [$fwdStp1, $fwdStp2, $fwdStp3,$revStp1, $revStp2, $revStp3,$fwdDon, $revDon,$fwdAcc, $revAcc];
      }

      # not used because slam-default doesn't fatten
      if ($self->fatten_value) {
        push(@{$exonbounds}, ['fatten', undef, undef, $self->fatten_value]) if(defined($self->fatten_value));
      }
      $self->{_exonbounds}=$exonbounds;
    }elsif ($exonbounds) {
      $self->{_exonbounds}=$exonbounds;
    }
  return $self->{_exonbounds};
}




sub fatten_value {
    my ($self, $fatten_value) = @_;
    $self->{_fatten_value} = $fatten_value if ($fatten_value);
    return $self->{_fatten_value};
}


sub fwdonly{
  my ($self,$fwdonly) = @_;

  $self->{_fwdonly} = $fwdonly if ($fwdonly);
  return $self->{_fwdonly};
}



  # Reads in approximate alignment from supplied filename

  sub read {
    my($self, $file) = @_;

    my($i,@fields,$aatLenY,$aatLenZ);
    $i = 0;
    $self->{'aat'} = [];

    open(IN,"<$file") || die "Could not open file: $file\n";
    while (<IN>) {
      @fields = split;
      die "aat->read(): bad format in line ".($i+1).", should have 3 fields\n" if(scalar(@fields) != 3);
      die "aat->read(): expected first entry of line ".($i+1)." to be $i.\n" if($fields[0] != $i);
      push(@{$self->{'aat'}},[$fields[1],$fields[2]]);
      $self->{'weight'} += ($fields[2] - $fields[1] + 1);
      $i++;
    }
    close(IN);
    $aatLenY = scalar(@{$self->{'aat'}});
    $aatLenZ = 1 + $self->{'aat'}->[$aatLenY-1][1];
    if (defined($self->{seqY}) && ($aatLenY != scalar(@{$self->{seqY}}))) {
      printf STDERR ("Unequal lengths for aatSeqY (%d) and seqY (%d).\n",$aatLenY,scalar(@{$self->{'seqY'}}));
    }
    if (defined($self->{seqZ}) && ($aatLenZ != scalar(@{$self->{seqZ}}))) {
      printf STDERR ("Unequal lengths for aatSeqZ (%d) and seqZ (%d).\n",$aatLenZ,scalar(@{$self->{'seqZ'}}));
    }
    return;
  }


# Writes sequence to a file (filename will be created by Bio::EnsEMBL::Pipeline::RunnableI->get_tmp_filename
# if $self->workdir is set, than $workdir will be used, otherwise /tmp/ will be used
# Remember: Bio::EnsEMBL::Pipeline::RunnableI must be in your path
# write returns the name and path of the written file

sub write {
  my ($self) = @_;

  my $n = scalar(@{$self->{'aat'}});
  my $aatfile = Bio::EnsEMBL::Pipeline::RunnableI->get_tmp_file($self->workdir,"approxAlignOutput","aat");
  print "ApproxAlign: $aatfile written\n";

  open(OUT,">$aatfile") || die "Could not write to file $aatfile\n";
  for (my $i=0; $i < $n; $i++) {
    print OUT $i . " " . ${$self->{'aat'}}[$i][0] . " " . ${$self->{'aat'}}[$i][1] . "\n";
  }
  close(OUT);
  $self->aatfile($aatfile);
  return;
}



sub aatfile {
  my $self = shift;

  $self->{'aatfile'} = shift if @_ ;
  return $self->{'aatfile'}
}


1;
