# EnsEMBL module for parsin blastz output
#
# Cared for by Abel Ureta-Vidal <abel@ebi.ac.uk>
#
# You may distribute this module under the same terms as perl itself
#

=pod

=head1 NAME

Bio::EnsEMBL::Pipeline::Tools::Blastz - Ensembl specific blastz output parser

=head1 SYNOPSIS

open F,"blastz_output_file";

my $fh = \*F;

my $BlastzParser = new Bio::EnsEMBL::Pipeline::Tools::Blastz(-fh => $fh);
 
or

my $blastz_output_file = "blastz_output_file";

my $BlastzParser = new Bio::EnsEMBL::Pipeline::Tools::Blastz(-file => $blastz_output_file);

while (defined (my $alignment = $BlastzParser->nextAlignment)) {
  print $alignment->percent_id," ",$alignment->score," ",$alignment->cigar_string,"\n";
  print $alignment->start," ",$alignment->end," ",$alignment->strand,"\n";
  print $alignment->hstart," ",$alignment->hend," ",$alignment->hstrand,"\n";
}

close F;

The constructor only need a filehandle opened on a blastz output file.
nextAlignment method return a Bio::EnsEMBL::DnaDnaAlignFeature object
corresponding to the next HSP-like alignment.

=head1 CONTACT

Ensembl development mailing list <ensembl-dev@ebi.ac.uk>
Abel Ureta-Vidal <abel@ebi.ac.uk>

=head1 APPENDIX

The rest of the documentation deals wtih each of the object methods.
Internal methods are usually preceded by a _

=cut 


package Bio::EnsEMBL::Pipeline::Tools::Blastz;

use vars qw(@ISA);
use Bio::EnsEMBL::DnaDnaAlignFeature;
use strict;

@ISA = qw(Bio::EnsEMBL::Root);

sub new {
  my ($class,@args) = @_;
  my $self = bless {}, $class;

  $self->{'_fh'} = undef; # filehandle on results file
  $self->{'_file'} = undef; # path for a results file
  $self->{'_eof'} = 0; # indicate if end of file and fh closed
  $self->{'_parsing_initialized'} = 0;
  $self->{'_command_line'} = "";
  $self->{'_matrix'} = "";
  $self->{'_options'} = "";
  
  
  my ($fh,$file) = $self->_rearrange([qw(FH FILE)], @args);

  if ((defined $fh && defined $file) ||
      !(defined $fh || defined $file)){ 
    $self->throw("Must pass in either fh or file argument");
  }
  if (defined $fh) {
    $self->{'_fh'} = $fh;
  } else {
    $self->file($file);
    open F, $self->file;
    $self->{'_fh'} = \*F;
  }
  $self->_initialize;
  return $self;
}

sub _initialize {
  my ($self) = @_;

  return undef if ($self->eof);

  my $fh = $self->fh;

  while (defined (my $line = <$fh>)) {
    next if ($line =~ /^\#:lav$/);
    last if ($line =~ /^\}$/);
    
# d stanza
    if ($line =~ /^d\s\{$/) {
      my $command_line = <$fh>;
      chomp $line;
      $command_line =~ s/^\s+\"//;
      $self->command_line($command_line);
      next;
    }
    if ($line =~ /^.*,.*$/) {
      $line =~ s/\"//g;
      $self->options($self->options.$line);
    } else {
      $self->matrix($self->matrix.$line);
    }
  }
  return $self->_parsing_initialized(1);
}

=head2 nextAlignmemt

 Args        : none

 Example     : $alignment = $Blastz->nextAligment

 Descritpion : return the next HSP-like alignment

 Returntype  : array of Bio::EnsEMBL::DnaDnaAlignFeature

 Exceptions  : none

 Caller      : general

=cut

sub nextAlignment {
  my ($self) = @_;

  return undef if ($self->eof);

  my $fh = $self->fh;

  while (defined (my $line = <$fh>)) {
    next if ($line =~ /^\#:lav$/);
    if ($line =~ /^\#:eof$/) {
      close $self->fh;
      $self->eof(1);
      return undef;
    }

# s stanza : get there sequence length and strand
    if ($line =~ /^s\s+\{$/) {
      # on query
      $line = <$fh>;
      if ($line =~ /^\s*\"\S+\"\s+(\d+)\s+(\d+)\s+(\d+)\s+\d+$/) {
	my($start,$end,$strand) = ($1,$2,$3);
	$self->length($end-$start+1);
	$self->strand(1) if ($strand == 0);
	$self->strand(-1) if ($strand == 1);
      }
      # on database
      $line = <$fh>;
      if ($line =~ /^\s*\"\S+\"\s+(\d+)\s+(\d+)\s+(\d+)\s+\d+$/) {
	my($hstart,$hend,$hstrand) = ($1,$2,$3);
	$self->hlength($hend-$hstart+1);
	$self->hstrand(1) if ($hstrand == 0);
	$self->hstrand(-1) if ($hstrand == 1);
      }

      <$fh>; # skip } line
      next;
    }

# h stanza : get there seqname and hseqname
    if ($line =~ /^h\s+\{$/) {

      # on query
      $line = <$fh>;
      if ($line =~ /^\s+\">(\S+)\s*.*\"$/) {
	my $seqname = $1;

	$self->seqname($seqname);
      }

      # on database
      $line = <$fh>;
      if ($line =~ /^\s+\">(\S+)\s*.*\"$/) {
	my $hseqname = $1;
	$self->hseqname($hseqname);
      }

      <$fh>; # skip } line
      next;
    }
    
# a stanza : get there a alignment, with score, percent_id and positions
    if ($line =~ /^a\s+\{$/) {
      my ($score,$sum_match_bases,$sum_block_length) ;

      my @feature_pairs;

      while (defined ($line = <$fh>)) {
	last if ($line =~ /^\}$/);
	if ($line =~ /^\s+s\s+(\d+)$/) {
	  $score = $1;
	  next;
	}
	next if ($line =~ /^\s+b\s+\d+\s+\d+$/);
	next if ($line =~ /^\s+e\s+\d+\s+\d+$/);
	if ($line =~ /^\s+l\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)$/) {
	   my ($start,$hstart,$end,$hend,$percid) = ($1,$2,$3,$4,$5);

	   my $block_length = $end - $start + 1;
	   $sum_match_bases += $percid*$block_length/100;
	   $sum_block_length += $block_length;

	   if ($self->strand == -1) {
	     $start = $self->length - $end + 1;
	     $end = $start + $block_length - 1;
	   }
	   if ($self->hstrand == -1) {
	     $hstart = $self->hlength - $hend + 1;
	     $hend = $hstart + $block_length - 1;
	   }

	   my $feature_pair = new Bio::EnsEMBL::FeaturePair;
	   $feature_pair->seqname($self->seqname);
	   $feature_pair->start($start);
	   $feature_pair->end($end);
	   $feature_pair->strand($self->strand);
	   $feature_pair->hseqname($self->hseqname);
	   $feature_pair->hstart($hstart);
	   $feature_pair->hend($hend);
	   $feature_pair->hstrand($self->hstrand);
	   $feature_pair->score($score);
	   push @feature_pairs,$feature_pair;
	}
      }

      # calculating the average of percentage identity over the whole HSP-like
      # not including indels,it probably should...
      my $average_pecent_id = int($sum_match_bases/$sum_block_length*100);

      foreach my $feature_pair (@feature_pairs) {
	$feature_pair->percent_id($average_pecent_id);
      }

      my $alignment = new Bio::EnsEMBL::DnaDnaAlignFeature(-features => \@feature_pairs);
      return $alignment;
    }
  }
}

=head2 fh

 Arg [1]     : filehandle $filehandle (optional)

 Example     : $Blastz->fh($filehandle)

 Descritpion : get/set the filehandle value and return it

 Returntype  : filehandle

 Exceptions  : thrown if $filehandle is not the GLOB reference

 Caller      : general

=cut

sub fh {
  my ($self,$value) = @_;
  if(defined $value) {
    if (ref($value) eq "GLOB") {
      $self->{'_fh'} = $value;
    } else {
      $self->throw("value for fh method should be a filehandle\n");
    }
  }
  return $self->{'_fh'};
}

=head2 file

 Arg [1]     : string $filename_path (optional)

 Example     : $Blastz->file($filename_path)

 Descritpion : get/set the filename_path value and return it

 Returntype  : string

 Exceptions  : thrown if $filename_path is not found

 Caller      : general

=cut

sub file {
  my ($self,$value) = @_;
  if(defined $value) {
    if (-e $value) {
      $self->{'_file'} = $value;
    } else {
      $self->throw("file $value not found\n");
    }
  }
  return $self->{'_file'};
}

sub eof {
  my ($self,$value) = @_;
  if(defined $value) {
    $self->{'_eof'} = $value;
  }
  return $self->{'_eof'};
}

=head2 file

 Arg [1]     : string $commandline (optional)
               command line used to obtain the blastz output which is parsed
  
 Example     : $Blastz->commandline($commandline)

 Descritpion : get/set the commandline value and return it

 Returntype  : string

 Exceptions  : none

 Caller      : general

=cut

sub command_line {
  my ($self,$value) = @_;
  if (defined $value) {
    $self->{'_command_line'} = $value;
  }
  return $self->{'_command_line'};
}

=head2 matrix

 Arg [1]     : string $matrix (optional)
               matrix used to obtain the blastz output which is parsed
  
 Example     : $Blastz->matrix($matrix)

 Descritpion : get/set the matrix value and return it

 Returntype  : string

 Exceptions  : none

 Caller      : general

=cut

sub matrix {
  my ($self,$value) = @_;
  if (defined $value) {
    $self->{'_matrix'} = $value;
  }
  return $self->{'_matrix'};
}

=head2 options

 Arg [1]     : string $options (optional)
               options used to obtain the blastz output which is parsed
  
 Example     : $Blastz->options($options)

 Descritpion : get/set the options value and return it

 Returntype  : string

 Exceptions  : none

 Caller      : general

=cut

sub options {
  my ($self,$value) = @_;
  if (defined $value) {
    $self->{'_options'} = $value;
  }
  return $self->{'_options'};
}

sub _parsing_initialized {
  my ($self,$value) = @_;
  if (defined $value) {
    $self->{'_parsing_initialized'} = $value;
  }
  return $self->{'_parsing_initialized'};
  
}

=head2 seqname

 Arg [1]     : string $seqname (optional)
               name of the query sequence
  
 Example     : $Blastz->seqname($seqname)

 Descritpion : get/set the seqname value and return it

 Returntype  : string

 Exceptions  : none

 Caller      : general

=cut

sub seqname {
  my ($self,$value) = @_;
  if (defined $value) {
    $self->{'_seqname'} = $value;
  }
  return $self->{'_seqname'};
}

=head2 hseqname

 Arg [1]     : string $hseqname (optional)
               name of the database sequence
  
 Example     : $Blastz->hseqname($hseqname)

 Descritpion : get/set the hseqname value and return it

 Returntype  : string

 Exceptions  : none

 Caller      : general

=cut

sub hseqname {
  my ($self,$value) = @_;
  if (defined $value) {
    $self->{'_hseqname'} = $value;
  }
  return $self->{'_hseqname'};
}

=head2 length

 Arg [1]     : int $length (optional)
               sequence length of the query sequence
  
 Example     : $Blastz->length($length)

 Descritpion : get/set the length value and return it

 Returntype  : int

 Exceptions  : none

 Caller      : general

=cut

sub length {
  my ($self,$value) = @_;
  if (defined $value) {
    $self->{'_length'} = $value;
  }
  return $self->{'_length'};
}

=head2 hlength

 Arg [1]     : int $length (optional)
               sequence length of the database sequence
  
 Example     : $Blastz->hlength($length)

 Descritpion : get/set the hlength value and return it

 Returntype  : int

 Exceptions  : none

 Caller      : general

=cut

sub hlength {
  my ($self,$value) = @_;
  if (defined $value) {
    $self->{'_hlength'} = $value;
  }
  return $self->{'_hlength'};
}

=head2 strand

 Arg [1]     : int $strand (optional)
               strand of the query sequence
  
 Example     : $Blastz->strand($strand)

 Descritpion : get/set the strand value and return it

 Returntype  : int

 Exceptions  : none

 Caller      : general

=cut

sub strand {
  my ($self,$value) = @_;
  if (defined $value) {
    $self->{'_strand'} = $value;
  }
  return $self->{'_strand'};
}

=head2 hstrand

 Arg [1]     : int $strand (optional)
               strand of the query sequence
  
 Example     : $Blastz->hstrand($strand)

 Descritpion : get/set the hstrand value and return it

 Returntype  : int

 Exceptions  : none

 Caller      : general

=cut

sub hstrand {
  my ($self,$value) = @_;
  if (defined $value) {
    $self->{'_hstrand'} = $value;
  }
  return $self->{'_hstrand'};
}

1;
