#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::Runnable::AxtFilter

=head1 SYNOPSIS

    my $obj = Bio::EnsEMBL::Pipeline::Runnable::AxtFilter->new
    (
     -features        => $features,
     -best            => 1,
     -subset          => $matrix_file);
    
    $obj->run

    my @newfeatures = $obj->output;


=head1 DESCRIPTION

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

# Let the code begin...

    package Bio::EnsEMBL::Pipeline::Runnable::AxtFilter;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Pipeline::RunnableI;
use Bio::EnsEMBL::DnaDnaAlignFeature;


@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI);

sub new {
  my ($class,@args) = @_;

  my $self = $class->SUPER::new(@args);
  
  my ($features, 
      $best, 
      $subset_matrix,
      $subset_cutoff,
      $query,
      $target_lengths,
      $target_nib_dir,
      $fa2nib,
      $subset_axt,
      $axt_best,
      $lav2axt) = $self->_rearrange([qw(FEATURES
                                        BEST
                                        SUBSETMATRIX
                                        SUBSETCUTOFF
                                        QUERY
                                        TARGETLENGTHS
                                        TARGETNIB
                                        FA2NIB
                                        SUBSETAXT
                                        AXTBEST
                                        LAVTOAXT)],
                                    @args);


  $self->throw("You must supply a reference to an array of features with -features\n") 
      if not defined $features;
  $self->throw("You must supply a query sequence\n") 
      if not defined $query;
  $self->throw("You must supply a directory of nib files for the targets with -targetnib")
      if not defined $target_nib_dir;
  $self->throw("You must filter with either -best 1 or -subset matrix_file or both\n") 
      if not $best and not -e $subset_matrix;

  $subset_cutoff = 3400 if not defined $subset_cutoff;

  $self->faToNib($fa2nib) if defined $fa2nib;
  $self->subsetAxt($subset_axt) if defined $subset_axt;
  $self->axtBest($axt_best) if defined $axt_best;
  $self->lavToAxt($lav2axt) if defined $lav2axt;

  $self->best($best);
  $self->subset_matrix($subset_matrix);
  $self->subset_cutoff($subset_cutoff);
  $self->query($query);
  $self->features($features);
  $self->target_lengths($target_lengths);
  $self->target_nib_dir($target_nib_dir);
  
  return $self;
}





=head2 run

  Title   : run
  Usage   : $self->run()
  Function: 
  Returns : none
  Args    : 

=cut

sub run {
  my ($self) = @_;

  my $contig_name = $self->query->name;

  my $work_dir = $self->workdir . "/" . $contig_name . ".$$." . "axt_tmp_dir";
  my $seq_file = "$work_dir/$contig_name.fa";
  my $query_nib_dir = "$work_dir/$contig_name.nib";
  my $lav_file = "$work_dir/$contig_name.lav";
  my $axt_file = "$work_dir/$contig_name.axt";
  mkdir $work_dir;

  my $fh;

  ##############################
  # write features in lav format
  ############################## 
  open $fh, ">$lav_file" or 
      $self->throw("could not open lav file '$lav_file' for writing\n");
  $self->_write_lav($fh);
  close($fh);

  #################################
  # write the query in nob format 
  # for use by lavToAxt;
  #################################
  my $seqio = Bio::SeqIO->new(-format => 'fasta',
                              -file   => ">$seq_file");
  $seqio->write_seq($self->query);
  $seqio->close;
  mkdir $query_nib_dir;
  system($self->faToNib, $seq_file, "$query_nib_dir/$contig_name.nib") 
      and $self->throw("Could not convert fasta file $seq_file to nib");

  ##############################
  # convert the lav file to axt
  ##############################
  system($self->lavToAxt, $lav_file, $query_nib_dir, $self->target_nib_dir, $axt_file)
      and $self->throw("Could not convert $lav_file to Axt format");
   
  ##############################
  # Filter with axtBest
  ##############################
  if ($self->best) {
    my $tmp_axt_file = "$work_dir/$contig_name.$$.res.axt";
    system($self->axtBest, $axt_file, $contig_name, $tmp_axt_file)
        and $self->throw("Something went wrong with axtBest\n");
    rename $tmp_axt_file, $axt_file;
  }

  
  ##############################
  # Filter with subsetAxt
  ##############################
  if ($self->subset_matrix) {
    my $tmp_axt_file = "$work_dir/$contig_name.$$.res.axt";
    system($self->subsetAxt, $axt_file, $tmp_axt_file, $self->subset_matrix, $self->subset_cutoff)
        and $self->throw("Something went wrong with subsetAxt\n");
    rename $tmp_axt_file, $axt_file;  
  }

    
  #################################
  # Read features from the Axt file
  ##################################
  open $fh, $axt_file or $self->throw("Something went badly wrong; could not read result axt file\n");
  my @out = $self->_read_Axt($fh);
  close($fh);

  $self->output(@out);
  
  ##########
  # clean up
  ##########
  unlink $lav_file, $axt_file, $seq_file, "$query_nib_dir/$contig_name.nib";
  rmdir $query_nib_dir;
  rmdir $work_dir;

  return 1;
}




sub _write_lav {  
  my ($self, $fh) = @_;

  my (%features);  
  foreach my $feat (sort { $a->start <=> $b->start} @{$self->features}) {
    push @{$features{$feat->hseqname}{$feat->strand}{$feat->hstrand}}, $feat;
  }
  
  my $query_length = $self->query->length;
  my $query_name   = $self->query->name;
  
  foreach my $target (sort { my($a1)= $a =~ /^(\S+)\./; 
                             my($b1)= $b =~ /^(\S+)\./; 
                             $a1 cmp $b1 
                           } keys %features) {

    print $fh "#:lav\n";
    print $fh "d {\n   \"generated by Runnable/AxtFilter.pm\"\n}\n";

    foreach my $qstrand (keys %{$features{$target}}) {
      foreach my $tstrand (keys %{$features{$target}{$qstrand}}) {
        
        my $query_strand = ($qstrand == 1) ? 0 : 1;
        my $target_strand = ($tstrand == 1) ? 0 : 1;
        
        my $target_length;
        if ($self->target_lengths and exists($self->target_lengths->{$target})) {
          $target_length = $self->target_lengths->{$target};
        } elsif ($target =~ /^\S+\.(\d+)\-(\d+)/) {
          $target_length = $2 - $1 + 1;
        } else {
          $self->throw("Error writing lav entry for $target; bould not get sequence length\n");
        }

        print $fh "#:lav\n";
        print $fh "s {\n";
        print $fh "   \"$query_name\" 1 $query_length $query_strand 1\n";
        print $fh "   \"$target\" 1 $target_length $target_strand 1\n";
        print $fh "}\n";
        
        print $fh "h {\n";
        print $fh "   \">$query_name";
        if ($query_strand) {
          print $fh " (reverse complement)";
        }
        print $fh "\"\n   \">$target";
        if ($target_strand) {
          print $fh " (reverse complement)";
        }
        print $fh "\"\n}\n";
	
        foreach my $reg (@{$features{$target}{$qstrand}{$tstrand}}) {
          my $qstart = $query_strand ?  $query_length - $reg->end + 1 : $reg->start; 
          my $qend = $query_strand ?  $query_length - $reg->start + 1 : $reg->end; 
          my $tstart = $target_strand ? $target_length - $reg->hend + 1 : $reg->hstart; 
          my $tend = $target_strand ? $target_length - $reg->hstart + 1 : $reg->hend; 
          
          printf $fh "a {\n   s %d\n", $reg->score;
          print $fh "   b $qstart $tstart\n"; 
          print $fh "   e $qend $tend\n";
          
          foreach my $seg ($reg->ungapped_features) {
            my $qstartl = $query_strand ?  $query_length - $seg->end + 1 : $seg->start; 
            my $qendl = $query_strand ?  $query_length - $seg->start + 1 : $seg->end; 
            
            my $tstartl = $target_strand ? $target_length - $seg->hend + 1 : $seg->hstart; 
            my $tendl = $target_strand ? $target_length - $seg->hstart + 1 : $seg->hend; 
            
            printf $fh "   l $qstartl $tstartl $qendl $tendl %d\n", $seg->percent_id;
            
          }
          print $fh "}\n";
        }
        
        print $fh "x {\n   n 0\n}\n"; 
      }
    }
    print $fh "m {\n   n 0\n}\n#:eof\n";
  }
}


sub _read_Axt {
  my ($self, $fh) = @_;
    
  my ($query_name, %to_keep, @out_list);

  my $query_length = $self->query->length;
  
  while(<$fh>) {
    /^\d+\s+(\S+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\d+)/ and do {
      
      # query strand is always + in Axt format

      my ($cname, $cst, $cen, $hname, $hst, $hen, $cstr, $hstr, $score) = 
          ($1, $2, $3, $4, $5, $6, 1, $7, $8);
      
      if ($hstr eq "+") {
        $hstr = 1;
      } else {
        $hstr = -1;
      }
      
      if (not $query_name) {
        $query_name = $cname;
      } elsif ($cname ne $query_name) {
        $self->throw("read_Axt error: Two different contig names in same file [$query_name, $cname]\n");
      }
      
      my $target_length;
      if ($self->target_lengths and exists($self->target_lengths->{$hname})) {
        $target_length = $self->target_lengths->{$hname};
      } elsif ($hname =~ /^\S+\.(\d+)\-(\d+)/) {
        $target_length = $2 - $1 + 1;
      } else {
        $self->throw("Error reading Axt; cannot determine length of $hname\n");
      }

      my $seq1 = <$fh>; chomp $seq1; my @dna1 = split //, $seq1;
      my $seq2 = <$fh>; chomp $seq2; my @dna2 = split //, $seq2;
      
      my ($local_cst, $local_hst) = ($cst, $hst);
      
      my ($in_match, $total_matches, $total_aligned) = (0,0,0);
      my @regs;
      
      for (my $i = 0; $i < @dna1; $i++) {
        if ($dna1[$i] ne "-" and $dna2[$i] ne "-") {
          # match
          if (not $in_match) {
            push @regs, { 
              cst => $local_cst,
              cen => $local_cst,
              hst => $local_hst,
              hen => $local_hst,
              len => 1,
              match => ($dna1[$i] eq $dna2[$i]) ? 1 : 0,
            };
          }
          else {
            $regs[-1]->{cen}++;
            $regs[-1]->{hen}++;
            $regs[-1]->{len}++;
            
            if ($dna1[$i] eq $dna2[$i]) {
              $regs[-1]->{match}++;
            }
          }
          
          if ($dna1[$i] eq $dna2[$i]) {
            $total_matches++;
          }
          $total_aligned++;
          
          $in_match = 1;
        }
        elsif ($dna1[$i] eq "-" and $dna2[$i] ne "-") {
          if ($in_match) {
            $local_hst = $regs[-1]->{hen} + 2;
            $local_cst = $regs[-1]->{cen} + 1;
          }
          else {
            $local_hst++;
          }
          $in_match = 0;
        }				  
        elsif ($dna2[$i] eq "-" and $dna1[$i] ne "-") {
          if ($in_match) {
            $local_cst = $regs[-1]->{cen} + 2;
            $local_hst = $regs[-1]->{hen} + 1;
          }
          else {
            $local_cst++;
          }
          $in_match = 0;
        }
        # else pair of gaps; ignore;
      }

      # Now process the regions
      my $av_perc_id = ($total_matches / $total_aligned) * 100;
      
      my @subreg_list;
      foreach my $reg (@regs) {
        my ($local_cst, $local_cen, $local_hst, $local_hen) = 
            ($reg->{cst}, $reg->{cen}, $reg->{hst}, $reg->{hen});
		
        # need to flip the co-ords of the target match back		
        if ($cstr < 0) {
          my $cst_save = $local_cst;
          $local_cst = $query_length - $local_cen + 1;
          $local_cen = $query_length - $cst_save + 1;
        }
        if ($hstr < 0) {
          my $hst_save = $local_hst;
          $local_hst = $target_length - $local_hen + 1;
          $local_hen = $target_length - $hst_save + 1;
        }

        #create featurepair
        my $fp = new Bio::EnsEMBL::FeaturePair->new();
        $fp->seqname($cname);
        $fp->start($local_cst);
        $fp->end($local_cen);
        $fp->strand($cstr);
        $fp->hseqname($hname);
        $fp->hstart($local_hst);
        $fp->hend($local_hen);
        $fp->hstrand($hstr);
        $fp->score($score);
        $fp->percent_id($av_perc_id);
        
        push @subreg_list, $fp;
      }
	  
      # Occasionally, Axt tools return an alignment with a characteristic that the
      # Ensembl core code finds dodgy (e.g. a insert followed immediately by a delete).
      # Any such weirdness will cause an Exception in the below, so catch it and discard
      # the alignment
      my $df;
      eval {
        $df = new Bio::EnsEMBL::DnaDnaAlignFeature(-features => \@subreg_list);
      };
      push @out_list, $df unless $@;
    }
  }
  
  return @out_list;
}



=head2 query

    Title   :   query
    Usage   :   $self->query($primaryseq)
    Function:  
    Returns :   
    Args    :   

=cut

sub query {
  my ($self,$arg) = @_;
  
  if (defined($arg)) {
    $self->{'_query'} = $arg;
  }
  
  if (!exists($self->{'_query'})) {
    return undef;
  } else {
    return $self->{'_query'};
  }
}



=head2 target_lengths

    Title   :   target_lengths
    Usage   :   $self->target_length($hash_ref)
    Function:  
    Returns :   
    Args    :   

=cut

sub target_lengths {
  my ($self,$val) = @_;
  
  if (defined $val) {
    $self->{'_target_length_hash'} = $val;
  }
  
  return $self->{'_target_length_hash'};
}




=head2 target_nib_dir

    Title   :   target_nib_dir
    Usage   :   $self->query($primaryseq)
    Function:  
    Returns :   
    Args    :   

=cut

sub target_nib_dir {
  my ($self,$dir) = @_;
  
  if (defined $dir) {
    $self->{'_target_nib_dir'} = $dir;
  }
  
  return $self->{'_target_nib_dir'};
}



sub _query_nib_dir {
  my ($self,$dir) = @_;
  
  if (defined $dir) {
    $self->{'_query_nib_dir'} = $dir;
  }
  
  return $self->{'_query_nib_dir'};
}



=head2 features

    Title   :   features
    Usage   :   $self->features($features)
    Function:   ref. to the array fo features to be projected
    Returns :   
    Args    :   

=cut

sub features {
  my ($self,$arg) = @_;
  
  if (defined($arg)) {
    $self->{'_features'} = $arg;
  }
  
  if (!exists($self->{'_features'})) {
    $self->{'_features'} = [];
  }    

  return $self->{'_features'};
}



=head2 best

    Title   :   best
    Usage   :   $self->best(1)
    Function:   Binary flag that determines whether the given feature set is 
      filtered using axtBest (so that the query sequence is coverred at each 
      position by at most 1 target sequence)
    Returns :   
    Args    :   

=cut

sub best {
  my ($self,$arg) = @_;
  
  if (defined($arg)) {
    $self->{'_best'} = $arg;
  }
  
  if (!exists($self->{'_best'})) {
    $self->{'_best'} = 0;
  }    

  return $self->{'_best'};
}


=head2 subset_matrix

    Title   :   subset_matrix
    Usage   :   $self->subset_matrix($matrix_file)
    Function:   
      When this is set to a valid blastz matrix file, 
      the Runnable will rescore and trim the features with
      the matriz using subsetAxt
    Returns :   
    Args    :   

=cut

sub subset_matrix {
  my ($self,$arg) = @_;
  
  if (defined($arg)) {
    $self->{'_subset_matrix'} = $arg;
  }
  
  if (!exists($self->{'_subset_matrix'})) {
    $self->{'_subset_matrix'} = 0;
  }    

  return $self->{'_subset_matrix'};
}


=head2 subset_cutoff

    Title   :   subset_cutoff
    Usage   :   $self->subset_cutoff($score)
    Function:   
      The threshold to use with subsetAxt
    Returns :   
    Args    :   

=cut

sub subset_cutoff {
  my ($self,$arg) = @_;
  
  if (defined($arg)) {
    $self->{'_subset_cutoff'} = $arg;
  }
  
  if (!exists($self->{'_subset_cutoff'})) {
    $self->{'_subset_cutoff'} = 0;
  }    

  return $self->{'_subset_cutoff'};
}

=head2 faToNib

    Title   :   faToNib
    Usage   :   $self->faToNib("/usr/local/bin/faToNib");
    Function:   
    Returns :   
    Args    :   

=cut

sub faToNib {
  my ($self,$arg) = @_;
  
  if (defined($arg)) {
    $self->{'_faToNib'} = $arg;
  }
  
  if (!exists($self->{'_faToNib'})) {
    $self->{'_faToNib'} = "/usr/local/ensembl/bin/faToNib";
  }    

  return $self->{'_faToNib'};
}



=head2 subsetAxt

    Title   :   subsetAxt
    Usage   :   $self->subsetAxt("/usr/local/bin/subsetAxt");
    Function:   
    Returns :   
    Args    :   

=cut

sub subsetAxt {
  my ($self,$arg) = @_;
  
  if (defined($arg)) {
    $self->{'_subsetAxt'} = $arg;
  }
  
  if (!exists($self->{'_subsetAxt'})) {
    $self->{'_subsetAxt'} = "/usr/local/ensembl/bin/subsetAxt";
  }    

  return $self->{'_subsetAxt'};
}



=head2 axtBest

    Title   :   axtBest
    Usage   :   $self->axtBest("/usr/local/bin/axtBest");
    Function:   
    Returns :   
    Args    :   

=cut

sub axtBest {
  my ($self,$arg) = @_;
  
  if (defined($arg)) {
    $self->{'_axtBest'} = $arg;
  }
  
  if (!exists($self->{'_axtBest'})) {
    $self->{'_axtBest'} = "/usr/local/ensembl/bin/axtBest";
  }    

  return $self->{'_axtBest'};
}


=head2 lavToAxt

    Title   :   lavToAxt
    Usage   :   $self->lavToAxt("/usr/local/bin/lavToAxt");
    Function:   
    Returns :   
    Args    :   

=cut

sub lavToAxt {
  my ($self,$arg) = @_;
  
  if (defined($arg)) {
    $self->{'_lavToAxt'} = $arg;
  }
  
  if (!exists($self->{'_axtBest'})) {
    $self->{'_lavToAxt'} = "/usr/local/ensembl/bin/lavToAxt";
  }    

  return $self->{'_lavToAxt'};
}





=head2 output

    Title   :   merge
    Usage   :   $self->merge($output)
    Function:
    Returns :   
    Args    :   

=cut

sub output {
  my ($self,@feats) = @_; 
  
  if (@feats) {
    $self->{'_output'} = [@feats];
  }
  
  if (!exists($self->{'_output'})) {
    $self->{'_output'} = [];
  }    

  return @{$self->{'_output'}};
}



1;
