
#
# Ensembl module for Bio::EnsEMBL::Pipeline::RunnableDB::CrossCloneMap.pm
#
# Cared for by Ewan Birney <birney@ebi.ac.uk>
#
# Copyright GRL and EBI
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Pipeline::RunnableDB::CrossSNPMap.pm - DESCRIPTION of Object

=head1 SYNOPSIS

Give standard usage here

=head1 DESCRIPTION

Describe the object here

  my $crossmap = Bio::EnsEMBL::Pipeline::RunnableDB::CrossSNPMap->new(-crossdb=>$snpdb,
								      -db      =>$db,
								      -score =>$score,
								      -minmatch =>$minmatch,
								      -masklevel =>$masklevel,
								      -debug => $debug,
								      -first =>$first
								     );
$crossmap->fetch_input;
$crossmap->run;
$crossmap->write_output;

=head1 CONTACT

Ensembl - ensembl-dev@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::Pipeline::RunnableDB::CrossSNPMap;
use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Root;
use Bio::EnsEMBL::Pipeline::RunnableDBI;
use Bio::EnsEMBL::Pipeline::Runnable::CrossMatch;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDBI Bio::EnsEMBL::Root);

sub new {
  my($class,@args) = @_;

  my $self = {};
 
  $self->{'_final'} = [];
  $self->{'_crossmatch'} = [];
  bless $self,$class;
  $self->id_generation('start');
  my ($snpdb,$db,$score,$minmatch,$masklevel,$debug,$first) = $self->_rearrange([qw(SNPDB DB SCORE MINMATCH MASKLEVEL DEBUG FIRST )],@args);

  if ((!$snpdb) || (!$snpdb->isa('Bio::EnsEMBL::ExternalData::SNPSQL::DBAdaptor'))) {
      $self->throw("You need to provide a snp database adaptor to run CrossSNPMap!");
  }
  print "this is first $first\n";
  $self->snp_dbobj($snpdb);
  $self->dbobj($db);
  $self->_score($score);
  $self->_minmatch($minmatch);
  $self->_masklevel($masklevel);
  $self->_debug($debug);
  $self->_first($first);
  $self->_minmatch($minmatch);
  $self->_masklevel($masklevel);
  return $self;

}



=head2 id_generation;

 Title   : id_generation; 
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub id_generation{
   my ($self,$start) = @_;

   if ($start) {
       $self->{'id'}=0; 
   }
   my $c = $self->{'id'};
   $self->{'id'}=$c+1;
   return $c;
}

=head2 check_clone_in_GP

 Title   : check_clone_in_GP
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub fetch_input{
  
  my ($self,$start_refnum,$end_refnum,@contigs) = @_;
  
  my @infos = $self->snp_dbobj->get_snp_info_between_two_refsnpids($start_refnum,$end_refnum);
  my (%all_version, %same_version, @diff_version, $no_clone);
  INFO : foreach my $info (@infos) {
    $all_version{$info->snpid}=$info;
    if ($info->score >2) {
      next INFO;
    }
    if ($self->_debug>1) {
      print STDERR "refsnpid is ",$info->snpid, " acc and version is ",$info->seqname, " strand is ",$info->strand,"\n";
    }
    $self->_acc($info->acc);
    $no_clone = "";
    my ($newclone, $new_version);
    eval {
      $newclone = $self->dbobj->get_Clone($info->acc);
    };
    if ($newclone) {
      print STDERR "here is id ",$newclone->id," and id2 is ",$info->acc," and version ", $newclone->embl_version, " ver2 is ",$info->version, "\n";
      if ($newclone->id == $info->acc and $newclone->embl_version == $info->version) {
	$same_version{$info->snpid}=$info;
      }
      elsif (@contigs or $newclone->id == $info->acc and $newclone->embl_version != $info->version) {
	#if it is the first attempt, goto next record to see wether we can pick up same clone with same version. If it is not first attempt, we will map it
	if ($self->_first) {
	  next INFO; 
	}
	
	my ($snp_seq,$sub_snp_seq,$snp_pos,@bases,$clone_seq,$snp_clone_seq);
	my $observed = $info->alleles;
	$observed =~ s/\/.*//;
	#print STDERR "seq5 is $seq5 and seq3 is $seq3\n";
	my $seq5 = $info->upStreamSeq;
	my $seq3 = $info->dnStreamSeq;
	if (!$seq5 or !$seq3) {
	  return ("no flank seqs");
	}
	else {
	  $snp_pos = (length $seq5) +1;
	  $snp_seq = $seq5.$observed.$seq3;
	  $snp_seq =~ s/\s+//g;
	}
	@bases = split "", $observed ;
	my $num_bases = @bases;
	$self->_snp_pos($snp_pos);
	$self->_num_bases($num_bases);
	
	my $snp_seqobj = Bio::PrimarySeq->new( -display_id => $info->snpid, -seq => $snp_seq);
	
	if ($newclone) {$new_version = $newclone->embl_version;}
	$self->_version($new_version);
	if (!@contigs) {
	  @contigs = $newclone->get_all_Contigs();
	}
	
	my @newseq;my %id_hash;my $clone_seqobj;
	foreach my $contig(@contigs ) {
	  $contig->id =~ /\w+\.(\d+)\.\d+\..*/;
	  my $version = $1;
	  my $id = $self->id_generation;
	  $self->_version($version);
	  $id_hash{$id}=$contig->id;
	  
	  my $seq = Bio::PrimarySeq->new( -display_id => "$id", -seq => $contig->seq);
	  if ($self->_debug>1) {
	    print STDERR "Created seq with id $id for ".$id_hash{$id}." contig id ".$contig->id."\n";
	  }
	  push(@newseq,$seq);
	}
	
	if ($self->_debug >1) {
	  print STDERR "Creating crossmatches for clone ".$self->_acc."\n";
	}
	foreach my $newseq ( @newseq ) {
	  my $cross = Bio::EnsEMBL::Pipeline::Runnable::CrossMatch->new( 
									-nocopy => 1,
									-seq1 => $snp_seqobj,
									-seq2 => $newseq,
									-score => $self->_score(),
									-minmatch => $self->_minmatch,
									-masklevel => $self->_masklevel,
								       );
	  
	  push(@{$self->{'_crossmatch'}},$cross);
	  
	  $self->idhash(\%id_hash);
	}
      }
    }
    else {
      print "no clone ",$info->acc," in database\n";
      $no_clone=$info;
    }
    $self->all_version(\%all_version);
    $self->same_version(\%same_version);
  }
  if ($no_clone) {
    return ($no_clone);
  }
}

=head2 idhash
  
  Title   : idhash
  Usage   : $obj->idhash($newval)
 Function: Getset for idhash value
  Returns : value of idhash
  Args    : newvalue (optional)
  
  
=cut
  
sub idhash{
  my $obj = shift;
  if( @_ ) {
    my $value = shift;
    $obj->{'idhash'} = $value;
  }
  return $obj->{'idhash'};
  
}

sub same_version {
  my $obj = shift;
  if( @_ ) {
    my $value = shift;
    $obj->{'_same_version'} = $value;
  }
  return $obj->{'_same_version'};
  
}
sub all_version {
  my $obj = shift;
  if( @_ ) {
    my $value = shift;
    $obj->{'_all_version'} = $value;
  }
  return $obj->{'_all_version'};
  
}

sub diff_version {
  my $obj = shift;
  if( @_ ) {
    my $value = shift;
    $obj->{'_diff_version'} = $value;
  }
  return $obj->{'_diff_version'};
  
}


=head2 run

 Title   : run
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub run{
  my ($self) = @_;
   
  #The feature pair array will represent the matrix of hits between inputs
  my @fp;
  #Run all the crossmatch runnables created in fetch_input
  if ($self->_debug>1) {
    print STDERR "Running crossmatches for clone ".$self->_acc."\n";
  }
  foreach my $crossmatch (@{$self->{'_crossmatch'}}) {
    $crossmatch->run();
    #Push all the feature pairs into the array
    push (@fp,$crossmatch->output);
  }
  
  #Sort the array by score
  my @sorted= sort { $b->score <=> $a->score} @fp;
  
  my $initial=1;
  my %looks_ok=();
  
  print STDERR "Analysing output for clone ".$self->_acc."\n";
  #The juicy part: look at the matrix, and sort out the mapping
 FP: foreach my $fp (@sorted) {
    if ($self->_debug>1) {
      print "Got feature pair between contig ".$fp->seqname." (".$fp->start."-".$fp->end.") and contig ".$fp->hseqname." (".$fp->hstart."-".$fp->hend.") with score ".$fp->score."\n";
    }
    #Take the first one as correct, because it hsa the highest score...
    if ($initial) {
      #print STDERR "Pushing it to final (first fp)\n";
      push (@{$looks_ok{$fp->seqname}},$fp);
      $initial=undef;
      next FP;
    }
       
    #Check if this feature pair is consistent with the rest
    
    #If seqname already in final map, check other matches
    
    if ($looks_ok{$fp->seqname}) {
      my @ha_match=@{$looks_ok{$fp->seqname}};
      #print STDERR "Contig already in final map, checking...\n";
      #sort other matches by start
      my @s_matches= sort { $a->start <=> $b->start} @ha_match;
      my $first=$s_matches[0];
      
      #Speed up by eliminating the two most obvious cases...
      #If the match is before the first match on this contig, add to final
      if ($fp->end < $first->start){
	#print STDERR "Pushing it to final (before first)\n";
	push (@{$looks_ok{$fp->seqname}},$fp);
	next FP;
      }
      #If the match is after the last match on this contig, add to final 
      my $size=scalar(@s_matches);
      if ($fp->start > $s_matches[$size-1]->end) {
	#print STDERR "Pushing it to final (after last)\n";
	push (@{$looks_ok{$fp->seqname}},$fp);
	next FP;
      } 
      
      #Loop through all the other matches and check if $fp does not
      #overlap any of them
      my $add=1;
      foreach my $match (@s_matches) {
	#If fp start or end are contained in any match, overlapping...
	if ((($fp->start > $match->start) && ($fp->start < $match->end)) || (($fp->end < $match->end) && ($fp->end > $match->start))) {
	  $add=0;
	}
      }
      if ($add) {
	#print STDERR "Pushing it to final (no overlaps)\n";
	push (@{$looks_ok{$fp->seqname}},$fp);
	next FP;
      } 
      #In any other case, do not add this feature pair
      #print STDERR "End of checks for this contig, not added to final\n";
    }
    
    #Else, just add to final map
    else {
      #print STDERR "This seqname was not found earlier, add to final map\n";
      push (@{$looks_ok{$fp->seqname}},$fp);
    }
  }
  my @final=();
  foreach my $seqname (keys %looks_ok) {
    foreach my $fp (@{$looks_ok{$seqname}}) {
      #print "seqname is ",$fp->seqname,"\n";
      my $ref = $self->idhash;
      my %id_hash = %$ref;
      my $backid = $id_hash{$fp->hseqname};
      #my $sn=$self->_acc.".".$backid;
      my $sn=$fp->seqname;
      $backid = $id_hash{$fp->hseqname}; 
      my $hsn;
      if ($backid) {
	$hsn=$backid;}
      else {$hsn=$self->_acc}
      
      $sn =~ s/\_/\./g;
      $hsn =~ s/\_/\./g;
      $fp->seqname($sn);
      $fp->hseqname($hsn);
      push (@final,$fp);
    }
  }
  %looks_ok=();
  
  #foreach my $fp (@final) {
  #print STDERR "In final, got ".$fp->seqname."-".$fp->hseqname."\n";
  
  #}
  push(@{$self->{'_final'}},@final);
}

	
sub sort_numeric {$a<=>$b}
=head2 write_output
    
    Title   : write_output
    Usage   :
  Function:
    Example :
    Returns : 
    Args    :
    
    
=cut
    
sub write_output{
  my ($self) = @_;
  
  my @fp=@{$self->{'_final'}};
  
  my %final_lists;
  
  foreach my $fp ( @fp ) {
    my $seqname = $fp->seqname;
    my $hseqname = $fp->hseqname;
    my $score = $fp->score;
    my $version = $self->_version;
    
    if ($self->_debug>1) {
      print "In write_out ",$fp->seqname," ",$fp->hseqname," ",$fp->score," ",$self->_snp_pos," ",$fp->start," ",$fp->end," ",$fp->hstart," ",$fp->hend,"\n";
    }
    my ($strand1,$new_start);
    my $snp_pos = $self->_snp_pos;

    if ($snp_pos>=$fp->start and $snp_pos<=$fp->end or $snp_pos>=$fp->end and $snp_pos<=$fp->start) {
      
      my $num_bases = $self->_num_bases;
      if ($fp->strand ==-1 or $fp->hstrand ==-1) {
	$strand1=-1;
	$new_start = $fp->hend-($snp_pos-$fp->start);
      }
      else {
	$strand1=1;
	$new_start = $fp->hstart+($snp_pos-$fp->start);
      }
      my $new_end = $new_start+$num_bases-1;
      push (@{$final_lists{$snp_pos}}, "$seqname\t$hseqname\t$score\t$version\t$new_start\t$new_end\t\tsnp_map\t$strand1");
      if ($self->_debug>1) {
	print $fp->seqname,"\t",$fp->hseqname,"\t",$fp->score,"\t",$self->_version,"\t$new_start\t$new_end\tsnp_map\t$strand1\n";
      }
    }
    else {
      if ($self->_debug>1) {
	print "snp_pos $snp_pos is not in matching range ",$fp->start, " ", $fp->end, "\n";
      }
    }
  }
  my $final_line;
  foreach my $snp_pos (keys % final_lists) {
    my (@sorted_lists, @a_fields, @b_fields,@last_words);
    if (@{$final_lists{$snp_pos}}>1) {
      @sorted_lists = sort {
	@a_fields = split /\t/, $a; 
	@b_fields = split /\t/, $b; 
	$a_fields[2] <=> $b_fields[2];
      } @{$final_lists{$snp_pos}};
      #print "again @sorted_lists\n";
      my $last_line = pop @sorted_lists;
      my $last_line2 = pop @sorted_lists;
      my @last_line = split /\t/, $last_line;
      my @last_line2 = split /\t/, $last_line2;
      if ($last_line and $last_line[2] != $last_line2[2]) {
      	print "this is last line $last_line\n";
	@last_words = split /\t/, $last_line;
	if ($last_words[1] =~ /\w+\.(\d+)\.\d+\..*/) {
	  $last_words[3] = $1;
	}
	$final_line = join "\t", @last_words[0,1,3..8];
      }
      else {
	print "last two hits for $snp_pos have same scores\n";
      }
    }
    else {
      my $last_line = pop @{$final_lists{$snp_pos}};
      if ($last_line) {
	#print "this is last line $last_line\n";
	@last_words = split /\t/, $last_line;
	if ($last_words[1] =~ /\w+\.(\d+)\.\d+\..*/) {
	  $last_words[3] = $1;
	}
	$final_line = join "\t", @last_words[0,1,3..8];
      }
    }
  }
  return ($final_line);
}

=head2 dbobj

 Title   : dbobj
 Usage   : $obj->dbobj($newval)
 Function: 
 Example : 
 Returns : value of dbobj
 Args    : newvalue (optional)


=cut

sub dbobj{
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'_dbobj'} = $value;
    }
    return $self->{'_dbobj'};

}

sub snp_dbobj{
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'_snp_dbobj'} = $value;
    }
    return $self->{'_snp_dbobj'};

}

=head2 _score

 Title   : _score
 Usage   : $obj->_score($newval)
 Function: 
 Example : 
 Returns : value of score
 Args    : newvalue (optional)


=cut

sub _score{
  my ($self,$value) = @_;
  if( defined $value) {
    $self->{'_score'} = $value;
  }
  return $self->{'_score'};
  
}

=head2 _debug

 Title   : _debug
 Usage   : $obj->_debug($newval)
 Function: 
 Example : 
 Returns : value of _debug
 Args    : newvalue (optional)


=cut


sub _debug{
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'_debug'} = $value;
    }
    return $self->{'_debug'};

}

=head2 _minmatch

 Title   : _minmatch
 Usage   : $obj->_minmatch($newval)
 Function: 
 Example : 
 Returns : value of _minmatch
 Args    : newvalue (optional)


=cut


sub _minmatch{
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'_debug'} = $value;
    }
    return $self->{'_minmatch'};

}

=head2 _masklevel

 Title   : _masklevel
 Usage   : $obj->_debug($newval)
 Function: 
 Example : 
 Returns : value of _masklevel
 Args    : newvalue (optional)


=cut


sub _masklevel{
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'_masklevel'} = $value;
    }
    return $self->{'_masklevel'};

}


=head2 _first

 Title   : _first
 Usage   : $obj->_first($newval)
 Function: 
 Example : 
 Returns : value of _first
 Args    : newvalue (optional)


=cut

sub _first{
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'_first'} = $value;
    }
    return $self->{'_first'};
}
 

=head2 _acc

 Title   : _acc
 Usage   : $obj->_acc($newval)
 Function: 
 Example : 
 Returns : value of _acc
 Args    : newvalue (optional)


=cut

sub _acc{
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'_acc'} = $value;
    }
    return $self->{'_acc'};

}

=head2 _version

 Title   : _debug
 Usage   : $obj->_version($newval)
 Function: 
 Example : 
 Returns : value of _version
 Args    : newvalue (optional)


=cut

sub _version{
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'_version'} = $value;
    }
    return $self->{'_version'};

}

=head2 _snp_pos

 Title   : _snp_pos
 Usage   : $obj->_snp_pos($newval)
 Function: 
 Example : 
 Returns : value of _snp_pos
 Args    : newvalue (optional)


=cut

sub _snp_pos{
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'_snp_pos'} = $value;
    }
    return $self->{'_snp_pos'};
}

=head2 _num_bases

 Title   : _num_bases
 Usage   : $obj->_debug($newval)
 Function: 
 Example : 
 Returns : value of _num_bases
 Args    : newvalue (optional)


=cut

sub _num_bases{
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'_num_bases'} = $value;
    }
    return $self->{'_num_bases'};

}

1;
#$refsnpid,$snpclass,$snptype,$observed,$seq5,$seq3,$acc,$version,$start,$end,$strand
