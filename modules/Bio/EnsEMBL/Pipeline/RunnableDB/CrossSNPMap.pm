
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

use Bio::Root::RootI;
use Bio::EnsEMBL::Pipeline::RunnableDBI;
use Bio::EnsEMBL::Pipeline::Runnable::CrossMatch;
#use Data::Dumper;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDBI Bio::Root::RootI);

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
  if ($first) {
    print "this is the first round through--checking wether we have a hit in current emsembl contig\n";
  }
  else {
    print "this is the second round through--finding a hit from current emsembl contig\n";
  }
  $self->snp_dbobj($snpdb);
  $self->dbobj($db);
  $self->_score($score);
  $self->_minmatch($minmatch);
  $self->_masklevel($masklevel);
  $self->_debug($debug);
  $self->_first($first);
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

sub fetch_input{
  
  my ($self,@args) = @_;
  my ($inseq,$start_refnum,$end_refnum,$contigs,$mouse) = $self->_rearrange([qw(INSEQ 
										START_REFNUM 
										END_REFNUM 
										CONTIGS 
										MOUSE)],@args);
  my @contigs=@$contigs if $contigs;
  print "inseq is $inseq start_num is $start_refnum and end is $end_refnum and contigs is @contigs and mouse is $mouse\n";

  my (@infos,$inseq_name,$inseq_seq,%id_hash,@newseq,$clone_seqobj);
  
  if ($start_refnum =~ /^\d+$/) {
    if ($self->_first) {
      @infos = $self->snp_dbobj->get_snp_info_between_two_internalids($start_refnum,$end_refnum,$mouse);
    }
    else {
      @infos = $self->snp_dbobj->get_snp_info_by_refsnpid($start_refnum,$mouse);
    }
    my (%all_version, %same_version, @diff_version, $no_clone, $no_clone_count); 
    my $infos_count = @infos;
    print "infos_count is $infos_count\n";
    INFO : foreach my $info (@infos) {
      if ($self->_first and !$same_version{$info->snpid} or !$self->_first) {
	$all_version{$info->snpid}=$info;
	if ($info->score >2) { ##get rid off mapweight >2##
	  if ($self->_debug>1) {
	    print STDERR "refsnpid ",$info->snpid, "has  mapweight >2\n";
	  }
	  next INFO;
	}
	if ($self->_debug>1) {
	  #print Dumper($info);
	  print STDERR "refsnpid is ",$info->snpid, " acc and version is ",$info->seqname, " strand is ",$info->strand,"\n";
	}
	$self->_acc($info->acc);
	$no_clone = "";
	my ($newclone, $new_version);
	eval {
	  $newclone = $self->dbobj->get_Clone($info->acc);
	};
	if ($newclone and !@contigs) { ##mainly for first round check to see clone and version are in ensembl
	  print STDERR "here is id ",$newclone->id," and id2 is ",$info->acc," and version ", $newclone->embl_version, " ver2 is ",$info->version, "\n";
	  if ($newclone->id == $info->acc and $newclone->embl_version == $info->version) {
	    print "same acc and version for ",$info->snpid, "\n";
	    $same_version{$info->snpid}=$info;
	    next INFO; 
	  }
	  elsif ($self->_first) { ##if it is the first round, just go through##
	    next INFO; 
	  }
	  else { ##mainly for we have same clone, diff version
	    $self->_store_snp_seqobj($info);
	    $new_version = $newclone->embl_version;
	    $self->_version($new_version);
	    
	    @contigs = $newclone->get_all_Contigs();
	  }
	}
	elsif (@contigs) { ###for having $refsnpid and @contigs (for normal mouse mapping)
	  $self->_store_snp_seqobj($info);
	}
	else {
	  print ++$no_clone_count," no clone ",$info->acc," or ssaha hits (contigs) in database\n";
	  $no_clone=$info;
	}
      }
    }
    $self->all_version(\%all_version);
    $self->same_version(\%same_version);
    if ($no_clone) {
      return ($no_clone);
    }
  } ##end of if start_num =~ /^\d+$/
  elsif ($inseq and @contigs) { ###write for hgbase snp or no refsnpid###
    ($inseq_name,$inseq_seq) = split /\_/, $inseq;
    
    my $inseq_length = length($inseq_seq);
    my $num_bases = ($inseq_seq =~ tr/[WACTG-\_]//);
    my ($fseq) = split /[WACTG-\_]+/,$inseq_seq;
    print "THIS FSEQ $fseq and num_bases is $num_bases\n";
    my $snp_pos = length($fseq)+1;
    $self->_snp_pos($snp_pos);
    $self->_num_bases($num_bases);
    my $snp_seqobj = Bio::PrimarySeq->new( -display_id =>$inseq_name , -seq => $inseq_seq);
    $self->_snp_seqobj($snp_seqobj);
    $self->_snp_seq_length($inseq_length);
    
  }
  
  foreach my $contig(@contigs ) {
    my $version;
    $contig->id =~ /\w+\.(\d+)\.\d+\..*/;
    $version = $1;
    if (!$version) {
      $version="0";
    }
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
    my $score; my $cross;
    #if ($self->_snp_seq_length) {
    #  $score = $self->_snp_seq_length-int($self->_snp_seq_length*0.05);###equvalent to 95%
    #}
    #else {
      $score = $self->_score;
    #}
    #print STDERR "The minscore is $score\n"; 
    print "newseq is $newseq and snp_seqobj is ",$self->_snp_seqobj,"\n";
    if ($newseq and $self->_snp_seqobj) {
      my $cross = Bio::EnsEMBL::Pipeline::Runnable::CrossMatch->new( 
								    -nocopy => 1,
								    -seq1 => $self->_snp_seqobj,
								    -seq2 => $newseq,
								    -score => $score,
								    -minmatch => $self->_minmatch,
								    -masklevel => $self->_masklevel,
								   );
      
      push(@{$self->{'_crossmatch'}},$cross);
      
      $self->idhash(\%id_hash);
    }
  }
}

sub _store_snp_seqobj {

  my ($obj,$info) = @_;
  my ($snp_seq,$snp_seq_length,$sub_snp_seq,$snp_pos,@bases,$clone_seq,$snp_clone_seq);
  my $observed = $info->alleles;
  $observed =~ s/\/.*//;
  my $seq5 = $info->upStreamSeq;
  my $seq3 = $info->dnStreamSeq;
  #print STDERR "in _store_snp_seqobj seq5 is $seq5 and seq3 is $seq3\n";
  if (!$seq5 and !$seq3) {
    return ("no flank seqs");
  }
  else {
    $snp_pos = (length $seq5) +1;
    $snp_seq = $seq5.$observed.$seq3;
    $snp_seq =~ s/\s+|-//g;
    $snp_seq_length=length($snp_seq);
  }
  @bases = split "", $observed ;
  my $num_bases = @bases;
  $obj->_snp_pos($snp_pos);
  $obj->_num_bases($num_bases);
  $obj->_snp_seq_length($snp_seq_length);

  my $snp_seqobj = Bio::PrimarySeq->new( -display_id => $info->snpid, -seq => $snp_seq);
  $obj->_snp_seqobj($snp_seqobj);
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
   
  my (@fp,@final);
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
  
  print STDERR "Analysing output for clone ".$self->_acc."\n";
  
 FP: foreach my $fp (@sorted) {
    if ($self->_debug>1) {
      print "Got feature pair between contig ".$fp->seqname." (".$fp->start."-".$fp->end.") and contig ".$fp->hseqname." (".$fp->hstart."-".$fp->hend.") with score ".$fp->score."\n";
    }
    my $ref = $self->idhash;
    my %id_hash = %$ref;
    my $backid = $id_hash{$fp->hseqname};
    #my $sn=$self->_acc.".".$backid;
    my $sn=$fp->seqname;
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
  if (!@fp) {return;}

  my (%final_lists, @snp_fps, $final_line, $hit_line);

  ##To find matches that cover snp_pos
  @snp_fps=$self->find_fps_cover_snp("0",@fp);
  
  ##if no match cover snp_pos, try take 1 more bp on either side of the matches
  if (@snp_fps==0) {
    @snp_fps=$self->find_fps_cover_snp("1",@fp);
  }
  if (@snp_fps==0 or @snp_fps >2 ) {
    if ($self->_debug>1) {
	print "snp_pos ",$self->_snp_pos," is not in matching range or more than 3 hits with same score\n";
      }
    return;
  }
  else {
   foreach my $fp (@snp_fps) {
     
     print "score for snp_fps lists ",$fp->score,"\n";
     
     if ($fp->hseqname =~ /\w+\.(\d+)\.\d+\..*/) {
       $self->_version($1);
     }
     else {
       $self->_version ("0");
     }
     
     $final_line = $fp->seqname."\t".$fp->hseqname."\t".$self->_version."\t".$fp->start."\t".$fp->end."\t\tsnp_map\t".$fp->strand,"\t".$fp->score;
     
     if (@snp_fps==2) {
       print "last two hits for snp_pos have same scores\n";
       print "TWO_HITS $final_line\n";
       $hit_line .="TWO_HITS $final_line\n";
      }
     else {
       $hit_line = $final_line;
     }
   }
   chomp($hit_line);
   return ($hit_line);####return if having 1 or 2 hits#####
 }
}


sub find_fps_cover_snp{

  my ($self, $num, @fp) = @_;

  my @snp_fps;
  
  foreach my $fp ( @fp ) {
    
    my $snp_pos = $self->_snp_pos;

    if ($self->_debug>1) {
      print "In write_out ",$fp->seqname," ",$fp->hseqname," ",$fp->score," snp_pos is ",$self->_snp_pos," ",$fp->start," ",$fp->end," ",$fp->hstart," ",$fp->hend,"\n";
    }
    
    if (($snp_pos>=($fp->start-$num) and $snp_pos<=($fp->end+$num) or $snp_pos>=($fp->end-$num) and $snp_pos<=($fp->start+$num)) and $fp->start != $fp->end) {
      $self->record_snp($fp);
      push (@snp_fps, $fp);
    }
  }
  return (@snp_fps);
}

sub record_snp{

  my ($self, $fp) = @_;
  
  my $snp_pos = $self->_snp_pos;
  my ($new_strand, $new_start);
  my $num_bases = $self->_num_bases;
  if ($fp->strand ==-1 or $fp->hstrand ==-1) {
    $new_strand=-1;
    $new_start = $fp->hstart+($fp->end-($snp_pos+$num_bases-1));
  }
  else {
    $new_strand=1;
    $new_start = $fp->hstart+($snp_pos-$fp->start);
  }
  my $new_end = $new_start+$num_bases-1;
  
  $fp->start($new_start);
  $fp->end($new_end);
  $fp->strand($new_strand);
  
  if ($self->_debug>1) {
    print $fp->seqname,"\t",$fp->hseqname,"\t",$fp->score,"\t",$self->_version,"\t$new_start\t$new_end\t\tsnp_map\t$new_strand\n";
  }
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
      $self->{'_minmatch'} = $value;
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

sub _snp_seq_length{
  my ($self,$value) = @_;
  if( defined $value) {
    $self->{'_snp_seq_length'} = $value;
  }
  return $self->{'_snp_seq_length'};
}

sub _snp_seqobj{
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'_snp_seqobj'} = $value;
    }
    return $self->{'_snp_seqobj'};
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
