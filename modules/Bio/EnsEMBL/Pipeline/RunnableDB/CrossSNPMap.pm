
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
								     );
$crossmap->fetch_input;
$crossmap->seq_map;
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
use Bio::EnsEMBL::Pipeline::Runnable::CrossMatch;
#use Data::Dumper;

@ISA = qw(Bio::EnsEMBL::Root);

sub new {
  my($class,@args) = @_;

  my $self = {};
 
  $self->{'_final'} = [];
  $self->{'_crossmatch'} = [];
  bless $self,$class;
  $self->id_generation('start');
  $self->id_generation1('start');
  my ($snpdb,$db,$score,$minmatch,$masklevel,$debug) = $self->_rearrange([qw(SNPDB 
										  DB 
										  SCORE 
										  MINMATCH 
										  MASKLEVEL 
										  DEBUG )],@args);

#  if ((!$snpdb) || (!$snpdb->isa('Bio::EnsEMBL::ExternalData::SNPSQL::SNPAdaptor'))) {
#      $self->throw("You need to provide a snp database adaptor to run CrossSNPMap!");
#  }
  
  $self->snp_dbobj($snpdb);
  $self->dbobj($db);
  $self->_score($score);
  $self->_minmatch($minmatch);
  $self->_masklevel($masklevel);
  $self->_debug($debug);
  $self->_masklevel($masklevel);
  return $self;
}



=head2 id_generation;

 Title   : id_generation; 
 Usage   : my $id = $self->id_generation
 Function: generate ids by auto increment 
 Example : my $id = $self->id_generation
 Returns : id
 Args    : no


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

sub id_generation1{
   my ($self,$start) = @_;

   if ($start) {
       $self->{'id1'}=0; 
   }
   my $c = $self->{'id1'};
   $self->{'id1'}=$c+1;
   return $c;
}

=head2 check_data_for_acc_version

 Title   : check_data_for_acc_version
 Usage   : $crossmap->check_data_for_acc_version($start_intnum,$end_intnum);
           $start_intnum is internal_id in RefSNP table
 Function: check to see acc and version are same in snpdb and core contig
 Example : $crossmap->check_data_for_acc_version($start_intnum,$end_intnum);
 Returns : %same_version
 Args    : start and end internal_ids


=cut

sub check_data_for_acc_version {
  
  my ($self,$start_intnum,$end_intnum) = @_;
  
  my (@infos) = $self->snp_dbobj->get_snp_info_between_two_internalids($start_intnum,$end_intnum);
    
  my (%same_version); 
  my $infos_count = @infos;
  print "infos_count is $infos_count\n";
  INFO : foreach my $info (@infos) {
    if (!$same_version{$info->snpid}) {
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
            
      my ($newclone, $new_version);
      eval {
	$newclone = $self->dbobj->get_CloneAdaptor->fetch_by_name($info->acc);
      };
      if ($newclone) { 
	print STDERR "here is id ",$newclone->id," and id2 is ",$info->acc," and version ", $newclone->embl_version, " ver2 is ",$info->version, "\n";
	if ($newclone->id == $info->acc and $newclone->embl_version == $info->version) {
	  print "same acc and version for ",$info->snpid, "\n";
	  $same_version{$info->snpid}=$info;
	  next INFO; 
	}
      }
    }
  }
  return(\%same_version);
}

=head2 fetch_input

 Title   : fetch_input
 Usage   : $crossmap->fetch_input($refsnpid,$contigs,$mouse);
 Function: generate snp_seqobj and objecs for cross_match runs
 Example : see usage
 Returns : no
 Args    : refsnpid and matching contigs_objs
 
=cut

sub fetch_input{
  
  my ($self,$refsnpid,$contigs,$mouse) = @_;
  
  print "refsnpid is $refsnpid contigs is @$contigs\n";
  
  my $info;
  if ($refsnpid =~ /^\d+$/) {
    ($info) = $self->snp_dbobj->get_snp_info_by_refsnpid($refsnpid,$mouse);
  }
  $self->_store_snp_seqobj($info);
  
  $self->_pre_run($contigs);
}

=head2 _pre_run

 Title   : _pre_run
 Usage   : internal method
 Function: generate objecs for cross_match runs
 Example : Called internally to the module by the constructor $self->_pre_run($contigs);
 Returns : no
 Args    : contigs_objs
 
=cut

sub _pre_run {

  my ($self,$contigs) = @_;
  
  my @contigs=@$contigs if $contigs;
  
  my (@infos,$inseq_name,$inseq_seq,%id_hash,@newseq,$clone_seqobj);

  foreach my $contig(@contigs ) {
    my ($version,$name);
    
    if ($contig->isa("Bio::EnsEMBL::RawContig" or "Bio::EnsEMBL::Slice")) {
      $name = $contig->name;
    }
    else  {
      $name = $contig->display_id;
    }
    #if ($contig->id =~ /\w+\.(\d+)\.\d+\..*/) { ####using old_schema
    if ($name =~ /\w+\.(\d+)\.\d+\..*/) {
      $version = $1;
    }
    if (!$version) {
      $version="0";
    }
    my $id = $self->id_generation; 
    print "name is $name\n";
    $self->_version($version);
    $id_hash{$id}=$name;
    $self->_acc($name);
    #$id_hash{$id}=$contig->id; ####using old_schema
    my $seq = Bio::PrimarySeq->new( -display_id => "$id", -seq => $contig->seq);
    if ($self->_debug>1) {
      print STDERR "Created seq with id $id for ".$id_hash{$id}." contig id ".$contig->name."\n";
      #print STDERR "Created seq with id $id for ".$id_hash{$id}." contig id ".$contig->id."\n"; ####using old_schema
    }
    push(@newseq,$seq);
  }
  
  if ($self->_debug >1) {
    print STDERR "Creating crossmatches for contig ".$self->_acc."\n";
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
    print "newseq is ",$newseq->id," and length is ",length($newseq->seq())," and snp_seqobj is ",$self->_snp_seqobj->id,"\n";
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

=head2 seq_map

 Title   : seq_map
 Usage   : $crossmap->seq_map;
 Function: generate snp_seqobj and objecs for cross_match runs
 Example : see usage
 Returns : no
 Args    : input_seqs and matching contigs_objs
 
=cut
 
###for hgbase sequence mapping

sub seq_map {
  
  my ($self,$inseq_name,$inseq_seq,$contigs) = @_;
  $inseq_seq = $$inseq_seq;
  ####one "-" is left to tell snp_pos###
  $inseq_seq =~ s/\-+/-/;
  if ($inseq_seq =~ /\(/ and $inseq_seq =~ /\)/) {
    $inseq_seq =~ s/\(|\)//g;
    $self->_big_gap(1);
  }
  elsif ($inseq_seq =~ /LARGEDELETION/) {
    $inseq_seq =~ s/LARGEDELETION/-/;
    $self->_big_gap(1);
  }
  my $inseq_length = length($inseq_seq); 
  print "inseq length is $inseq_length\n";
  my $num_bases = ($inseq_seq =~ tr/[WACTG\_]//);
  my ($fseq) = split /[WACTG\-\_]+/,$inseq_seq;
  #print "THIS FSEQ $fseq and num_bases is $num_bases\n";
  my $snp_pos = length($fseq)+1;
  $self->_snp_pos($snp_pos);
  $self->_num_bases($num_bases);
  my (%id_hash1,$id1,$display_id);
  $display_id = $inseq_name;
  if (length( $inseq_name) >=15) {
    $id1 = $self->id_generation1;
    $id_hash1{$id1}=$inseq_name;
    $self->idhash1(\%id_hash1);
    $display_id = $id1;
  }
  my $snp_seqobj = Bio::PrimarySeq->new( -display_id => $display_id , -seq => $inseq_seq);
  $self->_snp_seqobj($snp_seqobj);
  $self->_snp_seq_length($inseq_length);
  
  $self->_pre_run($contigs);
}

=head2 _store_snp_seqobj 

 Title   : _store_snp_seqobj
 Usage   : internal mrthod
 Function: save snp_objecs for cross_match runs
 Example : Called internally to the module by the constructor $self->_store_snp_seqobj($info);
 Returns : no
 Args    : snp_objs
 
=cut
  
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
    $snp_seq =~ s/\s+//g;
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

=head2 idhash1
  
  Title   : idhash1
  Usage   : $obj->idhash1($newval)
 Function: Getset for idhash1 value
  Returns : value of idhash1
  Args    : newvalue (optional)
  
  
=cut
  
sub idhash1{
  my $obj = shift;
  if( @_ ) {
    my $value = shift;
    $obj->{'idhash1'} = $value;
  }
  return $obj->{'idhash1'};
  
}

sub same_version {
  my $obj = shift;
  if( @_ ) {
    my $value = shift;
    $obj->{'_same_version'} = $value;
  }
  return $obj->{'_same_version'};
  
}

sub _big_gap {

  my $obj = shift;  
  if (@_) {
    my $value = shift;
    $obj->{'_big_gap'} = $value;
  }
 return $obj->{'_big_gap'};

}

 
=head2 run

 Title   : run
 Usage   : $crossmap->run;
 Function: run cross_match and generate @fp
 Example : $crossmap->run;
 Returns : no
 Args    : no


=cut

sub run{
  my ($self) = @_;
   
  my (@fp,@final);
  #Run all the crossmatch runnables created in _pre_run
  if ($self->_debug>1) {
    print STDERR "Running crossmatches for contig ".$self->_acc."\n";
  }
  foreach my $crossmatch (@{$self->{'_crossmatch'}}) {
    $crossmatch->run();
    #Push all the feature pairs into the array
    push (@fp,$crossmatch->output);
  }
  
  #Sort the array by score
  my @sorted= sort { $b->score <=> $a->score} @fp;
  
  print STDERR "Analysing output for clone ".$self->_acc."\n".@sorted."\n";
  
 FP: foreach my $fp (@sorted) {
    if ($self->_debug>1) {
      print STDERR "Got feature pair between contig ".$fp->seqname." (".$fp->start."-".$fp->end.") and contig ".$fp->hseqname." (".$fp->hstart."-".$fp->hend.") with score ".$fp->score."\n";
    }
    my $ref = $self->idhash;
    my $ref1 = $self->idhash1;
    my %id_hash = %$ref;
    my %id_hash1 = %$ref1 if $ref1;
    my $backid = $id_hash{$fp->hseqname};
    my $backid1 = $id_hash1{$fp->seqname};
    #my $sn=$self->_acc.".".$backid;
    my ($sn, $hsn);
    
    if ($backid) {
      $hsn=$backid;
    }
    else {
      $hsn=$self->_acc;
    }
    if ($backid1) {
      $sn = $backid1;
    }
    else {
      $sn = $fp->seqname;
    }
    #####don't change input_id###
    #$sn =~ s/\_/\./g;
    #$hsn =~ s/\_/\./g;
    $fp->seqname($sn);
    $fp->hseqname($hsn);
    push (@final,$fp);
  }
  push(@{$self->{'_final'}},@final);
}


sub sort_numeric {$a<=>$b}

=head2 write_output
  
 Title   : write_output
 Usage   : $crossmap->write_output
 Function:   generate output
 Example : see usage
 Returns : output to be read into the Hit table
 Args    : no
  
  
=cut
  
=head2 output_features
  
 Title   : output_features
 Usage   : $crossmap->output_features
 Function: output features
 Example : see usage
 Returns : output features
 Args    : no
  
  
=cut
 

sub output_features {
  my ($self) = @_;
  return @{$self->{'_final'}};
}

sub write_output{ 
  my ($self) = @_;
  
  my @fp=@{$self->{'_final'}};
  
  if (!@fp) {
    print "no hits\n";
    return;
  }

  
  #$self->write_feature_to_file(@fp);
		   
  my (%final_lists, @snp_fps, $snp_type, $snp_class,$final_line, $hit_line);

  ##To find matches that cover snp_pos
    
  @snp_fps=$self->_find_fps_cover_snp("0",@fp); 
  if (@snp_fps >0) {
    if ($self->_num_bases ==1) {
      ###looks like a real snp
      $snp_type = "snp";
      $snp_class = "exact";
    }
    else {
      $snp_type = "in-del";
      $snp_class = "between" if ($self->_num_bases >1);
      $snp_class = "range" if ($self->_num_bases ==0);
      my @redo_fps;
      if ((scalar @snp_fps ==1 or $snp_fps[0]->hseqname ne $snp_fps[1]->hseqname and $snp_fps[0]->score == $snp_fps[1]->score) and  scalar @fp >1) {###if it is a in-del, we need at least two hits
	foreach my $fp (@fp) {
	  ####only to pick up the ones that are not exactly cover the snp
	  foreach my $snp_fps (@snp_fps) {
	    if ($fp ne $snp_fps and $fp->hseqname eq $snp_fps->hseqname ) {
	      push @redo_fps, $fp;
	    }
	  }
	}
	@redo_fps = sort {$b->score <=> $a->score} @redo_fps;
	my @return_fps = $self->_find_fps_cover_snp("10",@redo_fps);
	push @snp_fps, @return_fps if (@return_fps);
      }
    }
  }

  ##if no match cover snp_pos, try take 10 more bp on either side of the matches
  elsif (@snp_fps ==0)  {
    @snp_fps=$self->_find_fps_cover_snp("10",@fp);
  }

  if (@snp_fps==0 or @snp_fps >10 ) {
    if ($self->_debug>0) {
      print "snp_pos ",$self->_snp_pos," is not in matching range or more than 10 hits with same score\n";
    }
    return;
  }
  else {
    my (%hit_name);
    foreach my $fp (@snp_fps) {
      push @{$hit_name{$fp->hseqname}}, $fp;
    }
    my $num_hit = keys %hit_name;	    
    foreach my $hit (keys %hit_name) {
      my @fp = @{$hit_name{$hit}};
      print "score for snp_fps lists ",$fp[0]->score,"\n";
      
      if ($fp[0]->hseqname =~ /\w+\.(\d+)\.\d+\..*/) {
	$self->_version($1);
      }
      else {
	$self->_version ("0");
      }
      
      ### for real snps or micrositalite having 1 hit###
      if (@fp ==1) {
      #foreach my $fp (@fp) { ####changed for HAP where we need more hit, not just first hit for real snp
	if ($fp[0]->start == $fp[0]->end) {
	  $snp_class = "exact";
	  
	  if ($self->_num_bases ==1) {
	    $snp_type = "snp";
	  }
	  elsif ($self->_num_bases == 0) {
	    $snp_type = "mix";   
	  }
	}
	elsif ($fp[0]->start < $fp[0]->end and $self->_big_gap ==1) {
	  $snp_type = "microsat";
	  $snp_class = "range";
	}
	my $qual = $fp[0]->score/length($self->_snp_seqobj->seq);
	$final_line = $fp[0]->seqname."\t".$fp[0]->hseqname."\t".$self->_version."\t".$fp[0]->start."\t".$fp[0]->end."\t$snp_class\t$snp_type\t".$fp[0]->strand."\t".$fp[0]->score."\t$qual";
      #}
      }
      ### for in-dels having more hits
      else {
	my ($higher_score_fp,@fp_diff_score);
	$snp_type = "in-del";
	if (scalar @fp >2) { ## use one higher hit plus a low hit
	  
	  @fp_diff_score = sort {$b->score<=>$a->score} @fp;
	  if ($fp_diff_score[0]->score != $fp_diff_score[1]->score) {
	    $higher_score_fp = $fp_diff_score[0];
	    shift (@fp_diff_score);
	    @fp_diff_score = sort {$a->hstart<=>$b->hstart} @fp_diff_score;
	    @fp = "";
	    
	    ####always cover wider range??? if 3 hits
	    if ($higher_score_fp->hstart > $fp_diff_score[1]->hstart) {
	      @fp = ($fp_diff_score[0], $higher_score_fp);
	      print "higher score at end ",$fp[0]->hstart,"\t",$fp[1]->hend,"\n";
	    }
	    elsif ($higher_score_fp->hstart < $fp_diff_score[0]->hstart) {
	      @fp = ($higher_score_fp, $fp_diff_score[1]);
	      print "higher score at beginning ",$fp[0]->hstart,"\t",$fp[1]->hend,"\n";
	    }
	  }
	  else { ##in case top two hits are having same score###
	    @fp = sort {$a->hstart<=>$b->hstart} @fp[0,1];
	    print "here 1\n";
	  }
	}
	else { ## incase only having two hits
	  @fp = sort {$a->hstart<=>$b->hstart} @fp;
	  print "here 2\n";
	}
	my ($start, $end, $mid);

	if ($fp[0]->strand == $fp[1]->strand) {
	  if ($fp[0]->hend+1 == $fp[1]->hstart) {
	    $start = $fp[0]->hend;
	    $end = $fp[1]->hstart;
	    $snp_class = "between";
	  }
	  else {
	    $start = $fp[0]->start;
	    $end = $fp[1]->start;
	  }
	  print "here 3 ",$fp[0]->strand, "\t",$fp[0]->start,"\t",$fp[1]->start,"\n";
	}
	else {
	  print "ERROR : strand are not same for ",$fp[0]->seqname,"\n";
	}

	if ($start == $end) {
	  $snp_class = "exact";
	  $snp_type = "snp";
	}
	elsif ($end < $start) {
	  print "end<start\n";
	  $end = $start;
	  $start = $end -1;
	  $snp_class = "between";
	}
	elsif ($start < $end) {
	  $snp_class = "range" if (!$snp_class);
	}
	
	my $qual = $fp[0]->score/length($self->_snp_seqobj->seq);
	my $warn = "big_gap";

	#####if without start and end, it means two hits are on different strand, it's uncommon
	$final_line = $fp[0]->seqname."\t".$fp[0]->hseqname."\t".$self->_version."\t$start\t$end\t$snp_class\t$snp_type\t".$fp[0]->strand."\t".$fp[0]->score."\t$qual";####added qual

	if ($end-$start >50 and $self->_big_gap) {
	  $final_line .="\t$warn";
	}
	elsif (($end-$start >50 and !$self->_big_gap) or (!$end and !$start)) {
	  my @fp = sort {$b->score<=>$a->score} @fp;
	  my $start = $fp[0]->start;
	  my $end = $fp[0]->end;
	  
	  if ($start == $end) {
	    $snp_class = "exact";
	    $snp_type = "snp";
	  }
	  elsif ($start > $end) {
	    $end = $start;
	    $start = $end -1;
	    $snp_class = "between";
	  }
	  elsif ($start < $end) {
	    $snp_class = "range";
	  }
	  my $qual = $fp[0]->score/length($self->_snp_seqobj->seq);
	  $final_line = $fp[0]->seqname."\t".$fp[0]->hseqname."\t".$self->_version."\t$start\t$end\t$snp_class\t$snp_type\t".$fp[0]->strand."\t".$fp[0]->score."\t$qual";
	  print "end-start>50, but not microsat and largedeletion, higher score chosen\n";
	}
      }
      
      print "final_line is $final_line\n";
      
      if ($num_hit ==1) { 
	$hit_line .= "$final_line\n" if ($final_line);
      }
      elsif ($num_hit >1) {
	$hit_line .= "MORE_HITS $final_line\n" if ($final_line);
      }
    }
    chomp($hit_line);
    return ($hit_line);####return if having 1 to 3 hits#####
  }
}

=head2 _find_fps_cover_snp

 Title   : _find_fps_cover_snp
 Usage   : $self->_find_fps_cover_snp("10",@fp);
 Function: find which fp cover snp region
 Example : see usage
 Returns : highest(s) fps
 Args    : "10" @fp 1 is the num to be added or take away from start and end


=cut

sub _find_fps_cover_snp{

  my ($self, $num, @fp) = @_;

  my (@snp_fps, @cover_snps, @snps_est, %done, %rec_hseqname);
  
  if (@fp) {
    my $first_hname = $fp[0]->hseqname;
    my $first_score = $fp[0]->score;
    
    foreach my $fp ( @fp ) {
      
      my $snp_pos = $self->_snp_pos;
      
      #if ($self->_debug>0) {
      print "In write_out ",$fp->seqname," ",$fp->hseqname," ",$fp->score," snp_pos is ",$self->_snp_pos," ",$fp->start," ",$fp->end," ",$fp->hstart," ",$fp->hend," ",$fp->strand,"\n";
      #}
      
      if (($snp_pos>=($fp->start-$num) and $snp_pos<=($fp->end+$num) or $snp_pos>=($fp->end-$num) and $snp_pos<=($fp->start+$num)) and $fp->start != $fp->end) {
	push @cover_snps, $fp;
	
	##the part that is part of highest score sequence, if not, the score is same as the highest core or position is very close to the highest
	if ($fp->score == $first_score or ($snp_fps[0] and (abs($fp->end - $self->_snp_pos) < 10 or abs($fp->start - $self->_snp_pos) < 10) and $rec_hseqname{$fp->hseqname})) {##allow two diff hseqname with lower score for each of them
	#if ($fp->score == $first_score or $fp->score > $first_score -1000 or ($snp_fps[0] and (abs($fp->end - $self->_snp_pos) < 10 or abs($fp->start - $self->_snp_pos) < 10) and $rec_hseqname{$fp->hseqname})) {##allow two diff hseqname with lower score for each of them######this one is for HAP mapping where right position is not always highest score
	  
	  $fp = $self->_record_snp($fp);
	  push (@snp_fps, $fp); 
	  $rec_hseqname{$fp->hseqname}=1;
	}
      }
    }
  }
  ######for est (for all now), the cover_snp match may not as hight score as others
  if (scalar @snp_fps == 0 and scalar @cover_snps >0 and $num >0) {
    my @cover_snps = sort {$b->score<=>$a->score} @cover_snps;
    my $first_score = $cover_snps[0]->score;
    foreach my $fp (@cover_snps) {
      if ($fp->score == $first_score) {
	$self->_record_snp($fp);
	push @snps_est, $fp;
      }
    }
    return @snps_est;
  }
  else {
    return (@snp_fps);
  }
}

=head2 _record_snp

 Title   : _record_snp
 Usage   : $self->_record_snp($fp)
 Function: to store fp that cover snp region
 Example : see usage
 Returns : no
 Args    : $fp


=cut

sub _record_snp{

  my ($self, $fp) = @_;

  if ($self->_debug>1) {
    print "in record_snp first ",$fp->seqname,"\t",$fp->hseqname,"\t",$fp->score,"\t",$self->_version,"\t",$self->_num_bases, "\t",$fp->start, "\t",$fp->end,"\t",$fp->hstart,"\t",$fp->hend,"\t",$fp->strand,"\n";
  }
  my $snp_pos = $self->_snp_pos;
  my ($new_strand, $new_start,$new_end);
  my $num_bases = $self->_num_bases;
  if ($fp->strand ==-1 or $fp->hstrand ==-1) {
    $new_strand=-1;
    if ($num_bases ==0) {
      $new_start = $fp->hstart+($fp->end-$snp_pos);
    }
    else {
      $new_start = $fp->hstart+($fp->end-($snp_pos+$num_bases-1));
    }
  }
  else {
    $new_strand=1;
    $new_start = $fp->hstart+($snp_pos-$fp->start);
  }
  if ($num_bases >=1) {
    $new_end = $new_start+$num_bases-1;
  }
  elsif ($num_bases ==0) {
    $new_end = $new_start;
  }
  
  $fp->start($new_start);
  $fp->end($new_end);
  $fp->strand($new_strand);
  
  print "in record_snp last ",$fp->seqname,"\t",$fp->hseqname,"\t",$fp->score,"\t",$self->_version,"\t$new_start\t$new_end\t\tsnp_map\t$new_strand\n";
  return ($fp);
  
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

sub _est{
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'_est'} = $value;
    }
    return $self->{'_est'};
 }

1;
#$refsnpid,$snpclass,$snptype,$observed,$seq5,$seq3,$acc,$version,$start,$end,$strand
