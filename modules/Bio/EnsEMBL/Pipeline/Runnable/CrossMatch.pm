#
# Ensembl module for Bio::EnsEMBL::Pipeline::Runnable::CrossMatch
#
# Cared for by Ewan Birney <birney@ebi.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Pipeline::Runnable::CrossMatch - CrossMatch featurepair generator

=head1 SYNOPSIS

    $rb = Bio::EnsEMBL::Pipeline::Runnable::CrossMatch->new( -seq1 => $seq1,-seq2 =>$seq2);
    $rb->run();
   
    @fp = $rb->output;



=head1 DESCRIPTION

This runs cross match (using James Gilberts/Tim Hubbard system) and provides
feature pair output

=head1 AUTHOR - Ewan Birney

This module is part of the Ensembl project http://www.ensembl.org

Email birney@ebi.ac.uk

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::Pipeline::Runnable::CrossMatch;
use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::RootI

use Bio::EnsEMBL::Pipeline::RunnableI;
use Bio::EnsEMBL::FeatureFactory;
use Bio::SeqIO;
# no need to 'use CrossMatch' as we embed the relevant pacakges here

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI);

# new() is written here 

sub new {
  my($class,@args) = @_;
  
  my $self = {};
  bless $self,$class;
  
  $self->{'_fp_array'} =[];
  my ($seq1,$seq2,$workdir,$score,$minmatch,$masklevel,$debug) = $self->_rearrange([qw(SEQ1 SEQ2 WORKDIR SCORE MINMATCH MASKLEVEL)],@args);
  if( !defined $seq1 || !defined $seq2 ) {
      $self->throw("Must pass in both seq1 and seq1 args");
  }

  $self->seq1($seq1);
  $self->seq2($seq2);
  if( $workdir) { 
      $self->workdir($workdir); 
  } else {
      $self->workdir("/tmp");
  }
  if($score) { 
      $self->score($score); 
  } else {
      $self->score(5000);
  }
  if( defined $minmatch ) {
    $self->minmatch($minmatch);
  }
  if( defined $masklevel ) {
    $self->masklevel($masklevel);
  }
  
# set stuff in self from @args
  return $self;
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
   my ($self,@args) = @_;

   # dump sequences to work directory
   my $file1 = "crossmatch1.".$$;
   my $file2 = "crossmatch2.".$$;
   
   open(F,">".$self->workdir."/$file1") || $self->throw("Yikes. Cannot make $file1 $!");
   my $seqout = Bio::SeqIO->new(-fh => \*F,-format => 'fasta');
   $seqout->write_seq($self->seq1);
   close(F);

   open(F,">".$self->workdir."/$file2") || $self->throw("Yikes. Cannot make $file2 $!");
   my $seqout2 = Bio::SeqIO->new(-fh => \*F,-format => 'fasta');
   $seqout2->write_seq($self->seq2);
   close(F);
   
   #build crossmatch factory
   my $cmf = CrossMatch::Factory->new($self->score,$self->minmatch,$self->masklevel);
   $cmf->dir($self->workdir);
   $cmf->alignments(1);

   print STDERR "About to run crossmatch from runnable...\n";
   #run crossmatch factory
   my $cm;
   eval {
       $cm = $cmf->crossMatch($file1,$file2);
   };
   print STDERR "Finished crossmatch in runnable...\n";

   if( $@ ) {
       unlink($file1);
       unlink($file2);
       $self->throw("Unable to run cross match! (probably no crossmatch executable).\n$@");
   }
   if (!$cm) {
     #print "0 match entries for ",$self->seq1->id,"\n";
     return;
   }
   #process alignments above score
   foreach my $fp ( $cm->fp($self->score) ) {
       my ($seq1,$seq2,$score,$start,$end,$hstart,$hend) = split(/:/,$fp);

       my ($strand,$hstrand,$swap);

       if( $start > $end ) {
	   $strand = -1;
	   $swap = $start;
	   $start = $end;
	   $end = $swap;
       } else {
	   $strand = 1;
       }

       if( $hstart > $hend ) {
	   $hstrand = -1;
	   $swap = $hstart;
	   $hstart = $hend;
	   $hend = $swap;
       } else {
	   $hstrand = 1;
       }


       $fp = Bio::EnsEMBL::FeatureFactory->new_feature_pair();
       #print STDERR "Processing FP with $start-$end to $hstart-$hend\n";

       $fp->start($start);
       $fp->end($end);
       $fp->strand($strand);
       $fp->seqname($self->seq1->id);
       $fp->hstart($hstart);
       $fp->hend($hend);
       $fp->hstrand($hstrand);
       $fp->hseqname($self->seq2->id);
       $fp->score($score);
       
       $self->add_fp($fp);
   }

   #unlink files

   unlink($file1);
   unlink($file2);

   # cool!


}

=head2 output

 Title   : output
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub output{
   my ($self,@args) = @_;

   return @{$self->{'_fp_array'}};
}



=head2 add_fp

 Title   : add_fp
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub add_fp{
   my ($self,@args) = @_;

   push(@{$self->{'_fp_array'}},@args);
}


=head2 workdir

 Title   : workdir
 Usage   : $obj->workdir($newval)
 Function: 
 Example : 
 Returns : value of workdir
 Args    : newvalue (optional)


=cut

sub workdir{
   my ($obj,$value) = @_;
   if( defined $value) {
       $obj->{'workdir'} = $value;
   }
   return $obj->{'workdir'};

}

=head2 score

 Title   : score
 Usage   : $obj->score($newval)
 Function: could set and return the score value
 Example : 
 Returns : value of minscore option used by crossmatch
 Args    : newvalue (optional)


=cut

sub score{
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'score'} = $value;
    }
    return $obj->{'score'};
}


=head2 minmatch

 Title   : minmatch
 Usage   : $obj->minmatch($newval)
 Function: could set and return the minmatch value
 Example : 
 Returns : value of minmatch option used by crossmatch
 Args    : newvalue (optional)


=cut

sub minmatch {
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'minmatch'} = $value;
    }
    return $obj->{'minmatch'};

}

=head2 masklevel

 Title   : masklevel
 Usage   : $obj->masklevel($newval)
 Function: could set and return the masklevel value
 Example : 
 Returns : value of masklevel option used by crossmatch
 Args    : newvalue (optional)


=cut

sub masklevel {
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'masklevel'} = $value;
    }
    return $obj->{'masklevel'};

}

=head2 seq1

 Title   : seq1
 Usage   : $obj->seq1($newval)
 Function: 
 Example : 
 Returns : value of seq1
 Args    : newvalue (optional)


=cut

sub seq1{
   my ($obj,$value) = @_;
   if( defined $value) {
      if( !ref $value || !$value->isa('Bio::PrimarySeqI') ) {
	  $obj->throw("$value is not a PrimarySeqI object. Cannot throw");
      }
      my $selfseq = Bio::PrimarySeq->new( -display_id => $value->id , -seq => $value->seq);
      if( $selfseq->length == 0 ) {
	  $obj->throw("attempting to crossmatch seemingly 0 length sequence!");
      }
      $obj->{'seq1'} = $selfseq;
    }
    return $obj->{'seq1'};

}


=head2 seq2

 Title   : seq2
 Usage   : $obj->seq2($newval)
 Function: 
 Example : 
 Returns : value of seq2
 Args    : newvalue (optional)


=cut

sub seq2{
   my ($obj,$value) = @_;
   if( defined $value) {
      if( !ref $value || !$value->isa('Bio::PrimarySeqI') ) {
	  $obj->throw("$value is not a PrimarySeqI object. Cannot throw");
      }
      my $selfseq = Bio::PrimarySeq->new( -display_id => $value->id , -seq => $value->seq);
      if( $selfseq->length == 0 ) {
	  $obj->throw("attempting to crossmatch seemingly 0 length sequence!");
      }
      $obj->{'seq2'} = $selfseq;

    }
    return $obj->{'seq2'};

}



=head1 Embedded CrossMatch module

To prevent another outside dependency, embedding the crossmatch module
in here directly.

=cut


package CrossMatch::Factory;

require 5.004;

use strict;
use Carp;
use Cwd;

# Makes a new CrossMatch factory object
sub new {
    # order here is important, as we used to only specify min score
    my ($pkg,$minscore,$minmatch,$masklevel) = @_;
    my( @dir );

    # Default to cwd if no directories supplied
    if (@_) {
	@dir = _dirs(@_);
    } else {
	$dir[0] = cwd();
    }

    $minmatch = $minscore unless (defined $minmatch); # $minmatch = $minscore compatatible with cvs version 1.7
    $minscore = 30 unless (defined $minscore); 
    $masklevel = 101 unless (defined $masklevel);

    return bless {
	dir => [ @dir ],
	minMatch  => $minmatch,
	minScore  => $minscore,
	maskLevel => $masklevel, # 101, Show all matches
	alignments => 0,  # don't keep/parse alignments by default
	_extn => ['']
	}, $pkg;
}

# Accessor methods for cross_match's minmatch, minscore and masklevel parameters
BEGIN {
    foreach my $func (qw( minMatch minScore maskLevel )) {
        no strict 'refs';
        
        *$func = sub {
            my( $matcher, $num ) = @_;
            
            if ($num) {
	        $num =~ /^\d+$/ or croak "non-integer '$num' supplied to $func";
	        $matcher->{$func} = $num;
            } else {
	        return $matcher->{$func};
            }
        }
    }
}

sub dir {
    my $matcher = shift;

    if (@_) {
	$matcher->{'dir'} = [ _dirs(@_) ];
    } else {
	return @{$matcher->{'dir'}};
    }
}
sub _dirs {
    my @dirs = @_;
    foreach (@dirs) {
	# Strip trailing "/" from name
	s{/$}{};
    }
    return @dirs;
}

# flag to switch alignment collection on
sub alignments {
    my $matcher = shift;
    $matcher->{'alignments'} = 1;
}

sub extn {
    my $matcher = shift;
    my( @extn );

    if (@_) {
        my( @values ) = @_;
	foreach my $extn (@values) {
	    # Add a leading dot
	    $extn =~ s|^\.?|\.| unless $extn eq '';
	    push( @extn, $extn );
	}
	$matcher->{'_extn'} = [ @extn ];
    } else {
	# Return all extensions
	return @{$matcher->{'_extn'}};
    }
}

sub findFile {
    my $matcher = shift;
    my $name = shift;

    foreach my $dir ($matcher->dir()) {
	foreach my $file (map "$dir/$name$_", ($matcher->extn()) ) {
	    # Check for a readable file of this name
	    return $file if -r $file;
	}
    }
    # Return undef if we didn't find a matching filename
    return;
}

sub _save_hit{
    my($match,$rares,$raalign)=@_;

    # only save if something to save
    if(@$rares){

	# parsing to deal different numbers of items
	# if 13 then second sequence is complement
	# if 12 then its not, so insert 'W'
	my @nres;
	if (@$rares == 13) {
            @nres=( @$rares[0..9, 12, 11, 10]       );
	} else {
            @nres=( @$rares[0..7], 'W', @$rares[8..11] );
	}

	# alignment is stored in a directional fashion: x->y maps to a->b where y>x
	my $raw_align;
	if(@$raalign){
	    # reverse if required
	    my $st1=$nres[5];
	    my $en1=$nres[6];
	    my $st2;
	    my $en2;
	    my $dirl;
	    if($$raalign[2] eq 'C'){
		$$raalign[3]='C';
		$$raalign[0]=reverse($$raalign[0]);
		$$raalign[1]=reverse($$raalign[1]);
		$dirl='C';
		$st2=$nres[11];
		$en2=$nres[10];
	    }else{
		$st2=$nres[10];
		$en2=$nres[11];
	    }

	    # check length
	    my $l1=length($$raalign[0]);
	    my $l2=length($$raalign[1]);
	    if($l1!=$l2){
		print "Lengths of alignment are different [$l1,$l2]\n$$raalign[0]\n$$raalign[1]\n";
		print join(',',@nres)."\n";
		die "failed";
	    }

	    # walk along sequence in blocks
	    $raw_align="$st1:$dirl$st2";
	    my $seq1=$$raalign[0];
	    my $seq2=$$raalign[1];
	    {
		my($s1a,$s1b,$s1c,$s2a,$s2b,$s2c);
		if($seq1=~/^([^\-]+)(\-+)(\S+)$/){
		    ($s1a,$s1b,$s1c)=($1,$2,$3);
		}else{
		    $s1a=$seq1;
		    $s1b=$s1c='';
		}
		if($seq2=~/^([^\-]+)(\-+)(\S+)$/){
		    ($s2a,$s2b,$s2c)=($1,$2,$3);
		}else{
		    $s2a=$seq2;
		    $s2b=$s2c='';
		}
                #print STDERR "heads are ",substr($s1a,0,10),";",substr($s1b,0,10),":",substr($s2a,0,10),";",substr($s2b,0,10),"\n";

		# escape if no more gaps
		next if(length($s1c)==0 && length($s2c)==0);
		# do shortest first
		my $lab;
                if( length($s1a.$s1b) == 0 ||
			length($s2a.$s2b) == 0 ) {
			print STDERR "Dodgy alignment processing catch! Bugging out\n";
			next;
		}
		if(length($s1a.$s1b)<length($s2a.$s2b)){
		    # update seq1
		    $lab=length($s1a.$s1b);
		    $st1+=length($s1a);
                    #print STDERR "st1 is $st1 with $lab\n";
         	    $seq1=$s1c;
		    #print STDERR "New head is ",substr($seq1,0,10),"\n";
		    # update seq2
		    $seq2=~s/^\S{$lab}//;
                    #print STDERR "new $lab.. seq2 head is ",substr($seq2,0,10),"\n";
		}else{
		    # update seq2
		    $lab=length($s2a);
		    $seq2=$s2c;
		    # update seq1
		    my $l2ab=length($s2a.$s2b);
		    $seq1=~s/^\S{$l2ab}//;
		    $st1+=$l2ab;
		}
		if($dirl eq 'C'){
		    $st2-=$lab;
		}else{
		    $st2+=$lab;
		}
		$raw_align.=",$st1:$dirl$st2";
		redo;
	    }
	    $raw_align.=",$en1:$dirl$en2";

	}
	$match->hit(@nres,$raw_align);

	# clear data
	@$rares=();
	@$raalign=();
    }
}
    
sub crossMatch {
    my $matcher = shift;
    my( @seq ) = (@_)[0,1];
    
    defined $seq[1]
	or croak qq(not enough arguments to crossMatch: "@_");
    
    # Find files with these names
    foreach (@seq) {
	my $name = $_;
	$_ = $matcher->findFile( $_ )
	    or croak qq(could not find a file for "$name");
    }
    

    ### Do cross_match on the found files ###
    my $minM = $matcher->minMatch();
    my $minS = $matcher->minScore();
    my $mask = $matcher->maskLevel();
    
    # Make new match object, saving parameters
    my @res;
    my @align;
    my $ialign;
    my $match = CrossMatch->new(@seq, $minM, $minS, $mask);
    print STDERR "Options used by crossmatch -minmatch $minM -minscore $minS -masklevel $mask\n";
    my $cross_command="/work2/elia/bin/cross_match -minmatch $minM -minscore $minS -masklevel $mask";
    if($matcher->{'alignments'}){
	$cross_command.=" -alignments";
    }
    $cross_command.=" @seq 2>/dev/null |";
    print STDERR "opening cross match pipe\n";
    open( CROSS_MATCH, $cross_command )
        or croak "Can't open pipe '$cross_command' : $!";
    while (<CROSS_MATCH>) {
	# process alignment lines if requested
        #print STDERR "Processing....$_";

	if($matcher->{'alignments'}){
	    if(/^(\w*)\s+(\S+)\s+(\d+)\s+(\S+)\s+(\d+)$/){
		if($2 eq $res[4] && $ialign==1){
		    $align[0].=$4;
		    $align[2]=$1;
		    $ialign=2;
		}elsif(($2 eq $res[8] || $2 eq $res[9]) && $ialign==2){
		    $align[1].=$4;
		    $align[3]=$1;
		    $ialign=1;
		}else{
		    die "alignment parsing error in Crossmatch.pm\n  $_";
		}
	    }
	}

	# this is used to exclude all except alignment summary lines
	next unless /\(\d+\).*\(\d+\)/;

        # Remove parentheses and asterisks
	tr/()\*//d;
        
	# save previous hit, with alignments
	&_save_hit($match,\@res,\@align);
	
	@res = split;
	$ialign=1;
    }
    print STDERR "About to process...\n";
    &_save_hit($match,\@res,\@align);
    print STDERR "saved hits \n";
        
    # Check exit status of cross_match
    unless (close CROSS_MATCH) {
        my $error = $! ? "Error from cross_match: $!"
                       : "Error: cross_match exited status $?";
	croak "$error\nCommand: '$cross_command'";
    }
    
    # Remove cross_match log file
    unlink "$seq[0].log";

    # Return Match object
    return $match;
}

# Pre-rearranged fields
# 0    1    2    3    4        5     6     7       8 9        10      11    12

# 98   0.00 0.00 0.00 130N4    1     104   84065   W 92M18    68800   68903 0
# 98   0.00 0.00 0.00 92M18    68800 68903 0       W 130N4    1       104   84065
# 8251 0.00 0.00 0.00 130N4    1     84169 0       W 130N4    1       84169 0
# 103  0.00 0.00 0.00 CFAT5    20771 20874 (0)     W A1280    1       104   (22149) | W Bs
# 103  0.00 0.00 0.00 CFAT5    20771 20874 (0)     C A1280.RC (0)     22253 22150   | C Ar Br
# 103  0.00 0.00 0.00 CFAT5.RC 1     104   (20770) C A1280    (22149) 104   1       | C As Bs
# 26355  1.37 0.42 0.36  Em:AP000350    133977 162328 (0)  C Em:AP000352   (125738) 28369     1

# 120 16.53 0.00 0.00  bK363A12    32474 32715 (50)  C cE129H9   (4343) 32827 32586  
# 100  0.00 0.00 0.00  bK363A12    32666 32765 (0)    cE129H9        1   100 (37070) *
#  * indicates that there is a higher-scoring match whose domain partly includes the domain of this match.

package CrossMatch;

require 5.004;

use strict;
use Carp;

sub new {
    my $pkg = shift;

    return bless {
	aFile     => $_[0],
	bFile     => $_[1],
	minM      => $_[2],
	minS      => $_[3],
	mask      => $_[4],
	hit       => [],
	active    => [],
        raw_align => [],
	}, $pkg;
}

# Create access methods for fields which record
# parameters given to cross_match
BEGIN {
    foreach my $field (qw( aFile bFile minM minS mask)) {
        no strict 'refs';
        
        *$field = sub {
            my( $match ) = @_;
            
            return $match->{$field};
        }
    }
}

# Create access methods to access the data from the matches
BEGIN {
    my( %fields );
    {
        my $i = 0;
        %fields = map {$_, $i++} qw( score pSub pDel pIns
                                     aName aStart aEnd aRemain
                                     strand
                                     bName bStart bEnd bRemain raw_align);
    }

    foreach my $func (keys %fields) {
        no strict 'refs';
        
        my $i = $fields{ $func };
        
        *$func = sub {
            my( $match, @rows ) = @_;
        
            if (wantarray) {
	        # Extract the requested values
	        unless (@rows) {
	            # Get a list of all the row indices
	            @rows = @{$match->{'active'}};
	        }
                
	        # Make a vertical slice through the hits
                return map $match->{'hit'}[$_][$i], @rows;
            } else {
	        # Return just one value in scalar context
	        if (defined $rows[0]) {
	            # Get field from row requested
	            return $match->{'hit'}[$rows[0]][$i];
	        } else {
	            # Just get value from first row in active list
	            my $row = $match->{'active'}[0];
	            return $match->{'hit'}[$row][$i];
	        }
            }
        }
    }
}

# Functions which provide access to data fields by sequence name
BEGIN {
    my %funcNames = (
		     start  => ['aStart', 'bStart'],
		     end    => ['aEnd', 'bEnd'],
		     remain => ['aRemain', 'bRemain']
		     );
    
    foreach my $call (keys %funcNames) {
        no strict 'refs';
        
        *$call = sub {
            my $match = shift;
            my $name = shift;
            my( $a, $b, $i, @newActive, $func );

            # Look through aName and bName columns for name
            foreach (@{$match->{'active'}}) {
	        # Compare $name to aName field
	        if ($name eq $match->{'hit'}[$_][4]) {
	            $a = 1; push( @newActive, $_ );
	        }
	        # Compare $name to bName field
	        if ($name eq $match->{'hit'}[$_][9]) {
	            $b = 1;
	            # Add to active list unless already have it
                    my $last = $newActive[$#newActive];
	            unless (defined($last) and $last == $_) {
		        push( @newActive, $_ );
	            }
	        }
            }

            # Save the new active list
            #warn "New active list is: (", join(' ', map "'$_'", @newActive), ")\n";
            $match->{'active'} = [ @newActive ];

            # Did we find the sequence name in the "a" or "b" column?
            if ($a) {
	        if ($b) {
	            carp qq($name matches in both "a" and "b" columns -- choosing "a");
	        }
	        $i = 0;
            } elsif ($b) {
	        $i = 1;
            }

            # Call function if found match to $name, or else return nothing
            if (defined $i) {
	        $func = $funcNames{$call}[$i];
	        $match->$func(@_);
            } else {
	        warn qq("$name" not found);
	        return;
            }
        }
    }
}

# Store data from a hit
sub hit {
    my $match = shift;

    if (@_ == 14) {
	push( @{$match->{'hit'}}, [@_]);
	push( @{$match->{'active'}}, $#{$match->{'hit'}} );
    } else {
	confess "Bad number of elements (", scalar @_, ") in '@_'";
    }
}

# Give list of active indices, or list of indices
# for all hits if no filter has been applied
sub list {
    my $match = shift;
    return @{$match->{'active'}};
}

# Restore active indices to list all hits
sub unfilter {
    my $match = shift;
    $match->{'active'} = [ (0..$#{$match->{'hit'}}) ];
    return (0..$#{$match->{'hit'}});
}

# Takes a ref to a subroutine as its argument
# which returns true or false for the line
# given as an argument
sub filter {
    my $match = shift;
    my $sub = shift;

    my( @true );
    foreach my $row ($match->list) {
	if ( &$sub($match, $row) ) {
	    push @true, $row;
	}
    }
    
    # Save the filter
    $match->{'active'} = [ @true ];
}

sub endHits {
    my $match = shift;
    my $row = shift;
    my $end_hit;
    
    $end_hit++ if $match->aStart($row)  == 1;
    $end_hit++ if $match->aRemain($row) == 0;
    $end_hit++ if $match->bStart($row)  == 1;
    $end_hit++ if $match->bRemain($row) == 0;
    
    if ($end_hit > 1) {
	return 1;
    } else {
	return;
    }
}

# Return number of hits
sub count {
    my $match = shift;
    
    return scalar @{$match->{'active'}};
}

sub seqlength {
    my( $match, $name ) = @_;
    return $match->end($name) + $match->remain($name)
}

# calls translation
# by default a2b returns
# if request by row
#  number
# or if request by name
#  name:number of the best scoring match in scalar mode
# or
#  name:number:score of each match in array mode
sub a2b {
    my( $match, $row, $index ) = @_;

    # row could be an number or a string
    # only need to look through 'a' since a2b

    # Look through aName and bName columns for name
    my @newActive;
    foreach (@{$match->{'active'}}) {
	# Compare $name to aName field
	if ($row eq $match->{'hit'}[$_][4]) {
	    push( @newActive, $_ );
	}
    }

    my $flag_no_name;
    if(!scalar(@newActive)){
	if($row=~/^\d+$/){
	    @newActive=($row);
	    $flag_no_name=1;
	}else{
	    print "WARN: Crossmatch->a2b: $row not found\n";
	    return;
	}
    }

    my @hits;
    my $maxscore;
    my $besthit;
    foreach my $row (@newActive){
	my $bname=$match->bName($row);
	my $score=$match->score($row);
	my $num=&_trans_a2b($match->raw_align($row), $index);
	my $hit;
	if($flag_no_name){
	    $hit=$num;
	}else{
	    $hit="$bname:$num";
	}
	if($score>$maxscore){
	    $maxscore=$score;
	    $besthit=$hit;
	}
	push(@hits,"$hit:$score");
    }

    # return all values if array, else best one only
    if (wantarray) {
	return @hits;
    }else{
	return $besthit;
    }
}

sub _trans_a2b{
    my($align,$n)=@_;
    my($st1p,$st2p);
    $st1p=$st2p=0;
    foreach my $pair (split(',',$align)){
	my $dirl;
	my($st1,$st2)=split(':',$pair);
	if($st2=~/C(\d+)/){
	    $st2=$1;
	    $dirl='C';
	}

	if($n<$st1 && $n>=$st1p){

	    if($st1p==0){
		# if lower than start then no match
		return -1;
	    }

	    my $n2;
	    if($dirl){
		$n2=$st2p-$n+$st1p;
		if($n2<=$st2){
		    # undefined as in a gap
		    return -2;
		}else{
		    return $n2;
		}
	    }else{
		$n2=$st2p+$n-$st1p;
		if($n2>=$st2){
		    # undefined as in a gap
		    return -2;
		}else{
		    return $n2;
		}
	    }
	}

	# copy to 'previous'
	($st1p,$st2p)=($st1,$st2);
    }
    # if higher than end then no match
    return -3;
}

# output's full featurepair list
sub fp {
    my( $match, $minscore ) = @_;
    # loop over all rows
    my @hits;
    foreach my $row (@{$match->{'active'}}) {
	my $score=$match->score($row);
	next if($score<$minscore);
	my $aname=$match->aName($row);
	my $bname=$match->bName($row);
	my @sted=&_trans_fp($match->raw_align($row));
	foreach my $sted2 (@sted){
	    push(@hits,join(":",$aname,$bname,$score,@$sted2));
	}
    }
    return @hits;
}

sub _trans_fp{
    my($align)=@_;
    my($st1p,$st2p);
    $st1p=$st2p=0;
    my @hits;
    foreach my $pair (split(',',$align)){
	my $dirl;
	my($st1,$st2)=split(':',$pair);
	if($st2=~/C(\d+)/){
	    $st2=$1;
	    $dirl='C';
	}

	# only output if previous saved
	if($st1p){
	    # calculate ends
	    my $l1=$st1-$st1p;
	    my $l;
	    if($dirl){
		$l=$st2p-$st2;
	    }else{
		$l=$st2-$st2p;
	    }
	    # if exact match, must be end, so need full length
	    if($l1==$l){
		$l++;
	    }elsif($l1<$l){
		$l=$l1;
	    }
	    $l1=$l-1;
	    my $ed1p=$st1p+$l1;
	    my $ed2p;
	    if($dirl){
		$ed2p=$st2p-$l1;
	    }else{
		$ed2p=$st2p+$l1;
	    }
	    push(@hits,[$st1p,$ed1p,$st2p,$ed2p]);
	}
	($st1p,$st2p)=($st1,$st2);
    }
    # if higher than end then no match
    return @hits;
}

1;


__END__

=head1 NAME - CrossMatch.pm

=head1 DESCRIPTION

Module to provide object-oriented access to Phil Green's B<cross_match>
Smith-Waterman alignment program.

=head1 SYNOPSIS

	use CrossMatch;
	
	# Create a factory object
	$matcher = CrossMatch::Factory->new( '/nfs/disk2001/this_dir',
	                                     '/home/jgrg/that_dir' );
	$matcher->minMatch( 80 );

        # process full crossmatch alignments
        $matcher->alignments;

	# Match two fasta fomat sequence files, generating
        # a match object
	$cm = $matcher->crossMatch( 'dJ334P19.seq', 'cB49C12.aa' );

=head1 FUNCTIONS

=over 4

=item new

Create a new CrossMatch factory object.  Optional arguments to new is a
list of directories to search (see B<dir> below).

=item crossMatch

Do a cross_match on the two B<fasta formatted> sequence files supplied
as agruments, returning a B<CrossMatch> object.  The B<CrossMatch>
object records the parameters used in the match, and the match data
returned my cross_match (if any).

I<Note:> The two fasta files may each contain multiple sequences.

=item dir

Change the list of directories in which B<crossMatch> searches for
files.  The list of directories supplied to B<dir> completely replaces
the existing list.  Directories are searched in the order supplied.

This list defaults to the current working directory when the
B<CrossMatch::Factory> object is created.

=item extn

A convenience function which allows you to supply a list of filename
extensions to be considered by B<crossMatch> when finding a file to
pass to cross_match.  For example:

	$matcher->extn( 'seq', 'aa' );
	$cm = $matcher->crossMatch( 'dJ334P19', 'cB49C12' );

B<crossMatch> will look for files "dJ334P19.seq", "dJ334P19.aa".  If
no B<extn>s had been set, then only "dJ334P19" would have been
searched for.

=item minMatch

Set the B<minmatch> parameter passed to cross_match.  Defaults to 50.

=item minScore

Set the B<minscore> parameter passed to cross_match.  Defaults to 50.

=item maskLevel

Set the B<masklevel> parameter passed to
cross_match.  Defaults to 101, which displays all
overlapping matches.

=item alignments

Causes the full crossmatch alignments to be parsed and stored in the
Crossmatch object.  These can be accessed by the a2b and fp methods.

=back

=head1 CrossMatch

Data stored in B<CrossMatch> objects, generated by a
B<CrossMatch::Factory>, can be retrieved with a variety of queries.
All the matches are stored internally in an array, in the order in
which they were generated by cross_match.  You can apply filters and
sorts (sorts are not yet implemented) to the object, which re-order or
hide matches.  These operations save their changes by altering the
active list, which is simply an array of indices which refer to
elements in the array of matches.  The data access funtions all
operate through this active list, and the active list can be reset to
show all the matches found by calling the B<unfilter> method.

=over 4

=item SYNOPSIS

	my $numberOfHits = $cm->count;
	foreach my $hit ($cm->list) {
	    print $cm->score($hit);
        }

	$firstName = $cm->aName;
	@allscores = $cm->score();
	@someSores = $cm->score(0,3,4);


=item FUNCTIONS

=over 4

=item count

Returns a count of the number of hits in the active list

=item list

Returns a list of array indices for the current list.

=item filter

	$cm->filter(\&CrossMatch::endHits);

Takes a reference to a subroutine (or an anonymous subroutine), and
alters the active list to contain only those hits for which the
subroutine returns true.  The subroutine should expect to be passed a
B<CrossMatch> object, and an integer corresponding to the index of a
match.  Applying a series of filters to the same object removes
successively more objects from the active list.

=item unfilter

Resets the active list to include all the hits, in the order in which
they were generated by cross_match.  Returns a list of the active
indices.

=item FILTERS

=over 4

=item endHits

A filter which returns true if one of the sequnces in a hit has its
end in the hit.

=back

=back

=item PARAMETERS

The parameters used by cross_match when doing the match can be
retrieved from the resulting Match object.  The two sequences matched
are labelled I<a> and I<b>, where I<a> is the sequence from the first
file passed to B<crossMatch>, and I<b>, the sequence from the second
file.  For example:

	$path_to_a = $cm->aFile;

retrieves the path of the file supplied as the first argument to
cross_match.

=over 4

=item aFile bFile

The full paths to the files used for sequences I<a> and I<b>.

=item minM

The minmatch parameter used.

=item minS

The minscore parameter used.

=item mask

The masklevel parameter used.

=back

=item DATA FIELDS

Syntax:

	$field = $match->FIELD();
	$field = $match->FIELD(INDEX);
	@fields = $match->FIELD();
	@fields = $match->FIELD(LIST);

Examples:

	$firstScore = $match->score(0);
	@aEnds = $match->aEnd();

These methods provide access to all the data fields for each hit found
by cross_match.  Each of the fields is described below.  In scalar
context, a single data field is returned, which is the field from the
first row in the active list if no I<INDEX> is given.  In array
context it returns a list with the I<FIELD> from all the hits in the
active list when called without arguments, or from the hits specified
by a I<LIST> of indices.

=over 4

=item score

The Smith-Waterman score of a hit, adjusted for the complexity of the
matching sequence.

=item pSub

The percent substitutions in the hit.

=item pDel

The percent deletions in the hit (sequence I<a> relative to sequence I<b>).

=item pIns

The percent insertions in the hit (sequence I<a> relative to sequence I<b>).

=item aName bName

ID of the first sequence (sequence I<a>) and second seqeunce (I<b>) in
the match, respecitvely.

=item aStart bStart

Start of hit in I<a>, or I<b>, respectively.

=item aEnd bEnd

End of hit in I<a>, or I<b>, respectively.

=item aRemain

Number of bases in I<a> after the end of the hit.  ("0" if the hit
extends to the end of I<a>.)

=item bRemain

The equivalent for sequence I<b> of aRemain, but note that if the
strand of I<b> which matches is the reverse strand, then this is the
number of bases left in I<b> after the hit, back towards the beginning
of I<b>.

=item strand

The strand on sequence I<b> which was matched.  This is "B<W>" for the
forward (or parallel or B<Watson>) strand, and "B<C>" for the reverse
(or anti-parallel or B<Crick>) strand.

=item a2b

Returns coordinate in sequence I<b> equivalent to value in sequence
I<a> passed to method.  If no corresponding base, returns I<-1>.

=item fp

Returns an array of strings contining full list of ungapped alignment
fragment coordinates, filtered by score value passed to method.

=back

=item QUERYING BY SEQUENCE NAME

These methods allow access to the B<start>, B<end>, and B<remain>
fields for those occasions where you know the name of your sequence,
but don't necessarily know if it is sequence I<a> or I<b>.  Syntax and
behaviour is the same as for the I<DATA FIELDS> functions, but the
first argument to the function is the name of the sequence you want to
get data about.

I<Note:> These methods perform a filtering operation, reducing the
active list to only those hits which match the name given (in either
the I<a> or I<b> columns).  You'll therefore need to B<unfilter> your
match if you need data from hits which have now become hidden.

If the sequence name is found in both columns I<a> and I<b>, then a
warning is printed, and the I<a> column is chosen.

For example, suppose you have matched two fasta files containing only
the sequences cN75H12 and bK1109B5 (in that order), then the following
calls retrieve the same data:

	@ends = $match->aEnd();
	@ends = $match->end('cN75H12');

	$start = $match->bStart(0);
	$start = $match->start('bK1109B5', 0);

A warning is printed to STDERR if a match contains hits from the name
supplied in both columns, and only hits from the a column are
returned.

=over 4

=item start end remain

Access the aStart or bStart, aEnd or bEnd, and aRemain or bRemain
fields, depending upon wether the name supplied matches the aName or
bName fields respectively.

=back

=head1 BUGS


=head1 AUTHOR

B<James Gilbert> email jgrg@sanger.ac.uk

B<Tim Hubbard> email th@sanger.ac.uk (alignment processing extensions)

=cut





