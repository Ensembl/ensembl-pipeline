# Author: Marc Sohrmann (ms2@sanger.ac.uk)
#
# You may distribute this module under the same terms as perl itself

=pod

=head1 NAME

  Bio::EnsEMBL::Pipeline::Runnable::Waba

=head1 SYNOPSIS

  my $seqstream = Bio::SeqIO->new ( -file => $clonefile,
                                    -fmt => 'Fasta',
                                  );
  $seq = $seqstream->next_seq;

  my $waba = Bio::EnsEMBL::Pipeline::Runnable::Waba->new ( -QUERY => $seq);
  $waba->workdir ($workdir);
  $waba->run;
  my @results = $waba->output;

=head1 DESCRIPTION

  Waba takes a Bio::Seq (or Bio::PrimarySeq) object
  and runs waba against a fasta database (Jim Kents aligner,
  initially written for elegans-briggsae comparison)
  The resulting output file is parsed to produce a set of features.

=head1 CONTACT

  B<ensembl-dev@ebi.ac.uk>

=head1 APPENDIX

  The rest of the documentation details each of the object methods.
  Internal methods are usually preceded with a _.

=cut

package Bio::EnsEMBL::Pipeline::Runnable::Waba;

use vars qw(@ISA);
use strict;
$| = 1;

use Bio::EnsEMBL::Root;
use Bio::EnsEMBL::Pipeline::RunnableI;
use Bio::EnsEMBL::SeqFeature;
use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::DnaDnaAlignFeature;
use Bio::SearchIO;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI);


=head2 new

    Arg [1]    : Bio::Seq $query
    Arg [2]    : char $waba_binary
    Arg [3]    : Bio::EnsEMBL::Analysis $analysis
    Description: initialises Waba object
    Returns    : a Waba object
    Exceptions : none
    Caller     : general

=cut

sub new {
    my ($class, @args) = @_;

    my $self = $class->SUPER::new (@_);

    $self->{'_fplist'} = [];              # an array of Bio::SeqFeatures
    $self->{'_sequence'}  = undef;        # location of Bio::Seq object
    $self->{'_waba'}      = undef;        # location of waba executable
    $self->{'_database'}  = undef;        # name of database
    $self->{'_workdir'}   = undef;        # location of tmp directory
    $self->{'_filename'}  = undef;        # file to store Bio::Seq object
    $self->{'_results'}   = undef;        # file to store results of waba run
    $self->{'_protected'} = [];           # a list of files protected from deletion

    my ($query, $waba, $database) = $self->_rearrange([qw(
	QUERY
        WABA
	DATABASE
    )], @args);

    $waba = $args[1]->program_file if ($args[1]->program_file);
    undef $waba unless (-e $waba);

    $database = $args[1]->db_file;

    $waba ||= 'waba';
    $self->query   ($query)    if ($query);
    $self->database($database) if ($database);
    $self->waba($self->find_executable($waba));

    return $self;
}


=head2 query

    Arg [1]    : Bio::SeqI $query
    Description: accessor for query sequence
    Returntype : Bio::SeqI
    Exceptions : query not a Bio::PrimarySeqI or Bio::SeqI
    Caller     : general

=cut

sub query {
    my ($self, $seq) = @_;
    if ($seq) {
	($seq->isa ("Bio::PrimarySeqI") || $seq->isa ("Bio::SeqI"))
	    || $self->throw("Input isn't a Bio::SeqI or Bio::PrimarySeqI");
	$self->{'_sequence'} = $seq;
	$self->filename ($self->query->id.".$$.seq");
	$self->results ($self->filename.".out");
	$self->file ($self->results);
    }
    return $self->{'_sequence'};
}


=head2 waba

    Arg [1]    : string $waba
    Description: accessor for waba binary
    Returntype : string
    Exceptions : binary not found and executable
    Caller     : general

=cut

sub waba {
    my ($self, $location) = @_;

    if ($location) {
        unless (-x $location) {
            $self->throw ("program not found at $location");
	}
        $self->{'_waba'} = $location;
    }
    return $self->{'_waba'};
}


=head2 database

    Arg [1]    : string $db
    Description: accessor for target database
    Returntype : string
    Exceptions : none
    Caller     : general

=cut

sub database {
    my $self = shift;
    if (@_) {
        $self->{'_database'} = shift;
    }
    return $self->{'_database'};
}


=head2 run

    Args       : none
    Description: runs the waba program and creates set of features
    Returntype : none
    Exceptions : none
    Caller     : general

=cut

sub run {
    my ($self, $dir) = @_;

    # check clone
    my $seq = $self->query || $self->throw("Clone required for Program\n");

    # set directory if provided
    $self->workdir ('/tmp') unless ($self->workdir($dir));
    $self->checkdir;

    # reset filename and results as necessary (adding the directory path)
    my $tmp = $self->workdir;
    my $input = $tmp."/".$self->filename;
    $self->filename ($input);
    $tmp .= "/".$self->results;
    $self->results ($tmp);

    # write sequence to file
    $self->writefile;

    # run program
    $self->run_program;

    # parse output
    $self->parse_results;
    $self->deletefiles;
}


=head2 run_program

    Args       : none
    Description: makes the system call to program
    Returntype : none
    Exceptions : system call unsuccessful
    Caller     : general

=cut

sub run_program {
    my ($self) = @_;

    # write the files containing the subject file lists

    my $path = $self->filename;
    open (QP, ">$path.query") || $self->throw("CANNOT create small program input file");
    print QP "$path\n";
    close QP;

    # write the files containing the query file lists
    # check database is a fasta file (first char '>') or a list of files

    my $db = $self->database;
    open DB, "< $db";
    my $char = getc DB;
    close DB;

    if ($char eq '>') {
	$db = $self->workdir . "/waba.$$.db";
        open DB, "> $db" || $self->throw("cannot create large program input file");
	print DB $self->database, "\n";
        close DB;
    }

    my $cmd = $self->waba . " all $path.query $db " . $self->results .
     "> /dev/null";

    $self->throw("Error running ".$self->waba." on ".$self->filename)
        unless ((system($cmd)) == 0);

    unlink "$path.query";
    unlink "$db" if $db ne $self->database;
}

=head2 output

    Arg [none] :
    Description: returns output of running dust
    Returntype : @{Bio::EnsEMBL::RepeatFeature}
    Exceptions : none
    Caller     : general

=cut

sub output {
    my ($self) = @_;

    return @{$self->{'_flist'}};
}



=head2 parse_results

 Title    :  parse_results
 Usage    :  $self->parse_results ($filename)
 Function :  parses program output to give a set of features
 Example  :
 Returns  : 
 Args     : filename (optional, can be filename, filehandle or pipe, not implemented)
 Throws   :

=cut

sub parse_results {
    my $self = shift;
    my $passed_results = shift;
    $self->results($passed_results) if $passed_results;
    my $filehandle;
    my $resfile = $self->results;
    
    if (-e $resfile) {
        # it's a filename
        if (-z $self->results) {  
	    print STDERR $self->program." didn't find anything\n";
	    return;
        }       
        else {
            open (OUT, "<$resfile") or $self->throw ("Error opening $resfile");
            $filehandle = \*OUT;
      }
    }
    else {
        # it'a a filehandle
        $filehandle = $resfile;
    }
    
    # parse

    my $switch = 0;
    my ($contig_match, $percent_id);
    my ($contig, $start, $end, $strand, $contig_file);
    my ($hid, $hstart, $hend);
    my ($seq, $hseq, $statepath);
    while (<$filehandle>) {
        chomp;
        if (/^(\S+)\s+align\s+(\S+)%\s+of\s+\d+\s+(\S+)\s+(\S+)\:(\d+)\-(\d+)\s+([+\-])\s+(\S+)\:(\d+)\-(\d+)/) {
            $contig_match = $1;
            $percent_id = $2;
            $contig_file = $3;
            $contig = $4;
            $start = $5;
            $end = $6;
            $strand = $7; $strand = ($strand eq "+") ? 1 : -1;
            $hid = $8;
            $hstart = $9;
            $hend = $10; 
            $switch = 1;
            # get rid of the path in the target name
            ($hid) = $hid =~ /\.([^\.]+)$/;
        }
        elsif ($switch == 1) {
            $seq = $_;
            $switch++; 
        }
        elsif ($switch == 2) {
            $hseq = $_;            
            $switch++;
        }
        elsif ($switch == 3  || eof) {
            $statepath = $_;
            $switch = 0;
            my @features = ();

            my @chars = split (//, $seq);
            my @hchars = split (//, $hseq);
            my @states = split (//, $statepath);
   
            # there is something weird going on with the coordinates in waba:
            # the start coordinate of both query and target have to be incremented
            # by 1 in order to correspond to the sequence coordinates (the end
            # coordinates seem fine, however).

            # query sequence can be sense or antisense
            my $substart = my $subend = ($strand == 1) ? $start+1 : $end;
            # target sequence always sense
            my $hsubstart = my $hsubend = $hstart+1;
            # set the increment value (depending on the strand information)
            my $inc = ($strand == 1) ? 1 : -1;
            # get the state of the first position
            my $current_state;
            if ($states[0] =~ /^[123]$/) {
                $current_state = "coding";
            }
            elsif ($states[0] =~ /^H$/) {
                $current_state = "high";
            }
            elsif ($states[0] =~ /^L$/) {
                $current_state = "low";
    	    }
            elsif ($states[0] =~ /^Q$/) {
                $current_state = "query";
	    }
            elsif ($states[0] =~ /^T$/) {
                $current_state = "target";
	    } 
            # loop over the remaining positions
            my @query_match;
            push (@query_match, $chars[0]);
            my @target_match;
            push (@target_match, $hchars[0]);
            my $count = 1;
            my $total_match = 0;
            my $total_mismatch = 0;
            my $total_block_score = 0;
            my $total_query_insert = 0;
            my $total_target_insert = 0;
            while ($count <= @chars) {
                if ($states[$count] =~ /^[123]$/ && $current_state eq "coding") {
                    $subend = ($chars[$count] ne '-') ? $subend+$inc : $subend;
                    $hsubend = ($hchars[$count] ne '-') ? $hsubend+1 : $hsubend;
                    push (@query_match, $chars[$count]);
                    push (@target_match, $hchars[$count]);
   	        }
                elsif ($states[$count] =~ /^H$/ && $current_state eq "high") {
                    $subend = ($chars[$count] ne '-') ? $subend+$inc : $subend;
                    $hsubend = ($hchars[$count] ne '-') ? $hsubend+1 : $hsubend;
                    push (@query_match, $chars[$count]);
                    push (@target_match, $hchars[$count]);
                }
                elsif ($states[$count] =~ /^L$/ && $current_state eq "low") {
                    $subend = ($chars[$count] ne '-') ? $subend+$inc : $subend;
                    $hsubend = ($hchars[$count] ne '-') ? $hsubend+1 : $hsubend;
                    push (@query_match, $chars[$count]);
                    push (@target_match, $hchars[$count]);
                }    
                elsif ($states[$count] =~ /^Q$/ && $current_state eq "query") {
                    $hsubend = ($hchars[$count] ne '-') ? $hsubend+1 : $hsubend;
                    push (@query_match, $chars[$count]);
                    push (@target_match, $hchars[$count]);
                }
                elsif ($states[$count] =~ /^T$/ && $current_state eq "target") {
                    $subend = ($chars[$count] ne '-') ? $subend+$inc : $subend;
                    push (@query_match, $chars[$count]);
                    push (@target_match, $hchars[$count]);
                }
                else {
                    # query gap
                    if ($substart == $subend) {
                        my $query_insert = $hsubend-$hsubstart+1;
                        # adjust $total_block_score;
                        $total_block_score += 8+4*log($query_insert);
                        $total_query_insert += $query_insert;
                    }
                    # target gap
                    elsif ($hsubstart == $hsubend) {
                        my $target_insert = ($subend > $substart) ? $subend-$substart+1 : $substart-$subend+1; 
                        # adjust $total_block_score;
                        $total_block_score += 8+4*log($target_insert);
                        $total_target_insert += $target_insert;
    		    }
                    # coding or high/low conservation
                    else {
                        # parse the subalignment
                        my ($align_coordinates, $cigar_string, $match, $mismatch, $qNumInsert, $qBaseInsert, $tNumInsert, $tBaseInsert, $block_score) = split_match ($substart, $hsubstart, \@query_match, \@target_match, $inc);
                        # increment "total" values
                        $total_match += $match;
                        $total_mismatch += $mismatch;
                        $total_block_score += $block_score;
                        $total_query_insert += $qBaseInsert;
                        $total_target_insert += $tBaseInsert;
                        # calculate score and percent identity
                        my $score = sprintf("%.2f", 2*($match-$mismatch)-$block_score);
                        my $pid = int(100*($match/($match+$mismatch+$qBaseInsert+$tBaseInsert)));
                        # adjust coordinates depending on strand 
			my $tmp_substart;
			my $tmp_subend;
                        if ($strand == -1){
                            my $tmp = $substart;
                            $tmp_substart = $subend;
                            $tmp_subend = $tmp;
		        }
			else {
			    $tmp_substart = $substart;
			    $tmp_subend = $subend;
			}
                        # build the featurepair
                        my %feature;
                        #$feature{name} = $self->clone->id;
                        $feature{name} = $self->clone_name;
                        $feature{percent_id} = $pid;
                        $feature{start} = $tmp_substart;
                        $feature{end} = $tmp_subend;
                        $feature{score} = $score;
                        $feature{strand} = $strand;
                        $feature{hname} = $hid;
                        $feature{hstart} = $hsubstart;
                        $feature{hend} = $hsubend;
                        $feature{hstrand} = 1;
                        ($feature{program}) = "waba";#$self->program =~ /([^\/]+)$/;
                        ($feature{db}) = $self->database =~ /([^\/]+)$/;
                        ($feature{logic_name}) = "waba";#$self->program =~ /([^\/]+)$/;
                        $feature{state} = $current_state;
                        $feature{cigar} = $cigar_string;
			$feature{p_value} = 0;#not actually used - req for seqfeature
                        my $fp = $self->create_feature (\%feature);
                        push (@features, $fp);
                        if ($count == @chars) {
                            last;
		        }
	       	    }
                    @query_match = $chars[$count];
                    @target_match = $hchars[$count];
                    if ($states[$count] =~ /^[123]$/) {
                        $current_state = "coding";
         	    } 
                    elsif ($states[$count] =~ /^H$/) {
                        $current_state = "high";
		    }
                    elsif ($states[$count] =~ /^L$/) {
                        $current_state = "low";
	  	    }
                    elsif ($states[$count] =~ /^Q$/) {
                        $current_state = "query";
		    }
                    elsif ($states[$count] =~ /^T$/) {
                        $current_state = "target";
		    }
                    if ($chars[$count] ne '-') {$substart = $subend+$inc;}
                    else {$substart = $subend;}

                    if ($hchars[$count] ne '-') {$hsubstart = $hsubend+1;}
                    else {$hsubstart = $hsubend;}
                    $subend = $substart;
                    $hsubend = $hsubstart;
   	        }
    	        $count++;
	    }
            my $total_score = sprintf ("%.2f", 2*($total_match-$total_mismatch)-$total_block_score);
            my $total_pid = int(100*($total_match/($total_match+$total_mismatch+$total_query_insert+$total_target_insert)));
            my %feature;
            #$feature{name} = $self->clone->id;
            $feature{name} = $self->clone_name;
            $feature{score} = $total_score;
            $feature{percent_id} = $total_pid;
            $feature{start} = $start+1;
            $feature{end} = $end;
            $feature{strand} = $strand;
            $feature{hname} = $hid;
            $feature{hstart} = $hstart+1;
            $feature{hend} = $hend;
            $feature{hstrand} = 1;
            ($feature{program}) = "waba";#$self->program =~ /([^\/]+)$/;
            ($feature{db}) = $self->database =~ /([^\/]+)$/;
            ($feature{logic_name}) = "waba";#$self->program =~ /([^\/]+)$/;
	    $feature{p_value} = 0;#not actually used - req for seqfeature
            my $fp = $self->create_feature (\%feature);
            unshift (@features, $fp);
            push (@{$self->{'_flist'}}, \@features);
        }
    }
    close $filehandle;   
}






=head2 create_feature

 Title    : create_feature
 Usage    : $self->create_feature ($feature)
 Function : creates a Bio::EnsEMBL::FeaturePair object from %feature,
            and pushes it onto @{$self->{'_flist'}}
 Example  :
 Returns  :
 Args     :
 Throws   :

=cut

sub create_feature {
    my ($self, $feat) = @_;

    # create analysis object (will end up in the analysisprocess table)
    my $analysis = Bio::EnsEMBL::Analysis->new ();
    $analysis->db ($feat->{db});
    $analysis->program ($feat->{program});
    $analysis->gff_source ($feat->{source});
    $analysis->gff_feature ($feat->{primary});
    $analysis->logic_name ($feat->{logic_name});

    # create featurepair object
    my $feature1 = Bio::EnsEMBL::SeqFeature->new ();
    $feature1->seqname ($feat->{name});
    $feature1->start ($feat->{start});
    $feature1->end ($feat->{end});
    $feature1->strand ($feat->{strand});
    $feature1->score ($feat->{score});
    $feature1->p_value ($feat->{p_value});
    $feature1->percent_id ($feat->{percent_id});
    $feature1->source_tag ($feat->{source});
    $feature1->primary_tag ($feat->{primary});
    $feature1->analysis ($analysis);

    if ($feat->{cigar}) {
        $feature1->add_tag_value ('cigar', $feat->{cigar});
    }
    if ($feat->{state}) {
        $feature1->add_tag_value ('state', $feat->{state});
    }

    my $feature2 = Bio::EnsEMBL::SeqFeature->new ();
    $feature2->seqname ($feat->{hname});
    $feature2->start ($feat->{hstart});
    $feature2->end ($feat->{hend});
    $feature2->strand ($feat->{hstrand});
    $feature2->score ($feat->{score});
    $feature2->p_value ($feat->{p_value});
    $feature2->percent_id ($feat->{percent_id});
    $feature2->source_tag ($feat->{source});
    $feature2->primary_tag ($feat->{primary});
    $feature2->analysis ($analysis);

    my $featurepair = Bio::EnsEMBL::FeaturePair->new ();
    $featurepair->feature1 ($feature1);
    $featurepair->feature2 ($feature2);

    if ($featurepair) {
      $featurepair->feature1->validate_prot_feature; 

      if ($feat->{cigar}) {
        $featurepair->feature1->add_tag_value ('cigar', $feat->{cigar});
      }
      if ($feat->{state}) {
        $featurepair->feature1->add_tag_value ('state', $feat->{state});
      }
      $featurepair->feature2->validate_prot_feature(2);
      return $featurepair;
    }
    return undef;
}

sub split_match {
    my ($qstart, $hstart, $qchars_ref, $hchars_ref, $inc) = @_;
    my @qchars = @$qchars_ref;
    my @hchars = @$hchars_ref;
    my @align_coordinates = ();
    my $cigar_string = "";

    # define some variable used to calculate the score
    my $match = 0;
    my $mismatch = 0;
    my $qNumInsert = 0;
    my $qBaseInsert = 0;
    my $tNumInsert = 0;
    my $tBaseInsert = 0;
    my $block_score = 0;
    my $insert = 0;

    # set the ungapped subalignment end
    my $qend   = $qstart;
    my $hend   = $hstart;

    # counter for the bases in the alignment    
    my $count = 0;
    # flag for ungapped subalignment
    my $found = 0;
    
    # loop over all positions
    while ($count <= $#qchars) {
        # We have hit an ungapped region. Increase the query and hit counters.
        # and set the subalignment flag
        if ($qchars[$count] ne '-' && $hchars[$count] ne '-') {
            $qend += $inc;
            $hend += 1;
            $found = 1;
            if ($qchars[$count] eq $hchars[$count]) {
                $match++;
	    }
            else {
                $mismatch++;
	    }
            if ($insert > 0) {
                $block_score += 8+4*log($insert);
                $insert = 0;
	    }
        # We have hit a gapped region
        } else {
            # If the subalignment flag is set,
            # push the coordinates on the array, and reset the start and end variables
            if ($found == 1) {
                push (@align_coordinates, [$qstart, $hstart]);
                if ($qchars[$count] eq '-') {
                    $qNumInsert++;
		}
                if ($hchars[$count] eq '-') {
                    $tNumInsert++;
		}
            }
            $insert++;
            # We're in a gapped region. We need to increment the sequence that
            # doesn't have the gap in it to keep the coordinates correct.
            # We also need to reset the current end coordinates.
            if ($qchars[$count] ne '-') {
                $qstart = $qend+$inc;
            } else {
                $qstart = $qend;
                $qBaseInsert++;
            }
            if ($hchars[$count] ne '-') {
                $hstart = $hend+1;
            } else {
                $hstart = $hend;
                $tBaseInsert++;
            }
            $qend = $qstart;
            $hend = $hstart;
            $found = 0;
        }
        $count++;
    }                        
    # Remember the last feature
    if ($found == 1) {
        push (@align_coordinates, [$qstart, $hstart]);
    }
    # make the cigar string
    my $last = $#align_coordinates;
    if ($last >= 0) {
        for (my $i = 0 ; $i < $last ; $i++) {
            my $ql = ($align_coordinates[$i+1]->[0])-($align_coordinates[$i]->[0]);
            $ql = ($ql > 0) ? $ql : -$ql;
            my $tl = ($align_coordinates[$i+1]->[1])-($align_coordinates[$i]->[1]);
            $tl = ($tl > 0) ? $tl : -$tl;
            my $length = ($ql > $tl) ? $ql-$tl : $tl-$ql;
            $cigar_string .= $align_coordinates[$i]->[0].",".$align_coordinates[$i]->[1].",$length:";
        }
        # add the final block
        $cigar_string .= $align_coordinates[$#align_coordinates]->[0].",".$align_coordinates[$#align_coordinates]->[1].",0";
    }

    return (\@align_coordinates, $cigar_string, $match, $mismatch, $qNumInsert, $qBaseInsert, $tNumInsert, $tBaseInsert, $block_score);
}


sub clone_name {
  my $self = shift;
  return $self->{_sequence}->display_id;
}



1;
