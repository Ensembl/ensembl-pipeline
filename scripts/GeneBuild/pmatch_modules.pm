# holds a pmatch protein hit - simply the name(idenitifier) of the protein
# and a list of ContigHit objects

package ProteinHit;
use Bio::Root::Object;

@ISA = qw(Bio::Root::Object);

=head2 new

 Title   : new
 Usage   :
 Function:
 Example :
 Returns : 
 Args    : constructor


=cut

sub new {
  my ($class, @args) = @_;
  my $self = bless {}, $class;

  my ($id) = $self->_rearrange(['ID'], @args);

  $self->throw("No id") unless defined $id;
  $self->id($id);
  
  return $self;

}

=head2 id

 Title   : id
 Usage   :
 Function: get/set for contig id
 Example :
 Returns : 
 Args    :


=cut

sub id {
  my ($self,$id) = @_;
  if ($id) {
    $self->{'id'} = $id;
  }
  return $self->{'id'};
}

=head2 add_ContigHit

 Title   : add_ContigHit
 Usage   :
 Function: adds a ContigHit into $self->{_contig_hits}
 Example :
 Returns : 
 Args    :


=cut

sub add_ContigHit {
  my ($self,$hit) = @_;
  $self->throw('No contig hit') unless defined $hit;
  $self->throw('$contig is not a ContigHit') unless $hit->isa("ContigHit");
  $self->{_contig_hits}{$hit->id()} = $hit;
}

=head2 each_ContigHit

 Title   : each_ContigHit
 Usage   :
 Function: returns all entries in $self->{_contig_hits}
 Example :
 Returns : 
 Args    :


=cut

sub each_ContigHit {
  my ($self) = @_;
  return values %{$self->{_contig_hits}};
}


=head2 get_ContigHit

 Title   : get_ContigHit
 Usage   :
 Function: returns entries for a particular contig in $self->{_contig_hits}
 Example :
 Returns : 
 Args    :


=cut

sub get_ContigHit {
  my ($self,$contig) = @_;
  return ($self->{_contig_hits}{$contig}) if defined $contig;
  return undef;
}

1;

# holds a MergedHit - produced by munging together CoordinatePairs
# MergedHit knows which contig and protein it is pairing, strand, overall 
# coverage and details of component CoordinatePairs

package MergedHit;
use Bio::Root::Object;

@ISA = qw(Bio::Root::Object);

=head2 new

 Title   : new
 Usage   :
 Function: constructor
 Example :
 Returns : 
 Args    : 


=cut

sub new {
  my ($class, @args) = @_;
  my $self = bless {}, $class;

  my ($query, $target, $strand,$coverage) = $self->_rearrange(['QUERY',
							       'TARGET',
							       'STRAND',
							       'COVERAGE'],@args);

  $self->throw("No query") unless defined $query;
  $self->query($query);

  $self->throw("No target") unless defined $target;
  $self->target($target);

  $self->throw("No strand") unless defined $strand;
  $self->strand($strand);

  $self->throw("No coverage") unless defined $coverage;
  $self->coverage($coverage);

  $self->{'_coord_pairs'} = [];

  return $self;

}

=head2 query

 Title   : query
 Usage   :
 Function: get/set for query (fpccontig name)
 Example :
 Returns : 
 Args    : 


=cut

sub query {
  my ($self,$arg) = @_;
  if ($arg) {
    $self->{'query'} = $arg;
  }
  return $self->{'query'};
}

=head2 target

 Title   : target
 Usage   :
 Function: get/set for target (protein name)
 Example :
 Returns : 
 Args    : 


=cut

sub target {
  my ($self,$arg) = @_;
  if ($arg) {
    $self->{'target'} = $arg;
  }
  return $self->{'target'};
}

=head2 coverage

 Title   : coverage
 Usage   :
 Function: get/set for coverage (% of target covered by this hit
 Example :
 Returns : 
 Args    : 


=cut

sub coverage {
  my ($self,$arg) = @_;
  if ($arg) {
    $self->{'coverage'} = $arg;
  }
  return $self->{'coverage'};
}

=head2 strand

 Title   : strand
 Usage   :
 Function: get/set for (query) strand
 Example :
 Returns : 
 Args    : 


=cut

sub strand {
  my ($self,$arg) = @_;
  if ($arg) {
    $self->{'strand'} = $arg;
  }
  return $self->{'strand'};
}

=head2 add_CoordPair

 Title   : add_CoordPair
 Usage   :
 Function: adds a CoordPair to the supporting data for this MergedHit
 Example :
 Returns : 
 Args    : 


=cut


sub add_CoordPair {
  my ($self,$pair) = @_;
  $self->throw('No coord pair') unless defined $pair;
  $self->throw('$pair is not a CoordPair') unless $pair->isa("CoordPair");


  push(@{$self->{_coord_pairs}},$pair);
  
  # need to do some sorting?
  if ($self->strand == 1) {
    @{$self->{_coord_pairs}} = sort {$a->qstart <=> $b->qstart} @{$self->{_coord_pairs}};
  }
  else {
    @{$self->{_coord_pairs}} = sort {$b->qstart <=> $a->qstart} @{$self->{_coord_pairs}};
  }
  
  $self->tstart(@{$self->{_coord_pairs}}[0]->tstart);
  $self->qstart(@{$self->{_coord_pairs}}[0]->qstart);

  $self->tend(@{$self->{_coord_pairs}}[scalar(@{$self->{_coord_pairs}})-1]->tend);
  $self->qend(@{$self->{_coord_pairs}}[scalar(@{$self->{_coord_pairs}})-1]->qend);
  
}

=head2 subsume_MergedHit

 Title   : subsume_MergedHit
 Usage   :
 Function: Incorporate all the pairs from another MergedHit into this one
 Example :
 Returns : 
 Args    : 


=cut
sub subsume_MergedHit {
  my ($self,$hit) = @_;
  $self->throw('No hit') unless defined $hit;
  $self->throw('$hit is not a MergedHit') unless $hit->isa("MergedHit");

  push(@{$self->{_coord_pairs}},$hit->each_CoordPair);
  
  # need to do some sorting?
  if ($self->strand == 1) {
    @{$self->{_coord_pairs}} = sort {$a->qstart <=> $b->qstart} @{$self->{_coord_pairs}};
  }
  else {
    @{$self->{_coord_pairs}} = sort {$b->qstart <=> $a->qstart} @{$self->{_coord_pairs}};
  }
  
  $self->tstart(@{$self->{_coord_pairs}}[0]->tstart);
  $self->qstart(@{$self->{_coord_pairs}}[0]->qstart);

  $self->tend(@{$self->{_coord_pairs}}[scalar(@{$self->{_coord_pairs}})-1]->tend);
  $self->qend(@{$self->{_coord_pairs}}[scalar(@{$self->{_coord_pairs}})-1]->qend);
  
}

=head2 each_CoordPair

 Title   : each_CoordPair
 Usage   :
 Function: returns the CoordPairs contibuting to this MergedHit
 Example :
 Returns : 
 Args    : 


=cut

sub each_CoordPair {
  my ($self) = @_;
  # sorted by qstart, based on strand

  return @{$self->{_coord_pairs}};
}

=head2 tstart

 Title   : tstart
 Usage   :
 Function: returns overall start of hit in protein coords
 Example :
 Returns : 
 Args    : 


=cut

sub tstart {
  my ($self,$arg) = @_;
  if ($arg) {
    $self->{'tstart'} = $arg;
  }
  return $self->{'tstart'};
}

=head2 tend

 Title   : tend
 Usage   :
 Function: returns overall end of hit in protein coords
 Example :
 Returns : 
 Args    : 


=cut

sub tend {
  my ($self,$arg) = @_;
  if ($arg) {
    $self->{'tend'} = $arg;
  }
  return $self->{'tend'};
}

=head2 qstart

 Title   : qstart
 Usage   :
 Function: returns overall start of hit in fpc contig coords
 Example :
 Returns : 
 Args    : 


=cut

sub qstart {
  my ($self,$arg) = @_;
  if ($arg) {
    $self->{'qstart'} = $arg;
  }
  return $self->{'qstart'};
}

=head2 qend

 Title   : qend
 Usage   :
 Function: returns overall end of hit in fpc contig coords
 Example :
 Returns : 
 Args    : 


=cut

sub qend {
  my ($self,$arg) = @_;
  if ($arg) {
    $self->{'qend'} = $arg;
  }
  return $self->{'qend'};
}

1;

# holds a pmatch contig hit - simply the name(identifier) of the contig
# and a list of start-end positions

package ContigHit;
use Bio::Root::Object;

@ISA = qw(Bio::Root::Object);


=head2 new

 Title   : new
 Usage   :
 Function: constructor
 Example :
 Returns : 
 Args    : 


=cut

sub new {
  my ($class, @args) = @_;
  my $self = bless {}, $class;

  my ($id) = $self->_rearrange(['ID'], @args);

  $self->throw("No id") unless defined $id;
  $self->id($id);
 
  $self->{'_forward_pairs'} = [];
  $self->{'_reverse_pairs'} = [];

  return $self;

}

=head2 id

 Title   : id
 Usage   :
 Function: get/set for contig id
 Example :
 Returns : 
 Args    : 


=cut

sub id {
  my ($self,$id) = @_;
  if ($id) {
    $self->{'id'} = $id;
  }
  return $self->{'id'};
}

=head2 add_CoordPair

 Title   : add_CoordPair
 Usage   :
 Function: adds a CoordPair to the list making up this hit
 Example :
 Returns : 
 Args    : 


=cut

sub add_CoordPair {
  my ($self,$pair) = @_;
  $self->throw('No coord pair') unless defined $pair;
  $self->throw('$pair is not a CoordPair') unless $pair->isa("CoordPair");
  if($pair->strand == 1) {
    push(@{$self->{_forward_pairs}},$pair);
  }
  else {
    push(@{$self->{_reverse_pairs}},$pair);
  }
}

=head2 each_ForwardPair

 Title   : each_ForwardPair
 Usage   :
 Function: returns CoordPairs represeting hits between a prtein and the forward strand of the contig
 Example :
 Returns : 
 Args    : 


=cut

sub each_ForwardPair {
  my ($self) = @_;
  return @{$self->{_forward_pairs}};
}

=head2 each_ReversePair

 Title   : each_Reverseair
 Usage   :
 Function: returns CoordPairs representing hits between a protein and the reverse strand of the contig
 Example :
 Returns : 
 Args    : 


=cut

sub each_ReversePair {
  my ($self) = @_;
  return @{$self->{_reverse_pairs}};
}
1;

# holds a Coordinate Pair; name of the protein and contig that are hit (as a cross check) and 
# the start and end positions of the 

package CoordPair;
use Bio::Root::Object;

@ISA = qw(Bio::Root::Object);

=head2 new

 Title   : new
 Usage   :
 Function: constructor
 Example :
 Returns : 
 Args    : 


=cut

sub new {
  my ($class, @args) = @_;
  my $self = bless {}, $class;

  my ($query, $target, $qstart, $qend, $tstart, $tend, $percent, $strand) = $self->_rearrange(['QUERY',
										     'TARGET',
										     'QSTART',
										     'QEND',
										     'TSTART',
										     'TEND', 
										     'PERCENT',
										     'STRAND'],@args);

  $self->throw("No query") unless defined $query;
  $self->query($query);

  $self->throw("No target") unless defined $target;
  $self->target($target);

  $self->throw("No qstart") unless defined $qstart;
  $self->qstart($qstart);

  $self->throw("No qend") unless defined $qend;
  $self->qend($qend);

  $self->throw("No tstart") unless defined $tstart;
  $self->tstart($tstart);

  $self->throw("No tend") unless defined $tend;
  $self->tend($tend);

  $self->throw("No percent") unless defined $percent;
  $self->percent($percent);

  $self->throw("No strand") unless defined $strand;
  $self->strand($strand);

  return $self;

}

=head2 query

 Title   : query
 Usage   :
 Function: get/set for query (contig name)
 Example :
 Returns : 
 Args    : 


=cut

sub query {
  my ($self,$arg) = @_;
  if ($arg) {
    $self->{'query'} = $arg;
  }
  return $self->{'query'};
}

=head2 target

 Title   : target
 Usage   :
 Function: get/set for target (protein name)
 Example :
 Returns : 
 Args    : 


=cut

sub target {
  my ($self,$arg) = @_;
  if ($arg) {
    $self->{'target'} = $arg;
  }
  return $self->{'target'};
}

sub qstart {
  my ($self,$arg) = @_;
  if ($arg) {
    $self->{'qstart'} = $arg;
  }
  return $self->{'qstart'};
}

sub qend {
  my ($self,$arg) = @_;
  if ($arg) {
    $self->{'qend'} = $arg;
  }
  return $self->{'qend'};
}

sub tstart {
  my ($self,$arg) = @_;
  if ($arg) {
    $self->{'tstart'} = $arg;
  }
  return $self->{'tstart'};
}

sub tend {
  my ($self,$arg) = @_;
  if ($arg) {
    $self->{'tend'} = $arg;
  }
  return $self->{'tend'};
}

sub percent {
  my ($self,$arg) = @_;
  if ($arg) {
    $self->{'percent'} = $arg;
  }
  return $self->{'percent'};
}

sub strand {
  my ($self,$arg) = @_;
  if ($arg) {
    $self->{'strand'} = $arg;
  }
  return $self->{'strand'};
}



1;

# holds the object that runs pmatch and does the first round of filtering
package First_PMF;
use Bio::Root::Object;

@ISA = qw(Bio::Root::Object);

sub new {
  my ($class, @args) = @_;
  my $self = bless {}, $class;

  my ($plengths, $protfile, $fpcfile, $pmatch, $tmpdir, $maxintronlen) = 
                 $self->_rearrange(['PLENGTHS','PROTFILE','FPCFILE','PMATCH','TMPDIR','MAXINTRONLEN'], @args);
  $self->throw("No protlengths data") unless defined $plengths;
  $self->plengths($plengths);

  $self->throw("No fpcfile") unless defined $fpcfile;
  $self->fpcfile($fpcfile);

  $self->throw("No protfile") unless defined $protfile;
  $self->protfile($protfile);

  $self->throw("No pmatch executable") unless defined $pmatch;
  $self->pmatch($pmatch);

  $self->throw("No tmpdir") unless defined $tmpdir;
  $self->tmpdir($tmpdir);

  $self->throw("No maximum intron length") unless defined $maxintronlen;
  $self->maxintronlen($maxintronlen);

  my %proteins = ();
  $self->{_proteins} = \%proteins; 

  return $self;

}

=head2 make_coord_pair

 Title   : make_coord_pair
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub make_coord_pair {
  my ($self) = @_;
  my @cols = split();
  # clean up the refseq IDs
  my $protid = $cols[5];
  # alter this as necessary 
  if($protid =~ /\w+\|\w+\|\w+\|([\w\.]+)\|/) {
    $protid = $1;
  }
  
  # sort out strand hmmm if strand = -1 then shouldn't we switch qstart & qend? Currently don't ...
  my $strand = 1;
  if ( $cols[3] < $cols[2]  ) { $strand=-1; }
  
  my $cp = new CoordPair(
			 -query  => $cols[1],
			 -target => $protid,
			 -qstart => $cols[2],
			 -qend   => $cols[3],
			 -tstart => $cols[6],
			 -tend   => $cols[7],
			 -percent=> $cols[8],
			 -strand => $strand,
			);

  # where to put the CoordPair?
  # find the relevant ProteinHit, or make a new one
  my $proteins = $self->{_proteins};
  my $ph = $$proteins{$protid};

  if (!defined $ph) {
     #make a new one and add it into %proteins
    $ph = new ProteinHit(-id=>$protid);
    $$proteins{$protid} = $ph;
  }
  
  # now find the relevant ContigHit, or make a new one
  my $ch = $ph->get_ContigHit($cols[1]);
  if (!defined $ch) {
     # make a new one and add it into $ph
    $ch = new ContigHit(-id => $cols[1]);
    $ph->add_ContigHit($ch);
  }
  
  # now add the CoordPair
  $ch->add_CoordPair($cp);

}

=head2 merge_hits

 Title   : merge_hits
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub merge_hits {
  my ($self) = @_;
  my @merged;

  # merge the hits together.
  my $protref = $self->{_proteins};

  foreach my $p_hit(values %$protref) {
    my @allhits = $self->make_mergelist($p_hit);
    
    my @chosen = $self->prune(@allhits);
#    print STDERR "\nNo hits good enough for " . $p_hit->id() . "\n"
#      unless scalar(@chosen);
    push(@merged,@chosen);
  }

  return @merged;
  
}

=head2 make_mergelist

 Title   : make_mergelist
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub make_mergelist {
  my ($self, $ph) = @_;

  my @hits = ();
  foreach my $ch($ph->each_ContigHit()) {
    # forward strand
    my @cps = $ch->each_ForwardPair();
    # sort!
    @cps = sort {$a->qstart <=> $b->qstart} @cps;
    #reverse strand
    my @mps = $ch->each_ReversePair();
    @mps = sort {$b->qstart <=> $a->qstart} @mps;
    push (@cps,@mps);
    
    # deal with forward & reverse separately?
    my $first = shift(@cps);
    my $mh = $self->new_merged_hit($first);
    
    push(@hits,$mh);

    my $prev = $hits[$#hits];
    my $prev_cp = $first;

    CP: 
    foreach my $cp(@cps) {

      # need to compare with the last entry in @merged

      # first the strand
      my $strand = 1;
      if ($cp->qend() < $cp->qstart()) { $strand = -1; }

      # does this CoordPair extend the current hit?
      # need a fudge factor - pmatch could overlap them by 1 or 2 ... or 3
      if( $strand == $prev->strand &&
	 ( ($cp->tstart >=  $prev_cp->tend) ||
	 (($prev_cp->tend -  $cp->tstart) <= 3)) &&
	  # Steve's fix - Added qstart/qend condition (*strand should allow it to work on
	  #     either strand)
	  # no overlap currently allowed
          $cp->qstart*$strand >= $prev_cp->qend*$strand 
	)
	{
	  #extend existing MergedHit
	  my $coverage = $cp->tend - $cp->tstart + 1;
	  $coverage += $prev->coverage();

	  # compensate for any overlap 
	  my $diff = $prev_cp->tend -  $cp->tstart;
	  if($diff >=0) {
	    $diff++;
	    $coverage -= $diff;
	  }

	  $prev->coverage($coverage);
	  $prev->add_CoordPair($cp);
          $prev_cp = $cp;
	}
      else {
	# make a new MergedHit
	my $mh = $self->new_merged_hit($cp);
	push(@hits,$mh);	

        $prev = $hits[$#hits];
        $prev_cp = $cp;
      }
    }
  }

  # my $times = join " ", times;
  # print "Before extend " . $times . "\n";
  # extend those merged hits

  # print "Number of hits = " . scalar(@hits) . "\n";

  @hits = $self->extend_hits(@hits);

  # $times = join " ", times;
  # print "After extend " . $times . "\n";
  # print "Number of hits = " . scalar(@hits) . "\n";

  # sort out coverage
  my $plengths = $self->plengths;
  my $protlen = $$plengths{$ph->id};
  warn "No length for " . $ph->id . "\n" unless $protlen;
  return unless $protlen; 
  foreach my $hit(@hits) {
    my $percent = $hit->coverage;
    $percent *= 100;
    $percent /= $protlen;
    $percent=sprintf("%.1f",$percent);
    $hit->coverage($percent);
  }

  return @hits;
}

=head2 new_merged_hit

 Title   : new_merged_hit
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub new_merged_hit {
  my ($self,$cp) = @_;
  my $coverage = $cp->tend - $cp->tstart + 1;
  my $mh = new MergedHit( -query    =>  $cp->query(),
			  -target   =>  $cp->target(),
			  -strand   =>  $cp->strand(),
			  -coverage =>  $coverage,
			);
  $mh->add_CoordPair($cp);
  return $mh;
}

sub plengths {
  my ($self, $plengths) = @_;
  
  # $plengths is a hash reference
  if(defined($plengths)){
    if (ref($plengths) eq "HASH") {
      $self->{_plengths} = $plengths;
    } 
    else {
      $self->throw("[$plengths] is not a hash ref.");
    }
  }
     
  return $self->{_plengths};

}

sub tmpdir {
  my ($self, $tmpdir) = @_;
  
  if(defined($tmpdir)){
    $self->throw("Cannot open [$tmpdir]") unless opendir(DIR, $tmpdir);
    closedir(DIR);
    $self->{_tmpdir} = $tmpdir;
  } 
  
  return $self->{_tmpdir};

}

sub protfile {
  my ($self, $protfile) = @_;
  
  if(defined($protfile)){
    if (-e $protfile) {
      $self->{_protfile} = $protfile;
    } 
    else {
      $self->throw("[$protfile] does not exist.");
    }
  }
  return $self->{_protfile};
  
}

sub maxintronlen {
  my ($self, $maxintronlen) = @_;
  
  if(defined($maxintronlen)){
    if ($maxintronlen =~ /\d+/) {
      $self->{_maxintronlen} = $maxintronlen;
    } 
    else {
      $self->throw("[$maxintronlen] is not numeric.");
    }
  }
  return $self->{_maxintronlen};
  
}

sub pmatch {
  my ($self, $pmatch) = @_;
  
  if(defined($pmatch)){
    if (-e $pmatch) {
      $self->{_pmatch} = $pmatch;
    } 
    else {
      $self->throw("[$pmatch] does not exist.");
    }
  }
  return $self->{_pmatch};
  
}

sub fpcfile {
  my ($self, $fpcfile) = @_;
  
  if(defined($fpcfile)){
    if (-e $fpcfile) {
      $self->{_fpcfile} = $fpcfile;
    } 
    else {
      $self->throw("[$fpcfile] does not exist.");
    }
  }  
  return $self->{_fpcfile};
  
}

sub run {
  my ($self) = @_;
  
  my $protfile = $self->protfile;
  my $fpcfile  = $self->fpcfile;
  my $pmatch   = $self->pmatch;
  print STDERR "$pmatch $protfile $fpcfile\n";
  
  # run pmatch  
  open (PM, "$pmatch -D $protfile $fpcfile |" );

  # parse results
  while(<PM>) {
    next unless /\S+/;
    $self->make_coord_pair($_);
  }

  close PM;

  # merge hits
  my @hits = $self->merge_hits;

  # store results
  if (!defined($self->{_output})) {
    $self->{_output} = [];
  } 
  
  push(@{$self->{_output}},@hits);
}

sub output{
  my ($self) = @_;
  if (!defined($self->{_output})) {
    $self->{_output} = [];
  } 
  return @{$self->{_output}};
}

sub extend_hits {
  my ($self, @hits) = @_;

  # we want to do essentially what we did to create the merged hits but we can skip over
  # intervening hits if we need to.
  my @fhits;
  my @rhits;
  my @newhits;

  foreach my $mh(@hits){
    if($mh->strand == 1){
      push (@fhits, $mh);
    }
    else {
      push (@rhits, $mh);
    }
  }

  @fhits = sort {$a->qstart <=> $b->qstart} @fhits;
  @rhits = sort {$b->qstart <=> $a->qstart} @rhits;


  while(scalar(@fhits)){
    my $hit = shift(@fhits);
    push (@newhits, $hit);

    # can we link it up to a subsequent hit?
    foreach my $sh(@fhits){
      die ("Argh!") unless $hit->strand == $sh->strand;
      last if ($hit->qend+$self->{_maxintronlen} < $sh->qstart);
      if ($sh->tstart > $hit->tend &&
	  $sh->tend   > $hit->tend &&
	  abs($sh->tstart - $hit->tend) <= 3 &&
	  # qstart/qend condition - no overlap currently allowed
	  $sh->qstart >= $hit->qend
	 ) {
	# add the coord pairs from $sh into $hit
        $hit->subsume_MergedHit($sh);
      }
    }
  }

  
# same for rhits
  while(scalar(@rhits)){
    my $hit = shift(@rhits);
    push (@newhits, $hit);
    

    # can we link it up to a subsequent hit?
    foreach my $sh(@rhits){
      die ("Argh!") unless $hit->strand == $sh->strand;
# hmmmm On minus strand qstart is currently > qend
# ie $sh->qend <= $sh->qstart <= $hit->qend <= $hit->qstart
#      last if ($hit->qstart-$self->{_maxintronlen} > $sh->qend);
      last if ($hit->qend-$self->{_maxintronlen} > $sh->qstart);
      if ($sh->tstart > $hit->tend &&
	  $sh->tend   > $hit->tend &&
	  abs($sh->tstart - $hit->tend) <= 3 &&
	  # qstart/qend condition - no overlap currently allowed
	  $hit->qend >= $sh->qstart
	 ) {
	# add the coord pairs from $sh into $hit
        $hit->subsume_MergedHit($sh);
      }
    }
  }

  # return extended hits
  return (@newhits);  
}

=head2 prune

 Title   : prune
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub prune {
  my ($self, @all) = @_;
  my @chosen = ();
  # reject hits with < 25% coverage
  my $lower_threshold = 25;

  # sort by descending order of coverage
  @all = sort {$b->coverage <=> $a->coverage} @all;

  my $first = shift(@all);
  if($first->coverage < $lower_threshold){
    # we're done
    return @chosen;
  }
  
  push (@chosen,$first);

  # don't select any hits that have coverage less than 2% below that of the first hit, be it 100 or 99.9 or ...
  my $curr_pc = $first->coverage() - 2;

 PRUNE:
  foreach my $hit(@all) {
    last PRUNE if $hit->coverage < $lower_threshold;
    last PRUNE if $hit->coverage < $curr_pc;
    push (@chosen,$hit);
  }

  return @chosen;

}

1;

# holds the object that does the second round of filtering
package Second_PMF;
use Bio::Root::Object;

@ISA = qw(Bio::Root::Object);

sub new {
  my ($class, @args) = @_;
  my $self = bless {}, $class;

  my ($phits) = $self->_rearrange(['PHITS'], @args);

  $self->throw("No pmatch hits data") unless defined $phits;
  $self->phits($phits);

  my %proteins = (); # hash of arrays of MergedHits, indexed by protin name
  $self->{_proteins} = \%proteins;

  return $self;

}

sub phits {
  my ($self, $phits) = @_;
  
  # $phits is an array reference
  if(defined($phits)){
    if (ref($phits) eq "ARRAY") {
      $self->{_phits} = $phits;
    } 
    else {
      $self->throw("[$phits] is not an array ref.");
    }
  }
     
  return $self->{_phits};

}

sub run {
  my ($self) = @_;

  # group hits by protein

  my %prots = %{$self->{_proteins}}; # just makes it a bit easier to follow

  foreach my $hit(@{$self->phits}){
    # print the details of all the constituent coord pairs separated by white space 
    push (@{$prots{$hit->target}}, $hit);
  }

  $self->{_proteins} = \%prots;

  # prune and store the hits
  $self->prune_hits;

}

sub prune_hits {
  my ($self) = @_;  
  my %prots = %{$self->{_proteins}}; # just makes it a bit easier to follow  

  PROTEIN:
  foreach my $p(keys %prots){
    my @chosen = ();
    my @allhits = @{$prots{$p}};
    
    # sort by descending order of coverage
    @allhits = sort {$b->coverage <=> $a->coverage} @allhits;
    
    my $first = shift(@allhits);
    
    # don't select any hits that have coverage less than 2% below that of the first hit, be it 100 or 99.9 or ...
    my $curr_pc = $first->coverage() - 2; 
    
    # lower bound threshold - reject anything with < 25% coverage
    my $lower_threshold = 25;
    next PROTEIN if $first->coverage < $lower_threshold;
    
    push (@chosen,$first) unless $first->coverage < $lower_threshold;
  PRUNE:
    foreach my $hit(@allhits) {
      
      last PRUNE if $hit->coverage() < $curr_pc;
      last PRUNE if $hit->coverage() < $lower_threshold;
      push (@chosen,$hit);
    }
    
    push(@{$self->{_output}},@chosen);
  }
  
}

sub output {
  my ($self) = @_;
  if (!defined($self->{_output})) {
    $self->{_output} = [];
  } 
  return @{$self->{_output}};
}

1;
