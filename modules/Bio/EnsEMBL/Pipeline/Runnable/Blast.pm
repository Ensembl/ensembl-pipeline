# Cared for by Michele Clamp  <michele@sanger.ac.uk>
#
# Copyright Michele Clamp
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::Runnable::Blast

=head1 SYNOPSIS

  # To run a blast job from scratch do the following.

  my $query = new Bio::Seq(-file   => 'somefile.fa',
                           -format => 'fasta');

  my $blast =  Bio::EnsEMBL::Pipeline::Runnable::Blast->new 
    ('-query'     => $query,
     '-program'   => 'blastp',
     '-database'  => 'swir',
     '-threshold' => 1e-6,
     '-filter'    => $filter,
     '-options'   => 'V=1000000');

  $blast->run();

  @featurepairs = $blast->output();

  foreach my $fp (@featurepairs) {
      print $fp->gffstring . "\n";
  }

  # Additionally if you have blast runs lying around that need parsing
  # you can do

  open(BLAST,"<blast.out");

  my @featurepairs = Bio::EnsEMBL::Pipeline::Runnable::Blast->parse_results(\*BLAST);

  close(BLAST);

=head1 DESCRIPTION

Blast takes a Bio::Seq (or Bio::PrimarySeq) object and runs blast with, the
output is parsed by BPLite and stored as Bio::EnsEMBL::FeaturePairs. 

The output features can be filtered by probability using the -threshold option.
Other options can be passed to the blast program using the -options method

=head1 CONTACT

Describe contact details here

=head1 APPENDIX


=cut

package Bio::EnsEMBL::Pipeline::Runnable::Blast;


use vars qw(@ISA);
use strict;
# Object preamble

use FileHandle;
use Bio::EnsEMBL::Pipeline::RunnableI;
use Bio::EnsEMBL::DnaDnaAlignFeature;
use Bio::EnsEMBL::DnaPepAlignFeature;
use Bio::EnsEMBL::PepDnaAlignFeature;
use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::SeqFeature;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::Pipeline::Runnable::FeatureFilter;
use Bio::EnsEMBL::Root;
use Bio::EnsEMBL::Pipeline::Tools::BPlite;
use Bio::EnsEMBL::Pipeline::Config::Blast;




@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI);


my %FASTA_HEADER;
my %BLAST_FLAVOUR;
my %REFILTER;

foreach my $db (@$DB_CONFIG) {
    my ($name, $header, $flavour, $refilter) = ($db->{'name'}, $db->{'header'}, $db->{'flavour'}, $db->{'refilter'});
    
    if($db && $name){
      $FASTA_HEADER{$name} = $header; 
      $BLAST_FLAVOUR{$name} = $flavour; 
      $REFILTER{$name} = $refilter;
    }else{
      my($p, $f, $l) = caller;
      warn("either db ".$db." or name ".$name." isn't defined so can't work $f:$l\n");
    }
}



=head2 new

    Title   :   new
    Usage   :   my obj =  Bio::EnsEMBL::Pipeline::Runnable::Blast->new 
    (-query    => $seq,
     -program  => 'blastp',
     -database => 'swir',
     -threshold => 1e-6,
     -filter    => $filter,
     -options   => 'V=1000000');

    Function:   Initialises Blast object
    Returns :   a Blast Object
    Args    :   A Bio::Seq object
                The blast executable (-BLAST) and database (-DB).

=cut

sub new {
    my ($class,@args) = @_;

    my $self = $class->SUPER::new(@args);    

    $self->{'_query'}     = undef;     # location of Bio::Seq object
    $self->{'_program'}   = undef;     # location of Blast
    $self->{'_database'}  = undef;     # name of database
    $self->{'_threshold'} = undef;     # Threshold for hit filterting
    $self->{'_threshold_type'} = undef;
    $self->{'_options'}   = undef;     # arguments for blast
    $self->{'_filter'}    = 1;         # Do we filter features?
    $self->{'_fplist'}    = [];        # an array of feature pairs (the output)

    $self->{'_workdir'}   = undef;     # location of temp directory
    $self->{'_filename'}  = undef;     # file to store Bio::Seq object
    $self->{'_results'}   = undef;     # file to store results of analysis
    $self->{'_prune'}     = 1;         # Don't allow hits to the same sequence
                                       # in the same region
    $self->{'_hardprune'} = 0;         # Don't force a hard limit on depth
                                       # of coverage of query
    $self->{'_coverage'}  = 10;        # Control parameter for coverage, prune
                                       # and hardprune
    $self->{'_ungapped'}  = undef;     # Do we create gapped features or not
    $self->{'_blast_re'}  = undef;

    #print STDERR "BLAST args : @args\n";
    # Now parse the input options and store them in the object
    my( $query, $program, $database, $threshold, $threshold_type, $filter,$coverage,$prune,$hardprune,$ungapped,$options) = 
            $self->_rearrange([qw(QUERY 
                                  PROGRAM 
                                  DATABASE 
                                  THRESHOLD
                                  THRESHOLD_TYPE
                                  FILTER 
                                  COVERAGE
                                  PRUNE
				  HARDPRUNE
                                  UNGAPPED
                                  OPTIONS)], 
                              @args);

    if ($query) {
      $self->query($query);
    } else {
      $self->throw("No query sequence input.");
    }

    $self->program($self->find_executable($program));
    if ($database) {
      $self->database($database);
    } else {
      $self->throw("No database input");
    }
    
    if ($options) {
      $self->options($options);
    } else {
      $self->options(' -cpus=1 ');  
    }
    
    if (defined($threshold)) {
      $self->threshold($threshold);
    }

    if (defined($threshold_type)) {
      $self->threshold_type($threshold_type);
    }

    if (defined($filter)) {
        $self->filter($filter);
    }
    
    if (defined($prune)) {
      $self->prune($prune);
    }
    if (defined($hardprune)) {
      $self->hardprune($hardprune);
    }
    if (defined($coverage)) {
      $self->coverage($coverage);
    }
    #print STDERR "setting ungapped = ".$ungapped."\n";
    if (defined($ungapped)) {
      $self->ungapped($ungapped);
    }
    return $self; # success - we hope!
}


###########
# Analysis methods
##########

=head2 run

    Title   :  run
    Usage   :   $obj->run()
    Function:   Runs blast and BPLite and creates array of feature pairs
    Returns :   none
    Args    :   none

=cut

sub run {
    my ($self, $dir) = @_;

    my $seq = $self->query || $self->throw("Query seq required for Blast\n");
   
    $self->workdir('/tmp') unless ($self->workdir($dir));
    $self->checkdir();
   
    #write sequence to file
    $self->writefile(); 
    $self->run_analysis();
    
    #parse output and create features
    $self->parse_results;
   
    $self->deletefiles();
   
    return 1;
}

sub databases {
  my ($self,@dbs) = @_;

  if (!defined($self->{_databases})) {
     $self->{_databases} = [];
  }
  if (@dbs) {
     push(@{$self->{_databases}},@dbs);
  }
  return @{$self->{_databases}};
}

=head2 run_analysis

    Title   :   run_analysis
    Usage   :   $obj->run_analysis
    Function:   Runs the blast query
    Returns :   nothing
    Args    :   none

=cut

sub run_analysis {
    my ($self) = @_;

    # This routine expands the database name into $db-1 etc for
    # split databases

    my @databases = $self->fetch_databases;
    
    $self->databases(@databases);

    foreach my $database (@databases) {
        my $db = $database;
        $db =~ s/.*\///;
        print STDERR "\n".$database."\n";
        #allow system call to adapt to using ncbi blastall. defaults to WU blast.       

        my $command = $self->program ;
        my $blastype = "";
        my $filename = $self->filename;

        if ($BLAST_FLAVOUR{$self->database} eq 'ncbi') {
            $command .= " -d $database -i $filename ";
        } else {
            $command .= " $database $filename ";
        }
        $command .= ' -gi '.$self->options. ' > '.$self->results . ".$db ";

        # Add the result file to our clean-up list.
        $self->file($self->results . ".$db");

        #$self->throw("Failed during blast run: $command". ($?/256) . " ". $!) unless (system ($command) == 0);
      }
  
}

=head2 fetch_databases

    Title   :   fetch_databases
    Usage   :   $obj->fetch_databases
    Function:   For split databases this
                method checks whether the database in $self->database
                exists.  If it doesn\'t it checks whether the
                database has been split into $database-1, $database-2 etc.
    Returns :   Array of database names
    Args    :   none

=cut


sub fetch_databases {
    my ($self) = @_;
    
    my @databases;

    my $dbname = $self->database; 
    #print "fetching databases for ".$dbname."\n";
    $dbname =~ s/\s//g;

    # prepend the environment variable $BLASTDB if
    # database name is not an absoloute path

    unless ($dbname =~ m!^/!) {
      $dbname = $ENV{BLASTDB} . "/" . $dbname;
    }

    # If the expanded database name exists put this in
    # the database array.
    #
    # If it doesn't exist then see if $database-1,$database-2 exist
    # and put them in the database array
    print STDERR "Checking if ".$dbname." exists\n";
    if (-f $dbname) {
      push(@databases,$dbname);
    } else {
      my $count = 1;
      while (-f $dbname . "-$count") 
        {
          push(@databases,$dbname . "-$count"); 	 
          $count++; 	 
        }
    }
    

    if (scalar(@databases) == 0) {
      $self->throw("No databases exist for " . $dbname);
    }

    return @databases;

}


sub get_parsers {
  my ($self)  = @_;

  my @parsers;

  foreach my $db ($self->databases) {
    $db =~ s/.*\///;

    my $fh = new FileHandle;
    $fh->open("<" . $self->results . ".$db");
    
    my $parser = new Bio::EnsEMBL::Pipeline::Tools::BPlite ('-fh' => $fh);
    
    push(@parsers,$parser);
  } 

  return @parsers;
}


=head2 parse_results

    Title   :   parse_results
    Usage   :   $obj->parse_results
    Function:   Parses the blast results and stores the output in
                an array of feature pairs.  The output will be filtered
                by some threshold.
                If we input a filehandle to this method the results 
                will be read from this file rather than the filename 
                stored in the object.  This is handy for processing
                blast runs that have been run previously.
    Returns :   @Bio::EnsEMBL::FeaturePair
    Args    :   Optional filehandle. 
                 

=cut

sub parse_results {
  my ($self,$fh) = @_;

  my %ids;
  my $count = 0;
  my @parsers;
 
  if (defined($fh)) {
    my $parser = new Bio::EnsEMBL::Pipeline::Tools::BPlite(-fh => $fh);
    push(@parsers,$parser);
  } else {
    @parsers = $self->get_parsers;
  }

  if ($self->filter) {
    my @hsps;
    foreach my $parser (@parsers) {
      while (my $sbjct = $parser->nextSbjct) {
        while (my $hsp = $sbjct->nextHSP) {
          push(@hsps,$hsp);
        }
      }
    }
    %ids = $self->filter_hits(@hsps);
    
  }
 
  #print STDERR "Ids " . keys(%ids) . "\n";

  @parsers = ();
  seek $fh, 0, 0 if ref $fh;

  if (defined($fh)) {
    my $parser = new Bio::EnsEMBL::Pipeline::Tools::BPlite(-fh => $fh);
    push(@parsers,$parser);
  } else {
    @parsers = $self->get_parsers;
  }
  my $db = $self->database;
  
  my $re = $self->get_regex($db);
  if(!$re){
    $self->throw("no regex defined for ".$db);
  }
  #print STDERR "have ".$re." regular expression\n";
  foreach my $parser (@parsers) {
    # print STDERR "New parser\n";
  NAME: while  ( my $sbjct =$parser->nextSbjct) {
      
    my $fasta_header = $sbjct->name ;     
    my ($name) = $fasta_header =~ /$re/;
    unless ($name) {
        $self->throw("Error getting a valid accession from \"" .
        $fasta_header .
        "\"; check your blast config and / or blast headers");
    }

    # print STDERR "Name " . $fasta_header . "\n";
     if (($self->filter == 1) && !defined($ids{$fasta_header})) {
      next NAME;
    }

    #print "Parsing name $name\n";
    
  HSP: while (my $hsp = $sbjct->nextHSP) {
      
      if ($self->threshold_type eq "PID") {
	next HSP
	  if defined $self->threshold and ($hsp->percent < $self->threshold);
      } elsif ($self->threshold_type eq "SCORE") {
	next HSP
	  if defined($self->threshold) and ($hsp->score < $self->threshold);
      } elsif ($self->threshold_type eq "PVALUE") {
	next HSP
	  if defined($self->threshold) and ($hsp->P > $self->threshold);
      }
      # Each HSP is a gapped alignment.
      # This method splits the gapped alignment into
      # ungapped pieces
      #print STDERR "HSP " . $hsp->P . "\n";
      $count++;
      $self->split_HSP($hsp,$name);
        
    }
  }
 }


# Alternate feature filter. If option not present in config, should default to FeatureFilter -prune

  if ($REFILTER{$self->database}){
    # re-filter, with pruning - rewrotee to use a local select_feature function instead of FeatureFilter 
    my @selected_features = $self->select_features($self->output);
    $self->output(@selected_features);
  } else {
    # re-filter, with pruning
    my @allfeatures = $self->output;
    if ($self->threshold_type eq "PID") {
      @allfeatures = sort {$b->percent_id <=> $a->percent_id} @allfeatures;
    } elsif ($self->threshold_type eq "SCORE") {
      @allfeatures = sort {$b->score <=> $a->score} @allfeatures;
    } else {
      @allfeatures = sort {$a->p_value <=> $b->p_value} @allfeatures;
    }
    if ($self->filter) {
      my $search = new Bio::EnsEMBL::Pipeline::Runnable::FeatureFilter(
                                         -prune     => $self->prune,
                                         -hardprune => $self->hardprune,
                                         -coverage  => $self->coverage);
      my @pruned = $search->run(@allfeatures);

      #print STDERR "dbg ", scalar(@allfeatures), " ", scalar(@pruned), "\n";
      $self->output(@pruned);
    }
  }
  #print "have parsed resultz\n";
  return $self->output;

}

sub prune {
  my ($self,$arg) = @_;

  if (defined($arg)) {
    $self->{_prune} = $arg;
  }
  return $self->{_prune};
}

sub hardprune {
  my ($self,$arg) = @_;

  if (defined($arg)) {
    $self->{_hardprune} = $arg;
  }
  return $self->{_hardprune};
}

sub coverage {
  my($self,$arg) = @_;

  if (defined($arg)) {
    $self->{_coverage} = $arg;
  }
  return $self->{_coverage};
}

sub filter_hits {
  my ($self,@hsps) = @_;

  my %ids;

  my @features;

  
    
 HSP: foreach my $hsp (@hsps) {
      
      my $name = $hsp->subject->seqname ;

      if ($self->threshold_type eq "PID") {
        next HSP
	  if defined $self->threshold and ($hsp->percent < $self->threshold);
      } elsif ($self->threshold_type eq "SCORE") {
        next HSP
	  if defined $self->threshold and ($hsp->score < $self->threshold);
      } elsif ($self->threshold_type eq "PVALUE") {
        next HSP
	  if defined $self->threshold and ($hsp->P > $self->threshold);
      } 
      
      my $qstart = $hsp->query->start();
      my $hstart = $hsp->subject->start();
      
      my $qend   = $hsp->query->end();
      my $hend   = $hsp->subject->end();      
      
      my $qstrand = $hsp->query->strand(),
      my $hstrand = $hsp->subject->strand();

      my $score  = $hsp->score;

      my $feature1 = new Bio::EnsEMBL::SeqFeature();
      $feature1->start($qstart);
      $feature1->end  ($qend);
      $feature1->strand($qstrand);
      $feature1->score($score);

      

      my $feature2 = new Bio::EnsEMBL::SeqFeature();
      $feature2->start  ($hstart);
      $feature2->end    ($hend);
      $feature2->strand ($hstrand);
      $feature2->score  ($score);
      $feature2->seqname($name);

        my $fp = new Bio::EnsEMBL::FeaturePair(-feature1 => $feature1,
                                               -feature2 => $feature2);
      $fp->p_value($hsp->P);
      $fp->percent_id($hsp->percent);

      push(@features,$fp);
    }

  if ($self->threshold_type eq "PID") {
    @features = sort {$b->percent_id <=> $a->percent_id} @features;
  } elsif ($self->threshold_type eq "SCORE") {
    @features = sort {$b->score <=> $a->score} @features;
  } elsif ($self->threshold_type eq "PVALUE") {
    @features = sort {$a->p_value <=> $b->p_value} @features;
  } 
  
  my $search = new Bio::EnsEMBL::Pipeline::Runnable::FeatureFilter(-coverage => $self->coverage);
  
  my @newfeatures = $search->run(@features);
  
  foreach my $f (@newfeatures) {
    my $id = $f->hseqname;
    
    $ids{$id} = 1;
  }
  
  return %ids;
}
    
=head2 split_HSP

    Title   :   split_HSP
    Usage   :   $obj->split_HSP
    Function:   Takes a gapped blast HSP 
                and turns it into an array of ungapped feature pairs.
    Returns :   Nothing
    Args    :   BPlite::HSP,name string
                 

=cut

sub split_HSP {
    my ($self,$hsp,$name) = @_;
    # First of all some jiggery pokery to find out what sort of alignment
    # we have - (dna-dna, dna-pep, pep-dna etc).
    # We also work out which strand each sequence is on.
    # 
    # Both these variables are needed to work out the increment
    # values for the query and the hit sequences.  As we split
    # up the gapped alignment we need to increment the coordinates.
    #
    # We have to take care with dna-pep alignments to increment the dna
    # sequences by 3 as there are 3 bases in a codon.
    #
    # For dna-dna alignments this will be 1 or -1 depending on
    # the strand.  
    #
    # For pep-dna alignments and vice-versa the increments will be +-3 
    # for the dna sequence and +- 1 for the peptide sequence. 

    # Find out which Blast program generated the results
    my $source = $self->program;             
    $source =~ s/\/.*\/(.*)/$1/;
   
    # Use the Blast analysis source to determine what the sequences are
    my ($qtype,  $htype)   = $self->_findTypes($source);
    
    my $qstrand = $hsp->query->strand(),
    my $hstrand = $hsp->subject->strand();
   
    my ($qinc,   $hinc)    = $self->_findIncrements($hsp,$qstrand,$hstrand,$qtype,$htype);
   
    #print STDERR "Alignment q : " . $hsp->query->start . "\t" . $hsp->query->end . "\t" . $hsp->querySeq . "\n";
    #print STDERR "Alignment s : " . $hsp->subject->start . "\t" . $hsp->subject->end . "\t" . $hsp->sbjctSeq . "\n";
#    print STDERR "types (increments) $qtype ($qinc) : $htype ($hinc) Strands : $qstrand $hstrand $name\n";

    #if ($qtype eq "dna" && $htype eq "dna") {exit()};
    
    # We split the alignment strings into arrays of one char each.  
    # We then loop over this array and when we come to a gap
    # in either the query sequence or the hit sequence we make 
    # a new feature pair. 
    #

    # Before making a new feature pair we check whether we have
    # something to make a feature pair out of - i.e. we aren't just on
    # the 2nd or greater base in a gapped region.  
    #
    # As we loop over the array we need to keep track of both the
    # query coordinates and the hit coordinates.  We track the last
    # start of a feature pair and also the current end of a feature
    # pair.  The feature pair start is reset every time we hit a
    # gapped position and the feature pair end is reset every time we
    # hit an aligned region.

    my @gap;


    my @qchars = split(//,$hsp->querySeq);  # split alignment into array of char
    my @hchars = split(//,$hsp->sbjctSeq);  # ditto for hit sequence
    
    my $qstart = $hsp->query->start();                # Start off the feature pair start
    my $hstart = $hsp->subject->start();              # ditto

    my $qend   = $hsp->query->start();                  # Set the feature pair end also
    my $hend   = $hsp->subject->start();                # ditto

   if ($qstrand == -1) {
      $qstart = $hsp->query->end;
      $qend   = $hsp->query->end;
    }
    if ($hstrand == -1) {
      $hstart = $hsp->subject->end;
      $hend   = $hsp->subject->end;
    }

    #print STDERR "hstart ". $hstart . "\t" . " hend " . $hend . " qstart " . $qstart . " qend " . $qend . "\n";
    
    my $count = 0;                                # counter for the bases in the alignment
    my $found = 0;                                # flag saying whether we have a feature pair
    

    my @tmpf;

    my $analysis = new Bio::EnsEMBL::Analysis(-db              => $self->database,
                                              -db_version      => 1,
                                              -program         => $source,
                                              -program_version => 1,
                                              -gff_source      => $source,
                                              -gff_feature     => 'similarity',
                                              -logic_name      => 'blast');
    
    # Here goes...
  
    while ($count <= $#qchars) {
        # We have hit an ungapped region.  Increase the query and hit counters.
        # and flag that we have a feature pair.
     
        if ($qchars[$count] ne '-' &&
            $hchars[$count] ne '-') {

            $qend += $qinc;
            $hend += $hinc;
            
            $found = 1;
        } else {

            # We have hit a gapped region.  If the feature pair flag is set ($found)
            # then make a feature pair, store it and reset the start and end variables.

            if ($found == 1) {

	      #print STDERR "hstart ". $hstart . "\t" . " hend " . $hend . " qstart " . $qstart . " qend " . $qend . "\n";

	      my $fp = $self->_convert2FeaturePair($qstart,$qend,$qstrand,$hstart,$hend,$hstrand,$qinc,$hinc,$hsp,$name, $analysis);
	      #print "Found " . $fp->gffstring . "\n";                
	      push(@tmpf,$fp);
	      #$self->growfplist($fp);                                         
            }
        
            # We're in a gapped region.  We need to increment the sequence that
            # doesn't have the gap in it to keep the coordinates correct.
            # We also need to reset the current end coordinates.

            if ($qchars[$count] ne '-') {
                $qstart = $qend   + $qinc;
            } else {
                $qstart = $qend;
            }
            if ($hchars[$count] ne '-') {
                $hstart = $hend   + $hinc;
            } else {
                $hstart = $hend;
            }
            
            $qend = $qstart;
            $hend = $hstart;

            $found = 0;
        }
        $count++;
    }                        
    
    # Remember the last feature
    if ($found == 1) {

        my $fp = $self->_convert2FeaturePair($qstart,$qend,$qstrand,$hstart,$hend,$hstrand,$qinc,$hinc,$hsp,$name, $analysis);
        #print "Found " . $fp->gffstring . "\n";
        push(@tmpf,$fp);
#       $self->growfplist($fp);                                         
    }
    #print STDERR "ungapped = ".$self->ungapped."\n";
    if ($self->ungapped) {
      foreach my $f (@tmpf) {
        $self->warn("can't store feature pairs this will fail\n");
        $self->growfplist($f);                                         
      } 
    } else {
      # Which type of feature do we want?
      my $fp;
      
     
      $qinc = abs( $qinc );
      $hinc = abs( $hinc );

      if( $qinc == 3 && $hinc == 1 ) {
        $fp = Bio::EnsEMBL::DnaPepAlignFeature->new(-features => \@tmpf);
      } elsif( $qinc == 1 && $hinc == 3 ) {
        $fp = Bio::EnsEMBL::PepDnaAlignFeature->new(-features => \@tmpf);
      } elsif( $qinc == 1 && $hinc == 1 ) {
        $fp = Bio::EnsEMBL::DnaDnaAlignFeature->new(-features => \@tmpf);
      } else {
        $self->throw( "Hardcoded values wrong?? " );
      }

      # helps debugging subsequent steps
      $fp->{'qseq'} = $hsp->querySeq();
      $fp->{'sseq'} = $hsp->sbjctSeq();

      $self->growfplist($fp);
    }
}

=head2 _convert2FeaturePair

    Title   :   convert2FeaturePair
    Usage   :   $obj->convert2FeaturePair
    Function:   Internal function taking a set of coords and
                converts them into a feature pair
    Returns :   Bio::EnsEMBL::FeaturePair
    Args    :   int,int,int,int,int,int,int,int,BPlite::HSP
                 

=cut

sub _convert2FeaturePair {
    my ($self,$qstart,$qend,$qstrand,$hstart,$hend,$hstrand,$qinc,$hinc,$hsp,$name,$analysis) = @_;

    # The actual end of the alignment is the previous character.
    
    my $tmpqend = $qend; $tmpqend -= $qinc;
    my $tmphend = $hend; $tmphend -= $hinc;
    
    my $tmpqstart = $qstart;
    my $tmphstart = $hstart;
    
    # This is for dna-pep alignments.  The actual end base
    # will be +- 2 bases further on.
    if (abs($qinc) > 1) {
        $tmpqend += $qstrand * 2;
    }
    if (abs($hinc) > 1) {
        $tmphend += $hstrand * 2;
    }
    
    # Make sure start is always < end
    if ($tmpqstart > $tmpqend) {
        my $tmp    = $tmpqstart;
        $tmpqstart = $tmpqend;
        $tmpqend   = $tmp;
    }
    if ($tmphstart > $tmphend) {
        my $tmp    = $tmphstart;
        $tmphstart = $tmphend;
        $tmphend   = $tmp;
    }
    
    # print STDERR "Creating feature pair " . $tmpqstart . "\t" . $tmpqend . "\t" . $qstrand . "\t" . $tmphstart . "\t" . $tmphend . "\t" . $hstrand . "\t" . $name . "\n";

    my $fp = $self->_makeFeaturePair($tmpqstart,$tmpqend,$qstrand,$tmphstart,$tmphend,$hstrand,$hsp->score,
                                     $hsp->percent,$hsp->P,$name,$analysis);

    return $fp;
}

=head2 _makeFeaturePair

    Title   :   _makeFeaturePair
    Usage   :   $obj->_makeFeaturePair
    Function:   Internal function that makes feature pairs
    Returns :   Bio::EnsEMBL::FeaturePair
    Args    :   int,int,int,int,int,int,int,int,BPlite::HSP
                 

=cut

sub _makeFeaturePair {
    my ($self,$qstart,$qend,$qstrand,$hstart,$hend,$hstrand,$score,$pid,$evalue,$name,$analysis)  = @_;

    my $source = $self->program;             
    $source =~ s/\/.*\/(.*)/$1/;

    my $feature1 = new Bio::EnsEMBL::SeqFeature(-seqname     => $self->query->id,
                                                -start       => $qstart,
                                                -end         => $qend,
                                                -strand      => $qstrand,
                                                -analysis    => $analysis,
                                                -score       => $score);
        
    $feature1->percent_id($pid);
    $feature1->p_value($evalue);


    my $feature2 = new Bio::EnsEMBL::SeqFeature(-seqname => $name,
                                                -start   => $hstart,
                                                -end     => $hend,
                                                -strand  => $hstrand,
                                                -analysis => $analysis,
                                                -score    => $score);
    
    my $fp = new Bio::EnsEMBL::FeaturePair(-feature1 => $feature1,
                                           -feature2 => $feature2);
   
    $feature2->percent_id($pid);
    $feature2->p_value($evalue);
 
    return $fp;
}

sub _findIncrements {
    my ($self,$hsp,$qstrand,$hstrand,$qtype,$htype) = @_;

    my $qinc   = 1 * $qstrand;
    my $hinc   = 1 * $hstrand;
    
    if ($qtype eq 'dna' && $htype eq 'pep') {
        $qinc = 3 * $qinc;
    } 
    if ($qtype eq 'pep' && $htype eq 'dna') {
        $hinc = 3 * $hinc;
    }

    return ($qinc,$hinc);
}


sub _findTypes {
    my ($self,$blast_type) = @_;

    my $query;
    my $sbjct;

    # Determine the type of sequences used in the Blast analysis
    # based on the following table:
    #
    # Analysis  Query  Subject
    #  blastp    pep     pep
    #  blastn    dna     dna
    #  blastx    dna     pep
    # tblastn    pep     dna
    # tblastx    dna     dna

    if ($blast_type =~ /tblastn/) {
      $query = 'pep';
      $sbjct = 'dna';
    }
    elsif ( $blast_type =~ /tblastx/ || $blast_type =~ /blastn/) {
      $query = 'dna';
      $sbjct = 'dna';
    }
    elsif ( $blast_type =~ /blastp/ ) {
      $query = 'pep';
      $sbjct = 'pep';
    }
    elsif ( $blast_type =~ /blastx/ ) {
      $query = 'dna';
      $sbjct = 'pep';
    }
    else {
      $self->throw("Unknown Blast analysis:" . $blast_type . "\n");
    }
    
    return ($query,$sbjct);
}

sub select_features {
 
        my ($self,@features) = @_;
 
        @features= sort {
                $a->strand<=> $b->strand
                       ||
                $a->start<=> $b->start
        } @features;
 
        my @selected_features;
 
        my $best_hit = $features[0];
 
        foreach my $feat (@features){
                if ($feat->overlaps($best_hit,'strong')) {
                        if ($feat->score > $best_hit->score) {
                                $best_hit = $feat;
                        }
                        }else {
                                push (@selected_features,$best_hit);
                                $best_hit = $feat;
                        }
        }
         
        my @newfeatures;        
        FEAT: foreach my $feat (@features){
                foreach my $sf (@selected_features){
                        if (($sf->hseqname eq $feat->hseqname)&&($sf->score == $feat->score)){
                                push (@newfeatures, $feat);
                                next FEAT;        
                        }        
                }
        }
        

        return @newfeatures;
}

##############
# input/output methods
#############

=head2 output

    Title   :   output
    Usage   :   obj->output()
    Function:   Returns an array of feature pairs
    Returns :   Returns an array of feature pairs
    Args    :   none

=cut

sub output {
    my ($self, @arg) = @_;

    if (@arg) {
      @{$self->{'_fplist'}} = @arg;
    }

    if (!defined($self->{'_fplist'})) {
       $self->{'_fplist'} = [];
    }
    return @{$self->{'_fplist'}};
}


#################
# get/set methods 
#################

sub query {
    my ($self, $seq) = @_;
    if ($seq) {
      unless ($seq->isa("Bio::PrimarySeqI") || $seq->isa("Bio::Seq")) {
        $self->throw("Input isn't a Bio::Seq or Bio::PrimarySeq");
      }

      $self->{'_query'} = $seq ;

      $self->filename($seq->id.".$$.seq");
      $self->results($self->filename.".blast.out");

      # Add file to list for later cleanup.
    }
    return $self->{'_query'};
}

=head2 program

    Title   :   program
    Usage   :   $obj->program('/usr/local/pubseq/bin/blastn');
    Function:   Get/set method for the location of the blast flavour
    Returns :   string
    Args    :   string

=cut

sub program {
  my ($self, $location) = @_;
  
  if ($location) {
    $self->throw("executable not found at $location: $!\n")     unless (-e $location && -x $location);
    $self->{'_program'} = $location ;
  }
  return $self->{'_program'};
}

=head2 database

    Title   :   database
    Usage   :   $obj->database('dbEST');
    Function:   Get/set method for the location of database
    Args    :   none

=cut

sub database {
  my ($self, $db) = @_;
  my $re_string;
    
  if (defined($db)) {
    $self->{'_database'} = $db;
    $db =~ s!.*/!!;       
    #print STDERR "have ".keys(%$FASTA_HEADER)." blast regexs defined\n";
  }
  
  return $self->{'_database'};
}

=head2 options

    Title   :   options
    Usage   :   $obj->options(' -I ');
    Function:   Get/set method for blast arguments
    Args    :   File path (optional)

=cut

sub options {
  my ($self, $args) = @_;
  
  if (defined($args)) {
    $self->{'_options'} = $args ;
  }
  return $self->{'_options'};
}

sub ungapped {
  my ($self,$arg) = @_;

  if (defined($arg)) {
    $self->{_ungapped} = $arg;
  }
  return $self->{_ungapped};
}



sub filter {
    my ($self,$args) = @_;

    if (defined($args)) {
        if ($args != 0 && $args != 1) {
            $self->throw("Filter option must be 0 or 1");
        }
        $self->{'_filter'} = $args;
    }
    return $self->{'_filter'};
}

sub get_threshold_types {
  my ($self) = @_;

  return ("PID","PVALUE","SCORE");
}

sub threshold_type {
  my ($self,$type) = @_;

  my %allowed = map { $_, 1 } $self->get_threshold_types;
  
  if ($type) {
    unless (defined $allowed{$type}) {
      $self->throw("Unallowed threshold type $type");
    }
    $self->{'_threshold_type'} = $type;
  }

  return $self->{'_threshold_type'};
}


sub blast_re {
  my ($self, $re_string) = @_;
  my $re;

  if (defined $re_string) {
    eval {
      $re = $re_string;
    };
    if ($@) {
      $self->throw("Illegal RE string $re_string");
    }
    $self->{'_blast_re'} = $re;
  }

  return $self->{'_blast_re'};
}


sub add_regex{
  my ($self, $name, $re_string) = @_;
  #print STDERR "adding regex ".$re_string." for ".$name."\n";
  $FASTA_HEADER{$name} = $re_string;
  
}

sub get_regex{
  my ($self, $name) = @_;
  #print STDERR "getting the regex for ".$name."\n";
  if($name =~/\/tmp\//){
    $name =~ s/\/tmp\///g;
  }
  return $FASTA_HEADER{$name};
}
1;

