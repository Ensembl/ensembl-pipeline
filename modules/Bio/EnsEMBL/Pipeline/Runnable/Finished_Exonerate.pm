### Bio::EnsEMBL::Pipeline::Runnable::Finished_Exonerate

package Bio::EnsEMBL::Pipeline::Runnable::Finished_Exonerate;

use strict;
use Bio::EnsEMBL::Pipeline::RunnableI;
use Tie::File;
use Fcntl;
use Bio::EnsEMBL::Pipeline::SeqFetcher::FileIndex;

use vars qw(@ISA $verbose);
@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI);

$verbose = 1;

$ENV{BLASTDB} = '/data/blastdb/Finished';

sub new {
	my ( $class, @args ) = @_;
	my $self = bless {}, $class;

	my ( $query_type, $query_db, $query_seqs, $query, $target, $analysis, $exo_options ) = $self->_rearrange(
		[
			qw(
			  QUERY_TYPE
			  QUERY_DB
			  QUERY_SEQS
			  QUERY
			  TARGET
			  ANALYSIS
			  EXO_OPTIONS
			  )
		],
		@args
	);
	$self->analysis($analysis);
	if ( defined($target) ) {
		$self->target($target);
	}
	if ( defined($query_db) ) {
		$self->query_db($query_db);
	}
	if ( defined($query) ) {
		$self->query($query);
	}
	if (defined($query_seqs)) {
	    if(ref($query_seqs) ne "ARRAY"){
	    	throw("You must supply an array reference with -query_seqs");
	    }
    	$self->query_seqs($query_seqs);
  	}
	if ( defined($exo_options) ) {
		$self->exo_options($exo_options);
	}
	if ( defined($query_type) ) {
    	$self->query_type($query_type);
  	} else{
    	$self->query_type('dna');
  	}
  	# These are what drives how we gather up the output
  	my $basic_options = "--showsugar false --showvulgar false --showalignment false --ryo \"RESULT: %S %pi %ql %tl %g %V\\n\" ";
	if ($self->options){
    	$basic_options .= $self->options;
  	}
  	$self->options($basic_options);
  	my $seqfetcher = $self->_make_seqfetcher($self->analysis->db_file);
    $self->seqfetcher($seqfetcher);


	return $self;
}

sub run {
	my ($self) = @_;
	my @query;

	if ( $self->query || $self->query_seqs ) {

		my $seq =  $self->query ? $self->query : $self->query_seqs;

		# Write query sequence(s) to file if necessary
		my $query_file = $self->write_seq_file( $seq );

		# register the file for deletion
		$self->files_to_delete($query_file);
		push @query, $query_file;

	} elsif ( $self->query_db ) {

		push @query, @{ $self->fetch_databases };
	}
	else {
		$self->throw('No query provided');
	}

	# Write target sequence to file
	my $target_file = $self->create_filename( 'target_seq', 'fa' );
	$self->write_seq_file( $self->target, $target_file );

	# register the file for deletion
	$self->files_to_delete($target_file);

	my $start = time;

	foreach my $q (@query) {
		# Build exonerate command
		my $command = join(
			'  ',
			(
				$self->analysis->program,
				$self->options, $self->exo_options,
				' --querytype ',  $self->query_type,
				' --targettype ', $self->target_type,
				' --query ',      $q,
				' --target ',     $target_file
			)
		);

		# Execute command and parse results
		print STDOUT "Running exonerate : $command\n";

		my $exo_fh;
		open( $exo_fh, "$command |" )
		  or $self->throw("Error opening exonerate command <$command>\n $? $!");

		$self->add_output($self->parse_results( $exo_fh ));

		close($exo_fh) or $self->throw("Error closing exonerate command: $? $!");

	}

	print STDOUT "=======================> Ran exonerate on @query  ".(time-$start)." secs (".time.")\n";

	$self->delete_files;

	return 1;
}

sub files_to_delete {
  my ($self, $file) = @_;
  if(!$self->{'del_list'}){
    $self->{'del_list'} = {};
  }
  if($file){
    $self->{'del_list'}->{$file} = 1;
  }
  return $self->{'del_list'};
}

sub files_to_protect {
  my ($self, $file) = @_;
  if(!$self->{'protect_list'}){
    $self->{'protect_list'} = {};
  }
  if($file){
    $self->{'protect_list'}->{$file} = 1;
  }
  return $self->{'protect_list'};
}

sub delete_files {
  my ($self, $filehash, $protected_hash) = @_;
  if(!$filehash){
    $filehash = $self->files_to_delete;
  }
  if(!$protected_hash){
    $protected_hash = $self->files_to_protect;
  }
  foreach my $name (keys(%$filehash)){
    if(!$protected_hash->{$name}){
      unlink $name;
    }
  }
  my @protected = keys(%$protected_hash);
  return \@protected;
}



sub fetch_databases {

	my ($self) = @_;
	my @databases;
	my $db_names = $self->query_db;
	$db_names =~ s/\s//g;

	if ( not exists $self->{databases} ) {
		$self->{databases} = [];
	}

	foreach my $dbname ( split( ",", $db_names ) )
	{    # allows the use of a comma separated list in $self->database
		    # prepend the environment variable $BLASTDB if
		    # database name is not an absoloute path
		unless ( $dbname =~ m!^/! ) {
			$dbname = $ENV{BLASTDB} . "/" . $dbname;
		}

		# If the expanded database name exists put this in
		# the database array.
		#
		# If it doesn't exist then see if $database-1,$database-2 exist
		# and put them in the database array
		if ( -f $dbname ) {
			push( @databases, $dbname );
		}
		else {
			my $count = 1;
			my $db_filename;
			while ( -f ( $db_filename = "${dbname}-${count}" ) ) {
				push( @databases, $db_filename );
				$count++;
			}
			$! = undef
			  ; # to stop pollution as it will be "No such file or directory" after while loop above.
		}
	}
	if ( scalar(@databases) == 0 ) {
		$self->throw( "No databases exist for " . $db_names );
	} else {
		foreach my $db_name (@databases){
			$self->get_db_version($db_name) if $db_name =~ /emnew_/;
		}
		$self->get_db_version($databases[0]);
		push @{$self->{databases}}, @databases;
	}

	return $self->{databases};

}

sub get_db_version {
	my ( $self, $db ) = @_;
	my $debug_this = 0;    # this just shows debug info.
	my $force_dbi  = 0;    # this will force a dbi call SLOW!!!!!!
	unless ( $self->{'_db_version_searched'} ) {
		if ($db) {
			$BlastableVersion::debug = $debug_this;
			$self->warn("BlastableVersion is cvs revision $BlastableVersion::revision \n") if $debug_this;

			my $ver = eval {
				my $blast_ver = BlastableVersion->new();
				$blast_ver->force_dbi($force_dbi);    # if set will be SLOW.
				$blast_ver->get_version($db);
				$blast_ver;
			};
			$self->throw("I failed to get a BlastableVersion for $db") if $@;

			my $dbv  = $ver->version();
			my $sgv  = $ver->sanger_version();
			my $name = $ver->name();
			my $date = $ver->date();
			unless ($dbv) {
				$self->throw(  "I know nothing about $db I tried to find out:\n"
					  . " - name <"
					  . $name . ">\n"
					  . " - date <"
					  . $date . ">\n"
					  . " - version <"
					  . $dbv . ">\n"
					  . " - sanger_version <"
					  . $sgv
					  . ">\n" );
			}
			$self->{'_db_version_searched'} = $dbv;
		}
		else {
			$self->throw(  "You've asked about what I searched, but I don't know."
				  . " It's not set. I need to be called with a database filename first"
			);

			# The code probably got here because of a problem with the MLDBM
			# cache file on the machine this was running on.
			# the cache file is stored @ /var/tmp/blast_versions
			# try <rm -f /var/tmp/blast_versions>
		}
	}
	return $self->{'_db_version_searched'};
}

sub target {
	my ( $self, $seq ) = @_;
	if ($seq) {
		unless ( $seq->isa("Bio::PrimarySeqI") || $seq->isa("Bio::SeqI") ) {
			$self->throw("target seq must be a Bio::SeqI or Bio::PrimarySeqI");
		}
		$self->{_target} = $seq;
	}
	return $self->{_target};
}

sub query {
  my $self = shift;
  my $slice = shift;
  if($slice){
    throw("Must pass Runnable::query a Bio::PrimarySeqI not a ".
          $slice) unless($slice->isa('Bio::PrimarySeqI'));
    $self->{'query'} = $slice;
  }
  return $self->{'query'};
}

sub query_seqs {
  my ($self, $seqs) = @_;
  if ($seqs){
    unless ($seqs->[0]->isa("Bio::PrimarySeqI") || $seqs->[0]->isa("Bio::SeqI")){
      throw("query seq must be a Bio::SeqI or Bio::PrimarySeqI");
    }
    $self->{_query_seqs} = $seqs;
  }
  return $self->{_query_seqs};
}

sub query_db {
	my ( $self, $db ) = @_;
	if ($db) {
		$self->{_query_db} = $db;
	}
	return $self->{_query_db};
}

sub options {
  my $self = shift;
  $self->{'options'} = shift if(@_);
  return $self->{'options'} || '';
}

sub exo_options {
  my $self = shift;
  $self->{'exo_options'} = shift if(@_);
  return $self->{'exo_options'} || '';
}

sub query_type {
  my ($self, $mytype) = @_;
  if (defined($mytype) ){
    my $type = lc($mytype);
    unless( $type eq 'dna' || $type eq 'protein' ){
      $self->throw("not the right query type: $type");
    }
    $self->{_query_type} = $type;
  }
  return $self->{_query_type};
}

sub target_type {
  my ($self) = @_;

  # the target type has to be DNA, because we are making transcripts

  return 'dna';
}

sub analysis{
  my $self = shift;
  my $analysis = shift;
  if($analysis){
    $self->throw("Must pass RunnableDB:analysis a Bio::EnsEMBL::Analysis".
          "not a ".$analysis) unless($analysis->isa
                                     ('Bio::EnsEMBL::Analysis'));
    $self->{'analysis'} = $analysis;
  }
  return $self->{'analysis'};
}

sub seqfetcher {
    my ( $self, $seqfetcher ) = @_;

    if ($seqfetcher) {
        $self->{'_seqfetcher'} = $seqfetcher;
    }
    return $self->{'_seqfetcher'};
}

sub _make_seqfetcher {
    my ( $self, $fasta_file ) = @_;

    my $seqfetcher;
    if ( $fasta_file){
        $seqfetcher = Bio::EnsEMBL::Pipeline::SeqFetcher::FileIndex->new( -seqfile => $fasta_file);
    }

    return $seqfetcher;
}

sub add_output {
    my ( $self, $feat ) = @_;

    my $ana = $self->analysis;
    foreach my $f (@$feat) {
        $f->analysis($ana);
    }

    $self->{'_output'} ||= [];

    push ( @{ $self->{'_output'} }, @$feat );

}


sub create_filename{
  my ($self, $stem, $ext, $dir) = @_;
  if(!$dir){
    $dir = '/tmp';
  }
  $stem = '' if(!$stem);
  $ext = '' if(!$ext);
  $self->throw($dir." doesn't exist Runnable:create_filename") unless(-d $dir);
  my $num = int(rand(100000));
  my $file = $dir."/".$stem.".".$$.".".$num.".".$ext;
  while(-e $file){
    $num = int(rand(100000));
    $file = $dir."/".$stem.".".$$.".".$num.".".$ext;
  }
  return $file;
}

sub write_seq_file{
  my ($self, $seq, $filename, $format) = @_;
  $format = 'fasta' if(!$format);
  my @seqs;
  if(ref($seq) eq "ARRAY"){
    @seqs = @$seq;
    $self->throw("Seqs need to be Bio::PrimarySeqI object not a ".$seqs[0])
      unless($seqs[0]->isa('Bio::PrimarySeqI'));
  }else{
    $self->throw("Need a Bio::PrimarySeqI object not a ".$seq)
      if(!$seq || !$seq->isa('Bio::PrimarySeqI'));
    @seqs = ($seq);
  }
  $filename = $self->create_filename('seq', 'fa', '/tmp')
    if(!$filename);
  my $seqout = Bio::SeqIO->new(
                               -file => ">".$filename,
                               -format => $format,
                              );
  foreach my $seq(@seqs){
    eval{
      $seqout->write_seq($seq);
    };
    if($@){
      $self->throw("FAILED to write $seq to $filename SequenceUtils:write_seq_file $@");
    }
  }
  return $filename;
}

sub parse_results {
	my ($self, $fh) = @_;
	my %strand_lookup = ( '+' => 1, '-' => -1, '.' => 1 );

	my @exon_sup_features;

	TRANSCRIPT:
	while (<$fh>){
		print STDERR $_ if $self->verbose;
		next unless /^RESULT:/;
		chomp;

		my (
			$tag,      $q_id,             $q_start,
			$q_end,    $q_strand,         $t_id,
			$t_start,  $t_end,            $t_strand,
			$score,    $perc_id,          $q_length,
			$t_length, $gene_orientation, @align_components
		  )
		  = split ;

		$t_strand         = $strand_lookup{$t_strand};
		$q_strand         = $strand_lookup{$q_strand};
		$gene_orientation = $strand_lookup{$gene_orientation};

		# Read vulgar information and extract exon regions.
		my $exons = $self->_parse_vulgar_block(
			$t_start,  $t_end,    $t_strand,
			$t_length, $q_start,  $q_end,
			$q_strand, $q_length, \@align_components
		);

		# now we have extracted the exons we can use
		# the gene_orienation to flip the strands if necessary
		if ( $gene_orientation == -1 ) {
			$t_strand *= -1;
			$q_strand *= -1;
		}

		my $covered_count = 0;
		my $coverage      = 1;
		if ($coverage) {
			foreach my $exon (@$exons) {
				foreach my $sf ( @{ $exon->{sf} } ) {
					$covered_count += $sf->{query_end} - $sf->{query_start} + 1;
				}
			}
		}
		else {
			$covered_count = abs( $q_end - $q_start );
		}

		$coverage = sprintf( "%.2f", 100 * $covered_count / $q_length );

		# Build FeaturePairs for each region of query aligned to a single
		# Exon.  Create a DnaDnaAlignFeature from these FeaturePairs and then
		# attach this to our Exon.
		my $transcript = Bio::EnsEMBL::Transcript->new();

		my (@tran_feature_pairs);

		foreach my $proto_exon (@$exons) {

			# Build our exon and set its key values.
			my $exon = Bio::EnsEMBL::Exon->new();

			$exon->seqname($t_id);
			$exon->start( $proto_exon->{exon_start} );
			$exon->end( $proto_exon->{exon_end} );
			$exon->phase( $proto_exon->{phase} );
			$exon->end_phase( $proto_exon->{end_phase} ) if $proto_exon->{end_phase};
			$exon->strand($t_strand);
			my @exon_feature_pairs;
			foreach my $sf ( @{ $proto_exon->{sf} } ) {
				my $ft1 = Bio::EnsEMBL::SeqFeature->new(
					-seqname    => $t_id,
					-start      => $sf->{target_start},
					-end        => $sf->{target_end},
					-strand     => $t_strand,
					-score      => $score,
					-percent_id => $perc_id);

				my $ft2 = Bio::EnsEMBL::SeqFeature->new(
					-seqname   => $q_id,
					-start     => $sf->{query_start},
					-end       => $sf->{query_end},
					-strand    => $q_strand,
					-coverage  => $coverage);

				my $feature_pair = Bio::EnsEMBL::FeaturePair->new(
					-feature1    => $ft1,
					-feature2      => $ft2,
				);

				push @exon_feature_pairs, $feature_pair;
				push @tran_feature_pairs, $feature_pair;
			}

			if (@exon_feature_pairs) {

				# Use our feature pairs for this exon to create a single
				# supporting feature (with cigar line).
				my $supp_feature;

				eval {
					if ( $self->query_type eq 'protein' )
					{
						$supp_feature =
						  Bio::EnsEMBL::DnaPepAlignFeature->new(
							-features => \@exon_feature_pairs );
					}
					else {
						$supp_feature =
						  Bio::EnsEMBL::DnaDnaAlignFeature->new(
							-features => \@exon_feature_pairs );
					}
				};
				if ($@) {
					$self->warn($@);
					return 0;
				}
				$exon->add_supporting_features($supp_feature);
				push @exon_sup_features,$supp_feature;
			}

			$transcript->add_Exon($exon);
		}


# Comment out the following peace of code below because it produces
# the "Insertion Deletion in the same alignment region" warning
# when creating AlignFeature objects from protein alignment with
# slip codons.


		# Create a single supporting feature for the whole transcript
#		my $t_supp_feat;
#		eval {
#			if ( $self->query_type eq 'protein' )
#			{
#				$t_supp_feat =
#				  Bio::EnsEMBL::DnaPepAlignFeature->new(
#					-features => \@tran_feature_pairs );
#			}
#			else {
#				$t_supp_feat =
#				  Bio::EnsEMBL::DnaDnaAlignFeature->new(
#					-features => \@tran_feature_pairs );
#			}
#		};
#		if ($@) {
#			$self->warn("Could not create Transcript supporting feature");
#		}
#		else {
#			#$transcript->add_supporting_features($t_supp_feat);
#		}

	}

	return \@exon_sup_features;

}

sub _parse_vulgar_block {
	my (
		$self,          $target_start,  $target_end,
		$target_strand, $target_length, $query_start,
		$query_end,     $query_strand,  $query_length,
		$vulgar_components
	  )
	  = @_;

	# This method works along the length of a vulgar line
	# exon-by-exon.  Matches that comprise an exon are
	# grouped and an array of 'proto-exons' is returned.
	# Coordinates from the vulgar line are extrapolated
	# to actual genomic/query coordinates.

	my @exons;
	my $exon_number = 0;

	# We sometimes need to increment all our start coordinates. Exonerate
	# has a coordinate scheme that counts _between_ nucleotides at the start.
	# However, for reverse strand matches

	my ( $query_in_forward_coords, $target_in_forward_coords );
	my ( $cumulative_query_coord,  $cumulative_target_coord );

	if ( $target_start > $target_end ) {
		warn(
"For target, start and end are in thew wrong order for a reverse strand match"
		  )
		  if $target_strand != -1;
		$cumulative_target_coord  = $target_start;
		$target_in_forward_coords = 1;
	}
	else {
		$cumulative_target_coord  = $target_start + 1;
		$target_in_forward_coords = 0;
	}
	if ( $query_start > $query_end ) {
		warn(
"For query, start and end are in thew wrong order for a reverse strand match"
		  )
		  if $query_strand != -1;
		$cumulative_query_coord  = $query_start;
		$query_in_forward_coords = 1;
	}
	else {
		$cumulative_query_coord  = $query_start + 1;
		$query_in_forward_coords = 0;
	}

	while (@$vulgar_components) {
		$self->throw(  "Something funny has happened to the input vulgar string."
			  . "  Expecting components in multiples of three, but only have ["
			  . scalar @$vulgar_components
			  . "] items left to process." )
		  unless scalar @$vulgar_components >= 3;

		my $type                = shift @$vulgar_components;
		my $query_match_length  = shift @$vulgar_components;
		my $target_match_length = shift @$vulgar_components;

		$self->throw(  "Vulgar string does not start with a match.  Was not "
			  . "expecting this." )
		  if ( scalar @exons == 0 ) && $type ne 'M' && $type ne 'S';

		if ( $type eq 'M' or $type eq 'S' ) {
			my %hash;

			if ( $target_strand == -1 ) {
				if ($target_in_forward_coords) {
					$hash{target_start} =
					  $cumulative_target_coord - ( $target_match_length - 1 );
					$hash{target_end} = $cumulative_target_coord;
				}
				else {
					$hash{target_end} =
					  $target_length - ( $cumulative_target_coord - 1 );
					$hash{target_start} =
					  $hash{target_end} - ( $target_match_length - 1 );
				}
			}
			else {
				$hash{target_start} = $cumulative_target_coord;
				$hash{target_end}   =
				  $cumulative_target_coord + ( $target_match_length - 1 );
			}

			if ( $query_strand == -1 ) {
				if ($query_in_forward_coords) {
					$hash{query_start} =
					  $cumulative_query_coord - ( $query_match_length - 1 );
					$hash{query_end} = $cumulative_query_coord;
				}
				else {
					$hash{query_end} =
					  $query_length - ( $cumulative_query_coord - 1 );
					$hash{query_start} =
					  $hash{query_end} - ( $query_match_length - 1 );
				}
			}
			else {
				$hash{query_start} = $cumulative_query_coord;
				$hash{query_end}   =
				  $cumulative_query_coord + ( $query_match_length - 1 );
			}

			if ( $type eq 'S' ) {
				$hash{incomplete_codon} = $target_match_length;
			}

			# there is nothing to add if this is the last state of the exon
			$exons[$exon_number]->{gap_end} = 0;
			push @{ $exons[$exon_number]->{sf} }, \%hash;
		}
		elsif ( $type eq "G" ) {
			if ( exists( $exons[$exon_number]->{sf} ) ) {

		 # this is the gap in the middle of an exon, or at the end. Assume it is
		 # at the end, and then reset if we see another match state in this exon
				$exons[$exon_number]->{gap_end} = $target_match_length;
			}
			else {

				# this is a gap at the start of an exon;
				$exons[$exon_number]->{gap_start} = $target_match_length;
			}
		}
		elsif ($type eq "I"
			or $type eq "F" )
		{

	   # in protein mode, any insertion on the genomic side should be treated as
	   # an intron to ensure that the result translates. However, we allow for
	   # codon insertions in the genomic sequence with respect to the protein.
	   # This introduces the possibility of in-frame stops, but I don't
	   # think "introning over" these insertions is appropriate here.

# if we see a gap/intron immediately after an intron, the current exon is "empty"
			if ( $exons[$exon_number] ) {
				$exon_number++;
			}
		}

		if ( $target_in_forward_coords and $target_strand == -1 ) {
			$cumulative_target_coord -= $target_match_length;
		}
		else {
			$cumulative_target_coord += $target_match_length;
		}
		if ( $query_in_forward_coords and $query_strand == -1 ) {
			$cumulative_query_coord -= $query_match_length;
		}
		else {
			$cumulative_query_coord += $query_match_length;
		}

	}

	for ( my $i = 0 ; $i < @exons ; $i++ ) {
		my $ex    = $exons[$i];
		my $ex_sf = $ex->{sf};

		if ( $target_strand == -1 ) {
			$ex->{exon_start} = $ex_sf->[-1]->{target_start};
			$ex->{exon_end}   = $ex_sf->[0]->{target_end};

			if ( exists $ex->{gap_start} ) {
				$ex->{exon_end} += $ex->{gap_start};
			}
			if ( exists $ex->{gap_end} ) {
				$ex->{exon_start} -= $ex->{gap_end};
			}

		}
		else {
			$ex->{exon_start} = $ex_sf->[0]->{target_start};
			$ex->{exon_end}   = $ex_sf->[-1]->{target_end};

			if ( exists $ex->{gap_start} ) {
				$ex->{exon_start} -= $ex->{gap_start};
			}
			if ( exists $ex->{gap_end} ) {
				$ex->{exon_end} += $ex->{gap_end};
			}
		}

		# filter out the incomplete supporting features
		my @sfs;
		foreach my $f (@$ex_sf) {
			if ( not exists( $f->{incomplete_codon} ) ) {
				push @sfs, $f;
			}
		}
		$ex->{sf} = \@sfs;
	}

	return \@exons;
}



1;

=head1 NAME - Bio::EnsEMBL::Pipeline::Runnable::Finished_Exonerate

=head1 AUTHOR

Mustapha Larbaoui B<email> ml6@sanger.ac.uk
