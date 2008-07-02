package Bio::EnsEMBL::Pipeline::SeqFetcher::Finished_Dfetch;

use strict;
use DBI;
use Bio::EnsEMBL::Pipeline::SeqFetcher::xdget;
use Bio::EnsEMBL::Analysis::Config::General;
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Exception qw(throw warning info);


my $uniprot_archive = 'uniprot_archive';
my $verbose = 0;

sub new {
    my ( $class, @args ) = @_;
    my $self = bless {}, $class;

	my ( $analysis , $type , $light , $mm_host , $mm_port , $mm_user ) = rearrange(
        [  'ANALYSIS' , 'TYPE' , 'LIGHT' , 'MM_HOST', 'MM_PORT','MM_USER' ], @args );
    $self->analysis($analysis) if $analysis;
    $self->is_light_fetch($light) if $light;
    throw("Must pass a molecule type, either protein or dna ".
          "not a ".$type) unless($type eq 'protein' || $type eq 'dna');
    $self->type($type);
    $self->mm_host($mm_host || 'cbi3');
    $self->mm_port($mm_port || 3306);
    $self->mm_user($mm_user || 'genero');

    return $self;
}


sub write_descriptions {
    my( $self, $dbobj, $id_list, $chunk_size ) = @_;

    $chunk_size ||= 1000;
	my ($descriptions,$failed) = $self->fetch_descriptions($id_list, $chunk_size);
    my $sth = $self->prepare_hit_desc_sth($dbobj);
	for my $accession (keys %$descriptions) {
		my @desc = map { $descriptions->{$accession}{$_}; } ('description','length','taxonid','database');
		eval{
			$sth->execute(
				$accession,
				@desc
			);
		};
		if ($@) {
            print STDERR "Unable to fetch description for $accession [$@]\n";
        }
	}
    return $failed;
}

sub fetch_descriptions {
    my( $self, $id_list, $chunk_size ) = @_;

	my $descriptions = {};	# results are stored here
	$chunk_size ||= 1000;
	my $failed = [];
	my $f;
	my $fetch_all = 1;
	my $query;
	my $dbs;

	# Get sequence length and description from local XDF database
	# if an Analysis object is set
	if($self->analysis) {
		$self->is_light_fetch(1);
		print STDOUT "Fetching ".@$id_list." sequences from XDF database\n" if $verbose;
		SEQ :
		foreach (@$id_list) {
			my $seq = $self->get_Seq_by_acc($_);
			if(!$seq){
				push @$failed,$_;
				next SEQ;
			}
			$descriptions->{$_}{description} = $seq->description();
			$descriptions->{$_}{length} = $seq->length();
		}
	}


	# Fetch data from M&M database
	# full sql query to get all sequence description: length, data class, description, taxon. id
	my $full_query = qq{
		SELECT e.accession_version, e.sequence_length,e.data_class, d.description, t.ncbi_tax_id
		FROM entry e, description d, taxonomy t
		WHERE e.entry_id = d.entry_id
		AND e.entry_id = t.entry_id
		AND e.accession_version IN ('};
	my $full_arch_query = qq{
		SELECT e.accession_version, e.sequence_length,MAX(r.data_class) AS data_class, d.description, t.ncbi_tax_id
		FROM entry e, db_release r,description d, taxonomy t
		WHERE e.entry_id = d.entry_id
		AND e.entry_id = t.entry_id
		AND e.entry_id = r.entry_id
		AND e.accession_version IN ('};
	my $iso_full_query = qq{
		SELECT i_e.accession_version, i_e.sequence_length,MAX(r.data_class) AS data_class, d.description, t.ncbi_tax_id
		FROM entry i_e, entry  p_e, db_release r, description d, taxonomy t, isoform i
		WHERE i_e.entry_id = d.entry_id
		AND i_e.entry_id = i.isoform_entry_id
		AND p_e.entry_id = i.parent_entry_id
		AND p_e.entry_id = t.entry_id
		AND p_e.entry_id = r.entry_id
		AND i_e.accession_version IN ('};
	# light sql query to get the data class and taxon. id only
	my $light_query = qq{
		SELECT e.accession_version, e.data_class, t.ncbi_tax_id
		FROM entry e, taxonomy t
		WHERE e.entry_id = t.entry_id
		AND e.accession_version IN ('};
	my $light_arch_query = qq{
		SELECT e.accession_version, MAX(r.data_class) AS data_class, t.ncbi_tax_id
		FROM entry e, db_release r, taxonomy t
		WHERE e.entry_id = t.entry_id
		AND e.entry_id = r.entry_id
		AND e.accession_version IN ('};
	my $iso_light_query = qq{
		SELECT i_e.accession_version, MAX(r.data_class) AS data_class, t.ncbi_tax_id
		FROM entry i_e, entry  p_e, db_release r, taxonomy t, isoform i
		WHERE i_e.entry_id = i.isoform_entry_id
		AND p_e.entry_id = i.parent_entry_id
		AND p_e.entry_id = t.entry_id
		AND p_e.entry_id = r.entry_id
		AND i_e.accession_version IN ('};

	if($self->type eq 'protein') {
		# get isoform ids
		my $ids_iso = [];
		my $ids = [];
		foreach (@$id_list) {
			if(/\-\d+/) {
				push @$ids_iso, $_ ;
			} else {
				push @$ids, $_ ;
			}
		}
		# fetch protein descriptions from uniprot db (query most recent first)
		$dbs = $self->get_dbnames_like("uniprot%");
		$query = $self->is_light_fetch() ? $light_query : $full_query ;
		$f = $ids;
		while(my $db = shift @$dbs)
		{
			$f = $self->fetch_mm_data($db,$query,0,$f, $chunk_size,$descriptions);
			last unless @$f;
		}
		# fetch remaining protein description from uniprot archive
		$query = $self->is_light_fetch() ? $light_arch_query : $full_arch_query ;
		$f = $self->fetch_mm_data($uniprot_archive,$query,1,$f, $chunk_size,$descriptions) if @$f;

		push @$failed, @$f;

		# fetch isoform protein descriptions from uniprot archive db
		$query = $self->is_light_fetch() ? $iso_light_query : $iso_full_query ;
		$f = $self->fetch_mm_data($uniprot_archive,$query,1,$ids_iso, $chunk_size,$descriptions);
		push @$failed, @$f;


	} else {
		$dbs = $self->get_dbnames_like("embl%");
		$query = $self->is_light_fetch() ? $light_query : $full_query ;
		$f = $id_list;
		while(my $db = shift @$dbs)
		{
			$f = $self->fetch_mm_data($db,$query,0,$f, $chunk_size,$descriptions);
			last unless @$f;
		}
		push @$failed, @$f;
	}

	# remove the failed id from the description hash
	foreach (@$failed) {
		delete $descriptions->{$_};
	}

	print STDOUT "Failed to fetch ".@$failed." ids\n" if @$failed && $verbose;


	return ($descriptions, $failed);
}

sub fetch_mm_data {
	my ($self,$dbname,$sql,$group,$ids, $chunk_size,$descriptions) = @_;

	my $dsn_mole = "DBI:mysql:host=".$self->mm_host.":port=".$self->mm_port.":".$dbname;
	my $dbh = DBI->connect( $dsn_mole, $self->mm_user, '',{ 'RaiseError' => 1 } );
	my $failed = [];

	for (my $i = 0; $i < @$ids; $i += $chunk_size) {
        my $j = $i + $chunk_size - 1;
        # Set second index to last element if we're off the end of the array
        $j = $#$ids if $#$ids < $j;
        # Take a slice from the array
        my $chunk = [@$ids[$i..$j]];
        my $query = $sql . join("','",@$chunk)."')";
        $query .= " GROUP BY accession_version" if $group;
        my $sth;

        # create a list of failed ids
        my %failed_list = map {$_ , 1 } @$chunk;

		print STDOUT "Fetching data from M&M database $dbname for ".@$chunk." ids\n" if $verbose;

        $sth = $dbh->prepare($query);
        $sth->execute() or die "Couldn't execute statement: " . $sth->errstr;
        while(my $hash = $sth->fetchrow_hashref()) {

			my $id = $hash->{'accession_version'};
			# remove id from failed list
			delete $failed_list{$id};

			my $hit_db = $hash->{'data_class'};
			if($dbname =~ /uniprot/i) {
				if($hit_db eq 'STD') {
					#print STDOUT "$hit_db => Swissprot\n";
					$hit_db = 'Swissprot';
				} elsif ($hit_db eq 'PRE') {
					#print STDOUT "$hit_db => TrEMBL\n";
					$hit_db = 'TrEMBL';
				}
			} else {
				#print STDOUT "$hit_db => EMBL\n";
				$hit_db = 'EMBL';
			}
			$descriptions->{$id}{database} = $hit_db;
			$descriptions->{$id}{taxonid} = $hash->{'ncbi_tax_id'};

			if(!$self->is_light_fetch()) {
				$descriptions->{$id}{description} = $hash->{'description'};
		    	$descriptions->{$id}{length} = $hash->{'sequence_length'};
			}
		}
		$sth->finish();

        push @$failed, keys %failed_list;
    }

    return $failed;
}


sub get_dbnames_like {
	my ($self, $db) = @_;

	my $ini_db = 'mm_ini';
	my @dbs;

	my $dsn_mole_ini = "DBI:mysql:host=".$self->mm_host.":port=".$self->mm_port.":".$ini_db;
	my $dbh_ini = DBI->connect( $dsn_mole_ini, $self->mm_user, '',{ 'RaiseError' => 1 } );

	my $query = "SELECT database_name
				 FROM ini
				 WHERE available = 'yes'
				 AND database_category LIKE \'$db\'
				 ORDER BY installation DESC";
	my $ini_sth = $dbh_ini->prepare($query);
	$ini_sth->execute() or die "Couldn't execute statement: " . $ini_sth->errstr;

	while(my $hash = $ini_sth->fetchrow_hashref()) {
		push @dbs, $hash->{'database_name'};
	}
	$ini_sth->finish();

	return \@dbs;
}


sub  get_Seq_by_acc {
  my ($self, $acc) = @_;

  if (!$self->analysis) {
    throw("Analysis object not set");
  }

  if (!defined($acc)) {
    throw("No accession input");
  }

  my $seq;
  my @seqfetchers = $self->_seqfetcher;

  foreach my $seqfetcher (@seqfetchers){
    eval{
      $seq = $seqfetcher->get_Seq_by_acc($acc);
    };
    if ( $@ ){
      warning("problem fetching sequence for $acc");
    }

    if ( defined $seq ){
      last;
    }
  }

  return $seq;
}


sub prepare_hit_desc_sth {
    my( $self, $dbobj ) = @_;

    my $sth = $dbobj->prepare(qq{
        REPLACE INTO hit_description (hit_name
              , hit_description
              , hit_length
              , hit_taxon
              , hit_db)
        VALUES (?,TRIM(?),?,?,?)
        });
    return $self->{'_hit_desc_sth'} = $sth;
}


sub analysis{
  my $self = shift;
  my $analysis = shift;
  if($analysis){
    throw("Must pass Bio::EnsEMBL::Pipeline::SeqFetcher::Finished_DFetch a Bio::EnsEMBL::Analysis".
          "not a ".$analysis) unless($analysis->isa('Bio::EnsEMBL::Analysis'));
    $self->{'analysis'} = $analysis;

    $self->db_version($analysis->db_version) if $analysis->db_version;

	throw("db_file or db must be defined in Analysis object") unless ($analysis->db || $analysis->db_file);
	my $db_files = $analysis->db || $analysis->db_file;

	# Create the xdget seqfetcher
	$db_files =~ s/\s//g;
	my @databases;

	foreach my $database ( split( ",", $db_files ) ) {
		if ( $database !~ /^\// ){
		      $database = $BLAST_DIR . "/" . $database;
		}

		if ( -f $database ) {
			push( @databases, $database );
		} else {
			my $count = 1;
			my $db_filename;
			while ( -f ( $db_filename = "${database}-${count}" ) ) {
				push( @databases, $db_filename );
				$count++;
			}
		}
	}

	my $seqfetcher = Bio::EnsEMBL::Pipeline::SeqFetcher::xdget->new(
		-db => \@databases,
        -executable => '/software/farm/bin/xdget');
	$self->_seqfetcher($seqfetcher);

  }

  return $self->{'analysis'};
}

sub _seqfetcher{
  my ($self, $fetcher) = @_;
  if ( $fetcher ){
    push( @{ $self->{'_seqfetcher'} }, $fetcher);
  }
  return @{ $self->{'_seqfetcher'} };
}

sub type {
    my ( $self, $type ) = @_;
    if ($type) {
        $self->{'_type'} = $type;
    }
    return $self->{'_type'};
}


sub is_light_fetch {
	my ( $self, $light) = @_;
	if($light){
		$self->{'light'} = $light;
	}

	return $self->{'light'};
}

sub db_version {
    my ( $self, $db_version ) = @_;
    if ($db_version) {
        $self->{'_db_version'} = $db_version;
    }
    return $self->{'_db_version'};
}

sub mm_host {
    my ( $self, $mm_host ) = @_;
    if ($mm_host) {
        $self->{'_mm_host'} = $mm_host;
    }
    return $self->{'_mm_host'};
}

sub mm_port {
    my ( $self, $mm_port ) = @_;
    if ($mm_port) {
        $self->{'_mm_port'} = $mm_port;
    }
    return $self->{'_mm_port'};
}

sub mm_user {
    my ( $self, $mm_user ) = @_;
    if ($mm_user) {
        $self->{'_mm_user'} = $mm_user;
    }
    return $self->{'_mm_user'};
}

1;