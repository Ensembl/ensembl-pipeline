package Bio::EnsEMBL::Pipeline::SeqFetcher::Finished_Dfetch;

use strict;
use DBI;
use Bio::DB::Flat::OBDAIndex;
use Bio::EnsEMBL::Analysis::Config::General;
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Exception qw(throw warning info);


my $uniprot_db_name = 'uniprot_archive';

sub new {
    my ( $class, @args ) = @_;
    my $self = bless {}, $class;

	my ( $analysis , $type , $mm_host , $mm_port , $mm_user ) = rearrange(
        [  'ANALYSIS' , 'TYPE' , 'MM_HOST', 'MM_PORT','MM_USER' ], @args );
    $self->analysis($analysis) if $analysis;
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

    $chunk_size ||= 200;
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
	my $failed = [];
	my $f;
	my $fetch_all = 1;
	my $query;

	# Get sequence length and description from local OBDA index
	# if an Analysis object is set
	if($self->analysis) {
		$fetch_all = 0;
		print STDOUT "Fetching ".@$id_list." sequences from OBDA index\n";
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
	my $iso_full_query = "	SELECT i_e.accession_version, i_e.sequence_length,p_e.data_class, d.description, t.ncbi_tax_id
							FROM entry i_e, entry  p_e, description d, taxonomy t, isoform i
							WHERE i_e.entry_id = d.entry_id
							AND i_e.entry_id = i.isoform_entry_id
							AND p_e.entry_id = i.parent_entry_id
							AND p_e.entry_id = t.entry_id
							AND i_e.accession_version IN ('";

	my $iso_light_query = " SELECT i_e.accession_version, p_e.data_class, t.ncbi_tax_id
							FROM entry i_e, entry  p_e, taxonomy t, isoform i
							WHERE i_e.entry_id = i.isoform_entry_id
							AND p_e.entry_id = i.parent_entry_id
							AND p_e.entry_id = t.entry_id
							AND i_e.accession_version IN ('";


	my $light_query = "  SELECT e.accession_version, e.data_class, t.ncbi_tax_id
						 FROM entry e, description d, taxonomy t
						 WHERE e.entry_id = d.entry_id
						 AND e.entry_id = t.entry_id
						 AND e.accession_version IN ('";

	my $full_query =  "	 SELECT e.accession_version, e.sequence_length,e.data_class, d.description, t.ncbi_tax_id
						 FROM entry e, description d, taxonomy t
						 WHERE e.entry_id = d.entry_id
						 AND e.entry_id = t.entry_id
						 AND e.accession_version IN ('";

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

		my $dbh = $self->get_dbh_mm_protein;
		$query = $fetch_all ? $full_query : $light_query;
		$f = $self->fetch_mm_data($dbh,$query,$fetch_all,$ids, $chunk_size,$descriptions);
		push @$failed, @$f;

		$query = $fetch_all ? $iso_full_query : $iso_light_query;
		$f = $self->fetch_mm_data($dbh,$query,$fetch_all,$ids_iso, $chunk_size,$descriptions);
		push @$failed, @$f;


	} else {
		my $dbh = $self->get_dbh_mm_dna();
		$query = $fetch_all ? $full_query : $light_query;
		$f = $id_list;
		while(my $con = shift @$dbh)
		{
			$f = $self->fetch_mm_data($con,$query,$fetch_all,$f, $chunk_size,$descriptions);
			last unless @$f;
		}
		push @$failed, @$f;
	}

	# remove the failed id from the description hash
	foreach (@$failed) {
		delete $descriptions->{$_};
	}

	print STDOUT "Failed to fetch ".@$failed." ids\n[@$failed]\n" if @$failed;


	return ($descriptions, $failed);
}

sub fetch_mm_data {
	my ($self,$dbh,$sql,$fetch_all,$ids, $chunk_size,$descriptions) = @_;
	my $failed = [];
	for (my $i = 0; $i < @$ids; $i += $chunk_size) {
        my $j = $i + $chunk_size - 1;
        # Set second index to last element if we're off the end of the array
        $j = $#$ids if $#$ids < $j;
        # Take a slice from the array
        my $chunk = [@$ids[$i..$j]];
        my $sql = $sql.join("','",@$chunk)."')";
        my $sth;

        # create a list of failed ids
        my %failed_list = map {$_ , 1 } @$chunk;

		print STDOUT "Fetching data from M&M database for ".@$chunk." ids\n";

        $sth = $dbh->prepare($sql);
        $sth->execute() or die "Couldn't execute statement: " . $sth->errstr;
        while(my $hash = $sth->fetchrow_hashref()) {

			my $id = $hash->{'accession_version'};
			# remove id from failed list
			delete $failed_list{$id};

			my $hit_db = $hash->{'data_class'};
			if($hit_db eq 'STD') {
				#print STDOUT "$hit_db => Swissprot\n";
				$hit_db = 'Swissprot';
			} elsif ($hit_db eq 'PRE') {
				#print STDOUT "$hit_db => TrEMBL\n";
				$hit_db = 'TrEMBL';
			} else {
				#print STDOUT "$hit_db => EMBL\n";
				$hit_db = 'EMBL';
			}
			$descriptions->{$id}{database} = $hit_db;
			$descriptions->{$id}{taxonid} = $hash->{'ncbi_tax_id'};

			if($fetch_all) {
				$descriptions->{$id}{description} = $hash->{'description'};
		    	$descriptions->{$id}{length} = $hash->{'sequence_length'};
			}
		}
		$sth->finish();

        push @$failed, keys %failed_list;
    }

    return $failed;
}

sub get_dbh_mm_protein {
	my ($self) = @_;
	my $dsn_mole = "DBI:mysql:host=".$self->mm_host.":port=".$self->mm_port.":".$uniprot_db_name;
	my $dbh = DBI->connect( $dsn_mole, $self->mm_user, '',{ 'RaiseError' => 1 } );
	return $dbh;
}

sub get_dbh_mm_dna {
	my ($self) = @_;

	my $ini_db = 'mm_ini';
	my @dbh;
	my @dbs;

	my $dsn_mole_ini = "DBI:mysql:host=".$self->mm_host.":port=".$self->mm_port.":".$ini_db;
	my $dbh_ini = DBI->connect( $dsn_mole_ini, $self->mm_user, '',{ 'RaiseError' => 1 } );

	my $query = "SELECT database_name
				 FROM ini
				 WHERE available = 'yes'
				 AND database_category LIKE 'embl%'
				 ORDER BY installation DESC";
	my $ini_sth = $dbh_ini->prepare($query);
	$ini_sth->execute() or die "Couldn't execute statement: " . $ini_sth->errstr;

	while(my $hash = $ini_sth->fetchrow_hashref()) {
		my $dsn_mole_embl = "DBI:mysql:host=".$self->mm_host.":port=".$self->mm_port.":".$hash->{'database_name'};
		push @dbh, DBI->connect( $dsn_mole_embl, $self->mm_user, '',{ 'RaiseError' => 1 } );
	}
	$ini_sth->finish();

	return \@dbh;
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
  my $have_secondary;
  foreach my $seqfetcher (@seqfetchers){
    $have_secondary = 1 if($seqfetcher->secondary_namespaces);
    eval{
      # note that this only works for OBDAIndex.pm
      $seq = $seqfetcher->get_Seq_by_id($acc);
    };
    if ( $@ ){
      warning("problem fetching sequence for $acc");
    }

    if ( defined $seq ){
      last;
    }
  }

  if(!defined $seq){
    my ($p, $f, $l) = caller;
    warning("OBDAIndexSeqFetcher: could not find sequence for primary key $acc $f:$l\n") if(!$have_secondary);

	FETCHER:
		foreach my $seqfetcher ( $self->_seqfetcher ){
			my @secondary_namespaces = $seqfetcher->secondary_namespaces;
			foreach my $name ( @secondary_namespaces ){
				my @seqs;
				eval{
					# this returns potentially an array of Bio::Seq
					@seqs = $seqfetcher->get_Seq_by_secondary($name,$acc);
				};
				if ( $@ ){
					warning("problem fetching sequence for secondary key $acc $@");
				}
				if ( @seqs > 1 ){
					warning("Multiple sequences (".scalar(@seqs).") for the same secondary accession $acc\n");
					next;
				}

				if ( defined $seqs[0] ){
					$seqs[0]->display_id( $acc );
					$seq = $seqs[0];
					last FETCHER;
				}
			}
		}
		unless ($seq){
			warning("could not find sequence for secondary key $acc");
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

	# Create the OBDA seqfetcher objects from the db index
	$db_files =~ s/\s//g;
	foreach my $database ( split( ",", $db_files ) ){
		if ( $database !~ /^\// ){
		      $database = $BLAST_DIR . "/" . $database;
		}
		if ( $database =~/(\S+)\/$/ ){
			$database = $1;
		}
		my @path = split /\//, $database;
		my $db_name = pop( @path );
		if ( $db_name =~/(\S+)\.fa/){
			$db_name = $1;
		}
		throw("Cannot define db_name") unless ( $db_name );
		my $index_dir = join '/', @path;
		throw("Cannot define index_dir") unless ( $index_dir );
		my $format = 'FASTA';
		#print STDOUT "$index_dir => $db_name\n";
		my $OBDAfetcher = new Bio::DB::Flat::OBDAIndex(
			-index_dir => $index_dir,
			-dbname    => $db_name,
			-format    => $format
		);

	    $self->_seqfetcher($OBDAfetcher);
	}
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