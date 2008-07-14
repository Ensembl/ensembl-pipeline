package Bio::EnsEMBL::Pipeline::Tools::MM_Taxonomy;

use strict;
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use DBI;

# Use this module to fetch information from the Mole&Mushroom taxonomy tree

sub new {
    my ( $class, @args ) = @_;
    my $self = bless {}, $class;

	my ( $mm_host , $mm_port , $mm_user, $mm_name ) = rearrange(
        [  'MM_HOST', 'MM_PORT','MM_USER', 'MM_NAME' ], @args );
    $self->mm_host($mm_host || 'cbi3');
    $self->mm_port($mm_port || 3306);
    $self->mm_user($mm_user || 'genero');
    $self->mm_name($mm_name || 'mm_ini');
    $self->init_connection();

    return $self;
}

sub init_connection {
	my ( $self ) = @_;
	my $dsn_mm_ini = "DBI:mysql:host=".$self->mm_host.":".$self->mm_name;
	my $dbh_mm_ini = DBI->connect( $dsn_mm_ini, $self->mm_user, '',{ 'RaiseError' => 1 } );
	my $mushroom_db;
	my $query = "SELECT database_name FROM ini
				 WHERE database_category = 'mushroom'
				 AND current = 'yes'
				 AND available = 'yes'";
	my $ini_sth = $dbh_mm_ini->prepare($query);
	$ini_sth->execute() or die "Couldn't execute statement: " . $ini_sth->errstr;
	while(my $hash = $ini_sth->fetchrow_hashref()) {
		$mushroom_db = $hash->{'database_name'};
	}
	$ini_sth->finish();
	my $dsn_mushroom = "DBI:mysql:host=".$self->mm_host.":$mushroom_db";
	my $dbh_mushroom = DBI->connect( $dsn_mushroom, $self->mm_user, '',{ 'RaiseError' => 1 } );

	$self->dbh($dbh_mushroom);
}

sub get_parent_id {
	my ( $self, $taxon_id ) = @_;
	my $sql = ' SELECT parent_id
				FROM taxonomy
				WHERE ncbi_tax_id=? ';
	my $sth = $self->dbh->prepare($sql);
	$sth->execute($taxon_id) or die "Couldn't execute statement: " . $sth->errstr;;
	my ($parent_id) = $sth->fetchrow_array();

	return $parent_id;
}

# Get all taxon id upstream of the given node id

sub get_all_parent_id {
	my ( $self, $taxon_id ) = @_;
	my $ids = [];
	while(my $id = $self->get_parent_id($taxon_id)) {
		push @$ids, $id;
		$taxon_id = $id;
	}

	return $ids;
}

sub get_children_id {
	my ( $self, $taxon_id ) = @_;
	my $ids = [];
	my $sql =' SELECT ncbi_tax_id
				FROM taxonomy
				WHERE parent_id=? ';
	my $sth = $self->dbh->prepare($sql);
	$sth->execute($taxon_id) or die "Couldn't execute statement: " . $sth->errstr;;
	while (my ($id) = $sth->fetchrow_array()) { push @$ids, $id };

	return $ids;
}

# Get all taxon id downstream of the given node id

sub get_all_children_id {
	my ( $self, $taxon_id ) = @_;
	my $ids = [];
	for my $id (@{$self->get_children_id($taxon_id)}) {
		push @$ids, $id;
		push @$ids, @{$self->get_all_children_id($id)};
	}

	return $ids;

}

# this method does the same thing in a single query but is slower

sub get_all_children {
	my ( $self, $taxon_id ) = @_;
	my $ids = [];
	my $sql =' SELECT  t.ncbi_tax_id
			FROM taxonomy p, taxonomy t
			WHERE p.ncbi_tax_id = ?
			AND p.lft < t.lft
			AND p.rgt > t.rgt';
	my $sth = $self->dbh->prepare($sql);
	$sth->execute($taxon_id) or die "Couldn't execute statement: " . $sth->errstr;;
	while (my ($id) = $sth->fetchrow_array()) { push @$ids, $id };

	return $ids;
}

sub get_taxon_id {
	my ( $self, $taxon_id ) = @_;
	my $sql = ' SELECT n.name, t.rank, t.parent_id
				FROM taxonomy_name n, taxonomy t
				WHERE t.ncbi_tax_id= ?
				AND n.name_type="scientific name"
				AND n.ncbi_tax_id = t.ncbi_tax_id;';
	my $sth = $self->dbh->prepare($sql);
	$sth->execute($taxon_id) or die "Couldn't execute statement: " . $sth->errstr;;
	my $hash = $sth->fetchrow_hashref();

	return $hash;
}

sub print_taxon {
	my ( $self, $tax_arr ) = @_;
	print STDOUT "name\trank\tparent_id\n" if @$tax_arr;
	for (@$tax_arr) {
		my $hash = $self->get_taxon_id($_);
		print STDOUT $hash->{name}."\t".$hash->{rank}."\t".$hash->{parent_id}."\n";
	}
}


sub dbh {
	my ( $self, $dbh ) = @_;
	if ($dbh) {
        $self->{'_dbh'} = $dbh;
    }
    return $self->{'_dbh'};
}

sub mm_name {
    my ( $self, $mm_name ) = @_;
    if ($mm_name) {
        $self->{'_mm_name'} = $mm_name;
    }
    return $self->{'_mm_name'};
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


