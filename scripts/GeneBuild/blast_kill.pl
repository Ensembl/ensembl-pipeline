#!/usr/local/bin/perl -w

=head1 NAME

  blast_kill.pl

=head1 SYNOPSIS

  blast_kill.pl

=head1 DESCRIPTION

blast_kill.pl removes blast hits to proteins whose IDs are in the
"kill list" of Swissprot IDs to ignore. This should prevent
their use as supporting evidence during the similarity gene
build.

=head1 OPTIONS

Options are to be set in GeneBuild config modules
The important ones for this script are:
    GeneBuild::Databases::GB_DBHOST
    GeneBuild::Databases::GB_DBNAME
    GeneBuild::Databases::GB_DBUSER
    GeneBuild::Databases::GB_DBPASS
    GeneBuild::Scripts::GB_KILL_LIST

=cut

use strict;
use DBI;

use Bio::EnsEMBL::Pipeline::Config::GeneBuild::Databases qw (
							     GB_DBHOST
							     GB_DBNAME
							     GB_DBUSER
							     GB_DBPASS
							     );

use Bio::EnsEMBL::Pipeline::Config::GeneBuild::Scripts qw (
							   GB_KILL_LIST
							  );
my $host = $GB_DBHOST;
my $name = $GB_DBNAME;
my $user = $GB_DBUSER;
my $pass = $GB_DBPASS;
my $kill_list = $GB_KILL_LIST;

sub get_kill_list {
  my @kill_list_arr = ();
  open KILL_LIST_FH, $kill_list or die "can't open $kill_list";
  while (<KILL_LIST_FH>) {
    my @element = split;
    if (scalar(@element) == 0) {        # blank or empty line
      next;
    }
    push @kill_list_arr, $element[0];
  }
  close KILL_LIST_FH or die "file error for $kill_list";
  return \@kill_list_arr;
}


my $dbh = DBI->connect(
                        "DBI:mysql:host=$host;database=$name",
                           $user, $pass, {RaiseError => 1}
                      );

my $kill_list_arr_ref = get_kill_list();
foreach my $id_to_kill (@$kill_list_arr_ref) {
  my $sth = $dbh->prepare("DELETE FROM supporting_feature WHERE hid = \"$id_to_kill\"");
  $sth->execute;
}
$dbh->disconnect;
