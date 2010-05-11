#! /usr/local/ensembl/bin/perl -w

use strict;
use Getopt::Long 'GetOptions';
use DBI;

my ($dbhost, $dbuser, $dbpass, $dbport, $dbname, $help,$bkp_dir);

$bkp_dir = "/nfs/anacode/protein_pipeline/backups";

GetOptions('dbhost=s' => \$dbhost,
		   'dbuser=s' => \$dbuser,
		   'dbpass=s' => \$dbpass,
		   'dbport=s' => \$dbport,
		   'dbname=s' => \$dbname,
		   'bkp_dir=s'=> \$bkp_dir,
		  );    # plus default options

exec('perldoc', $0) if !($dbhost && $dbuser && $dbpass && $dbport && $dbname);

my $dbh = DBI->connect("DBI:mysql:$dbname:$dbhost:$dbport", $dbuser, $dbpass, {RaiseError => 1, AutoCommit => 0})
        || die "cannot connect to $dbname, $DBI::errstr";

$dbh->debug();

#----- note that the protein_feature table is of type MyISAM, so transaction won't work (ie, backup before doing)
#      Transaction works for InnoDB table

# back up protein_feature table
my $pf_bk = "$bkp_dir/protein_feature_$dbname.$$";

my $chk = system("mysqldump --host=$dbhost --user=$dbuser --pass=$dbpass --port=$dbport $dbname protein_feature > $pf_bk");
if ( $chk != 0 ){
  die "protein_feature table backup failed - exited!";
}
else {
  print "protein_feature table backup successful: [$pf_bk]\n\n";
}

eval {
  # ----- sql query to spot duplications in protein_feature table populated by
  #       ensembl protein annotation pipeline

  my $sql = "select translation_id, seq_start, seq_end, hit_start, hit_end, hit_name, score, protein_feature_id, analysis_id from protein_feature";
  my $sth = $dbh->prepare($sql);
  $sth->execute;

  my %translation_feature;
  my %dup_analysis;
  my @to_delete;

  while (my $row = $sth->fetchrow_arrayref) {
    my $dupl = $row->[0].'_'.$row->[1].'_'.$row->[2].'_'.$row->[3].'_'.$row->[4].'_'.$row->[5].'_'.$row->[6];
    #warn $dupl;
    push(@{$translation_feature{$dupl}},$row->[7]);
    push(@{$dup_analysis{$dupl}},$row->[8]);
  }

  if ( !%translation_feature ) {
    print "Not getting any data! Somethins is wrong\n";
    $sth->finish;
    $dbh->disconnect;
    exit(0)
  }

  else {

    foreach ( keys %translation_feature ) {
      if ( scalar @{$translation_feature{$_}} > 1){
        my @dupl = sort { $a<=>$b } @{$translation_feature{$_}};
        #warn "Analysis: @{$dup_analysis{$_}}\n";

        shift @dupl;  # keep the oldest one
        push(@to_delete, @dupl);  # kick the rest for trashing
      }
    }
  }

  if ( @to_delete ){
    my $to_delete = join(',', @to_delete); # put commas betw. protein_feature_ids for IN () in SQL

    $sql = "DELETE FROM protein_feature WHERE protein_feature_id IN ($to_delete)";
    $dbh->do($sql);
    $dbh->commit;
    print scalar @to_delete, " duplicated protein_feature id deleted successfully\n";
#	warn $to_delete
  }

  else {
    print "No duplication found.\n";
  }
};

if ($@) {
  print "$@\n";  # just print out error, as transaction cannot be done with MyISam tables
}

$dbh->disconnect;



__END__

=head1 NAME delete_duplicates_from_protfeat_table.pl

=head1 DESCRIPTION

Look for protein feature ids in protein feature table of vega databases after running Ensembl protein annotation pipeline.

=head1 AUTHOR

Chao-Kung Chen B<email> ck1@sanger.ac.uk

Modified by ml6@sanger.ac.uk
