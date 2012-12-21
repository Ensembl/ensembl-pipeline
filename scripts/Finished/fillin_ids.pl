#!/usr/bin/env perl
# $Source: /tmp/ENSCOPY-ENSEMBL-PIPELINE/scripts/Finished/fillin_ids.pl,v $
# $Revision: 1.2 $
#
# fill ids into input_id_analysis table
# 17.10.2002 Kerstin Jekosch [kj2@sanger.ac.uk]
# changed 14.10.2003 Mario Caccamo [mc2@sanger.ac.uk]
# changed 20.10.2003 to allow for partly filled in tables by Kerstin
# Added to cvs the 11st of May 2010 by ml6
# This script is used during the preparation of the Vega database before running the protein pipeline

use warnings ;
use strict;
use DBI;
use Getopt::Long 'GetOptions';
use vars qw($dbname $dbhost $dbuser $dbpass $dbport $fillid $type $old $translation $file $contig);

GetOptions(	'dbname=s'	  => \$dbname,
			'dbhost=s'	  => \$dbhost,
			'dbuser=s'	  => \$dbuser,
			'dbpass=s'    => \$dbpass,
			'dbport=s'    => \$dbport,
			'id=s'        => \$fillid,
			'type=s'      => \$type,
			'old'         => \$old,
			'translation' => \$translation,
			'file:s'      => \$file,
			'contig'      => \$contig,
		  )
	or die "Error processing command line";
	
$dbpass = ''      unless $dbpass;
$dbuser = 'ensro' unless $dbuser;
$fillid = 3       unless $fillid;
$type = 'CONTIG'  unless $type;
$type = 'TRANSLATIONID' if ($translation);
$type = 'FILENAME' if ($file);

my @fillids = split /,/, $fillid;

#######################
# connect to database #
#######################

my $dbh=DBI->connect("DBI:mysql:$dbname:$dbhost:$dbport", $dbuser, $dbpass, {RaiseError => 1}) 
	|| die "cannot connect to db, $DBI::errstr";


$dbh->debug();

##########################
# get ids for the clones #
##########################

my @ids;
my $sth1;
if ($file) {
    my @a = glob("$file/chunk*");
    @ids  = map /\/(chunk\..*)/, @a;
}

elsif ($translation) {
    $sth1 = $dbh->prepare ( q{SELECT translation_id FROM translation } );
}
elsif ($contig) {
    $sth1 = $dbh->prepare ( q{SELECT name FROM contig } );
}


unless ($file) {
    $sth1->execute;
    while (my $row = $sth1->fetchrow_arrayref) {
	    my @f = @$row;
	    my $name = $f[0];
	    push @ids, $name;	
    }
    $sth1->finish;
}

#######################################
# insert into input_id_analysis table #
#######################################

# check whether ids already present
my $sth2 = $dbh->prepare ( q{SELECT input_id FROM input_id_analysis WHERE analysis_id = ?} );
$sth2->execute($fillid);
my %not;
while (my $row = $sth2->fetchrow_arrayref) {
    my @f = @$row;
    my $there = $f[0];
    $not{$there}++;
    print "let out $there\n";
}    
$sth2->finish;

my $sth3;
foreach my $name (@ids) {
    foreach my $id (@fillids) {
        if ($old) { # old schema
            die "can't deal with translation and old schema\n" if ($translation);
            $sth3 = $dbh->prepare ( q{INSERT IGNORE INTO input_id_analysis VALUES (?,?,NOW(),0) } ); 
        }
        else { # new schema
            if ($translation) {
                $sth3 = $dbh->prepare ( q{INSERT IGNORE INTO input_id_analysis VALUES (?,'TRANSLATIONID',?,NOW(),'','',0) } ); 
	    }
            elsif ($contig) {
                $sth3 = $dbh->prepare ( q{INSERT IGNORE INTO input_id_analysis VALUES (?,'CONTIG',?,NOW(),'','',0) } ); 
	    }
            elsif ($file) {
                $sth3 = $dbh->prepare ( q{INSERT IGNORE INTO input_id_analysis VALUES (?,'FILENAME',?,NOW(),'','',0) } );
            }
        }
        #print "execute fillin for $name, $id\n" unless (exists $not{$name});
#        print "INSERT INTO input_id_analysis VALUES ($name,'TRANSLATIONID',$id,NOW(),'','',0) \n" if ($translation);
#        print "INSERT INTO input_id_analysis VALUES ('$name','FILENAME',$id,NOW(),'','',0) \n" if ($file);
#        print "INSERT INTO input_id_analysis VALUES ($name,'CONTIG',$id,NOW(),'','',0) \n" if ($contig);
        $sth3->execute($name,$id) unless (exists $not{$name});
    }
}
$sth3->finish;

############################################
# insert into input_id_type_analysis table #
############################################

unless ($old || $translation || $file) { # new schema
    my $sth4 = $dbh->prepare ( q{SELECT analysis_id FROM analysis } );
    $sth4->execute;
    while (my ($analysis_id) = $sth4->fetchrow) {
        my $sth5 = $dbh->prepare ( q{INSERT IGNORE INTO input_id_type_analysis VALUES (?,?)});
        $sth5->execute($analysis_id,$type);
        $sth5->finish;
    }
    $sth4->finish;
}



