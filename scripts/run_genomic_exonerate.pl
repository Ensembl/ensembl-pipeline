#!/usr/local/bin/perl

use strict;
use Bio::EnsEMBL::Pipeline::RunnableDB::ExonerateGenomic;
use Bio::EnsEMBL::DBSQL::DBAdaptor;

$| = 1;

my $dnadb=Bio::EnsEMBL::DBSQL::DBAdaptor->new(-dbname=>'homo_sapiens_core_110',-host=>'ensrv3',-user=>'ensro');

my $mousedb=Bio::EnsEMBL::DBSQL::DBAdaptor->new(-dbname=>'mouse',-host=>'ensrv3',-user=>'ensadmin');

my $newdb=Bio::EnsEMBL::DBSQL::DBAdaptor->new(-dbname=>'exonerate_genomic',-host=>'ecs1d',-user=>'ensadmin', -dnadb=>$dnadb);

my $estgen = Bio::EnsEMBL::Pipeline::RunnableDB::ExonerateGenomic->new(-querydb => $dnadb, -targetdb => $mousedb, -dbobj => $newdb);
$estgen->fetch_input();
$estgen->run;



