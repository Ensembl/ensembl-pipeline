#!/bin/sh

BLASTDB="/data/blastdb/Ensembl"
BLASTMAT="/usr/local/ensembl/data/blastmat"
BLASTFILTER="/usr/local/ensembl/bin"

export BLASTDB BLASTMAT BLASTFILTER

root_dir="/ecs4/work5/finished/production"

PERL5LIB="$root_dir/ensembl-pipeline/modules\
:$root_dir/ensembl/modules\
:$root_dir/bioperl-0.7.2\
:$root_dir/bioperl-1.2\
:$root_dir/PerlModules"

export PERL5LIB

exec ./Finished_RuleManager.pl \
-dbhost otterpipe1 \
-dbport 3302 \
-dbuser ottadmin \
-dbpass lutralutra \
-start_from SubmitContig \
$@
