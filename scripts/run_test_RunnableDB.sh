#!/bin/sh

root_dir="/ecs4/work5/finished/production"

PERL5LIB="$root_dir/ensembl-pipeline/modules\
:$root_dir/ensembl/modules\
:$root_dir/bioperl-0.7.2\
:$root_dir/bioperl-1.2\
:$root_dir/PerlModules"

export PERL5LIB

exec ./test_RunnableDB \
-dbhost otterpipe1 \
-dbname zebrafish_finished \
-dbport 3302 \
-dbuser ottadmin \
-dbpass lutralutra \
-input_id AL672185.7.1.47697 \
-logic_name refseq_zebrafish \
-verbose \
$@
