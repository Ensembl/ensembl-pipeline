# Copyright EMBL-EBI 2000
# Author: Arne Stabenau
# Creation: 11.07.2000

# configuration information
# give useful keynames to things

# if states are involved maybe good to have statenames in key.
# parameters should be avail after compile time so
# that you can make changes at runtime.

BEGIN {
package main;

%pipeConf = ( 
    'nfstmp.dir' => '/work1/scp/out',
    'DBI.driver' => 'mysql',
    'dbhost'      => 'ecs1a',
    'dbname'      => 'simon_oct07',
    'dbuser'      => 'ensadmin',
    'queue'       => 'acarichunky',
    'batchsize'   => 50,
);
}

1;
