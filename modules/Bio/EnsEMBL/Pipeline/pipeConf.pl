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
    'nfstmp.dir'  => '/nfs/disk100/humpub3/pipeline',
    'DBI.driver'  => 'mysql',
    'dbhost'      => 'ensrv3',
    'dbname'      => 'chr20',
    'dbuser'      => 'ensadmin',
    'queue'       => 'ultra_blast_farm',
    'batchsize'   => 1,
);
}

1;
