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
    'ControlDB.name' => 'analysis_test',
    'ControlDB.host' => 'ensrv1.sanger.ac.uk',
    'ControlDB.user' => 'ensadmin',
    'EnsEMBL.name' => 'analysis_test',
    'EnsEMBL.host' => 'ensrv1.sanger.ac.uk',
    'EnsEMBL.user' => 'ensadmin',
    'EnsEMBL.pass' => '',
    'DBI.driver' => 'mysql',
    'Contig.Masked.repeatdata' => '/nfs/some/stupid/file',
);
}

1;
