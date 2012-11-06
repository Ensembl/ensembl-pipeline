# $Source: /tmp/ENSCOPY-ENSEMBL-PIPELINE/sql/patch_08_2012.sql,v $
# $Revision: 1.1 $
# This patch needs to be applied to reference databases if they used an ensembl-pipeline 
# checkout before August 2012. An info column has been added to the job table.

ALTER TABLE job ADD info varchar(200) DEFAULT '' ;

