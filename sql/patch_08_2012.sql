# This patch needs to be applied to reference databases if they used an ensembl-pipeline 
# checkout before August 2012. An info column has been added to the job table.

ALTER TABLE job_status ADD info varchar(200) DEFAULT '' ;

