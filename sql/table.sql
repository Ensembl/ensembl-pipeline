# TABLE job - Stores information about all jobs in the pipeline system
# job_id        - job internal ID
# input_id      - name (e.g. accession/Ensembl ID) of input
# taskname      - name of the task which created this job
# submission_id - the submission identifier provided by the submission system
#                 e.g. LSF submission ID or local PID
# array_index   - the index in the job array that was created, probably only
#                 useful for LSF (if no job array then NULL)
# parameter     - parameters passed on to the module from the runner script
#                 used instead of parameters in the analysis table or 
#                 command line arguments
# module        - the module which this job runs
# *_file        - files created to contain job output/error

CREATE TABLE job (
  job_id            int unsigned NOT NULL auto_increment,
  taskname          varchar(150) NOT NULL,
  input_id          varchar(150) NOT NULL,
  submission_id     int unsigned,
  job_name          varchar(255),
  array_index       mediumint unsigned,
  parameters        varchar(255),
  module            varchar(255),
  stdout_file       varchar(255),
  stderr_file       varchar(255),
  retry_count       tinyint unsigned NOT NULL,

  PRIMARY KEY (job_id),
  UNIQUE      (taskname, input_id),
  KEY         (job_name(15), array_index),
  KEY         (input_id(50)),
  KEY         (taskname(15))
);



# TABLE job_status  
# contains a history of the states of all jobs
# the current status of a job is the one with the 
# highest sequence_num 

CREATE TABLE job_status (
  job_id            int(10) unsigned NOT NULL,
  status            enum('CREATED', 'SUBMITTED', 'READING', 'WRITING',
                         'RUNNING', 'SUCCESSFUL', 'FATAL', 'KILLED', 'FAILED',
                         'RETRIED') NOT NULL,

  time              datetime NOT NULL,
  sequence_num      int unsigned NOT NULL auto_increment, 

  KEY (job_id),
  KEY (sequence_num)
);



CREATE TABLE config (
  header          varchar(255),
  key_name        varchar(255),
  value           varchar(255)
);

