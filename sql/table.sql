
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
  job_id            int(10) unsigned NOT NULL auto_increment,
  taskname          varchar(40) NOT NULL,
  input_id          varchar(40) NOT NULL,
  submission_id     varchar(40),
  array_index       varchar(255),
  parameter         varchar(255),
  module            varchar(255),
  stdout_file       varchar(100),
  stderr_file       varchar(100),

  PRIMARY KEY (job_id),
  KEY         (input_id),
  KEY         (taskname)
);



# TABLE job_status - contains a history of the states of all jobs
#                    the current status of a job is the one with the most
#                    recent timestamp
#
# job_id     - job internal ID
# status     - text string (e.g. 'CREATED' , 'RUNNING')


CREATE TABLE job_status (
  job_id   int(10) unsigned NOT NULL,
  status   enum('CREATED', 'SUBMITTED', 'READING', 'WRITING',
                'RUNNING', 'SUCCESFUL', 'FATAL', 'KILLED') NOT NULL,
  time     datetime NOT NULL,

  KEY (job_id),
  KEY (status)
);



CREATE TABLE config (
  header          varchar(255),
  key_name        varchar(255),
  value           varchar(255)
);

