CREATE TABLE job (
  job_id            int(10) unsigned DEFAULT '0' NOT NULL auto_increment,
  input_id          varchar(40) NOT NULL,
  analysis_id       smallint(5) unsigned NOT NULL,
  submission_id     mediumint(10) unsigned NOT NULL,
  stdout_file       varchar(100) NOT NULL,
  stderr_file       varchar(100) NOT NULL,
  retry_count       tinyint(2) unsigned default 0,

  PRIMARY KEY (job_id),
  KEY         (input_id),
  KEY         (analysis_id)
);

# job_id        - job internal ID
# input_id      - name (e.g. accession/Ensembl ID) of input
# analysis_id   - internal ID of analysis (analysis table)
# submission_id - ID of job in LSF
# *_file        - files created to contain job output/error
# retry_count   - number of times job restarted

# ?? what is job.objectfile - do we need/use it?




CREATE TABLE job_status (
  job_id            int(10) unsigned NOT NULL,
  status            varchar(40) DEFAULT 'CREATED' NOT NULL,
  time              datetime NOT NULL,
  is_current        enum('n', 'y') DEFAULT 'n',

  KEY (job_id),
  KEY (status),
  KEY (is_current)
);

# job 'history' table - tracks each state of a job in its 'life'
# one line per job/status
#
# job_id     - job internal ID
# status     - text string (e.g. 'CREATED' , 'RUNNING')
# is_current - whether this status is the current status


CREATE TABLE config (
  header          varchar(255),
  key_name        varchar(255),
  value           varchar(255)
);

