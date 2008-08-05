--
-- Tables structures for efg analysis pipeline
--

DROP TABLE IF EXISTS `meta`;
CREATE TABLE meta (

  meta_id                     INT NOT NULL AUTO_INCREMENT,
  species_id                  INT UNSIGNED DEFAULT 1,
  meta_key                    VARCHAR(40) NOT NULL,
  meta_value                  VARCHAR(255) NOT NULL,

  PRIMARY   KEY (meta_id),
  UNIQUE    KEY species_key_value_idx (species_id, meta_key, meta_value),
            KEY species_value_idx (species_id, meta_value)

) COLLATE=latin1_swedish_ci TYPE=MyISAM;

DROP TABLE IF EXISTS `analysis`;
CREATE TABLE `analysis` (
  `analysis_id` int(10) unsigned NOT NULL auto_increment,
  `created` datetime NOT NULL default '0000-00-00 00:00:00',
  `logic_name` varchar(100) NOT NULL default '',
  `db` varchar(120) default NULL,
  `db_version` varchar(40) default NULL,
  `db_file` varchar(120) default NULL,
  `program` varchar(80) default NULL,
  `program_version` varchar(40) default NULL,
  `program_file` varchar(80) default NULL,
  `parameters` varchar(255) default NULL,
  `module` varchar(80) default NULL,
  `module_version` varchar(40) default NULL,
  `gff_source` varchar(40) default NULL,
  `gff_feature` varchar(40) default NULL,
  PRIMARY KEY  (`analysis_id`),
  UNIQUE KEY `logic_name` (`logic_name`),
  KEY `logic_name_idx` (`logic_name`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

DROP TABLE IF EXISTS `analysis_description`;
CREATE TABLE analysis_description (

  analysis_id	          SMALLINT UNSIGNED NOT NULL,
  description             TEXT,
  display_label           VARCHAR(255),
  displayable             BOOLEAN NOT NULL DEFAULT 1,
  web_data                TEXT,

  UNIQUE KEY analysis_idx (analysis_id)
  
) COLLATE=latin1_swedish_ci TYPE=MyISAM;

DROP TABLE IF EXISTS `job`;
CREATE TABLE job (
  job_id            int(10) unsigned NOT NULL auto_increment,
  input_id          varchar(100) NOT NULL,
  analysis_id       smallint(5) unsigned NOT NULL,
  submission_id     mediumint(10) unsigned NOT NULL,
  stdout_file       varchar(200) NOT NULL,
  stderr_file       varchar(200) NOT NULL,
  retry_count       tinyint(2) unsigned default 0,
  temp_dir          varchar(100) DEFAULT '',
  exec_host         varchar(40) DEFAULT '',       
  PRIMARY KEY       (job_id),
  KEY               (input_id),
  KEY               (analysis_id)
);

# job_id        - job internal ID
# input_id      - name (e.g. accession/Ensembl ID) of input
# analysis_id   - internal ID of analysis (analysis table)
# submission_id - ID of job in LSF
# *_file        - files created to contain job output/error
# retry_count   - number of times job restarted

# ?? what is job.objectfile - do we need/use it?



DROP TABLE IF EXISTS `job_status`;
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



# rule tables: rule_goal and rule_conditions
#
# an analysis A depends on a number of rules R1 .. RN. it can
# be run when all rules with rule_goal.analysis = A are satisfied.
# in order to satisfy a rule R, at least one of the conditions
# rule_conditions.condition must have been met for that rule.
#
# rule_goal.analysis and rule_conditions.contitions relate to
# entries in table 'analysis'
# currently, rule_goal.analysis is an analysis internal ID,
# rule_conditions.condition is an analysis logic_name

DROP TABLE IF EXISTS `rule_goal`;
CREATE TABLE rule_goal (
  rule_id           smallint(10) unsigned NOT NULL auto_increment,
  goal              varchar(40),

  PRIMARY KEY (rule_id)
);

DROP TABLE IF EXISTS `rule_conditions`;
CREATE TABLE rule_conditions (
  rule_id           smallint(10) unsigned NOT NULL,
  rule_condition        varchar(40),

  KEY (rule_id)
);

# rule_id     - rule internal ID
# goal - internal ID (if +ve integer) or logic_name of analysis
# condition   - a literal, such as analysis.logic_name

# need also to include rules like 'run job X on all contigs'


DROP TABLE IF EXISTS `input_id_analysis`;
CREATE TABLE input_id_analysis (
  input_id          varchar(100) not null,
  input_id_type     varchar(40) not null,
  analysis_id       smallint(10) unsigned NOT NULL,
  created           datetime NOT NULL,
  runhost           varchar(20) NOT NULL,
  db_version        varchar(40) NOT NULL,
  result            smallint(10) unsigned NOT NULL,

  PRIMARY KEY       (analysis_id, input_id),
  KEY input_created (input_id, created),
  KEY input_type_idx (input_id_type, input_id)
);

# pipeline 'history' table - records each job performed in the
# pipeline with the time it was completed
#
# input_id/analysis_id       - see table 'job'
# timestamp                  - when this job was completed

# ?? need an extra column to hold exit status
# e.g. for repeat masked seq need to know if seq is all repeat
# (prevent e.g. blast running) - and store failures here that
# are understood (e.g. genscan producing garbage)
# and if, say, no genscan predictions, don\'t bother BlastGenscanXXX

DROP TABLE IF EXISTS `input_id_type_analysis`;
CREATE table input_id_type_analysis (
  analysis_id       smallint(10) unsigned NOT NULL,
  input_id_type     varchar(40) NOT NULL
);
