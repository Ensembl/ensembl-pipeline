#
# Table structure for table 'current_status'
#
CREATE TABLE current_status (
  jobId             int(10) unsigned DEFAULT '0' NOT NULL,
  status            varchar(40) DEFAULT '' NOT NULL,

  PRIMARY KEY (jobId),
  KEY status_index (status)
);

#
# Table structure for table 'job'
#
CREATE TABLE job (
  jobId             int(10) unsigned DEFAULT '0' NOT NULL auto_increment,
  input_id          varchar(40) DEFAULT '' NOT NULL,
  class             enum ("clone","contig","vc","gene") not null,
  analysisId        int(10) unsigned DEFAULT '0' NOT NULL,
  LSF_id            int(10) unsigned DEFAULT '0',
  stdout_file       varchar(100) DEFAULT '' NOT NULL,
  stderr_file       varchar(100) DEFAULT '' NOT NULL,
  object_file       varchar(100) DEFAULT '' NOT NULL,
  retry_count       int default 0,

  PRIMARY KEY (jobId),
  KEY input_index (input_id),
  KEY analysis_index (analysisId)
);

#
# Table structure for table 'jobstatus'
#
CREATE TABLE jobstatus (
  jobId             int(10) unsigned DEFAULT '0' NOT NULL,
  status            varchar(40) DEFAULT 'CREATED' NOT NULL,
  time              datetime DEFAULT '0000-00-00 00:00:00' NOT NULL,

  KEY (jobId),
  KEY status_index (status)
);

CREATE TABLE RuleGoal (
  ruleId            int unsigned default '0' not null auto_increment,
  goalAnalysisId    int unsigned,
 
  PRIMARY KEY (ruleID)
);

CREATE TABLE RuleConditions (
  ruleId            int not null,
  conditionLiteral  varchar(20)
);

CREATE TABLE InputIdAnalysis (
  inputId           varchar(40) not null,
  class             enum ("clone","contig","vc","gene") not null,
  analysisId        int not null,
  created           datetime not null,

  PRIMARY KEY (analysisId, inputId, class),
  KEY inputIdx      (inputId, created),
  KEY inputId_class (inputId, class)
);

CREATE TABLE VoidInputIdAnalysis (
  inputId           varchar(40) DEFAULT '' NOT NULL,
  analysisId        int(11) DEFAULT '0' NOT NULL,
  exception         varchar(40),

  PRIMARY KEY (inputId,analysisId),
  KEY exc (exception)
);

#
# Table structure for table 'exon_external'
#
CREATE TABLE exon_external (
  internal_id       int(10) unsigned default NULL auto_increment,
  external_id       varchar(40) NOT NULL default '',

  PRIMARY KEY (internal_id)
);

#
# Table structure for table 'gene_external'
#
CREATE TABLE gene_external (
  internal_id       int(10) unsigned default NULL auto_increment,
  external_id       varchar(40) NOT NULL default '',

  PRIMARY KEY (internal_id)
);

#
# Table structure for table 'transcript_external'
#
CREATE TABLE transcript_external (
  internal_id       int(10) unsigned default NULL auto_increment,
  external_id       varchar(40) NOT NULL default '',

  PRIMARY KEY (internal_id)
);

#
# Table structure for table 'translation_external'
#
CREATE TABLE translation_external (
  internal_id       int(10) unsigned default NULL auto_increment,
  external_id       varchar(40) NOT NULL default '',

  PRIMARY KEY (internal_id)
);
