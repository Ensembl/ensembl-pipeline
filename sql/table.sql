#
# Table structure for table 'analysisprocess'
#
CREATE TABLE analysisprocess (
  analysisId int(10) unsigned DEFAULT '0' NOT NULL auto_increment,
  created datetime DEFAULT '0000-00-00 00:00:00' NOT NULL,
  logic_name varchar(40) not null,
  db varchar(40),
  db_version varchar(40),
  db_file varchar(80),
  program varchar(40),
  program_version varchar(40),
  program_file varchar(40),
  parameters varchar(80),
  module varchar(80),
  module_version varchar(40),
  gff_source varchar(40),
  gff_feature varchar(40),
  PRIMARY KEY (analysisId)
);

#
# Table structure for table 'current_status'
#
CREATE TABLE current_status (
  jobId int(10) unsigned DEFAULT '0' NOT NULL,
  status varchar(40) DEFAULT '' NOT NULL,
  PRIMARY KEY (jobId),
  KEY status_index (status)
);

#
# Table structure for table 'exon_pair'
#
CREATE TABLE exon_pair (
  exon1_id varchar(40) DEFAULT '' NOT NULL,
  exon2_id varchar(40) DEFAULT '' NOT NULL,
  created datetime DEFAULT '0000-00-00 00:00:00' NOT NULL,
  type varchar(40) DEFAULT '' NOT NULL,
  exon1_version varchar(40) DEFAULT '0' NOT NULL,
  exon2_version varchar(40) DEFAULT '0' NOT NULL,
  PRIMARY KEY (exon1_id,exon2_id,exon1_version,exon2_version,type),
  KEY exon1_index (exon1_id),
  KEY exon2_index (exon2_id)
);

#
# Table structure for table 'job'
#
CREATE TABLE job (
  jobId int(10) unsigned DEFAULT '0' NOT NULL auto_increment,
  input_id varchar(40) DEFAULT '' NOT NULL,
  analysisId int(10) unsigned DEFAULT '0' NOT NULL,
  LSF_id int(10) unsigned DEFAULT '0',
  stdout_file varchar(100) DEFAULT '' NOT NULL,
  stderr_file varchar(100) DEFAULT '' NOT NULL,
  object_file varchar(100) DEFAULT '' NOT NULL,
  retry_count int default 0,

  PRIMARY KEY (jobId),
  KEY input_index (input_id),
  KEY analysis_index (analysisId)
);

#
# Table structure for table 'jobstatus'
#
CREATE TABLE jobstatus (
  jobId int(10) unsigned DEFAULT '0' NOT NULL,
  status varchar(40) DEFAULT 'CREATED' NOT NULL,
  time datetime DEFAULT '0000-00-00 00:00:00' NOT NULL
);

#
# Table structure for table 'timclone'
#
CREATE TABLE timclone (
  disk_id varchar(40) DEFAULT '' NOT NULL,
  clone_group set('SU','EU','SF','EF') DEFAULT '' NOT NULL,
  chromosome varchar(10) DEFAULT 'unk' NOT NULL,
  last_check datetime DEFAULT '0000-00-00 00:00:00' NOT NULL,
  created datetime DEFAULT '0000-00-00 00:00:00' NOT NULL,
  dna_update_state int(10) unsigned DEFAULT '0' NOT NULL,
  update_state int(10) unsigned DEFAULT '2' NOT NULL,
  update_label varchar(40) DEFAULT 'dna_read' NOT NULL,
  update_date datetime DEFAULT '0000-00-00 00:00:00' NOT NULL,
  modified datetime DEFAULT '0000-00-00 00:00:00' NOT NULL,
  internal_lock int(10) DEFAULT '1' NOT NULL,
  external_lock int(10) DEFAULT '1' NOT NULL,
  PRIMARY KEY (disk_id)
);

CREATE TABLE RuleGoal (
  ruleId int unsigned default '0' not null auto_increment,
  goalAnalysisId int unsigned,
 
  PRIMARY KEY (	ruleID )
);

CREATE TABLE RuleConditions (
  ruleId int not null,
  conditionLiteral varchar(20)
);

CREATE TABLE InputIdAnalysis (
  inputId varchar(20) not null,
  class enum( "clone","contig","vc","gene" ) not null,
  analysisId int not null,
  created datetime not null,

  PRIMARY KEY ( analysisId, inputId, class ),
  KEY inputIdx( inputId, created )
);


