# MySQL dump 5.13
#
# Host: localhost    Database: analysis
#--------------------------------------------------------
# Server version	3.22.22

#
# Table structure for table 'analysis'
#
CREATE TABLE analysis (
  id int(10) unsigned DEFAULT '0' NOT NULL auto_increment,
  created datetime DEFAULT '0000-00-00 00:00:00' NOT NULL,
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
  gff_feature varchar(40) DEFAULT '' NOT NULL,
  PRIMARY KEY (id)
);

#
# Table structure for table 'analysis_list'
#
CREATE TABLE analysis_list (
  set_id int(10) DEFAULT '0' NOT NULL,
  analysis_id int(10) DEFAULT '0' NOT NULL
);

#
# Table structure for table 'analysis_set'
#
CREATE TABLE analysis_set (
  id int(10) DEFAULT '0' NOT NULL auto_increment,
  created datetime DEFAULT '0000-00-00 00:00:00' NOT NULL,
  PRIMARY KEY (id)
);

#
# Table structure for table 'clone'
#
CREATE TABLE clone (
       disk_id varchar(40) NOT NULL,
       clone_group set("SU","EU","SF","EF") NOT NULL,
       chromosome varchar(10) DEFAULT 'unk' NOT NULL,
       last_check datetime NOT NULL,
       created datetime NOT NULL,
       dna_update_state int(10) unsigned DEFAULT '0' NOT NULL,
       update_state int(10) unsigned DEFAULT '2' NOT NULL,
       update_label varchar(40) DEFAULT 'dna_read' NOT NULL,
       update_date datetime DEFAULT 'now()' NOT NULL,
       modified datetime NOT NULL,
       internal_lock int(10) DEFAULT '1' NOT NULL,
       external_lock int(10) DEFAULT '1' NOT NULL,
       PRIMARY KEY (disk_id)
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
  PRIMARY KEY (exon1_id,exon2_id,exon1_version,exon2_version,type)
);

#
# Table structure for table 'job'
#
CREATE TABLE job (
  id int(10) unsigned DEFAULT '0' NOT NULL auto_increment,
  input_id int(10) unsigned DEFAULT '0' NOT NULL,
  analysis int(10) unsigned DEFAULT '0' NOT NULL,
  LSF_id int(10) unsigned DEFAULT '0',
  machine varchar(40) DEFAULT '',
  object mediumtext,
  queue varchar(40) DEFAULT '',
  PRIMARY KEY (id)
);

#
# Table structure for table 'jobstatus'
#
CREATE TABLE jobstatus (
  id int(10) unsigned DEFAULT '0' NOT NULL,
  status varchar(40) DEFAULT 'CREATED' NOT NULL,
  time datetime DEFAULT '0000-00-00 00:00:00' NOT NULL
);

