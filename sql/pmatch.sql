# MySQL dump 8.10
#
# Host: ecs1e    Database: mouse_sanger_Oct01
#--------------------------------------------------------
# Server version	3.23.25-beta

#
# Table structure for table 'pmatch_feature'
#

CREATE TABLE pmatch_feature (
  feature_internal_id int(10) unsigned NOT NULL auto_increment,
  protein_internal_id int(10) unsigned DEFAULT '0' NOT NULL,
  chr_name varchar(40) DEFAULT '' NOT NULL,
  start int(10) unsigned DEFAULT '0' NOT NULL,
  end int(10) unsigned DEFAULT '0' NOT NULL,
  coverage double(16,4) DEFAULT '0.0000' NOT NULL,
  PRIMARY KEY (feature_internal_id),
  UNIQUE protein_internal_id (protein_internal_id,chr_name,start,end),
  KEY protein_internal_id_2 (protein_internal_id)
);

# MySQL dump 8.10
#
# Host: ecs1e    Database: mouse_sanger_Oct01
#--------------------------------------------------------
# Server version	3.23.25-beta

#
# Table structure for table 'protein'
#

CREATE TABLE protein (
  protein_internal_id int(10) unsigned NOT NULL auto_increment,
  protein_id varchar(40) DEFAULT '' NOT NULL,
  cdna_id varchar(40) DEFAULT '' NOT NULL,
  PRIMARY KEY (protein_internal_id),
  UNIQUE protein_id (protein_id)
);

