# MySQL dump 5.13
#
# Host: localhost    Database: analysis_test
#--------------------------------------------------------
# Server version	3.23.12c-alpha

#
# Table structure for table 'object_state'
#
CREATE TABLE object_state (
  object_state_id int(11) NOT NULL auto_increment,
  object_id varchar(40) DEFAULT '' NOT NULL,
  state_description_id int(11) DEFAULT '0' NOT NULL,
  last_change datetime DEFAULT '0000-00-00 00:00:00' NOT NULL,
  inTransition enum('true','false') DEFAULT 'true' NOT NULL,
  PRIMARY KEY (object_state_id),
  KEY stateIdx (state_description_id),
  KEY transitionIdx (inTransition,object_id)
);

#
# Table structure for table 'transition_log'
#
CREATE TABLE transition_log (
  object_state_id int(11) DEFAULT '0' NOT NULL,
  object_id varchar(40) DEFAULT '' NOT NULL,
  state_description_id int(11) DEFAULT '0' NOT NULL,
  start_transit datetime DEFAULT '0000-00-00 00:00:00' NOT NULL,
  finish_transit datetime,
  PRIMARY KEY (object_state_id,start_transit),
  KEY stateIdx (state_description_id)
);

#
# Table structure for table 'object_error'
#
CREATE TABLE object_error (
  object_state_id int(11) DEFAULT '0' NOT NULL,
  error_mesg mediumtext DEFAULT '' NOT NULL,
  PRIMARY KEY (object_state_id)
);

#
# Table structure for table 'state_description'
#
CREATE TABLE state_description (
  state_description_id int(11) NOT NULL auto_increment,
  state_nickname varchar(40) DEFAULT '' NOT NULL,
  object_class varchar(40) DEFAULT '' NOT NULL,
  endState enum('true','false') DEFAULT 'true' NOT NULL,
  transition_module varchar(80),
  needsMonitoring enum('true','false') DEFAULT 'true' NOT NULL,
  PRIMARY KEY (state_description_id),
  KEY monitorIdx (needsMonitoring)
);

