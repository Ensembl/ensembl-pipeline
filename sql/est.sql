#
# Table structure for table 'est'
#
CREATE TABLE est (
  est_internal_id   int(10) unsigned DEFAULT '0' NOT NULL auto_increment,
  est_id            varchar(40) DEFAULT '' NOT NULL,
  est_length        int(10) unsigned DEFAULT '0',

  PRIMARY KEY (est_internal_id),
  UNIQUE est_id(est_id)
);

