create table flag (
flag_id            int(11) unsigned not null auto_increment,
ensembl_id  varchar(20) not null,
table_name  varchar(20) NOT NULL,
analysis_id  int(10) unsigned NOT NULL,
  PRIMARY KEY       (flag_id),
  KEY (ensembl_id),
  KEY (analysis_id)
);
