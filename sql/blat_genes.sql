#
# Table structure for table 'blat_gene_dump'
#

CREATE TABLE blat_gene_dump (
  blat_gene_dump_id int(10) unsigned NOT NULL auto_increment,
  chr_name varchar(40) NOT NULL default '',
  gene_start int(10) NOT NULL default '0',
  gene_end int(10) NOT NULL default '0',
  dumped_gene mediumtext,
  type varchar(40) NOT NULL default '',
  mapped_status int(10) unsigned default NULL,
  PRIMARY KEY  (blat_gene_dump_id)
) TYPE=MyISAM;

