
CREATE TABLE denormalised_gene_dump (
  denorm_gene_dump_id int(10) unsigned NOT NULL auto_increment,
  dumped_gene mediumtext,
  PRIMARY KEY  (denorm_gene_dump_id)
) MAX_ROWS = 4500000  AVG_ROW_LENGTH = 3500 TYPE=MyISAM;


CREATE TABLE denormalised_gene (
  denorm_gene_id int(10) unsigned NOT NULL auto_increment,
  denorm_gene_dump_id int(10) unsigned NOT NULL,
  chr_name varchar(40) NOT NULL default '',
  gene_start int(10) NOT NULL default '0',
  gene_end int(10) NOT NULL default '0',
  type varchar(40) NOT NULL default '',
  PRIMARY KEY  (denorm_gene_id),
  KEY chr_gene_idx (chr_name, gene_start, gene_end)
);


CREATE TABLE gene_dumpedgene (
  denorm_gene_dump_id int(10) unsigned NOT NULL,
  gene_id int(10) unsigned NOT NULL,
  PRIMARY KEY  (denorm_gene_dump_id, gene_id),
  KEY denorm_gene (denorm_gene_dump_id),
  KEY gene (gene_id)
);

