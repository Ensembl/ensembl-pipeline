# Table structure for table 'hit_description'
# Table for the Finished Analysis Runnables to store hit information in.
# hit_name         joins to hit_name of dna_align_feature & protein_align_feature tables
# hit_length       sequence length.
# hit_description  DE line from EMBL/TrEMBL/Swissprot format file
# hit_taxon        Taxon ID for sequence
# hit_db column open to change
#

CREATE TABLE hit_description (
  hit_name varchar(40) DEFAULT '' NOT NULL,
  hit_length int(10) unsigned,
  hit_description text,
  hit_taxon int(10) unsigned,
  hit_db enum('EMBL','Swissprot','TrEMBL','Pfam') DEFAULT 'EMBL' NOT NULL,
  PRIMARY KEY (hit_name),
  KEY hit_db (hit_db)
);

# To change from previous version [patch] does not affect column values.
#ALTER TABLE hit_description MODIFY hit_db ENUM('EMBL','Swissprot','TrEMBL','Pfam') DEFAULT 'EMBL' NOT NULL;