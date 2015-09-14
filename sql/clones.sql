DROP TABLE clones IF EXISTS;
CREATE TABLE clones (
clones_id INT(11) UNSIGNED NOT NULL AUTO_INCREMENT,
trace_id  BIGINT UNSIGNED NOT NULL,
clone_id VARCHAR(20) NOT NULL,
direction ENUM('R', 'F'),
library VARCHAR(20) NOT NULL,
insert_size  MEDIUMINT UNSIGNED NOT NULL,
insert_stdev MEDIUMINT UNSIGNED NOT NULL,
trace_name   BIGINT UNSIGNED NOT NULL,
length       SMALLINT UNSIGNED NOT NULL,
clip_qleft   SMALLINT UNSIGNED NOT NULL,
clip_qright  SMALLINT UNSIGNED NOT NULL,
clip_vleft   SMALLINT UNSIGNED NOT NULL,
clip_vright  SMALLINT UNSIGNED NOT NULL,
analysis_id  SMALLINT UNSIGNED NOT NULL,
  PRIMARY KEY       (clones_id),
  KEY (trace_id),
  KEY (analysis_id),
  KEY (clones_id)
) ENGINE MyISAM;

# trace_id     id of the sequence, stored in dna_align_feature table by ExonerateAlignFeature
# value        clone name, library,...  CH243-100A1:F:CH243:184000:36800:1098268172037:1001
# analysis_id  Just to add more flexibility...
