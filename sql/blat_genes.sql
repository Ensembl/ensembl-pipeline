CREATE TABLE blat_exon (
  exon_id       int unsigned NOT NULL auto_increment,
  chr_id     int(10) unsigned NOT NULL,            # foreign key, contig:internal_id
  chr_start  int(10) unsigned NOT NULL,                     # start of exon within contig
  chr_end    int(10) unsigned NOT NULL,                     # end of exon within specified contig
  chr_strand tinyint(2) NOT NULL,                  # 1 or -1 depending on the strand of the exon

  phase         tinyint(2) NOT NULL,
  end_phase     tinyint(2) NOT NULL,
  sticky_rank   tinyint DEFAULT '1' NOT NULL,         # see note above
  
  PRIMARY KEY ( exon_id, sticky_rank),
  KEY chr_idx (chr_id, chr_start )
);


CREATE TABLE blat_exon_transcript (
  exon_id          INT unsigned NOT NULL, # foreign key exon:exon_id
  transcript_id    INT unsigned NOT NULL, # foregin key transcript:transcript_id
  rank          int(10) NOT NULL,         # Indicates the 5' to 3' position of the exon
                                          # within the transcript ie rank of 1 means
                                          # the exon is the 5' most within this transcript
  PRIMARY KEY (exon_id,transcript_id,rank),
  KEY transcript (transcript_id),
  KEY exon ( exon_id )
);


CREATE TABLE blat_dna_align_feature (
  dna_align_feature_id int unsigned not null auto_increment,
  chr_id varchar(10) NOT NULL,
  chr_start int(10) unsigned NOT NULL,
  chr_end int(10) unsigned NOT NULL,
  chr_strand tinyint(1) NOT NULL,
  hit_start int NOT NULL,
  hit_end int NOT NULL,
  hit_strand tinyint(1) NOT NULL,
  hit_name varchar(40) NOT NULL,
  analysis_id int(10) unsigned NOT NULL,

#  What scoring do we need ?

  score double,
  evalue double,
  perc_ident float,
  cigar_line text,

  PRIMARY KEY ( dna_align_feature_id ),
  KEY hit_idx( hit_name ),
  KEY ctg_idx( chr_id, chr_start ),
  KEY ana_idx( analysis_id )
) MAX_ROWS=100000000 AVG_ROW_LENGTH=80;


CREATE TABLE blat_gene (
  gene_id   int unsigned NOT NULL auto_increment,
  type VARCHAR(40) NOT NULL,
  analysis_id int,
  transcript_count int NOT NULL,
  display_xref_id int unsigned NOT NULL,

  PRIMARY KEY (gene_id),
  KEY xref_id_index ( display_xref_id )
);


CREATE TABLE blat_supporting_feature (
  exon_id int(11) DEFAULT '0' NOT NULL,
  feature_type enum('dna_align_feature','protein_align_feature'),
  feature_id int(11) DEFAULT '0' NOT NULL,
  UNIQUE all_idx (exon_id,feature_type,feature_id),
  KEY feature_idx (feature_type,feature_id)
) MAX_ROWS=100000000 AVG_ROW_LENGTH=80;

 

CREATE TABLE blat_transcript (
  transcript_id    INT UNSIGNED NOT NULL auto_increment,  
  gene_id          INT UNSIGNED NOT NULL,          # foreign key gene:gene_id
  translation_id   INT UNSIGNED NOT NULL,          # foreign key translation:translation_id
  exon_count int NOT NULL,
  display_xref_id int unsigned NOT NULL,

  PRIMARY KEY (transcript_id),
  KEY gene_index (gene_id),
  KEY translation_index ( translation_id ),
  KEY xref_id_index ( display_xref_id )

);


