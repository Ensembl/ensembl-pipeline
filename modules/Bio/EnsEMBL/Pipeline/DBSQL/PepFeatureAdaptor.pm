package Bio::EnsEMBL::Pipeline::DBSQL::PepFeatureAdaptor;


use Bio::EnsEMBL::DBSQL::BaseAdaptor;


@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);



sub store {
  my ($self,@features) = @_;

  foreach my $f (@features) {
    my $query = "insert into feature values(null,'" . 
      $f->seqname  . "'," . 
      $f->start . ","  .
      $f->end   . ",'" . 
      $f->hseqname  . "'," . 
      $f->hstart . "," . 
      $f->hend   . "," .
      $f->score  ." ,'".
      $f->p_value ."',".
      $f->percent_id    .",'" .
      $f->qtrans ."',".
      $f->qpos   .",".
      $f->qres   .",'".
      $f->qchr   ."','".
      $f->htrans ."',".
      $f->hpos   .",".
      $f->hres   .",'".
      $f->hchr   ."',".
      $f->qgenenum .",".
      $f->hgenenum . ")";

    my $sth = $self->prepare($query);
    my $res = $sth->execute;
  }

}

sub get_all_query_gene_id {
  my ($self) = @_;

  my $query = "select distinct qgene from feature";

  my $sth = $self->prepare($query);

  my $res = $sth->execute;

  my @genes;

  while (my $row = $sth->fetchrow_hashref) {
    my $gene = $row->{qgene};
    push(@genes,$gene);
  }

  return @genes;
}



sub fetch_PepFeatures_by_hchr {
  my ($self,$hchr) = @_;
 
  my $query = "select * from feature where hchr = '$hchr'";

  my $sth = $self->prepare($query);

  my $res = $sth->execute;

  my @features;

  while (my $row = $sth->fetchrow_hashref) {
    my $pep = $self->_create_PepFeature_from_hash($row);
    push(@features,$pep);
  }

  return @features;
  
}

sub fetch_PepFeatures_by_query_gene_id {
  my ($self,$gene) = @_;
 
  my $query = "select * from feature where qgene = '$gene'";

  my $sth = $self->prepare($query);

  my $res = $sth->execute;

  my @features;

  while (my $row = $sth->fetchrow_hashref) {
    my $pep = $self->_create_PepFeature_from_hash($row);
    push(@features,$pep);
  }

  return @features;
  
}

sub fetch_PepFeatures_by_hchr_and_pvalue {
  my ($self,$hchr,$pvalue) = @_;
 
  my $query = "select * from feature where hchr = '$hchr' and 1*evalue < $pvalue";

  my $sth = $self->prepare($query);

  my $res = $sth->execute;

  my @features;

  while (my $row = $sth->fetchrow_hashref) {
    my $pep = $self->_create_PepFeature_from_hash($row);
    push(@features,$pep);
  }

  return @features;
  
}

 
sub fetch_PepFeatures_below_pvalue {
  my ($self,$pvalue) = @_;

  my $query = "select * from feature where 1*evalue < $pvalue";

  my $sth = $self->prepare($query);

  my $res = $sth->execute;

  my @features;

  while (my $row = $sth->fetchrow_hashref) {
    my $pep = $self->_create_PepFeature_from_hash($row);
    push(@features,$pep);
  }

  return @features;
  
}




sub _create_PepFeature_from_hash {
  my ($self,$row) = @_;

  my $feat1 = new Bio::EnsEMBL::SeqFeature(-seqname     => $row->{qgene},
					   -start       => $row->{qstart},
					   -end         => $row->{qend},
					   -source_tag  => 'pepfeature',
					   -primary_tag => 'similarity',
					   -score       => $row->{score},
					   -p_value     => $row->{evalue},
					   -percent_id  => $row->{perc_id});

  my $feat2 = new Bio::EnsEMBL::SeqFeature(-seqname     => $row->{hgene},
					   -start       => $row->{hstart},
					   -end         => $row->{hend},
					   -source_tag  => 'pepfeature',
					   -primary_tag => 'similarity',
					   -score       => $row->{score},
					   -p_value     => $row->{evalue},
					   -percent_id  => $row->{perc_id});

  my $pepfeature = new Bio::EnsEMBL::Pipeline::PepFeature(-feature1 => $feat1,
							  -feature2  => $feat2,
							  -qgenenum  => $row->{qgenenum},
							  -qtrans    => $row->{qtrans},
							  -qpos      => $row->{qpos},
							  -qres      => $row->{qres},
							  -qchr      => $row->{qchr},
							  -hgenenum  => $row->{hgenenum},
							  -htrans    => $row->{htrans},
							  -hpos      => $row->{hpos},
							  -hres      => $row->{hres},
							  -hchr      => $row->{hchr});
  
  return $pepfeature;
}

  
