use lib 't/pipeline';
use strict;
use warnings;

BEGIN { $| = 1;
	use Test ;
	plan tests => 15;
}


use TestUtils qw(debug test_getter_setter);

use Bio::EnsEMBL::Pipeline::IDSet;

my $list_1 = [1, 2, 3, 4, 5, 6];
my $list_2 = [2, 4, 6, 8, 10, 12];



my $one = Bio::EnsEMBL::Pipeline::IDSet->new(
					     -ID_LIST => $list_1,
				            );

ok($one);					

my $two = Bio::EnsEMBL::Pipeline::IDSet->new(
					     -ID_LIST => $list_2,
				            );



my $intersect = $one->and($two);

ok($intersect);

my @array = @{$intersect->ID_list};

ok(@array);

my @correct = (2, 4, 6);

@correct = sort{$a <=> $b} @correct;
@array = sort{$a <=> $b} @array;
my $value = 1;
for(my $i = 0 ;$i==2;$i++){
   
   if($array[$i] != $correct[$i]){
     $value = 0;
     last;		 
   }	
}

ok($value);

my $union = $one->or($two);


@array = @{$union->ID_list};

ok(@array);

@correct = (1, 2, 3, 4, 5, 6, 8, 10, 12);

@correct = sort{$a <=> $b} @correct;
@array = sort{$a <=> $b} @array;
$value = 1;
for(my $i = 0 ;$i==8;$i++){
   if($array[$i] != $correct[$i]){
     $value = 0;
     last;		 
   }	
}

ok($value);


my $difference = $one->not($two);

@array = @{$difference->ID_list};

ok(@array);

@correct = (1, 3, 5);

@correct = sort{$a <=> $b} @correct;
@array = sort{$a <=> $b} @array;
$value = 1;
for(my $i = 0 ;$i==2;$i++){
   if($array[$i] != $correct[$i]){
     $value = 0;
     last;		 
   }	
}

ok($value);


my $xdiff = $one->xor($two);

@array = @{$xdiff->ID_list};

ok(@array);

@correct = (1, 3, 5, 8, 10, 12);

@correct = sort{$a <=> $b} @correct;
@array = sort{$a <=> $b} @array;
$value = 1;
for(my $i = 0 ;$i==5;$i++){
   if($array[$i] != $correct[$i]){
     $value = 0;
     last;		 
   }	
}

ok($value);

my $element = 14;

ok($two->add_element($element));

@array = @{$two->ID_list};

$value = 0;
foreach my $a(@array){
   if($a == 14){
     $value = 1;	 
   }  
}

ok($value);

ok($two->delete_element($element));

@array = @{$two->ID_list};

$value = 1;
foreach my $a(@array){
   if($a == 14){
     $value = 0;	 
   }  
}

ok($value);

my $number = $one->count;
ok($number==6);
