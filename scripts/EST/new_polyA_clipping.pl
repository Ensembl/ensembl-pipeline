#!/usr/local/ensembl/bin/perl -w
use strict;

use Bio::Seq;
use Bio::SeqIO;

#script to parse a fasta file and identify sequences with polyA/T tails/heads
#these are then clipped and stored in a file 

#This script is an alternative to polyA_clipping.pl written to enable clipping of those sequences
#which have short non-polyA/T strings at the end of the polyA/T string as well as those which end in A/T strings

#clipping:
#	the non-A/T sequences at the ends must be <=10bp (set by $buffer)  
#	the polyA/T string must be >4bp to be removed
#	it only clips polyA tails or polyT heads using a sliding window of 3 bp
#	the clipping is recursive but only clips one end of a sequence
#	the head/tail is only clipped if the polyA/T string is longer than the non-polyA/T string at the end of the sequence

#perl new_polyA_clipping.pl sequences.fasta polyat_clipped.out

my $data = $ARGV[0];
my $seqoutfile = $ARGV[1];

my $seqout = new Bio::SeqIO(-file => ">$seqoutfile", "-format" => "Fasta");
my $cdna = new Bio::Seq;

my $buffer = 10; #length of region allowed to mismatch at ends of sequence
my $window_size = 3; 
my $clipped_seq = "";

my $end_region = $buffer + $window_size;

local $/ = "\n>";

open(DATA, "<$data") or die ("Can't read $data $! \n");

while(<DATA>){ 
	#have a sequence:
	
	s/>//g;
	
	my ($name, $seq);
	$clipped_seq = undef;
	
	if ($_=~/^([\w\.]+)\s+([\w\s]+)/m){
		$name = $1;
		my $tmp = $2;
		($seq = $tmp)=~s/\s//g; 
	}

	
	if (length $seq >= $end_region){

		#check for poly-A/T tails/heads - which should have mostly been removed...
	
		#decide whether looking for polyT head or polyA tail:
		my $head = substr($seq, 0, $end_region); 
		my $tail = substr($seq, -$end_region); 

		my $t_count = 0;
		my $a_count = 0;
		
		while ($head=~/([Tt]+)/g){ #will match greedily
			if (length $1 > $t_count){
				$t_count = length $1;
			}
		}
		while ($tail=~/([Aa]+)/g){ #will match greedily
			if (length $1 > $a_count){
				$a_count = length $1;
			}
		}
	
		#decide whether to clip heads, tails or unsure... set appropriate flag
		my $clip_end = "";

		if (($a_count > $t_count) && ($a_count > 4)){ #call the subroutine to trim the polyA tail:
			
			$clipped_seq = clip_polya($seq, "tail");
	
		}elsif (($a_count < $t_count) && ($t_count > 4)){ #call the subroutine to trim the polyT head:
			
			$clipped_seq = clip_polya($seq, "head");
		
		}else{
			if ($a_count > 4){ #only do for ones which appear to have a head/tail:
				#tied - not sure which to do -try both and choose afterwards
				my $clipped_head = clip_polya($seq, "head");
				if (!defined $clipped_head){
					$clipped_head = $seq;
				}
			
				my $clipped_tail = clip_polya($seq, "tail");
				if (!defined $clipped_tail){
					$clipped_tail = $seq;
				}
						
				#choose one which clipped the most:
				if (length $clipped_head < length $clipped_tail){
					$clipped_seq = $clipped_head;
				}elsif (length $clipped_tail < length $clipped_head){
					$clipped_seq = $clipped_tail;
				}else{ #still can't tell, leave as original seq
					$clipped_seq = $seq;
				}
			}
		}
	}
	#else{
	#	print "$seq too short to prune\n";
	#}
		
	if (!defined $clipped_seq){
		$clipped_seq = $seq;
	}
	
	eval{
		$cdna->seq($clipped_seq);
 		$cdna->display_id($name);
   		$cdna->desc("");
 	};
	$seqout->write_seq($cdna);
	#print "$name: $seq\n$name: $clipped_seq\n\n";

}
close DATA;

sub clip_polya{
	
	my $seq = shift @_;
	my $end = shift @_;
	
	my @seq = split//, $seq;
	my $length = length $seq;
	my $a_count = 0;
	my $t_count = 0;
	
	
	if ($end eq "tail"){
		my $tail = substr($seq, -$end_region); 
		while ($tail=~/([Aa]+)/g){ #will match greedily
			if (length $1 > $a_count){
				$a_count = length $1;
			}
		}
	}elsif ($end eq "head"){
		my $head = substr($seq, 0, $end_region); 
		while ($head=~/([Tt]+)/g){ #will match greedily
			if (length $1 > $t_count){
				$t_count = length $1;
			}
		}
	}
	
	if($a_count > 4){
		#looking only for polyA tail - use moving window looking for strings of 2/3 As:
	
		#moving through seq starting from end - allow for buffer region:
		for (my $i = ($length - 1); $i > ($length - $buffer); $i--){

			my $match = 0;
			for (my $j = $i; $j > ($i - $window_size); $j--){ #check a window 
				if ($seq[$j] eq 'A'){ 
					$match++;
				}
			}

			if ($match > 1){ #if (2+)/3 = A - looks like a polyA tail 	
				
				#in a polyA region - want to see how far this extends:
				my $pos = $i;
				
				while($pos > ($window_size - 1)){
					
					#move the window one position:
					if (ord($seq[$pos]) == 65){ #65 = decimal for 'A'
						$match--;
					}
					if (ord($seq[$pos - $window_size]) == 65){
						$match++;
					}
					
					if ($match < 2){ 
						#at end of the polyA region:
						
						#find length of polyA region: polya_len = (seq length - non-tail length - post-tail buffer length)
						my $polya_len = $length - ($pos - ($window_size - 1)) - ($length  - ($i + 1)); 
						
						#test to see if polyA string > post polyA region:
						if ($polya_len > ($length - $i)){ 
						
							#we now want to clip end of sequence
							#identify last non-A in this window:
							my $len;
							for (my $j = ($pos - $window_size); $j <= $pos; $j++){
								if (ord($seq[$j]) != 65){
									$len = ($j + 1);
								}else{
									last;
								}
							}
							
							$clipped_seq = substr($seq, 0, $len);
							
							#now, it might be that the sequence look something like ....AAAAACGAAAAA
							#in which case the newly clipped ....AAAAACG can be reexamined
							$clipped_seq = clip_polya($clipped_seq, $end);
							
						}							
						$pos = 0; #break out of while loop
					}
					$pos--;						
				}
				last; #move onto a new sequence
			}
		}
	}elsif($t_count > 4){

		#looking only for polyT head - use moving window looking for strings of 2/3 Ts:
	
		#moving through seq from front:
		for (my $i = 0; $i <= $buffer; $i++){
			my $match = 0;
			for (my $j = $i; $j < ($i + $window_size); $j++){ #check a window 
				if ($seq[$j] eq 'T'){ 
					$match++;
				}
			}
			
			
			if ($match > 1){ #if (2+)/3 = T - looks like a polyT head 

				my $pos = $i;

				#in a polyT region - want to see how far this extends:
				while($pos < ($length - $window_size)){
					
					#move the window one position
					if (ord($seq[$pos]) == 84){ #eq 'T'
						$match--;
					}
					if (ord($seq[$pos + $window_size]) == 84){
						$match++;
					}
					
					if ($match < 2){ 
						
						#at end of polyT region:
						
						#find length of polyT region: polyt_len = (head length - pre-head buffer length)
						my $polyt_len = ($pos + ($window_size - 1)) - ($i - 1);  
						
						#test to see if polyT string > pre polyT region:
						if ($polyt_len > $i){ 
						
							#we now want to clip front of sequence:
							#identify first non-T in this window:
							my $len;
							for (my $j = ($pos + ($window_size)); $j >= $pos; $j--){
								if (ord($seq[$j]) != 84){
									$len = $j;
								}else{
									last;
								}
							}

							$clipped_seq = substr($seq, $len);

							#now, it might be that the sequence look something like TTTTTCGTTTTT...
							#in which case the newly clipped CGTTTTT... can be reexamined
							$clipped_seq = clip_polya($clipped_seq, $end);
						}							
						$pos = $length; #break out of while loop
					}
					$pos++;						
				}
				last; #move onto a new sequence
			}
		}		
	}else{
		#broken!
	}
	return $clipped_seq;
}
