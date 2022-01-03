#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;                       # Erweiterte Verarbeitung von Befehlszeilenoptionen 
use autodie qw(:all);
use Data::Dumper;                   # stringified perl data structures, suitable for both printing and eval
use Algorithm::NeedlemanWunsch;
use Bio::SeqIO;
use Bio::Perl;
use Bio::Seq;
use Bio::Tools::CodonTable;
use Time::Duration;
use Parallel::ForkManager;
use List::MoreUtils qw(uniq);
use Sort::Key::Natural qw(natsort);
use List::Util qw(any);
use Array::Utils qw(:all);
use List::MoreUtils qw(first_index);
use POSIX;

#------------------------------------------------------------------------
#Definitions
#

my $extension = 35;
my $query;
my %model_informations;
my $locus;
my $distance_first = 10000; 	#Can be defined by the user in dependence of the intron size expected for the first (few) exons.   
my $distance_last = 10000;
my $ID;
my $core = 8;	
my $genomeIn;
my $querysequenceIn;
my $blastoutIn;
my $mode;

my $path_to_maxentscan_3 = "/scratch/franziska/bin/maxentscan/fordownload/score3.pl";
my $path_to_maxentscan_5 = "/scratch/franziska/bin/maxentscan/fordownload/score5.pl";

GetOptions("genome=s" => \$genomeIn, "querysequence=s" => \$querysequenceIn, "locus=s" => \$locus, "ID=s" => \$ID);
#------------------------------------------------------------------------
#Main

if (-d "MaxEntScan/"){
	`rm -r MaxEntScan`;
}

mkdir "MaxEntScan/";


my $strand_last ="";
my $strand;
my $start;
my $end;
my $first = `grep -P '$locus\t' SearchQuery/${ID}.blastout|sort -n -k 7,8| awk '\$11 <= 0.000001 {print \$0}'| awk '\$3 > 50 {print \$0}'| grep 'e-'|head -n 1`;
my @first = split "\t", $first;

if ($first[8]<$first[9]){
$strand = "+";
}
elsif ($first[8] > $first[9]){
$strand = "-";
};
my $last = `grep -P '$locus\t' SearchQuery/${ID}.blastout|sort -n -k 7,8| awk '\$11 <= 0.000001 {print \$0}'| awk '\$3 > 50 {print \$0}'| grep 'e-'|tail -n 1`;
my @last = split "\t", $last;
if ($last[8] < $last[9]){
$strand_last = "+";
}
elsif($last[8] > $last[9]){
$strand_last = "-";
}

if ($strand ne $strand_last){
print "WARNING: The paralog $ID has high scoring blast/hmm-hits on different strands in the query genome. Check this region with the help of the SearchQuery/${ID}.blastout. This might be a second paralog on the same contig or an assembly mistake. In the first case, split this contig in two parts so that each paralog is situated on a differently named contig and restart the pipeline. In the second case consider taking a different genome version/assembly.\n";
}
elsif ($strand eq $strand_last){

	if ($strand eq "+"){
  		$start = $first[8] - $distance_first;
  		$end = $last[9]+ $distance_last;
  		`fastacmd -d $genomeIn -s $locus -L $start,$end > SearchQuery/${ID}_locus.fa`;	
  		`formatdb -i SearchQuery/${ID}_locus.fa -p F -o`
	}
	elsif($strand eq "-"){
  		$end = $first[9] + $distance_first;
  		$start = $last[8] - $distance_last;
  		`fastacmd -d $genomeIn -S 2 -s $locus -L $start,$end > SearchQuery/${ID}_locus.fa`;
  		`formatdb -i SearchQuery/${ID}_locus.fa -p F -o`
	};

#nochmal blast auf locus
}
my @call = "blat SearchQuery/${ID}_locus.fa SearchQuery/${ID}.fa SearchQuery/${ID}_locus.blast8 -q=prot -t=dnax -out=blast8";
system (@call) == 0 or die "blat against query locus failed. $!\n";
print "Blat of query sequences against query locus done:\n";



open (my $blatout_handle, "<", "SearchQuery/${ID}_locus.blast8") or die "Can't open Blatout, $!";
my @a;
my @b;
my $a_align;
my $b_align;

#data from blatout
while (defined (my $line = <$blatout_handle>)){
	chomp $line;	
	#print Dumper ($line);
	next if ($line =~ m/^\#/ || $line =~ m/^$/);  # skip comments and empty lines
	if ($line =~ m/${locus}/){
		my @line_array = split /\t/, $line;
		if ($line_array[2] > 90.00){
			if ($line_array[10] < 0.1 || $line_array[2] == 100.00){

				my %data = search_data($line);

				my $genome_seqIO = Bio::SeqIO->new(-file   => "SearchQuery/".${ID}."_locus.fa", -format => "Fasta");
				$a = $genome_seqIO->next_seq;
				$a = $a -> seq;

				my @a_short;

				if ($data{strand} eq "+"){				
					if (($data{start}) - $extension < 0){
						$data{extension_start} = ($data{start} - $extension);
					}
					else{
						$data{extension_start} = $extension +1;
					}
					
					my $extended_start = $data{start} - $extension -1;

					if (((length $a) - $data{end}) < $extension){
						$data{extension_end} = ((length $a) - $data{end});
					}
					else{
						$data{extension_end} = $extension -1;
					}
					my $extended_end = $data{end} + $data{extension_end};	

					@a = split //, $a;
					for(my $i = $extended_start;$i <= $extended_end;$i++) {
						push (@a_short, $a[$i]);
					}
				}
				elsif($data{strand} eq "-"){				
					if (($data{end}) - $extension < 0){
						$data{extension_end} = ($data{end} - $extension);
					}
					else{
						$data{extension_end} = $extension +1;
					}
					
					my $extended_start = $data{end} -$extension -1;

					if (((length $a) - $data{start}) < $extension){
						$data{extension_start} = ((length $a) - $data{start});
					}
					else{
						$data{extension_start} = $extension -1;
					}
					my $extended_end = $data{start} + $data{extension_start};	
					@a = split //, $a;
					for(my $i = $extended_start;$i <= $extended_end;$i++) {
						push (@a_short, $a[$i]);
					}
					foreach (@a_short){
						if ($_ eq 'A'){$_ =~ s/A/T/;}
						elsif ($_ eq 'T'){$_ =~ s/T/A/;}
						elsif ($_ eq 'G'){$_ =~ s/G/C/;}
						elsif ($_ eq 'C'){$_ =~ s/C/G/;}
					}
					@a_short = reverse(@a_short);
				}
		
				$b = $querysequenceIn;
				my @b = split //, $b;
				$data{genome_part} = join ("",@a_short);
				$data{query} = $b;

				my $key_for_hash = $data{paralogID}."_".$data{start_query};
				$model_informations{$key_for_hash} = \%data;
			}
		}
	}					
}
close $blatout_handle;

if( !%model_informations){
	exit;
}

# for my $i (keys %model_informations) {
# 	if ($model_informations{$i}->{'bias'} < 90){
# 		delete $model_informations{$i};
# 	}
# }
my @exons;

my $pm = new Parallel::ForkManager($core);

for my $i (keys %model_informations) {
	push (@exons, $model_informations{$i}->{'start_query'});
}	
my @sorted_exons = sort {$a <=> $b} @exons;
my $lastexon=$sorted_exons[-1];

#%model_informations = find_best_exon_mapping(\@unique_exons, \%model_informations);

for my $i (keys %model_informations) {
	#print Dumper ($model_informations{$i});
	#$pm->start and next;
	my @a_short = split //, $model_informations{$i}->{'genome_part'};

	#search stopcodon next to last exon
	if ($model_informations{$i}->{'start_query'} == $lastexon){
		my ($stopcodon, $stop_nt) = find_stopcodon(@a_short);	
		$model_informations{$i}->{'stop_nt'}  = $stop_nt;		
		$model_informations{$i}->{'modulo_stop'} = $model_informations{$i}->{'stop_nt'} % 3;

		if ($stopcodon == 1 && ($model_informations{$i}->{'modulo_stop'} == 0)){
			$model_informations{$i}->{'score_5splicesite'} = 0;
			$model_informations{$i}->{'stopcodon'} = 1;
			$model_informations{$i}->{'startcodon'} = 0;
		}

		else {		
			until ($model_informations{$i}->{'modulo_stop'} == 0 || $model_informations{$i}->{'stop_nt'} >= $extension -4){
				#print Dumper ($model_informations{$i}->{'stop_nt'});
				($stopcodon, $stop_nt) = find_stopcodon_after_modulo(\@a_short,$model_informations{$i}->{'stop_nt'});	
				$model_informations{$i}->{'stop_nt'}  = $stop_nt;		
				$model_informations{$i}->{'modulo_stop'} = $model_informations{$i}->{'stop_nt'} % 3;
			}		
			if ($stopcodon == 1 && ($model_informations{$i}->{'modulo_stop'} == 0)){
				$model_informations{$i}->{'score_5splicesite'} = 0;
				$model_informations{$i}->{'stopcodon'} = 1;
				$model_informations{$i}->{'startcodon'} = 0;
			}
			else{
				$model_informations{$i}->{'stopcodon'} = 0;			
				$model_informations{$i}->{'startcodon'} = 0;
				print "no stopcodon: $model_informations{$i}\n";
				print Dumper ($model_informations{$i})
			}
		}

		$model_informations{$i} = best_three_splicesite($model_informations{$i}, 100);
	}

	#search startcodon in exon 0
	elsif ($model_informations{$i}->{'start_query'} eq "1"){
		my ($startcodon, $start_nt) = find_startcodon(\@a_short, $model_informations{$i}->{'start_query'});
		$model_informations{$i}->{'start_nt'}  = $start_nt;		
		$model_informations{$i}->{'modulo_start'} = $model_informations{$i}->{'start_nt'} % 3;
		if ($startcodon == 1 and $model_informations{$i}->{'modulo_start'} == 0){
			$model_informations{$i}->{'score_3splicesite'} = 0;
			$model_informations{$i}->{'startcodon'} = 1;
			$model_informations{$i}->{'stopcodon'} = 0;
		}
		else {
			until ($model_informations{$i}->{'modulo_start'} == 0 || $model_informations{$i}->{'start_nt'} >= $extension -4){
				my ($startcodon, $start_nt) = find_startcodon_after_modulo(\@a_short, $model_informations{$i}->{'start_nt'});
				$model_informations{$i}->{'start_nt'}  = $start_nt;		
				$model_informations{$i}->{'modulo_start'} = $model_informations{$i}->{'start_nt'} % 3;
			}
			if ($startcodon == 1 and $model_informations{$i}->{'modulo_start'} == 0){
				$model_informations{$i}->{'score_3splicesite'} = 0;
				$model_informations{$i}->{'startcodon'} = 1;
				$model_informations{$i}->{'stopcodon'} = 0;
			}			
			else{
				$model_informations{$i}->{'startcodon'} = 0;
				$model_informations{$i}->{'score_3splicesite'} = -333;
				$model_informations{$i}->{'stopcodon'} = 0;
				$model_informations{$i}->{'three_nt'}  = 0;
				print "no startcodon: $model_informations{$i}->{'start_query'}\n";
			}
		}

		$model_informations{$i} = best_five_splicesite($model_informations{$i}, 100);
	}

	# schauen, ob GT....AG (canonical) im Intron ist, GC-AG (noncanonical)
	# mammalian non-canonical splice sites, then the 99.24% of splice site pairs should be GT-AG, 0.69% GC-AG, 0.05% AT-AC and finally only 0.02% could consist of other types of non-canonical splice sites
	# (Analysis of canonical and non-canonical splice sites in mammalian genomes; M. Burset; 2000)

	else{	
		#print Dumper ("möp",$model_informations{$i});

		$model_informations{$i} = best_five_splicesite($model_informations{$i}, 100);

		$model_informations{$i} = best_three_splicesite($model_informations{$i}, 100);		
		#print Dumper ("mäp",$model_informations{$i});
	}
		#$pm->finish
}
#$pm->wait_all_children;
#maxentscan output positiv
# -> wie gross ist $nt?
#	-> <3, dann schauen, ob gesplittete as
for my $i (keys %model_informations) {
	#print Dumper ($i);
	if($model_informations{$i}->{'score_3splicesite'} == 0){
		$model_informations{$i}->{'three_nt_extension'} = 0;
		$model_informations{$i}->{'splitted_AS3'} = 0;
	}
	else{
	#if ($model_informations{$i}->{'score_3splicesite'} > 0){
	#		if ((int $data{three_nt} / 3) > 0); ganzzahlige division
		if ($model_informations{$i}->{'three_nt'} > 0 || $model_informations{$i}->{'three_nt'} < 0){
			#print Dumper ($i, $model_informations{$i}->{'three_nt'});
			if ($model_informations{$i}->{'three_nt'} > 0){
				$model_informations{$i}->{'splitted_AS3'} = $model_informations{$i}->{'three_nt'} % 3;
			}
			else {
				#print Dumper ($i, $model_informations{$i}->{'three_nt'});
				$model_informations{$i}->{'splitted_AS3'} = ($model_informations{$i}->{'three_nt'} +3) % 3;
				#print Dumper ($model_informations{$i}->{'splitted_AS3'});
			}

			if ($model_informations{$i}->{'splitted_AS3'} > 0 || $model_informations{$i}->{'splitted_AS3'} < 0){
			# in fasta header 
				$model_informations{$i}->{'splitted_for_header3'} = "3: ".$model_informations{$i}->{'splitted_AS3'};
				# print Dumper ("splitted_for_header: ".$model_informations{$i}->{'splitted_for_header3'});
			}
		#	-> durch 3 teilbar, query verlaengern
			if ($model_informations{$i}->{'strand'} eq "+"){
				$model_informations{$i}->{'three_nt_extension'} = POSIX::floor($model_informations{$i}->{'three_nt'} / 3) * 3;
				my $extension = (int $model_informations{$i}->{'three_nt'} / 3) * 3;
				$model_informations{$i}->{'start'} = $model_informations{$i}->{'start'} - $extension;
			}
			elsif($model_informations{$i}->{'strand'} eq "-"){
				$model_informations{$i}->{'three_nt_extension'} = POSIX::floor($model_informations{$i}->{'three_nt'} / 3) * 3;
				my $extension = (int $model_informations{$i}->{'three_nt'} / 3) * 3;
				$model_informations{$i}->{'start'} = $model_informations{$i}->{'start'} + $extension;
			}
		}
		else{
			$model_informations{$i}->{'splitted_AS3'} = 0;
			$model_informations{$i}->{'three_nt_extension'} = 0;			
		}
	}
	if ($model_informations{$i}->{'score_5splicesite'} == 0){
		$model_informations{$i}->{'five_nt_extension'} = 0;
		$model_informations{$i}->{'splitted_AS5'} = 0;
	}
	else{
		if ($model_informations{$i}->{'five_nt'} > 0 || $model_informations{$i}->{'five_nt'} < 0){			
			if ($model_informations{$i}->{'five_nt'} > 0){
				$model_informations{$i}->{'splitted_AS5'} = $model_informations{$i}->{'five_nt'} % 3;
			}
			else {
				$model_informations{$i}->{'splitted_AS5'} = ($model_informations{$i}->{'five_nt'} +3) % 3;
			}

			if ($model_informations{$i}->{'splitted_AS5'} > 0 || $model_informations{$i}->{'splitted_AS5'} < 0){
			# in fasta header 
				$model_informations{$i}->{'splitted_for_header5'} = "5: ".$model_informations{$i}->{'splitted_AS5'};
				# print Dumper ("splitted_for_header: ".$model_informations{$i}->{'splitted_for_header5'});
					#	-> durch 3 teilbar, query verlaengern
			}
			if ($model_informations{$i}->{'strand'} eq "+"){
				$model_informations{$i}->{'five_nt_extension'} = POSIX::floor($model_informations{$i}->{'five_nt'} / 3) * 3;
				my $extension = (int $model_informations{$i}->{'five_nt'} / 3) * 3;
				$model_informations{$i}->{'end'} = $model_informations{$i}->{'end'} + $extension;
			}
			elsif($model_informations{$i}->{'strand'} eq "-"){
				$model_informations{$i}->{'five_nt_extension'} = POSIX::floor($model_informations{$i}->{'five_nt'} / 3) * 3;
				my $extension = (int $model_informations{$i}->{'five_nt'} / 3) * 3;
				$model_informations{$i}->{'end'} = $model_informations{$i}->{'end'} - $extension;
			}
		}
		else{
			$model_informations{$i}->{'splitted_AS5'} = 0;
			$model_informations{$i}->{'five_nt_extension'} = 0;
		}
	}
}


my @keys2 = (natsort keys %model_informations);
#print Dumper (@keys2);
for my $i (@keys2) {
	#print Dumper (%model_informations{$i});		
	my $idx = first_index {$_ =~ m/$i/} @keys2;
	my $model_compare = "";
	if (defined $keys2[$idx - 1] && ($idx -1) >= 0){
		$model_compare = $keys2[$idx - 1];
		#print Dumper($model_compare);
	}
	else{
		next;
	}

	#print Dumper ("i= ".$i);
	if ($model_informations{$model_compare}->{'splitted_AS5'} + $model_informations{$i}->{'splitted_AS3'} == 3 || $model_informations{$model_compare}->{'splitted_AS5'} + $model_informations{$i}->{'splitted_AS3'} == 0){
			#print Dumper ("splitted AS3");
	}
	else{
		my %possible_matches;
		my $scores_ref_three;
		my $ref_five_scores = $model_informations{$model_compare}->{'five_scores'};
		#print Dumper ("key", %$ref_five_scores);
		for my $entry (keys %$ref_five_scores){
			#print Dumper ("entry",$entry);
			for (my $n = 0;$n <= $model_informations{$i}->{'number_three_scores'}-1;$n++){

				$scores_ref_three = $model_informations{$i}->{'three_scores'}; 
				my $next_best_three_score = (sort { $a <=> $b } (keys %$scores_ref_three))[-$n];
				#print Dumper ($next_best_three_score);

				$model_informations{$i}->{'score_3splicesite'} = $next_best_three_score;
				$model_informations{$i}->{'three_nt'} = $model_informations{$i}->{'three_scores'}->{$next_best_three_score};

				if ($model_informations{$i}->{'three_nt'} > 0){
					$model_informations{$i}->{'splitted_AS3'} = $model_informations{$i}->{'three_nt'} % 3;
					$model_informations{$i}->{'splitted_for_header3'} = "3: ".$model_informations{$i}->{'splitted_AS3'};
				}
				else{
					$model_informations{$i}->{'splitted_AS3'} = ($model_informations{$i}->{'three_nt'} +3) % 3;
					$model_informations{$i}->{'splitted_for_header3'} = "3: ".$model_informations{$i}->{'splitted_AS3'};
				}

				$model_informations{$model_compare}->{'score_5splicesite'} = $entry;
				$model_informations{$model_compare}->{'five_nt'} = $$ref_five_scores{$entry};
				#print Dumper ($model_informations{$model_compare}->{'score_5splicesite'}, $model_informations{$model_compare}->{'five_nt'} );		
				if ($model_informations{$model_compare}->{'five_nt'} > 0){				
					$model_informations{$model_compare}->{'splitted_AS5'} = $model_informations{$model_compare}->{'five_nt'} % 3;
					$model_informations{$model_compare}->{'splitted_for_header5'} = "5: ".$model_informations{$model_compare}->{'splitted_AS5'};
				}
				else{
					$model_informations{$model_compare}->{'splitted_AS5'} = ($model_informations{$model_compare}->{'five_nt'}+3) % 3;
					$model_informations{$model_compare}->{'splitted_for_header5'} = "5: ".$model_informations{$model_compare}->{'splitted_AS5'};

				}
				#print Dumper ($model_informations{$i}->{'splitted_AS3'} ,$model_informations{$model_compare}->{'splitted_AS5'} );
				if ($model_informations{$model_compare}->{'splitted_AS5'} + $model_informations{$i}->{'splitted_AS3'} == 3 || $model_informations{$model_compare}->{'splitted_AS5'} + $model_informations{$i}->{'splitted_AS3'} == 0){
					my %temp_match;
					$temp_match{$model_informations{$i}->{'score_3splicesite'}} = $model_informations{$i}->{'three_nt'};
					$temp_match{$model_informations{$model_compare}->{'score_5splicesite'}} = $model_informations{$model_compare}->{'five_nt'};
					my $key = ($model_informations{$i}->{'score_3splicesite'} + $model_informations{$model_compare}->{'score_5splicesite'});
					$possible_matches{$key} = \%temp_match;
				}
			}
		}
		#print Dumper (%possible_matches);
		if( !%possible_matches){
			next;
		}

		my $best_match = (sort { $a <=> $b } (keys %possible_matches))[-1];

		my $common = "";
		my $common_three = "";

		foreach (keys %$ref_five_scores) {
			if (exists $possible_matches{$best_match}->{$_}){
				$common = $_;
				}
		}
		foreach (keys %$scores_ref_three) {
			if (exists $possible_matches{$best_match}->{$_}){
				$common_three = $_;
				}
		}
		#print Dumper ($common_three);		

		$model_informations{$model_compare}->{'score_5splicesite'} = $common;
		$model_informations{$model_compare}->{'five_nt'} = $$ref_five_scores{$common};	
		if ($model_informations{$model_compare}->{'five_nt'} > 0){				
			$model_informations{$model_compare}->{'splitted_AS5'} = $model_informations{$model_compare}->{'five_nt'} % 3;
			$model_informations{$model_compare}->{'splitted_for_header5'} = "5: ".$model_informations{$model_compare}->{'splitted_AS5'};
		}
		else{
			$model_informations{$model_compare}->{'splitted_AS5'} = ($model_informations{$model_compare}->{'five_nt'}+3) % 3;
			$model_informations{$model_compare}->{'splitted_for_header5'} = "5: ".$model_informations{$model_compare}->{'splitted_AS5'};

		}
		#print Dumper ($model_informations{$model_compare}->{'score_5splicesite'}, $model_informations{$model_compare}->{'five_nt'} );



		$model_informations{$i}->{'score_3splicesite'} = $common_three;
		$model_informations{$i}->{'three_nt'} = $model_informations{$i}->{'three_scores'}->{$common_three};				
		if ($model_informations{$i}->{'three_nt'} > 0){
			$model_informations{$i}->{'splitted_AS3'} = $model_informations{$i}->{'three_nt'} % 3;
			$model_informations{$i}->{'splitted_for_header3'} = "3: ".$model_informations{$i}->{'splitted_AS3'};
		}
		else{
			$model_informations{$i}->{'splitted_AS3'} = ($model_informations{$i}->{'three_nt'} +3) % 3;
			$model_informations{$i}->{'splitted_for_header3'} = "3: ".$model_informations{$i}->{'splitted_AS3'};
		}
	}
}

# fasta output for every paralog 
my $dna;
my @targetsequence;

#open (my $OUT, ">> SearchQuery/$model_informations{$entry}->{'paralogID'}_exon$idx");
#open (my $OUT, ">> HomologousExons/paralog_exon_alignment.fa");
my $myCodonTable   = Bio::Tools::CodonTable->new();

for my $exon (@sorted_exons){
	#print Dumper ($exon); 
	my $entry = $ID."_".$exon;
	my ($new_three_nt, $new_five_nt);
	#print Dumper($model_informations{$entry}->{'score_5splicesite'},$model_informations{$entry}->{'score_3splicesite'});
	if ($model_informations{$entry}->{'score_5splicesite'} >= -6 && $model_informations{$entry}->{'score_3splicesite'} >= -6){
		if ($model_informations{$entry}->{'start_query'} == 1){
			$new_three_nt = $model_informations{$entry}->{'start_nt'} + 3;
			$new_five_nt = $model_informations{$entry}->{'five_nt_extension'};
		}

		elsif($model_informations{$entry}->{'start_query'} == $lastexon){
			$new_three_nt = $model_informations{$entry}->{'three_nt_extension'};
			$new_five_nt = $model_informations{$entry}->{'stop_nt'};
		}

		else{
			$new_three_nt = $model_informations{$entry}->{'three_nt_extension'};
			$new_five_nt = $model_informations{$entry}->{'five_nt_extension'};
		}
		#print Dumper ($new_three_nt, $new_five_nt);
		my @exonsequence = sequence_final_exon($model_informations{$entry}->{'genome_part'}, $new_three_nt, $new_five_nt);
		my $dna_exonsequence= join ("",@exonsequence);
		#print Dumper ($dna_exonsequence);
		my @triplets = unpack '(a3)*', $dna_exonsequence;
		my @protein_aa_exon;
		for my $codon (@triplets){
			my $aa_exon = $myCodonTable->translate($codon);
	    	push (@protein_aa_exon, $aa_exon);
		}
		my $protein_exon = join ("",@protein_aa_exon);

		my $idx = first_index {$_ =~ m/$model_informations{$entry}->{'start_query'}/} @sorted_exons;
		open (my $OUT, ">> SearchQuery/$ID.inf");
		print $OUT ">$model_informations{$entry}->{'paralogID'}_exon$idx\n$protein_exon\n";
	}
}

for my $exon (@sorted_exons){
#print Dumper ($exon); 
	my $entry = $ID."_".$exon;
	my ($new_three_nt, $new_five_nt);
	if ($model_informations{$entry}->{'score_5splicesite'} >= -5 && $model_informations{$entry}->{'score_3splicesite'} >= -5){
		if ($model_informations{$entry}->{'start_query'} == 1){
			$new_three_nt = $model_informations{$entry}->{'start_nt'} + 3;
			$new_five_nt = $model_informations{$entry}->{'five_nt'};
		}

		elsif($model_informations{$entry}->{'start_query'} == $lastexon){
			$new_three_nt = $model_informations{$entry}->{'three_nt'};
			$new_five_nt = $model_informations{$entry}->{'stop_nt'};
		}

		else{
			$new_three_nt = $model_informations{$entry}->{'three_nt'};
			$new_five_nt = $model_informations{$entry}->{'five_nt'};
		}

		my @exonsequence = sequence_final_exon($model_informations{$entry}->{'genome_part'}, $new_three_nt, $new_five_nt);
		push (@targetsequence, @exonsequence);
	}
	else {
		print Dumper("ne", $exon);
	}
}

$dna= join ("",@targetsequence);
#print Dumper ($dna);
my @triplets = unpack '(a3)*', $dna;
my @protein_aa;
for my $codon (@triplets){
	my $aa = $myCodonTable->translate($codon);
    push (@protein_aa, $aa);
}
my $protein = join ("",@protein_aa);
open (my $OUTQUERY, ">> HomologousExons/paralogs.fa");	
print $OUTQUERY ">$ID\n$protein\n";

print Dumper ($protein);
#------------------------------------------------------------------------
# SUBROUTINES
#

#find start, end, paralog name from blastout

sub search_data{
	my ($blastoutIn_line) = @_;

	my @line_array = split "\t", $blastoutIn_line;

	my $paralogID = $line_array[0];
	my $locus = $line_array[1];
	my $bias = $line_array[2];
	my $mismatch = $line_array[4];
	my $gap = $line_array[5];
	my $start_query = $line_array[6];
	my $end_query = $line_array[7];
	my $start = $line_array[8];
	my $end = $line_array[9];
	my $evalue = $line_array[10];
	my $strand;

	if ($line_array[8] < $line_array[9]){
		$strand = "+";
	}
	elsif ($line_array[8] > $line_array[9]){
		$strand = "-";
	};

	my %line_record = (
		start_query => $start_query,
    	end_query => $end_query,
    	start => $start,
    	end => $end,
    	strand => $strand,
    	paralogID => $paralogID,
    	locus => $locus,
    	evalue => $evalue,
    	bias => $bias,
    	mismatch => $mismatch,
    	gap => $gap,
	);
}

sub best_three_splicesite{	
	my ($model_informations, $nt_last) = @_ ;
	#print Dumper ("nt_last",$nt_last);
	my @a_short = split //, $model_informations->{genome_part};

	my ($three_splicesite, $three_nt) = three_splice_site(\@a_short, $model_informations->{mismatch}, $model_informations->{start_query});

	$model_informations->{three_nt} = $three_nt;
	$model_informations->{three_splicesite} = $three_splicesite;
	
	my %three_scores;			

	if ($three_splicesite == 1){
		my @testsequence_three_splicesite = sequence_three_splicesite(\@a_short, $three_nt);	

		if (scalar @testsequence_three_splicesite == 23){
   			$model_informations->{testsequence_three_splicesite} = join ("",@testsequence_three_splicesite);
			open (testsequences_three_splicesite_file, '>>', "MaxEntScan/Testsequences_three_splicesite.txt") or die "$!";
			print testsequences_three_splicesite_file join ("",@testsequence_three_splicesite);
		}
		else {
			print "unvollstaendig 3: $model_informations->{model}\n";
		}

		if(-e "MaxEntScan/Testsequences_three_splicesite.txt"){
			my $output_3splicesite = qx(perl $path_to_maxentscan_3 "MaxEntScan/Testsequences_three_splicesite.txt");
			#print Dumper ($output_3splicesite);
			my @array_3splicesite = split /\t/, $output_3splicesite;
			if(defined $array_3splicesite[1]){
				chomp ($array_3splicesite[1]);
			}
			$model_informations->{score_3splicesite} = $array_3splicesite[1];	
			my $key_for_hash = $model_informations->{score_3splicesite};
			$three_scores{$key_for_hash} = $model_informations->{three_nt};
			unlink "MaxEntScan/Testsequences_three_splicesite.txt";
		}
		$three_splicesite = 0;	
	}

	if($three_splicesite == 0 || $model_informations->{score_3splicesite} < 0){
		$model_informations->{three_nt} = 0;
		my ($three_splicesite, $three_nt) = three_splice_site_minus_three(@a_short);

		$model_informations->{three_nt} = $three_nt;
		$model_informations->{three_splicesite} = $three_splicesite;

		if ($three_splicesite == 1){
			my @testsequence_three_splicesite = sequence_three_splicesite(\@a_short, $three_nt);	

			if (scalar @testsequence_three_splicesite == 23){
	   			$model_informations->{testsequence_three_splicesite} = join ("",@testsequence_three_splicesite);
				open (testsequences_three_splicesite_file, '>>', "MaxEntScan/Testsequences_three_splicesite.txt") or die "$!";
				print testsequences_three_splicesite_file join ("",@testsequence_three_splicesite);
			}
			else {
				print "unvollstaendig 3: $model_informations->{model}\n";
			}

			if(-e "MaxEntScan/Testsequences_three_splicesite.txt"){
				my $output_3splicesite = qx(perl $path_to_maxentscan_3 "MaxEntScan/Testsequences_three_splicesite.txt");
				#print Dumper ($output_3splicesite);
				my @array_3splicesite = split /\t/, $output_3splicesite;
				if(defined $array_3splicesite[1]){
					chomp ($array_3splicesite[1]);
				}
				$model_informations->{score_3splicesite} = $array_3splicesite[1];
				my $key_for_hash = $model_informations->{score_3splicesite};
				$three_scores{$key_for_hash} = $model_informations->{three_nt};
				unlink "MaxEntScan/Testsequences_three_splicesite.txt";
			}
		}
		$three_splicesite = 0;
	}

	if($three_splicesite == 0 || $model_informations->{score_3splicesite} < 0){
		$model_informations->{score_3splicesite} = -333;
		until ($model_informations->{score_3splicesite} > 0 || $model_informations->{three_nt} >= $extension -22){
			my @a_short = split //, $model_informations->{genome_part};
			my ($three_splicesite, $three_nt) = three_splice_site_after_negativ_score(\@a_short, $model_informations->{three_nt});
			$model_informations->{three_nt} = $three_nt;
			$model_informations->{three_splicesite} = $three_splicesite;
			if ($three_splicesite == 1){
				my @testsequence_three_splicesite = sequence_three_splicesite(\@a_short, $model_informations->{three_nt});
				if (scalar @testsequence_three_splicesite == 23){
					$model_informations->{testsequence_three_splicesite} = join ("",@testsequence_three_splicesite);
					open (testsequences_three_splicesite_file, '>>', "MaxEntScan/Testsequences_three_splicesite.txt") or die "$!";
					print testsequences_three_splicesite_file join ("",@testsequence_three_splicesite);
					my $output_3splicesite = qx(perl $path_to_maxentscan_3 "MaxEntScan/Testsequences_three_splicesite.txt");
					#print Dumper ($output_3splicesite);
					my @array_3splicesite = split /\t/, $output_3splicesite;
					chomp ($array_3splicesite[1]);
					$model_informations->{score_3splicesite} = $array_3splicesite[1];
					my $key_for_hash = $model_informations->{score_3splicesite};
					$three_scores{$key_for_hash} = $model_informations->{three_nt};
					if(-e "MaxEntScan/Testsequences_three_splicesite.txt"){
						unlink "MaxEntScan/Testsequences_three_splicesite.txt";
					}
				}
			}				
			else{
			 	$model_informations->{score_3splicesite} = -333;
			}
		}
	}
	#print Dumper(%three_scores);
	if( !%three_scores){
		return ($model_informations);
	}
	#print Dumper ($model_informations);
	#print Dumper(%three_scores);
	my $best_three_score = (sort { $a <=> $b } (keys %three_scores))[-1];
	#print Dumper ($best_three_score);
	$model_informations->{number_three_scores} = %three_scores;
	$model_informations->{three_scores} = \%three_scores;
	$model_informations->{score_3splicesite} = $best_three_score;
	$model_informations->{three_nt} = $three_scores{$best_three_score};

	# my $best_three_score = (sort { $a <=> $b } (keys %three_scores))[-1];
	# if ($three_scores{$best_three_score} == $nt_last){
	# 	$best_three_score = (sort { $a <=> $b } (keys %three_scores))[-2];	
	# }
	# #print Dumper ($best_three_score);
	# $model_informations->{score_3splicesite} = $best_three_score;
	# $model_informations->{three_nt} = $three_scores{$best_three_score};  
	# if ($model_informations->{score_3splicesite} < -3){
	# 		$model_informations->{score_3splicesite} = 0;
	# } 
	return ($model_informations);
}

sub best_five_splicesite{    
	my ($model_informations, $nt_last) = @_ ;
	#print Dumper ($model_informations, $nt_last);
	my @a_short = split //, $model_informations->{genome_part};
	# print Dumper (@a_short);
	$model_informations->{stopcodon} = 0;
	$model_informations->{startcodon} = 0;

	my ($five_splicesite, $five_nt) = five_splice_site_canonical(\@a_short, $model_informations->{mismatch}, $model_informations->{gap});
	#print Dumper($five_splicesite);
	$model_informations->{five_nt} = $five_nt;
	$model_informations->{five_splicesite} = $five_splicesite;
	#print Dumper($model_informations->{five_nt} );
	my %five_scores;

	if ($five_splicesite == 1){
		#print Dumper ("1");
		my @testsequence_five_splicesite = sequence_five_splicesite(\@a_short, $five_nt); 
		if (defined $testsequence_five_splicesite[8]){
			$model_informations->{testsequence_five_splicesite} = join ("",@testsequence_five_splicesite);
			open (testsequences_five_splicesite_file, '>>', "MaxEntScan/Testsequences_five_splicesite.txt") or die "$!";
			print testsequences_five_splicesite_file join ("",@testsequence_five_splicesite);
		}
		else {
			print "unvollstaendig 5: %model_informations{model}\n";
		}
		if(-e "MaxEntScan/Testsequences_five_splicesite.txt"){
			my $output_5splicesite = qx(perl $path_to_maxentscan_5 "MaxEntScan/Testsequences_five_splicesite.txt");
			my @array_5splicesite = split /\t/, $output_5splicesite;
			#print Dumper($output_5splicesite);
			if (defined $array_5splicesite[1]){
				chomp ($array_5splicesite[1]);
			}
			$model_informations->{score_5splicesite} = $array_5splicesite[1];	
			my $key_for_hash = $model_informations->{score_5splicesite};
			$five_scores{$key_for_hash} = $model_informations->{five_nt};
			unlink "MaxEntScan/Testsequences_five_splicesite.txt";
		}		
	}
	$five_splicesite = 0;
	if ($five_splicesite == 0 || $model_informations->{score_5splicesite} < 0 ){
		$model_informations->{score_5splicesite} = -555;
		until ($model_informations->{score_5splicesite} > 0 || $model_informations->{five_nt} > $extension -4){			
			#print Dumper ("3");
			my ($five_splicesite, $five_nt) = five_splice_site_canonical_after_negative_score(\@a_short, $model_informations->{five_nt});
			$model_informations->{five_nt} = $five_nt;
			$model_informations->{five_splicesite} = $five_splicesite;
			if ($five_splicesite == 1){
				my @testsequence_five_splicesite = sequence_five_splicesite(\@a_short, $model_informations->{five_nt});
				if (scalar @testsequence_five_splicesite == 9){
					$model_informations->{testsequence_five_splicesite} = join ("",@testsequence_five_splicesite);
					open (testsequences_five_splicesite_file, '>>', "MaxEntScan/Testsequences_five_splicesite.txt") or die "$!";
					print testsequences_five_splicesite_file join ("",@testsequence_five_splicesite);
					my $output_5splicesite = qx(perl $path_to_maxentscan_5 "MaxEntScan/Testsequences_five_splicesite.txt");
					my @array_5splicesite = split /\t/, $output_5splicesite;
					#print Dumper($output_5splicesite);
					chomp ($array_5splicesite[1]);
					$model_informations->{score_5splicesite} = $array_5splicesite[1];			
					my $key_for_hash = $model_informations->{score_5splicesite};
					$five_scores{$key_for_hash} = $model_informations->{five_nt};
					if(-e "MaxEntScan/Testsequences_five_splicesite.txt"){
						unlink "MaxEntScan/Testsequences_five_splicesite.txt";
					}
				}
			}
		}
	}

	if ($five_splicesite == 0 || $model_informations->{score_5splicesite} < 0 ){
		#print Dumper ("4");
		$model_informations->{five_nt} = -3;
		until ($model_informations->{score_5splicesite} > 0 || $model_informations->{five_nt} > $extension -4){
			my ($five_splicesite, $five_nt) = five_splice_site_noncanonical(\@a_short,$model_informations->{five_nt});
			$model_informations->{five_nt} = $five_nt;
			$model_informations->{five_splicesite} = $five_splicesite;
			if ($five_splicesite == 1){
				my @testsequence_five_splicesite = sequence_five_splicesite(\@a_short, $model_informations->{five_nt});
				if (scalar @testsequence_five_splicesite == 9){
					$model_informations->{testsequence_five_splicesite} = join ("",@testsequence_five_splicesite);
					open (testsequences_five_splicesite_file, '>>', "MaxEntScan/Testsequences_five_splicesite.txt") or die "$!";
					print testsequences_five_splicesite_file join ("",@testsequence_five_splicesite);
					my $output_5splicesite = qx(perl $path_to_maxentscan_5 "MaxEntScan/Testsequences_five_splicesite.txt");
					my @array_5splicesite = split /\t/, $output_5splicesite;
					chomp ($array_5splicesite[1]);
					$model_informations->{score_5splicesite} = $array_5splicesite[1];					
					my $key_for_hash = $model_informations->{score_5splicesite};
					$five_scores{$key_for_hash} = $model_informations->{five_nt};
					if(-e "MaxEntScan/Testsequences_five_splicesite.txt"){
						unlink "MaxEntScan/Testsequences_five_splicesite.txt";
					}
				}
			}
			else{
			 	$model_informations->{score_5splicesite} = -555;
			}
		}
	}			

	if( !%five_scores){
		return ($model_informations);
	}
	#print Dumper ($model_informations);
	#print Dumper(%five_scores);
	my $best_five_score = (sort { $a <=> $b } (keys %five_scores))[-1];
	#print Dumper ($best_five_score);
	$model_informations->{number_five_scores} = %five_scores;
	$model_informations->{five_scores} = \%five_scores;
	$model_informations->{score_5splicesite} = $best_five_score;
	$model_informations->{five_nt} = $five_scores{$best_five_score};



	# print Dumper(%five_scores);
	# my $best_five_score = (sort { $a <=> $b } (keys %five_scores))[-1];	
	# if ($five_scores{$best_five_score} == $nt_last){
	# 	$best_five_score = (sort { $a <=> $b } (keys %five_scores))[-2];
	# }
	# print Dumper ($best_five_score);
	# $model_informations->{score_5splicesite} = $best_five_score;
	# $model_informations->{five_nt} = $five_scores{$best_five_score}; 		
	# if ($model_informations->{score_5splicesite} < -3){
	# 		$model_informations->{score_5splicesite} = 0;
	# 	}
	# #print Dumper ("ups",$model_informations);
    return ($model_informations);
}

sub find_stopcodon {
    my (@genome_target) = @_;  

    my $stopcodon; #1 if stopcodon is found, else 0
    my $stop_nt = 0; #nt between stopcodon and query exon

	for(my $i = ((scalar @genome_target) - $extension);$i <= (scalar @genome_target) - 3;$i++){	
        if (($genome_target[$i] eq "T" and $genome_target[$i+1] eq "A" and $genome_target[$i+2] eq "G") or ($genome_target[$i] eq "T" and $genome_target[$i+1] eq "G" and $genome_target[$i+2] eq "A") or ($genome_target[$i] eq "T" and $genome_target[$i+1] eq "A" and $genome_target[$i+2] eq "A")) {     

        	$stopcodon = 1;
        	return ($stopcodon, $stop_nt);
        } 
        else {
        	$stopcodon = 0;
        	$stop_nt++;
        }
    }
    return ($stopcodon, $stop_nt);
}

sub find_stopcodon_after_modulo {
    my ($genome_target_ref, $stop_nt) = @_; 
    my @genome_target = @$genome_target_ref; 
	my $stopcodon = 0;

	for(my $i = ((scalar @genome_target) - $extension + $stop_nt + 1);$i <= (scalar @genome_target) - 3;$i++){	
        if (($genome_target[$i] eq "T" and $genome_target[$i+1] eq "A" and $genome_target[$i+2] eq "G") or ($genome_target[$i] eq "T" and $genome_target[$i+1] eq "G" and $genome_target[$i+2] eq "A") or ($genome_target[$i] eq "T" and $genome_target[$i+1] eq "A" and $genome_target[$i+2] eq "A")) {     

        	$stopcodon = 1;
        	$stop_nt++;
        	return ($stopcodon, $stop_nt);
        } 
        else {
        	$stopcodon = 0;
        	$stop_nt++;
        }
    }
    return ($stopcodon, $stop_nt);
}
sub find_startcodon {
    my ($genome_target_ref, $qstarts) = @_;   
	my @genome_target = @$genome_target_ref;
	my $qstart = (split ",", $qstarts)[0];

    my $startcodon; #1 if startcodon is found, else 0
    my $start_nt = -6 + ($qstart *3);; #nt between startcodon and query exon

	for(my $i = $extension + 2 - ($qstart*3) + 3;$i >= 1;$i--){
		#print Dumper ($genome_target[$i]);
        if ($genome_target[$i] eq "G" and $genome_target[$i-1] eq "T" and $genome_target[$i-2] eq "A"){        
        	$startcodon = 1;    
        	return ($startcodon, $start_nt);
        } 
        else {
        	$startcodon = 0;
        	$start_nt++;
        }
    }
    return ($startcodon, $start_nt);
}

sub find_startcodon_after_modulo {
    my ($genome_target_ref, $start_nt) = @_;   
	my @genome_target = @$genome_target_ref;
    my $startcodon = 0; #1 if startcodon is found, else 0


	for(my $i = $extension +2  - $start_nt -1;$i >= 1;$i--){
		#print Dumper ($genome_target[$i]);
        if ($genome_target[$i] eq "G" and $genome_target[$i-1] eq "T" and $genome_target[$i-2] eq "A"){        
        	$startcodon = 1;    
        	return ($startcodon, $start_nt);
        	$start_nt++;
        } 
        else {
        	$startcodon = 0;
        	$start_nt++;
        }
    }
    return ($startcodon, $start_nt);
}
#find 3' splice site in intron AG
sub three_splice_site_block{    
	my ($genome_target_ref, $start) = @_;   
	my @genome_target = @$genome_target_ref;
    my $three_splicesite; #1 if ss is found, else 0
 	my $three_nt = 0;


	for(my $i = (scalar @genome_target) - $extension - ($start * 3) ;$i >= 16;$i--){
		#print Dumper ($genome_target[$i]);
        if ($genome_target[$i] eq "G" and $genome_target[$i-1] eq "A"){        
        	$three_splicesite = 1;    
        	return ($three_splicesite, $three_nt);
        } 
        else {
        	$three_splicesite = 0;
        	$three_nt++;
        }
    }
    return ($three_splicesite, $three_nt);
}

sub three_splice_site {
    my ($genome_target_ref, $mismatch, $qstarts) = @_;   
	my @genome_target = @$genome_target_ref;
	my $qstart = (split ",", $qstarts)[0];
	my $three_nt;
	#print Dumper ($qstart);
    my $three_splicesite = 0; #1 if ss is found, else 0
    if ($mismatch < 5){
	$three_nt = 0 + ($qstart *3);
		for(my $i = $extension -1 - ($qstart*3) ;$i >= 20;$i--){			
			#print Dumper ($genome_target[$i]);
	        if ($genome_target[$i] eq "G" and $genome_target[$i-1] eq "A"){        
	        	$three_splicesite = 1;    
	        	return ($three_splicesite, $three_nt);
	        } 
	        else {
	        	$three_splicesite = 0;
	        	$three_nt++;
	        }
	    }
	}
	else{
		$three_nt = 0 + ($qstart *3) - ($mismatch *3);		
		for(my $i = $extension -1 - ($qstart*3) + ($mismatch *3); $i >= 16;$i--){			
			#print Dumper ($genome_target[$i]);
	        if ($genome_target[$i] eq "G" and $genome_target[$i-1] eq "A"){        
	        	$three_splicesite = 1;    
	        	return ($three_splicesite, $three_nt);
	        } 
	        else {
	        	$three_splicesite = 0;
	        	$three_nt++;
	        }
	    }
	}
    return ($three_splicesite, $three_nt);
}

sub three_splice_site_minus_three {
    my (@genome_target) = @_;   
    my $three_splicesite; #1 if ss is found, else 0
 	my $three_nt = - 9;
	for(my $i = $extension +8 ;$i >= 20;$i--){
		#print Dumper ($genome_target[$i]);
        if ($genome_target[$i] eq "G" and $genome_target[$i-1] eq "A"){        
        	$three_splicesite = 1;    
        	return ($three_splicesite, $three_nt);
        } 
        else {
        	$three_splicesite = 0;
        	$three_nt++;
        }
    }
    return ($three_splicesite, $three_nt);
}

sub three_splice_site_after_negativ_score {
    my ($genome_target_ref, $three_nt) = @_;   
    my @genome_target = @$genome_target_ref;
    my $three_splicesite = 0; #1 if ss is found, else 0

	for(my $i = ($extension -3 - $three_nt);$i >= 20;$i--){
        if ($genome_target[$i] eq "G" and $genome_target[$i-1] eq "A"){        
        	$three_splicesite = 1;    
        	return ($three_splicesite, ($three_nt + 2));
        } 
        else {
        	$three_splicesite = 0;
        	$three_nt++;
        }
    }
    return ($three_splicesite, $three_nt);
}

sub three_splice_site_intron {
    my (@genome_target) = @_;   
    my $three_splicesite = 0; #1 if ss is found, else 0
	my $three_nt = 0;
	#print Dumper (@genome_target);
	for(my $i = ((scalar @genome_target) -1);$i >= 2;$i--){
		#print Dumper ($genome_target[$i]);
        if ($genome_target[$i] eq "G" and $genome_target[$i-1] eq "A"){  
        	#print Dumper ($genome_target[$i]);      
        	$three_splicesite = 1;    
        	return ($three_splicesite, $three_nt);
        } 
        else {
        	$three_splicesite = 0;
        	$three_nt++;
        }
    }
    return ($three_splicesite, $three_nt);
}
sub five_splice_site_canonical_block{
    my ($genome_target_ref, $start) = @_;  
    #print Dumper (@$genome_target_ref);
    my @genome_target = @$genome_target_ref;
    my $five_splicesite; #1 if splicesite is found, else 0
    my $five_nt = 0; #nt between splice site and query exon

	for(my $i = ($extension + ($start * 3) -3 );$i <= (scalar @genome_target) - $extension -2;$i++){	
		#print Dumper ($genome_target[$i]);
        if ($genome_target[$i] eq "G" and $genome_target[$i+1] eq "T"){     

        	$five_splicesite = 1;    
        	return ($five_splicesite, $five_nt -3 );
        } 
        else {
        	$five_splicesite = 0;
        	$five_nt++;
        }
    }
    return ($five_splicesite, $five_nt);
}
sub five_splice_site_canonical {
    my ($genome_target_ref, $mismatch, $gap) = @_;  
    #print Dumper ($genome_target_ref);
    my @genome_target = @$genome_target_ref;
    my $five_splicesite; #1 if splicesite is found, else 0
    my $five_nt = -3; #nt between splice site and query exon

	for(my $i = ((scalar @genome_target) - $extension -3 );$i <= (scalar @genome_target) - 2;$i++){	
		#print Dumper ($genome_target[$i]);
        if ($genome_target[$i] eq "G" and $genome_target[$i+1] eq "T"){     

        	$five_splicesite = 1;    
        	return ($five_splicesite, $five_nt);
        } 
        else {
        	$five_splicesite = 0;
        	$five_nt++;
        }
    }
    return ($five_splicesite, $five_nt);
}

sub five_splice_site_canonical_after_negative_score {
    my ($genome_target_ref, $five_nt) = @_;  
    my @genome_target = @$genome_target_ref;
    my $five_splicesite = 0; #1 if splicesite is found, else 0

	for(my $i = ((scalar @genome_target) - $extension + 2 + $five_nt);$i <= (scalar @genome_target) - 2;$i++){	
		# print Dumper ($five_nt);
        if ($genome_target[$i] eq "G" and $genome_target[$i+1] eq "T"){     
        	# print Dumper ("yes");
        	$five_splicesite = 1;    
        	return ($five_splicesite, ($five_nt + 2));
        } 
        else {
        	# print Dumper ("no");
        	$five_splicesite = 0;
        	$five_nt++;
        }
    }
    return ($five_splicesite, $five_nt);
}

sub five_splice_site_noncanonical {   
	my ($genome_target_ref, $five_nt) = @_;  
    my @genome_target = @$genome_target_ref;
    my $five_splicesite; #1 if splicesite is found, else 0
    #my $five_nt = -3; #nt between splice site and query exon

	for(my $i = ((scalar @genome_target) - $extension +2 + $five_nt);$i <= (scalar @genome_target) - 2;$i++){
	#for(my $i = ((scalar @genome_target) - $extension -3);$i <= (scalar @genome_target) - 2;$i++){	
        if ($genome_target[$i] eq "G" and $genome_target[$i+1] eq "C"){     

        	$five_splicesite = 1;    
        	return ($five_splicesite, $five_nt + 2);
        } 
        else {
        	$five_splicesite = 0;
        	$five_nt++;
        }
    }
    return ($five_splicesite, $five_nt);
}

sub five_splice_site_canonical_intron {
    my (@genome_target) = @_;  
    my $five_splicesite; #1 if splicesite is found, else 0
    my $five_nt = 0; #nt between splice site and query exon

	for(my $i = 0;$i <= (scalar @genome_target) -3;$i++){	
		#print Dumper ($genome_target[$i]);
        if ($genome_target[$i] eq "G" and $genome_target[$i+1] eq "T"){     

        	$five_splicesite = 1;    
        	return ($five_splicesite, $five_nt);
        } 
        else {
        	$five_splicesite = 0;
        	$five_nt++;
        }
    }
    return ($five_splicesite, $five_nt);
}

sub sequence_three_splicesite{						# score3ss: sequence 23 bases long. [20 intron][3 exon]
	my ($genome_target_ref, $three_nt) = @_;      
	my @genome_target = @$genome_target_ref;
	my @testsequence_three_splicesite;
	for(my $i = ($extension -20 - $three_nt);$i <= ($extension -20 - $three_nt + 22);$i++) {
		# print Dumper ($i);
		# print Dumper ($genome_target[$i]);
		push (@testsequence_three_splicesite, $genome_target[$i]);
		}

		# print Dumper (@testsequence_three_splicesite);
	return (@testsequence_three_splicesite);
}

sub sequence_five_splicesite{						# score5ss: sequence 9 bases long. [3 exon][6 intron] 
	my ($genome_target_ref, $five_nt) = @_;      
	my @genome_target = @$genome_target_ref;
	my @testsequence_five_splicesite;	
	# print Dumper (scalar @genome_target);
	# print Dumper ($five_nt);
	for(my $i = (scalar @genome_target - $extension -3 + $five_nt);$i <= (scalar @genome_target - $extension + 5 + $five_nt);$i++) {
		# print Dumper ($i);
		# print Dumper ($genome_target[$i]);
		push (@testsequence_five_splicesite, $genome_target[$i]) if defined $genome_target[$i];
		}
	return (@testsequence_five_splicesite);
}

sub sequence_final_exon{
	my ($genome_target, $three_nt, $five_nt) = @_;      
	my @genome_target_array = split //, $genome_target;
#print Dumper ($genome_target, $three_nt, $five_nt);
	my @sequence_final_exon;
	for(my $i = ($extension - $three_nt);$i <= (scalar @genome_target_array - $extension -1  + $five_nt);$i++) {
		#print Dumper ($i);
		#print Dumper ($genome_target_array[$i]);
		push (@sequence_final_exon, $genome_target_array[$i]) if defined $genome_target_array[$i];
		}

		#print Dumper (@sequence_final_exon);
	return (@sequence_final_exon);
}