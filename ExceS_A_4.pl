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
my $target;
my $strandIn;
my $core = 8;
my $extension = 35;
my $query;
my %model_informations;
my $myCodonTable = Bio::Tools::CodonTable->new();
my $paralog;
my $homologousexons =  "HomologousExons/paralog_exon_alignment.fa";
my $mode = "fasta";

my $path_to_maxentscan_3 = "/scratch/franziska/bin/maxentscan/fordownload/score3.pl";
my $path_to_maxentscan_5 = "/scratch/franziska/bin/maxentscan/fordownload/score5.pl";

GetOptions("target=s" => \$target, "query=s" => \$query, "strand=s" => \$strandIn, "paralog=s" => \$paralog, "mode=s" => \$mode, "homologousexons=s" => \$homologousexons);

#------------------------------------------------------------------------
#Main
`rm -r MaxEntScan`;
mkdir "MaxEntScan/";

my @paralog = split "_", $paralog;
my ($type, $querypara) = split "-", $paralog[0];

my %exonsequence;    
my $seqio = Bio::SeqIO->new(-file => "$homologousexons", -format => "fasta");
while(my$seqobj = $seqio->next_seq) {
    my $id  = $seqobj->display_id;
    my $seq = $seqobj->seq;
    $exonsequence{$id} = $seq;
}

my $blat = "SearchTarget/Target_vs_query.blatout_${paralog}_merged";
my @call_blat = "blat $target $homologousexons $blat -noHead -q=prot -t=dnax -out=pslx -maxIntron=20000";
system (@call_blat) == 0 or die "blat against contig genome failed. $!\n";

open (my $blatout_handle, "<", $blat) or die "Can't open Blatout, $!";
my @a;

my @b;
my $a_align;
my $b_align;

#data from blatout
while (defined (my $line = <$blatout_handle>)){
	chomp $line;	
	#print Dumper ($line);
	next if ($line =~ m/^\#/ || $line =~ m/^$/);  # skip comments and empty lines
	if ($line =~ m/${querypara}/){
		#print Dumper ($line);
		my @line_split = split /\t/, $line;
		if ($line_split[8] eq "+${strandIn}"){
			my %data = search_data($line);

			my $genome_seqIO = Bio::SeqIO->new(-file   => $target, -format => "Fasta");
			$a = $genome_seqIO->next_seq;
			$a = $a -> seq;

			my @a_short;			

			if (($data{start}) - $extension < 0){
				$data{extension_start} = ($extension - $data{start});
			}
			else{
				$data{extension_start} = $extension;
			}
			
			if ($data{exonID} == 0){
				if ($data{start_query} > 10){
					$data{extension_start} = $extension + ($data{start_query} * 3)
				}
			}


			my $extended_start = $data{start} -$data{extension_start} ;

			#print Dumper ((length $a),$data{end});	
			if (((length $a) - $data{end}) < $extension){
				$data{extension_end} = ((length $a) - $data{end} -1);
			}
			else{
				$data{extension_end} = $extension -1;
			}
			#print Dumper ($data{extension_end});
			my $extended_end = $data{end} + $data{extension_end};


			if (($data{end_query} < ($data{query_size} - 10)) && (((length $a) - $data{end}) > $extension)){
				my $end_extension = (($data{query_size} - $data{end_query}) *3);
					
				if ($data{strand} eq "+" ){
					$extended_end += $end_extension;
				}
				else {
					$extended_start -= $end_extension;
				}
			}
			#print Dumper ($data{start},$data{end}, $extended_start, $extended_end);
			if ($data{strand} eq "+"){		
				@a = split //, $a;
				for(my $i = $extended_start;$i <= $extended_end;$i++) {
					#print Dumper ($i,$a[$i]);
					push (@a_short, $a[$i]);
				}
			}
			elsif($data{strand} eq "-"){
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
			#print Dumper ((scalar @a));
			$b = $exonsequence{$data{model}};
			my @b = split //, $b;
			$data{genome_part} = join ("",@a_short);
			$data{query} = $b;
			$data{blat} = "yes";
			$data{three_nt_noIntron} = 0;
			my $key_for_hash = $data{model}."_".$data{start};
			$model_informations{$key_for_hash} = \%data;
		}
	}					
}
close $blatout_handle;

#print Dumper(%model_informations);

if( !%model_informations){
	exit;
}
for my $i (keys %model_informations) {
	if($model_informations{$i}->{'exonID'} > 0){
		if ($model_informations{$i}->{'start_query'} > 12 || $model_informations{$i}->{'end_query'} < $model_informations{$i}->{'query_size'} -10){
			delete $model_informations{$i};
		}
	}
}

#my $pm = new Parallel::ForkManager($core);

if( !%model_informations){
	exit;
}
my @exons;

`grep \'>${querypara}_\' $homologousexons > exons.txt`;

open (my $exon_handle, "< exons.txt") or die "Can't open exons.txt, $!";
while (defined (my $line = <$exon_handle>)){
	chomp $line;
	my $line_2 = (split "_", $line)[1];
	$line_2 =~ s/exon//;
	push (@exons, $line_2);
}
my @sorted_exons = sort {$a <=> $b} @exons;
my @unique_exons = uniq @sorted_exons;
my $lastexon=$sorted_exons[-1];
#print Dumper(@sorted_exons);
my @exons_found;
for my $i (keys %model_informations) {
	push (@exons_found, $model_informations{$i}->{'exonID'});
}	

my @sorted_exons_found = sort {$a <=> $b} @exons_found;
#print Dumper(@sorted_exons_found);
my @keys = (natsort keys %model_informations);
#print Dumper(@keys);

my $model;
my $start1;
my $end1;
my $last_start;
my $last_end;
my $last_model;

if (scalar(@keys) == 1){
	$model = $keys[0];
	$model_informations{$model}->{'five_splicesite'} = 0;
	$model_informations{$model}->{'three_splicesite'} = 0;
}

for (my $i = 1; $i < (scalar @keys);$i++) {
	if ($strandIn eq "+"){		
		$model = $keys[$i];
		#print Dumper ($keys[$i]);
		$start1 = $model_informations{$model}->{'start'};
		$end1 = $model_informations{$model}->{'end'};
		$last_model = $keys[$i-1];
		$last_start = $model_informations{$last_model}->{'start'};
		$last_end = $model_informations{$last_model}->{'end'};
		if (($start1 - $last_end) <= 3){			
			#print Dumper("yes");
			$model_informations{$model}->{'three_splicesite'} = 1;
			$model_informations{$model}->{'three_nt'} = 0;
			$model_informations{$model}->{'three_nt_noIntron'} = ($start1 - $last_end);
			$model_informations{$model}->{'score_3splicesite'} = 100;
			$model_informations{$model}->{'five_splicesite'} = 0;
			$model_informations{$last_model}->{'five_splicesite'} = 1;
			$model_informations{$last_model}->{'five_nt'} = 0;
			$model_informations{$last_model}->{'score_5splicesite'} = 100;
			$model_informations{$last_model}->{'three_splicesite'} = 0;
		}
		if (!$model_informations{$last_model}->{'three_splicesite'}){$model_informations{$last_model}->{'three_splicesite'} = 0;}
		if (!$model_informations{$last_model}->{'five_splicesite'}){ $model_informations{$last_model}->{'five_splicesite'} = 0;}
		if (!$model_informations{$model}->{'five_splicesite'}){ $model_informations{$model}->{'five_splicesite'} = 0;}
		if (!$model_informations{$model}->{'three_splicesite'}){ $model_informations{$model}->{'three_splicesite'} = 0;}

	}	
	elsif ($strandIn eq "-"){
		$model = $keys[$i];
		#print Dumper ($keys[$i]);
		$start1 = $model_informations{$model}->{'start'};
		$end1 = $model_informations{$model}->{'end'};
		$last_model = $keys[$i-1];
		$last_start = $model_informations{$last_model}->{'start'};
		$last_end = $model_informations{$last_model}->{'end'};
		if (($last_start - $end1) <= 3){
			#print Dumper("yes");
			$model_informations{$model}->{'three_splicesite'} = 1;
			$model_informations{$model}->{'three_nt'} = 0;
			$model_informations{$model}->{'three_nt_noIntron'} = ($last_start - $end1);
			$model_informations{$model}->{'score_3splicesite'} = 100;
			$model_informations{$model}->{'five_splicesite'} = 0;
			$model_informations{$last_model}->{'five_splicesite'} = 1;
			$model_informations{$last_model}->{'five_nt'} = 0;
			$model_informations{$last_model}->{'score_5splicesite'} = 100;
			$model_informations{$last_model}->{'three_splicesite'} = 0;
		}
		if (!$model_informations{$last_model}->{'three_splicesite'}){$model_informations{$last_model}->{'three_splicesite'} = 0;}
		if (!$model_informations{$last_model}->{'five_splicesite'}){ $model_informations{$last_model}->{'five_splicesite'} = 0;}
		if (!$model_informations{$model}->{'five_splicesite'}){ $model_informations{$model}->{'five_splicesite'} = 0;}
		if (!$model_informations{$model}->{'three_splicesite'}){ $model_informations{$model}->{'three_splicesite'} = 0;}
	}
}

my @diff = array_diff(@unique_exons, @sorted_exons_found);
#print Dumper(@diff,($#diff));
my @missing_exons; 
my $j = 0;
my %missing_exon_hash;


if (defined $diff[0]){
	for my $i (0..$#diff){
		#print Dumper("$i: ".$i); 
		push (@missing_exons, $diff[$i]);

		if ($diff[$i] == 0){
			$missing_exon_hash{$j}->{'start_exon'} = 0;
			$missing_exon_hash{$j}->{'end_exon'} = 0;
		}
		elsif($diff[$i] == $diff[$i-1]+1){
			$missing_exon_hash{$j}->{'end_exon'} = $diff[$i];
		}
		else{
			$j = $i;
			$missing_exon_hash{$j}->{'start_exon'} = $diff[$i];
			$missing_exon_hash{$j}->{'end_exon'} = $diff[$i];
		}
	}
}

#print Dumper(%missing_exon_hash);

for my $i (keys %missing_exon_hash){
	#print Dumper($i, @keys);
	my $start;
	my $end;
	my $paralogID;
	if ($strandIn eq "+"){
		#print Dumper($missing_exon_hash{$i}->{'start_exon'});
		if ($missing_exon_hash{$i}->{'start_exon'} == 0){
			$start = 0;
		}
		else{		
			my $searched_end_exon = $missing_exon_hash{$i}->{'start_exon'}-1;
			#print Dumper($searched_end_exon);
			for( @keys ){
					if( $_ =~ m/exon${searched_end_exon}/){
					$start = $model_informations{$_}->{'end'} -50;
					#print Dumper($start);
					$paralogID = $_;
					last;
					}
			}
		}
		if ($missing_exon_hash{$i}->{'end_exon'} == $lastexon){
			$end = (scalar @a);
			#print Dumper($end);
		}
		else{						
			my $searched_start_exon = $missing_exon_hash{$i}->{'end_exon'}+1;
			for( @keys ){
				if( $_ =~ m/exon${searched_start_exon}/){
					$end = $model_informations{$_}->{'start'} +50;
					#print Dumper($end);
					$paralogID = $_;
					last;
				}
			}	
		}
	}	
	elsif($strandIn eq "-"){
		#print Dumper($missing_exon_hash{$i}->{'end_exon'});
		if ($missing_exon_hash{$i}->{'end_exon'} == $lastexon){
			$start = 0;
		}
		else{		
			my $searched_end_exon = $missing_exon_hash{$i}->{'end_exon'}+1;
			#print Dumper($searched_end_exon);
			for( @keys ){
					if( $_ =~ m/exon${searched_end_exon}/){
					$start = $model_informations{$_}->{'end'} -50;
					#print Dumper($start);
					$paralogID = $_;
					last;
					}
			}
		}

		if ($missing_exon_hash{$i}->{'start_exon'} == 0){
			$end = (scalar @a);
			#print Dumper($end);
		}
		else{						
			my $searched_start_exon = $missing_exon_hash{$i}->{'start_exon'}-1;
			#print Dumper ($searched_start_exon);
			for( @keys ){
				#print Dumper ("muh",$_ );
				if( $_ =~ m/exon${searched_start_exon}/){
					#print Dumper ("mäh");
					$end = $model_informations{$_}->{'start'} +50 ;
					#print Dumper($end);
					$paralogID = $_;
					last;
				}
			}

			if (!defined $end){
				$searched_start_exon = $missing_exon_hash{$i}->{'start_exon'}-2;
				if( $_ =~ m/exon${searched_start_exon}/){
					$end = $model_informations{$_}->{'start'} +50 ;
					#print Dumper($end);
					$paralogID = $_;
					last;
				}
			}	
		}
	}	

	if($start < 0){
		$start = 0;
	}

	if($end > (length $a)){
		$end = (length $a);
	}

	unless (-s "TargetSequences/${paralog}_${start}_${end}.fa"){ 
			`formatdb -i $target -p F -o`;
			`fastacmd -d $target -s ${paralog}_merged -L ${start},${end} | sed \'s/>lcl|/>/g\' > TargetSequences/${paralog}_${start}_${end}.fa`;
			`formatdb -i TargetSequences/${paralog}_${start}_${end}.fa -p F -o`;	
	}
	#print Dumper($paralogID,$model_informations{$paralogID});

	for my $missing_exon (@missing_exons){
		`grep -A1 \'$model_informations{$paralogID}->{'paralogID'}_exon${missing_exon}\' $homologousexons --no-group-separator >> HomologousExons/missing_exons_splitaligner_${paralog}_merged.fa`;
	}

	my $blat3 = "SearchTarget/Target_vs_query.blastout_${paralog}_merged";
	my @call_blast = "blastall -p tblastn -d TargetSequences/${paralog}_${start}_${end}.fa -i HomologousExons/missing_exons_splitaligner_${paralog}_merged.fa -e 0.01 -a 10 -m8 -G 11 -E 1 -C F > $blat3";
	system (@call_blast) == 0 or die "blat against contig genome failed. $!\n";
	#print "@call_blast\n";

	`rm -r HomologousExons/missing_exons_splitaligner_${paralog}_merged.fa`;

	if (-e $blat3){
		open (my $blatout_handle, "<", $blat3); #or die "Can't open Blatout, $!";
		my @a;
		my @b;
		my $a_align;
		my $b_align;

		#data from blatout
		while (defined (my $line = <$blatout_handle>)){
			chomp $line;	
			next if ($line =~ m/^\#/ || $line =~ m/^$/);  # skip comments and empty lines
			my %data;
			my $counter = 0;
			if (($line =~ m/${querypara}/) && $counter == 0){
				my @line_split = split /\t/, $line;	
				if ($line_split[2] > 50.00){	
					%data = search_data_blast($line);
					$counter = 1;
				}
				else{
					next;
				}
				#print Dumper ($line, $data{start},$data{end});
				#print "read sequences for alignment\n";
				my $genome_seqIO = Bio::SeqIO->new(-file   => "TargetSequences/${paralog}_${start}_${end}.fa", -format => "Fasta");
				$a = $genome_seqIO->next_seq;
				$a = $a -> seq;
				#print Dumper ($line, $data{start},$data{end}, (length $a));
				my @a_short;
				my $extended_start;
				my $extended_end;

				if ($data{strand} eq "+"){
					if (($data{start}) - 1 - $extension < 0){
						#print Dumper ("1");
					$data{extension_start} = ($data{start} - 1);
					}
					else{
						#print Dumper ("2");
						$data{extension_start} = $extension;
					}
					
					$extended_start = $data{start} -1 -$extension;

					if (((length $a) - $data{end}) -1 < $extension){
						#print Dumper ("3");
						$data{extension_end} = ((length $a) - $data{end} -1);
					}
					else{
						#print Dumper ("4");
						$data{extension_end} = $extension -1;
					}
					$extended_end = $data{end} + $data{extension_end} ;


					if (($data{end_query} < ($data{query_size} - 10)) && (((length $a) - $data{end} -1) > $extension)){
						#print Dumper ("5");
						my $end_extension = (($data{query_size} - $data{end_query}) *3);	
						$extended_end += $end_extension;
					}
				}
				else{
					if (($data{end}) - 1 - $extension < 0){
						#print Dumper ("1");
						$data{extension_end} = ($data{end} - 1);
					}
					else{
						#print Dumper ("2");
						$data{extension_end} = $extension -1;
					}
					
					$extended_start = $data{end} -2 - $data{extension_end};
					#print Dumper($extended_start, $data{extension_end} );
					if (((length $a) - $data{start}) -1 < $extension){
						#print Dumper ("3");
						$data{extension_start} = ((length $a) - $data{start}-1);
					}
					else{
						#print Dumper ("4");
						$data{extension_start} = $extension;
					}
					$extended_end = $data{start} -1 + $data{extension_start};
					#print Dumper($extended_end, $data{extension_start} );

					if (($data{end_query} < ($data{query_size} - 10)) && (((length $a) - $data{start} -1) > $extension)){
						#print Dumper ("5");
						my $end_extension = (($data{query_size} - $data{end_query}) *3);
						$extended_start -= $end_extension;
					}
				}

				if ($data{strand} eq "+"){		
					@a = split //, $a;
					#print Dumper ("6", $a);
					for(my $i = $extended_start;$i <= $extended_end;$i++) {
						#print Dumper ($i, $a[$i]);
						push (@a_short, $a[$i]);
					}
				}
				elsif($data{strand} eq "-"){
					#print Dumper ("6", $a);
					@a = split //, $a;
					for(my $i = $extended_start;$i <= $extended_end;$i++) {
						#print Dumper ($i, $a[$i]);
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

				$b = $exonsequence{$data{model}};
				my @b = split //, $b;
				$data{genome_part} = join ("",@a_short);
				#print Dumper ($data{genome_part});
				$data{query} = $b;
				$data{three_splicesite} = 0,
				$data{five_splicesite} = 0,
				$data{blat} = "no";
				$data{three_nt_noIntron} = 0;
				my $key_for_hash = $data{model}."_".$data{start};
				#print Dumper(%data);
				$model_informations{$key_for_hash} = \%data;				
			}
		}
	}
}

#%model_informations = find_best_exon_mapping(\@unique_exons, \%model_informations);

for my $i (keys %model_informations) {
	#print Dumper ($model_informations{$i});
	#$pm->start and next;
	my @a_short = split //, $model_informations{$i}->{'genome_part'};

	#search stopcodon next to last exon
	if ($model_informations{$i}->{'exonID'} == $lastexon){
		$model_informations{$i}->{'stopcodon'} = 0;
		my ($stopcodon, $stop_nt) = find_stopcodon(\@a_short, $model_informations{$i}->{'extension_end'});	
		$model_informations{$i}->{'stop_nt'}  = $stop_nt;		
		$model_informations{$i}->{'modulo_stop'} = $model_informations{$i}->{'stop_nt'} % 3;

		if ($stopcodon == 1 && ($model_informations{$i}->{'modulo_stop'} == 0)){
			$model_informations{$i}->{'score_5splicesite'} = 0;
			$model_informations{$i}->{'stopcodon'} = 1;
			$model_informations{$i}->{'startcodon'} = 0;	
		}

		else {		
			until ($model_informations{$i}->{'modulo_stop'} == 0 || $model_informations{$i}->{'stop_nt'} >= $model_informations{$i}->{'extension_end'} -5){
				#print Dumper ($model_informations{$i}->{'stop_nt'});
				($stopcodon, $stop_nt) = find_stopcodon_after_modulo(\@a_short,$model_informations{$i}->{'stop_nt'}, $model_informations{$i}->{'extension_end'});	
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
				print "no stopcodon: $model_informations{$i}->{'model'}\n";
				delete $model_informations{$i};
				next;
			}
		}

		$model_informations{$i} = best_three_splicesite($model_informations{$i});
	}

	#search startcodon in exon 0
	elsif ($model_informations{$i}->{'exonID'} eq "0"){
		$model_informations{$i}->{'startcodon'} = 0;
		my ($startcodon, $start_nt) = find_startcodon(\@a_short, $model_informations{$i}->{'start_query'}, $model_informations{$i}->{'extension_start'});
		$model_informations{$i}->{'start_nt'}  = $start_nt;	
		#print Dumper ($model_informations{$i}->{'start_nt'});
		$model_informations{$i}->{'modulo_start'} = $model_informations{$i}->{'start_nt'} % 3;

		if ($startcodon == 1 and $model_informations{$i}->{'modulo_start'} == 0){
			$model_informations{$i}->{'score_3splicesite'} = 0;
			$model_informations{$i}->{'startcodon'} = 1;
			$model_informations{$i}->{'stopcodon'} = 0;
		}
		else {
			until ($model_informations{$i}->{'modulo_start'} == 0 || $model_informations{$i}->{'start_nt'} >= $model_informations{$i}->{'extension_start'} -4){
				my ($startcodon, $start_nt) = find_startcodon_after_modulo(\@a_short, $model_informations{$i}->{'start_nt'}, $model_informations{$i}->{'extension_start'});
				$model_informations{$i}->{'start_nt'}  = $start_nt;		
				$model_informations{$i}->{'modulo_start'} = $model_informations{$i}->{'start_nt'} % 3;
				#print Dumper ($model_informations{$i});
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
				print "no startcodon: $model_informations{$i}->{'model'}\n";
				delete $model_informations{$i};
				next;
				#print Dumper ($model_informations{$i});
			}
		}

		$model_informations{$i} = best_five_splicesite($model_informations{$i});
	}

	# schauen, ob GT....AG (canonical) im Intron ist, GC-AG (noncanonical)
	# mammalian non-canonical splice sites, then the 99.24% of splice site pairs should be GT-AG, 0.69% GC-AG, 0.05% AT-AC and finally only 0.02% could consist of other types of non-canonical splice sites
	# (Analysis of canonical and non-canonical splice sites in mammalian genomes; M. Burset; 2000)

	else{	
		#print Dumper ("möp",$model_informations{$i});

		$model_informations{$i} = best_five_splicesite($model_informations{$i});
		#print Dumper ("mäp");
		$model_informations{$i} = best_three_splicesite($model_informations{$i});		
		#print Dumper ("mäp",$model_informations{$i});
	}
		#$pm->finish
}
#$pm->wait_all_children;

#maxentscan output positiv
# -> wie gross ist $nt?
#	-> <3, dann schauen, ob gesplittete as
#print Dumper (%model_informations);

my $model_informations_ref = \%model_informations;
#print Dumper(%model_informations);
%model_informations = splitted_as($model_informations_ref);

#print Dumper (%model_informations);
%model_informations = matching_neighbours($model_informations_ref);

%model_informations = splitted_as($model_informations_ref);


# fasta output for every paralog 
my $dna;
my @targetsequence;
my @sorted_exons_found2;
my @exons_found2;
			
for my $entry (natsort keys %model_informations){		
	if (defined ($model_informations{$entry}->{'score_3splicesite'} and $model_informations{$entry}->{'score_5splicesite'} and $model_informations{$entry}->{'stopcodon'} and $model_informations{$entry}->{'startcodon'})){
		if (($model_informations{$entry}->{'score_3splicesite'} > -10 and $model_informations{$entry}->{'score_5splicesite'} > -10) or ($model_informations{$entry}->{'score_3splicesite'} > -10 and $model_informations{$entry}->{'stopcodon'} == 1) or ($model_informations{$entry}->{'startcodon'} == 1 and $model_informations{$entry}->{'score_5splicesite'} > -10)){
			#print Dumper($model_informations{$entry}->{'exonID'});
			push (@exons_found2, $model_informations{$entry}->{'exonID'});
			@sorted_exons_found2 = sort {$a <=> $b} @exons_found2;
		}
	}
}

#print Dumper (@sorted_exons_found2);
#print Dumper(%model_informations);

my ($new_three_nt, $new_five_nt);
for my $entry (natsort keys %model_informations){
	if (defined ($model_informations{$entry}->{'score_3splicesite'} and $model_informations{$entry}->{'score_5splicesite'} and $model_informations{$entry}->{'stopcodon'} and $model_informations{$entry}->{'startcodon'})){
		if (($model_informations{$entry}->{'score_3splicesite'} > -10 and $model_informations{$entry}->{'score_5splicesite'} > -10) or ($model_informations{$entry}->{'score_3splicesite'} > -10 and $model_informations{$entry}->{'stopcodon'} == 1) or ($model_informations{$entry}->{'startcodon'} == 1 and $model_informations{$entry}->{'score_5splicesite'} > -10)){
			
			#print Dumper($model_informations{$entry});
			if ($model_informations{$entry}->{'exonID'} == 0){
				if (any {$_ == ($model_informations{$entry}->{'exonID'} +1)} @sorted_exons_found){
					$new_three_nt = $model_informations{$entry}->{'start_nt'} + 3;
					$new_five_nt = $model_informations{$entry}->{'five_nt'};
					#print Dumper("1");
				}
				else{
					$new_three_nt = $model_informations{$entry}->{'start_nt'} + 3;
					$new_five_nt = $model_informations{$entry}->{'five_nt_extension'};
					#print Dumper("2");
				}
			}
			elsif($model_informations{$entry}->{'exonID'} == $lastexon){
				if (any {$_ == ($model_informations{$entry}->{'exonID'} -1)} @sorted_exons_found2) {
					$new_three_nt = $model_informations{$entry}->{'three_nt'} + $model_informations{$model}->{'three_nt_noIntron'};
					$new_five_nt = $model_informations{$entry}->{'stop_nt'};
					#print Dumper("3");
				}
				else{
					$new_three_nt = $model_informations{$entry}->{'three_nt_extension'} + $model_informations{$model}->{'three_nt_noIntron'};
					$new_five_nt = $model_informations{$entry}->{'stop_nt'};
					#print Dumper("4");
				}
			}
			else{
				if ((any {$_ == ($model_informations{$entry}->{'exonID'} +1)} @sorted_exons_found2) && (any {$_ == ($model_informations{$entry}->{'exonID'} -1)} @sorted_exons_found2)) {
					$new_three_nt = $model_informations{$entry}->{'three_nt'} + $model_informations{$model}->{'three_nt_noIntron'};
					$new_five_nt = $model_informations{$entry}->{'five_nt'};
					#print Dumper("5");
				}
				elsif(any {$_ == ($model_informations{$entry}->{'exonID'} +1)} @sorted_exons_found2){
					$new_three_nt = $model_informations{$entry}->{'three_nt_extension'} + $model_informations{$model}->{'three_nt_noIntron'};
					$new_five_nt = $model_informations{$entry}->{'five_nt'};
					#print Dumper("6");
				}
				elsif(any {$_ == ($model_informations{$entry}->{'exonID'} -1)} @sorted_exons_found2){
					$new_three_nt = $model_informations{$entry}->{'three_nt'} + $model_informations{$model}->{'three_nt_noIntron'};
					$new_five_nt = $model_informations{$entry}->{'five_nt_extension'};
					#print Dumper("7");
				}
				else{
					$new_three_nt = $model_informations{$entry}->{'three_nt_extension'} + $model_informations{$model}->{'three_nt_noIntron'};
					$new_five_nt = $model_informations{$entry}->{'five_nt_extension'};
					#print Dumper("8");
				}
			}
			if ($model_informations{$entry}->{'blocks'} > 1){
				my @new_genome_part;
				my @blocks = split /,/, $model_informations{$entry}->{'blocksizes'};
				my @qstarts = split /,/, $model_informations{$entry}->{'qstarts'};

				if ($model_informations{$entry}->{'blocks'} == 2){
					# my $block_1 = $blocks[0];
					# my $block_2 = $blocks[1];
					# my $delete_aa_end_1 = ($model_informations{$entry}->{'extension_start'} + ($block_1 * 3));
					# my $delete_aa_start_2 = ((length $model_informations{$entry}->{'genome_part'}) - $model_informations{$entry}->{'extension_end'} - ($block_2 * 3));	
					# my $gap = abs($delete_aa_start_2 - $delete_aa_end_1);
					# print Dumper ($gap,$delete_aa_end_1, $delete_aa_start_2);


					my $block_1 = $blocks[0];
					my $block_2 = $blocks[1];
					my $delete_aa_end_1 = ($model_informations{$entry}->{'extension_start'} + ($block_1 * 3));
					my $delete_aa_start_2 = ((length $model_informations{$entry}->{'genome_part'}) - $model_informations{$entry}->{'extension_end'} - ($block_2 * 3));	
					my $gap = abs($delete_aa_start_2 - $delete_aa_end_1 -1);
					#print Dumper ($blocks[0],$blocks[1],$model_informations{$entry}->{'extension_start'},$model_informations{$entry}->{'extension_end'});
					#print Dumper ($gap,$delete_aa_end_1, $delete_aa_start_2);
					#$VAR1 = -22;
					#$VAR2 = 170;
					#$VAR3 = 148;

					if ($gap <=3){
					substr($model_informations{$entry}->{'genome_part'}, ($delete_aa_end_1 +1), $gap) = "";
					}
					else{
						my $testsequence = substr($model_informations{$entry}->{'genome_part'}, ($delete_aa_end_1), $gap + 3);
						#print Dumper ($testsequence);
						my @testsequence_array = split //, $testsequence;
						my $three_nt = 0;
						my $three_splicesite = 0;
						my $score_3splicesite = 0;

						until ($score_3splicesite > 0 || $three_nt >= $gap){
							
							($three_splicesite, $three_nt) = three_splice_site_intron(\@testsequence_array, $three_nt);
							if ($three_splicesite == 1){
								my $testsequence_three_splicesite = substr($model_informations{$entry}->{'genome_part'}, ($delete_aa_end_1 + $gap + 3 - $three_nt - 20), 23);	
								#print Dumper ($testsequence_three_splicesite);
								my @testsequence_three_splicesite = split //, $testsequence_three_splicesite;
								if (scalar @testsequence_three_splicesite == 23){
									open (testsequences_three_splicesite_file, '>>', "MaxEntScan/Testsequences_three_splicesite.txt") or die "$!";
									print testsequences_three_splicesite_file join ("",@testsequence_three_splicesite);
								}

								if(-e "MaxEntScan/Testsequences_three_splicesite.txt"){
									my $output_3splicesite = qx(perl $path_to_maxentscan_3 "MaxEntScan/Testsequences_three_splicesite.txt");
									#print Dumper ($output_3splicesite);
									my @array_3splicesite = split /\t/, $output_3splicesite;
									if(defined $array_3splicesite[1]){
										chomp ($array_3splicesite[1]);
									}
									$score_3splicesite = $array_3splicesite[1];	
									# my $key_for_hash = $model_informations->{score_3splicesite};
									# $three_scores{$key_for_hash} = $model_informations->{three_nt};
									unlink "MaxEntScan/Testsequences_three_splicesite.txt";
								}
								
								if ($score_3splicesite < 0){
									$three_nt += 2;
								}	
							}
						}

						if ($score_3splicesite < 0){
							$model_informations{$entry}->{'genome_part'} = $model_informations{$entry}->{'genome_part'};
						}
						else{
							my ($five_splicesite, $five_nt) = five_splice_site_canonical_intron(@testsequence_array);
							#print Dumper ($three_splicesite,$five_splicesite, $three_nt, $five_nt);
							if ($three_splicesite == 1 && $five_splicesite == 1){
								substr($model_informations{$entry}->{'genome_part'}, ($delete_aa_end_1 +1 + $five_nt), $gap + 3 -$three_nt - $five_nt) = "";
							}
						}
					}
				}
				if ($model_informations{$entry}->{'blocks'} >2){
					#print Dumper ("muh");	
					my $firstblock_sequence;
					my $lastblock_sequence;
					my $newquery;
					#my @testsequence_array = split //, $testsequence;
					my @a_short = split //, $model_informations{$entry}->{'genome_part'};
					my ($five_splicesite, $five_nt) = five_splice_site_canonical_block(\@a_short, $blocks[0]);
					#print Dumper ($five_splicesite);
					if ($five_splicesite == 1){
						$firstblock_sequence = substr($model_informations{$entry}->{'genome_part'}, 0, ($model_informations{$entry}->{'extension_start'} + ($blocks[0] * 3)+ $five_nt));				
						my ($three_splicesite, $three_nt) = three_splice_site_block(\@a_short, $blocks[-1]);
						#print Dumper ($three_splicesite);	
						if ($three_splicesite == 1){
							$lastblock_sequence = substr($model_informations{$entry}->{'genome_part'}, ((length $model_informations{$entry}->{'genome_part'}) - $model_informations{$entry}->{'extension_end'} - ($blocks[-1] * 3) - $three_nt),(length $model_informations{$entry}->{'genome_part'}));

							my $three_nt_extension_aa = POSIX::floor($three_nt / 3);
							my $five_nt_extension_aa = POSIX::floor($five_nt / 3);

							$newquery = substr($model_informations{$entry}->{'query'}, ($blocks[0] + $five_nt_extension_aa),((length $model_informations{$entry}->{'query'}) - $qstarts[-1] - $three_nt_extension_aa));
							#print Dumper((length $model_informations{$entry}->{'query'}), $blocks[-1],$newquery);
							open (NEWQUERY, ">SearchTarget/${entry}_new.fa");	
							print NEWQUERY ">${entry}\n";
							print NEWQUERY "${newquery}\n";
							open (OUTTARGET, ">TargetSequences/${paralog}_${entry}_block.fa");	
							print OUTTARGET ">${entry}\n";
							print OUTTARGET "$model_informations{$entry}->{'genome_part'}\n";


							my $blat = "SearchTarget/Target_vs_query.blatout_${entry}";
							#my @call_blat = "blat TargetSequences/${contigIn}_${entry}_block.fa SearchTarget/${entry}_new.fa $blat -noHead -q=prot -t=dnax -out=pslx -maxIntron=90";
							#system (@call_blat) == 0 or die "blat against contig genome failed. $!\n";
							`blat TargetSequences/${paralog}_${entry}_block.fa SearchTarget/${entry}_new.fa $blat -noHead -q=prot -t=dnax -out=pslx -maxIntron=90`;
							open (my $blatout_handle, "<", $blat) or die "Can't open Blatout, $!";

							while (defined (my $line = <$blatout_handle>)){
								chomp $line;	
								my @line_array = split/\t/, $line;
								if ($line_array[17] == 1){
									my $midblock_testsequence = substr($model_informations{$entry}->{'genome_part'}, ($line_array[15] -  $extension), ($line_array[16] +  $extension));
									#print Dumper ($midblock_testsequence);
									my @midblock_testsequence = split //, $midblock_testsequence;
									my ($five_splicesite, $five_nt) = five_splice_site_canonical(\@midblock_testsequence, 0, 0,$extension);
									my ($three_splicesite, $three_nt) = three_splice_site(\@midblock_testsequence, 0, -1,$extension);
									#print Dumper ($three_splicesite,$five_splicesite);
									if ($three_splicesite == 1 && $five_splicesite == 1){
										if(($line_array[16] + $five_nt) < ((length $model_informations{$entry}->{'genome_part'}) - (length $lastblock_sequence))){
											my $midblock_sequence = substr($model_informations{$entry}->{'genome_part'}, ($line_array[15] - $three_nt -3), ($line_array[16] + $five_nt));
											$model_informations{$entry}->{'genome_part'} = ${firstblock_sequence}.${midblock_sequence}.${lastblock_sequence};
											#print Dumper (${firstblock_sequence},${midblock_sequence},${lastblock_sequence},$model_informations{$entry}->{'genome_part'});
										}
										else{
											$lastblock_sequence = substr($model_informations{$entry}->{'genome_part'}, ($line_array[15] - $three_nt -3),(length $model_informations{$entry}->{'genome_part'}));
											$model_informations{$entry}->{'genome_part'} = ${firstblock_sequence}.${lastblock_sequence};
										}
									}
									elsif ($three_splicesite == 0 && $five_splicesite == 1){
										$firstblock_sequence = substr($model_informations{$entry}->{'genome_part'}, 0,($line_array[16] + $five_nt));
										$model_informations{$entry}->{'genome_part'} = ${firstblock_sequence}.${lastblock_sequence};
									}
									elsif ($three_splicesite == 1 && $five_splicesite == 0){
										$lastblock_sequence = substr($model_informations{$entry}->{'genome_part'}, ($line_array[15] - $three_nt -3),(length $model_informations{$entry}->{'genome_part'}));
										$model_informations{$entry}->{'genome_part'} = ${firstblock_sequence}.${lastblock_sequence};
									}
									# else{
									# }
								}
								else{
									if ($line_array[17] == 2){
										my $newtarget = $model_informations{$entry}->{'genome_part'};
										my @line_array = split/\t/, $line;
										my $tstart = $line_array[15];
										my $tend = $line_array[16];
										my @blocks = split /,/, $line_array[18];
										my @qstarts = split /,/, $line_array[19];
										my $block_1 = $blocks[0];
										my $block_2 = $blocks[1];
										my $delete_aa_end_1 = ($tstart + ($block_1 * 3));
										my $delete_aa_start_2 = ($tend - ($block_2 * 3));	
										my $gap = $delete_aa_start_2 - $delete_aa_end_1;
										#print Dumper ($gap);
										if ($gap <=3){
											#print Dumper("1");
											substr($newtarget, ($delete_aa_end_1 +1), $gap) = "";
											$model_informations{$entry}->{'genome_part'} = ${firstblock_sequence}.${newtarget}.${lastblock_sequence};
										}
										else{#print Dumper("2");
											my $testsequence = substr($newtarget, ($delete_aa_end_1 +1), $gap);
											#print Dumper ($testsequence);
											my @testsequence_array = split //, $testsequence;
											my ($three_splicesite, $three_nt) = three_splice_site_intron(@testsequence_array);
											my ($five_splicesite, $five_nt) = five_splice_site_canonical_intron(@testsequence_array);
											#print Dumper ($three_splicesite, $five_splicesite);

											if ($three_splicesite == 1 && $five_splicesite == 1){
												substr($newtarget, ($delete_aa_end_1 +1 + $five_nt), $gap -$three_nt) = "";
												my $midblock_sequence = substr($newtarget, ($tstart - $three_nt), ($tend + $five_nt));
											}
											else{
												my $midblock_testsequence = substr($newtarget, ($line_array[15] -  $extension), ($line_array[16] +  $extension));
												#print Dumper ($midblock_testsequence);
												my @midblock_testsequence = split //, $midblock_testsequence;
												my ($five_splicesite, $five_nt) = five_splice_site_canonical(\@midblock_testsequence, 0, 0,$extension);
												my ($three_splicesite, $three_nt) = three_splice_site(\@midblock_testsequence, 0, 0,$extension);
												#print Dumper ($three_splicesite, $five_splicesite, $three_nt,$five_nt);
												my $midblock_sequence = substr($newtarget, ($line_array[15] - $three_nt - 3), ($line_array[16] + $five_nt));
												#print Dumper ($firstblock_sequence,$midblock_sequence,$lastblock_sequence);
												$model_informations{$entry}->{'genome_part'} = ${firstblock_sequence}.${midblock_sequence}.${lastblock_sequence};
											}
										}
										#print Dumper ($model_informations{$entry}->{'genome_part'});
									}
								}
							}
						}
						else{
							next;
						}
					}
					else{
						next;
					}
				}
			}

			my ($new_three_nt_2, $new_five_nt_2);
			#print Dumper($model_informations{$entry}->{'score_5splicesite'},$model_informations{$entry}->{'score_3splicesite'});
			if ($model_informations{$entry}->{'exonID'} == 0){
				$new_three_nt_2 = $model_informations{$entry}->{'start_nt'} + 3;
				$new_five_nt_2 = $model_informations{$entry}->{'five_nt_extension'};
			}

			elsif($model_informations{$entry}->{'exonID'} == $lastexon){
				$new_three_nt_2 = $model_informations{$entry}->{'three_nt_extension'};
				$new_five_nt_2 = $model_informations{$entry}->{'stop_nt'};
			}

			else{
				$new_three_nt_2 = $model_informations{$entry}->{'three_nt_extension'};
				$new_five_nt_2 = $model_informations{$entry}->{'five_nt_extension'};
			}
			#print Dumper ($new_three_nt, $new_five_nt);
			my @exonsequence = sequence_final_exon($model_informations{$entry}->{'genome_part'}, $new_three_nt_2, $new_five_nt_2, $model_informations{$entry}->{'extension_start'}, $model_informations{$entry}->{'extension_end'});
			my $dna_exonsequence= join ("",@exonsequence);
			print Dumper ($dna_exonsequence);
			my @triplets = unpack '(a3)*', $dna_exonsequence;
			my @protein_aa_exon;
			for my $codon (@triplets){
				my $aa_exon = $myCodonTable->translate($codon);
		    	push (@protein_aa_exon, $aa_exon);
			}
			my $protein_exon = join ("",@protein_aa_exon);
			print Dumper($protein_exon);
			open (my $OUTTARGET, ">> SearchTarget/${paralog}_merged_target_exons");
			print $OUTTARGET ">${paralog}_$model_informations{$entry}->{'exonID'}_merged_final\n$protein_exon\n";

			#print Dumper($model_informations{$entry}->{'exonID'},$model_informations{$entry}->{'genome_part'}, $new_three_nt, $new_five_nt);
			my @exonsequence_2 = sequence_final_exon($model_informations{$entry}->{'genome_part'}, $new_three_nt, $new_five_nt, $model_informations{$entry}->{'extension_start'}, $model_informations{$entry}->{'extension_end'});
			push (@targetsequence, @exonsequence_2);
		}
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

if(length $protein > 1){	
	open (my $OUTTARGET, "> SearchTarget/${paralog}_target.fa");	
	print $OUTTARGET ">${paralog}_merged_final\n$protein\n";
}
print Dumper ($protein);
#------------------------------------------------------------------------
# SUBROUTINES
#

#find start, end, paralog name from blastout
sub search_data{
	my ($blastoutIn_line) = @_;	
	my @line_array = split "\t", $blastoutIn_line;

	my $match = $line_array[0];
	my $mismatch= $line_array[1];	
	my $gap = $line_array[6];	
	my $strand = $line_array[8];
	my $model = $line_array[9];
	my $targetname = $line_array[13];
	my $start_query = $line_array[11];
	my $end_query = $line_array[12];
	my $blocks = $line_array[17];
	my $blocksizes = $line_array[18];
	my $tstarts = $line_array[20];
	my $qstarts = $line_array[19];
	my $start = $line_array[15];
	$strand =~ s/^.//;
	my $end = $line_array[16];
	my $query_size = $line_array[10];

	my @modelsplit = split "_",$model;
	my $exonID = $modelsplit[1];
	$exonID =~ s/exon//;
	my $paralogID = $modelsplit[0];

	my %line_record = (		
		start_query => $start_query,
    	end_query => $end_query,
    	start => $start,
    	end => $end,
    	strand => $strand,
    	model => $model,
    	paralogID => $paralogID,
    	targetname => $targetname,
		match => $match,
    	exonID => $exonID, 	
    	mismatch => $mismatch,
    	gap => $gap,	
    	blocks => $blocks,
		blocksizes => $blocksizes,
		tstarts => $tstarts,
		qstarts => $qstarts,
		query_size => $query_size,
	);
	return %line_record;
}

sub search_data_blast{
	my ($blastoutIn_line) = @_;

	my @line_array = split "\t", $blastoutIn_line;
	my $model = $line_array[0];
	my $targetname = $line_array[1];
	my $bias = $line_array[2];	
	my $mismatch = $line_array[4];
	my $gap = $line_array[5];	
	my $start_query = ($line_array[6] -1);
	my $end_query = $line_array[7];
	my $start = $line_array[8];
	my $end = $line_array[9];
	my $evalue = $line_array[10];
	my $bitscore = $line_array[11];
	my $strand;	
	my $blocks = 1;
	my $blocksizes = 0;
	my $tstarts = 0;
	my $qstarts = 0;	
	my $query_size = $line_array[3];

	if ($line_array[8]<$line_array[9]){
		$strand = "+";
	}
	elsif ($line_array[8] > $line_array[9]){
		$strand = "-";
	};

	my @modelsplit = split "_",$model;
	my $exonID = $modelsplit[1];
	$exonID =~ s/exon//;
	my $paralogID = $modelsplit[0];

	my %line_record = (		
		start_query => $start_query,
    	end_query => $end_query,
    	start => $start,
    	end => $end,
    	strand => $strand,
    	model => $model,
    	paralogID => $paralogID,
    	targetname => $targetname,
    	evalue => $evalue,
    	bitscore => $bitscore,
    	bias => $bias,
    	exonID => $exonID, 	
    	mismatch => $mismatch,
    	gap => $gap,    	
    	blocks => $blocks,
		blocksizes => $blocksizes,
		tstarts => $tstarts,
		qstarts => $qstarts,
		query_size => $query_size,

	);
	return %line_record;
}

sub find_best_exon_mapping{
    my ($exons_ref, $model_informations_ref) = @_;
    my %number; 
    my $remember;
    for my $elem (@$exons_ref){
    	for my $entry (sort {$a cmp $b} keys %$model_informations_ref){
    		if ($$model_informations_ref{$entry}->{'exonID'} == $elem){
    			$number{$elem}++;
    			if($number{$elem} == 2){
    				if ($$model_informations_ref{$entry}->{'bitscore'} <= $$model_informations_ref{$remember}->{'bitscore'}){
    					delete $$model_informations_ref{$entry};
    				}
    				elsif ($$model_informations_ref{$remember}->{'bitscore'} <= $$model_informations_ref{$entry}->{'bitscore'}){
    					delete ($$model_informations_ref{$remember});
    				}
    			}
    			else{
    				$remember = $entry; 
    			}
    		}
    	}
    }
    return (%$model_informations_ref);
}


sub best_three_splicesite{	
	my ($model_informations) = @_ ;
	#print Dumper ($model_informations);
	my @a_short = split //, $model_informations->{genome_part};

	if ($model_informations->{three_splicesite} == 0){
		#print Dumper("1");
		my ($three_splicesite, $three_nt) = three_splice_site(\@a_short, $model_informations->{mismatch}, $model_informations->{start_query}, $model_informations->{extension_start});

		$model_informations->{three_nt} = $three_nt;
		$model_informations->{three_splicesite} = $three_splicesite;
		
		my %three_scores;			

		if ($three_splicesite == 1){
			my @testsequence_three_splicesite = sequence_three_splicesite(\@a_short, $three_nt, $model_informations->{extension_start});	

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
			#print Dumper("2");
			$model_informations->{three_nt} = 0;
			my ($three_splicesite, $three_nt) = three_splice_site_minus_three(\@a_short, $model_informations->{extension_start});
			#print Dumper ($three_splicesite, $three_nt);
			$model_informations->{three_nt} = $three_nt;
			$model_informations->{three_splicesite} = $three_splicesite;

			if ($three_splicesite == 1){
				my @testsequence_three_splicesite = sequence_three_splicesite(\@a_short, $three_nt, $model_informations->{extension_start});	

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
			#print Dumper("3");
			$model_informations->{score_3splicesite} = -333;
			until ($model_informations->{score_3splicesite} > 0 || $model_informations->{three_nt} >= $model_informations->{extension_start}-18){
				my @a_short = split //, $model_informations->{genome_part};
				my ($three_splicesite, $three_nt) = three_splice_site_after_negativ_score(\@a_short, $model_informations->{three_nt}, $model_informations->{extension_start});
				$model_informations->{three_nt} = $three_nt;
				$model_informations->{three_splicesite} = $three_splicesite;
				if ($three_splicesite == 1){
					my @testsequence_three_splicesite = sequence_three_splicesite(\@a_short, $model_informations->{three_nt}, $model_informations->{extension_start});
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
					#print Dumper("muh");
				 	$model_informations->{score_3splicesite} = -333;
				}
			}
		}
		if( !%three_scores){
			$model_informations->{number_three_scores} = 0;
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
	}    
	return ($model_informations);
}

sub best_five_splicesite{    
	my ($model_informations) = @_ ;
	#print Dumper ($model_informations);
	my @a_short = split //, $model_informations->{genome_part};
	# print Dumper (@a_short);
	$model_informations->{stopcodon} = 0;
	$model_informations->{startcodon} = 0;
	
	if ($model_informations->{five_splicesite} == 0){		
		my %five_scores;	
		my ($five_splicesite, $five_nt) = five_splice_site_canonical(\@a_short, $model_informations->{mismatch}, $model_informations->{gap}, $model_informations->{extension_end});
		#print Dumper($five_splicesite);
		$model_informations->{five_nt} = $five_nt;
		$model_informations->{five_splicesite} = $five_splicesite;
		#print Dumper($model_informations->{five_nt} );
		if ($five_splicesite == 1){
			#print Dumper ("1");
			my @testsequence_five_splicesite = sequence_five_splicesite(\@a_short, $five_nt, $model_informations->{extension_end}); 
			if (defined $testsequence_five_splicesite[8]){
				$model_informations->{testsequence_five_splicesite} = join ("",@testsequence_five_splicesite);
				open (testsequences_five_splicesite_file, '>>', "MaxEntScan/Testsequences_five_splicesite.txt") or die "$!";
				print testsequences_five_splicesite_file join ("",@testsequence_five_splicesite);
			}
			else {
				print "unvollstaendig 5: $model_informations->{model}\n";
				$model_informations->{five_splicesite} = 0;
				$model_informations = new_splicesite_after_neighbouring_mismatch_three($model_informations);
				$model_informations = new_splicesite_after_neighbouring_mismatch_five($model_informations);
				#print Dumper ($model_informations);
				return ($model_informations);
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
		if ($five_splicesite == 0 || $model_informations->{score_5splicesite} < 0 ){
			#print Dumper ($model_informations->{five_nt});
			$model_informations->{score_5splicesite} = -555;
			until ($model_informations->{score_5splicesite} > 0 || $model_informations->{five_nt} >= $model_informations->{extension_end} -4){			
				#print Dumper ("3",$model_informations->{five_nt} );
				my ($five_splicesite, $five_nt) = five_splice_site_canonical_after_negative_score(\@a_short, $model_informations->{five_nt}, $model_informations->{extension_end});
				$model_informations->{five_nt} = $five_nt;
				$model_informations->{five_splicesite} = $five_splicesite;
				if ($five_splicesite == 1){
					my @testsequence_five_splicesite = sequence_five_splicesite(\@a_short, $model_informations->{five_nt}, $model_informations->{extension_end});
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
			until ($model_informations->{score_5splicesite} > 0 || $model_informations->{five_nt} >= $model_informations->{extension_end} -4){
				#print Dumper ($model_informations->{five_nt});
				my ($five_splicesite, $five_nt) = five_splice_site_noncanonical(\@a_short,$model_informations->{five_nt}, $model_informations->{extension_end});
				$model_informations->{five_nt} = $five_nt;
				$model_informations->{five_splicesite} = $five_splicesite;
				if ($five_splicesite == 1){
					my @testsequence_five_splicesite = sequence_five_splicesite(\@a_short, $model_informations->{five_nt}, $model_informations->{extension_end});
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
				else{#print Dumper ("ups");
				 	$model_informations->{score_5splicesite} = -555;
				}
			}
		}
		if( !%five_scores){
			$model_informations->{number_five_scores} = 0;
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
		#print Dumper ($model_informations);
	}		
	#print Dumper ("ups",$model_informations);

    return ($model_informations);	
}

sub find_stopcodon {
    my ($genome_target_ref, $extension) = @_;  
	my @genome_target = @$genome_target_ref; 
    my $stopcodon = 0; #1 if stopcodon is found, else 0
    my $stop_nt = 0; #nt between stopcodon and query exon

	for(my $i = ((scalar @genome_target) - $extension -1);$i <= (scalar @genome_target) - 3;$i++){	
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
    my ($genome_target_ref, $stop_nt, $extension) = @_; 
    my @genome_target = @$genome_target_ref; 
	my $stopcodon = 0;

	for(my $i = ((scalar @genome_target) - $extension + $stop_nt);$i <= (scalar @genome_target) - 3;$i++){	
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
    my ($genome_target_ref, $qstarts, $extension) = @_;   
	my @genome_target = @$genome_target_ref;
	my $qstart = (split ",", $qstarts)[0];
#print Dumper ($qstart);
    my $startcodon = 0; #1 if startcodon is found, else 0
    my $start_nt = -3 + ($qstart *3);; #nt between startcodon and query exon
	#print Dumper ($start_nt);
	for(my $i = $extension + 2 - ($qstart*3);$i >= 1;$i--){
		#print Dumper ($genome_target[$i]);
        if ($genome_target[$i] eq "G" and $genome_target[$i-1] eq "T" and $genome_target[$i-2] eq "A"){        
        	$startcodon = 1;    
        	return ($startcodon, $start_nt);
        } 
        else {
        	$startcodon = 0;
        	$start_nt++;
        	#print Dumper ($start_nt);
        }
    }
    return ($startcodon, $start_nt);
}

sub find_startcodon_after_modulo {
    my ($genome_target_ref, $start_nt, $extension) = @_;   
	my @genome_target = @$genome_target_ref;
    my $startcodon = 0; #1 if startcodon is found, else 0

#print Dumper ($start_nt);
	for(my $i = $extension - $start_nt -2;$i >= 2;$i--){
		#print Dumper ($genome_target[$i]);
        if ($genome_target[$i] eq "G" and $genome_target[$i-1] eq "T" and $genome_target[$i-2] eq "A"){        
        	$startcodon = 1;    #print Dumper ($start_nt);
        	return ($startcodon, $start_nt);
        	#$start_nt++;
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
    my ($genome_target_ref, $mismatch, $qstarts, $extension) = @_;   
	my @genome_target = @$genome_target_ref;
	my $qstart = (split ",", $qstarts)[0];
	my $three_nt;
	#print Dumper ($qstart);
    my $three_splicesite = 0; #1 if ss is found, else 0
    if ($mismatch < 5 || $mismatch > 15){
	$three_nt = 0 + ($qstart *3);
		for(my $i = $extension - 1 -($qstart*3) ;$i >= 16;$i--){			
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
		for(my $i = $extension - 1 -($qstart*3) + ($mismatch *3); $i >= 16;$i--){			
			#print Dumper ($genome_target[$i], $i, $three_nt);
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
    my ($genome_target_ref, $extension) = @_;  
    my @genome_target = @$genome_target_ref; 
    my $three_splicesite = 0; #1 if ss is found, else 0
 	my $three_nt = - 9;
	for(my $i = $extension +8 ;$i >= 16;$i--){
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
    my ($genome_target_ref, $three_nt, $extension) = @_;   
    my @genome_target = @$genome_target_ref;
    my $three_splicesite = 0; #1 if ss is found, else 0

	for(my $i = ($extension -3 - $three_nt);$i >= 16;$i--){
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
    my ($genome_target_ref, $three_nt) = @_;   
    my @genome_target = @$genome_target_ref;
    my $three_splicesite = 0; #1 if ss is found, else 0
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
    my ($genome_target_ref, $mismatch, $gap, $extension) = @_;  
    #print Dumper ($genome_target_ref);
    my @genome_target = @$genome_target_ref;
    my $five_splicesite; #1 if splicesite is found, else 0
    my $five_nt = -3; #nt between splice site and query exon

	for(my $i = ((scalar @genome_target) - $extension -4 );$i <= (scalar @genome_target) - 2;$i++){	
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
sub five_splice_site_canonical_minus_three {
    my ($genome_target_ref, $mismatch, $gap, $extension) = @_;  
    #print Dumper ($genome_target_ref);
    my @genome_target = @$genome_target_ref;
    my $five_splicesite; #1 if splicesite is found, else 0
    my $five_nt = -9; #nt between splice site and query exon

	for(my $i = ((scalar @genome_target) - $extension -9 );$i <= (scalar @genome_target) - 2;$i++){	
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
    my ($genome_target_ref, $five_nt, $extension) = @_;  
    my @genome_target = @$genome_target_ref;
    my $five_splicesite = 0; #1 if splicesite is found, else 0

	for(my $i = ((scalar @genome_target) - $extension +1 + $five_nt);$i <= (scalar @genome_target) - 2;$i++){	
		#print Dumper ($five_nt);
		#print Dumper ($i, $genome_target[$i]);
        if ($genome_target[$i] eq "G" and $genome_target[$i+1] eq "T"){     
        	#print Dumper ("yes");
        	$five_splicesite = 1;    
        	return ($five_splicesite, ($five_nt + 2));
        } 
        else {
        	#print Dumper ("no");
        	$five_splicesite = 0;
        	$five_nt++;
        }
    }
    return ($five_splicesite, $five_nt);
}

sub five_splice_site_noncanonical {   
    my ($genome_target_ref, $five_nt, $extension) = @_;  
    my @genome_target = @$genome_target_ref;
    my $five_splicesite; #1 if splicesite is found, else 0
    #my $five_nt = -3; #nt between splice site and query exon

	for(my $i = ((scalar @genome_target) - $extension +3 + $five_nt);$i <= (scalar @genome_target) - 2;$i++){
	#for(my $i = ((scalar @genome_target) - $extension -3);$i <= (scalar @genome_target) - 2;$i++){	
        if ($genome_target[$i] eq "G" and $genome_target[$i+1] eq "C"){     

        	$five_splicesite = 1;    
        	return ($five_splicesite, ($five_nt + 2));
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
	my ($genome_target_ref, $three_nt, $extension) = @_;      
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
	my ($genome_target_ref, $five_nt, $extension) = @_;      
	my @genome_target = @$genome_target_ref;
	my @testsequence_five_splicesite;	
	# print Dumper (scalar @genome_target);
	# print Dumper ($five_nt);
	for(my $i = (scalar @genome_target - $extension -4 + $five_nt);$i <= (scalar @genome_target - $extension + 4 + $five_nt);$i++) {
		# print Dumper ($i);
		# print Dumper ($genome_target[$i]);
		push (@testsequence_five_splicesite, $genome_target[$i]) if defined $genome_target[$i];
		}
	return (@testsequence_five_splicesite);
}

sub splitted_as{
	my ($model_informations_ref) = @_ ;
	my %model_informations = %$model_informations_ref;

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
				else{
					$model_informations{$i}->{'splitted_for_header3'} = "3: 0";
				}
			#	-> durch 3 teilbar, query verlaengern
				if ($model_informations{$i}->{'strand'} eq "+"){
					$model_informations{$i}->{'three_nt_extension'} = POSIX::floor($model_informations{$i}->{'three_nt'} / 3) * 3;
					my $extension = (int $model_informations{$i}->{'three_nt'} / 3) * 3;
					#$model_informations{$i}->{'start'} = $model_informations{$i}->{'start'} - $model_informations{$i}->{'extension_start'};
				}
				elsif($model_informations{$i}->{'strand'} eq "-"){
					if ($model_informations{$i}->{'blat'} eq "yes"){
						$model_informations{$i}->{'three_nt_extension'} = POSIX::floor($model_informations{$i}->{'three_nt'} / 3) * 3;
						my $extension = (int $model_informations{$i}->{'three_nt'} / 3) * 3;
						#$model_informations{$i}->{'end'} = $model_informations{$i}->{'end'} + $model_informations{$i}->{'extension_end'};
					}
					else{
						$model_informations{$i}->{'three_nt_extension'} = POSIX::floor($model_informations{$i}->{'three_nt'} / 3) * 3;
						my $extension = (int $model_informations{$i}->{'three_nt'} / 3) * 3;
						#$model_informations{$i}->{'start'} = $model_informations{$i}->{'start'} + $model_informations{$i}->{'extension_start'};
					}
				}
			}
			else{
				$model_informations{$i}->{'splitted_AS3'} = 0;
				$model_informations{$i}->{'three_nt_extension'} = 0;			
			}
		}
		#print Dumper $model_informations{$i}->{'exonID'};
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
					#$model_informations{$i}->{'end'} = $model_informations{$i}->{'end'} + $model_informations{$i}->{'extension_end'} +1;
				}
				elsif($model_informations{$i}->{'strand'} eq "-"){
					if ($model_informations{$i}->{'blat'} eq "yes"){
						$model_informations{$i}->{'five_nt_extension'} = POSIX::floor($model_informations{$i}->{'five_nt'} / 3) * 3;
						my $extension = (int $model_informations{$i}->{'five_nt'} / 3) * 3;
						#$model_informations{$i}->{'start'} = $model_informations{$i}->{'start'} - $model_informations{$i}->{'extension_start'} +1;
					}
					else{
						$model_informations{$i}->{'five_nt_extension'} = POSIX::floor($model_informations{$i}->{'five_nt'} / 3) * 3;
						my $extension = (int $model_informations{$i}->{'five_nt'} / 3) * 3;
						#$model_informations{$i}->{'end'} = $model_informations{$i}->{'end'} - $model_informations{$i}->{'extension_end'} +1;
					}
				}
			}
			else{
				$model_informations{$i}->{'splitted_AS5'} = 0;
				$model_informations{$i}->{'five_nt_extension'} = 0;
			}
		}
	}
	return (%model_informations);
}


sub matching_neighbours{
	my ($model_informations_ref) = @_ ;
	my %model_informations = %$model_informations_ref;
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
				$model_informations{$i} = new_splicesite_after_neighbouring_mismatch_three($model_informations{$i});
				$model_informations{$model_compare} = new_splicesite_after_neighbouring_mismatch_five($model_informations{$model_compare});

				#matching_neighbours
				$ref_five_scores = $model_informations{$model_compare}->{'five_scores'};
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
				
				if( !%possible_matches){

					#looking, if neighbouring exons are matching without intron
					#print Dumper ("voegh nevuhrmtgchimrui");
					my $border_end;
					my $border_start;
					my $nt_start;
					my $nt_end;

					if ($model_informations{$i}->{'blat'} eq "yes" && $model_informations{$model_compare}->{'blat'} eq "yes"){
						print Dumper ("wichtig");
					}
					elsif($model_informations{$i}->{'blat'} eq "yes" && $model_informations{$model_compare}->{'blat'} eq "no"){										
						if ($model_informations{$i}->{'strand'} eq "-"){
							$border_end = $model_informations{$i}->{'end'};
							$border_start = $model_informations{$model_compare}->{'end'} - $model_informations{$i}->{'extension_end'} + $model_informations{$i}->{'end'};
							$nt_start = ((length $model_informations{$model_compare}->{'query'}) - $model_informations{$model_compare}->{'query_size'})*3;
							$nt_end = $model_informations{$i}->{'start_query'} *3;				
							if ($border_end - $border_start - $nt_start - $nt_end < 21)	{
								#print Dumper ($border_end, $border_start, $nt_start, $nt_end);

								$model_informations{$i}->{'three_splicesite'} = 1;
								$model_informations{$i}->{'three_nt'} = ($border_start - $border_end +2);
								$model_informations{$i}->{'score_3splicesite'} = 100;
								$model_informations{$i}->{'five_splicesite'} = 0;
								$model_informations{$model_compare}->{'five_splicesite'} = 1;
								$model_informations{$model_compare}->{'five_nt'} = 0;
								$model_informations{$model_compare}->{'score_5splicesite'} = 100;
								$model_informations{$model_compare}->{'three_splicesite'} = 0;
							}
						}
						else{
							print Dumper ("wichtig");
						}

					}
					elsif($model_informations{$i}->{'blat'} eq "no" && $model_informations{$model_compare}->{'blat'} eq "yes"){	
						print Dumper ("wichtig");	

					}
					else{
						print Dumper ("wichtig");
					}
					next;
				}
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
	return (%model_informations);
}

sub new_splicesite_after_neighbouring_mismatch_three{
	my ($model_informations) = @_ ;
	#print Dumper("huhuhuhuhuhuhuhuhuhuhuhuu");
	
	#$model_informations->{start_query} = $model_informations->{start_query} +5;

	#vergrößerung extension	
	my $extension_2 = 53;	
	my @a_short;	
	if($model_informations->{blat} eq "yes"){
		my $genome_seqIO = Bio::SeqIO->new(-file   => $target, -format => "Fasta");
		$a = $genome_seqIO->next_seq;
		$a = $a -> seq;	
	}
	else{
		my $location = (split ":",$model_informations->{targetname})[-1];
		my $start = (split "-", $location)[0];
		my $end = (split "-", $location)[1];

		my $genome_seqIO = Bio::SeqIO->new(-file   => "TargetSequences/${paralog}_${start}_${end}.fa", -format => "Fasta");
		$a = $genome_seqIO->next_seq;
		$a = $a -> seq;
	}	

	my $extended_start;
	my $extended_end;

	if($model_informations->{blat} eq "no"){

		if ($model_informations->{strand} eq "+"){
			if (($model_informations->{start}) - 1 - $extension_2 < 0){
				#print Dumper ("1");
			$model_informations->{extension_start} = ($model_informations->{start} - $extension_2);
			}
			else{
				#print Dumper ("2");
				$model_informations->{extension_start} = $extension_2;
			}
			
			$extended_start = $model_informations->{start} -1 -$extension_2;

			if (((length $a) - $model_informations->{end}) -1 < $extension_2){
				#print Dumper ("3");
				$model_informations->{extension_end} = ((length $a) - $model_informations->{end} -1);
			}
			else{
				#print Dumper ("4");
				$model_informations->{extension_end} = $extension_2 -1;
			}
			$extended_end = $model_informations->{end} + $model_informations->{extension_end} ;


			if ($model_informations->{end_query} < ($model_informations->{query_size} - 10)){
				#print Dumper ("5");
				my $end_extension = (($model_informations->{query_size} - $model_informations->{end_query}) *3);	
				$extended_end += $end_extension;
			}
		}
		else{
			if (($model_informations->{end}) - 1 - $extension_2 < 0){
				#print Dumper ("1");
			$model_informations->{extension_end} = ($model_informations->{end} - $extension_2);
			}
			else{
				#print Dumper ("2");
				$model_informations->{extension_end} = $extension_2 -1;
			}
			
			$extended_start = $model_informations->{end} +1 -$extension_2;

			if (((length $a) - $model_informations->{start}) -1 < $extension_2){
				#print Dumper ("3");
				$model_informations->{extension_start} = ((length $a) - $model_informations->{start} -1);
			}
			else{
				#print Dumper ("4");
				$model_informations->{extension_start} = $extension_2;
			}
			$extended_end = $model_informations->{start} -1 + $model_informations->{extension_start};

			if ($model_informations->{end_query} < ($model_informations->{query_size} - 10)){
				#print Dumper ("5");
				my $end_extension = (($model_informations->{query_size} - $model_informations->{end_query}) *3);
				$extended_start -= $end_extension;
			}
		}

		if ($model_informations->{strand} eq "+"){		
			@a = split //, $a;
			for(my $i = $extended_start;$i <= $extended_end;$i++) {
				push (@a_short, $a[$i]);
			}
		}
		elsif($model_informations->{strand} eq "-"){
			#print Dumper ("6", $a);
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
	}
	else{
		if (($model_informations->{start}) - $extension_2 < 0){
				$model_informations->{extension_start} = ($model_informations->{start} - $extension_2);
			}
		else{
			$model_informations->{extension_start} = $extension_2;
		}
		
		my $extended_start = $model_informations->{start} - $model_informations->{extension_start};

		if (((length $a) - $model_informations->{end}) < $extension_2){
			#print Dumper("mäp");
			$model_informations->{extension_end} = ((length $a) - $model_informations->{end} -1);
		}
		else{
			$model_informations->{extension_end} = $extension_2 -1;
		}
		my $extended_end = $model_informations->{end} + $model_informations->{extension_end};


		if ($model_informations->{end_query} < ($model_informations->{query_size} - 10)){
			my $end_extension = (($model_informations->{query_size} - $model_informations->{end_query}) *3);
				#print Dumper("kleiner");
			if ($model_informations->{strand} eq "+" ){
				$extended_end += $end_extension;
			}
			else {
				$extended_start -= $end_extension;
			}
		}
		#print Dumper ($model_informations->{start},$model_informations->{end}, $extended_start, $extended_end);
		if ($model_informations->{strand} eq "+"){		
			@a = split //, $a;
			for(my $i = $extended_start;$i <= $extended_end;$i++) {
				push (@a_short, $a[$i]);
			}
		}
		elsif($model_informations->{strand} eq "-"){
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
	}

	$model_informations->{genome_part} = join ("",@a_short);
#print Dumper ($model_informations->{genome_part});
	#print Dumper("1");
	my ($three_splicesite, $three_nt);
	my %three_scores;

	($three_splicesite, $three_nt) = three_splice_site(\@a_short, $model_informations->{mismatch}, $model_informations->{start_query}, $model_informations->{extension_start});

	$model_informations->{three_nt} = $three_nt;
	$model_informations->{three_splicesite} = $three_splicesite;

	if ($three_splicesite == 1){
		my @testsequence_three_splicesite = sequence_three_splicesite(\@a_short, $three_nt, $model_informations->{extension_start});	

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
		#print Dumper("2");
		$model_informations->{three_nt} = 0;
		my ($three_splicesite, $three_nt) = three_splice_site_minus_three(\@a_short, $model_informations->{extension_start});

		$model_informations->{three_nt} = $three_nt;
		$model_informations->{three_splicesite} = $three_splicesite;

		if ($three_splicesite == 1){
			my @testsequence_three_splicesite = sequence_three_splicesite(\@a_short, $three_nt, $model_informations->{extension_start});	

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

	until ($model_informations->{extension_start} - $model_informations->{three_nt} <= 25){		
		#print Dumper ($model_informations->{extension_start} - $model_informations->{three_nt});

		#print Dumper("3");
		my @a_short = split //, $model_informations->{genome_part};
		my ($three_splicesite, $three_nt) = three_splice_site_after_negativ_score(\@a_short, $model_informations->{three_nt}, $model_informations->{extension_start});
		$model_informations->{three_nt} = $three_nt;
		$model_informations->{three_splicesite} = $three_splicesite;
		if ($three_splicesite == 1){
			my @testsequence_three_splicesite = sequence_three_splicesite(\@a_short, $model_informations->{three_nt}, $model_informations->{extension_start});
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
			#print Dumper("muh");
		 	$model_informations->{score_3splicesite} = -333;
		}
		$model_informations->{three_nt} = $model_informations->{three_nt} + 2;
	}

	if( !%three_scores){
		$model_informations->{number_three_scores} = 0;
		return ($model_informations);
	}
	#print Dumper(%three_scores);
	my $best_three_score = (sort { $a <=> $b } (keys %three_scores))[-1];
	#print Dumper ($best_three_score);
	$model_informations->{number_three_scores} = %three_scores;
	$model_informations->{three_scores} = \%three_scores;
	$model_informations->{score_3splicesite} = $best_three_score;
	$model_informations->{three_nt} = $three_scores{$best_three_score};
	#print Dumper ($model_informations);
	return ($model_informations);
}

sub new_splicesite_after_neighbouring_mismatch_five{	
	my ($model_informations) = @_ ;
	my @a_short = split //, $model_informations->{genome_part};
	$model_informations->{stopcodon} = 0;
	$model_informations->{startcodon} = 0;
	
	if ($model_informations->{five_splicesite} == 0){		
		my %five_scores;	
		my ($five_splicesite, $five_nt) = five_splice_site_canonical(\@a_short, $model_informations->{mismatch}, $model_informations->{gap}, $model_informations->{extension_end});
		$model_informations->{five_nt} = $five_nt;
		$model_informations->{five_splicesite} = $five_splicesite;

		if ($five_splicesite == 1){
			my @testsequence_five_splicesite = sequence_five_splicesite(\@a_short, $five_nt, $model_informations->{extension_end}); 
			if (defined $testsequence_five_splicesite[8]){
				$model_informations->{testsequence_five_splicesite} = join ("",@testsequence_five_splicesite);
				open (testsequences_five_splicesite_file, '>>', "MaxEntScan/Testsequences_five_splicesite.txt") or die "$!";
				print testsequences_five_splicesite_file join ("",@testsequence_five_splicesite);
			}
			else {
				print "unvollstaendig 5: %model_informations{model}\n";
				$model_informations->{five_splicesite} = 5;
			}
			if(-e "MaxEntScan/Testsequences_five_splicesite.txt"){
				my $output_5splicesite = qx(perl $path_to_maxentscan_5 "MaxEntScan/Testsequences_five_splicesite.txt");
				my @array_5splicesite = split /\t/, $output_5splicesite;

				if (defined $array_5splicesite[1]){
					chomp ($array_5splicesite[1]);
				}
				$model_informations->{score_5splicesite} = $array_5splicesite[1];							
				my $key_for_hash = $model_informations->{score_5splicesite};
				$five_scores{$key_for_hash} = $model_informations->{five_nt};
				unlink "MaxEntScan/Testsequences_five_splicesite.txt";
			}		
		}
		if ($five_splicesite == 0 || $model_informations->{score_5splicesite} < 0 ){
			$model_informations->{score_5splicesite} = -555;
			$model_informations->{five_nt} = -15;
			until ($model_informations->{score_5splicesite} > 0 || $model_informations->{five_nt} > $model_informations->{extension_end} -5){			
				my ($five_splicesite, $five_nt) = five_splice_site_canonical_after_negative_score(\@a_short, $model_informations->{five_nt}, $model_informations->{extension_end});
				$model_informations->{five_nt} = $five_nt;
				$model_informations->{five_splicesite} = $five_splicesite;
				if ($five_splicesite == 1){
					my @testsequence_five_splicesite = sequence_five_splicesite(\@a_short, $model_informations->{five_nt}, $model_informations->{extension_end});
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
			}
		}
		if ($five_splicesite == 0 || $model_informations->{score_5splicesite} < 0 ){
			$model_informations->{five_nt} = -3;
			until ($model_informations->{score_5splicesite} > 0 || $model_informations->{five_nt} >= $model_informations->{extension_end} -4){
				my ($five_splicesite, $five_nt) = five_splice_site_noncanonical(\@a_short,$model_informations->{five_nt}, $model_informations->{extension_end});
				$model_informations->{five_nt} = $five_nt;
				$model_informations->{five_splicesite} = $five_splicesite;
				if ($five_splicesite == 1){
					my @testsequence_five_splicesite = sequence_five_splicesite(\@a_short, $model_informations->{five_nt}, $model_informations->{extension_end});
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
		my $best_five_score = (sort { $a <=> $b } (keys %five_scores))[-1];
		$model_informations->{number_five_scores} = %five_scores;
		$model_informations->{five_scores} = \%five_scores;
		$model_informations->{score_5splicesite} = $best_five_score;
		$model_informations->{five_nt} = $five_scores{$best_five_score};
	}

	return ($model_informations);
}

sub sequence_final_exon{
	my ($genome_target, $three_nt, $five_nt,$extension_start, $extension_end) = @_;      
	my @genome_target_array = split //, $genome_target;
#print Dumper ($genome_target, $three_nt, $five_nt,$extension_start, $extension_end);
	my @sequence_final_exon;
	for(my $i = ($extension_start - $three_nt);$i <= (scalar @genome_target_array - $extension_end -2 + $five_nt);$i++) {
		#print Dumper ($i);
		#print Dumper ($genome_target_array[$i]);
		push (@sequence_final_exon, $genome_target_array[$i]) if defined $genome_target_array[$i];
		}

		#print Dumper (@sequence_final_exon);
	return (@sequence_final_exon);
}