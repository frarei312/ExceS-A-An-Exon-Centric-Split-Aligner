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
use Time::Duration;

#------------------------------------------------------------------------
#Definitions
#

my $genomeIn;
my $blastoutIn;
my $mode;
my %model_informations;
my $core = 8;
my $extension = 35;


my $path_to_maxentscan_3 = "/scratch/franziska/bin/maxentscan/fordownload/score3.pl";
my $path_to_maxentscan_5 = "/scratch/franziska/bin/maxentscan/fordownload/score5.pl";

GetOptions("genome=s" => \$genomeIn, "blastout=s" => \$blastoutIn, "mode=s" => \$mode);

#------------------------------------------------------------------------
#Main
#
my $start_time = time();

unless (-d "MaxEntScan/"){
	mkdir "MaxEntScan/";
}

create_TargetSequences($blastoutIn,$genomeIn);
my %exonsequence;

my $seqio = Bio::SeqIO->new(-file => "HomologousExons/paralog_exon_alignment_final.fa", -format => "fasta");
while(my$seqobj = $seqio->next_seq) {
    my $id  = $seqobj->display_id;
    my $seq = $seqobj->seq;
    $exonsequence{$id} = $seq;
}

open (my $blastout_handle, "<", $blastoutIn) or die "Can't open Blastout, $!";
my @a;
my @b;
my $a_align;
my $b_align;

#data from blastout

while (defined (my $line = <$blastout_handle>)){
	chomp $line;	
	#print Dumper ($line);
	next if ($line =~ m/^\#/ || $line =~ m/^$/);  # skip comments and empty lines
	my %data;
	
	my @line_split = split /\t/, $line;	
	if ($line_split[2] > 50.00){
		%data = search_data_blast($line);
		#print Dumper ($start);
	}
	else{
		next;
	}

	my $genome_seqIO = Bio::SeqIO->new(-file   => "TargetSequences/".$data{targetname}.".fa", -format => "Fasta");
	$a = $genome_seqIO->next_seq;
	$a = $a -> seq;

	my @a_short;			

	if (($data{start}) - $extension < 0){
		$data{extension_start} = ($extension - $data{start} +1);
	}
	else{
		$data{extension_start} = $extension +1;
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
	$data{three_splicesite} = 0,
	$data{five_splicesite} = 0,
	$data{blat} = "no";	
	$data{stopcodon} = 0;
	$data{startcodon} = 0;
	my $key_for_hash = $data{model}."_".$data{start};
	$model_informations{$key_for_hash} = \%data;			
		
}
close $blastout_handle;

#print "1\n";
#print "Runtime ", duration(time() - $start_time), ".\n";

my @exons;

for my $i (keys %model_informations) {
	push (@exons, $model_informations{$i}->{'exonID'});
}	
my @sorted_exons = sort {$a <=> $b} @exons;

my $lastexon=$sorted_exons[-1];

for my $i (keys %model_informations) {
	#print Dumper ($model_informations{$i}->{'model'});
	my @a_short = split //, $model_informations{$i}->{'genome_part'};

	#search stopcodon next to last exon
	if ($model_informations{$i}->{'exonID'} == $lastexon){
		#print Dumper ($model_informations{$i});
		$model_informations{$i}->{'stopcodon'} = 0;
		my ($stopcodon, $stop_nt) = find_stopcodon(\@a_short, $model_informations{$i}->{'extension_end'});	
		$model_informations{$i}->{'stop_nt'}  = $stop_nt;		
		$model_informations{$i}->{'modulo_stop'} = $model_informations{$i}->{'stop_nt'} % 3;

		if ($stopcodon == 1 && ($model_informations{$i}->{'modulo_stop'} == 0)){
			#print Dumper("yes");
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
		if ($model_informations{$i}->{'score_3splicesite'} == -333){
			$model_informations{$i} = extend_search_radius($model_informations{$i});
			$model_informations{$i} = best_three_splicesite($model_informations{$i});
		}
		#print Dumper ($model_informations{$i});
	}
	#search startcodon in exon 0
	elsif ($model_informations{$i}->{'exonID'} eq "1"){
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
			until ($model_informations{$i}->{'modulo_start'} == 0 || $model_informations{$i}->{'start_nt'} >= $model_informations{$i}->{'extension_start'} - 8){
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

	#print Dumper ("möp",$model_informations{$i});
	$model_informations{$i} = best_five_splicesite($model_informations{$i});
	#print Dumper ("mäp");
	$model_informations{$i} = best_three_splicesite($model_informations{$i});		
		#print Dumper ("mäp",$model_informations{$i});	

	if ($model_informations{$i}-> {'score_3splicesite'} < 0 || $model_informations{$i}-> {'score_5splicesite'} < 0){
		$model_informations{$i} = extend_search_radius($model_informations{$i});
		if ($model_informations{$i}-> {'score_3splicesite'} < 0){
			$model_informations{$i} = best_three_splicesite($model_informations{$i});
		}	
		if ($model_informations{$i}-> {'score_5splicesite'} < 0){
			$model_informations{$i} = best_five_splicesite($model_informations{$i});
		}
		if ($model_informations{$i}-> {'score_3splicesite'} < 0 || $model_informations{$i}-> {'score_5splicesite'} < 0){
			$model_informations{$i} = extend_search_radius($model_informations{$i});
			if ($model_informations{$i}-> {'score_3splicesite'} < 0){
				$model_informations{$i} = best_three_splicesite($model_informations{$i});
			}	
			if ($model_informations{$i}-> {'score_5splicesite'} < 0){
				$model_informations{$i} = best_five_splicesite($model_informations{$i});
			}
		}
	}
		#$pm->finish
}

#print "2";
#print "Runtime ", duration(time() - $start_time), ".\n";

my $model_informations_ref = \%model_informations;
#print Dumper(%model_informations);
%model_informations = splitted_as($model_informations_ref);
			
for my $i (keys %model_informations) {	
	$model_informations{$i}->{'start'} = $model_informations{$i}->{'start'} - $model_informations{$i}->{'three_nt_extension'};
	$model_informations{$i}->{'end'} = $model_informations{$i}->{'end'} + $model_informations{$i}->{'five_nt_extension'};
	if ($model_informations{$i}->{'strand'} eq "-")
	{
		my $help = $model_informations{$i}->{'start'};
		$model_informations{$i}->{'start'} = $model_informations{$i}->{'end'};
		$model_informations{$i}->{'end'} = $help;
	}
}


#create ExonMatchSolver1.input
#first line: #TargetName     Strandinformation       Model   E-Value Bitscore        Bias    ParalogID       ExonID
if ($mode eq "target"){
	open (ExonMatchSolver1_input_file, '>>', "ExonMatchSolver1.input") or die "$!";
	print ExonMatchSolver1_input_file "#TargetName\tStrandinformation\tModel\tE-Value\tBitscore\tBias\tParalogID\tExonID\n";

	for my $i (keys %model_informations) {	
		#print Dumper (%model_informations{$i});
		# print Dumper ($model_informations{$i}->{'score_3splicesite'});
		# print Dumper ($model_informations{$i}->{'score_5splicesite'});
		# print Dumper ($model_informations{$i}->{'startcodon'});
		# print Dumper ($model_informations{$i}->{'stopcodon'});
		if (defined ($model_informations{$i}->{'score_3splicesite'} and $model_informations{$i}->{'score_5splicesite'} and $model_informations{$i}->{'stopcodon'} and $model_informations{$i}->{'startcodon'})){
			if (($model_informations{$i}->{'score_3splicesite'} > 0 and $model_informations{$i}->{'score_5splicesite'} > 0) or ($model_informations{$i}->{'score_3splicesite'} > 0 and $model_informations{$i}->{'stopcodon'} == 1) or ($model_informations{$i}->{'startcodon'} == 1 and $model_informations{$i}->{'score_5splicesite'} > 0)){
				print ExonMatchSolver1_input_file "$model_informations{$i}->{'targetname'}\t$model_informations{$i}->{'start'}_$model_informations{$i}->{'end'}\t$model_informations{$i}->{'model'}\t$model_informations{$i}->{'evalue'}\t$model_informations{$i}->{'bitscore'}\t$model_informations{$i}->{'bias'}\t$model_informations{$i}->{'paralogID'}\t$model_informations{$i}->{'exonID'}\n";
				#print ExonMatchSolver1_input_file "$model_informations{$i}->{'targetname'}\t0_100\tabc\t0\t5\t0\t1\t4\n";
		}
		}
	}
	close ExonMatchSolver1_input_file;

	#print "4";
	#print "Runtime ", duration(time() - $start_time), ".\n";
}
# fasta output for every paralog 


#------------------------------------------------------------------------
# Aufgaben
#

#output: splice sites (canonical (noncanonical)), exons, gesplittete As

#A: ueberlappende Isoforme gene locus aufheben

#A: uebersetzung As -> nt
#A: blastoutput uebersetzt?


#------------------------------------------------------------------------
# SUBROUTINES
#

#find start, end, paralog name from blastout
sub search_data_blast{
	my ($blastoutIn_line) = @_;
	chomp $blastoutIn_line;

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
	$bitscore =~ s/ //;
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

	if ($strand eq "-"){
		$start = $line_array[9];
		$end = $line_array[8];
	}
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

#create a file for each target sequence
sub create_TargetSequences{
	my ($blastoutIn,$genomeIn) = @_;
	unless (-d "TargetSequences/"){
		mkdir "TargetSequences/";
	}

	my %wantedSequences;

	open (my $blastout_handle, "<", $blastoutIn) or die "Can't open Blastout, $!";
	while (defined (my $line = <$blastout_handle>)){
		chomp $line;
		my @line_array = split "\t", $line;

		$wantedSequences{$line_array[1]}=1; 
	}

	foreach my $keys (keys %wantedSequences){				
		unless (-s "TargetSequences/${keys}.fa"){
			`fastacmd -d $genomeIn -s $keys| sed \'s/>lcl|/>/g\' > TargetSequences/${keys}.fa`;
			`formatdb -i  TargetSequences/${keys}.fa -p F -o`
		}	
	}
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
			$model_informations->{testsequence_three_splicesite} = join ("",@testsequence_three_splicesite);
			if ($model_informations->{testsequence_three_splicesite} =~ m/N+/){
				#print Dumper("mist");
				$model_informations->{score_3splicesite} = -300;	
				my $key_for_hash = $model_informations->{score_3splicesite};
				$three_scores{$key_for_hash} = $model_informations->{three_nt};
			}					
			elsif (scalar @testsequence_three_splicesite == 23){
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
				$model_informations->{testsequence_three_splicesite} = join ("",@testsequence_three_splicesite);
				if ($model_informations->{testsequence_three_splicesite} =~ m/N+/){
					#print Dumper("mist");
					$model_informations->{score_3splicesite} = -300;	
					my $key_for_hash = $model_informations->{score_3splicesite};
					$three_scores{$key_for_hash} = $model_informations->{three_nt};
				}	
				elsif (scalar @testsequence_three_splicesite == 23){
		   			#print Dumper ($model_informations->{testsequence_three_splicesite});
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
					$model_informations->{testsequence_three_splicesite} = join ("",@testsequence_three_splicesite);

					if ($model_informations->{testsequence_three_splicesite} =~ m/N+/){
						#print Dumper("mist");
						$model_informations->{score_3splicesite} = -300;	
						my $key_for_hash = $model_informations->{score_3splicesite};
						$three_scores{$key_for_hash} = $model_informations->{three_nt};
					}
					elsif (scalar @testsequence_three_splicesite == 23){
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
	#$model_informations->{startcodon} = 0;

	if ($model_informations->{five_splicesite} == 0){		
		my %five_scores;	
		my ($five_splicesite, $five_nt) = five_splice_site_canonical(\@a_short, $model_informations->{mismatch}, $model_informations->{gap}, $model_informations->{extension_end});
		#print Dumper($five_splicesite);
		$model_informations->{five_nt} = $five_nt;
		$model_informations->{five_splicesite} = $five_splicesite;
		#print Dumper($model_informations->{five_nt});
		if ($five_splicesite == 1){
			#print Dumper ("1");
			my @testsequence_five_splicesite = sequence_five_splicesite(\@a_short, $five_nt, $model_informations->{extension_end}); 
			if (defined $testsequence_five_splicesite[8]){
				$model_informations->{testsequence_five_splicesite} = join ("",@testsequence_five_splicesite);
				if ($model_informations->{testsequence_five_splicesite} =~ m/N+/){
					$model_informations->{score_5splicesite} = -500;
				}						
				else{
					open (testsequences_five_splicesite_file, '>>', "MaxEntScan/Testsequences_five_splicesite.txt") or die "$!";
					print testsequences_five_splicesite_file join ("",@testsequence_five_splicesite);
				}
			}
			else {
				print "unvollstaendig 5: $model_informations->{model}\n";
				#print Dumper ($model_informations->{extension_end});
				#print Dumper ($model_informations->{extension_start});
				$model_informations->{five_splicesite} = 0;
				$model_informations->{score_5splicesite} = -555;
				unless($model_informations->{extension_end} < 35){
					$model_informations = extend_search_radius($model_informations);
					$model_informations = best_five_splicesite($model_informations);
				}
				return ($model_informations);
			}
			if(-e "MaxEntScan/Testsequences_five_splicesite.txt"){
				my $output_5splicesite = qx(perl $path_to_maxentscan_5 "MaxEntScan/Testsequences_five_splicesite.txt");
				my @array_5splicesite = split /\t/, $output_5splicesite;
				#print Dumper($output_5splicesite);
				if (defined $array_5splicesite[1]){
					chomp ($array_5splicesite[1]);
					$model_informations->{score_5splicesite} = $array_5splicesite[1];							
					my $key_for_hash = $model_informations->{score_5splicesite};
					$five_scores{$key_for_hash} = $model_informations->{five_nt};
				}
				unlink "MaxEntScan/Testsequences_five_splicesite.txt";
			}		
		}
		if ($five_splicesite == 0 || $model_informations->{score_5splicesite} < 0 ){
			#print Dumper ($model_informations->{five_nt});
			$model_informations->{score_5splicesite} = -555;
			until ($model_informations->{score_5splicesite} > 0 || $model_informations->{five_nt} >= $model_informations->{extension_end} -4){			
				print Dumper ("3",$model_informations->{five_nt} );
				my ($five_splicesite, $five_nt) = five_splice_site_canonical_after_negative_score(\@a_short, $model_informations->{five_nt}, $model_informations->{extension_end});
				$model_informations->{five_nt} = $five_nt;
				$model_informations->{five_splicesite} = $five_splicesite;
				if ($five_splicesite == 1){
					my @testsequence_five_splicesite = sequence_five_splicesite(\@a_short, $model_informations->{five_nt}, $model_informations->{extension_end});
					if (scalar @testsequence_five_splicesite == 9){
						$model_informations->{testsequence_five_splicesite} = join ("",@testsequence_five_splicesite);
						if ($model_informations->{testsequence_five_splicesite} =~ m/N+/){
							$model_informations->{score_5splicesite} = -500;
						}				
						else{
							open (testsequences_five_splicesite_file, '>>', "MaxEntScan/Testsequences_five_splicesite.txt") or die "$!";
							print testsequences_five_splicesite_file join ("",@testsequence_five_splicesite);
						}

						my $output_5splicesite = qx(perl $path_to_maxentscan_5 "MaxEntScan/Testsequences_five_splicesite.txt");
						my @array_5splicesite = split /\t/, $output_5splicesite;
						#print Dumper($output_5splicesite);				
						if (defined $array_5splicesite[1]){
							chomp ($array_5splicesite[1]);
							$model_informations->{score_5splicesite} = $array_5splicesite[1];
							my $key_for_hash = $model_informations->{score_5splicesite};
							$five_scores{$key_for_hash} = $model_informations->{five_nt};
						}
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
						if ($model_informations->{testsequence_five_splicesite} =~ m/N+/){
							$model_informations->{score_5splicesite} = -500;
						}
						else{
							open (testsequences_five_splicesite_file, '>>', "MaxEntScan/Testsequences_five_splicesite.txt") or die "$!";
							print testsequences_five_splicesite_file join ("",@testsequence_five_splicesite);
						}
						my $output_5splicesite = qx(perl $path_to_maxentscan_5 "MaxEntScan/Testsequences_five_splicesite.txt");
						my @array_5splicesite = split /\t/, $output_5splicesite;				
						if (defined $array_5splicesite[1]){
							#print Dumper($output_5splicesite);
							chomp ($array_5splicesite[1]);
							$model_informations->{score_5splicesite} = $array_5splicesite[1];
							my $key_for_hash = $model_informations->{score_5splicesite};
							$five_scores{$key_for_hash} = $model_informations->{five_nt};
						}
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
		#print Dumper ($i);
        if (($genome_target[$i] eq "T" and $genome_target[$i+1] eq "A" and $genome_target[$i+2] eq "G") or ($genome_target[$i] eq "T" and $genome_target[$i+1] eq "G" and $genome_target[$i+2] eq "A") or ($genome_target[$i] eq "T" and $genome_target[$i+1] eq "A" and $genome_target[$i+2] eq "A")) {     
			#print Dumper ($genome_target[$i], $stop_nt);
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
        	$startcodon = 1;    
        	#print Dumper ($start_nt);
        	return ($startcodon, ($start_nt + 1));
        	#$start_nt++;
        } 
        else {
        	$startcodon = 0;
        	$start_nt++;
        }
    }
    return ($startcodon, $start_nt);
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
    my $five_splicesite = 0; #1 if splicesite is found, else 0
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

	for(my $i = ((scalar @genome_target) - $extension +1 + $five_nt);$i <= (scalar @genome_target) - 2;$i++){
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

sub extend_search_radius{
	my ($model_informations) = @_ ;
	my $extension_2 = $model_informations->{extension_start} + 18;
	if ($model_informations->{start_query} >= 9 || $model_informations->{query_size} - $model_informations->{end_query} >= 9){
		$extension_2 += 12;
	}
	my @a_short;

		my $genome_seqIO = Bio::SeqIO->new(-file   => "TargetSequences/".$model_informations->{targetname}.".fa", -format => "Fasta");
		$a = $genome_seqIO->next_seq;
		$a = $a -> seq;

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
			#print Dumper ($i, $a[$i]);
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
	$model_informations->{genome_part} = join ("",@a_short);
	return ($model_informations);
}