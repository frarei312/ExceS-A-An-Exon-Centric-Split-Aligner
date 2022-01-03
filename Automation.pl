#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use autodie qw(:all);
use Data::Dumper;
use Bio::SeqIO;
use Bio::Perl;
use Bio::Seq;
use List::Util qw(any);
use Bio::AlignIO;
use File::Slurp qw(read_file);

#------------------------------------------------------------------------
#Definitions
#

#Set this variables to your path once and save the script! 					#
#
my $pathtotargetgenomeIn = '';
my $pathtoquerygenomeIn = "";	
my $EMSpipeline = "";			
#
#Set this variable to the name of the query once and save the script! 					#
#
my $query_name ="";

#------------------------------------------------------------------------

my $paraMax = "";
my $querygenomeIn = "";
my $paralogsIn = "";
my $targetgenomesIn = "";
my $firsttargetgenomeout = "";
my $m = "";
my $manualCuration = "";
my $round = 1;
my %paralogexonsequence;
my %model_informations_results;

GetOptions("Paralogs=s" => \$paralogsIn, "TargetGenomes=s" => \$targetgenomesIn, "paraMax=i" => \$paraMax, "m=s" => \$m, "manualCuration=s" => \$manualCuration, "queryName=s" => \$query_name);

#------------------------------------------------------------------------
#Main
#


#paralognames in array

open (my $paralogshandle, "<", $paralogsIn) or die "Can't open ParalogsIn, $!";
my @paralogs = qx(grep ">" $paralogsIn);
my %paralogsequence;
my @notusedparalogssave;

my $seqio = Bio::SeqIO->new(-file => "$paralogsIn", -format => "fasta");
while(my$seqobj = $seqio->next_seq) {
    my $id  = $seqobj->display_id;
    my $seq = $seqobj->seq;
    $paralogsequence{$id} = $seq;
}
foreach my $paralog(@paralogs){
	chomp $paralog;
	my @paralog_split = (split/\s/, $paralog);
	$paralog_split[0] =~ s/>//;
	push @notusedparalogssave, $paralog_split[0];
	open (my $paralogoutput,">>", "${paralog_split[0]}.fa") or die "Can't open Paralogoutput, $!";
	print $paralogoutput "${paralog}_${query_name}\n";
	print $paralogoutput "$paralogsequence{$paralog_split[0]}\n";
}

my $out1;
my $out2;
my $usedlines = 0;
my $genome;

#outfile definieren
open (my $genome_names, "<" , "$targetgenomesIn") or die "Cannot open directory: $!";
while (defined (my $line = <$genome_names>)){
	chomp $line;
	$genome = $line;
	($out1, $out2) = (split/\./, $line)[0, 1];	
	my @notusedparalogs = @notusedparalogssave;

	if ($usedlines == 0){
		$firsttargetgenomeout = "EMS_${out1}_${out2}";
	}

	unless ($usedlines > 4){
		unless (-e "EMS_${out1}_${out2}"){
			print "perl $EMSpipeline -i $paralogsIn -o EMS_${out1}_${out2}  -mode fasta -target $pathtotargetgenomeIn/$line -query $pathtoquerygenomeIn -paraMax $paraMax\n";

			system("perl $EMSpipeline -i $paralogsIn -o EMS_${out1}_${out2}  -mode fasta -target $pathtotargetgenomeIn/$line -query $pathtoquerygenomeIn -paraMax $paraMax");


			unless (-e "Paralog_exon_alignments/"){
				mkdir "Paralog_exon_alignments/";

				my $query;
				open (EXONFILE_HANDLE, "<" , "EMS_${out1}_${out2}/HomologousExons/paralog_exon_alignment.fa") or die "Cannot open directory: $!";
				while (defined (my $line = <EXONFILE_HANDLE>)){
					chomp $line;
					if($line =~ m/>/){
						$line =~ s/>//;
						$query = "${line}_${query_name}";
						open (EXONOUT,">> Paralog_exon_alignments/${line}.fa") or die "Can't open Paralog_exon_alignments/${line}, $!";
						print EXONOUT ">${line}_${query_name}\n";
					}
					else{
						print EXONOUT "$line\n";		
						close (EXONOUT);
					}
					$paralogexonsequence{$query} = $line;
				}
			}

			if (-e "EMS_${out1}_${out2}/Summary.fa"){
				my %model_informations;
				my $seqio = Bio::SeqIO->new(-file => "EMS_${out1}_${out2}/Summary.fa", -format => "fasta");
				while(my $seqobj = $seqio->next_seq) {
				    my $id  = $seqobj->display_id;		    					
		   			$id =~ s/^SearchTarget\///;
					$id =~ s/_final//;
					$id =~ s/_initial//;
					my @id_array = split/\_/, $id;
					$id_array[0] =~ s/^[0-9]+\-//;
				    my $seq = $seqobj->seq;

				    if (exists $model_informations{$id_array[0]}){
				    	if ((length $seq) > (length $model_informations{$id_array[0]}->{'sequence'})){
				    		my %id = (		
								id => $id,
									sequence => $seq,
							);
							@notusedparalogs = grep {$_ ne $id_array[0]} @notusedparalogs;
							my $key_for_hash = $id_array[0];
							$model_informations{$key_for_hash} = \%id;
				    	}
				    }
				    else{					
				    	my %id = (		
							id => $id,
							sequence => $seq,
						);
						@notusedparalogs = grep {$_ ne $id_array[0]} @notusedparalogs;
						my $key_for_hash = $id_array[0];
						$model_informations{$key_for_hash} = \%id;			
						#print Dumper (%model_informations);
				    }
				}
				for my $entry (keys %model_informations){
					#print Dumper ($entry);
						divide_exons($model_informations{$entry});						
						my $key_for_hash = $model_informations{$entry}->{id}."_${out1}_${out2}";
						$model_informations_results{$key_for_hash} = $model_informations{$entry}->{sequence};
				}
			}#print Dumper (%model_informations_results);
		}
		$usedlines++;
	}

# second step
# find orthologs in other related species
#

	else{
		unless (-e "EMS_${out1}_${out2}"){
			print "perl $EMSpipeline -i $paralogsIn -m $m -o EMS_${out1}_${out2} -mode alignment -user ${firsttargetgenomeout}/HomologousExons/paralog_exon_alignment_final.fa -target $pathtotargetgenomeIn/$line -paraMax $paraMax\n";

			system("perl $EMSpipeline -i $paralogsIn -m $m -o EMS_${out1}_${out2} -mode alignment -user ${firsttargetgenomeout}/HomologousExons/paralog_exon_alignment_final.fa -target $pathtotargetgenomeIn/$line -paraMax $paraMax");

			if (-e "EMS_${out1}_${out2}/Summary.fa"){
				my %model_informations;
				my $seqio = Bio::SeqIO->new(-file => "EMS_${out1}_${out2}/Summary.fa", -format => "fasta");
				while(my $seqobj = $seqio->next_seq) {
				    my $id  = $seqobj->display_id;		    					
		   			$id =~ s/^SearchTarget\///;
					$id =~ s/_final//;
					$id =~ s/_initial//;
					my @id_array = split/\_/, $id;
					$id_array[0] =~ s/^[0-9]+\-//;
				    my $seq = $seqobj->seq;

				    if (exists $model_informations{$id_array[0]}){
				    	if ((length $seq) > (length $model_informations{$id_array[0]}->{'sequence'})){
				    		my %id = (		
								id => $id,
									sequence => $seq,
							);
							@notusedparalogs = grep {$_ ne $id_array[0]} @notusedparalogs;
							my $key_for_hash = $id_array[0];
							$model_informations{$key_for_hash} = \%id;
				    	}
				    }
				    else{					
				    	my %id = (		
							id => $id,
							sequence => $seq,
						);
						@notusedparalogs = grep {$_ ne $id_array[0]} @notusedparalogs;
						my $key_for_hash = $id_array[0];
						$model_informations{$key_for_hash} = \%id;
				    }
				}
				for my $entry (keys %model_informations){
						divide_exons($model_informations{$entry});
						my $key_for_hash = $model_informations{$entry}->{id}."_${out1}_${out2}";
						$model_informations_results{$key_for_hash} = $model_informations{$entry}->{sequence};
				}
				#print Dumper (%model_informations_results);
			}
		}
	}

	if (scalar @notusedparalogs > 0){
		open (my $searchagain,">>", "SearchAgain") or die "Can't open SearchAgain, $!";
		#print $searchagain "$genome\n";
	}
}

if ($manualCuration eq "yes"){
		print "Please curate manually and restart.\n";
			#print Dumper (%model_informations_results);
		for my $entry (keys %model_informations_results){
			my $name = (split /_/, $entry)[0];
			$name =~ s/^[0-9]+-//;
			#print Dumper ($name);
			open (my $paralog_result,">>", "${name}.fa") or die "Can't open Paralogoutput, $!";
			print $paralog_result ">$entry\n$model_informations_results{$entry}\n";
		}
		exit;
}

else{
	if (-e "SearchAgain"){
		my $file_string_new = read_file("SearchAgain");
		my $file_string_old = "";

		while ($file_string_new ne $file_string_old) {
			my ($firstout1, $firstout2, $firstout3) = (split/_/, $firsttargetgenomeout)[1, 2, 3];			
			open (my $SearchAgain, "<" , "SearchAgain") or die "Cannot open directory: $!";
			while (defined (my $line = <$SearchAgain>)){
				chomp $line;
				($out1, $out2) = (split/\./, $line)[0, 1];

				if ($line =~ m/${firstout1}_${firstout2}.${firstout3}/){
					#print Dumper ("first");
				}
				else{
					`rm -r EMS_${out1}_${out2}`;
				
				open (my $searchagain_round,">>", "SearchAgain_${round}") or die "Can't open SearchAgain, $!";
				print $searchagain_round "$line\n";
				}	
			}

			`rm SearchAgain`;
			open (my $newround, "<" , "SearchAgain_${round}") or die "Cannot open directory: $!";
			$round++;
			while (defined (my $line = <$newround>)){
				chomp $line;
				$genome = $line;
				($out1, $out2) = (split/\./, $line)[0, 1];	
				my @notusedparalogs = @notusedparalogssave;
				unless (-e "EMS_${out1}_${out2}"){
					print "perl $EMSpipeline -i $paralogsIn -m $m -o EMS_${out1}_${out2} -mode alignment -user ${firsttargetgenomeout}/HomologousExons/paralog_exon_alignment_final.fa -target $pathtotargetgenomeIn/$line -paraMax $paraMax\n";

					system("perl $EMSpipeline -i $paralogsIn -m $m -o EMS_${out1}_${out2} -mode alignment -user ${firsttargetgenomeout}/HomologousExons/paralog_exon_alignment_final.fa -target $pathtotargetgenomeIn/$line -paraMax $paraMax");

					if (-e "EMS_${out1}_${out2}/Summary.fa"){
						my %model_informations;
						my $seqio = Bio::SeqIO->new(-file => "EMS_${out1}_${out2}/Summary.fa", -format => "fasta");
						while(my $seqobj = $seqio->next_seq) {
						    my $id  = $seqobj->display_id;		    					
				   			$id =~ s/^SearchTarget\///;
		 					$id =~ s/_final//;
							$id =~ s/_initial//;
							my @id_array = split/\_/, $id;
							$id_array[0] =~ s/^[0-9]+\-//;
						    my $seq = $seqobj->seq;

						    if (exists $model_informations{$id_array[0]}){
						    	if ((length $seq) > (length $model_informations{$id_array[0]}->{'sequence'})){
						    		my %id = (		
										id => $id,
		    								sequence => $seq,
		    						);
		    						@notusedparalogs = grep {$_ ne $id_array[0]} @notusedparalogs;
									my $key_for_hash = $id_array[0];
									$model_informations{$key_for_hash} = \%id;
						    	}
						    }
						    else{					
						    	my %id = (		
									id => $id,
		    						sequence => $seq,
		    					);
		    					@notusedparalogs = grep {$_ ne $id_array[0]} @notusedparalogs;
								my $key_for_hash = $id_array[0];
								$model_informations{$key_for_hash} = \%id;
						    }
						}
						for my $entry (keys %model_informations){
								divide_exons($model_informations{$entry});						
								my $key_for_hash = $model_informations{$entry}->{id}."_${out1}_${out2}";
								$model_informations_results{$key_for_hash} = $model_informations{$entry}->{sequence};
						}
					}
				}
				if (scalar @notusedparalogs > 0){
					open (my $searchagain,">>", "SearchAgain") or die "Can't open SearchAgain, $!";
					print $searchagain "$genome\n";
				}
				#print Dumper (%model_informations_results);
			}			
			if (-e "SearchAgain"){
				$file_string_old = $file_string_new;
				$file_string_new = read_file("SearchAgain")
			}	
			else{
				last;
			}
		}
	}				
	for my $entry (keys %model_informations_results){
		my $name = (split /_/, $entry)[0];
		$name =~ s/^[0-9]+-//;
		my $sequence = $model_informations_results{$entry};
		$sequence =~ s/\*/X/;
		open (my $paralog_result,">>", "${name}.fa") or die "Can't open Paralogoutput, $!";
		print $paralog_result ">$entry\n$sequence\n";
	}

	print "Search completed.\n";
}
#------------------------------------------------------------------------
# SUBROUTINES
#	


sub divide_exons{

	my ($model_informations_ref) = @_;
	my %model_informations_handle = %$model_informations_ref;
	my $paralogname = $model_informations_handle{id};
	my @paralogname_array = split/\_/, $paralogname;
	#print Dumper ($paralogname_array[0]);
	my $paralog_search = $paralogname_array[0];
	$paralog_search =~ s/^[0-9]+\-//;

	my $exonnumber;
	my @ExonFiles = glob ("EMS_${out1}_${out2}/SearchTarget/".${paralogname_array}[0]."*_exons");
	foreach my $file(@ExonFiles){
		open(EXONS, $file) or die "Could not open $file\n";
		while (defined (my $line = <EXONS>)){
			chomp $line;
			#print Dumper ($line);
			#if($line =~ m/${paralog_search}/){
				#print Dumper("yes");
				if($line =~ m/>/){
					my @line = split/\_/, $line;
					$exonnumber = $line[1];
					#print Dumper($exonnumber);
					open (EXONFILES, ">> Paralog_exon_alignments/${paralog_search}_${exonnumber}.fa") or warn "Could not read ${paralog_search}_${exonnumber}.fa\n";
					my $speciesname = "${out1}_${out2}";			
					$line =~ s/SearchTarget\///;
				 	$line =~ s/_final//;
					$line =~ s/_initial//;
					print EXONFILES "${line}_${speciesname}\n";
				}

				else{				   		
					$line =~ s/\*/X/;	
					print EXONFILES "$line\n";
					close(EXONFILES);
					qx(muscle3.8.31_i86linux64 -in Paralog_exon_alignments/${paralog_search}_${exonnumber}.fa -out Paralog_exon_alignments/${paralog_search}_${exonnumber}.fasta -quiet);
				}
			#}
		}
		close(EXONS);
	}	
}