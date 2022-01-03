#!/usr/bin/perl -w
use strict;
use Bio::Perl;
use Bio::SeqIO;
use Bio::Seq;
use Data::Dumper;
use Bio::DB::Fasta;
use Statistics::Standard_Normal qw/z_to_pct pct_to_z/;
use Getopt::Long;
use List::Util;
use FindBin qw($Bin);
use Cwd;
use DateTime;
use Parallel::ForkManager;

#------------------------------------------------------------------------------
#DEPENDENCIES
#perl-modules (see above)
#if mode = fasta: blastall, formatdb, clustalw;
#if mode = alignment: getORF, hmmer (hmmbuild, hmmsearch);
#for all modi: bl2seq, fastacmd;
#linux command line: grep, awk, sed, wc;

# -----------------------------------------------------------------------------
# GLOBALS

#################################################################################################
#Set this variables to your path once and save the script! 					#
#
#Path to ExonMatchSolver5.jar									#	
my $Bin = "";											#

#Path to procompart
my $path = "";
#
#Path to sensitive HMMER (alignment-mode only):							#
my $hmmer_path = "";										#
#
#Path to a scratch with sufficient file space for a translated genome (alignment-mode only) 	#
my $space = "";											#
#	
#Path to ExceS_A_1:							#
my $ExceS_A_1 = "";		#
#
#Path to ExceS_A_2:							#
my $ExceS_A_2 = "";
#Path to ExceS_A_3:							#
my $ExceS_A_3 = "";
#Path to ExceS_A_4:							#
my $ExceS_A_4 = "";

my $distance_first = 1500000; 	#Can be defined by the user in dependence of the intron size expected for the first (few) exons.   
my $distance_last = 1500000; #Can be defined by the user in dependence of the intron size expected for the last (few) exon.
#################################################################################################


use vars qw ($help $genome_target $genome_query $input $outfolder);
my $core = 10;			#Numbers of cores available for the pipeline
my $intervall = "";		#Option handed over to ExonMatchSolver
my $max_solution = "";		#Option handed over to ExonMatchSolver
my $WGD = "no";
my $input2 = "";
my $mode = "fasta";
my $exonerate = "no";
my $length_cutoff = undef;      #Minimal amino acid length of TCEs to be included in the homology assignment. 
my $models = "";
my $z_cutoff = 3;		
my $nonDuplicated = "";
my $maxparalogs = "";

# -----------------------------------------------------------------------------
# OPTIONS

GetOptions (
    "i=s"          => \$input,
    "user=s"       => \$input2,
    "m=s"	       => \$models,
    "mode=s"       => \$mode,
    "o=s"          => \$outfolder,
    "WGD=s"        => \$WGD,
    "noWGD=s"      => \$nonDuplicated,
    "paraMax=i"    => \$maxparalogs,	
    "query=s"      => \$genome_query,
    "target=s"     => \$genome_target,
    "dFirst=i"     => \$distance_first,
    "dLast=i"      => \$distance_last,
    "scipio_opt=s" => \$scipio_options,
    "l=i"          => \$length_cutoff,
    "z=f"          => \$z_cutoff,
    "c=i"	       => \$core,
    "s=s"	       => \$exonerate,
    "int=i"	       => \$intervall,
    "max=i"        => \$max_solution,
    "help"         => \$help,
    "h"            => \$help);
usage() if ($help || !$input || !$outfolder|| !$genome_target|| !$mode|| !$maxparalogs);


# -----------------------------------------------------------------------------
# MAIN

my %loci;
my $strand;
my $start;
my $end;
my $genome_query = $genome_query;
my $genome_target = $genome_target;
my $genome_target_original = $genome_target;
my $paralog;
my $exon;
my %missing_scores;
my @numle;
my $e = 1;
my %single_exons;
my $myInput = $input;
my $myInput2 = $input2;
my $myOutFolder = $outfolder;
my @exons;
my @sorted_exons;
my $blast;
my $hmm;
my @IDs;
my $cutpoints;

#Input: protein sequence or nucleotide sequence from Swiss-Prot/ RefSeq in fasta-format
#Test existence of target and query genome, build blast database if that does not exist
#Create folder with name of choice (-o), copy input file (-i) in this folder

die "The folder \"$myOutFolder\" already exists. Please change the name of the outputfolder (-o).\n" 
    if -e $myOutFolder;
printHeader();
mkdir "$myOutFolder/";

die "The target genome \"$genome_target\" does not exist. $!\n"
    unless -e $genome_target;
#if ($mode eq "fasta"|| $mode eq "user"){
    my @buildtarget = "formatdb -o -i $genome_target -p F";
    unless (-e "${genome_target}.nhr" && -e "${genome_target}.nin" && -e "${genome_target}.nsd" && -e "${genome_target}.nsi" && -e "${genome_target}.nsq"){
	system (@buildtarget) == 0 or die "Formatdb with target genome failed. $!\n";
    };
#};

my $paralog_num = `grep '>' $myInput| wc -l`;
my @numbers = (1..$paralog_num);

`cp $myInput $myOutFolder/.`;
chdir "$myOutFolder/";
my $dir = getcwd();
mkdir "SearchQuery/";
mkdir "HomologousExons/";
mkdir "SingleExons/";

my $in  = new Bio::SeqIO(-file  => "$myInput");
while (my $seq = $in->next_seq) {
    push @IDs, $seq->display_id;
};

#IDs of genes of interest are stored in @IDs
foreach my $ID(@IDs){
    my $out = Bio::SeqIO->new(-file => ">SearchQuery/${ID}.fa",
			      -format=>'fasta'
	);
    my $in  = new Bio::SeqIO(-file  => "$myInput");                                
    while (my $seq = $in->next_seq) {
	if ($seq->display_id =~ /^$ID$/) {
	    $out->write_seq($seq);
	};         	
    };                                
}

##################In fasta mode, the gene (family) of interest has first to be identified in the query genome and homologous TCEs have to be identified
##################
if ($mode eq "fasta"){
    die "The query genome \"$genome_query\" does not exist. $!\n"
	unless -e $genome_query;
    #print "$genome_query\n";
    my @buildquery = "formatdb -i ${genome_query} -p F -o";
    unless (-e "${genome_query}.nhr" && -e "${genome_query}.nin" && -e "${genome_query}.nsd" && -e "${genome_query}.nsi" && -e "${genome_query}.nsq"){
	system (@buildquery) == 0 or die "Formatdb with query genome failed. $!\n";
    }

    #Blast the query sequences against the query genome
    foreach my $ID(@IDs){
	my @call = "blastall -p tblastn -d $genome_query -i SearchQuery/${ID}.fa -a $core -m 8 -G 12 -E 1 -o SearchQuery/${ID}.blastout";
	system (@call) == 0 or die "blastall against query genome failed. $!\n";
    };
    #print "Blast of query sequences against query genome done:\n";

    #Extract the corresponding locus from the blast output
    foreach my $ID(@IDs){
	my $locus = `sort -k11,11g SearchQuery/${ID}.blastout| awk '{if (\$3>90) print \$2}'| head -n1`;
	chomp $locus;
	#print "Paralog ${ID} is located on locus $locus.\n";
	$loci{$ID} = $locus;
    }

	foreach my $ID(@IDs){
		print Dumper ("perl $ExceS_A_3 --genome $genome_query --locus $loci{$ID} --ID $ID --querysequence SearchQuery/${ID}.fa");
		qx(perl $ExceS_A_3 --genome $genome_query --locus $loci{$ID} --ID $ID --querysequence SearchQuery/${ID}.fa);
	}
    #Retrieve the intron-exon structure for each paralog ($count), the accumulated TCE-length are stored in single arrays within a hash ($hash{$count}) (one array per paralog). This is done by the sub get_exons.
    #Lengths are extracted from the .inf-file if prosplign is used and summed up
    #They can be used as start and stopp to extract the subsequence of the original amino acid sequence

    my $Dozscore = 1;
    if (scalar(@numbers) eq "1"){
	$Dozscore = 0;
    }

	my %hash2;
	foreach my $ID(@IDs){
		$hash2{$ID} = `grep ">" SearchQuery/${ID}.inf | tail -1 | cut -d"_" -f2 | sed \'s/exon//g\'`;
	}
	foreach my $count(sort {$a cmp $b} keys %hash2){
		`cat SearchQuery/${count}.inf >> HomologousExons/paralog_exon_alignment.fa`;
	}


    #Pairwise scores are calculated for every TCE pair using Clustalw 
    `clustalw2 -INFILE=HomologousExons/paralog_exon_alignment.fa > HomologousExons/pairwisescores.clw`;
    my $switch = 0;
    my @scores;
    my @array;
    my $value;
    my @new_array;
    my %hash;
    my $x;
    my $y;
    my @length;
    my $length;
    my $length2;
    my @length2;
    my $switch3 = 0;
    my @id;
    my %length;
    my %scores;
    my @otherscores;

    #Parse the clustalw-input file
    #Store the sequence-id as hash key and the aa-length[0] and full sequence[1] name as array in the hash %length
    open (IN, "< HomologousExons/pairwisescores.clw");
    while (my $line = <IN>){
		chomp $line;
		if ($line =~ /Sequence\s(\d)+:/){
		    my @length = split (/\s+/,$line);
		    $_ = $length[1];
		    s/://g;
		    $length[1] = $_;
		    $length{$length[1]} = [int($length[3]), $length[2]]; 		
		};
		if ($line =~ /Aligned\.\sScore/){
		    @scores = split ":", $line;
		    $_ = $scores[0];
		    s/Sequences\s\(//g;
		    $scores[0] = $_;
		    $_ = $scores[1];
		    s/\)\sAligned.\sScore//g;
		    $scores[1] = $_;
		    $_ = $scores[2];
		    s/(\s)+(\d+)/$2/g;
		    $scores[2] =$_;	
		    foreach my $ID(@IDs){
				if ($length{$scores[0]}[1] =~ /$ID/ && $length{$scores[1]}[1] =~ /$ID/){
			   	 push @{$scores{$ID}}, [$scores[0], $scores[1], $scores[2]];
				};	   
			}	
		push @otherscores, [$scores[0], $scores[1], $scores[2]];	
		}
    }
    close IN;    
    #print Dumper (@scores); 
   	#print Dumper (@otherscores);
    my @maximum;
    for my $seq(keys %length){
	push @maximum, $length{$seq}[0];
    } 

    my @sortedmax = sort {$a <=> $b} @maximum;
    my $maxLength = $sortedmax[-1];
    unless (defined $length_cutoff){
	$length_cutoff = ($maxLength / 10);
	#print "Length-cutoff: $length_cutoff\n";
    };

    #Store the first id as a key in the hash %zscores, in this hash, store the second id as key with the score as value 
    my %zscores;

    for my $paralog(keys %scores){
	foreach my $element(0..$#{$scores{$paralog}}){
	    foreach my $number1($scores{$paralog}[$element][0]){	
		foreach my $number2($scores{$paralog}[$element][1]){
		    foreach my $all(0..$#{$scores{$paralog}}){
			if ($scores{$paralog}[$all][0] == $number1){
			    $zscores{$paralog}{$number1}{$scores{$paralog}[$all][1]}=[$scores{$paralog}[$all][2],]; 		
			};
			if ($scores{$paralog}[$all][1] == $number1){
			    $zscores{$paralog}{$number1}{$scores{$paralog}[$all][0]}=[$scores{$paralog}[$all][2],]; 
			}
			if ($scores{$paralog}[$all][0] == $number2){
			    $zscores{$paralog}{$number2}{$scores{$paralog}[$all][1]}=[$scores{$paralog}[$all][2],]; 		
			};
			if ($scores{$paralog}[$all][1] == $number2){
			    $zscores{$paralog}{$number2}{$scores{$paralog}[$all][0]}=[$scores{$paralog}[$all][2],]; 
			}
		    }
		}
	    }
	}
    };

    my %otherzscores;
    foreach my $element(0..$#otherscores){
	foreach my $number($otherscores[$element][0]){	
	    foreach my $all(0..$#otherscores){
		if ($otherscores[$all][0] == $number){
		    $otherzscores{$number}{$otherscores[$all][1]}=[$otherscores[$all][2],]; 
		};
		if ($otherscores[$all][1] == $number){
		    $otherzscores{$number}{$otherscores[$all][0]}=[$otherscores[$all][2],]; 
		}
	    }
	}
    }
	#print Dumper (%otherzscores);
    #Parse %zscores
    #Store the first id as a key of the hash %statistics and calculate mean and standard deviation as values of %statistics 
    my %statistics;
    my @average;
    my @number;
    for my $para(sort {$a cmp $b} keys %zscores){
    	#print Dumper ($para);
	for my $number(sort {$a <=> $b}keys %{$zscores{$para}}){
	    for my $pair(keys %{$zscores{$para}{$number}}){
		push @average, $zscores{$para}{$number}{$pair}[0];
	    };
	    $statistics{$number}=[average(@average),stdev(@average)];	
	    #print Dumper ($statistics{$number});
	}
    }
    #print Dumper (%statistics);
    #Calculate the zscores of every id-pair with the mean and average parsed from %statistics, add as value to %otherzscores
    my @array2;
    my @array3;

    for my $seq1(sort {$a <=> $b}keys %otherzscores){
    	#print Dumper ($seq1);
	for my $seq2(sort {$a <=> $b}keys %{$otherzscores{$seq1}}){
		    	#print Dumper ($seq2);
		#if ($statistics{$seq1}[1] == 0){
		#	$statistics{$seq1}[1] = 0.0000000001;
		#}
		#print Dumper ($otherzscores{$seq1}{$seq2}[0]);
	    $otherzscores{$seq1}{$seq2}[1]= (($otherzscores{$seq1}{$seq2}[0] - $statistics{$seq1}[0])/$statistics{$seq1}[1]);
	    #print Dumper ($otherzscores{$seq1}{$seq2}[1]);
	    if ($otherzscores{$seq1}{$seq2}[1] > $z_cutoff){
		if ($length{$seq1}[0] > $length_cutoff && $length{$seq2}[0] > $length_cutoff){	
		    push @array2, [$seq1, $seq2, $otherzscores{$seq1}{$seq2}[0]];		
		}	
	    }
	}
    }
    #print Dumper (@array2);
    #Filter for only the pairs that are reciprocally included in the set (from both z-score filterings)  
    for my $i (0..$#array2){
	for my $j (0..$#array2){
	    if ($array2[$i][0] == $array2[$j][1] && $array2[$i][1] == $array2[$j][0]){
		push @array3, [$array2[$i][0], $array2[$i][1], $array2[$i][2]];
	    }
	}
    }
    #print Dumper (@array3);
    my @sorted = sort{$a->[0] <=> $b->[0]}@array3;
    #print Dumper (@sorted);
    #Assign the homologous TCEs into groups (hash in hash with groups as keys in %hash, in the groups are the ids keys with occurences as values)
    #Groups have letters as names
    my $switch2 = 0;
    for $x (0..$#sorted){
	$switch2 = 0;
	foreach my $exon(keys %hash){
	    if (exists $hash{$exon}{$sorted[$x][0]}){
		$switch2 = 1;		
		$hash{$exon}{$sorted[$x][0]}++;
		$hash{$exon}{$sorted[$x][1]}++;
	    }elsif (exists $hash{$exon}{$sorted[$x][1]}){
		$switch2 = 1;		
		$hash{$exon}{$sorted[$x][0]}++;
		$hash{$exon}{$sorted[$x][1]}++;	
	    }elsif(!exists $hash{$exon}{$sorted[$x][0]} && !exists $hash{$exon}{$sorted[$x][1]}){
		next;
	    };
	};
	if ($switch2 == 0){
	    $hash{$e}{$sorted[$x][0]}++;
	    $hash{$e}{$sorted[$x][1]}++;
	    $e++;
	};
	
    };
    my $string;
    my %newhash;
    for my $exon(sort keys %hash){
	$string = "";
	for my $seq (sort {$a <=> $b} keys %{$hash{$exon}}){
	    $string .= ":".$seq;
	}
	unless (exists $newhash{$string}){
	    $newhash{$string} = $exon;
	}
    }
    #print Dumper (%newhash);
    my %finalhash;
    my $counter = 1;
    for my $seq (sort {$newhash{$a} <=> $newhash{$b}} keys %newhash){
	$finalhash{$counter} = $seq; 
	$counter ++;
    }

    open (EXONS, "< HomologousExons/paralog_exon_alignment.fa");
    open (OUT, "> HomologousExons/paralog_exon_alignment_final.fa");

    my @exons;    	
    #print Dumper (%finalhash);
    for my $groups(sort keys %finalhash){
	#print Dumper ($groups);
	push @exons, $groups;
	my @sequences = split ":", $finalhash{$groups};
	if ($#sequences-1 > $paralog_num){
	    #print "WARNING: Number of paralogs assigned to one homologous TCE group is greater than number of paralogs. Check assignment of homologous TCEs. Increase length_cutoff (-l).\n";
	}; 
	foreach my $i (1..$#sequences){
		#print Dumper ("i= ".$i);
	    while(<EXONS>){
	    	#print Dumper ("EXON: ".$_);
		chomp;
		if ($switch3 == 1){		    
			#print Dumper ($_);
		    print OUT "$_\n";
		    $switch3 = 0;	
		    last;				
		};
		if ($_ =~ /^>$length{$sequences[$i]}[1]$/){
			#print Dumper ("length: ".$_);
		    s/_exon\d+/_exon$groups/g;	
		    #print Dumper ($_);
		    print OUT "$_\n";
		    $switch3=1;					
		}
		else{
		    $switch3=0;				
		}
	    }
	    seek EXONS, 0, 0;	
	} 
    };
    close OUT;
    close EXONS;

    @sorted_exons = sort{$a <=> $b}@exons;

    if ($Dozscore == 0) {
	`cp HomologousExons/paralog_exon_alignment.fa HomologousExons/paralog_exon_alignment_final.fa`;
	#print "I skip z-score filtering for identification of homologous exons as just one paralog was detected\n";
	open (IN, "< HomologousExons/paralog_exon_alignment_final.fa");
	while (<IN>){
	    chomp $_;
	    if ($_ =~ /(.*)_exon(\d+)/){ 
		push @exons, $2;
	    };
	};
	my @sort_exons = sort{$a <=> $b}@exons;
	@sorted_exons = uniq (@sort_exons);
    }
    else {
    	#print "Homologous TCEs between paralogs identified.\n";
    }
}

elsif ($mode eq "user"|| $mode eq "alignment"){	
	$myInput2 =~ s/_final\.fa//;
	#print Dumper ($myInput2);
    `cp ../${myInput2}_final.fa $dir/HomologousExons/paralog_exon_alignment_final.fa`;
    `cp ../${myInput2}.fa $dir/HomologousExons/paralog_exon_alignment.fa`;
    open (IN, "< HomologousExons/paralog_exon_alignment_final.fa") or die "Please provide an additional input file with -user <file> in alignment and user mode.\n";
    while (<IN>){
	if ($_ =~ /(.*)_exon(\d+)/){ 
	    push @exons, $2;
	};
    };
    my @sort_exons = sort{$a <=> $b}@exons;
    @sorted_exons = uniq (@sort_exons);
}

else{die "Specify input file with -user <file> in alignment or user mode\n."};


################# Initial homology search starts here
#################


if ($mode eq "alignment"){
    # translate the genomic sequence into 6 ORFs, stored in ${space}${myOutFolder}_ORF.aas (amino-acid sequence) with getorf
    chdir $space or die "Couldn't access $space\n";
    system "awk '{print \$1}' $genome_target > tmp1";
    unless (-e "${space}${myOutFolder}_ORF.aas" && -e $genome_target) {
	my $ret = system "getorf -minsize 9 -reverse Y -stdout Y -sequence tmp1 -outseq stdout|sed 's/ - /-/g'| awk '{print \$1, \$3, \$2}' > tmp2";
	if ($ret != 0){
	    die;
	};
	open (TMP, "< tmp2") or die "Couldn't open tmp2, $!";
	open (DEST, ">> ${space}${myOutFolder}_ORF.aas") or die "Couldn't open ${space}${myOutFolder}_ORF.aas, $!";
	while (<TMP>){
	    chomp;
	    s/\(//;
	    s/_(\d+)\s\s\[/ $1\(FORWARD\)\[/;
	    s/\sREVERSE/\(REVERSE\)/;
	    s/_(\d+)\(REVERSE\) / $1\(REVERSE\)/;
	    print DEST "$_\n";
	}
	close TMP;
	close DEST;
    };
    chdir $dir;

    mkdir "hmmModels/";
    my @models = glob "${models}*.fasta";
    if (scalar @models < 1) {
	die "Please provide models in alignment mode with option -m\n";
    }

    #Build paralog- and exon-specifc pHMM and perform homology search
    foreach my $alignments(@models){
	    	#print Dumper($alignments);
		$alignments =~ s/${models}//g;
		#print Dumper($alignments);
		$alignments =~ s/\.fasta//g;
		#print Dumper($alignments);
		`hmmbuild --amino ${dir}/hmmModels/${alignments}.model ${models}${alignments}.fasta`;
    };
    mkdir "SearchTarget/";  
    my @files = glob "hmmModels/*.model";
    #print "Progress of hmm search against the target genome:\n";
    foreach my $file(@files){
		$_ = $file;
		s/hmmModels\///;
		s/\.model//;
		$file = $_;
		#print "$file\n";
		system "hmmsearch --cpu $core --domtblout SearchTarget/${file}_domtableout -E 0.01 -o SearchTarget/${file}_out hmmModels/${file}.model ${space}${myOutFolder}_ORF.aas >> SearchTarget/Target_vs_query.hmmsearchout";
	};

    #print "Hmmsearch of query models against target genome done\n";

    ## cat all output hits to one file all_exons_out_$spec, format this file to be a list with the hits only, substitute with spaces and - with tabs, add columns with paralogID and exonID, save a copy in the file data in ExonMatchSolver/ 
    $hmm = "SearchTarget/hmmsearch_tableout";
    system "cat SearchTarget/*_domtableout| sed \'s/-\\([0-9]\\+\\)]/_\\1]/g\' > $hmm";   #anders
    system ("cat SearchTarget/*_out > SearchTarget/hmmsearch.out");			# warum? wird nie verwendet
    ($hmm, $cutpoints) = CompartGenomeHMM($hmm);
    #print "Dollar hmm is $hmm\n";
    #print Dumper("cutpoints");
    #print Dumper \%$cutpoints;

	if (-e $hmm){
	    open (IN, "< $hmm") or die "Couldn't open $hmm, $!";
	    open (HITS, "> ExonMatchSolver1.input");
	    print HITS "#TargetName\tStrandinformation\tModel\tE-Value\tBitscore\tBias\tParalogID\tExonID\n";
	    while (<IN>){
	    	#print Dumper ($_);
		unless ($_ =~ /#/) {
		    chomp;
		    my @hits = split /[\t\s]+/, $_;   				#anders
		    my @name = split /_/, $hits[3]; 
		    #print Dumper($hits[0]);
		    $hits[30] = $name[0];
		    #print Dumper($hits[30]);
		    #$hits[30] =~ s/Paralog//g;
		    $hits[31] = $name[1];	
		    #print Dumper($hits[31]);						#anders
		    $hits[31] =~ s/exon//g;
		    #print Dumper($hits[31]);
		    print HITS "$hits[0]\t$hits[22]\t$hits[3]\t$hits[6]\t$hits[7]\t$hits[8]\t$hits[30]\t$hits[31]\n";				#anders
		}
	    }
	    close HITS;
	    close IN;
	}
}


elsif($mode eq "fasta" || $mode eq "user"){
    #Run blast with the single TCEs
    #Reformat the blast-output to fit the ExonMatchSolver's needs, Blast-Hits outputted to Target_vs_query.blastout  
    mkdir "SearchTarget/";  
    #V3.2 modification: Stricter Evalue in Blast and hmmsearch filtering -e 0.001
    #print "$genome_target\t$blast\n";
    my @buildtarget = "formatdb -i ${genome_target} -p F -o";
    #print "@buildtarget\n";
    #print "$genome_target\n";
    unless (-e "${genome_target}.nhr" && -e "${genome_target}.nin" && -e "${genome_target}.nsd" && -e "${genome_target}.nsi" && -e "${genome_target}.nsq"){
	system (@buildtarget) == 0 or die "Formatdb with target genome failed. $!. You need a permission to write files in the target genome directory\n";
    };

    my $blast = "SearchTarget/Target_vs_query.blastout";
    my @call = "blastall -p tblastn -d $genome_target -i HomologousExons/paralog_exon_alignment_final.fa -e 0.01 -a $core -m8 -G 11 -E 1 -C F > $blast";
    system (@call) == 0 or die "blastall against target genome failed. $!\n";
    #print "Blast of query sequences against target genome done.\n";

    ($blast, $cutpoints) = CompartGenome($blast);

    @buildtarget = "formatdb -o -i $genome_target -p F";
    unless (-e "${genome_target}.nhr" && -e "${genome_target}.nin" && -e "${genome_target}.nsd" && -e "${genome_target}.nsi" && -e "${genome_target}.nsq"){
	system (@buildtarget) == 0 or die "Formatdb with target genome failed. $!\n";
    };
	print Dumper ("perl $ExceS_A_1 -genome ${genome_target} -blastout $blast -mode target");
	qx(perl $ExceS_A_1 -genome ${genome_target} -blastout $blast -mode target);
    
};
################Run ExonMatchSolver first time: purpose is to extract those exons that have not been scored in the first run
################Invoke ExonMatchSolver a first time with parameters given, print output and invokation options to ExonMatchSolver.out
{ 
	if (-e "$dir/ExonMatchSolver1.input"){
	    chdir "$Bin/" or die "Couldn't access $Bin/, $!";
	    my @exonSolver = "java -jar ExonMatchSolver5.1-12.10.jar -i $dir/ExonMatchSolver1.input -p $maxparalogs -u -t 1200 -c 40 >> $dir/ExonMatchSolver1.out"; 
	    # my @exonSolver = "java -jar ExonMatchSolver4.jar -i $dir/ExonMatchSolver1.input -p $maxparalogs $intervall $max_solution -u >> $dir/ExonMatchSolver1.out";
	    #print Dumper (@exonSolver);
	    if (-e "$dir/ExonMatchSolver1.out"){
		unlink "$dir/ExonMatchSolver1.out";
	    };
	    open (SOLUTION, ">> $dir/ExonMatchSolver1.out") or die "Can't open ExonMatchSolver1.out, $!";
	    print SOLUTION "ExonMatchSolver invoked with @exonSolver\n";
	    system (@exonSolver) == 0 or die "@exonSolver failed: $!";
	}
	else{
		print "No ExonMatchSolver input exists\n";
		exit;
	}
}
#print "ExonMatchSolver run for the first time.\n";

chdir $dir;

my $paranum = getParaNum("ExonMatchSolver1.out");
print "NUMBER OF PARALOGS IDENTIFIED IN TARGET GENOME: $paranum\n";

#Read first solution from ExonMatchSolver; this delivers a list of unscored TCEs (score is too low to be found with standard blast options), read this list only, put contig names as keys into a hash and unscored paralog-TCE names into an array/ reference this array
open (SOLUTION, "< ExonMatchSolver1.out") or die "Can't open ExonMatchSolver1.out, $!";
while (<SOLUTION>){
    chomp;
    if ($_ =~ /\#Solution/){
	last;
    }
    else{
		unless ($_ =~ /\bjava\b/ or $_ =~ /\#Unscored/ or $_ =~ /^$/){    	
		    my @exon_order = split "\t", $_;
		    #print Dumper(@exon_order);
		    if ($mode eq "user" || $mode eq "fasta"){
		    	$missing_scores{$exon_order[0]} = [@exon_order[1..$#exon_order]];
		    }
		    elsif($mode eq "alignment"){
		    	my @exon_order_new;
		    	for(my $i = 1;$i <= $#exon_order;$i++) {
		    		${exon_order[$i]} =~ s/-/_exon/;
		    		if (-e "hmmModels/${exon_order[$i]}.model"){
		    			push (@exon_order_new, $exon_order[$i]);
		    			#${exon_order[$i]} =~ s/_exon/-/;
		    		}
					$missing_scores{$exon_order[0]} = [@exon_order_new[0..$#exon_order_new]];
		    	}
			}
		}
    }
}
close SOLUTION;

#print Dumper (\%missing_scores);

#Rerun the homology search with more sensitive settings to retrieve scores for exons that were not scored until now
unless (-d "TargetSequences/"){
    mkdir "TargetSequences/";
};

if ($mode eq "user" || $mode eq "fasta"){
    #Extract genomic sequences of the contigs, stored as keys in %missing_scores
    foreach my $keys (keys %missing_scores){
	unless (-s "TargetSequences/${keys}.fa"){
	    #print "Here: $genome_target\t$keys\n";
	    `fastacmd -d $genome_target -s $keys| sed \'s/>lcl|/>/g\' > TargetSequences/${keys}.fa`;
	}
    }

    ###Extract the single exon sequences in single files that are unscored until now
    #Run blast a second time to find the missing scores
    if (-e "SearchTarget/Target_vs_query_missing_unsorted.blastout"){
	unlink "SearchTarget/Target_vs_query_missing_unsorted.blastout";
    };
    my $pm = new Parallel::ForkManager($core);

    #print Dumper (\%missing_scores);
    foreach my $contig (keys %missing_scores){
	my @exons_1 = @{$missing_scores{$contig}};
	foreach my $element (@exons_1){
	    $pm->start and next;
	    #@numle = split "-", $element;
	    #$numle[0] =~ s/(\d)\w/$1/;

	    `grep -A1 \'>$element\\s\*\$\' HomologousExons/paralog_exon_alignment_final.fa > SingleExons/${element}.fa`;
	    if (-s "SingleExons/${element}.fa"){	
		`bl2seq -p tblastn -i SingleExons/${element}.fa -j TargetSequences/${contig}.fa -G 11 -E 1 -D 1 -e 1 -F F| grep -v \'#\'|awk -v OFS='\t' '{print \$2, \$9"_"\$10, \$1, \$11, \$12, \$3}' >> SearchTarget/Target_vs_query_missing_unsorted.blastout`;
		
	    #`grep -A1 \'>${numle[0]}_exon${numle[1]}\\s\*\$\' HomologousExons/paralog_exon_alignment_final.fa > SingleExons/${numle[0]}_exon${numle[1]}.fa`;
	    #if (-s "SingleExons/${numle[0]}_exon${numle[1]}.fa"){	
		#`bl2seq -p tblastn -i SingleExons/${numle[0]}_exon${numle[1]}.fa -j TargetSequences/${contig}.fa -G 11 -E 1 -D 1 -e 1 -F F| grep -v \'#\'|awk -v OFS='\t' '{print \$2, \$9"_"\$10, \$1, \$11, \$12, \$3}' >> SearchTarget/Target_vs_query_missing_unsorted.blastout`;
		#print "bl2seq -p tblastn -i SingleExons/${numle[0]}_exon${numle[1]}.fa -j TargetSequences/${contig}.fa -G 11 -E 1 -D 1 -e 1 -F F| grep -v \'#\'|awk -v OFS='\t' '{print \$2, \$9\"_\"\$10, \$1, \$11, \$12, \$3}' >> SearchTarget/Target_vs_query_missing_unsorted.blastout\n";
	    }
	    else{
		unlink "SingleExons/${element}.fa"
	    }	
	    $pm->finish
	};
    };
    $pm->wait_all_children;
}

elsif ($mode eq "alignment"){

    {
	# the translated sequences (fasta) are stored in a file so that one sequence including the header is printed on one line for further processing		 
	unless (-e "${space}${myOutFolder}_ORF_grep"){
	    my $seq = "";		 
	    open (GENOME, "< ${space}${myOutFolder}_ORF.aas") or die "Couldn't open ${space}${myOutFolder}_ORF.aas, $!";
	    open (AAS, "> ${space}${myOutFolder}_ORF_grep") or die "Couldn't open ${space}${myOutFolder}_ORF_grep, $!";
	    while (<GENOME>){
		chomp;
		if (/^>/){
		    print AAS "$seq\n";
		    print AAS "$_ ";
		    $seq = ""; next;
		}
		chomp;	
		$seq .= $_;
	    }
	    print AAS "$seq\n";
	}
    };

    #print Dumper (\%missing_scores);
    ##Extract AAS sequences of the contigs, stored as keys in %missing_scores
    my %FILE;

    foreach my $keys (keys %missing_scores){
	unless (-s "TargetSequences/${keys}.aas"){
	    open ($FILE{$keys}, ">> TargetSequences/${keys}.aas") or die;
	}
    };
    open (AAS, "< ${space}${myOutFolder}_ORF_grep") or die "Couldn't open ${space}${myOutFolder}_ORF_grep, $!";
    while (<AAS>){
	my $name;
	chomp;
	if ($_ =~ /^>(\S+)\s/){
	    if (exists $missing_scores{$1}){
		#print "\$1 starts: $1\n";
		$name = $1;
		s/]\s/]\n/g;
		print {$FILE{$name}} "$_\n";
	    }
	}
	#close ($FILE{$name});
    };

    if (-e "SearchTarget/Target_vs_query_missing.hmmout"){
	unlink "SearchTarget/Target_vs_query_missing.hmmout";
    };
    foreach my $contig (keys %missing_scores){
	my @exons_1 = @{$missing_scores{$contig}};
	foreach my $element (@exons_1){
		#print Dumper ($element);
	    #@numle = split "-", $element;
	    #print Dumper("numle".@numle);
	    #Check this again
	    #$numle[0] =~ s/(\d)\w/$1/;	
	    #print Dumper("contig".$contig);
	    system "${hmmer_path}hmmsearch --cpu $core -E 1 --domE 1 --domtblout SearchTarget/${element}_${contig}_domtableout hmmModels/${element}.model TargetSequences/${contig}.aas >> SearchTarget/Target_vs_query_missing.hmmout";
	}
    }
    #if (-e "SearchTarget/*_*_*_domtableout"){
    system "cat SearchTarget/*_*_domtableout |grep -v '#'| sed -e \'s/\\s\\+/\\t/g\'| sed \'s/-\\([0-9]\\+\\)]/_\\1]/g\'| cat - $hmm > SearchTarget/hmmsearch_tableout_missingscores";
    #print "cat SearchTarget/*_*_*_domtableout |grep -v '#'| sed -e \'s/\\s\\+/\\t/g\'| cat - $hmm > SearchTarget/hmmsearch_tableout_missingscores\n";
	#}

    open (IN, "< SearchTarget/hmmsearch_tableout_missingscores") or die "Couldn't open SearchTarget/hmmsearch_tableout_missingscores, $!";
    open (HITS, "> ExonMatchSolver_missing.input");
    #print HITS "#TargetName\tStrandinformation\tModel\tE-Value\tBitscore\tBias\tParalogID\tExonID\n";
    while (<IN>){
	unless ($_ =~ /#/) {
	    chomp;
	    my @hits = split /[\s\t]+/, $_;
	    my @name = split /_/, $hits[3]; 
	    #$hits[20] = $name[0];
	    #$hits[20] =~ s/Paralog//g;
	    #$hits[21] = $name[1];
	    $hits[30] = $name[0];
	    #$hits[30] =~ s/Paralog//g;
	    $hits[31] = $name[1];
	    $hits[31] =~ s/exon//g;
	    my @coordinates = split /_/, $hits[22];
	    print HITS "$hits[0]\t$coordinates[0]_$coordinates[1]\t$hits[3]\t$hits[6]\t$hits[7]\t$hits[8]\t$hits[30]\t$hits[31]\n";
	}
    }
    close HITS;
    close IN;

};

#Check if there were additional hits retrieved
my $NoMissing = "0";

if ($mode eq "fasta" || $mode eq "user"){
    if (-e "SearchTarget/Target_vs_query_missing_unsorted.blastout") {
	`sort -k1,2 SearchTarget/Target_vs_query_missing_unsorted.blastout |uniq > SearchTarget/Target_vs_query_missing.blastout`;
    }
    else{
	#print "No exons missing from query paralog(s). Your query paralog(s) were detected in full length.\n";
	$NoMissing = "1";
	#print "ExonMatchSolver is not run for a second time.\n";
	`cp ExonMatchSolver1.out ExonMatchSolver2.out`;
	`cp ExonMatchSolver1.input ExonMatchSolver2.input`;
    };
};

#If additional hits were retrieved in more sensitive homology search, run ExonMatchSolver a second time
unless ($NoMissing == "1"){

    if ($mode eq "fasta" || $mode eq "user") {
	open (HITSTOTAL, "> ExonMatchSolver_missing.input") or die "Couldn't open ExonMatchSolver_missing.input, $!";
	open (TOTAL, "< SearchTarget/Target_vs_query_missing.blastout") or die "Couldn't open SearchTarget/Target_vs_query_total.blastout, $!";

	while (<TOTAL>){
	    chomp;
	    my @hits = split '[\s\t]+', $_;
	    my @name = split /_/, $hits[2]; 
	    $name[1] =~ s/exon//;
	    print HITSTOTAL "$_\t$name[0]\t$name[1]\n";
	}
	close HITSTOTAL;
	close TOTAL;

    }
    elsif($mode eq "alignment"){
	open (TOTAL, "< SearchTarget/hmmsearch_tableout_missingscores") or die "Couldn't open SearchTarget/hmmsearch_tableout_missingscores, $!";			#why?
    };

    if ($mode eq "alignment"){
	open (INPUT, "> ExonMatchSolver2.input") or die "Couldn't open , ExonMatchSolver2.input $!";
	print INPUT "#TargetName\tStrandinformation\tModel\tE-Value\tBitscore\tBias\tParalogID\tExonID\n";
	system "sort -r ExonMatchSolver_missing.input |uniq >> ExonMatchSolver2.input";
    }
    elsif($mode eq "user"|| $mode eq "fasta"){
	`cat ExonMatchSolver1.input ExonMatchSolver_missing.input > ExonMatchSolver2.input`; 
    }


    # invoke ExonMatchSolver a second time with parameters given, print output and invokation options in ExonMatchSolver2.out
    {
	chdir "$Bin/" or die "Couldn't access dir $Bin, $!";
	my @exonSolver = "java -jar ExonMatchSolver5.1-12.10.jar -i $dir/ExonMatchSolver2.input -p $paranum -u -t 1200 -c 40 >> $dir/ExonMatchSolver2.out";
	#my @exonSolver = "java -jar ExonMatchSolver4.jar -i $dir/ExonMatchSolver2.input $paranum $intervall $max_solution >> $dir/ExonMatchSolver2.out"; 
	if (-e "$dir/ExonMatchSolver2.out"){
	    unlink "$dir/ExonMatchSolver2.out"
	};
	open (SOLUTIONTOTAL, ">> $dir/ExonMatchSolver2.out") or die "Can't open ExonMatchSolver2.out, $!";
	print SOLUTIONTOTAL "ExonMatchSolver invoked with @exonSolver\n";
	system (@exonSolver) == 0 or die "system @exonSolver failed: $!";
    }
    #print "ExonMatchSolver run for the second time.\n";
    chdir $dir;

} #Close loop $NoMissing


#################### Parse output of second ExonMatchSolver run to identify scaffolds/fragments, where the paralog might be annotated completely
#################### Number of paralogs encoded in the target genome is retrieved from the output of the ExonMatchSolver run 2
chdir $dir;
open (SOLUTIONTOTAL, "< ExonMatchSolver2.out") or die "Can't open ExonMatchSolver2.out, $!";

my @lost_exons;
my %complete_contigs;
my %all_contigs;


{			
    my $switch = 0;
    my @paralogs;
    $paralogs[0] = 0;
    # filter output from ExonmatchSolver for contigs that might have the whole paralog annotated in solution 1 (first filter for contigs with more than one exon, sort alphabetically, filter for contigs that have TCEs 1, 2, 3, 4 or 5 and at the same time one of the last three TCEs, names of these contigs are given in @complete_contigs) 
    #print "@sorted_exons\n";
    #print "(One) of best ILP solutions: Scaffold\tParalog type\tExons on scaffold\n";

    open (SOLUTIONTOTAL, "< ExonMatchSolver2.out") or die "Can't open ExonMatchSolver2.out, $!";	
    if ($paranum eq "1") {
		#print Dumper("1");
	    $switch = 1;    
	};

    while (<SOLUTIONTOTAL>){
	chomp;
	#print Dumper ($_);
	if ($_ =~ /#paralogs/){
		#print Dumper("2");
	    @paralogs = split ': ', $_;
	    if ($paralogs[1] +1 eq "$paranum"){
		$switch = 1;
	    };
	};
	if ($_ =~ /Solution 1:/ && $switch eq "1") {
		#print Dumper("3");
	    $switch = 2;
	};   
	if ($_ =~ /\#paralogs: $paranum/ && $switch eq "2"){
		#print Dumper("4");
	    last;
	} 
	elsif ($switch eq "2" && $_ =~ m/paralog[0-9]+/){
		#print Dumper("5");
	    $_ =~ s/\(\d+\.\d+\)//g;
	    my @exon_order = split "\t", $_;
	    my $number = $exon_order[1];
	    $number=~ s/paralog//;
	    $exon_order[1] = $number . "-" . $exon_order[2];
	    my @exon_numbers = split " ", $exon_order[3];
	    if (@exon_numbers > 2 && scalar @sorted_exons < 10) {
			my @sorted_array = (@exon_order[0,1], sort {$a <=> $b} @exon_numbers[0..$#exon_numbers]);
			$all_contigs{$sorted_array[0]} = [$sorted_array[1], $sorted_array[2]];
				#print Dumper (@sorted_array);
				#print Dumper (@sorted_exons);		
			if ($sorted_array[2] <= $sorted_exons[0] || $sorted_array[2] <= $sorted_exons[1] || $sorted_array[2] <= $sorted_exons[2]|| $sorted_array[2] <= $sorted_exons[3]){
			    if ($sorted_array[-1] >= $sorted_exons[-3] || $sorted_array[-1] >= $sorted_exons[-2] || $sorted_array[-1] >= $sorted_exons[-1]){
					$complete_contigs{$sorted_array[0]} = [$sorted_array[1], $sorted_array[2]];
					foreach my $letter(@sorted_exons){
					    my $switch = 0;	
					    foreach $exon(@sorted_array){
						if ($letter eq $exon){
						    $switch = 1;
						}
					    }
					    if ($switch == 0){						
						push @lost_exons, [$sorted_array[0],$sorted_array[1],$letter];
					    }
					}
			    }
			};
			#print "@sorted_array\n";
			#print "@sorted_exons\n";
	    }		
	    elsif(@exon_numbers > 2){
	    	my @sorted_array = (@exon_order[0,1], sort {$a <=> $b} @exon_numbers[0..$#exon_numbers]);
			$all_contigs{$sorted_array[0]} = [$sorted_array[1], $sorted_array[2]];
			#print Dumper (@sorted_array);
			#print Dumper (@sorted_exons);
			if ($sorted_array[2] eq $sorted_exons[0] || $sorted_array[2] eq $sorted_exons[1] || $sorted_array[2] eq $sorted_exons[2]|| $sorted_array[2] eq $sorted_exons[3]|| $sorted_array[2] eq $sorted_exons[4] || $sorted_array[2] eq $sorted_exons[5] || $sorted_array[2] eq $sorted_exons[6] || $sorted_array[2] eq $sorted_exons[7]){
			    if ($sorted_array[-1] eq $sorted_exons[-6] || $sorted_array[-1] eq $sorted_exons[-5] || $sorted_array[-1] eq $sorted_exons[-4] || $sorted_array[-1] eq $sorted_exons[-3] || $sorted_array[-1] eq $sorted_exons[-2] || $sorted_array[-1] eq $sorted_exons[-1]){
				$complete_contigs{$sorted_array[0]} = [$sorted_array[1], $sorted_array[2]];
				foreach my $letter(@sorted_exons){
				    my $switch = 0;	
				    foreach $exon(@sorted_array){
					if ($letter eq $exon){
					    $switch = 1;
					}
				    }
				    if ($switch == 0){						
					push @lost_exons, [$sorted_array[0],$sorted_array[1],$letter];
				    }
				}
			    }
			}
		}			
		elsif(@exon_numbers = 2 && scalar @sorted_exons < 5){
	    	my @sorted_array = (@exon_order[0,1], sort {$a <=> $b} @exon_numbers[0..$#exon_numbers]);
			$all_contigs{$sorted_array[0]} = [$sorted_array[1], $sorted_array[2]];
			#print Dumper (@sorted_array);
			#print Dumper (@sorted_exons);
			if ($sorted_array[2] eq $sorted_exons[0] || $sorted_array[2] eq $sorted_exons[1] || $sorted_array[2] eq $sorted_exons[2]|| $sorted_array[2] eq $sorted_exons[3]|| $sorted_array[2] eq $sorted_exons[4] || $sorted_array[2] eq $sorted_exons[5]){
			    if ($sorted_array[-1] eq $sorted_exons[-4] || $sorted_array[-1] eq $sorted_exons[-3] || $sorted_array[-1] eq $sorted_exons[-2] || $sorted_array[-1] eq $sorted_exons[-1]){
				$complete_contigs{$sorted_array[0]} = [$sorted_array[1], $sorted_array[2]];
				foreach my $letter(@sorted_exons){
				    my $switch = 0;	
				    foreach $exon(@sorted_array){
					if ($letter eq $exon){
					    $switch = 1;
					}
				    }
				    if ($switch == 0){						
					push @lost_exons, [$sorted_array[0],$sorted_array[1],$letter];
				    }
				}
			    }
			}
		};
	}
    }
    close SOLUTIONTOTAL;


    if (keys %complete_contigs) {
	#print "Genomic fragments for which annotation of TCEs is (nearly) complete for one paralog are:\nGenomic fragment => [TargetParalogNumber-QueryParalog, first TCE detected]\n";
	#print Dumper (\%complete_contigs);
    }
    else{
	#print "No (nearly) complete gene copies detected on one genomic fragment.\n";
    }
}


#print "Annotation of nearly complete hits starts.\n";

#print Dumper(%complete_contigs);
runSplicedAligner("ExonMatchSolver2.input", \%complete_contigs);

chdir $dir;
if (-e "SearchTarget/Exonerate_additional_hits.blastout"){
    unlink "SearchTarget/Exonerate_additional_hits.blastout";
}

#Identify TCEs that are additionally found by the exonerate annotation
#How? bl2seq of TCEs that are missing in the putative complete annotations against the Exonerate-protein annotation 
#Filter via e-value (lower than 0.01) -> Make this more flexible?
if (@lost_exons) {
	#print Dumper(@lost_exons);
    #print "The following TCEs are missing in the (nearly) complete genomic fragments as proposed by the ExonMatchSolver. Its amino acid sequences are blasted against the Exonerate annotation of the respective paralog on the respective genomic fragments:\n";
    foreach my $lines(0..$#lost_exons){
	my $full_name = ${lost_exons[$lines][1]};
	${lost_exons[$lines][1]} =~ s/[0-9]\-//; 
	`grep -A1 '>${lost_exons[$lines][1]}_exon${lost_exons[$lines][2]}\\s\*\$' HomologousExons/paralog_exon_alignment_final.fa > SingleExons/${lost_exons[$lines][1]}_exon${lost_exons[$lines][2]}.fa`;
	if (-s "SingleExons/${lost_exons[$lines][1]}_exon${lost_exons[$lines][2]}.fa"){
	    #print "Paralog $lost_exons[$lines][1]\texon $lost_exons[$lines][2]\n";
	    my @bl2seq = "bl2seq -p blastp -i SingleExons/${lost_exons[$lines][1]}_exon${lost_exons[$lines][2]}.fa -j SearchTarget/${full_name}_target.fa -G 11 -E 1 -D 1 -e 0.01 -F F| grep -v '#' >> SearchTarget/Exonerate_additional_hits.blastout";
	    #print "@bl2seq\n";
	    system (@bl2seq) == 0 or warn "Bl2seq against Target-Annotation failed. $!\n";
	}
	else{
	    #print "The amino acid sequence of paralog ${lost_exons[$lines][1]}, exon${lost_exons[$lines][2]} does not exist. This TCE might be missing in this paralog (deleted) or it is too derived to be recognized as a homologous exon.\n";
	}
    };
};

#Add these Hits to the ExonMatchSolver input, set score to 0 and add the hits for every paralog
if (-s "SearchTarget/Exonerate_additional_hits.blastout"){
    open (PRO, "> $dir/SearchTarget/Exonerate_additional_hits_ExonMatchSolver.input") or warn "Couldn't open Exonerate_additional_hits_ExonMatchSolver.input, $!";
    open (ADD, "< SearchTarget/Exonerate_additional_hits.blastout") or warn "Couldn't open Exonerate_additional_hits_ExonMatchSolver.blastout, $!";
    while (<ADD>){
        chomp $_;
	foreach my $lines(0..$#lost_exons){
	    #print "${lost_exons[$lines][1]}_exon${lost_exons[$lines][2]}\n";
	    #${lost_exons[$lines][1]} =~ s/(\d)\w/$1/;
	    if ($_ =~ /${lost_exons[$lines][1]}_exon${lost_exons[$lines][2]}/){
		#print "$_\n";
		my $additional_hits = `grep '${lost_exons[$lines][0]}\t' ExonMatchSolver2.input| sort -nk 8,8| head -n 1`;
		chomp $additional_hits;
		#print "$additional_hits\n";
		my @add = split "\t", $additional_hits;
		#$_ = $additional_hits;
		my @names = split "_", $add[2];
		$names[1] =~ s/exon//g;
		foreach my $ID(@IDs){
		    print PRO $add[0] ."\t0_0\t" . ${ID} . "_exon" . ${lost_exons[$lines][2]} . "\t0\t0\t0\t" . ${ID} . "\t" . ${lost_exons[$lines][2]} . "\n";
		}
	    };
	};	
    };
    close PRO;
    close ADD;
}

#Check if a third run of EMS is necessary
if (-s "$dir/SearchTarget/Exonerate_additional_hits_ExonMatchSolver.input"){
    #print "Additional TCEs found by Exonerate.\n";
    system ("cat ExonMatchSolver2.input SearchTarget/Exonerate_additional_hits_ExonMatchSolver.input > ExonMatchSolver3.input") == 0 or warn "Concatination of input from first ExonMatchSolver-run and additional hits did not work.\n";
    
    #system ("cat ExonMatchSolver2.input SearchTarget/Exonerate_additional_hits_ExonMatchSolver.input|sort -r > ExonMatchSolver3.input") == 0 or warn "Concatination of input from first ExonMatchSolver-run and additional hits did not work.\n";
    #`rm SearchTarget/Exonerate_additional_hits_ExonMatchSolver.input`;

    #invoke ExonMatchSolver a third time with parameters given (hit list contains hits from two blast runs and additional Exonerate hits), print output and invokation options to $solution

    chdir "$Bin/" or die "Couldn't access dir $Bin, $!";
    my @exonSolver = "java -jar ExonMatchSolver5.1-12.10.jar -i $dir/ExonMatchSolver3.input -u -p $maxparalogs -t 1200 -c 40 >> $dir/ExonMatchSolver3.out"; 
    #my @exonSolver = "java -jar ExonMatchSolver4.jar -i $dir/ExonMatchSolver3.input $paranum $intervall $max_solution >> $dir/ExonMatchSolver3.out"; 
    if (-e "$dir/ExonMatchSolver3.out"){
	unlink "$dir/ExonMatchSolver3.out";
    };
    open (SOLUTIONCOMPLETE, ">> $dir/ExonMatchSolver3.out") or die "Can't open $dir/ExonMatchSolver3.out, $!";
    print SOLUTIONCOMPLETE "ExonMatchSolver invoked with @exonSolver\n";
    system (@exonSolver) == 0 or warn "system @exonSolver failed: $!";
    #print "ExonMatchSolver run for the third time.\n";
    chdir $dir;
}
else{
    #print "No additional TCEs found by Exonerate. ExonMatchSolver is not run a third time.\n";
    `cp ExonMatchSolver2.out ExonMatchSolver3.out`;
    `cp ExonMatchSolver2.input ExonMatchSolver3.input`;
};
close SOLUTIONCOMPLETE;


###################### Post-processing starts here with the final ExonMatchSolver input 
######################
#check for two paralogs on the same scaffold/contig/chromosome, search for Hits that were found by the same TCE twice on the same contig and had an e-value < 0, discrimination with e-value is necessary, otherwise all low scoring hits are shown as well which are probably not real TCEs 
{
    my $check = `cat ExonMatchSolver3.input| sort -g -k4,4| grep \'e-\'| awk \'{print \$1, \$3, \$7}\'| sort -k1,1|uniq -d| sort -n -k1,1| wc -l`;
    chomp $check;
    if ($check == 0){
	#print "Check for two paralogs on one unit: negative.\n";
    }
    elsif ($check > 0){
	#print "WARNING: Check for two hits of the following exons on one unit: positive. This can indicate many things e.g. a change in gene structure, an exon duplication or a complete gene duplication on the way to degradation or can just be a background hit. It is up to you to find this out! More then one blast-hit on the same contig with the same TCE were observed for:\n";
	system "cat ExonMatchSolver3.input| sort -g -k4,4| grep \'e-\'| awk \'{print \$1, \$3}\'| sort -k1,1|uniq -d| sort -n -k1,1"; 
	#print "This could also suggest a split exon.\n";
    };
};

#length normalization of single TCEs on contigs: Are they reliable?
{
    my $switch = 0;
    #my @sorted_array;
    #my @sorted_exons;
    my @paralogs;
    $paralogs[0] = 0;		
    
    # filter output from ExonmatchSolver for contigs that have only one TCE matched, extract bitscore for these hits
    open (SOLUTION, "< ExonMatchSolver3.out") or die "Can't open $dir/ExonMatchSolver3.out, $!";
    while (<SOLUTION>){
	chomp;
	#print Dumper($_);
	#print "Paralog number of final EMS solution is $paranum\n";
	if ($paranum eq "1") {
		#print Dumper("1");
	    $switch = 1;    
	};
	if ($_ =~ /#paralogs/){
	    @paralogs = split ': ', $_;
	    if ($paralogs[1] +1 eq "$paranum"){
		$switch = 1;
		#print "$paralogs[1]\n";
		#print "switch \n";
	    };
	};
	if ($_ =~ /Solution 1:/ && $switch eq "1") {
				#print Dumper("2");
	    $switch = 2;
	};
	#if($_ =~ /Solution 2:/ && $switch eq "2"){
	#print "WARNING: Several different optimal solutions were obtained by the last run of ExonMatchSolver, check $dir/ExonMatchSolver3.out.\n";
	#};
	if ($_ =~ /\#paralogs: $paranum/ && $switch eq "2"){
				#print Dumper("3");
	    last;
	}
	elsif ($switch eq "2" && $_ =~ m/paralog\s/){
	    #print "$_\n";
	    $_ =~ s/\(\d+\.\d+\)//g;
	    my @exon_order = split "\t", $_;
	    my @name = split /[\s]+/, $exon_order[1];
	    $name[3] =~ s/\)//g;
	    $exon_order[1] = $name[1] . "-" . $name[3];
	    #if ($_ =~ /\s(\d)+\s$/){
	    #unless ($_ =~ /ExonMatchsolver/ || $_ =~ /-/ || $_ =~ /#/ || $_ =~ /^$/){	
	    #my @exon_order = split "\t", $_;
	    #$exon_order[1] =~ s/[a-zA-Z]//g;
	    #$exon_order[1] =~ s/ //g;
	    #$exon_order[1] =~ s/\)//g;
	    #$exon_order[1] =~ s/\(/./g;
	    if (@exon_order < 4) {
		#print "$_\n";
		$single_exons{$exon_order[0]} = [$exon_order[1], $exon_order[2]];
		#print "Singles: $exon_order[0]\t$exon_order[1]\t$exon_order[2]\n";
	    }
	    #	}
	}
    }
    close SOLUTION;

    #If single exons exist, check whether they are relaible (by bitscore filtering) and take care of their placement
    my @reliable;
    if (keys %single_exons) {    
	#print "Single exons:\n";
	#print Dumper \%single_exons;

	if (-e "Single_exons"){
	    unlink "Single_exons";
	};
	if (-e "tmp"){
	    unlink "tmp";
	};
	open (TMP, ">> tmp");
	print TMP "Reliable exons, normalized Bitscore > 1.5:\n";

	#print: I stopped here with paralog name substitution

	#fix this part here, paralog name is wrong here again
	open (RESULT, ">> Single_exons");
	print RESULT "Single exons on contigs:\n#Contig\tParalog\tExon\tBitScore\tNormalizedBitScore\n";
	foreach my $contig(keys %single_exons){
	    my @data = @{$single_exons{$contig}};
	    my $letter = $data[1];
	    #$data[0] =~ s/\d+\.//g;
	    my ($type, $para) = split "-", $data[0];
	    #my $para = $data[0];
	    #print "Para: $para\n";
	    #print "Exon: $letter\n";
	    my $length = ();
	    my $score = qx(grep -P '$contig\t' $dir/ExonMatchSolver3.input| grep -v '#'| sort -n -k 5,5| tail -n 1|awk '{print \$5}'); 
	    chomp $score;
	    if ($mode eq "user"|| $mode eq "fasta"){
		$length = qx(grep -A1 '${para}_exon${letter}\$' $dir/HomologousExons/paralog_exon_alignment_final.fa| tail -n1| sed ':a;N;\$!ba;s/\\n//g'|wc -c);
		chomp $length;
	    }
	    elsif($mode eq "alignment"){
		$length = qx(grep 'LENG' hmmModels/${para}_exon${letter}.model| awk '{print \$2}');	
	    };
	    if ($score =~ /\d/){	
		my $normalized_score = $score/$length;	
		print RESULT "$contig\t$para\t$letter\t$score\t$normalized_score\n";
		if ($normalized_score > 1.5){	
		    print TMP "$contig\t$para\t$letter\t$score\t$normalized_score\n";
		}
	    }
	    else{
	     	#print "WARNING: The extraction of bitscore and length of single Exons failed. grep -P '$contig\t' $dir/ExonMatchSolver3.input| grep -v '#'| grep '${para}_'| sort -n -k 5,5| tail -n 1|awk '{print \$5}'\n"
	    }
	}
	close TMP;
	close RESULT;

	my $reliable_exons = `cat tmp| grep -v '^\$'| grep -v '#'|grep -v 'exons'|sort -k 2,3| awk '{print \$2"\t"\$3}'| uniq -d| grep -v -f - tmp| grep -v 'exons'`;
	my $non_explicit_exons = `cat tmp| grep -v '^\$'| grep -v '#'|grep -v 'exons'|sort -k 2,3| awk '{print \$2"\t"\$3}'| uniq -d| grep -f - tmp`;
	open (RESULT_2, "> Reliable_exons");
	print RESULT_2 "#Contig\tParalog\tExon\tBitScore\tNormalizedBitScore\n";
	print RESULT_2 "Reliable, explicit matching exons:\n$reliable_exons";
	print RESULT_2 "Reliable, non-explicit matching exons:\n$non_explicit_exons";

	#If there are reliable single exons, those have to be inserted into the other contigs
	#{
	open (RESULT_2, "< Reliable_exons");
	while (<RESULT_2>){
	    chomp;
	    unless ($_ =~ /exons/ || $_ =~ /#/){
		my @array = split "\t", $_;
		push @reliable, [$array[0], $array[1], $array[2]];
	    }
	}
	#print "@reliable\n";
	#print "Reliability check of single TCEs completed.\n";
	`rm tmp`;
	close RESULT_2;
    }
    else{
	#print "No single exons exist in final solution.\n";
    }

    #Assemble the contigs with more than one hit toegther with the reliable single hits for the optimal solution 1, they are stored in the hash %complete (inside is another hash with TCEs as keys, contigs as values)
    #Unreliable single TCEs are not further considered in the automated pipeline
    my %complete;
    #print Dumper \%complete;
    {
	my $switch = "0";
	my @sorted_array;
	my @sorted_exons;
	my @paralogs;
	$paralogs[0] = 0;

	open (EMSOUT, "< ExonMatchSolver3.out");
	while (<EMSOUT>){
	    chomp;
	    if ($paranum eq "1") {
	    	$switch = 1;    
	    };
	    if ($_ =~ /#paralogs/){
	    	@paralogs = split ': ', $_;
	        if ($paralogs[1]+1 eq "$paranum"){
		    $switch = "1";
		    #print "$paralogs[1] Switch is 1, $paranum\n";
		    #print "$paralogs[1]\n";
		    #print "switch \n";
	        };
	    };
	    if ($_ =~ /Solution 1:/ && $switch eq "1") {
	        $switch = "2";
		#print "$switch $_\n";
	    }
	    elsif($_ =~ /Solution 2:/ && $switch eq "1"){
		#print "WARNING: Several different optimal solutions were obtained by the last run of ExonMatchSolver, check $dir/ExonMatchSolver3.out.\n";
	    };
	    if ($_ =~ /\#paralogs: $paranum/ && $switch eq "2"){
		last;
	    } 
	    elsif ($switch eq "2" && $_ =~ m/paralog[0-9]+/){
	
		$_ =~ s/\(\d+\.\d+\)//g;

	   	my @exon_order = split "\t", $_;
	    my $number = $exon_order[1];
	    $number=~ s/paralog//;
	    $exon_order[1] = $number . "-" . $exon_order[2];
	    my @exon_numbers = split " ", $exon_order[3];
	    if (@exon_numbers > 2) {
			my @sorted_array = (@exon_order[0,1], sort {$a <=> $b} @exon_numbers[0..$#exon_numbers]);
		    #print "@sorted_array\n";
		    foreach my $i(0..$#reliable){		
			if ($sorted_array[1] =~ $reliable[$i][1]){
			    $complete{$exon_order[1]}{$reliable[$i][2]}=$reliable[$i][0];
			}	
		    };
		    foreach my $j(@sorted_array[2..$#sorted_array]){
			$complete{$exon_order[1]}{$j}=$sorted_array[0];
		    }
		}
		#}
	    }
	};
	seek EMSOUT, 0, 0;
    }


    #print "Solution of the last run of ExonMatchSolver including contigs with single, reliable TCEs\n";
    #print "TargetParalogNumber-QueryParalog => Exon => Scaffold\n";
   # print Dumper (\%complete);
    close EMSOUT;

    # Check for every paralog if the previous contig/scaffold is the same as the current one. If not, keep these entries in %merge (array in this array keeps paralog, exon, contig)
    my %merge; 
    my $i=0;
    $switch = 0;
    for my $paralog (keys %complete){
		my $previous_paralog="";
		my $previous_exon ="";	
		for my $exon (sort {$a <=> $b} keys %{$complete{$paralog}}){
			#print Dumper($exon, $complete{$paralog}{$exon}, $previous_exon);
		    if ($complete{$paralog}{$exon} ne $previous_exon) {
				if ($paralog ne '' && $previous_paralog ne '' && $previous_exon ne ''){
				    if (exists $merge{$paralog}{$previous_exon}){
					for my $i (0..$#{$merge{$paralog}{$previous_exon}}){
					    if ($previous_paralog eq $merge{$paralog}{$previous_exon}[$i]){
						$switch =1;					
					    }			
					}
					if ($switch != 1){						
					    push (@{$merge{$paralog}{$previous_exon}}, $previous_paralog);
					}
					$switch = 0;				
				    }
				    elsif(!exists $merge{$paralog}{$previous_exon}){
					$merge{$paralog}{$previous_exon} = [$previous_paralog,];
				   }		
				};
				if ($paralog ne '' && $exon ne '' && $complete{$paralog}{$exon} ne ''){
					#print Dumper($merge{$paralog}{$complete{$paralog}{$exon}});
				    if (exists $merge{$paralog}{$complete{$paralog}{$exon}}){
						for my $j (0..$#{$merge{$paralog}{$complete{$paralog}{$exon}}}){
						    if ($exon eq $merge{$paralog}{$complete{$paralog}{$exon}}[$j]){
								$switch =1;					
						    }			
						}
						if ($switch != 1){
						    push (@{$merge{$paralog}{$complete{$paralog}{$exon}}},$exon);
						}
						$switch = 0;
					}
					elsif(!exists $merge{$paralog}{$complete{$paralog}{$exon}}){
						$merge{$paralog}{$complete{$paralog}{$exon}} = [$exon,];			
				    }				
				};
				$previous_exon=	$complete{$paralog}{$exon};
				$previous_paralog= $exon;
			    }
			else{
				$previous_exon=$complete{$paralog}{$exon};
				$previous_paralog= $exon
		    };
		}
    }
    #print Dumper (\%merge);
    my $exonCounter =0;
    my %split;
    my %cat;
    my %test;

    foreach my $para(keys %merge){
	foreach my $contig(keys %{$merge{$para}}){
		#print Dumper($contig);
	    unless (-s "TargetSequences/${contig}.fa"){
		`fastacmd -d $genome_target -s $contig| sed \'s/>lcl|/>/g\' > TargetSequences/${contig}.fa`;
		#print "fastacmd -d $genome_target -s $contig| sed \'s/>lcl|/>/g\' > TargetSequences/${contig}.fa";
	    }
	}
    }
	#print Dumper (%merge);
    for my $scaff(keys %merge){
	my @moreExons;
	#simplest case: only one scaffold
	if (scalar (keys %{$merge{$scaff}}) != 1){
	    #for the cases with more than one scaffold
	    my @twoElements = ();
	    my @keys =();
	    my @exons1 = ();
	    my @sorted_exons1 = ();
	    my $switch = 0;
	    #for my $element (sort {$merge{$scaff}{$a} <=> $merge{$scaff}{$b}} keys %{$merge{$scaff}}){
	    #push @twoElements, $scaff;
	    #};
	    #my @uniqTwoElements = uniq (@twoElements);
	    #for my $scaffold(@uniqTwoElements){
	    for my $element (sort {$merge{$scaff}{$a}[0] <=> $merge{$scaff}{$b}[0]} keys %{$merge{$scaff}}){
		foreach my $j(0..$#{$merge{$scaff}{$element}}){		
		    push @exons1, $merge{$scaff}{$element}[$j]; 
		};
	    };
	    #print "@exons1\n";
	    @sorted_exons1 = sort {$a <=> $b} @exons1;
	    #print "SORTED_EXONS: @sorted_exons1\t@exons1\n";
	    foreach my $i(0..$#exons1){	
		if($sorted_exons1[$i] != $exons1[$i]){
		    $switch = 1;
		}
	    }	
	    if ($switch == 1){
		for my $element (keys %{$merge{$scaff}}){
		    foreach my $i (0..$#{$merge{$scaff}{$element}}){
			$split{$scaff}{$element}{$merge{$scaff}{$element}[$i]} = [],;
		    };
		}
	    }
	    elsif ($switch == 0){
		for my $element (keys %{$merge{$scaff}}){
		    foreach my $i (0..$#{$merge{$scaff}{$element}}){
			$cat{$scaff}{${merge{$scaff}{$element}[0]}} = [$element,];
		    }
		}
	    }
	    #}
	    delete $merge{$scaff};	
	}
	elsif(scalar (keys %{$merge{$scaff}}) == 1){
	    for my $element (keys %{$merge{$scaff}}){
		foreach my $contig(keys %complete_contigs){
		    #print "TEST: $contig\t$element\t$scaff\n"; 
		    if ($element eq $contig){
			delete $merge{$scaff};
		    }
		}
	    }	
	}
    }

    foreach my $para(keys %merge){
	open (OUTPROSPLIGN, "> SearchTarget/${para}_target.fa");

		foreach my $contig(keys %{$merge{$para}}){
		    my $querypara = (split "-", $para)[1];
		    if ($mode eq "fasta"|| $mode eq "user"){
				my $strand_last ="";

				my $in = "ExonMatchSolver3.input";
				my @output = ExtractFirstLast($in, $contig);
				$strand = $output[0];

				my @first;
				my @last;
				$first[0] = $output[1];
				$first[1] = $output[2];
				$last[0] = $output[3];
				$last[1] = $output[4];
				
				if ($strand eq "+"){
				    ${merge{$para}{$contig}[1]}= $first[0] - $distance_first;
				    ${merge{$para}{$contig}[2]} = $last[1]+ $distance_last;		
				}
				elsif($strand eq "-"){
				    ${merge{$para}{$contig}[2]} = $first[1] + $distance_first;
				    ${merge{$para}{$contig}[1]} = $last[0] - $distance_last;
				};
		    }
	    	elsif($mode eq "alignment"){

				my $in = "ExonMatchSolver3.input";
				#print Dumper($in, $contig);
				my @output = ExtractFirstLast($in, $contig);
				#print Dumper("4");
				#print Dumper(@output);
				if (! @output || ! $output[1]){
					next;
				}

				$strand = $output[0];
				#print "$input, $contig, @output\n";
				my @first;
				my @last;
				$first[0] = $output[1];
				$first[1] = $output[2];
				$last[0] = $output[3];
				$last[1] = $output[4];		
				
				if ($strand eq "+"){
				    ${complete_contigs{$contig}[2]}= $first[0] - $distance_first;
				    ${complete_contigs{$contig}[3]} = $last[1]+ $distance_last;		
				}
				elsif($strand eq "-"){
				    ${complete_contigs{$contig}[3]} = $first[1] + $distance_first;
				    ${complete_contigs{$contig}[2]} = $last[0] - $distance_last;
				};
			};
			my $para_original = $para;
		    $_ = $para_original;			
		    s/[0-9]\.//;
		    $querypara = $_;
		    my $querypara_old = (split "-", $para_original)[1];

			if ($mode eq "fasta"|| $mode eq "user"){
				print Dumper ("perl $ExceS_A_2 --contig $contig --query $para_original --strand $strand --start ${merge{$para}{$contig}[1]}--end ${merge{$para}{$contig}[2]} --mode 2");
				system("perl $ExceS_A_2 --contig $contig --query $para_original --strand $strand --start ${merge{$para}{$contig}[1]} --end ${merge{$para}{$contig}[2]} --mode 2");
			}


			elsif($mode eq "alignment"){
		   	 	print Dumper ("perl $ExceS_A_2 --contig $contig --query $para_original --strand $strand --start ${complete_contigs{$contig}[2]} --end ${complete_contigs{$contig}[3]} --mode 2");
				system("perl $ExceS_A_2 --contig $contig --query $para_original --strand $strand --start ${complete_contigs{$contig}[2]} --end ${complete_contigs{$contig}[3]} --mode 2");
			}
	    	#print Dumper ("perl $ExceS_A_2 --target TargetSequences/${contig}.fa --query SearchQuery/${querypara_old}.fa --strand $strand --paralog $para_original --mode 4");
			#system ("perl $ExceS_A_2 --target TargetSequences/${contig}.fa --query SearchQuery/${querypara_old}.fa --strand $strand --paralog $para_original --mode 4");

		}
    }

    #For the case %cat
    #Fill %cat with the following information: paralog => TCE => [contig, start, stop, strand, seq], start and stop are extracted from the first TCE on this contig
    #seq is always extracted so that it is on the positive strand
    #sort TCEs and cat sequences of the corresponding contigs to eachother
    {
	for my $para(keys %cat){
		#print Dumper($para);
	    my $para_original = $para;
	    #$para =~ s/(\d)\w/$1/;
	    #$_ = $para;			
	    #s/[0-9]\.//;
	    my $querypara = (split "-", $para)[1];
	    #my $querypara = $_;
	    #print "$para_original\t$querypara\n";

	    open (FINALCAT, "> TargetSequences/${para_original}_merged.fa");
	    open (OUTCAT, ">  SearchTarget/${para_original}_target_final.fa");	

	    my $fastacat = "";
	    #print Dumper (%{$cat{$para_original}});
	    for my $exon (sort {$a <=> $b} keys %{$cat{$para_original}}){
		my $seq = "";
		#$para_original = $para;
		my $querypara = (split "-", $para_original)[1];
		#$para =~ s/(\d)\w/$1/;	
		my $loci = `grep -P '$cat{$para_original}{$exon}[0].*_exon${exon}\t' ExonMatchSolver3.input| sort -g -k4,4| head -n1| awk '{print \$2}'`;
		#print "grep -P '$cat{$para_original}{$exon}[0].*_exon${exon}\t' ExonMatchSolver3.input| sort -g -k4,4| head -n1| awk '{print \$2}'\n";
		chomp $loci;
		#print Dumper ($loci);
		if ($mode eq "fasta"|| $mode eq "user"){		
		    my @temp = split (/_/, $loci);
		    foreach my $i (0..$#temp){
			push @{$cat{$para_original}{$exon}}, $temp[$i];
		    };				
		}
		elsif ($mode eq "alignment"){
		    my @temp = split (/_/, $loci);
		    $temp[0] =~ s/^.+\[//g;
		    $temp[1] =~ s/]//g;
		    #print Dumper ($temp[0]);print Dumper ($temp[1]);
		    foreach my $i (0..$#temp){
			push @{$cat{$para_original}{$exon}}, $temp[$i];
		    };
		}			
		if ($cat{$para_original}{$exon}[1] > $cat{$para_original}{$exon}[2]){
		    $cat{$para_original}{$exon}[3] = "-";
		}
		elsif($cat{$para_original}{$exon}[1] < $cat{$para_original}{$exon}[2]){
		    $cat{$para_original}{$exon}[3] = "+";
		};
		#}		
		#print Dumper ($cat{$para_original}{$exon});
		my $s = (split "-", $cat{$para_original}{$exon}[0])[0];
		if ($cat{$para_original}{$exon}[3] eq "-"){
		    $seq = `fastacmd -d $genome_target -s $s -S 2 | grep -v '>'|sed ':a;N;\$!ba;s/\\n//g'`
		}
		elsif($cat{$para_original}{$exon}[3] eq "+"){
		    $seq = `fastacmd -d $genome_target -s $s| grep -v '>'|sed ':a;N;\$!ba;s/\\n//g'`;

		};
		chomp $seq;
		#print "$seq\n";
		$cat{$para_original}{$exon}[4] = $seq;
		$fastacat .= $cat{$para_original}{$exon}[4];
	    };
	    $_ = $fastacat;
	    s/.{32000}/$&\n/g;
	    $fastacat = $_;
	    print FINALCAT ">${para_original}\n$fastacat\n";
	    close FINALCAT;

	    my $para_original = $para;
	    $_ = $para_original;			
	    s/[0-9]\.//;
	    $querypara = $_;
	    my $querypara_old = (split "-", $para_original)[1];
		
	    print Dumper ("perl $ExceS_A_2 --target TargetSequences/${para_original}_merged.fa --query ${querypara_old} --strand + --paralog $para_original --mode 4");
		system ("perl $ExceS_A_2 --target TargetSequences/${para_original}_merged.fa --query ${querypara_old} --strand + --paralog $para_original --mode 4");
    }
	}

    #For the case %split: 

    my $strand_overall_positive;
    my $strand_overall_negative;
    for my $para(keys %split){
    	#print Dumper (%split);
    	#print Dumper ($para);
	my $para_original = $para;
	my $querypara = (split "-", $para)[1];
	my $fastasplit = "";
	my $difference = 0;
	my $fastacat = "";	
	open (FINAL, "> TargetSequences/${para_original}_merged.fa");
	#if ($exonerate eq "no"){	
	open (OUT, ">  SearchTarget/${para_original}_target_final.fa");	
	#};
	$strand_overall_positive= 0;
	$strand_overall_negative= 0;
	my @exons = ();
	my @sorted_exons = ();
	#Extract maxcontig, the contig with the highest numbers of exons. In case two contigs have the same number of exons, the first one is kept. The strand of this contig is used as reference for the other contigs of this paralog. Maxcontig is also used to extract the locus where the single contigs should be inserted in maxcontig (triples). 
	my $maxsize =0;
	my $maxcontig =0; 
	my $maxcontig2 =0;
	my $switchMinus =0;
	my $switchPlus =0;		
	my $seq = ();
	my $majorstrand = 0;	
	my $scrambled = 0;
	my @counter = ();
	my %maxcontig;
	my @maxsorted_exons =();
	my %exons_minor;
	for my $contig2(keys %{$complete{$para_original}}){
	    $maxcontig{$complete{$para_original}{$contig2}}++;
	}
	#print Dumper (\%maxcontig);
	for my $contig0(keys %{$split{$para_original}}){
	    #print "$contig2\t$maxcontig{$contig2}\n"; 	
	    push @counter, $maxcontig{$contig0};
	}	
	my @sorted_counter = sort {$a <=> $b}@counter;		
	#print "$sorted_counter[-1]\n";
	my %rmaxcontig = reverse %maxcontig;		
	$maxcontig = $rmaxcontig{$sorted_counter[-1]};
	#print "MAX: $maxcontig\n";

	for my $pa(keys %split){
	    for my $sca (keys %{$split{$pa}}){
		if ($sca ne $maxcontig) {
		    for my $parax(keys %complete){
			for my $exon(keys %{$complete{$parax}}){
			    if ($complete{$parax}{$exon} eq $sca) {
				if (!exists $split{$parax}{$sca}{$exon}) {
				    $split{$parax}{$sca}{$exon} = [],;
				}
			    }
			}
		    }
		}
	    }
	}

	#Extract the strand information and loci of blast hits for the exons in %split, at the same time count the strand for all exons to infer whether they are all on the negative strand ($strand_overall_positive = 0) or on the postive strand ($strand_overall_negative = 0)
	for my $contig (keys %{$split{$para_original}}){
	    for my $exon (keys %{$split{$para_original}{$contig}}){
		push @exons, $exon;
		@sorted_exons = sort {$a <=> $b} @exons;	
		#print "@sorted_exons\n";
		#print "$para_original\t$para\n";		
		my $loci = `grep -P '$contig.*_exon$exon\t' ExonMatchSolver3.input| sort -g -k4,4| head -n1| awk '{print \$2}'`;
		chomp $loci;
		#print "$loci\n";
		if ($mode eq "user"|| $mode eq "fasta"){
		    my @temp = split (/_/, $loci);
		    foreach my $i (0..$#temp){
			$split{$para_original}{$contig}{$exon}[$i] = $temp[$i];
		    };				
		}elsif ($mode eq "alignment"){
		    my @temp = split (/_/, $loci);
		    $temp[0] =~ s/^.+\[//g;
		    $temp[1] =~ s/]//g;
		    foreach my $i (0..$#temp){
			$split{$para_original}{$contig}{$exon}[$i] = $temp[$i];
		    };					
		}
		if ($split{$para_original}{$contig}{$exon}[0] > $split{$para_original}{$contig}{$exon}[1]){
		    $split{$para_original}{$contig}{$exon}[2] = "-";
		    $strand_overall_negative ++;
		}elsif($split{$para_original}{$contig}{$exon}[0] < $split{$para_original}{$contig}{$exon}[1]){
		    $split{$para_original}{$contig}{$exon}[2] = "+";
		    $strand_overall_positive ++;			
		};
	    }
	}
	#print "$para: $strand_overall_positive\n";
	#print "$para: $strand_overall_negative\n";

	#For maxcontig, put all its exons in an array and sort it, extract the strand-information for maxcontig
	#For the other contigs, extract the genomic sequence of the blast hit + 10 nt upstream and downstream, store it in %split in the strand orientation of maxcontig, if $strand_overall_positive or negative are 0, don't change the strands at all 
	for my $contig3(keys %{$split{$para_original}}){
		#print Dumper ($contig3);
	    if($contig3 eq $maxcontig){
		    @exons = ();
		    for my $exon3(sort {$a <=> $b} keys %{$split{$para_original}{$maxcontig}}){	
				push @exons, $exon3;
		    };
		    #print Dumper(@exons);
		    @maxsorted_exons = sort {$a <=> $b} @exons;
		    #print "2: @maxsorted_exons\n";
		    for my $j (0..$#maxsorted_exons){
		    	#print Dumper ($j);
				if (exists $split{$para_original}{$maxcontig}{$maxsorted_exons[$j]}[2]){
				    if ($split{$para_original}{$maxcontig}{$maxsorted_exons[$j]}[2] eq "-"){
					$switchMinus ++;
				    }
				    elsif ($split{$para_original}{$maxcontig}{$maxsorted_exons[$j]}[2] eq "+"){
					$switchPlus ++;					
				    }
				}
		    };
		    #print Dumper ($switchMinus, $switchPlus);
		    if ($switchMinus >= 1 && $switchPlus >= 1){
			#print "WARNING: There are hits on the positive and negative strand of the same contig! This could be a scrambled locus. Please check: $maxcontig\n";			
		    };
		    if ($switchMinus > $switchPlus){
			$majorstrand = "-";				
		    }
		    elsif ($switchMinus < $switchPlus){
			$majorstrand = "+";				
		    }
		    else{
		    	#print Dumper ("mäh");
		    	@exons = ();
		    	@maxsorted_exons = ();
		    	$switchMinus = 0;
		    	$switchPlus = 0;
		    	for my $exon6(keys %maxcontig){
		    		if ($exon6 ne $maxcontig && $maxcontig{$exon6} == $sorted_counter[-1]){
				    	#print Dumper ($exon6);
					    for my $exon5(sort {$a <=> $b} keys %{$split{$para_original}{$exon6}}){	
							push @exons, $exon5;
					    };
					    #print Dumper(@exons);
					    @maxsorted_exons = sort {$a <=> $b} @exons;
					    #print "2: @maxsorted_exons\n";
					    for my $j (0..$#maxsorted_exons){
					    	#print Dumper ($j);
							if (exists $split{$para_original}{$exon6}{$maxsorted_exons[$j]}[2]){
							    if ($split{$para_original}{$exon6}{$maxsorted_exons[$j]}[2] eq "-"){
								$switchMinus ++;
							    }
							    elsif ($split{$para_original}{$exon6}{$maxsorted_exons[$j]}[2] eq "+"){
								$switchPlus ++;					
							    }
							}
					    };	    
					    if ($switchMinus >= 1 && $switchPlus >= 1){
							#print "WARNING: There are hits on the positive and negative strand of the same contig! This could be a scrambled locus. Please check: $maxcontig\n";	
						}		
						if ($switchMinus > $switchPlus){
							$majorstrand = "-";				
						}
						elsif ($switchMinus < $switchPlus){
							$majorstrand = "+";				
						}
					}
					$maxcontig = $exon6;
				}
			}
		    #print "MAJORSTRAND: $majorstrand of $para_original\n";
		    $switchMinus = 0;
		    $switchPlus = 0;
		    for my $j (0..$#maxsorted_exons-1){
				if ($majorstrand eq "+"){
				    if ($split{$para_original}{$maxcontig}{$maxsorted_exons[$j]}[0] > $split{$para_original}{$maxcontig}{$maxsorted_exons[$j+1]}[0]){			
					$scrambled = 1;						
				    };
				}
				elsif ($majorstrand eq "-"){
				    if ($split{$para_original}{$maxcontig}{$maxsorted_exons[$j]}[0] < $split{$para_original}{$maxcontig}{$maxsorted_exons[$j+1]}[0]){	
					$scrambled = 1;
				    };
				};
		    };
		    if ($scrambled == 1){
				#print "WARNING: There are hits in a scrambled order on the same contig! Please check: $maxcontig\n";				
				$scrambled = 0;	
			}		
		};

	    if ($contig3 ne $maxcontig){
			#print "CONTROL: $contig3\n";
			for my $exon3(sort {$a <=> $b} keys %{$split{$para_original}{$contig3}}){
			    $exons_minor{$exon3} = $contig3;
			    my $start = 0;
			    my $end = 0;
			    #print "$para\t$contig3\t$exon3\n";
			    #my $length = qx(grep -v \'>\' TargetSequences/${contig3}.fa |sed ':a;N;\$!ba;s/\\n//g'|wc -m);
			    if ($strand_overall_positive == 0 || $strand_overall_negative == 0){	
					if ($split{$para_original}{$contig3}{$exon3}[2] eq "-"){
				    	if ($mode eq "fasta"|| $mode eq "user"){
							$start = $split{$para_original}{$contig3}{$exon3}[1] - 30;
							$end = $split{$para_original}{$contig3}{$exon3}[0] + 30;
				    	}
				    	elsif ($mode eq "alignment"){
							$start = $split{$para_original}{$contig3}{$exon3}[1] - 30;
							$end = $split{$para_original}{$contig3}{$exon3}[0] + 30;
				    	};
				    	#print Dumper ("12");
				    	my @buildtarget = "formatdb -o -i TargetSequences/${contig3}.fa -p F";
    					unless (-e "TargetSequences/${contig3}.nhr" && -e "TargetSequences/${contig3}.nin" && -e "TargetSequences/${contig3}.nsd" && -e "TargetSequences/${contig3}.nsi" && -e "TargetSequences/${contig3}.nsq"){
						system (@buildtarget) == 0 or die "Formatdb with target genome failed. $!\n";
						}
				    	#print "fastacmd -d TargetSequences/${contig3}.fa -s $contig3 -L $start,$end| grep -v '>'|sed ':a;N;\$!ba;s/\\n//g'\n";					
				    	$seq = `fastacmd -d TargetSequences/${contig3}.fa -s $contig3 -L $start,$end 2> error| grep -v '>'|sed ':a;N;\$!ba;s/\\n//g'`;
				    	if ($seq eq ''){
				    		#print Dumper ("11");
				    		my @buildtarget = "formatdb -o -i TargetSequences/${contig3}.fa -p F";
    						unless (-e "TargetSequences/${contig3}.nhr" && -e "TargetSequences/${contig3}.nin" && -e "TargetSequences/${contig3}.nsd" && -e "TargetSequences/${contig3}.nsi" && -e "TargetSequences/${contig3}.nsq"){
							system (@buildtarget) == 0 or die "Formatdb with target genome failed. $!\n";
							}
							$seq = `fastacmd -d TargetSequences/${contig3}.fa -s $contig3 -L 0,0 | grep -v '>'|sed ':a;N;\$!ba;s/\\n//g'`;					
				    	}		
				    	chomp $seq;
					}
					elsif ($split{$para_original}{$contig3}{$exon3}[2] eq "+"){
				    	if ($mode eq "fasta"|| $mode eq "user"){
							$start = $split{$para_original}{$contig3}{$exon3}[0] - 30;
							$end = $split{$para_original}{$contig3}{$exon3}[1] + 30;
				    	}
				    	elsif ($mode eq "alignment"){
							$start = $split{$para_original}{$contig3}{$exon3}[0] - 30;
							$end = $split{$para_original}{$contig3}{$exon3}[1] + 30;
				    	}
				    	#print Dumper ("10");
				    	my @buildtarget = "formatdb -o -i TargetSequences/${contig3}.fa -p F";
    					unless (-e "TargetSequences/${contig3}.nhr" && -e "TargetSequences/${contig3}.nin" && -e "TargetSequences/${contig3}.nsd" && -e "TargetSequences/${contig3}.nsi" && -e "TargetSequences/${contig3}.nsq"){
						system (@buildtarget) == 0 or die "Formatdb with target genome failed. $!\n";
						}
				    	#print "fastacmd -d TargetSequences/${contig3}.fa -s $contig3 -L $start,$end | grep -v '>'|sed ':a;N;\$!ba;s/\\n//g'\n";
				    	$seq = `fastacmd -d TargetSequences/${contig3}.fa -s $contig3 -L $start,$end 2> error| grep -v '>'|sed ':a;N;\$!ba;s/\\n//g'`;
				    	if ($seq eq ''){
				    		my @buildtarget = "formatdb -o -i TargetSequences/${contig3}.fa -p F";
    						unless (-e "TargetSequences/${contig3}.nhr" && -e "TargetSequences/${contig3}.nin" && -e "TargetSequences/${contig3}.nsd" && -e "TargetSequences/${contig3}.nsi" && -e "TargetSequences/${contig3}.nsq"){
							system (@buildtarget) == 0 or die "Formatdb with target genome failed. $!\n";
							}
							$seq = `fastacmd -d TargetSequences/${contig3}.fa -s $contig3 -L 0,0| grep -v '>'|sed ':a;N;\$!ba;s/\\n//g'`;					
				    	}					
				    	chomp $seq;
					};			
			    }
			    elsif($majorstrand eq "+" && $strand_overall_negative != 0){
					if ($split{$para_original}{$contig3}{$exon3}[2] eq "-"){
				    	if ($mode eq "fasta"|| $mode eq "user"){
							$start = $split{$para_original}{$contig3}{$exon3}[1] - 30;
							$end = $split{$para_original}{$contig3}{$exon3}[0] + 30;
				    	}
				    	elsif ($mode eq "alignment"){
							#print "TEST:$split{$para_original}{$contig3}{$exon3}[1]\n";
							$start = $split{$para_original}{$contig3}{$exon3}[1] - 30;
							$end = $split{$para_original}{$contig3}{$exon3}[0] + 30;
				    	}			
				    	#print Dumper ("14");				    		
				    	my @buildtarget = "formatdb -o -i TargetSequences/${contig3}.fa -p F";
    					unless (-e "TargetSequences/${contig3}.nhr" && -e "TargetSequences/${contig3}.nin" && -e "TargetSequences/${contig3}.nsd" && -e "TargetSequences/${contig3}.nsi" && -e "TargetSequences/${contig3}.nsq"){
						system (@buildtarget) == 0 or die "Formatdb with target genome failed. $!\n";
						}
				    	$seq = `fastacmd -d TargetSequences/${contig3}.fa -s $contig3 -L $start,$end -S2 2> error| grep -v '>'|sed ':a;N;\$!ba;s/\\n//g'`;
				    	#print "fastacmd -d TargetSequences/${contig3}.fa -s $contig3 -L $start,$end -S2| grep -v '>'|sed ':a;N;\$!ba;s/\\n//g'\n";
				    	if ($seq eq ''){
				    		#print Dumper ("15");
				    		my @buildtarget = "formatdb -o -i TargetSequences/${contig3}.fa -p F";
    						unless (-e "TargetSequences/${contig3}.nhr" && -e "TargetSequences/${contig3}.nin" && -e "TargetSequences/${contig3}.nsd" && -e "TargetSequences/${contig3}.nsi" && -e "TargetSequences/${contig3}.nsq"){
							system (@buildtarget) == 0 or die "Formatdb with target genome failed. $!\n";
							}
							$seq = `fastacmd -d TargetSequences/${contig3}.fa -s $contig3 -L 0,0 -S2| grep -v '>'|sed ':a;N;\$!ba;s/\\n//g'`;					
				    	}					
				    	chomp $seq;
					}
					elsif ($split{$para_original}{$contig3}{$exon3}[2] eq "+"){
				    	if ($mode eq "fasta"|| $mode eq "user"){
							$start = $split{$para_original}{$contig3}{$exon3}[0] - 30;
							$end = $split{$para_original}{$contig3}{$exon3}[1] + 30;
				    	}
				    	elsif ($mode eq "alignment"){
							$start = $split{$para_original}{$contig3}{$exon3}[0] - 30;
							$end = $split{$para_original}{$contig3}{$exon3}[1] + 30;
				    	}
				    	#print Dumper ("16");
				    	my @buildtarget = "formatdb -o -i TargetSequences/${contig3}.fa -p F";
    					unless (-e "TargetSequences/${contig3}.nhr" && -e "TargetSequences/${contig3}.nin" && -e "TargetSequences/${contig3}.nsd" && -e "TargetSequences/${contig3}.nsi" && -e "TargetSequences/${contig3}.nsq"){
						system (@buildtarget) == 0 or die "Formatdb with target genome failed. $!\n";
						}
				    	#print "fastacmd -d TargetSequences/${contig3}.fa -s $contig3 -L $start,$end| grep -v '>'|sed ':a;N;\$!ba;s/\\n//g'";
				    	$seq = `fastacmd -d TargetSequences/${contig3}.fa -s $contig3 -L $start,$end 2> error| grep -v '>'|sed ':a;N;\$!ba;s/\\n//g'`;
				    	if ($seq eq ''){
				    		#print Dumper ("17");
				    		my @buildtarget = "formatdb -o -i TargetSequences/${contig3}.fa -p F";
    						unless (-e "TargetSequences/${contig3}.nhr" && -e "TargetSequences/${contig3}.nin" && -e "TargetSequences/${contig3}.nsd" && -e "TargetSequences/${contig3}.nsi" && -e "TargetSequences/${contig3}.nsq"){
							system (@buildtarget) == 0 or die "Formatdb with target genome failed. $!\n";
							}
							$seq = `fastacmd -d TargetSequences/${contig3}.fa -s $contig3 -L 0,0 | grep -v '>'|sed ':a;N;\$!ba;s/\\n//g'`;					
				    	}					
				    	chomp $seq;
					};			
			    }
			    elsif($majorstrand eq "-" && $strand_overall_positive != 0){
					if ($split{$para_original}{$contig3}{$exon3}[2] eq "-"){
				    	if ($mode eq "fasta"|| $mode eq "user"){
						$start = $split{$para_original}{$contig3}{$exon3}[1] - 30;
						$end = $split{$para_original}{$contig3}{$exon3}[0] + 30;
				    	}
				    	elsif ($mode eq "alignment"){
							$start = $split{$para_original}{$contig3}{$exon3}[1] - 30;
							$end = $split{$para_original}{$contig3}{$exon3}[0] + 30;
				    	}
				    	#print Dumper ("18");
				    	my @buildtarget = "formatdb -o -i TargetSequences/${contig3}.fa -p F";
    					unless (-e "TargetSequences/${contig3}.nhr" && -e "TargetSequences/${contig3}.nin" && -e "TargetSequences/${contig3}.nsd" && -e "TargetSequences/${contig3}.nsi" && -e "TargetSequences/${contig3}.nsq"){
						system (@buildtarget) == 0 or die "Formatdb with target genome failed. $!\n";
						}
				    	#print "fastacmd -d TargetSequences/${contig3}.fa -s $contig3 -L $start,$end| grep -v '>'|sed ':a;N;\$!ba;s/\\n//g'\n";					
				    	$seq = `fastacmd -d TargetSequences/${contig3}.fa -s $contig3 -L $start,$end 2> error| grep -v '>'|sed ':a;N;\$!ba;s/\\n//g'`;
				    	if ($seq eq ''){
				    		#print Dumper ("19");
				    		my @buildtarget = "formatdb -o -i TargetSequences/${contig3}.fa -p F";
    						unless (-e "TargetSequences/${contig3}.nhr" && -e "TargetSequences/${contig3}.nin" && -e "TargetSequences/${contig3}.nsd" && -e "TargetSequences/${contig3}.nsi" && -e "TargetSequences/${contig3}.nsq"){
							system (@buildtarget) == 0 or die "Formatdb with target genome failed. $!\n";
							}
							$seq = `fastacmd -d TargetSequences/${contig3}.fa -s $contig3 -L 0,0 | grep -v '>'|sed ':a;N;\$!ba;s/\\n//g'`;					
				    	}					
				    	chomp $seq;
					}
					elsif ($split{$para_original}{$contig3}{$exon3}[2] eq "+"){
				    	if ($mode eq "fasta"|| $mode eq "user"){
							$start = $split{$para_original}{$contig3}{$exon3}[0] - 30;
							$end = $split{$para_original}{$contig3}{$exon3}[1] + 30;
				    	}
				    	elsif ($mode eq "alignment"){
							$start = $split{$para_original}{$contig3}{$exon3}[0] - 30;
							$end = $split{$para_original}{$contig3}{$exon3}[1] + 30;
				    	}
				    	#print Dumper ("8");
				    		my @buildtarget = "formatdb -o -i TargetSequences/${contig3}.fa -p F";
    						unless (-e "TargetSequences/${contig3}.nhr" && -e "TargetSequences/${contig3}.nin" && -e "TargetSequences/${contig3}.nsd" && -e "TargetSequences/${contig3}.nsi" && -e "TargetSequences/${contig3}.nsq"){
							system (@buildtarget) == 0 or die "Formatdb with target genome failed. $!\n";
							}
				    	#print "fastacmd -d TargetSequences/${contig3}.fa -s $contig3 -L $start,$end -S2| grep -v '>'|sed ':a;N;\$!ba;s/\\n//g'\n";
				    	$seq = `fastacmd -d TargetSequences/${contig3}.fa -s $contig3 -L $start,$end -S2 2> error| grep -v '>'|sed ':a;N;\$!ba;s/\\n//g'`;
				    	if ($seq eq ''){
				    		#print Dumper ("7");
				    		my @buildtarget = "formatdb -o -i TargetSequences/${contig3}.fa -p F";
    						unless (-e "TargetSequences/${contig3}.nhr" && -e "TargetSequences/${contig3}.nin" && -e "TargetSequences/${contig3}.nsd" && -e "TargetSequences/${contig3}.nsi" && -e "TargetSequences/${contig3}.nsq"){
							system (@buildtarget) == 0 or die "Formatdb with target genome failed. $!\n";
							}
							$seq = `fastacmd -d TargetSequences/${contig3}.fa -s $contig3 -L 0,0 -S2 | grep -v '>'|sed ':a;N;\$!ba;s/\\n//g'`;					
				    	}					
				    	chomp $seq;		
					}			
			    };
			    $split{$para_original}{$contig3}{$exon3}[3] = $seq;
			};
	    };
	};
		#print Dumper (\%split);
		#For the other contigs, find the next exon upstream and downstream on maxcontig 
		#Extract the sequence between these two exons on maxcontig, search for more than three consecutive Ns and substitute the N-region with the putative new exon (stored in %split)
		#If the N-region is longer, fill up with Ns, if there are no Ns in the region, insert the sequence 10 nucleotides after the previous exon hit together with a N before and after
		#Run exonerate on this merged sequence
		#print "SORTED: @maxsorted_exons\n";
		#print "$majorstrand\n$maxcontig\n";
		if ($majorstrand eq "+"){
		    for my $contigMinorPlus(sort{$a <=> $b} keys %exons_minor){
			my $contig4 = $exons_minor{$contigMinorPlus};
			my $exon4 = $contigMinorPlus;
			#for my $contig4(sort {${$split{$para_original}}{$a} <=> ${$split{$para_original}}{$b}} keys %{$split{$para_original}}){		
			#print "$contig4\t$contigMinorPlus\n";		
			if ($contig4 ne $maxcontig){
			    my $find = "";				
			    my $addition;
			    my $replace;
			    #print "$contig4\t$maxcontig\n";
			    #for my $exon4(sort {$a <=> $b} keys %{$split{$para_original}{$contig4}}){
			    my @greater = ();
			    my @smaller= ();
			    my $after = ();
			    my $before = ();
			    my @greater_sorted =();
			    my @smaller_sorted=();
			    #print "$exon4\n";				
			    for my $i (0..$#maxsorted_exons){
				if ($exon4 > $maxsorted_exons[$i]){	
				    push @greater, $maxsorted_exons[$i];
				}elsif($exon4 < $maxsorted_exons[$i]){
				    push @smaller, $maxsorted_exons[$i];				
				}
			    };
			    #print "SMALLER/GREATER: @smaller\t@greater\n";
			    if ($fastasplit eq ""){						    		
			    	my @buildtarget = "formatdb -o -i TargetSequences/${maxcontig}.fa -p F";
	    			unless (-e "TargetSequences/${maxcontig}.nhr" && -e "TargetSequences/${maxcontig}.nin" && -e "TargetSequences/${maxcontig}.nsd" && -e "TargetSequences/${maxcontig}.nsi" && -e "TargetSequences/${maxcontig}.nsq"){
					system (@buildtarget) == 0 or die "Formatdb with target genome failed. $!\n";
					}		
					$fastasplit = `fastacmd -d TargetSequences/${maxcontig}.fa -s $maxcontig| grep -v '>'|sed ':a;N;\$!ba;s/\\n//g'`;
					chomp $fastasplit;
			    };
			    if (scalar(@smaller) == 0){
				my $fastasplit_old = $fastasplit;				
				$fastasplit = $fastasplit_old . $split{$para_original}{$contig4}{$exon4}[3];
				
			    }elsif (scalar(@greater) == 0){

				my $fastasplit_old = $fastasplit;				
				$fastasplit = $split{$para_original}{$contig4}{$exon4}[3] . $fastasplit_old;									
				$difference = (length ($fastasplit) - length ($fastasplit_old));			
			
			    }
			    elsif(scalar(@greater) != 0 && scalar(@smaller) != 0){
				@greater_sorted = sort {$a <=> $b} @greater;
				@smaller_sorted = sort {$a <=> $b} @smaller;
				$before = $greater_sorted[-1];
				$after = $smaller_sorted[0];	

				if ($split{$para_original}{$maxcontig}{$smaller_sorted[0]}[2] eq "+"){
				    #print "$maxcontig, $majorstrand\n";
				    #print "$split{$para_original}{$maxcontig}{$before}[1], $split{$para_original}{$maxcontig}{$after}[0]\t$difference\n";
				    $find = substr ($fastasplit, $split{$para_original}{$maxcontig}{$before}[1] - 1 + $difference, ($split{$para_original}{$maxcontig}{$after}[0] - $split{$para_original}{$maxcontig}{$before}[1] + 1 + $difference));			
				}elsif ($split{$para_original}{$maxcontig}{$smaller_sorted[0]}[2] eq "-"){
				    #print "$split{$para_original}{$maxcontig}{$after}[0], $split{$para_original}{$maxcontig}{$before}[1]\t$difference\n";						
				    $find = substr ($fastasplit, $split{$para_original}{$maxcontig}{$after}[0] - 1 + $difference, ($split{$para_original}{$maxcontig}{$before}[1] - $split{$para_original}{$maxcontig}{$after}[0] + 1 + $difference));			
				}

				my $seq = $split{$para_original}{$contig4}{$exon4}[3];	
				my $switch = 0;	
			
				if ($find =~ /((N){3,})/){
				    if (length $1 > length ($seq)){	
					$addition = "N" x (length ($1) - length ($seq) - 2);
				    }else {
					$addition = "N";	
				    };
				    $_ = $find;
				    s/((N){3,})/NN${seq}${addition}/;
				    $replace = $_;				
				}else{
				    $replace = substr ($find, 0, length($find)-10). "N" . $seq . "N" . substr ($find, length($find)-10, 10);
				};

				my $fastasplit_old = $fastasplit;
				my $pos = index($fastasplit,$find);
				substr ($fastasplit, $pos, length($find), $replace);
				$difference = length ($fastasplit) - length ($fastasplit_old);
			    };		
			}			
		    };	
		}
		elsif($majorstrand eq "-"){
		    for my $contigMinorMinus(sort{$b <=> $a} keys %exons_minor){
		    	#print Dumper (%exons_minor);
			my $contig4 = $exons_minor{$contigMinorMinus};
			my $exon4 = $contigMinorMinus;
			#print Dumper("contig 4 ${contig4}", "exon 4 ${exon4}");
		
			if ($contig4 ne $maxcontig){
			    my $find = "";				
			    my $addition;
			    my $replace;
			    my @greater = ();
			    my @smaller= ();
			    my $after = ();
			    my $before = ();
			    my @greater_sorted =();
			    my @smaller_sorted=();
				
			    for my $i (0..$#maxsorted_exons){
				if ($exon4 > $maxsorted_exons[$i]){	
				    push @greater, $maxsorted_exons[$i];
				}elsif($exon4 < $maxsorted_exons[$i]){
				    push @smaller, $maxsorted_exons[$i];				
				}
			    };

			    if ($fastasplit eq ""){				    	
			    	my @buildtarget = "formatdb -o -i TargetSequences/${maxcontig}.fa -p F";
	    			unless (-e "TargetSequences/${maxcontig}.nhr" && -e "TargetSequences/${maxcontig}.nin" && -e "TargetSequences/${maxcontig}.nsd" && -e "TargetSequences/${maxcontig}.nsi" && -e "TargetSequences/${maxcontig}.nsq"){
					system (@buildtarget) == 0 or die "Formatdb with target genome failed. $!\n";
					}		
					$fastasplit = `fastacmd -d TargetSequences/${maxcontig}.fa -s $maxcontig| grep -v '>'|sed ':a;N;\$!ba;s/\\n//g'`;
					chomp $fastasplit;
			    };
			    if (scalar(@smaller) == 0){
			    	#print Dumper ($fastasplit, $split{$para_original}{$contig4}{$exon4}[3]);
				my $fastasplit_old = $fastasplit;				
				$fastasplit = $split{$para_original}{$contig4}{$exon4}[3] . $fastasplit_old;
				#my $difference_old = $difference;				
				$difference = length ($fastasplit) - length ($fastasplit_old);	
			    }elsif (scalar(@greater) == 0){
				my $fastasplit_old = $fastasplit;				
				$fastasplit = $fastasplit_old . $split{$para_original}{$contig4}{$exon4}[3];						
			    }elsif(scalar(@greater) != 0 && scalar(@smaller) != 0){
				@greater_sorted = sort {$a <=> $b} @greater;
				@smaller_sorted = sort {$a <=> $b} @smaller;

				$before = $greater_sorted[-1];
				$after = $smaller_sorted[0];	

				if ($split{$para_original}{$maxcontig}{$smaller_sorted[0]}[2] eq "+"){

				    $find = substr ($fastasplit, $split{$para_original}{$maxcontig}{$before}[1] - 1 + $difference, ($split{$para_original}{$maxcontig}{$after}[0] - $split{$para_original}{$maxcontig}{$before}[1] + 1 + $difference));			
				}elsif ($split{$para_original}{$maxcontig}{$smaller_sorted[0]}[2] eq "-"){
				    $find = substr ($fastasplit, $split{$para_original}{$maxcontig}{$after}[0] - 1 + $difference, ($split{$para_original}{$maxcontig}{$before}[1] - $split{$para_original}{$maxcontig}{$after}[0] + 1 + $difference));			
				}

				my $seq = $split{$para_original}{$contig4}{$exon4}[3];	
				my $switch = 0;	
				if ($find =~ /((N){3,})/){
				    if (length $1 > length ($seq)){	
					$addition = "N" x (length ($1) - length ($seq) - 2);
				    }else {
					$addition = "N";	
				    };
				    $_ = $find;
				    if ($strand eq "+"){
					s/((N){3,})/$addition${seq}NN/;
					$replace = $_;				
				    }elsif($strand eq "-"){
					s/((N){3,})/NN${seq}${addition}/;
					$replace = $_;
				    }						
				}else{
				    $replace = substr ($find, 0, 10). "N" . $seq . "N" . substr ($find, 10, length($find));
				};
				
				my $fastasplit_old = $fastasplit;
				my $pos = index($fastasplit,$find);
				substr ($fastasplit, $pos, length($find), $replace);
				$difference = length ($fastasplit) - length ($fastasplit_old);
			    };	
			}	
		    }
		}

		$_ = $fastasplit;	
		s/.{32000}/$&\n/g;
		$fastasplit = $_;		
		print FINAL ">${para_original}_merged\n$fastasplit\n";
		close FINAL;

	    $_ = $para_original;			
	    s/[0-9]\.//;
	    my $querypara = $_;
	    my $querypara_old = (split "-", $para_original)[1];
	    if ($strand_overall_positive == 0){ #All exons of this paralog are on the negative strand			
	    	print Dumper ("perl $ExceS_A_2 --target TargetSequences/${para_original}_merged.fa --query ${querypara_old} --strand - --paralog $para_original --mode 4");
			system ("perl $ExceS_A_2 --target TargetSequences/${para_original}_merged.fa --query querypara_old} --strand - --paralog $para_original --mode 4");
	    }
	    elsif ($strand_overall_negative == 0){ #All exons of this paralog are on the positive strand		    	
	    	print Dumper ("perl $ExceS_A_2 --target TargetSequences/${para_original}_merged.fa --query ${querypara_old} --strand + --paralog $para_original --mode 4");
			system ("perl $ExceS_A_2 --target TargetSequences/${para_original}_merged.fa --query ${querypara_old} --strand + --paralog $para_original --mode 4");
		}
	    elsif ($majorstrand eq "-"){ #The exons on the major contig are on the negative strand, but at least one exon is on the positive strand		    	
	    	print Dumper ("perl $ExceS_A_2 --target TargetSequences/${para_original}_merged.fa --query ${querypara_old} --strand - --paralog $para_original --mode 4");
			system ("perl $ExceS_A_2 --target TargetSequences/${para_original}_merged.fa --query ${querypara_old} --strand - --paralog $para_original --mode 4");
	    }
	    elsif ($majorstrand eq "+"){ #The exons on the major contig are on the positive strand, but at least one exon is on the negative strand		    	
	    	print Dumper ("perl $ExceS_A_2 --target TargetSequences/${para_original}_merged.fa --query ${querypara_old} --strand + --paralog $para_original --mode 4");
			system ("perl $ExceS_A_2 --target TargetSequences/${para_original}_merged.fa --query ${querypara_old} --strand + --paralog $para_original --mode 4");
		}
    }

    if (-e "Summary.fa"){
	unlink "Summary.fa";
    };
    #if (-e "SearchTarget/*-*_target*.fa"){
    `cat SearchTarget/*-*_target*.fa >> Summary.fa`;
	#}
	#else{
		#print "No paralog found\n";
	#}
    print "ExonMatchSolver pipeline completed\n";
}

my @ExonerateFiles = glob ("SearchTarget/*_target.aln");
foreach my $file(@ExonerateFiles){
    `rm $file`;
};

my @ExonerateFiles2 = glob ("SearchTarget/*_target_final.aln");
foreach my $file(@ExonerateFiles2){
    `rm $file`;
};

#unless ($mode eq "alignment"){
#`rm SearchTarget/*.blastout`;
#}

#`rm -r SingleExons`;

if (glob("tmp*")){
    `rm tmp*`;
}

if (-e "error"){
    `rm error`; 
}
if (-e "error.log"){
    `rm error.log`; 
}
if (-e "Exonerate_help"){
    `rm Exonerate_help`; 
}


#------------------------------------------------------------------------------------------------------
#FUNCTIONS

sub usage {
    print STDERR "\nExonMatchSolver\nA pipeline assisting curation of gene annotations\n";
    print STDERR "usage: ExonMatchSolver.pl -i <fasta> -o <dir> -mode <fasta/alignment/user> -target <file> -paraMax <integer> -[OPTIONS]\n";
    print STDERR "\n";
    print STDERR "[INPUT]\n";
    print STDERR " -mode <alignment/fasta/user>	mode to be run: alignment (HMM are built), fasta or user\n"; 
    print STDERR " -i <file>    			input file (protein sequences, one per paralog of the protein family in the query species)\n";
    print STDERR " -o <dir>    			output directory\n";
    print STDERR " -target <file>     		absolute path to the target genome\n";
    #print STDERR " -WGD <yes/no>			WGD expected? Default: no\n";
    #print STDERR " -noWGD <integer,integer,...> 	specify paralogs for which no WGD is expected separated by commas only (no space), if WGD is set to yes i.e. the specified paralogs are unduplicated\n";
    print STDERR " -s <yes/no>                  	spliced alignmnet programm used: exonerate or ProSplign, Default: no (ProSplign)\n";		
    print STDERR " -query <file>     		absolute path to the query genome, required in fasta-mode\n";
    print STDERR " -user <file>                	input file with paralog-specific and TCE-individual protein sequences (.fa), required in user- and alignment-mode\n";
    print STDERR " -m <path>                   	input folder with paralog-specific and TCE-individual protein alignments (.stk, one alignment per paralog- and exon), required in alignment-mode\n";
    print STDERR " -l <integer>			length cutoff for the preparational step in fasta-mode\n";	
    print STDERR " -z <integer>			Z-score cutoff applied during preparational step in fasta-mode. Default: 3 (3 standard deviations).\n";
    print STDERR "[OPTIONS for BLAST]\n";
    print STDERR " -dFirst <integer>		nucleotide distance considered upstream of the first blast hit found\n";
    print STDERR " -dLast <integer>		nucleotide distance considered downstream of the last blast hit found\n";	
    print STDERR "[OPTIONS for the ILP]\n";
    print STDERR " -paraMax <integer>	        maximal number of paralogs expected to be located in the target genome; this corresponds to the highest number of paralogs EMS will run the ILP for. If you are in doubt about the number of paralogs encoded, please specify slightly more. We will identify the exact number for you.\n";
    print STDERR " -int <integer/fraction>	integer: round bitscore to a number dividable by this number; fraction: return all optimal solutions within this fraction\n";
    print STDERR " -max <integer>			maximal number of solution to be returned\n";
    print STDERR "[OPTIONS for SCIPIO]\n";
    print STDERR "-scipio_opt <string>           	see scipio documentation for options\n";
    print STDERR " -c <integer>                  	number of cores available to use\n";
    print STDERR " -h <file>    		 	this (useful) help message\n";
    print STDERR "[VERSION]\n";
    print STDERR " 24-10-2017\n";
    print STDERR "[BUGS]\n";
    print STDERR " Please report bugs to henrike\@bioinf.uni-leipzig.de\n";
    print STDERR "\n";
    exit(-1);
}

sub prettyTime{
    my @months = qw(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec);
    my @weekDays = qw(Sun Mon Tue Wed Thu Fri Sat Sun);
    my ($second, $minute, $hour, $dayOfMonth, $month,
	$yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings) = localtime();
    my $year = 1900 + $yearOffset;
    return "$hour:$minute:$second, $weekDays[$dayOfWeek] $months[$month] $dayOfMonth, $year";
}

sub printHeader{
    print "\n";
    print "  ============================================================\n";
    print "  =                       ExonMatchSolver                    =\n";
    print "  =     A pipeline assisting curation of gene annotations    =\n";
    print "  ============================================================\n";
    print "   startet at ".prettyTime()."\n";
    print "   -> input file:       $input\n";
    print "   -> output folder:    $outfolder\n";
    print "\n";
}

sub uniq {
    my %seen;
    return grep { !$seen{$_}++ } @_;
}

sub get_exons {
    my $length_total = 0;
    my $length = 0;
    my @array;
    push @array, 0;
    if ($exonerate eq "no"){
	while(<INF>){
	    chomp;
	    if ($_ =~ m/id:/g){
		my @coordinates = split "\t", $_; 
		if ($coordinates[2] - $coordinates[3] > 0){
		    my $length = $coordinates[2] - $coordinates[3] + 1;
		    $length_total = $length + $length_total;
		    push @array, $length_total;
		}elsif ($coordinates[3] - $coordinates[2] > 0){
		    my $length = $coordinates[3] - $coordinates[2] + 1;
		    $length_total = $length + $length_total;
		    push @array, $length_total;
		}
	    }
	}
    }elsif ($exonerate eq "yes"){
	while(<GFF>){
	    chomp;
	    if ($_ =~ m/\texon\t/g){
		my @coordinates = split "\t", $_; 
		if ($coordinates[3] - $coordinates[4] > 0){
		    my $length = $coordinates[3] - $coordinates[4] + 1;
		    $length_total = $length + $length_total;
		    push @array, $length_total;
		}elsif ($coordinates[4] - $coordinates[3] > 0){
		    my $length = $coordinates[4] - $coordinates[3] + 1;
		    $length_total = $length + $length_total;
		    push @array, $length_total;
		}
	    }
	}
    }
    return (@array);
}

sub average{
    my( $mean );    #   the mean of the values, calculated by this subroutine
    my( $n );      #   the number of values
    my( $sum );    #   the sum of the values
    my( $x );      #   set to each value

    $n = 0;
    $sum = 0;
    foreach $x ( @_ ){
        $sum += $x;
        $n++;
    }
    $mean = $sum / $n;
    return $mean;      
}

sub stdev{
    my( $n );               #   the number of values
    my( $stddev );          #   the standard deviation, calculated by this
    #   subroutine, from the variance
    my( $sum );             #   the sum of the values
    my( $sumOfSquares );    #   the sum of the squares of the values
    my( $variance );        #   the variance, calculated by this subroutine
    my( $x );               #   set to each value

    $n = 0;
    $sum = 0;
    $sumOfSquares = 0;

    foreach $x ( @_ ){
	$sum += $x;
	$n++;
	$sumOfSquares += $x * $x;
    }
    $variance = ( $sumOfSquares - ( ( $sum * $sum ) / $n ) ) / ( $n - 1 );
    $stddev = sqrt( $variance );
    return $stddev;
}

sub getParaNum {
    my ($input) =  @_;
    my $X1 = 0;
    my $Y1 = 0;
    my $X2;
    my $Y2;
    my @M;
    my @diff;
    $M[0] =0;
    my $count = 0;
    my $differ;
    my $parastop;
    my $max = 0;
    my @monotony;

    open (INPUT, "<", $input);

    #From x=#paralog and y = unweighted score, calculate the slope (written to @M) and the difference of two adjacent slopes (@dif)
    while (<INPUT>) {
	chomp;
	$_ =~ s/^\s//g;
	if ($_ =~ /#paralogs/) {
	    my @paralogs = split ':', $_;
	    $X2 = $paralogs[1];
	    #}elsif ($_ =~ /#Solution/) {
	    #my @solution = split '[: ]', $_;
	    #$solution = $solution[1];
	}elsif($_ =~ /#unweighted score/) {
	    my @unweighted= split ': ', $_;
	    $Y2 = $unweighted[1];
	}elsif($_ =~ /#Best score/){
	    $count++;
	    #print "$paralog\t$solution\t$unweighted\n";    
	    #$X2 = $paralog;
	    #$Y2 = $unweighted;
	    #print Dumper($Y2, $X2);
	    $M[$count] = ($Y1 - $Y2)/($X1 - $X2);
	    $differ = $M[$count - 1] - $M[$count];
	    push @diff, $differ; 
	};
    };
#print Dumper(@diff);
#print Dumper($#diff);
#print Dumper(scalar @diff);

	if($#diff == 0){
		$parastop = 1;
		return $parastop;

	}
	else{
	    #Loop through the slope differences to extract monotony, also calculate maximal value of slope differences
	    for my $slopediff (1..$#diff) {
	    	#print Dumper ($diff[$slopediff-1], $diff[$slopediff]);
		if ($diff[$slopediff-1] > $diff[$slopediff]) {
		    $monotony[$slopediff-1] = 0;  
		}elsif($diff[$slopediff-1] < $diff[$slopediff]){
		    $monotony[$slopediff-1] = 1;
		};
		if ($diff[$max] < $diff[$slopediff]) {
		    $max = $slopediff;
		    #print Dumper ($max);
		}
	    };
	}
    #print "Monotonie: @monotony\n";
    #print "Max:$max\t$diff[$max]\n";
    #print "@M\n";

    #Loop through @monotony array to identify the index, which is the first rising one coming from the last element (highest number of paralogs)
    for my $para (reverse 0..$#monotony){
	if ($monotony[$para] eq "1") {
	    $parastop = $para +1;
	    return $parastop;
	    #print "Para: $parastop\n";
	    #print "Max: $max\n";
	    exit;
	}   
    }
}

sub CompartGenome{
    #Check if there are several exons on the same contig wich overlap in their query coverage 
    #my $check = `sort -k1,1n $blast | awk '{print \$1"\t"\$2"\t"\$7"\t"\$8}'| sort -n | uniq -d| wc -l`;
    #chomp $check;
    my $blastinside =  shift;
    #print "$blastinside\n";
    my $check =0;
    open(BLAST, "< $blastinside") or die ("Could not read blast\n");
    my %blast;
    #my $i = 0;
    while (<BLAST>) {
	chomp $_;
	#$i ++;
	my @blast = split /[\s+]+/, $_; 
	my $key = $blast[0] . "_" . $blast[1]; 
	push @{$blast{$key}}, [$blast[6], $blast[7]];    
    };
    close BLAST;
    foreach my $key(keys %blast){
	if (scalar @{$blast{$key}} > 1) {         
            foreach my $locus (0 .. $#{$blast{$key}}){
		my $prevlocus = 0;
		if ($blast{$key}[$locus][1] - $blast{$key}[$prevlocus][0] > 5) {
		    $check ++;

		}
		$prevlocus = $locus;
	    };
	    
	};
    };


    #Enter loop of cutting the genome if there are exons that cover the same part of query on same contig
    if ($check > 0) {
	my $targetMod = "${genome_target}_${IDs[0]}Mod";
	if (-e $targetMod){
		`rm $targetMod`;
	}
	my $blastoutSynt = "compartments";
	#print "There might be two gene copies situated on the same genomic fragment on the target genome. Compartments will be identified.\n";
	my %compart;
	
	#Run blast of full length query seq against target genome in order to identify the compartments, call procompart to get compartement coordinates
	#print "`blastall -p tblastn -d $genome_target -i SearchQuery/${IDs[0]}.fa -e 0.01 -a 2 -m8 -G 11 -E 1 -C F| awk '{for(i=n;i<=NF;i++)\$(i-(n-1))=\$i;NF=NF-(n-1);print \"Paralog\\t\"\$0}' n=2 | ${path}procompart -max_extent 0 -t | sort -k2,2 -k4,4n > $blastoutSynt`\n";
	my @control = `blastall -p tblastn -d $genome_target -i SearchQuery/${IDs[0]}.fa -e 0.01 -a 2 -m8 -G 11 -E 1 -C F| awk '{for(i=n;i<=NF;i++)\$(i-(n-1))=\$i;NF=NF-(n-1);print "Paralog\\t"\$0}' n=2 | ${path}procompart -max_extent 0 -max_intron 100 -t | sort -k2,2 -k4,4n > $blastoutSynt`;
	system (@control) != 0 or die "blastall with query sequences against target genome failed. $!\n";    
	
	#Parse protcompart output and store in %compart
	my @loci;
	my %seen;
	open(COMPART, "< $blastoutSynt") or warn "Could not read $blastoutSynt\n";
        while (<COMPART>) {
	    chomp $_;
	    my @locus = split '\t', $_;
	    push (@{$seen{$locus[1]}}, $locus[0]);
	    $compart{$locus[0]} = [$locus[1], $locus[3], $locus[4]];
        };

	#Loop through every key of %seen, for every value in array -> record cutting points in %cutpoints
	#find cutpoint by looking at coordinates from previous and next array element
	my %cutpoints;
	foreach my $contigs (keys %seen){
	    if(@{$seen{$contigs}} > 1){	
		my $i =0;
		foreach my $locus (@{$seen{$contigs}}){
		    if ($i ne "0") {
			#print "$seen{$contigs}[$i-1]\n";
			$cutpoints{$seen{$contigs}[$i]} = [$contigs, $compart{$seen{$contigs}[$i-1]}[2], ];
		    }else{$cutpoints{$seen{$contigs}[$i]} = [$contigs, 1, ]};
		    if ($i ne  $#{$seen{$contigs}}) {
			$cutpoints{$seen{$contigs}[$i]}[2] = $compart{$seen{$contigs}[$i+1]}[1];
		    }else{$cutpoints{$seen{$contigs}[$i]}[2] = "last"};
		    $i ++;
		}
	    };
	}	
	
	#Loop through genome file and write the modified genome in new file. Modified genome is cut in compartments for genes that lie on the same contig.
	my $in  = Bio::SeqIO->new(-file => ${genome_target}, -format => 'Fasta');
	open(GENOME, ">$targetMod") or warn ("Could not open $targetMod for writing. $!\n");
	mkdir "TargetSequences/";
	while (my $seq = $in->next_seq) {
	    foreach my $contigs (keys %seen){
		if(@{$seen{$contigs}} > 1){
		    my $i =0;
		    if ($seq->display_id eq "$contigs") {
			my $last = $seq->length();
			foreach my $locus (@{$seen{$contigs}}){
			    if($cutpoints{$seen{$contigs}[$i]}[2] eq "last"){
				$cutpoints{$seen{$contigs}[$i]}[2] = $last;
			    };
			    print GENOME ">", $seq->display_id, "-", $locus, "\n", $seq->subseq($cutpoints{$seen{$contigs}[$i]}[1], $cutpoints{$seen{$contigs}[$i]}[2]), "\n";
			    open (LOCUS, "> TargetSequences/${contigs}-${locus}.fa");
			    print LOCUS ">", $seq->display_id, "-", $locus, "\n", $seq->subseq($cutpoints{$seen{$contigs}[$i]}[1], $cutpoints{$seen{$contigs}[$i]}[2]), "\n";
			    close LOCUS;
			    $i ++;
			};
		    }else{
			print GENOME ">", $seq->display_id, "\n", $seq->seq(), "\n"
		    };
		};
	    };
	}

	my $blastinsideMod = "${blastinside}Mod";

	#Modify orignial blastinput file
	open(BLASTMOD, "> $blastinsideMod");
	open(BLAST, "< $blastinside") or die ("Could not read $blastinside\n");
	while (<BLAST>) {
	    chomp $_;

	    my @blast = split /[\s]+/, $_;
	    my $finlocus;
	    my $switch =0;
	    foreach my $locus (reverse sort keys %cutpoints){
		#print "$blast[1] eq $cutpoints{$locus}[0] && $blast[8] <= $cutpoints{$locus}[2] && $blast[9] <= $cutpoints{$locus}[2] && $blast[8] >= $cutpoints{$locus}[1] && $blast[9] >= $cutpoints{$locus}[1]\n";
		if ($blast[1] eq $cutpoints{$locus}[0] && $blast[8] <= $cutpoints{$locus}[2] && $blast[9] <= $cutpoints{$locus}[2] && $blast[8] >= $cutpoints{$locus}[1] && $blast[9] >= $cutpoints{$locus}[1]) {
		    $blast[28] = $blast[8] - $cutpoints{$locus}[1] + 1;
		    $blast[29] = $blast[9] - $cutpoints{$locus}[1] + 1;
		    $finlocus = $locus;
		    $switch =1;
		};
	    };
	    if ($switch eq "0") {
		print BLASTMOD "$blast[0]\t$blast[1]\t$blast[2]\t$blast[3]\t$blast[4]\t$blast[5]\t$blast[6]\t$blast[7]\t$blast[8]\t$blast[9]\t$blast[10]\t$blast[11]\n";		
	    }elsif($switch eq "1"){
		print BLASTMOD "$blast[0]\t$blast[1]-${finlocus}\t$blast[2]\t$blast[3]\t$blast[4]\t$blast[5]\t$blast[6]\t$blast[7]\t$blast[28]\t$blast[29]\t$blast[10]\t$blast[11]\n";
	    }
	};
	close BLAST;
	close BLASTMOD;
	my $blastinsideModFinal = "${blastinside}Modfinal";
	`sort -n $blastinsideMod| uniq > $blastinsideModFinal`;
	return ($blastinsideModFinal, \%cutpoints);
    }else{
	my %cutpoints;
	return ($blastinside, \%cutpoints);
    }
};

sub CompartGenomeHMM{		
    my $targetAAS = "${space}${myOutFolder}_ORF.aas";				
    #Check if there are several exons on the same contig wich overlap in their query coverage in the hmmoutput file for single exons
    my $hmminside =  shift;
    my $check =0;
    open(BLAST, "< $hmminside") or die ("Could not read hmm output\n");
    my %blast;
    #my $i = 0;
    while (<BLAST>) {
	chomp $_;
	unless ($_ =~ /#/){
	    #print "$_\n";
	    my @blast = split '\s+', $_;
	    #print "$blast[3]\n";
	    my $key = $blast[3] . "_" . $blast[0]; 
	    #print Dumper ($key);
	    push @{$blast{$key}}, [$blast[15], $blast[16]];    
	}
    };
    close BLAST;
    #print Dumper (\%blast);
    foreach my $key(keys %blast){
	if (scalar @{$blast{$key}} > 1) {         
            foreach my $locus (0 .. $#{$blast{$key}}){
		my $prevlocus = 0;
		if ($blast{$key}[$locus][1] - $blast{$key}[$prevlocus][0] > 5) {
		    $check ++;
		    #print "Overlap: $blast{$key}[$locus][0], $blast{$key}[$locus][1], $blast{$key}[$prevlocus][0], $blast{$key}[$prevlocus][1] \n";
		}else{
		    #print"Case not covered\n";
		    #print "$key\t$locus: $blast{$key}[$locus][0], $blast{$key}[$locus][1], $blast{$key}[$prevlocus][0], $blast{$key}[$prevlocus][1] \n";      
		};
		$prevlocus = $locus;
	    };  
	}
	else{
		#print Dumper ("yes")
		return ($hmminside, \%blast);
	};
    };

    #Enter loop of cutting the genome if there are exons that cover the same part of query on same contig
    if ($check > 0) {
	my $targetMod = "${genome_target}Mod";
	my $blastoutSynt = "compartments";
	#print "There might be a tandem duplication. Compartments will be identified.\n";
	my %compart;
	my %hmm;
	
	#Run hmmsearch of full length query seq against target genome in order to identify the compartments, call procompart to get compartement coordinates
	`clustalo --force --outfmt=st -i $myInput -o ${myInput}.stk`;
	#system (@clustal) != 0 or die "Alignment of query sequences with clustal failed. $!\n";    
	`hmmbuild --amino ${myInput}.model ${myInput}.stk`;
	#system (@hmmbuild) != 0 or die "Hmmbuild of query sequence alignment (${myInput}.model) failed. $!\n";    
	my @control = `hmmsearch --domtblout ${myInput}_domtblout -E 0.01 -o ${myInput}_out ${myInput}.model $targetAAS`;
	system (@control) != 0 or die "hmmsearch with query sequences against target genome failed. $!\n";    
	#print "Hmmsearch with full-length query models done against target genome\n";
	
	#Reformat hmmsearch output to get a blast like output and run procompart with this output 
	open(HMM, "< ${myInput}_domtblout") or warn "${myInput}_domtblout could not be opened.\n";
	open(REFHMM, "> ${myInput}_domtblout_blast");
	while (<HMM>) {
	    chomp $_;
	    unless ($_ =~ /#/){
	    		    print "$_\n";
		my @hmm = split '\s+', $_;
		my @coordinates = split '[\[\]-]', $hmm[22];
		my $length = $hmm[16] - $hmm[15];
		if ($coordinates[0] =~ /FORWARD/) {        
		    my $targetstart = $coordinates[1] + 3*$hmm[19];
		    my $targetstopp =  $coordinates[1] + 3*$hmm[20];
		    print REFHMM "$hmm[3]\t$hmm[0]\t1\t$length\t1\t1\t$hmm[15]\t$hmm[16]\t$targetstart\t$targetstopp\t$hmm[6]\t$hmm[7]\n";
		}elsif($coordinates[0] =~ /REVERSE/){
		    my $targetstart = $coordinates[2] - 3*$hmm[19];
		    my $targetstopp =  $coordinates[2] - 3*$hmm[20];
		    print REFHMM "$hmm[3]\t$hmm[0]\t1\t$length\t1\t1\t$hmm[15]\t$hmm[16]\t$targetstart\t$targetstopp\t$hmm[6]\t$hmm[7]\n";
		};
	    }
	};
	`cat ${myInput}_domtblout_blast| ${path}procompart -max_extent 0 -max_intron 100 -t| sort -k2,2 -k4,4n > $blastoutSynt`;

	#Parse procompart output and store in %compart
	my @loci;
	my %seen;
	open(COMPART, "< $blastoutSynt") or warn "Could not read $blastoutSynt\n";
        while (<COMPART>) {
	    chomp $_;
	    my @locus = split '\t', $_;
	    #print Dumper (@locus);	    
	    push (@{$seen{$locus[1]}}, $locus[0]);
	    $compart{$locus[0]} = [$locus[1], $locus[3], $locus[4]];
	    #print Dumper($compart{$locus[0]});
        };

	`rm ${myInput}_domtblout_blast`;
	`rm ${myInput}_domtblout`;
	`rm ${myInput}.model`;
	`rm ${myInput}_out`;
	`rm ${myInput}.stk`;
	
	#Loop through every key of %seen, for every value in array -> record cutting points in %cutpoints
	#find cutpoint by looking at coordinates from previous and next array element
	my %cutpoints;
	#print Dumper (\%seen);
	foreach my $contigs (keys %seen){
	    if(@{$seen{$contigs}} > 1){	
		my $i =0;
		foreach my $locus (@{$seen{$contigs}}){
		    #print "$seen{$contigs}[$i]\n";

		    if ($i ne "0") {
				#print "$seen{$contigs}[$i-1]\n";
				$cutpoints{$seen{$contigs}[$i]} = [$contigs, $compart{$seen{$contigs}[$i-1]}[2], ];
		    }
		    else{
		    	$cutpoints{$seen{$contigs}[$i]} = [$contigs, 1, ];
		    }
		    
		    #print "Last element: $#{$seen{$contigs}}\n";


		    if ($i ne  $#{$seen{$contigs}}) {
				$cutpoints{$seen{$contigs}[$i]}[2] = $compart{$seen{$contigs}[$i+1]}[1];
		    }
		    else{
		    	$cutpoints{$seen{$contigs}[$i]}[2] = "last"};
		    	$i ++;
		}
	    };
	}	
	#print Dumper (\%cutpoints);
	#Loop through genome file and write the modified scaffolds in a new file.
	#Genome file, nucleotide sequence
	my $in  = Bio::SeqIO->new(-file => $genome_target, -format => 'Fasta');
	open(GENOME, ">$targetMod") or warn ("Could not open $targetMod for writing. $!\n");
	mkdir "TargetSequences/";
	while (my $seq = $in->next_seq) {
	    foreach my $contigs (keys %seen){
		if(@{$seen{$contigs}} > 1){	
		    my $i =0;
		    if ($seq->display_id eq "$contigs") {
			my $last = $seq->length();
			foreach my $locus (@{$seen{$contigs}}){
			    if($cutpoints{$seen{$contigs}[$i]}[2] eq "last"){
				$cutpoints{$seen{$contigs}[$i]}[2] = $last;
			    };
			    #print Dumper($cutpoints{$seen{$contigs}[$i]}[1], $cutpoints{$seen{$contigs}[$i]}[2] );
			    #if($cutpoints{$seen{$contigs}[$i]}[1] < 1){
			    #	$cutpoints{$seen{$contigs}[$i]}[1] = 1;
			    #}
			    print GENOME ">", $seq->display_id, "-", $locus, "\n", $seq->subseq($cutpoints{$seen{$contigs}[$i]}[1], $cutpoints{$seen{$contigs}[$i]}[2]), "\n";
			    open (LOCUS, "> TargetSequences/${contigs}-${locus}.fa");
			    print LOCUS ">", $seq->display_id, "-", $locus, "\n", $seq->subseq($cutpoints{$seen{$contigs}[$i]}[1], $cutpoints{$seen{$contigs}[$i]}[2]), "\n";
			    close LOCUS;
			    $i ++;
			    #print Dumper ("getorf -minsize 9 -reverse Y -stdout Y -sequence TargetSequences/${contigs}-${locus}.fa -outseq stdout|sed 's/ - /-/g'| awk '{print \$1, \$3, \$2}' > tmp2");
			    my $ret = system "getorf -minsize 9 -reverse Y -stdout Y -sequence TargetSequences/${contigs}-${locus}.fa -outseq stdout|sed 's/ - /-/g'| awk '{print \$1, \$3, \$2}' > tmp2";
			    if ($ret != 0){
				die;
			    };
			    open (TMP, "< tmp2") or die "Couldn't open tmp2, $!";
			    open (DEST, "> TargetSequences/${contigs}-${locus}.aas");
			    while (<TMP>){
				chomp;
				s/\(//;
				s/_(\d+)\s\s\[/ $1\(FORWARD\)\[/;
				s/\sREVERSE/\(REVERSE\)/;
				s/_(\d+)\(REVERSE\) / $1\(REVERSE\)/;
				print DEST "$_\n";
			    }
			    close TMP;
			    close DEST;
			    
			}
		    };
		};
	    };
	}	
	close GENOME;
	
	#Modify the hmmsearch original outputfile that was specific for exons/paralogs
	#print Dumper (\%cutpoints);
	my $blastinsideMod = "${hmminside}Mod";
	#print "$blastinsideMod\n";
	#Modify orignial blastinput file
	open(BLASTMOD, "> $blastinsideMod");
	#print "$hmminside\n";
	open(BLAST, "< $hmminside") or die ("Could not read $hmminside\n");
	while (<BLAST>) {
	    chomp $_;
	    #print "$_\n";
	    #print "Blast\n";
	    unless ($_ =~ /#/){
		my @hmm = split '\s+', $_;
		#print "$hmm[0]\n";
		my $targetstart;
		my $targetstopp;
		my @coordinates = split '[\[\]_]+', $hmm[22];
		my $length = $hmm[16] - $hmm[15];
		if ($coordinates[0] =~ /FORWARD/) {
		    $targetstart = $coordinates[1] + 3*$hmm[19];
		    $targetstopp =  $coordinates[1] + 3*$hmm[20];
		    print REFHMM "$hmm[3]\t$hmm[0]\t1\t$length\t1\t1\t$hmm[15]\t$hmm[16]\t$targetstart\t$targetstopp\t$hmm[6]\t$hmm[7]\n";
		}elsif($coordinates[0] =~ /REVERSE/){
		    #print "$coordinates[2]\n";
		    $targetstart = $coordinates[2] - 3*$hmm[19];
		    $targetstopp =  $coordinates[2] - 3*$hmm[20];
		    #print 
		    print REFHMM "$hmm[3]\t$hmm[0]\t1\t$length\t1\t1\t$hmm[15]\t$hmm[16]\t$targetstart\t$targetstopp\t$hmm[6]\t$hmm[7]\n";
		};
		my $finlocus;
		#if ($cutpoints{$hmm[0]}) {
		my $switch = 0;
		foreach my $locus (reverse sort keys %cutpoints){
		    if ($hmm[0] eq $cutpoints{$locus}[0] && $targetstart <= $cutpoints{$locus}[2] && $targetstopp <= $cutpoints{$locus}[2] && $targetstart >= $cutpoints{$locus}[1] && $targetstopp >= $cutpoints{$locus}[1]) {
			$finlocus = $locus;
			$switch = 1;
		    };       
		};
		if ($switch eq "0") {
		    print BLASTMOD "$hmm[0]\t$hmm[1]\t$hmm[2]\t$hmm[3]\t$hmm[4]\t$hmm[5]\t$hmm[6]\t$hmm[7]\t$hmm[8]\t$hmm[9]\t$hmm[10]\t$hmm[11]\t$hmm[12]\t$hmm[13]\t$hmm[14]\t$hmm[15]\t$hmm[16]\t$hmm[17]\t$hmm[18]\t$hmm[19]\t$hmm[20]\t$hmm[21]\t$hmm[22]\n";          
		}elsif($switch eq "1"){
		    my $newstart = $targetstart - $cutpoints{$finlocus}[1] + 1;
		    my $newstopp = $targetstopp - $cutpoints{$finlocus}[1] + 1;
		    my $strand;
		    if ($newstart < $newstopp) {
			$strand = "(FORWARD)";	
		    }else {$strand = "(REVERSE)";
		    };
		    print BLASTMOD "$hmm[0]-${finlocus}\t$hmm[1]\t$hmm[2]\t$hmm[3]\t$hmm[4]\t$hmm[5]\t$hmm[6]\t$hmm[7]\t$hmm[8]\t$hmm[9]\t$hmm[10]\t$hmm[11]\t$hmm[12]\t$hmm[13]\t$hmm[14]\t$hmm[15]\t$hmm[16]\t$hmm[17]\t$hmm[18]\t$hmm[19]\t$hmm[20]\t$hmm[21]\t${strand}[${newstart}_${newstopp}]\n";                        			
		}
	    };
	}
	close BLAST;
	close BLASTMOD;
	my $blastinsideModFinal = "${hmminside}Modfinal";
	`sort -n $blastinsideMod| uniq > $blastinsideModFinal`;
	return ($blastinsideModFinal, \%cutpoints);
    }
};

sub runSplicedAligner{
    my $input = $_[0];	
    #print Dumper ($input);
    my %locus =  %{$_[1]};
    #print Dumper (%locus);
    my $strand;

    #Extract missing contigs  
    foreach my $contig(keys %locus){
		unless (-s "TargetSequences/${contig}.fa"){
		    `fastacmd -d $genome_target -s $contig| sed \'s/>lcl|/>/g\' > TargetSequences/${contig}.fa`;		    
		    `formatdb -i TargetSequences/${contig}.fa -p F -o`;
		}
    };
    foreach my $contig(keys %locus){
    	#print Dumper($contig);
		$_ = ${locus{$contig}[0]};
		my ($type, $querypara) = split "-", $_;
		my @output = ExtractFirstLast($input, $contig, ${locus{$contig}[0]});
		if (! @output || ! $output[1]){
			next;
		}

		$strand = $output[0];

		my @first;
		my @last;
		$first[0] = $output[1];
		$first[1] = $output[2];
		$last[0] = $output[3];
		$last[1] = $output[4];

		if ($strand eq "+"){
		    ${locus{$contig}[2]} = $first[0] - $distance_first;
		    ${locus{$contig}[3]} = $last[1]+ $distance_last;		
		}elsif($strand eq "-"){
		    ${locus{$contig}[3]} = $first[1] + $distance_first;
		    ${locus{$contig}[2]} = $last[0] - $distance_last;
		};
		if (${locus{$contig}[2]} < 0){
			${locus{$contig}[2]} = 0;
		}

		my $genome_seqIO = Bio::SeqIO->new(-file   => "TargetSequences/${contig}.fa", -format => "Fasta");
		$a = $genome_seqIO->next_seq;
		$a = $a -> seq;

		if (${locus{$contig}[3]} > (length $a)){
			${locus{$contig}[3]}  = (length $a);
		}

		print Dumper ("perl $ExceS_A_2 --contig $contig --query ${locus{$contig}[0]} --strand $strand --start ${locus{$contig}[2]} --end ${locus{$contig}[3]} --mode 2");
		system("perl $ExceS_A_2 --contig $contig --query ${locus{$contig}[0]} --strand $strand --start ${locus{$contig}[2]} --end ${locus{$contig}[3]} --mode 2");
    };
}


sub ExtractFirstLast{
    my $input =  shift;
    my $contig =  shift;    

    #print "$input\n";
    #print "$contig\n";
    open (IN, "< $input") or die "Could not open $input.\n"; 
    my $strand;
    my %strand;
    my %hit;
    my %valid;
    $strand{"+"} = 0;
    $strand{"-"} = 0;
    my $first;
    my $first2;
    my $last2;
    my $last;
    my $counter =0;
    my $locusprev;
    my $exonprev;
    my $check = "0";

    while(<IN>){
	chomp $_;
	#print "$_\n";
	if ($_ =~ /^$contig/) {
	    #print "$_\n";
	    my @entry = split '[\t]', $_;
	    my @hit = split '_', $entry[1];
	    $hit[0] =~ s/^.*\[//g;
	    $hit[1] =~ s/]//g;
	    #print Dumper ($hit[0]);
	    #print Dumper ($hit[1]);
	    if ($hit[0] > $hit[1] && $entry[3] < 0.001) {
		$strand{"-"}++;    
	    }elsif($hit[0] < $hit[1] && $entry[3] < 0.001){
		$strand{"+"}++;
	    };
	    $check = "1";
	}
    }
    close IN;

    if ($check eq "0") {
	die "No $contig in $input\n";
    }

#print Dumper ($strand{"+"});
#print Dumper ($strand{"-"});
    if ($strand{"+"} > $strand{"-"}) {
	$strand = "+";    
    }elsif($strand{"+"} < $strand{"-"}) {
	$strand = "-";    
    }
    else{
    	next;
    }

    #print "$strand\n";

    open (IN, "< $input") or die "Could not open $input.\n"; 
    while(<IN>){
	chomp $_;
	#print Dumper($_);
	#print Dumper ($contig);
	if ($_ =~ /^$contig/) {                
	    my @entry = split '[\t]', $_;
	    my @hit = split '_', $entry[1];
	    $hit[0] =~ s/^.*\[//g;
	    #print Dumper ($hit[0]);
	    $hit[1] =~ s/]//g;
	   	#print Dumper ($hit[1]);
	   	#print Dumper ($entry[3]);
	    my $strandlocal;
	    if ($hit[0] > $hit[1] && $entry[3] < 0.001) {
		$strandlocal = "-";
	    }elsif($hit[1] > $hit[0] && $entry[3] < 0.001){
		$strandlocal = "+";    
	    }
	    else{
	    	#print Dumper ("MÄH");
	    	next;
	    }
	    if ($strandlocal eq $strand) {
		$hit{$entry[7]} = [$hit[0], $hit[1]];
	    }
	};
    }

    foreach my $exons (sort {$a<=>$b} keys %hit){
    	#print Dumper ($exons);
		my $locus = $hit{$exons}[0];
		my $exon = $exons;
		#print "$locus\n";
		if ($locusprev) {
		    #print Dumper ("$strand eq + && $locusprev < $locus");
		    if ($strand eq "+" && $locusprev <= $locus) {
			$valid{$exonprev} = [$hit{$exonprev}[0], $hit{$exonprev}[1]];
			$valid{$exons} = [$hit{$exons}[0], $hit{$exons}[1]];
    		}
	    	elsif($strand eq "-" && $locusprev >= $locus){
			$valid{$exonprev} = [$hit{$exonprev}[0], $hit{$exonprev}[1]];
			$valid{$exons} = [$hit{$exons}[0], $hit{$exons}[1]];
	    	}
	    	else{
	    	#print Dumper ("MuH");
	    	next;
	    	};
		}
		$locusprev = $locus;
		$exonprev = $exon;
		#print "$locusprev\n";
    }

    my @output;
    push @output, $strand;
#print Dumper(%valid);
    foreach my $exons (sort {$a<=>$b} keys %valid){
	if($counter eq "0"){
	    push @output, $valid{$exons}[0];   
	    push @output, $valid{$exons}[1];
	};
	$counter ++;
	$last = $valid{$exons}[0];
	$last2 = $valid{$exons}[1];	
	#print Dumper($last);
	#print Dumper($last2);
    }

    push @output, $last;   
    push @output, $last2;
	#print Dumper (@output);
    return (@output);
}