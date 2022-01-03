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

#use File::Compare::compare_text;
#------------------------------------------------------------------------
#Definitions
#
my $pathtotargetgenomeIn = '/scratch/franziska/genomes/worms';
#Path to ExonMatchSolver.pl		
my $EMSpipeline = "/scratch/franziska/ExonMatchSolver/ExonMatchSolver5/ExonMatchSolver5_fasta_4_0_0.pl";	

my $space = "/scratch/franziska/genomes/";


my $genome_query;
my $genome_target;
my $mode;
my $hmmModels;

		
my $paraMax = "";
my $querygenomeIn = "";
my $paralogsIn = "";
my $targetgenomesIn = "";

my $distance_first = 10000; 	#Can be defined by the user in dependence of the intron size expected for the first (few) exons.   
my $distance_last = 10000;

my $core = 8;	
my $genomeIn;
my $querysequenceIn;
my $blastoutIn;
my $mode;

GetOptions("genomeQuery=s" => \$genome_query, "genomeTarget=s" => \$genome_target, "querysequence=s" => \$querysequenceIn, "mode=s" => \$mode, "hmmModels=s" => \$hmmModels);

#------------------------------------------------------------------------
#Main


my @buildquery = "formatdb -i ${genome_query} -p F -o";
unless (-e "${genome_query}.nhr" && -e "${genome_query}.nin" && -e "${genome_query}.nsd" && -e "${genome_query}.nsi" && -e "${genome_query}.nsq"){
	system (@buildquery) == 0 or die "Formatdb with query genome failed. $!\n";
}

my @call = "blastall -p tblastn -d $genome_query -i $querysequenceIn -a $core -m 8 -G 12 -E 1 -o ${ID}.blastout";
system (@call) == 0 or die "blastall against query genome failed. $!\n";

my $ID = $querysequence;
$ID =~ s/.fa//;

my $locus = `sort -k11,11g ${ID}.blastout| awk '{if (\$3>90) print \$2}'| head -n1`;
chomp $locus;
$loci{$ID} = $locus;

print Dumper ($locus);

qx(perl $ExceS_A_3 --genome $genome_query --locus $locus --ID $ID --querysequence $querysequenceIn);

if ($mode eq "alignment"){
    # translate the genomic sequence into 6 ORFs, stored in ${space}${myOutFolder}_ORF.aas (amino-acid sequence) with getorf
    chdir $space or die "Couldn't access $space\n";
    system "awk '{print \$1}' $genome_target > tmp1";
    unless (-e "${space}${ID}_ORF.aas" && -e $genome_target) {
	my $ret = system "getorf -minsize 9 -reverse Y -stdout Y -sequence tmp1 -outseq stdout|sed 's/ - /-/g'| awk '{print \$1, \$3, \$2}' > tmp2";
	if ($ret != 0){
	    die;
	};
	open (TMP, "< tmp2") or die "Couldn't open tmp2, $!";
	open (DEST, ">> ${space}${ID}_ORF.aas") or die "Couldn't open ${space}${ID}_ORF.aas, $!";
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

    mkdir "SearchTarget/";  
    my @files = glob "$hmmModels/*";
    print "Progress of hmm search against the target genome:\n";
    foreach my $file(@files){
		$_ = $file;
		s/hmmModels\///;
		s/\.model//;
		$file = $_;
		print "$file\n";
		system "hmmsearch --cpu $core --domtblout SearchTarget/${file}_domtableout -E 0.01 -o SearchTarget/${file}_out hmmModels/${file}.model ${space}${ID}_ORF.aas >> SearchTarget/Target_vs_query.hmmsearchout";
	};

    print "Hmmsearch of query models against target genome done\n";

	$hmm = "SearchTarget/hmmsearch_tableout";
    system "cat SearchTarget/*_domtableout| sed \'s/-\\([0-9]\\+\\)]/_\\1]/g\' > $hmm";

	if (-e $hmm){
	    open (IN, "< $hmm") or die "Couldn't open $hmm, $!";
	    open (HITS, "> ExonMatchSolver1.input");
	    print HITS "#TargetName\tStrandinformation\tModel\tE-Value\tBitscore\tBias\tParalogID\tExonID\n";
	    while (<IN>){
			unless ($_ =~ /#/) {
			    chomp;
			    my @hits = split /[\t\s]+/, $_;   				
			    my @name = split /_/, $hits[3]; 
			    $hits[30] = $name[0];
			    $hits[31] = $name[1];	
			    $hits[31] =~ s/exon//g;

			    print HITS "$hits[0]\t$hits[22]\t$hits[3]\t$hits[6]\t$hits[7]\t$hits[8]\t$hits[30]\t$hits[31]\n";
			}
	    }
	    close HITS;
	    close IN;
	}
}

elsif($mode eq "fasta"){

	my @buildtarget = "formatdb -i ${genome_target} -p F -o";
	unless (-e "${genome_target}.nhr" && -e "${genome_target}.nin" && -e "${genome_target}.nsd" && -e "${genome_target}.nsi" && -e "${genome_target}.nsq"){
	system (@buildtarget) == 0 or die "Formatdb with target genome failed. $!. You need a permission to write files in the target genome directory\n";
	};

	my $blast = "SearchTarget/Target_vs_query.blastout";
	my @call = "blastall -p tblastn -d $genome_target -i HomologousExons/paralog_exon_alignment_final.fa -e 0.01 -a $core -m8 -G 11 -E 1 -C F > $blast";
	system (@call) == 0 or die "blastall against target genome failed. $!\n";

	@buildtarget = "formatdb -o -i $genome_target -p F";
	unless (-e "${genome_target}.nhr" && -e "${genome_target}.nin" && -e "${genome_target}.nsd" && -e "${genome_target}.nsi" && -e "${genome_target}.nsq"){
	system (@buildtarget) == 0 or die "Formatdb with target genome failed. $!\n";
	};
}

system("perl $ExceS_A_2 --contig $contig --query ${locus{$contig}[0]} --strand $strand --start ${locus{$contig}[2]} --end ${locus{$contig}[3]}");
