# ExceS-A-An-Exon-Centric-Split-Aligner

ExceS-A is a fast, exon-centric spliced aligner and is included in the ExonMatchSolver pipeline.

please read:

http://www.bioinf.uni-leipzig.de/~henrike/index_pdf.pdf                                                                                                                        

additional Download:

•Cplex: https://www.ibm.com/de-de/products/ilog-cplex-optimization-studio

•MaxEntScan

For running the whole pipeline, following programs are required(linux only):

•perl including the following packages: Bio::Perl; Bio::SeqIO; Bio::Seq;Data::Dumper; Bio::DB::Fasta; Statistics::Standard_Normal; Bio::AlignIO; Getopt::Long;List::Util; FindBin, Cwd; DateTime; File::Slurp; List::MoreUtils; Array::Utils; Sort::Key::Natural; Bio::Tools::CodonTable

•BLAST+ tool suite (v 2.2.26)

Install the ExonMatchSover5.jar in a directory of your choice. The path to this directory may be specified by changing $Bin in ExonMatchSolver5.pl. Additional programs will be required for the following modes:

Fasta-mode:

•clustalw (v.2.0.12)

Alignment-mode:

•getORF (v EMBOSS:6.6.0.0)

•HMMER (v 3.1b1)
To get a more sensitive HMMER, it is necessary to copy its source code in a different directory and change the source code in /src/p7_pipeline.c
Outcommand the following two lines of code with “/*” (see also Cryp-togenomicon). The path to the more sensitive hmmer may be specified by changing $hmmer_path in the perl-script.

if (pli->ddef->nregions == 0) return eslOK;

/* score passed threshold but there’s no discrete domains here */

if (pli->ddef->nenvelopes == 0) return eslOK;

/* rarer: region was found, stochastic clustered, no envelopes found */

# command ExonMatchSolver:

perl ExonMatchSolver5.pl -i [input paralog file] -o [output file] -mode fasta -target [target genome] -paraMax [number] -query [query genome]



# Automation:

Automatic generation of exon-level HMMs for protein homology search

command:
perl Automation.pl --Paralogs [input paralog file] --TargetGenomes [list] --paraMax [number] --m [/path to current directory/Paralog_exon_alignments/] --manualCuration [yes/no]
