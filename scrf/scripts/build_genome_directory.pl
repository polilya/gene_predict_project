#!/usr/bin/perl -w
use strict;

use File::Path;

my $UCSC_URL="http://hgdownload.cse.ucsc.edu/goldenPath";
my $TRAIN_FRAGMENT_SIZE = 1000000;
my $PREDICT_FRAGMENT_SIZE = 5000000;

my $usage = "$0 <genome name> <alignment directories>\n";

@ARGV >= 2 || die $usage;

my $genome_name = $ARGV[0];

#create genome directory
mkdir $genome_name;

#remove any old subdirectories we're about to replace
rmtree "$genome_name/chr_seq";
rmtree "$genome_name/align";
rmtree "$genome_name/estseq";
rmtree "$genome_name/genes";
rmtree "$genome_name/split";


#download and extract genomic sequence
mkdir "$genome_name/chr_seq";

my $command = "wget -P $genome_name/chr_seq $UCSC_URL/$genome_name/bigZips/chromFa.tar.gz";
print STDERR "*** Downloading genomic sequence: $command\n";
system($command);
if (! -e "$genome_name/chr_seq/chromFa.tar.gz") {
    #try .zip
    my $command = "wget -P $genome_name/chr_seq $UCSC_URL/$genome_name/bigZips/chromFa.zip";
    print STDERR "*** Downloading genomic sequence: $command\n";
    system($command);

    if (! -e "$genome_name/chr_seq/chromFa.zip") {
	fatal_error("Could not download chromosome sequence");
    }

    $command = "unzip $genome_name/chr_seq/chromFa.zip -d $genome_name/chr_seq";
    print STDERR "*** Extracting genomic sequence: $command\n";
    system($command);
}
else {
    $command = "tar -zxf $genome_name/chr_seq/chromFa.tar.gz -C $genome_name/chr_seq";
    print STDERR "*** Extracting genomic sequence: $command\n";
    system($command);
}


#download and extract RepeatMasker output
$command = "wget -P $genome_name/chr_seq $UCSC_URL/$genome_name/bigZips/chromOut.tar.gz";
print STDERR "*** Downloading RepeatMasker output: $command";
system($command);

if (! -e "$genome_name/chr_seq/chromOut.tar.gz") {
    #try .zip
    $command = "wget -P $genome_name/chr_seq $UCSC_URL/$genome_name/bigZips/chromOut.zip";
    print STDERR "*** Downloading RepeatMasker output: $command";
    system($command);

    if (! -e "$genome_name/chr_seq/chromOut.zip") {
	fatal_error("Could not download RepeaterMasker output");
    }

    $command = "unzip $genome_name/chr_seq/chromOut.zip -d $genome_name/chr_seq";
    print STDERR "*** Extracting RepeatMasker output: $command\n";
    system($command);
}
else {
    $command = "tar -zxf $genome_name/chr_seq/chromOut.tar.gz -C $genome_name/chr_seq";
    print STDERR "*** Extracting RepeatMasker output: $command\n";
    system($command);
}


#move everything to the top level of the chr_seq directory
opendir(CHR_SEQ, "$genome_name/chr_seq");
my @chr_dir_files = readdir(CHR_SEQ);
closedir(CHR_SEQ);
for (my $i=0; $i<@chr_dir_files; $i++) {
    if (!($chr_dir_files[$i] =~ /^\./) && -d "$genome_name/chr_seq/$chr_dir_files[$i]") {
	system ("mv $genome_name/chr_seq/$chr_dir_files[$i]/* $genome_name/chr_seq");
	system ("rmdir $genome_name/chr_seq/$chr_dir_files[$i]");
    }
}

#remove non-standard chromosome files
system ("rm $genome_name/chr_seq/chrM*");
system ("rm $genome_name/chr_seq/chrU*");
system ("rm $genome_name/chr_seq/chrNA*");
system ("rm $genome_name/chr_seq/chr*random*");
system ("rm $genome_name/chr_seq/chr*h*");

#build list of chromosomes
my @chromosomes;
opendir(CHR_SEQ, "$genome_name/chr_seq");
@chr_dir_files = readdir(CHR_SEQ);
closedir(CHR_SEQ);
for (my $i=0; $i<@chr_dir_files; $i++) {
    if ($chr_dir_files[$i] =~ /^chr(.+)\.fa$/) {
	push (@chromosomes, $1);
    }
}


#download genomic alignments
my $pairwise = 0;  # do we have one or more pairwise alignments (as opposed to a single multiple alignment)?
my @informants;
mkdir "$genome_name/align";
for (my $i=1; $i<@ARGV; $i++) {
    if ($ARGV[$i] =~ /^vs(.+)/) {
	$pairwise = 1;  #pairwise alignments
	mkdir "$genome_name/align/$ARGV[$i]";
	my $aligned_species = $1;
	substr($aligned_species, 0, 1, lc(substr($aligned_species, 0, 1)));
	push(@informants, $aligned_species);
	print STDERR "*** Downloading genomic alignments from $aligned_species\n";
	for (my $j=0; $j<@chromosomes; $j++) {
	    $command = "wget -P $genome_name/align/$ARGV[$i] $UCSC_URL/$genome_name/$ARGV[$i]/axtNet/chr$chromosomes[$j].$genome_name.$aligned_species.net.axt.gz";
	    print STDERR "\t$command\n";
	    system($command);
	    if (! -e "$genome_name/align/$ARGV[$i]/chr$chromosomes[$j].$genome_name.$aligned_species.net.axt.gz") {
		fatal_error("Could not download genomic alignment file");
	    }
	}
	$command = "gunzip $genome_name/align/$ARGV[$i]/*gz";
	print STDERR "*** Extracting genomic alignments for $aligned_species: $command\n";
	system($command);
    }
    else {
	# multiple alignment
	if (@ARGV != 2) {
	    fatal_error("At most one multiple alignment directory can be specified");
	}
	if (! ($ARGV[$i] =~ /^(\S+)\:(\S+)$/)) {
	    fatal_error("Wrong format for multiple alignment argument.  Should be \"directory:informant1,informant2,informant3,...,informantN\"");
	}
	my $alignment_dir = $1;
	my $informant_list = $2;
	@informants = split(/,/, $informant_list);
	print STDERR "Using following informants from multiple alignment: ";
	for (my $i=0; $i<@informants; $i++) {
	    print STDERR "$informants[$i] ";
	}
	print STDERR "\n";
	print STDERR "*** Downloading multiple alignment of genomic sequences\n";
	for (my $j=0; $j<@chromosomes; $j++) {
	    $command = "wget -P $genome_name/align $UCSC_URL/$genome_name/$alignment_dir/chr$chromosomes[$j].maf.gz";
	    print STDERR "\t$command\n";
	    system($command);
	    if (! -e "$genome_name/align/chr$chromosomes[$j].maf.gz") {
		# Try maf subdirectory
		$command = "wget -P $genome_name/align $UCSC_URL/$genome_name/$alignment_dir/maf/chr$chromosomes[$j].maf.gz";
		print STDERR "\t$command\n";
		system($command);
		if (! -e "$genome_name/align/chr$chromosomes[$j].maf.gz") {
		    fatal_error("Could not download genomic alignment file");
		}
	    }
	}
	$command = "gunzip $genome_name/align/*gz";
	print STDERR "*** Extracting multiple alignment: $command\n";
	system($command);
    }
}
my $species_list = "$genome_name ";
for (my $i=0; $i<@informants; $i++) {
    $species_list .= "$informants[$i] ";
}

#download EST alignments
mkdir "$genome_name/estseq";
$command = "wget -P $genome_name/estseq $UCSC_URL/$genome_name/database/all_est.txt.gz";
print STDERR "*** Downloading EST alignments: $command\n";
system ($command);
if (! -e "$genome_name/estseq/all_est.txt.gz") {
    fatal_error("Could not download EST alignments");
}
$command = "gunzip $genome_name/estseq/all_est.txt.gz";
print STDERR "*** Extracting EST alignments: $command\n";
system ($command);


#download gene annotations
my $mrnas = 0;  # were we forced to download GenBank mRNAs for our set of known genes?
mkdir "$genome_name/genes";
#first try for CCDS
$command = "wget -P $genome_name/genes $UCSC_URL/$genome_name/database/ccdsGene.txt.gz";
print STDERR "*** Trying CCDS genes: $command\n";
system ($command);
if (! -e "$genome_name/genes/ccdsGene.txt.gz") {
    #next try MGC
    #$command = "wget -P $genome_name/genes $UCSC_URL/$genome_name/database/mgcGenes.txt.gz";
    #print STDERR "*** Trying MGC genes: $command\n";
    #system ($command);
    if (! -e "$genome_name/genes/mgcGenes.txt.gz") {
	#next try RefSeq
	$command = "wget -P $genome_name/genes $UCSC_URL/$genome_name/database/refGene.txt.gz";
	print STDERR "*** Trying RefSeq genes: $command\n";
	system ($command);
	if (! -e "$genome_name/genes/refGene.txt.gz") {
	    #finally try mRNAs from GenBank
	    $mrnas = 1;
	    $command = "wget -P $genome_name/genes $UCSC_URL/$genome_name/database/all_mrna.txt.gz";
	    print STDERR "*** Trying GenBank mRNAs: $command\n";
	    system ($command);
	    if (! -e "$genome_name/genes/all_mrna.txt.gz") {
		fatal_error("Could not find any gene annotations\n");
	    }
	}
    }
}
$command = "gunzip $genome_name/genes/*gz";
print STDERR "*** Extracting gene annotations: $command\n";
system ($command);



#mask repeats, except for low-complexity and simple repeats
print STDERR "*** Masking repeats\n";
for (my $i=0; $i<@chromosomes; $i++) {
    $command = "mask_from_rm_output $genome_name/chr_seq/chr$chromosomes[$i].fa $genome_name/chr_seq/chr$chromosomes[$i].fa.out > $genome_name/chr_seq/chr$chromosomes[$i].fa.masked\n";
    #$command = "cp $genome_name/chr_seq/chr$chromosomes[$i].fa $genome_name/chr_seq/chr$chromosomes[$i].fa.masked";
    print STDERR "\t$command\n";
    system($command);
}


#create .align files
print STDERR " *** Creating .align files\n";
if ($pairwise) {
    for (my $i=0; $i<@chromosomes; $i++) {
	$command = "axts_to_align $genome_name/chr_seq/chr$chromosomes[$i].fa.masked ";
	for (my $j=0; $j<@informants; $j++) {
	    $command .= "$genome_name/align/".$ARGV[$j+1]."/chr$chromosomes[$i].$genome_name.$informants[$j].net.axt ";
	}
	$command .= "$species_list > $genome_name/align/chr$chromosomes[$i].align";
	print STDERR "\t$command\n";
	system($command);
	if (-s "$genome_name/align/chr$chromosomes[$i].align" < 1000) {
	    fatal_error("Size of .align file extremely small, looks like there was an error in creating it");
	}
    }
}
else {
    for (my $i=0; $i<@chromosomes; $i++) {
	$command = "maf_to_align $genome_name/chr_seq/chr$chromosomes[$i].fa.masked $genome_name/align/chr$chromosomes[$i].maf $species_list > $genome_name/align/chr$chromosomes[$i].align";
	print STDERR "\t$command\n";
	system($command);
	if (-s "$genome_name/align/chr$chromosomes[$i].align" < 1000) {
	    fatal_error("Size of .align file extremely small, looks like there was an error in creating it");
	}
    }
}

#create .estseq files
print STDERR "*** Creating EST sequence\n";

#separate by chromosome
for (my $i=0; $i<@chromosomes; $i++) {
    $command = "cat $genome_name/estseq/all_est.txt | grep -E \"chr$chromosomes[$i]\t\" > $genome_name/estseq/chr$chromosomes[$i]_est.txt";
    print STDERR "\t$command\n";
    system($command);
}

#run make_estseq
for (my $i=0; $i<@chromosomes; $i++) {
    $command = "make_estseq $genome_name/chr_seq/chr$chromosomes[$i].fa $genome_name/estseq/chr$chromosomes[$i]_est.txt > $genome_name/estseq/chr$chromosomes[$i].estseq";
    print STDERR "\t$command\n";
    system($command);
}



#process gene annotations
if ($mrnas) {
    fatal_error("Processing mRNA annotations not yet supported");
}
else {
    $command = "ucsc2gtf.rb $genome_name/genes/* > $genome_name/genes/genes.gtf";
    print STDERR "*** Converting UCSC annotations to GTF: $command\n";
    system($command);
}

#separate GTF by chromosome
for (my $i=0; $i<@chromosomes; $i++) {
    $command = "cat $genome_name/genes/genes.gtf | grep -E \"chr$chromosomes[$i]\t\" > $genome_name/genes/chr$chromosomes[$i].gtf.raw";
    print STDERR "\t$command\n";
    system($command);
}
unlink "$genome_name/genes/genes.gtf";

#label transcripts by overlap
for (my $i=0; $i<@chromosomes; $i++) {
    $command = "label_transcripts_by_overlap.pl $genome_name/genes/chr$chromosomes[$i].gtf.raw > $genome_name/genes/chr$chromosomes[$i].gtf";
    print STDERR "\t$command\n";
    system($command);
    unlink "$genome_name/genes/chr$chromosomes[$i].gtf.raw";
}

#clean annotations
mkdir "$genome_name/genes/clean";
$command = "clean_gtf.pl $genome_name/genes $genome_name/chr_seq $genome_name/genes/clean";
print STDERR "*** Cleaning annotations: $command\n";
system($command);


#get alignment statistics
#print STDERR " *** Computing alignment statistics\n";
#for (my $i=0; $i<@chromosomes; $i++) {
#    $command = "get_align_stats.pl $genome_name/align/chr$chromosomes[$i].align $genome_name/genes/clean/chr$chromosomes[$i].gtf 0.01 > $genome_name/align/chr$chromosomes[$i].stats.txt\n";
#    print STDERR "\t$command\n";
#    system($command);
#}


#split input files and annotations into fragments
mkdir "$genome_name/split";

$command = "split_align.pl $genome_name/align $genome_name/split $TRAIN_FRAGMENT_SIZE";
print STDERR "*** Splitting alignment into fragments: $command\n";
system($command);

$command = "split_estseq.pl $genome_name/estseq $genome_name/split $TRAIN_FRAGMENT_SIZE";
print STDERR "*** Splitting EST sequence into fragments: $command\n";
system($command);

$command = "split_gtf.pl $genome_name/genes/clean $genome_name/split $TRAIN_FRAGMENT_SIZE";
print STDERR "*** Splitting annotations into fragments: $command\n";
system($command);


#also create another split directory with longer fragments for prediction
mkdir "$genome_name/split_long";

$command = "split_align.pl $genome_name/align $genome_name/split_long $PREDICT_FRAGMENT_SIZE";
print STDERR "*** Splitting alignment into fragments: $command\n";
system($command);

$command = "split_estseq.pl $genome_name/estseq $genome_name/split_long $PREDICT_FRAGMENT_SIZE";
print STDERR "*** Splitting EST sequence into fragments: $command\n";
system($command);

$command = "split_gtf.pl $genome_name/genes/clean $genome_name/split_long $PREDICT_FRAGMENT_SIZE";
print STDERR "*** Splitting annotations into fragments: $command\n";
system($command);


#create list of files for CRF training, SVM training, and holdout
$command = "make_training_testing_lists.pl $genome_name/split 0";
print STDERR "*** Dividing fragments into training and holdout sets: $command\n";
system($command);

$command = "make_training_testing_lists.pl $genome_name/split_long 0";
print STDERR "*** Dividing fragments into training and holdout sets: $command\n";
system($command);


#generate paramter file
$command = "generate_parameters.pl $species_list > $genome_name/split/random_parameters.txt";
print STDERR "*** Generating parameter file: $command\n";
system($command);

#tar split directories for transport to cluster
$command = "tar -cf $genome_name.split.tar $genome_name/split $genome_name/split_long";
print STDERR "*** Archiving split directory: $command\n";
system($command);


sub fatal_error {
    print "Fatal error: $_[0]\n";
    exit;
}
