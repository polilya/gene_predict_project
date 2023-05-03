#!/usr/bin/perl -w
use strict;

my $usage = "$0 <species list>\n";

@ARGV > 0 || die $usage;

my $NUM_RANDOMS = 10000000;
my $AA_K = 0;
my $DNA_K = 6;
my $DNA_PAIR_K = 3;
my $DNA_KMER_ARRAY_K = 0;
my $DNA_KMER_PAIR_ARRAY_K = 0;

my $SVMS_ON = 1;

my $GC_DONOR_ON = 1;

my $EST_ON = 1;

my @implicit_sfs = ("Intron",
		    "CDS0",
		    "CDS1",
		    "CDS2");

my @dna_kmer_arrays = qw(StartCodon StopCodon DonorSplice DonorSpliceGC AcceptorSplice);
my @dna_kmer_array_lengths = (12, 6,  9, 9, 43);
my @dna_kmer_array_offsets = ( 8, 2,  5, 5, 41); 

my @svms;
my @svm_types;
my @svm_sample_rates;
my @svm_kernels;
my @svm_kernel_orders;
my @svm_bins;
my @svm_lengths;
my @svm_offsets; 

if ($GC_DONOR_ON) {
    @svms = qw(StartCodon StopCodon DonorSplice DonorSpliceGC AcceptorSplice);
    @svm_types = qw(Positional Positional Positional Positional Positional);
    @svm_sample_rates = (1.0, 1.0, 0.5, 1.0, 0.5);
    @svm_kernels = qw(Polynomial Polynomial Polynomial Polynomial Polynomial);
    @svm_kernel_orders = (2, 2, 2, 2, 2);
    @svm_bins = (5, 5, 10, 5, 10);
    @svm_lengths = (14, 9, 11, 11, 30);
    @svm_offsets = ( 8, 3,  3,  3, 27); 
}
else {
    @svms = qw(StartCodon StopCodon DonorSplice AcceptorSplice);
    @svm_types = qw(Positional Positional Positional Positional);
    @svm_sample_rates = (1.0, 1.0, 1.0, 1.0);
    @svm_kernels = qw(Polynomial Polynomial Polynomial Polynomial);
    @svm_kernel_orders = (2, 2, 2, 2);
    @svm_bins = (5, 5, 10, 10);
    @svm_lengths = (14, 9, 11, 30);
    @svm_offsets = ( 8, 3,  3, 27); 
}

#take list of species from the command line
my @species = @ARGV;
my @svm_species = @ARGV;

#create an array of many random numbers between -1 and 1
my @randoms;
for (my $i=0; $i<$NUM_RANDOMS; $i++) {
    $randoms[$i] = -0.02 * rand() + 0.01;
}

my $i = 0;

print "[States]\n\n";

printf "Intergenic\tIntergenic\tNo\t.\t%.6f\t.\n", $randoms[$i++];

printf "StartCodonCDS0\tCDS0\tYes\t+\t%.6f\t.\n", $randoms[$i++];
printf "StartCodonCDS1\tCDS1\tYes\t+\t%.6f\t.\n", $randoms[$i++];
printf "SingleCDS0\tCDS0\tYes\t+\t%.6f\t.\n", $randoms[$i++];
printf "SingleCDS1\tCDS1\tYes\t+\t%.6f\t.\n", $randoms[$i++];
printf "SingleCDS2\tCDS2\tYes\t+\t%.6f\tTAA,TAG,TGA\n", $randoms[$i++];
printf "InitialCDS0\tCDS0\tYes\t+\t%.6f\t.\n", $randoms[$i++];
printf "InitialCDS1\tCDS1\tYes\t+\t%.6f\t.\n", $randoms[$i++];
printf "InitialCDS2\tCDS2\tYes\t+\t%.6f\tTAA,TAG,TGA\n", $randoms[$i++];
printf "Intron0\tIntron\tNo\t+\t%.6f\t.\n", $randoms[$i++];
printf "Intron1\tIntron\tNo\t+\t%.6f\t.\n", $randoms[$i++];
printf "Intron1T\tIntron\tNo\t+\t%.6f\t.\n", $randoms[$i++];
printf "Intron2\tIntron\tNo\t+\t%.6f\t.\n", $randoms[$i++];
printf "Intron2TA\tIntron\tNo\t+\t%.6f\t.\n", $randoms[$i++];
printf "Intron2TG\tIntron\tNo\t+\t%.6f\t.\n", $randoms[$i++];
printf "InternalCDS0\tCDS0\tYes\t+\t%.6f\t.\n", $randoms[$i++];
printf "InternalCDS1\tCDS1\tYes\t+\t%.6f\t.\n", $randoms[$i++];
printf "InternalCDS2\tCDS2\tYes\t+\t%.6f\tTAA,TAG,TGA\n", $randoms[$i++];
printf "TerminalCDS0\tCDS0\tYes\t+\t%.6f\t.\n", $randoms[$i++];
printf "TerminalCDS1\tCDS1\tYes\t+\t%.6f\t.\n", $randoms[$i++];
printf "TerminalCDS2\tCDS2\tYes\t+\t%.6f\tTAA,TAG,TGA\n", $randoms[$i++];

printf "StartCodonCDS0-\tCDS0\tYes\t-\t%.6f\t.\n", $randoms[$i++];
printf "StartCodonCDS1-\tCDS1\tYes\t-\t%.6f\t.\n", $randoms[$i++];
printf "SingleCDS0-\tCDS0\tYes\t-\t%.6f\t.\n", $randoms[$i++];
printf "SingleCDS1-\tCDS1\tYes\t-\t%.6f\t.\n", $randoms[$i++];
printf "SingleCDS2-\tCDS2\tYes\t-\t%.6f\tTTA,CTA,TCA\n", $randoms[$i++];
printf "InitialCDS0-\tCDS0\tYes\t-\t%.6f\t.\n", $randoms[$i++];
printf "InitialCDS1-\tCDS1\tYes\t-\t%.6f\t.\n", $randoms[$i++];
printf "InitialCDS2-\tCDS2\tYes\t-\t%.6f\tTTA,CTA,TCA\n", $randoms[$i++];
printf "Intron0-\tIntron\tNo\t-\t%.6f\t.\n", $randoms[$i++];
printf "Intron0A-\tIntron\tNo\t-\t%.6f\t.\n", $randoms[$i++];
printf "Intron1-\tIntron\tNo\t-\t%.6f\t.\n", $randoms[$i++];
printf "Intron1TA-\tIntron\tNo\t-\t%.6f\t.\n", $randoms[$i++];
printf "Intron1CA-\tIntron\tNo\t-\t%.6f\t.\n", $randoms[$i++];
printf "Intron2-\tIntron\tNo\t-\t%.6f\t.\n", $randoms[$i++];
printf "InternalCDS0-\tCDS0\tYes\t-\t%.6f\t.\n", $randoms[$i++];
printf "InternalCDS1-\tCDS1\tYes\t-\t%.6f\t.\n", $randoms[$i++];
printf "InternalCDS2-\tCDS2\tYes\t-\t%.6f\tTTA,CTA,TCA\n", $randoms[$i++];
printf "TerminalCDS0-\tCDS0\tYes\t-\t%.6f\t.\n", $randoms[$i++];
printf "TerminalCDS1-\tCDS1\tYes\t-\t%.6f\t.\n", $randoms[$i++];
printf "TerminalCDS2-\tCDS2\tYes\t-\t%.6f\tTTA,CTA,TCA\n", $randoms[$i++];

print "\n\n";

print "[Transitions]\n\n";

printf "Intergenic\tIntergenic\t.\t\t%.6f\t.\t.\t.\t.\n", $randoms[$i++];

printf "Intergenic\tStartCodonCDS0\tStartCodon\t%.6f\t.\tATG\t.\t.\n", $randoms[$i++];
printf "StartCodonCDS0\tStartCodonCDS1\t.\t%.6f\t.\t.\t.\t.\n", $randoms[$i++];
printf "StartCodonCDS1\tSingleCDS2\t.\t%.6f\t.\t.\t.\t.\n", $randoms[$i++];
printf "StartCodonCDS1\tInitialCDS2\t.\t%.6f\t.\t.\t.\t.\n", $randoms[$i++];

printf "SingleCDS0\tSingleCDS1\t.\t\t%.6f\t.\t.\t.\t.\n", $randoms[$i++];
printf "SingleCDS1\tSingleCDS2\t.\t\t%.6f\t.\t.\t.\t.\n", $randoms[$i++];
printf "SingleCDS2\tSingleCDS0\t.\t\t%.6f\t.\t.\t.\t.\n", $randoms[$i++];
printf "SingleCDS2\tIntergenic\tStopCodon\t%.6f\t.\tTAA,TAG,TGA\t.\t.\n", $randoms[$i++];

printf "Intergenic\tInitialCDS0\tStartCodon\t%.6f\t.\tATG\t.\t.\n", $randoms[$i++];
printf "InitialCDS0\tInitialCDS1\t.\t\t%.6f\t.\t.\t.\t.\n", $randoms[$i++];
printf "InitialCDS1\tInitialCDS2\t.\t\t%.6f\t.\t.\t.\t.\n", $randoms[$i++];
printf "InitialCDS2\tInitialCDS0\t.\t\t%.6f\t.\t.\t.\t.\n", $randoms[$i++];

printf "InitialCDS0\tIntron1\t\tDonorSplice\t%.6f\t.\tGT\tT\t.\n", $randoms[$i++];
printf "InitialCDS0\tIntron1T\t\tDonorSplice\t%.6f\tT\tGT\t.\t.\n", $randoms[$i++];
printf "InitialCDS1\tIntron2\t\tDonorSplice\t%.6f\t.\tGT\tTA,TG\t.\n", $randoms[$i++];
printf "InitialCDS1\tIntron2TA\t\tDonorSplice\t%.6f\tTA\tGT\t.\t.\n", $randoms[$i++];
printf "InitialCDS1\tIntron2TG\t\tDonorSplice\t%.6f\tTG\tGT\t.\t.\n", $randoms[$i++];
printf "InitialCDS2\tIntron0\t\tDonorSplice\t%.6f\t.\tGT\t.\t.\n", $randoms[$i++];

if ($GC_DONOR_ON) {
    printf "InitialCDS0\tIntron1\t\tDonorSpliceGC\t%.6f\t.\tGC\tT\t.\n", $randoms[$i++];
    printf "InitialCDS0\tIntron1T\t\tDonorSpliceGC\t%.6f\tT\tGC\t.\t.\n", $randoms[$i++];
    printf "InitialCDS1\tIntron2\t\tDonorSpliceGC\t%.6f\t.\tGC\tTA,TG\t.\n", $randoms[$i++];
    printf "InitialCDS1\tIntron2TA\t\tDonorSpliceGC\t%.6f\tTA\tGC\t.\t.\n", $randoms[$i++];
    printf "InitialCDS1\tIntron2TG\t\tDonorSpliceGC\t%.6f\tTG\tGC\t.\t.\n", $randoms[$i++];
    printf "InitialCDS2\tIntron0\t\tDonorSpliceGC\t%.6f\t.\tGC\t.\t.\n", $randoms[$i++];
}

printf "Intron0\t\tIntron0\t\t.\t\t%.6f\t.\t.\t.\t.\n", $randoms[$i++];
printf "Intron1\t\tIntron1\t\t.\t\t%.6f\t.\t.\t.\t.\n", $randoms[$i++];
printf "Intron1T\t\tIntron1T\t\t.\t\t%.6f\t.\t.\t.\t.\n", $randoms[$i++];
printf "Intron2\t\tIntron2\t\t.\t\t%.6f\t.\t.\t.\t.\n", $randoms[$i++];
printf "Intron2TA\t\tIntron2TA\t\t.\t\t%.6f\t.\t.\t.\t.\n", $randoms[$i++];
printf "Intron2TG\t\tIntron2TG\t\t.\t\t%.6f\t.\t.\t.\t.\n", $randoms[$i++];
printf "Intron0\t\tInternalCDS0\tAcceptorSplice\t%.6f\tAG\t.\t.\t.\n", $randoms[$i++];
printf "Intron0\t\tTerminalCDS0\tAcceptorSplice\t%.6f\tAG\t.\t.\t.\n", $randoms[$i++];
printf "Intron1\t\tInternalCDS1\tAcceptorSplice\t%.6f\tAG\t.\t.\t.\n", $randoms[$i++];
printf "Intron1\t\tTerminalCDS1\tAcceptorSplice\t%.6f\tAG\t.\t.\t.\n", $randoms[$i++];
printf "Intron1T\t\tInternalCDS1\tAcceptorSplice\t%.6f\tAG\t.\t.\tAA,AG,GA\n", $randoms[$i++];
printf "Intron1T\t\tTerminalCDS1\tAcceptorSplice\t%.6f\tAG\t.\t.\tAA,AG,GA\n", $randoms[$i++];
printf "Intron2\t\tInternalCDS2\tAcceptorSplice\t%.6f\tAG\t.\t.\t.\n", $randoms[$i++];
printf "Intron2\t\tTerminalCDS2\tAcceptorSplice\t%.6f\tAG\t.\t.\t.\n", $randoms[$i++];
printf "Intron2TA\t\tInternalCDS2\tAcceptorSplice\t%.6f\tAG\t.\t.\tA,G\n", $randoms[$i++];
printf "Intron2TA\t\tTerminalCDS2\tAcceptorSplice\t%.6f\tAG\t.\t.\tA,G\n", $randoms[$i++];
printf "Intron2TG\t\tInternalCDS2\tAcceptorSplice\t%.6f\tAG\t.\t.\tA\n", $randoms[$i++];
printf "Intron2TG\t\tTerminalCDS2\tAcceptorSplice\t%.6f\tAG\t.\t.\tA\n", $randoms[$i++];

printf "InternalCDS0\tInternalCDS1\t.\t\t%.6f\t.\t.\t.\t.\n", $randoms[$i++];
printf "InternalCDS1\tInternalCDS2\t.\t\t%.6f\t.\t.\t.\t.\n", $randoms[$i++];
printf "InternalCDS2\tInternalCDS0\t.\t\t%.6f\t.\t.\t.\t.\n", $randoms[$i++];

printf "InternalCDS0\tIntron1\t\tDonorSplice\t%.6f\t.\tGT\tT\t.\n", $randoms[$i++];
printf "InternalCDS0\tIntron1T\t\tDonorSplice\t%.6f\tT\tGT\t.\t.\n", $randoms[$i++];
printf "InternalCDS1\tIntron2\t\tDonorSplice\t%.6f\t.\tGT\tTA,TG\t.\n", $randoms[$i++];
printf "InternalCDS1\tIntron2TA\t\tDonorSplice\t%.6f\tTA\tGT\t.\t.\n", $randoms[$i++];
printf "InternalCDS1\tIntron2TG\t\tDonorSplice\t%.6f\tTG\tGT\t.\t.\n", $randoms[$i++];
printf "InternalCDS2\tIntron0\t\tDonorSplice\t%.6f\t.\tGT\t.\t.\n", $randoms[$i++];

if ($GC_DONOR_ON) {
    printf "InternalCDS0\tIntron1\t\tDonorSpliceGC\t%.6f\t.\tGC\tT\t.\n", $randoms[$i++];
    printf "InternalCDS0\tIntron1T\t\tDonorSpliceGC\t%.6f\tT\tGC\t.\t.\n", $randoms[$i++];
    printf "InternalCDS1\tIntron2\t\tDonorSpliceGC\t%.6f\t.\tGC\tTA,TG\t.\n", $randoms[$i++];
    printf "InternalCDS1\tIntron2TA\t\tDonorSpliceGC\t%.6f\tTA\tGC\t.\t.\n", $randoms[$i++];
    printf "InternalCDS1\tIntron2TG\t\tDonorSpliceGC\t%.6f\tTG\tGC\t.\t.\n", $randoms[$i++];
    printf "InternalCDS2\tIntron0\t\tDonorSpliceGC\t%.6f\t.\tGC\t.\t.\n", $randoms[$i++];
}

printf "TerminalCDS0\tTerminalCDS1\t.\t\t%.6f\t.\t.\t.\t.\n", $randoms[$i++];
printf "TerminalCDS1\tTerminalCDS2\t.\t\t%.6f\t.\t.\t.\t.\n", $randoms[$i++];
printf "TerminalCDS2\tTerminalCDS0\t.\t\t%.6f\t.\t.\t.\t.\n", $randoms[$i++];
printf "TerminalCDS2\tIntergenic\tStopCodon\t%.6f\t.\tTAA,TAG,TGA\t.\t.\n", $randoms[$i++];

printf "Intergenic\tSingleCDS2-\tStopCodon\t%.6f\tTTA,CTA,TCA\t.\t.\t.\n", $randoms[$i++];
printf "SingleCDS2-\tSingleCDS1-\t.\t\t%.6f\t.\t.\t.\t.\n", $randoms[$i++];
printf "SingleCDS1-\tSingleCDS0-\t.\t\t%.6f\t.\t.\t.\t.\n", $randoms[$i++];
printf "SingleCDS0-\tSingleCDS2-\t.\t\t%.6f\t.\t.\t.\t.\n", $randoms[$i++];

printf "Intergenic\tTerminalCDS2-\tStopCodon\t%.6f\tTTA,CTA,TCA\t.\t.\t.\n", $randoms[$i++];
printf "TerminalCDS2-\tTerminalCDS1-\t.\t\t%.6f\t.\t.\t.\t.\n", $randoms[$i++];
printf "TerminalCDS1-\tTerminalCDS0-\t.\t\t%.6f\t.\t.\t.\t.\n", $randoms[$i++];
printf "TerminalCDS0-\tTerminalCDS2-\t.\t\t%.6f\t.\t.\t.\t.\n", $randoms[$i++];
printf "TerminalCDS2-\tIntron1-\tAcceptorSplice\t%.6f\t.\tCT\t.\t.\n", $randoms[$i++];
printf "TerminalCDS2-\tIntron1TA-\tAcceptorSplice\t%.6f\t.\tCT\tC,T\t.\n", $randoms[$i++];
printf "TerminalCDS2-\tIntron1CA-\tAcceptorSplice\t%.6f\t.\tCT\tT\t.\n", $randoms[$i++];
printf "TerminalCDS1-\tIntron0-\tAcceptorSplice\t%.6f\t.\tCT\t.\t.\n", $randoms[$i++];
printf "TerminalCDS1-\tIntron0A-\tAcceptorSplice\t%.6f\t.\tCT\tCT,TC,TT\t.\n", $randoms[$i++];
printf "TerminalCDS0-\tIntron2-\tAcceptorSplice\t%.6f\t.\tCT\t.\t.\n", $randoms[$i++];

printf "Intron2-\tIntron2-\t.\t\t%.6f\t.\t.\t.\t.\n", $randoms[$i++];
printf "Intron1-\tIntron1-\t.\t\t%.6f\t.\t.\t.\t.\n", $randoms[$i++];
printf "Intron1TA-\tIntron1TA-\t.\t\t%.6f\t.\t.\t.\t.\n", $randoms[$i++];
printf "Intron1CA-\tIntron1CA-\t.\t\t%.6f\t.\t.\t.\t.\n", $randoms[$i++];
printf "Intron0-\tIntron0-\t.\t\t%.6f\t.\t.\t.\t.\n", $randoms[$i++];
printf "Intron0A-\tIntron0A-\t.\t\t%.6f\t.\t.\t.\t.\n", $randoms[$i++];

printf "Intron2-\tInternalCDS2-\tDonorSplice\t%.6f\tAC\t.\t.\t.\n", $randoms[$i++];
printf "Intron2-\tInitialCDS2-\tDonorSplice\t%.6f\tAC\t.\t.\t.\n", $randoms[$i++];
printf "Intron1-\tInternalCDS1-\tDonorSplice\t%.6f\tAC\t.\t.\tTA,CA\n", $randoms[$i++];
printf "Intron1-\tInitialCDS1-\tDonorSplice\t%.6f\tAC\t.\t.\tTA,CA\n", $randoms[$i++];
printf "Intron1TA-\tInternalCDS1-\tDonorSplice\t%.6f\tAC\tTA\t.\t.\n", $randoms[$i++];
printf "Intron1TA-\tInitialCDS1-\tDonorSplice\t%.6f\tAC\tTA\t.\t.\n", $randoms[$i++];
printf "Intron1CA-\tInternalCDS1-\tDonorSplice\t%.6f\tAC\tCA\t.\t.\n", $randoms[$i++];
printf "Intron1CA-\tInitialCDS1-\tDonorSplice\t%.6f\tAC\tCA\t.\t.\n", $randoms[$i++];
printf "Intron0-\tInternalCDS0-\tDonorSplice\t%.6f\tAC\t.\t.\tA\n", $randoms[$i++];
printf "Intron0-\tInitialCDS0-\tDonorSplice\t%.6f\tAC\t.\t.\tA\n", $randoms[$i++];
printf "Intron0A-\tInternalCDS0-\tDonorSplice\t%.6f\tAC\tA\t.\t.\n", $randoms[$i++];
printf "Intron0A-\tInitialCDS0-\tDonorSplice\t%.6f\tAC\tA\t.\t.\n", $randoms[$i++];

if ($GC_DONOR_ON) {
    printf "Intron2-\tInternalCDS2-\tDonorSpliceGC\t%.6f\tGC\t.\t.\t.\n", $randoms[$i++];
    printf "Intron2-\tInitialCDS2-\tDonorSpliceGC\t%.6f\tGC\t.\t.\t.\n", $randoms[$i++];
    printf "Intron1-\tInternalCDS1-\tDonorSpliceGC\t%.6f\tGC\t.\t.\tTA,CA\n", $randoms[$i++];
    printf "Intron1-\tInitialCDS1-\tDonorSpliceGC\t%.6f\tGC\t.\t.\tTA,CA\n", $randoms[$i++];
    printf "Intron1TA-\tInternalCDS1-\tDonorSpliceGC\t%.6f\tGC\tTA\t.\t.\n", $randoms[$i++];
    printf "Intron1TA-\tInitialCDS1-\tDonorSpliceGC\t%.6f\tGC\tTA\t.\t.\n", $randoms[$i++];
    printf "Intron1CA-\tInternalCDS1-\tDonorSpliceGC\t%.6f\tGC\tCA\t.\t.\n", $randoms[$i++];
    printf "Intron1CA-\tInitialCDS1-\tDonorSpliceGC\t%.6f\tGC\tCA\t.\t.\n", $randoms[$i++];
    printf "Intron0-\tInternalCDS0-\tDonorSpliceGC\t%.6f\tGC\t.\t.\tA\n", $randoms[$i++];
    printf "Intron0-\tInitialCDS0-\tDonorSpliceGC\t%.6f\tGC\t.\t.\tA\n", $randoms[$i++];
    printf "Intron0A-\tInternalCDS0-\tDonorSpliceGC\t%.6f\tGC\tA\t.\t.\n", $randoms[$i++];
    printf "Intron0A-\tInitialCDS0-\tDonorSpliceGC\t%.6f\tGC\tA\t.\t.\n", $randoms[$i++];
}

printf "InternalCDS2-\tInternalCDS1-\t.\t\t%.6f\t.\t.\t.\t.\n", $randoms[$i++];
printf "InternalCDS1-\tInternalCDS0-\t.\t\t%.6f\t.\t.\t.\t.\n", $randoms[$i++];
printf "InternalCDS0-\tInternalCDS2-\t.\t\t%.6f\t.\t.\t.\t.\n", $randoms[$i++];
printf "InternalCDS2-\tIntron1-\tAcceptorSplice\t%.6f\t.\tCT\t.\t.\n", $randoms[$i++];
printf "InternalCDS2-\tIntron1TA-\tAcceptorSplice\t%.6f\t.\tCT\tC,T\t.\n", $randoms[$i++];
printf "InternalCDS2-\tIntron1CA-\tAcceptorSplice\t%.6f\t.\tCT\tT\t.\n", $randoms[$i++];
printf "InternalCDS1-\tIntron0-\tAcceptorSplice\t%.6f\t.\tCT\t.\t.\n", $randoms[$i++];
printf "InternalCDS1-\tIntron0A-\tAcceptorSplice\t%.6f\t.\tCT\tCT,TC,TT\t.\n", $randoms[$i++];
printf "InternalCDS0-\tIntron2-\tAcceptorSplice\t%.6f\t.\tCT\t.\t.\n", $randoms[$i++];

printf "InitialCDS2-\tInitialCDS1-\t.\t\t%.6f\t.\t.\t.\t.\n", $randoms[$i++];
printf "InitialCDS1-\tInitialCDS0-\t.\t\t%.6f\t.\t.\t.\t.\n", $randoms[$i++];
printf "InitialCDS0-\tInitialCDS2-\t.\t\t%.6f\t.\t.\t.\t.\n", $randoms[$i++];

printf "InitialCDS2-\tStartCodonCDS1-\t.\t%.6f\t.\t.\t.\t.\n", $randoms[$i++];
printf "SingleCDS2-\tStartCodonCDS1-\t.\t%.6f\t.\t.\t.\t.\n", $randoms[$i++];
printf "StartCodonCDS1-\tStartCodonCDS0-\t.\t%.6f\t.\t.\t.\t.\n", $randoms[$i++];
printf "StartCodonCDS0-\tIntergenic\tStartCodon\t%.6f\tCAT\t.\t.\t.\n", $randoms[$i++];

print "\n\n";

print "[Lengths]\n\n";

print "\n\n";

print "[MaskingFeatures]\n\n";
for (my $j=0; $j<@implicit_sfs; $j++) {
    printf "$implicit_sfs[$j]\n";
    printf "%.6f\n\n", $randoms[$i++]; 
}
print "\n\n";

my $num_dna_params = 4**$DNA_K + $DNA_K + 1;
print "[DNAKmerFeatures]\n\n";
if ($DNA_K > 0) {
    for (my $j=0; $j<@implicit_sfs; $j++) {
	printf "$implicit_sfs[$j]\t$DNA_K\t.\t$species[0]\tStandard\n";
	for (my $k=0; $k<$num_dna_params; $k++) {
	    printf "%.6f\n", $randoms[$i++];
	}
	printf "\n";
    }
}
print "\n\n";

my $num_dna_pair_params = 4**(2 * $DNA_PAIR_K) + (2**$DNA_PAIR_K - 1) + 3 + 1;
print "[DNAKmerPairFeatures]\n\n";
if ($DNA_PAIR_K > 0) {
    for (my $j=0; $j<@implicit_sfs; $j++) {
	for (my $a=0; $a == 0; $a++) {
	    for (my $b=$a+1; $b<@species; $b++) {
		printf "$implicit_sfs[$j]\t$DNA_PAIR_K\t.\t$species[$a]\t$species[$b]\tStandard\n";
		for (my $k=0; $k<$num_dna_pair_params; $k++) {
		    printf "%.6f\n", $randoms[$i++];
		}
		printf "\n";
	    }
	}
    }
}
print "\n\n";

print "[DNAKmerArrayFeatures]\n\n";
if ($DNA_KMER_ARRAY_K > 0) {
    for (my $j=0; $j<@dna_kmer_arrays; $j++) {
	printf "$dna_kmer_arrays[$j]\t$DNA_KMER_ARRAY_K\t$dna_kmer_array_lengths[$j]\t$dna_kmer_array_offsets[$j]\t$species[0]\tStandard\n";
	my $num_params = (4**$DNA_KMER_ARRAY_K + $DNA_KMER_ARRAY_K + 1) * $dna_kmer_array_lengths[$j]; 
	for (my $k=0; $k<$num_params; $k++) {
	    printf "%.6f\n", $randoms[$i++];
	}
	print "\n";
    }
}
print "\n\n";

print "[DNAKmerPairArrayFeatures]\n\n";
if ($DNA_KMER_PAIR_ARRAY_K > 0) {
    for (my $j=0; $j<@dna_kmer_arrays; $j++) {
	for (my $a=0; $a == 0; $a++) {
	    for (my $b=$a+1; $b<@species; $b++) {
		printf "$dna_kmer_arrays[$j]\t$DNA_KMER_PAIR_ARRAY_K\t$dna_kmer_array_lengths[$j]\t$dna_kmer_array_offsets[$j]\t$species[$a]\t$species[$b]\tStandard\n";
		my $num_params = $dna_kmer_array_lengths[$j] * (4**(2 * $DNA_KMER_PAIR_ARRAY_K) 
								+ (2**$DNA_KMER_PAIR_ARRAY_K - 1)
								+ 3 + 1);
		#my $num_params = $dna_kmer_array_lengths[$j] * 7**(2*$DNA_KMER_PAIR_ARRAY_K);
		for (my $k=0; $k<$num_params; $k++) {
		    printf "%.6f\n", $randoms[$i++];
		}
		print "\n";
	    }
	}
    }
}

print "[SVMFeatures]\n\n";
if ($SVMS_ON) {
for (my $j=0; $j<@svms; $j++) {
    print "$svms[$j]\t$svm_types[$j]\t$svm_kernels[$j]\t$svm_kernel_orders[$j]\t$svm_bins[$j]\t$svm_lengths[$j]\t$svm_offsets[$j]\t$svm_sample_rates[$j]\t";
    for (my $k=0; $k<@svm_species; $k++) {
	print "$svm_species[$k]";
	if ($k < @svm_species - 1) {
	    print ",";
	}
    }
    print "\n";
    for (my $k=0; $k<$svm_bins[$j]; $k++) {
	printf "%.6f\t0\n", $randoms[$i++];
    }
    print "\n";
}
}

print "[ESTPositionFeatures]\n\n";
if ($EST_ON) {
    for (my $j=0; $j<@implicit_sfs; $j++) {
	print "$implicit_sfs[$j]\n";
	for (my $k=0; $k<5; $k++) {
	    printf "%.6f\n", $randoms[$i++];
	}
	print "\n";
    }
    print "\n";
}

print "[ESTTransitionFeatures]\n\n";
if ($EST_ON) {
    for (my $j=0; $j<@dna_kmer_arrays; $j++) {
	print "$dna_kmer_arrays[$j]\n";
	for (my $k=0; $k<(5*5); $k++) {
	    printf "%.6f\n", $randoms[$i++];
	}
	print "\n";
    }
    print "\n";
}
