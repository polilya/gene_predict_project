#!/usr/bin/perl -w
use strict;
use FastaRandom;

my $EXTENSION = "estseq";

my $TRUE = 1;
my $FALSE = 0;

#It splits each file into fragments

my $usage = "split_estseq.pl <estseq directory> <output directory> <fragment size>\n";

@ARGV == 3 || die $usage;

my ($estseq_dir, $out_dir, $frag_size) = @ARGV;

opendir (ESTSEQ_DIR, $estseq_dir);
my @estseq_files = grep(/\.estseq$/, readdir(ESTSEQ_DIR));
closedir (ESTSEQ_DIR);

@estseq_files = sort {$a cmp $b} @estseq_files;
my $fr;

for (my $i = 0; $i < @estseq_files; $i++) {
    print STDERR "$estseq_files[$i]\n";
    $estseq_files[$i] =~ /chr(\S+)\.estseq$/;

    #create a new directory for this file
    my $chr_num = $1; 
    mkdir $out_dir."/".$chr_num;

    my $this_seq_file = $estseq_dir."/".$estseq_files[$i];
    $fr = FastaRandom->new($this_seq_file);
    
    my $fragment_num = 1;
    mkdir ("$out_dir/$chr_num/fragment_$fragment_num");
    my $fragment_name = "$out_dir/$chr_num/fragment_$fragment_num/$fragment_num.$EXTENSION";
    my $fragment_start = 1;
    my $fragment_stop = $frag_size;
    while (writeFragment ($this_seq_file, $fragment_num, $fragment_name, $fragment_start, $fragment_stop, "chr$chr_num-$fragment_num")) {
	$fragment_num++;
	mkdir ("$out_dir/$chr_num/fragment_$fragment_num");
	$fragment_name = "$out_dir/$chr_num/fragment_$fragment_num/$fragment_num.$EXTENSION";
	$fragment_start += $frag_size;
	$fragment_stop += $frag_size;
    }
}

sub writeFragment {
    my ($seq_file_name, $frag_num, $frag_name, $frag_start, $frag_stop, $header) = @_;
    #print ("Writing fragment $frag_name, $frag_start-$frag_stop\n");

    my $last_fragment = $FALSE;
    my $frag_length;

    open (FRAG_OUT, ">$frag_name");
    if ($frag_stop > $fr->{LineLength}) {
	$frag_stop = $fr->{LineLength};
	$last_fragment = $TRUE;
    }
    my $fragment_ref = $fr->getSeq($frag_start, $frag_stop);
    my $fragment = $$fragment_ref;
    print FRAG_OUT "$fragment\n";
    $frag_length = length($fragment);
    $frag_stop = $frag_start + $frag_length - 1;

    print STDERR "fragment_".$frag_num."\t".($frag_start)."\t".($frag_start+$frag_length-1)."\n";
    system ("gzip --fast --force $frag_name"); 

    if (! $last_fragment) {
	return 1;
    }
    else {
	return 0;  #we hit end of file, tell main function to stop
    }
}
