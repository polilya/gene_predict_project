#!/usr/bin/perl -w
use strict;
use MultiRandom;

my $EXTENSION = "align";

my $TRUE = 1;
my $FALSE = 0;

#This script takes as input a directory containing a set of FASTA files
#It splits each file into fragments

my $usage = "split_align.pl <align directory> <output directory> <fragment size>\n";

@ARGV == 3 || die $usage;

my ($align_dir, $out_dir, $frag_size) = @ARGV;

opendir (ALIGN_DIR, $align_dir);
my @align_files = grep(/\.align$/, readdir(ALIGN_DIR));
closedir (ALIGN_DIR);

@align_files = sort {$a cmp $b} @align_files;
my $mr;

for (my $i = 0; $i < @align_files; $i++) {
    print STDERR "$align_files[$i]\n";
    $align_files[$i] =~ /chr(\S+)\.align$/;

    #create a new directory for this file
    my $chr_num = $1; 
    mkdir $out_dir."/".$chr_num;
    my $this_seq_file = $align_dir."/".$align_files[$i];
    $mr = MultiRandom->new($this_seq_file);
    my $header = ">".join(" ",@{$mr->{SeqNames}});
    
    my $fragment_num = 1;
    mkdir ("$out_dir/$chr_num/fragment_$fragment_num");
    my $fragment_name = "$out_dir/$chr_num/fragment_$fragment_num/$fragment_num.$EXTENSION";
    my $fragment_start = 1;
    my $fragment_stop = $frag_size;
    while (writeFragment ($this_seq_file, $fragment_num, $fragment_name, $fragment_start, $fragment_stop, $header)) {
	$fragment_num++;
	mkdir ("$out_dir/$chr_num/fragment_$fragment_num");
	$fragment_name = "$out_dir/$chr_num/fragment_$fragment_num/$fragment_num.$EXTENSION";
	$fragment_start += $frag_size;
	$fragment_stop += $frag_size;
    }
}

sub writeFragment {
    my ($seq_file_name, $frag_num, $frag_name, $frag_start, $frag_stop, $header) = @_;
    my $last_fragment = $FALSE;
    my $frag_length;

    open (FRAG_OUT, ">$frag_name");
    print FRAG_OUT "$header\n";
    for (my $j=0; $j<$mr->{N}; $j++) {
	my $fragment_ref = $mr->getSeq($j, $frag_start, $frag_stop);
	print STDERR "Writing sequence $j\n";
	my $fragment = $$fragment_ref;
	if ($fragment =~ /(\S+)\n/) {
	    $last_fragment = $TRUE;
	    $fragment = $1;
	    print STDERR "Found end of chromosome, fragment length ".length($fragment)."\n";
	}
	print FRAG_OUT "$fragment\n";
	$frag_length = length($fragment);
	$frag_stop = $frag_start + $frag_length - 1;
    }

    print STDERR "fragment_".$frag_num."\t".($frag_start)."\t".($frag_start+$frag_length-1)."\n";
    system ("gzip --fast --force $frag_name"); 

    if (! $last_fragment) {
	return 1;
    }
    else {
	return 0;  #we hit end of file, tell main function to stop
    }
}
