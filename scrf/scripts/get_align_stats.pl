#!/usr/bin/perl -w
use strict;
use GTF_Parser_UTR_SU;
use MultiRandom;

my $usage = "get_align_stats.pl <align file> <GTF file> <noncoding sampling rate>\n";

@ARGV == 3 || die $usage;

my ($align_file, $gtf_file, $nc_sampling_rate) = @ARGV;

my $mr = MultiRandom->new($align_file);
my $n = $mr->{N};

my $CDS0 = 0;
my $CDS1 = 1;
my $CDS2 = 2;
my $NC = 3;
my @functional_types = qw (CDS0 CDS1 CDS2 NC); 

my $TOTAL = 0;
my $ALIGNED = 1;
my $MATCH = 2;

my @counts;
for (my $seq=0; $seq<$n-1; $seq++) {
    for (my $i=0; $i<4; $i++) {
	for (my $j=0; $j<3; $j++) {
	    $counts[$seq][$i][$j] = 0;
	}
    }
}

my @unsorted_genes = @{GTF_Parser_UTR_SU::parse_gtf($gtf_file)};
#sort genes by start coord
my @genes = sort by_start @unsorted_genes;

for (my $i=0; $i<@genes; $i++) {
    print STDERR "Gene $i\n";
    if ($i % 50 == 0 && $i>0) {
	outputStats();
    }

    my @this_gene = @{$genes[$i]};
    
    if ($i != @genes-1 && rand() < $nc_sampling_rate) {
	#get intergenic counts
	my $this_gene_stop = $this_gene[@this_gene - 1][$STOP];
	my $next_gene_start = $genes[$i+1][0][$START];
	#countNC ($this_gene_stop+7, $next_gene_start-7, "+");
    }
    
    for (my $j=0; $j<@this_gene; $j++) {
	my $type = $this_gene[$j][$TYPE];
	my $strand = $this_gene[$j][$STRAND];
	my $start = $this_gene[$j][$START];
	my $stop = $this_gene[$j][$STOP];    
	my $frame = $this_gene[$j][$FRAME];
	
	if (! (defined $genes[$i][$j][$TYPE])) {
	    print STDERR "Type undefined for $i, $j\n";
	}
	
	if ($type == $EINITU || $type == $EINITS) {
	    if ($strand eq "+") {
		countCDS ($frame, $start+6, $stop-3, $strand);
	    }
	    elsif ($strand eq "-") {
		countCDS ($frame, $start+3, $stop-6, $strand);
	    }
	}
	elsif ($type == $EXON) {
	    if ($strand eq "+") {
		countCDS ($frame, $start-3, $stop-3, $strand);
	    }
	    elsif ($strand eq "-") {
		countCDS ($frame, $start+3, $stop-3, $strand);
	    }
	}
	elsif ($type == $ETERM) {
	    if ($strand eq "+") {
		countCDS ($frame, $start+3, $stop, $strand);		   
	    }
	    elsif ($strand eq "-") {
		countCDS ($frame, $start, $stop-3, $strand);
	    }
	}
	elsif ($type == $ESNGL) {
	    if ($strand eq "+") {
		countCDS ($frame, $start+6, $stop, $strand);
	    }
	    elsif ($strand eq "-") {
		countCDS ($frame, $start, $stop-6, $strand);
	    }
	}
	elsif ($type == $INTRON && rand() < $nc_sampling_rate) {
	    if ($strand eq "+") {
		countNC ($start+6, $stop-39, $strand);
	    }
	    elsif ($strand eq "-") {
		countNC ($start+39, $stop-6, $strand);
	    }
	}
    }
}

outputStats();


sub outputStats {
    print "\t\t\t\tCDS1\tCDS2\tCDS3\tIntron\n";
    for (my $seq=0; $seq<$mr->{N}-1; $seq++) {
	my $seq_name = $mr->{SeqNames}->[$seq+1];
	print "$seq_name    % aligned\t";
	for (my $i=0; $i<4; $i++) {
	    my $percent_aligned = 100 * $counts[$seq][$i][1] / $counts[$seq][$i][0];
	    printf ("\t%.2lf", $percent_aligned);
	}
	print "\n";
	print "$seq_name    % identity\t";
	for (my $i=0; $i<4; $i++) {
	    if ($counts[$seq][$i][1] > 0) {
		my $percent_id = 100 * $counts[$seq][$i][2] / $counts[$seq][$i][1];
		printf ("\t%.2lf", $percent_id);
	    }
	}
	print "\n\n";
    }
}

sub countNC {
    my ($start, $stop, $strand) = @_;

    if ($stop < $start) {
	return;
    }

    my @seq_arr;

    if ($strand eq "+") {
	for (my $i=0; $i<$n; $i++) {
	    my $seq_ref = $mr->getSeq($i, $start, $stop);
	    my @seq_split = split (//, $$seq_ref);
	    $seq_arr[$i] = \@seq_split;
	}
    }
    elsif ($strand eq "-") {
	for (my $i=0; $i<$n; $i++) {
	    my $seq_ref = $mr->getRevCompSeq($i, $start, $stop);
	    my @seq_split = split (//, $$seq_ref);
	    $seq_arr[$i] = \@seq_split;
	}
    }
    
    for (my $i=0; $i<@{$seq_arr[0]}; $i++) {
	for (my $j=1; $j<$n; $j++) {
	    $counts[$j-1][$NC][$TOTAL]++;
	    if ($seq_arr[$j][$i] ne ".") {
		$counts[$j-1][$NC][$ALIGNED]++;
	    }
	    if ($seq_arr[$j][$i] eq $seq_arr[0][$i]) {
		$counts[$j-1][$NC][$MATCH]++;
	    }
	}
    }
}

sub countCDS {
    my ($frame, $start, $stop, $strand) = @_;

    if ($stop < $start) {
	return;
    }

    my @seq_arr;
    
    if ($strand eq "+") {
	for (my $i=0; $i<$n; $i++) {
	    my $seq_ref = $mr->getSeq($i, $start, $stop);
	    #print "$$seq_ref\n";
	    my @seq_split = split (//, $$seq_ref);
	    $seq_arr[$i] = \@seq_split;
	}
    }
    elsif ($strand eq "-") {
	for (my $i=0; $i<$n; $i++) {
	    my $seq_ref = $mr->getRevCompSeq($i, $start, $stop);
	    my @seq_split = split (//, $$seq_ref);
	    $seq_arr[$i] = \@seq_split;
	}
    }

    for (my $i=0; $i<@{$seq_arr[0]}; $i++) {
	my $type;
	if ($frame == 0) {
	    $type = $CDS0;
	}
	elsif ($frame == 1) {
	    $type = $CDS1;
	}
	elsif ($frame == 2) {
	    $type = $CDS2;
	}

	for (my $j=1; $j<$n; $j++) {
	    $counts[$j-1][$type][$TOTAL]++;
	    if ($seq_arr[$j][$i] ne ".") {
		$counts[$j-1][$type][$ALIGNED]++;
	    }
	    if ($seq_arr[$j][$i] eq $seq_arr[0][$i]) {
		$counts[$j-1][$type][$MATCH]++;
	    }
	}
	$frame = ($frame + 1) % 3;
    }
}

sub by_start {    
    my @gene_a = @$a;
    my @gene_b = @$b;

    return $gene_a[0][$START] <=> $gene_b[0][$START];
}
