#!/usr/bin/perl -w
use strict;
use GTF_Parser;

my $usage = "label_transcripts_by_overlap.pl <GTF file>\n";

@ARGV == 1 || die $usage;

my ($gtf_file) = @ARGV;

my %genes = %{GTF_Parser::parse_gtf($gtf_file)};
my %gtf_hash = %{GTF_Parser::get_gtf_hash($gtf_file)};

my %transcripts;
foreach my $gene (keys %genes) {
    foreach my $transcript (keys %{$genes{$gene}}) {
	$transcripts{$transcript} = $genes{$gene}{$transcript};
    }
}

my @transcript_ids = keys %transcripts;

geneLoop:
for (my $i=0; $i<@transcript_ids; $i++) { #transcript A
    for (my $j=0; $j<$i; $j++) { #transcript B (comes before A)
	my @a_exons = @{$transcripts{$transcript_ids[$i]}};
	my $a_start = $a_exons[0][$START];
	my $a_stop = $a_exons[@a_exons - 1][$STOP];
	my $a_strand = $a_exons[0][$STRAND];

	my @b_exons = @{$transcripts{$transcript_ids[$j]}};
	my $b_start = $b_exons[0][$START];
	my $b_stop = $b_exons[@b_exons - 1][$STOP];
	my $b_strand = $b_exons[0][$STRAND];
	
	if (($a_strand eq $b_strand &&
	     (($a_start <= $b_start && $a_stop >= $b_start) ||
	      ($a_start >= $b_start && $a_start <= $b_stop)))) {
	    #overlapping transcripts on the same strand
	    #print out gene A with gene ID set equal to B
	    print STDERR "$transcript_ids[$i] overlaps with $transcript_ids[$j]\n";

	    my $b_id = $transcript_ids[$j];
	    my $a_string = $gtf_hash{$transcript_ids[$i]};
	    $a_string =~ s/gene_id \"\S+\"/gene_id \"$b_id\"/g;
	    print $a_string;
	 
	    next geneLoop;
	}
    }
    #didn't find an overlap, so just print out gene A as is
    print $gtf_hash{$transcript_ids[$i]};
}
