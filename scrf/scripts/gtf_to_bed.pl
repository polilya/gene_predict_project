#!/usr/bin/perl -w
use strict;
use GTF_Parser;

my $OFFSET = 0;

my @SCORE_BINS = (0.9, 0.99, 1000);
my @COLORS = (13448960,16753920,2263842);


my $usage = "$0 <GTF directory>\n";

@ARGV == 1 || die $usage;

my ($gtf_dir) = @ARGV; 

print "track name=\"CONTRAST\" description=\"CONTRAST Gene Predictions\" type=\"coloredExon\"\n";

opendir(GTF_DIR, $gtf_dir);
my @gtf_files = grep(/\.gtf$/, readdir(GTF_DIR));
closedir(GTF_DIR);

for (my $filenum=0; $filenum<@gtf_files; $filenum++) {

    $gtf_files[$filenum] =~ /(.+)\.gtf/;
    my $sequence_name = $1;

    my %genes = %{GTF_Parser::parse_gtf("$gtf_dir/$gtf_files[$filenum]")};

    my $id = 1;

    my @gene_array = sort {
	my @a_transcripts = keys %{$genes{$a}};
	my @b_transcripts = keys %{$genes{$b}};
	$genes{$a}{$a_transcripts[0]}[0][$START] <=> $genes{$b}{$b_transcripts[0]}[0][$START];
    } (keys %genes); 
    
    for (my $i=0; $i<@gene_array; $i++) {
	foreach my $transcript (keys %{$genes{$gene_array[$i]}}) {
	    my @exons = @{$genes{$gene_array[$i]}{$transcript}};
	    my $strand = $exons[0][$STRAND];
	    
	    my $gene_start = $exons[0][$START] - 1 + $OFFSET;
	    my $gene_stop = $exons[@exons - 1][$STOP] + $OFFSET;  #weird UCSC half-open coordinates
	    if ($exons[0][$STRAND] eq "+") {
		$gene_stop += 3;  # stop codon included in UCSC gene displays
	    }
	    else {
		$gene_start -= 3; # stop codon included in UCSC gene displays
	    }
	    
	    my $num_exons = @exons + 0;
	    my @exon_starts;
	    my @exon_lengths;
	    my @exon_scores;
	    for (my $i=0; $i<@exons; $i++) {
		my $exon_start = $exons[$i][$START] -1 + $OFFSET - $gene_start;
		if ($exons[$i][$STRAND] eq "-" && ($exons[$i][$TYPE] == $ETERM || $exons[$i][$TYPE] == $ESNGL)) {
		    $exon_start -= 3;  # stop codon included in UCSC gene displays
		}
		push(@exon_starts, $exon_start);
		
		my $exon_length = $exons[$i][$STOP] - $exons[$i][$START] + 1;
		if ($exons[$i][$TYPE] == $ETERM || $exons[$i][$TYPE] == $ESNGL) {
		    $exon_length += 3;  # stop codon included in UCSC gene displays
		}
		push(@exon_lengths, $exon_length);
		push(@exon_scores, $exons[$i][$SCORE]);
	    }
	    print "$sequence_name\t$gene_start\t$gene_stop\tCONTRAST.$sequence_name.$id\t1000\t$strand\t$gene_start\t$gene_stop\t0\t$num_exons\t";
	    for (my $i=0; $i<@exon_lengths; $i++) {
		print "$exon_lengths[$i],";
	    }
	    print "\t";
	    for (my $i=0; $i<@exon_starts; $i++) {
		print "$exon_starts[$i],";
	    }
	    print "\t$num_exons\t";
	    for (my $i=0; $i<@exon_scores; $i++) {
		my $bin = 0;
		while ($exon_scores[$i] > $SCORE_BINS[$bin]) {
		    $bin++;
		}
		print "$COLORS[$bin],";
	    }
	    print "\n";
	    
	    $id++;
	}
    }
}
