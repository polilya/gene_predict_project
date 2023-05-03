use strict;
package GTF_Parser;

#Provides functions for parsing a GTF file

use base "Exporter";
our @EXPORT = qw($EINIT $EXON $ETERM $ESNGL $EA $ENC $EP $EPA $TS $TNC $TT $TST @type_names $TYPE $START $STOP $SCORE $STRAND $FOH $GENE_ID $TRANSCRIPT_ID);

our $EINIT = 0;
our $EXON = 1;
our $ETERM = 2;
our $ESNGL = 3;
our $EA = 4;
our $ENC = 5;
our $EP = 6;
our $EPA = 7;
our $TS = 8;
our $TNC = 9;
our $TT = 10;
our $TST = 11;

our @type_names = qw (Einit Exon Eterm Esngl Ea Enc Ep Epa Ta Tnc Tp Tpa);

our $TYPE = 0;
our $START = 1;
our $STOP = 2;
our $STRAND = 3;
our $FOH = 4;
our $GENE_ID = 5;
our $TRANSCRIPT_ID = 6;
our $SCORE = 7;

our $NO_FOH = 3;

sub get_gtf_hash {
  my ($gtf_file) = @_;

  my %gtf_hash;

  open (GTF_FILE, $gtf_file) || die "Couldn't open $gtf_file";
  my @gtf_lines = <GTF_FILE>;
  close (GTF_FILE);

  my $current_transcript_id = "";
  my $gene_lines;
  for (my $i=0; $i<@gtf_lines; $i++) {
    $gtf_lines[$i] =~ /\S+\t\S+\t(\S+)\t(\S+)\t(\S+)\t\S+\t(\S+)\t(\S+)\tgene_id \"(\S+)\"; transcript_id \"(\S+)\"/;  
    my $gene_id = $6;
    my $transcript_id = $7;
    if (exists $gtf_hash{$transcript_id}) {
      $gtf_hash{$transcript_id} .= $gtf_lines[$i]; 
    }
    else {
      $gtf_hash{$transcript_id} = $gtf_lines[$i];
    }
  }

  return \%gtf_hash;
}

sub parse_gtf {

    my ($gtf_file) = @_;
    
    open (GTF_FILE, $gtf_file) || die "Couldn't open $gtf_file";
    my @gtf_lines = <GTF_FILE>;
    close (GTF_FILE);
    
    my %genes;
    
    for (my $j=0; $j<@gtf_lines; $j++) {
	$gtf_lines[$j] =~ /\S+\t\S+\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\tgene_id \"(\S+)\"; transcript_id \"(\S+)\"/;
	my $feature = $1;
	my $start = $2;
	my $stop = $3;
	my $score = $4;
	my $strand = $5;
	my $foh = $6;
	my $gene_id = $7;
	my $transcript_id = $8;
	
	if ($feature eq "CDS" || $feature eq "5UTR" || $feature eq "3UTR") {
	    my $feature_num = 0;
	    if (exists $genes{$gene_id}{$transcript_id}) {
		$feature_num = @{$genes{$gene_id}{$transcript_id}} + 0;
	    }
	    
	    $genes{$gene_id}{$transcript_id}[$feature_num][$START] = $start;
	    $genes{$gene_id}{$transcript_id}[$feature_num][$STOP] = $stop;
	    $genes{$gene_id}{$transcript_id}[$feature_num][$SCORE] = $score;
	    $genes{$gene_id}{$transcript_id}[$feature_num][$STRAND] = $strand;
	    $genes{$gene_id}{$transcript_id}[$feature_num][$GENE_ID] = $gene_id;
	    $genes{$gene_id}{$transcript_id}[$feature_num][$TRANSCRIPT_ID] = $transcript_id;

	    if ($feature eq "CDS") {
		$genes{$gene_id}{$transcript_id}[$feature_num][$TYPE] = $EXON; #placeholder
	    }  
	    elsif ($feature eq "5UTR") {
		$genes{$gene_id}{$transcript_id}[$feature_num][$TYPE] = $EPA; #placeholder
	    }  
	    elsif ($feature eq "3UTR") {
		$genes{$gene_id}{$transcript_id}[$feature_num][$TYPE] = $TST; #placeholder
	    }  
	}
    }

    # sort all features in each transcript by start coordinate
    foreach my $gene (keys %genes) {
	foreach my $transcript (keys %{$genes{$gene}}) {
	    my @sorted_exons = sort by_start_coordinate @{$genes{$gene}{$transcript}};
	    $genes{$gene}{$transcript} = \@sorted_exons;
	}
    }

    # set feature types
    foreach my $gene (keys %genes) {
	foreach my $transcript (keys %{$genes{$gene}}) {

	    # first count number of 5' UTR, 3' UTR, and CDS features
	    my $num_5utr = 0;
	    my $num_3utr = 0;
	    my $num_cds = 0;
	    my @exons = @{$genes{$gene}{$transcript}};
	    for (my $i=0; $i<@exons; $i++) {
		if ($exons[$i][$TYPE] == $EXON) {
		    $num_cds++;
		}
		elsif ($exons[$i][$TYPE] == $EPA) {
		    $num_5utr++;
		}
		elsif ($exons[$i][$TYPE] == $TST) {
		    $num_3utr++;
		}

	    }

	    my $strand = $exons[0][$STRAND];
	    my $num_cds_seen = 0;
	    my $num_5utr_seen = 0;
	    my $num_3utr_seen = 0;

	    # now set types
	    for (my $i=0; $i<@exons; $i++) {
		if ($exons[$i][$TYPE] == $EXON) {
		    $num_cds_seen++;

		    if ($num_cds == 1) {
			$exons[$i][$TYPE] = $ESNGL;
		    }
		    elsif ( ($strand eq "+" && $num_cds_seen == 1) ||
			    ($strand eq "-" && $num_cds_seen == $num_cds) ) {
			$exons[$i][$TYPE] = $EINIT;
		    }
		    elsif ( ($strand eq "+" && $num_cds_seen == $num_cds) ||
			    ($strand eq "-" && $num_cds_seen == 1) ) {
			$exons[$i][$TYPE] = $ETERM;
		    }
		    else {
			$exons[$i][$TYPE] = $EXON;
		    }
		}
		elsif ($exons[$i][$TYPE] == $EPA) {
		    $num_5utr_seen++;

		    if ($num_5utr == 1) {
			$exons[$i][$TYPE] = $EPA;
		    }
		    elsif ( ($strand eq "+" && $num_5utr_seen == 1) ||
			    ($strand eq "-" && $num_5utr_seen == $num_5utr) ) {
			$exons[$i][$TYPE] = $EP;
		    }
		    elsif ( ($strand eq "+" && $num_5utr_seen == $num_5utr) ||
			    ($strand eq "-" && $num_5utr_seen == 1) ) {
			$exons[$i][$TYPE] = $EA;
		    }
		    else {
			$exons[$i][$TYPE] = $ENC;
		    }
		}
		elsif ($exons[$i][$TYPE] == $TST) {
		    $num_3utr_seen++;

		    if ($num_3utr == 1) {
			$exons[$i][$TYPE] = $TST;
		    }
		    elsif ( ($strand eq "+" && $num_3utr_seen == 1) ||
			    ($strand eq "-" && $num_3utr_seen == $num_3utr) ) {
			$exons[$i][$TYPE] = $TS;
		    }
		    elsif ( ($strand eq "+" && $num_3utr_seen == $num_3utr) ||
			    ($strand eq "-" && $num_3utr_seen == 1) ) {
			$exons[$i][$TYPE] = $TT;
		    }
		    else {
			$exons[$i][$TYPE] = $TNC;
		    }		    
		}
	    }
	}
    }

    #set all coding exon feature foh
    foreach my $gene (keys %genes) {
	#set_fohs ($genes[$i]);
    }

#   foreach my $gene (keys %genes) {
#	foreach my $transcript (keys %{$genes{$gene}}) {
#	    print STDERR "Transcript $transcript\n";
#	    my @exons = @{$genes{$gene}{$transcript}};
#	    for (my $i=0; $i<@exons; $i++) {
#		print STDERR "$type_names[$exons[$i][$TYPE]]\t$exons[$i][$START]\t$exons[$i][$STOP]\n";
#	    }
#	}
#      print STDERR "\n";
#    }
    
    return \%genes;
}

sub by_start_coordinate {
    return $a->[$START] <=> $b->[$START];
}

sub set_fohs {
    my ($gene_ref) = @_;
    my @gene = @{$gene_ref};

    if ($gene[0][$STRAND] eq "+") {
	my $foh = 0;
	for (my $i=0; $i<@gene; $i++) {
	    my $length = $gene[$i][$STOP] - $gene[$i][$START] + 1;

	    if ($gene[$i][$TYPE] == $EINIT ||
		$gene[$i][$TYPE] == $ESNGL) {
		$gene[$i][$FOH] = 0;
		$foh = ($foh + ($length % 3)) % 3;
	    }
	    elsif ($gene[$i][$TYPE] == $EXON ||
		   $gene[$i][$TYPE] == $ETERM) {
		$gene[$i][$FOH] = $foh;
		$foh = ($foh + ($length % 3)) % 3;
	    }
	    else {
		$gene[$i][$FOH] = $NO_FOH;
	    }
	}
    }
    elsif ($gene[0][$STRAND] eq "-") {
	my $foh = 0;
	for (my $i=@gene-1; $i>=0; $i--) {
	    my $length = $gene[$i][$STOP] - $gene[$i][$START] + 1;

	    if ($gene[$i][$TYPE] == $EINIT ||
		$gene[$i][$TYPE] == $ESNGL) {
		$gene[$i][$FOH] = 0;
		$foh = ($foh + ($length % 3)) % 3;
	    }
	    elsif ($gene[$i][$TYPE] == $EXON ||
		   $gene[$i][$TYPE] == $ETERM) {
		$gene[$i][$FOH] = $foh;
		$foh = ($foh + ($length % 3)) % 3;
	    }
	    else {
		$gene[$i][$FOH] = $NO_FOH;
	    }
	}
    }
}
