#!/usr/bin/perl -w
use strict;
use GTF_Parser;
use FastaRandom;

my $MULTIPLE_TRANSCRIPTS = 1;  # should we include multiple transcripts per gene?
my $CANONICAL = 0;
my $MAX_LENGTH = -1;  # maximum exon length, -1 means no maximum
my $MIN_LENGTH = -1;  # minimum exon length, -1 means no minimum

#Goes through all GTF files in the given directory and outputs
#a corresponding "clean" GTF file in the given output directory
#of all transcripts meeting the following conditions:
#
# Start codon is ATG
# Stop codon is TAA, TAG, or TGA
# Acceptor sites are all AG, Donor sites all GT or GC (canonical mode)
# Acceptor sites are either AG or AC, 
# Donor sites either GT, GC, or AT (non-canonical mode)
# No in-frame stop codons
# Correct labeling of exon frames
# Total transcript length is a multiple of 3
# No zero-length introns

my $usage = "clean_gtf <GTF directory> <FASTA directory> <output directory>\n";

@ARGV == 3 || die $usage;

my ($gtf_dir, $seq_dir, $out_dir) = @ARGV;

opendir (GTF_DIR, $gtf_dir);
my @gtf_files = grep(/\.gtf$/, readdir(GTF_DIR));
closedir (GTF_DIR);

my $transcripts_ok = 0;
my $transcripts_with_errors = 0;
my $genes_ok = 0;
my $genes_with_errors = 0;
my $plus_errors = 0;
my $minus_errors = 0;
my $short_einits = 0;

for (my $i=0; $i<@gtf_files; $i++) {
  $gtf_files[$i] =~ /(chr[XYRLh\d]+)/;
  print STDERR "$gtf_files[$i]\n";
  my $seq_file = "$1.fa";
  my %genes = %{GTF_Parser::parse_gtf("$gtf_dir/$gtf_files[$i]")};
  my %gtf_hash = %{GTF_Parser::get_gtf_hash("$gtf_dir/$gtf_files[$i]")};
  my $fr = FastaRandom->new("$seq_dir/$seq_file");
  my @clean_transcripts;
  my %clean_starts;
  my %clean_stops;
  my %clean_strands;

  foreach my $gene (keys %genes) {
      my $best_transcript = "";
      my $max_coding_length = 0;
      my $max_exons = -1;
      my $gene_strand;

    transcriptLoop:
      foreach my $transcript (keys %{$genes{$gene}}) {  
	  my $exons = 0;
	  my $errors = 0;
	  my $coding_sequence = "";

	  my @exons = @{$genes{$gene}{$transcript}};
	  $gene_strand = $exons[0][$STRAND];

	  my $transcript_start = $exons[0][$START];
	  my $transcript_stop = $exons[@exons - 1][$STOP];

	  for (my $i=0; $i<@exons; $i++) {

	      my $start = $exons[$i][$START];
	      my $stop = $exons[$i][$STOP];
	      my $strand = $exons[$i][$STRAND];
	      my $type = $exons[$i][$TYPE];

	      my $length = $stop - $start + 1;
	      if ($start > $stop) { # check for negative length
		  $errors++;
		  print STDERR "Transcript $transcript has feature of length $length\n";
		  $transcripts_with_errors++;
		  next transcriptLoop;  #just give up on this transcript
	      }

	      # check against max and min lengths
	      if (($length > $MAX_LENGTH && $MAX_LENGTH > 0) ||
		  ($length < $MIN_LENGTH && $MIN_LENGTH > 0)) {
		  print STDERR "Exon length of $length out of bounds: $MIN_LENGTH to $MAX_LENGTH\n";
		  $errors++;
	      }

	      if ($type == $EINIT) {
		  if ($length < 3) {
		      $short_einits++;
		      print STDERR "********** Found $short_einits short initial exons\n";
		  }
		  $errors += check_donor($fr, $start, $stop, $strand);
		  add_to_transcript($fr, $start, $stop, $strand, \$coding_sequence);
		  $exons++;
	      }
	      elsif ($type == $EXON) {
		  $errors += check_donor($fr, $start, $stop, $strand);
		  $errors += check_acceptor($fr, $start, $stop, $strand);
		  add_to_transcript($fr, $start, $stop, $strand, \$coding_sequence);
		  $exons++;
	      }
	      elsif ($type == $ETERM) {
		  $errors += check_acceptor($fr, $start, $stop, $strand);
		  if ($strand eq "+") {
		      add_to_transcript($fr, $start, $stop+3, $strand, \$coding_sequence);
		  }
		  else {
		      add_to_transcript($fr, $start-3, $stop, $strand, \$coding_sequence);
		  }
		  $exons++;
	      }
	      elsif ($type == $ESNGL) {
		  if ($strand eq "+") {
		      add_to_transcript($fr, $start, $stop+3, $strand, \$coding_sequence);
		  }
		  else {
		      add_to_transcript($fr, $start-3, $stop, $strand, \$coding_sequence);
		  }
		  $exons++;
	      }
	      elsif ($type == $EA) {
		  $errors += check_acceptor($fr, $start, $stop, $strand);
	      }
	      elsif ($type == $ENC) {
		  $errors += check_donor($fr, $start, $stop, $strand);
		  $errors += check_acceptor($fr, $start, $stop, $strand);
	      }
	      elsif ($type == $EP) {
		  $errors += check_donor($fr, $start, $stop, $strand);
	      }
	      elsif ($type == $EPA) {
		  #check nothing
	      }
	      elsif ($type == $TS) {
		  $errors += check_donor($fr, $start, $stop, $strand);
	      }
	      elsif ($type == $TNC) {
		  $errors += check_donor($fr, $start, $stop, $strand);
		  $errors += check_acceptor($fr, $start, $stop, $strand);
	      }
	      elsif ($type == $TT) {
		  $errors += check_acceptor($fr, $start, $stop, $strand);
	      }
	      elsif ($type == $TST) {
		  #check nothing
	      }
	      else {
		  die "Unknown exon type\n";
	      }
	  }

	  my $coding_length = length($coding_sequence);
	  if ($coding_length % 3 != 0) {
	      $errors++;
	      print STDERR "*****\tTranscript $transcript ($gene_strand) $transcript_start - $transcript_stop coding length $coding_length not a multiple of 3\n";
	      $transcripts_with_errors++;  
	      next transcriptLoop;  #give up
	  }
	  
	  #make sure transcript has valid start and stop codons
	  my $start_codon = substr(uc($coding_sequence), 0, 3);
	  my $stop_codon = substr(uc($coding_sequence), length($coding_sequence) - 3, 3);
	  if ($start_codon ne "ATG") {
	      print STDERR "Invalid start codon: $start_codon\n";
	      $errors++;
	  } 
	  if ($stop_codon ne "TAA" && $stop_codon ne "TAG" && $stop_codon ne "TGA") {
	      print STDERR "Invalid stop codon: $stop_codon\n";
	      $errors++;
	  }

	  #scan transcript for in-frame stop codons
	  my @coding_bases = split (//, uc($coding_sequence));
	  for (my $i=0; $i<@coding_bases - 3; $i += 3) {
	      my $codon = $coding_bases[$i].
		  $coding_bases[$i+1].$coding_bases[$i+2];
	      if ($codon eq "TAA" || $codon eq "TAG" || $codon eq "TGA") {
		  $errors++;
		  my $strand = $exons[0][$STRAND];
		  print STDERR "Transcript $transcript ($strand) has in-frame stop codon at position $i: $codon\n";

		  my $cumulative_length = 0;
		  for (my $j=0; $j<@exons; $j++) {
		      my $type = $exons[$j][$TYPE];
		      if ($type == $EINIT || $type == $EXON || $type == $ETERM || $type == $ESNGL) { 
			  $cumulative_length += $exons[$j][$STOP] - $exons[$j][$START] + 1;
			  print STDERR "$cumulative_length\n";
		      }
		  }
	      }
	  }

	  #check to make sure this gene does not overlap any transcripts
	  #we have already decided are clean
	  if (! $MULTIPLE_TRANSCRIPTS) {
	      for (my $i=0; $i<@clean_transcripts; $i++) {
		  my $clean_start = $clean_starts{$clean_transcripts[$i]};
		  my $clean_stop = $clean_stops{$clean_transcripts[$i]};
		  my $this_start = $exons[0][$START];
		  my $this_stop = $exons[$#exons][$STOP];
		  
		  if (($this_start <= $clean_start && $this_stop >= $clean_start) ||
		      ($this_start >= $clean_start && $this_start <= $clean_stop)) {
		      $errors++;
		      print STDERR "Transcript $transcript overlaps with transcript $clean_transcripts[$i]\n";
		      print STDERR "$this_start - $this_stop vs. $clean_start - $clean_stop\n";
		      print STDERR "$exons[0][$STRAND]\t$clean_strands{$clean_transcripts[$i]}\n";
		  }
	      }
	  }

	  if ($exons == 0) {
	      $errors++;
	      print STDERR "Transcript $transcript had zero coding exons\n";
	  } 

	  if ($errors == 0) {  #this transcript is OK
	      #if ($coding_length > $max_coding_length) {
		  #$max_coding_length = $coding_length;
		  #$best_transcript = $transcript;
	      #}
	      if ($exons > $max_exons) {
		  $max_exons = $exons;
		  $best_transcript = $transcript;
	      }
	      if ($MULTIPLE_TRANSCRIPTS) {
		  push(@clean_transcripts, ($transcript));
	      }
	      $clean_starts{$transcript} = $exons[0][$START];
	      $clean_stops{$transcript} = $exons[$#exons][$STOP];
	      $clean_strands{$transcript} = $exons[0][$STRAND];
	  }
	  else {
	      print STDERR "*****\tError on transcript $transcript ($gene_strand) $transcript_start - $transcript_stop\n";
	  }
      }
      if ($best_transcript ne "") {
	  $genes_ok++;
	  if (! $MULTIPLE_TRANSCRIPTS) {
	      push (@clean_transcripts, ($best_transcript));
	  }
      }
      else {
	  $genes_with_errors++;
	  if ($gene_strand eq "+") {
	      $plus_errors++;
	  }
	  else {
	      $minus_errors++;
	  }
	  print STDERR "$genes_ok genes OK, $genes_with_errors genes with errors\n";
	  print STDERR (100 * $genes_ok / ($genes_ok + $genes_with_errors))."% OK\n";
	  print STDERR "$plus_errors errors on + strand, $minus_errors errors on minus strand\n";
	  print STDERR "\n";
      }
  }

  #print out all the clean transcripts 
  my @sorted_clean_transcripts = sort {$clean_starts{$a} <=> $clean_starts{$b}} @clean_transcripts;

  open (my $out_fh, ">$out_dir/$gtf_files[$i]");
  for (my $i=0; $i<@sorted_clean_transcripts; $i++) {
      print $out_fh $gtf_hash{$sorted_clean_transcripts[$i]};
  }
  close ($out_fh);
}


sub add_to_transcript {
  my ($fr, $start, $stop, $strand, $transcript_ref) = @_;

  if ($strand eq "+") {
    my $exon_seq_ref = $fr->getSeq($start, $stop);
    $$transcript_ref = ($$transcript_ref).($$exon_seq_ref);
  }

  elsif ($strand eq "-") {
    my $exon_seq_ref = $fr->getRevCompSeq($start, $stop);
    $$transcript_ref = ($$exon_seq_ref).($$transcript_ref);
  }
}

sub check_donor {
  my ($fr, $start, $stop, $strand) = @_;
  
  my $donor_ref;
  if ($strand eq "+") {
    $donor_ref = $fr->getSeq($stop+1, $stop+2);
  }
  elsif ($strand eq "-") {
    $donor_ref = $fr->getRevCompSeq($start-2, $start-1);
  }
  if (! defined $donor_ref) {
      print STDERR "undefined donor site sequence\n";
      return 1;
  }
  if (($CANONICAL && uc($$donor_ref) ne "GT" && uc($$donor_ref) ne "GC") ||  
      (uc($$donor_ref) ne "GT" && uc($$donor_ref) ne "GC" && uc($$donor_ref) ne "AT")) {
    print STDERR "check_donor failed, site was $$donor_ref, strand $strand\n";
    return 1;
  }
  else {
    return 0;
  }
}

sub check_acceptor {
  my ($fr, $start, $stop, $strand) = @_;

  my $acceptor_ref;
  if ($strand eq "+") {
    $acceptor_ref = $fr->getSeq($start-2, $start-1);
  }
  elsif ($strand eq "-") {
    $acceptor_ref = $fr->getRevCompSeq($stop+1, $stop+2);
  }
  if (! defined $acceptor_ref) {
      print STDERR "undefined acceptor site sequence\n";
      return 1;
  }
  if ( ($CANONICAL && uc($$acceptor_ref) ne "AG") ||
       (uc($$acceptor_ref) ne "AG" &&
	uc($$acceptor_ref) ne "AC")) {
    print STDERR "check_acceptor failed, site was $$acceptor_ref, strand $strand\n";
    return 1;
  }
  else {
    return 0;
  }
}
