use strict;
package GTF_Parser_UTR_SU;

#Provides functions for parsing a GTF file
#into Zoe states according to the UTR model

use base "Exporter";
our @EXPORT = qw($EINITU $EINITS $EXON $ETERM $ESNGL $EA $ENC $EP $EPA $INTRON $INC $INTER @type_names $TYPE $START $STOP $STRAND $ISOCHORE $FRAME $ID $SCORE $NO_FRAME $UNDEFINED_ISOCHORE);

our $EINITU = 0;
our $EINITS = 1;
our $EXON = 2;
our $ETERM = 3;
our $ESNGL = 4;
our $EA = 5;
our $ENC = 6;
our $EP = 7;
our $EPA = 8;
our $INTRON = 9;
our $INC = 10;
our $INTER = 11;

our @type_names = qw (EinitU EinitS Exon Eterm Esngl Ea Enc Ep Epa Intron Inc);

our $TYPE = 0;
our $START = 1;
our $STOP = 2;
our $STRAND = 3;
our $ISOCHORE = 4;
our $FRAME = 5;
our $ID = 6;
our $SCORE = 7;

our $NO_FRAME = 3;
our $UNDEFINED_ISOCHORE = -2;

#returns a hash which maps gene IDs to strings composed
#of the portion of the GTF file describing that gene
sub get_gtf_hash {
  my ($gtf_file) = @_;

  my %gtf_hash;

  open (GTF_FILE, $gtf_file) || die "Couldn't open $gtf_file";
  my @gtf_lines = <GTF_FILE>;
  close (GTF_FILE);

  my $current_gene_id = "";
  my $gene_lines;
  for (my $i=0; $i<@gtf_lines; $i++) {
    $gtf_lines[$i] =~ /\S+\t\S+\t(\S+)\t(\S+)\t(\S+)\t\S+\t(\S+)\t(\S+)\tgene_id \"(\S+)\"/;  
    my $id = $6;
    if (exists $gtf_hash{$id}) {
      $gtf_hash{$id} .= $gtf_lines[$i]; 
    }
    else {
      $gtf_hash{$id} = $gtf_lines[$i];
    }
  }

  return \%gtf_hash;
}


sub parse_gtf {

  my ($gtf_file, $iso_file) = @_;

  open (GTF_FILE, $gtf_file);
  my @gtf_lines = <GTF_FILE>;
  close (GTF_FILE);

  my @isochore_lines;
  if (defined $iso_file) {
    open (ISO_FILE, $iso_file);
    @isochore_lines = <ISO_FILE>;
    close (ISO_FILE);
  }
  
  my @genes;
  my $current_gene_id = "";
  my $gene_num = -1;
  my $feature_num;
  
  my $isochore; #isochore current gene falls in
  my $num_cds; #number of CDS features in current gene
  my $num_utr; #number of 5UTR features in current gene
  my $current_cds; #number of current CDS we are processing
  my $current_utr; #number of current UTR we are processing
  for (my $j=0; $j<@gtf_lines; $j++) {
    $gtf_lines[$j] =~ /\S+\t\S+\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\tgene_id \"(\S+)\"/;
    my $feature = $1;
    my $start = $2;
    my $stop = $3;
    my $score = $4;
    my $strand = $5;
    my $frame = $6;
    my $id = $7;
    if ($id ne $current_gene_id) {
      #a new gene
      $gene_num++;
      $current_gene_id = $id;
      $feature_num = 0;
      $current_cds = 0;
      $current_utr = 0;
      
      my $gene_start = $start;
      my $gene_stop = $stop;
      
      #count the number of CDS and 5UTR features
      $num_cds = 0;
      $num_utr = 0;
      my $k=$j;
      while ($k<@gtf_lines && $gtf_lines[$k] =~ /$id/) {
	$gtf_lines[$k] =~ /\S+\t\S+\t(\S+)\t\S+\t(\S+)/;
	$gene_stop = $2;
	if ($1 eq "CDS") {
	  $num_cds++;
	}
	elsif ($1 eq "5UTR") {
	  $num_utr++;
	}
	$k++;
      }
      
      #determine isochore of this gene
      if (defined $iso_file) { 
	$isochore = lookup_isochore(\@isochore_lines, $gene_start, $gene_stop);
      }
      else {
	$isochore = $UNDEFINED_ISOCHORE;
      }
    }
    
    my $type;
    if ($feature eq "5UTR") {
      $current_utr++;
	
      #determine UTR type
      if ($num_utr == 1) {
	$type = $EPA;
      }
      elsif ($strand eq "+") {
	if ($current_utr == 1) {
	  $type = $EP;
	}
	elsif ($current_utr == $num_utr) {
	  $type = $EA;
	}
	else {
	  $type = $ENC;
	}
      }
      elsif ($strand eq "-") {
	if ($current_utr == 1) {
	  $type = $EA;
	}
	elsif ($current_utr == $num_utr) {
	  $type = $EP;
	}
	else {
	  $type = $ENC;
	}
      }
      $genes[$gene_num][$feature_num][$TYPE] = $type;
      $genes[$gene_num][$feature_num][$START] = $start;
      $genes[$gene_num][$feature_num][$STOP] = $stop;
      $genes[$gene_num][$feature_num][$STRAND] = $strand;
      $genes[$gene_num][$feature_num][$ISOCHORE] = $isochore;
      $genes[$gene_num][$feature_num][$FRAME] = $NO_FRAME;
      $genes[$gene_num][$feature_num][$SCORE] = $score;
      
      $feature_num++;

      #add Inc state, if appropriate
      if (($strand eq "+" && ($type == $EP || $type == $ENC)) ||
	  ($strand eq "-" && ($type == $EA || $type == $ENC))) {
	$gtf_lines[$j+1] =~ /\S+\t\S+\t\S+\t(\S+)/;
	my $next_start = $1;
	$genes[$gene_num][$feature_num][$TYPE] = $INC;
	$genes[$gene_num][$feature_num][$START] = $stop+1;
	$genes[$gene_num][$feature_num][$STOP] = $next_start-1;
	$genes[$gene_num][$feature_num][$STRAND] = $strand;
	$genes[$gene_num][$feature_num][$ISOCHORE] = $isochore;
	$genes[$gene_num][$feature_num][$FRAME] = $NO_FRAME;
	$genes[$gene_num][$feature_num][$SCORE] = ".";
	$feature_num++;
      }
    }
    elsif ($feature eq "CDS") {
      $current_cds++;
      
      #determine CDS type
      if ($num_cds == 1) {
	$type = $ESNGL;
      }
      elsif ($strand eq "+") {
	if ($current_cds == 1) {
	  if ($num_utr == 1) {
	    $type = $EINITU;
	  }
	  else {
	    $type = $EINITS;
	  }
	}
	elsif ($current_cds == $num_cds) {
	  $type = $ETERM;
	}
	  else {
	    $type = $EXON;
	  }
      }
      elsif ($strand eq "-") {
	if ($current_cds == 1) {
	  $type = $ETERM;
	  }
	elsif ($current_cds == $num_cds) {
	  if ($num_utr == 1) {
	    $type = $EINITU;
	  }
	  else {
	    $type = $EINITS;
	  }
	}
	else {
	  $type = $EXON;
	}
      }

      $genes[$gene_num][$feature_num][$TYPE] = $type;
      $genes[$gene_num][$feature_num][$START] = $start;
      $genes[$gene_num][$feature_num][$STOP] = $stop;      
      $genes[$gene_num][$feature_num][$STRAND] = $strand;
      $genes[$gene_num][$feature_num][$ISOCHORE] = $isochore;
      $genes[$gene_num][$feature_num][$FRAME] = $frame;
      $genes[$gene_num][$feature_num][$SCORE] = $score;
      $feature_num++;

      #add Intron state, if appropriate
      if (($strand eq "+" && ($type == $EINITS || $type == $EINITU ||
			      $type == $EXON)) ||
	  ($strand eq "-" && ($type == $ETERM || $type == $EXON))) {
	$gtf_lines[$j+1] =~ /\S+\t\S+\t\S+\t(\S+)/;
	my $next_start = $1;
	$genes[$gene_num][$feature_num][$TYPE] = $INTRON;
	$genes[$gene_num][$feature_num][$START] = $stop+1;
	$genes[$gene_num][$feature_num][$STOP] = $next_start-1;
	$genes[$gene_num][$feature_num][$STRAND] = $strand;
	$genes[$gene_num][$feature_num][$ISOCHORE] = $isochore;
	$genes[$gene_num][$feature_num][$FRAME] = $NO_FRAME;
	$genes[$gene_num][$feature_num][$SCORE] = ".";
	$feature_num++;
      }
    }
  }

  #set all exon frames
  for (my $i=0; $i<@genes; $i++) {
      set_frames ($genes[$i]);
  }

  return \@genes;
}

sub lookup_isochore {
  my ($isochore_lines_ref, $gene_start, $gene_stop) = @_;
  for (my $i=0; $i<@{$isochore_lines_ref}; $i++) {
    $isochore_lines_ref->[$i] =~ /(\S+)\t(\S+)\t\S+\t(\S+)/;
    my $frag_start = $1;
    my $frag_stop = $2;
    my $isochore;
    if ($3 eq "X") {
      $isochore = $UNDEFINED_ISOCHORE;
    }
    else {
      $isochore = $3 - 1;  #we index isochores 0, 1, 2, 3 not I, II, III, IV
    }
    if ($gene_start >= $frag_start && $gene_stop <= $frag_stop) {
      return $isochore;
    }
  }
  return $UNDEFINED_ISOCHORE;  #gene fell on the border between two 1Mb fragments
}

sub set_frames {
    my ($gene_ref) = @_;
    my @gene = @{$gene_ref};

    if ($gene[0][$STRAND] eq "+") {
	my $frame = 0;
	for (my $i=0; $i<@gene; $i++) {
	    my $length = $gene[$i][$STOP] - $gene[$i][$START] + 1;

	    if ($gene[$i][$TYPE] == $EINITS ||
		$gene[$i][$TYPE] == $EINITU ||
		$gene[$i][$TYPE] == $ESNGL) {
		$gene[$i][$FRAME] = 0;
		$frame = ($frame + ($length % 3)) % 3;
	    }
	    elsif ($gene[$i][$TYPE] == $EXON ||
		   $gene[$i][$TYPE] == $ETERM) {
		$gene[$i][$FRAME] = $frame;
		$frame = ($frame + ($length % 3)) % 3;
	    }
	    else {
		$gene[$i][$FRAME] = $NO_FRAME;
	    }
	}
    }
    elsif ($gene[0][$STRAND] eq "-") {
	my $frame = 0;
	for (my $i=@gene-1; $i>=0; $i--) {
	    my $length = $gene[$i][$STOP] - $gene[$i][$START] + 1;

	    if ($gene[$i][$TYPE] == $EINITS ||
		$gene[$i][$TYPE] == $EINITU ||
		$gene[$i][$TYPE] == $ESNGL) {
		$gene[$i][$FRAME] = 0;
		$frame = ($frame + ($length % 3)) % 3;
	    }
	    elsif ($gene[$i][$TYPE] == $EXON ||
		   $gene[$i][$TYPE] == $ETERM) {
		$gene[$i][$FRAME] = $frame;
		$frame = ($frame + ($length % 3)) % 3;
	    }
	    else {
		$gene[$i][$FRAME] = $NO_FRAME;
	    }
	}
    }
}
