use strict;
package FastaRandom;

#Provides fast random access into a single-sequence FASTA file 
#without having to load the entire sequence into memory

#***** This class assumes that every non-header line in the 
#      FASTA file is of the same length.  If this assumption does not
#      hold, the class will not work properly!

#***** This class does not understand multi-sequence FASTA files!

my %complement = (A => "T",
		  C => "G",
		  G => "C",
		  T => "A",
		  a => "t",
		  c => "g",
		  g => "c",
		  t => "a",
		  N => "N");

#Creates a new FastaRandom instance allowing random
#access into the given FASTA file
sub new {
  shift;
  my ($filename) = @_;
  open (my $fh, $filename) || die "FastaRandom constructor:  could not open Fasta file $filename\n";
  my $header_length = length(<$fh>);
  my $line_length = length(<$fh>) - 1;
  my $this = bless({FileHandle => $fh, HeaderLength => $header_length,
	      LineLength => $line_length});
  return $this;
}

#closes the filehandle
#should be invoked when we're done with the FastaRandom object
sub close {
  my $this = shift;
  close ($this->{FileHandle});
}

#returns a reference to the sequence defined by the given start and stop 
#coordinates
#note that the first base in the sequence is considered 1 (not 0), like in GTF
sub getSeq {
  my $this = shift;
  my ($start, $stop) = @_;
  
  if ($start > $stop) {
      my $seq = "";
      return \$seq;
  }

  #read in the sequence from disk
  my $startOffset = $this->{HeaderLength} + $start-1 + 
      int(($start-1)/$this->{LineLength});
  my $stopOffset = $this->{HeaderLength} + $stop-1 + 
      int(($stop-1)/$this->{LineLength});
  my $readLength = $stopOffset - $startOffset + 1;

  my $raw_sequence;
  sysseek ($this->{FileHandle}, $startOffset, 0);
  sysread ($this->{FileHandle}, $raw_sequence, $readLength);
  
  #remove the '\n' characters from the sequence
  my $seq = "";
  my @raw_arr = split (/\n/, $raw_sequence);
  for (my $i=0;$i<@raw_arr; $i++) {
    $seq .= $raw_arr[$i];
  }
  return \$seq;
}

sub getRevCompSeq {
  my $this = shift;
  my ($start, $stop) = @_;
  
  #first get the original sequence
  my $seq_ref = getSeq($this, $start, $stop);
  my @seq_arr = split(//, $$seq_ref);

  #now reverse complement it
  my $rev_comp_seq = "";
  my $j=0;
  for (my $i=$#seq_arr; $i>=0; $i--) {
    $rev_comp_seq .= $complement{$seq_arr[$i]};
  }
  return \$rev_comp_seq;
}

sub getRevSeq {
  my $this = shift;
  my ($start, $stop) = @_;
  
  #first get the original sequence
  my $seq_ref = getSeq($this, $start, $stop);
  my @seq_arr = split(//, $$seq_ref);

  #now reverse it
  my $rev_seq = "";
  my $j=0;
  for (my $i=$#seq_arr; $i>=0; $i--) {
    $rev_seq .= $seq_arr[$i];
  }
  return \$rev_seq;
}

return 1;  #Obligatory (Perl is weird)
