use strict;
package MultiRandom;

#Provides fast random access into a multiple sequence alignment file 
#without having to load the entire sequence into memory
#The format expected is the PHYLOSCAN input format

my %complement = ("A" => "T",
		  "C" => "G",
		  "G" => "C",
		  "T" => "A",
		  "a" => "t",
		  "c" => "g",
		  "g" => "c",
		  "t" => "a",
		  "_" => "_",
		  "." => ".",
		  "N" => "N");

sub new {
  shift;
  my ($filename) = @_;
  my $file_size = (-s $filename);
  open (my $fh, $filename) || die "MultiRandom constructor:  could not open file $filename\n";
  my $header = <$fh>;
  my $header_length = length($header);
  chomp $header;
  $header =~ s/>//;
  my @seq_names = split (/ /, $header);
  my $num_seq = @seq_names + 0;
  my $line_length = ($file_size - $header_length)/$num_seq;

  my $this = bless({FileHandle => $fh, HeaderLength => $header_length,
		    LineLength => $line_length, SeqNames => \@seq_names, 
		    N => $num_seq});
  return $this;
}

#closes the filehandle
#should be invoked when we're done with the FastaRandom object
sub close {
  my $this = shift;
  close ($this->{FileHandle});
}

#returns a reference to the sequence defined by the given start and stop 
#coordinates and sequence number
#note that the first base in the sequence is considered 1 (not 0), like in GTF
sub getSeq {
  my $this = shift;
  my ($seq_num, $start, $stop) = @_;
  
  #read in the sequence from disk
  my $startOffset = $this->{HeaderLength} + $seq_num * $this->{LineLength} + 
      $start - 1;
  my $stopOffset = $this->{HeaderLength} + $seq_num * $this->{LineLength} + 
      $stop - 1;
  my $readLength = $stopOffset - $startOffset + 1;

  my $sequence;
  sysseek ($this->{FileHandle}, $startOffset, 0);
  sysread ($this->{FileHandle}, $sequence, $readLength);
  
  return \$sequence;
}

sub getRevCompSeq {
    my $this = shift;
    my ($seq_num, $start, $stop) = @_;

    my $orig_seq_ref = getSeq($this, $seq_num, $start, $stop);
    my @orig_seq_arr = split(//, $$orig_seq_ref);

    my $rev_comp_seq = "";
    for (my $i=@orig_seq_arr - 1; $i>=0; $i--) {
	$rev_comp_seq .= $complement{$orig_seq_arr[$i]};
    }

    return \$rev_comp_seq;
}

return 1;
