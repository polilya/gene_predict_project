#!/usr/bin/perl -w
use strict;

my $usage = "split_gtf.pl <GTF directory> <output directory> <fragment size>\n"; 

@ARGV == 3 || die $usage;

my ($gtf_dir, $out_dir, $frag_size) = @ARGV;

opendir(GTF_DIR, $gtf_dir);
my @gtf_dir_contents = readdir(GTF_DIR);
closedir(GTF_DIR);

for (my $i=0; $i<@gtf_dir_contents; $i++) {
    if ($gtf_dir_contents[$i] =~ /chr(\S+)\.gtf/) {
	my $chr_name = $1;
	mkdir "$out_dir/$chr_name";

	my $gtf_file = "$gtf_dir/chr$chr_name.gtf";
	my %split_files;

	open (GTF_IN, $gtf_file);
	while (<GTF_IN>) {
	    $_ =~ /\S+\t\S+\t\S+\t(\S+)\t(\S+)/;
	    my $start = $1;
	    my $end = $2;
	    if (int($start / $frag_size) == int($end / $frag_size)) {
		#feature not split across fragments, go ahead
		my $frag_num = int($start / $frag_size);
		my $out_string = $_;
		
		#replace start/end coordinates with ones that make this
		#fragment begin at zero
		my $new_start = $start - $frag_num * $frag_size;
		my $new_end = $end - $frag_num * $frag_size;
		
		$out_string =~ s/$start/$new_start/;
		$out_string =~ s/$end/$new_end/;
		
		if (exists $split_files{$frag_num + 1}) {
		    $split_files{$frag_num + 1} .= $out_string;
		}
		else {
		    $split_files{$frag_num + 1} = $out_string;
		}
	    } 
	}
	close(GTF_IN);
	
	foreach my $frag_num (keys %split_files) {
	    mkdir "$out_dir/$chr_name/fragment_$frag_num";
	    
	    my $out_file = "$out_dir/$chr_name/fragment_$frag_num/$frag_num.gtf";
	    print STDERR "Writing $out_file\n";
	    
	    open (FRAG_OUT, ">$out_file");
	    print FRAG_OUT $split_files{$frag_num}; 
	    close (FRAG_OUT);
	}
	
	# make sure all fragments with no annotations have an empty GTF file
	opendir(SPLIT_DIR, "$out_dir/$chr_name");
	my @fragment_dirs = grep(/fragment/, readdir(SPLIT_DIR));
	for (my $j=0; $j<@fragment_dirs; $j++) {
	    if ($fragment_dirs[$j] =~ /fragment_(\d+)/) {
		my $frag_num = $1;
		my $command = "touch $out_dir/$chr_name/$fragment_dirs[$j]/$frag_num.gtf";
		print STDERR "Running: $command\n";
		system($command);
	    }
	}
	closedir(SPLIT_DIR);
    }
}
