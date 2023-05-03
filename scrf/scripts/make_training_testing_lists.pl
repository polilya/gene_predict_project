#!/usr/bin/perl -w
use strict;

my $CRF_FRACTION = 0.5;
my $SVM_FRACTION = 0.25;
my $HOLDOUT_FRACTION = 0.25;

my $usage = "make_training_testing_lists.pl <chr split directory> <number of cross validation sets>";

@ARGV == 2 || die $usage;

my ($split_dir, $num_cv_sets) = @ARGV;

opendir (SPLIT_DIR, $split_dir);
my @chr_dirs = readdir(SPLIT_DIR);
closedir (SPLIT_DIR);

if ($num_cv_sets == 0) {
    # no cross validation, just divide everything into three sets
    open (CRF_LIST, ">$split_dir/crf_list.txt");
    open (SVM_LIST, ">$split_dir/svm_list.txt");
    open (HOLDOUT_LIST, ">$split_dir/holdout_list.txt");
    open (ANN_LIST, ">$split_dir/ann_list.txt");
    open (PRED_LIST, ">$split_dir/pred_list.txt");

    for (my $i=0; $i<@chr_dirs; $i++) {
	print STDERR "Processing $split_dir/$chr_dirs[$i]\n";
	opendir (CHR_DIR, "$split_dir/$chr_dirs[$i]");
	my @frag_dirs = grep(/fragment_\d+/, readdir(CHR_DIR));
	for (my $j=0; $j<@frag_dirs; $j++) {
	    print STDERR "\tProcessing $frag_dirs[$j]: ";
	    $frag_dirs[$j] =~ /fragment_(\d+)/;
	    my $frag_num = $1;
	    my $align_file = "$chr_dirs[$i]/$frag_dirs[$j]/$frag_num.align.gz";
	    my $estseq_file = "$chr_dirs[$i]/$frag_dirs[$j]/$frag_num.estseq.gz";
	    my $gtf_file = "$chr_dirs[$i]/$frag_dirs[$j]/$frag_num.gtf";
	    my $random_val = rand();
	    if ($random_val < $CRF_FRACTION) {
		# add to CRF training list
		print STDERR "CRF training\n";
		print CRF_LIST "$align_file\t$estseq_file\t$gtf_file\n";
	    }
	    elsif ($random_val < $CRF_FRACTION + $SVM_FRACTION) {
		# add to SVM training list
		print STDERR "SVM training\n";
		print SVM_LIST "$align_file\t$estseq_file\t$gtf_file\n";
	    }
	    else {
		# add to holdout list
		print STDERR "Holdout\n";
		print HOLDOUT_LIST "$align_file\t$estseq_file\t$gtf_file\n";

		# add to annotation and prediction lists for Eval
		print ANN_LIST "$gtf_file\n";
		print PRED_LIST "$chr_dirs[$i]/$frag_dirs[$j]/$frag_num.pred.gtf\n";
	    }
	}
	closedir(CHR_DIR);
    }
    
    close (CRF_LIST);
    close (SVM_LIST);
    close (HOLDOUT_LIST);
    close (ANN_LIST);
    close (PRED_LIST);
    
    system("cat $split_dir/crf_list.txt $split_dir/svm_list.txt $split_dir/holdout_list.txt > $split_dir/all_list.txt");
}
else {
    for (my $cv=1; $cv<=$num_cv_sets; $cv++) {
	open (CRF_LIST, ">$split_dir/crf_list.no_set$cv.txt");
	open (SVM_LIST, ">$split_dir/svm_list.no_set$cv.txt");
	open (HOLDOUT_LIST, ">$split_dir/holdout_list.no_set$cv.txt");
	
	for (my $i=0; $i<@chr_dirs; $i++) {
	    print STDERR "Processing $split_dir/$chr_dirs[$i]: ";
	    opendir (CHR_DIR, "$split_dir/$chr_dirs[$i]");
	    my @frag_dirs = grep(/fragment_\d+/, readdir(CHR_DIR));
	    for (my $j=0; $j<@frag_dirs; $j++) {
		print STDERR "\tProcessing $frag_dirs[$j]\n";
		$frag_dirs[$j] =~ /fragment_(\d+)/;
		my $frag_num = $1;
		my $align_file = "$chr_dirs[$i]/$frag_dirs[$j]/$frag_num.align.gz";
		my $estseq_file = "$chr_dirs[$i]/$frag_dirs[$j]/$frag_num.estseq.gz";
		my $gtf_file = "$chr_dirs[$i]/$frag_dirs[$j]/$frag_num.no_set$cv.gtf";
		my $random_val = rand();
		if ($random_val < $CRF_FRACTION) {
		    # add to CRF training list
		    print STDERR "CRF training\n";
		    print CRF_LIST "$align_file\t$estseq_file\t$gtf_file\n";
		}
		elsif ($random_val < $CRF_FRACTION + $SVM_FRACTION) {
		    # add to SVM training list
		    print STDERR "SVM training\n";
		    print SVM_LIST "$align_file\t$estseq_file\t$gtf_file\n";
		}
		else {
		    # add to holdout list
		    print STDERR "Holdout\n";
		    print HOLDOUT_LIST "$align_file\t$estseq_file\t$gtf_file\n";
		}
	    }
	    closedir(CHR_DIR);
	}
	
	close (CRF_LIST);
	close (SVM_LIST);
	close (HOLDOUT_LIST);
    }
}
