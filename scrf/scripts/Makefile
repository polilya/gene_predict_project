CC=g++
CFLAGS = -O3 

all: maf_to_align make_estseq mask_from_rm_output axts_to_align

maf_to_align: maf_to_align.o AlignmentBlock.o
	$(CC) $(CFLAGS) maf_to_align.o AlignmentBlock.o -lboost_regex -lz -o maf_to_align

make_estseq: make_estseq.o
	$(CC) $(CFLAGS) make_estseq.o -o make_estseq

mask_from_rm_output: mask_from_rm_output.o
	$(CC) $(CFLAGS) mask_from_rm_output.o -o mask_from_rm_output

axts_to_align: axts_to_align.o
	$(CC) $(CFLAGS) axts_to_align.o -o axts_to_align

%.o:	%.cpp
	$(CC) -c $(CFLAGS) $<

clean:
	rm -f *.o core
	rm -f maf_to_align
	rm -f make_estseq
	rm -f mask_from_rm_output
	rm -f axts_to_align
