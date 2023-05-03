CC=ccache g++-3.4
CFLAGS = -g -pg -O3 -march=athlon64 --param large-function-growth=10000 --param inline-unit-growth=10000 --param max-inline-insns-single=10000 -Winline
LINKFLAGS = -lm -lz

SRCS = \
	AlignmentSequence.cpp\
	DNAKmerArrayFS.cpp\
	DNAKmerFS.cpp\
	DNAKmerPairArrayFS.cpp\
	DNAKmerPairFS.cpp\
	DPMatrix.cpp\
	ESTPositionFS.cpp\
	ESTTransitionFS.cpp\
	ESTSequence.cpp\
	FeatureSet.cpp\
	GTF.cpp\
	GTFEvaluation.cpp\
	GTFFeature.cpp\
	GTFTranscript.cpp\
	KmerIterator.cpp\
	LengthFS.cpp\
	Main.cpp\
	MaskingFeature.cpp\
	Minimizer.cpp\
	NumericFS.cpp\
	NumericSequence.cpp\
	Optimizer.cpp\
	RegularizationSelector.cpp\
	Scaler.cpp\
	Segment.cpp\
	Segmentation.cpp\
	SemiCRF.cpp\
	Sequence.cpp\
	SequenceFeatureSet.cpp\
	State.cpp\
	StateType.cpp\
	svm.cpp\
	SVMFS.cpp\
	SVMTrainer.cpp\
	TransitionFeature.cpp\
	TransitionType.cpp\
	WorkUnit.cpp\

OBJS = $(SRCS:%.cpp=%.o) 

all: main

main: $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) $(LINKFLAGS) -o scrf

intel:
	make all CC="icpc" CFLAGS="-O3 -xP -openmp -vec-report0"

intelmulti:
	make all CC="mpic++" LAMHCP="icpc" OMPI_MPICXX="icpc" CFLAGS="-O3 -xP -openmp -vec-report0 -DMULTI"

intelquad:
	make all CC="mpic++" LAMHCP="icpc" OMPI_MPICXX="icpc" CFLAGS="-O3 -xP -openmp -vec-report0 -DQUAD_PRECISION" LINKFLAGS="-lm -lqd -lz"

multi:
	make all CC="ccache mpic++" OMPI_MPICXX="g++" CFLAGS="-g -O3 -DMULTI -march=athlon64 --param large-function-growth=10000 --param inline-unit-growth=10000 --param max-inline-insns-single=10000 -Winline"

quad:
	make all CFLAGS="-DQUAD_PRECISION -O3 -march=athlon64 --param large-function-growth=10000 --param inline-unit-growth=10000 --param max-inline-insns-single=10000 -Winline" LINKFLAGS="-lm -lqd -lz"

debugmulti:
	make all CC="ccache mpiCC" CFLAGS="-g -ggdb -DMULTI -march=athlon64 --param large-function-growth=10000 --param inline-unit-growth=10000 --param max-inline-insns-single=10000 -Winline"

debug:
	make all CFLAGS="-g -ggdb -march=athlon64 --param large-function-growth=10000 --param inline-unit-growth=10000 --param max-inline-insns-single=10000 -Winline"


%.o:	%.cpp
	$(CC) -c $(CFLAGS) $<

clean:
	rm -f *.o core
	rm -f scrf
