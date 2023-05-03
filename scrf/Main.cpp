#include "SemiCRF.h"
#include "GTF.h"
#include "Optimizer.h"
#include "SVMTrainer.h"
#include <getopt.h>

#ifdef QUAD_PRECISION
#include <qd/fpu.h>
#endif

#ifdef MULTI
#include <mpi.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#include "MemoryDebug.h"

void printUsageAndExit() {

  string usage = "Usage: ";
  usage += "scrf [options] <parameter file>\n\n";
  usage += "Options:\n";
  usage += "-t --training-list <CRF training list>\t\ttrain CRF on the examples given in the list\n";
  usage += "-s --svm-training-list <SVM training list>\ttrain SVMs on the examples given in the list\n";
  usage += "-e --testing-list <testing list>\t\twhile training, evaluate performance at each iteration on the testing list\n";
  usage += "-p --predict <predict list>\t\tmake predictions for each file in the list\n";
  usage += "-i --initial-guess\t\t\t\tgenerate new HMM-style initial guess from training examples\n";
  usage += "-b --pseudocount <pseudocount>\t\t\tset pseudocount for initial guess (default 0.1)\n";
  usage += "-w --write-gtf <GTF output directory>\t\toutput GTF files of predictions on testing set\n";
  usage += "-r --regularize <regularizer>\t\t\tset the regularization coefficient.  Default is 0.\n";
  usage += "-g --gamma <gamma>\t\t\t\tset gamma parameter (sharpness of sigmoid)\n";
  usage += "-l --lambda <lambda>\t\t\t\tset lambda parameter (label weight)\n";
  usage += "-k --kappa <kappa>\t\t\t\tset kappa parameter (decoding weight)\n";
  usage += "-o --objective <objective>\t\t\tset optimization objective: ll, ea, a, ess, ss, esa, asa\n";
  usage += "-d --decoding <decoding algorithm>\t\tset decoding algorithm: viterbi, mea, or mesa\n";
  usage += "-f --feature-scaling\t\t\t\tscale features by the number of times they appear in the training data\n";
  usage += "-h --stochastic\t\t\t\t\tuse stochastic gradient descent";

  cerr << "Invalid arguments" << endl;
  cerr << usage << endl;
  exit(0);
}

int main(int argc, char** argv) {
#ifdef MULTI
  MPI_Init(&argc, &argv);
#endif

#ifdef _OPENMP
  omp_set_num_threads(2);
#endif

  srand(time(NULL));

  string parameterFile = "";
  string alignFile = "";
  string batchFile = "";
  string crfTrainingList = "";
  string svmTrainingList = "";
  string testingList = "";
  string predictList = "";
  string gtfOutputDir = "NULL";
  bool newInitialGuess = false;
  bool scaling = false;
  weight_t regularizationCoeff = 0;
  objective_t objective = LIKELIHOOD;
  decode_t decode = VITERBI;
  double pseudocount = 0.1;
  weight_t gamma = 5;
  weight_t lambda = 1;
  weight_t kappa = 1;
  bool stochastic = false;

  struct option options[] = {
    {"training-list", required_argument, NULL, 't'},
    {"svm-training-list", required_argument, NULL, 's'},
    {"testing-list", required_argument, NULL, 'e'},
    {"predict", required_argument, NULL, 'p'},
    {"initial-guess", no_argument, NULL, 'i'},
    {"pseudocount", required_argument, NULL, 'b'},
    {"write-gtf", required_argument, NULL, 'w'},
    {"regularize", required_argument, NULL, 'r'},
    {"gamma", required_argument, NULL, 'g'},
    {"lambda", required_argument, NULL, 'l'},
    {"kappa", required_argument, NULL, 'k'},
    {"objective", required_argument, NULL, 'o'},
    {"decoding", required_argument, NULL, 'd'},
    {"feature-scaling", no_argument, NULL, 'f'},
    {"stochastic", no_argument, NULL, 'h'}
  };

  int arg;
  while (true) {
    int option_index;
    arg = getopt_long(argc, argv, "t:s:e:p:ib:w:r:g:l:k:o:d:fh", options, &option_index);
    
    if (arg == -1) break;
    
    switch (arg) {
    case 't':
      crfTrainingList = optarg;
      break;
    case 's':
      svmTrainingList = optarg;
      break;
    case 'e':
      testingList = optarg;
      break;
    case 'p':
      predictList = optarg;
      break;
    case 'i':
      newInitialGuess = true;
      break;
    case 'b':
      pseudocount = atof(optarg);
      break;
    case 'w':
      gtfOutputDir = optarg;
      break;
    case 'r':
      regularizationCoeff = atof(optarg);
      break;
    case 'g':
      gamma = atof(optarg);
      break;
    case 'l':
      lambda = atof(optarg);
      break;
    case 'k':
      kappa = atof(optarg);
      break;
    case 'o':
      if (string(optarg) == "ll")
	objective = LIKELIHOOD;
      else if (string(optarg) == "ea")
	objective = EA;
      else if (string(optarg) == "a")
	objective = A;
      else if (string(optarg) == "ess")
	objective = ESS;
      else if (string(optarg) == "ss")
	objective = SS;
      else if (string(optarg) == "esa")
	objective = ESA;
      else if (string(optarg) == "asa")
	objective = ASA;
      else {
	cerr << "Unknown objective: use \"ll\", \"ea\", \"a\", \"ess\", \"ss\", \"asa\" or \"asa\"" << endl;
	exit(0);
      }
      break;
    case 'd':
      if (string(optarg) == "viterbi")
	decode = VITERBI;
      else if (string(optarg) == "mea")
	decode = MEA;
      else if (string(optarg) == "mesa")
	decode = MESA;
      else {
	cerr << "Unknown decoding algorithm: use \"viterbi\" or \"mea\" or \"mesa\"" << endl;
	exit(0);
      }
      break;
    case 'f':
      scaling = true;
      break;
    case 'h':
      stochastic = true;
      break;
    case '?':
      printUsageAndExit();
      break;
    }
  }
  
  if (optind != argc - 1) {
    printUsageAndExit();
  }

  parameterFile = argv[optind];
  SemiCRF* scrf = new SemiCRF(parameterFile);

#ifdef QUAD_PRECISION
  fpu_fix_start(NULL); //set round-to-double flag on CPU for QD library
#endif
    
  if (predictList != "") {
    //prediction mode
    vector<string> alignFilenames;
    vector<string> estseqFilenames;
    vector<string> gtfFilenames;
    string s;

    ifstream predictListStream(predictList.c_str());
    while (! predictListStream.eof()) {
      predictListStream >> s;
      if (! predictListStream.eof())
	alignFilenames.push_back(s);
      predictListStream >> s;
      if (! predictListStream.eof())
	estseqFilenames.push_back(s);
      predictListStream >> s;
      if (! predictListStream.eof())
	gtfFilenames.push_back(s);
    }

    assert(alignFilenames.size() == estseqFilenames.size());

    int id = 0;
    int numProcs = 1;

#ifdef MULTI
    MPI_Comm_rank(MPI_COMM_WORLD, &id);  // get our id
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs); // get number of process
#endif

    for (int i=id; i<alignFilenames.size(); i+= numProcs) {
      cerr << "Process " << id << " predicting on sequence " << i << ": " << alignFilenames[i] << ", " << estseqFilenames[i] << endl;
      set<string> seqNames;
      scrf->getNeededSequences(seqNames);
      AlignmentSequence alignment(alignFilenames[i], seqNames);
      NumericSequence* numSeq = scrf->createNumericSequence(&alignment);
      ESTSequence* estseq = NULL;
      if  (estseqFilenames[i] != "None")
	estseq = new ESTSequence(estseqFilenames[i]);
      scrf->setSequences(&alignment, numSeq, estseq);
      scrf->precomputeForDP(NULL);
      Segmentation prediction;
      
      if (decode == VITERBI)
	scrf->viterbi(prediction);
      else if (decode == MESA || decode == MEA) {
	//#pragma omp parallel sections
	{
	  //#pragma omp section
	  scrf->forward();
	  //#pragma omp section  
	  scrf->backward();
	}
	if (decode == MEA)
	  scrf->meaDecode(prediction, kappa);
	else if (decode == MESA)
	  scrf->mesaDecode(prediction, kappa);
      }       

      GTF predictionGTF(prediction, *scrf);

      string outputGTFFile = gtfFilenames[i];
      if (outputGTFFile.substr(outputGTFFile.length() - 3, 3) != "gtf")
	fatalError("GTF filename does not end with \"gtf\"");
      outputGTFFile.replace(outputGTFFile.length() - 3, 3, "pred.gtf");
      predictionGTF.write(outputGTFFile);
    }
  }
  else if (crfTrainingList != "") {
    //training mode
    if (gtfOutputDir != "NULL" && testingList == "")
      printUsageAndExit();

    //first train SVMs, if any
    if (svmTrainingList != "" && newInitialGuess) {
      SVMTrainer* svmTrainer = new SVMTrainer(scrf, svmTrainingList, crfTrainingList);
      svmTrainer->train();
      delete svmTrainer;
    }

    //delete scrf;
    //exit(0);

    //now optimize CRF
    Optimizer opt(scrf, crfTrainingList, testingList, objective, decode);
    opt.optimize(newInitialGuess, pseudocount, gtfOutputDir, scaling, gamma, 
		 lambda, kappa, stochastic, regularizationCoeff);
  }
  else
    printUsageAndExit();

#ifdef MULTI
  MPI_Finalize();
#endif

  return 0;
}
