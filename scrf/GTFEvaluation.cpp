#include "GTFEvaluation.h"
#include <set>
#include <iomanip>

void GTFEvaluation::init() {
  transcriptStats.resize(numStatTypes, 0);
  combinedExonStats.resize(numStatTypes, 0);
  cdsPositionalStats.resize(numStatTypes, 0);
  individualExonStats.resize(numStatTypes, vector<int>(numExonTypes, 0));
  exonTypePositionalStats.resize(numStatTypes, vector<int>(numExonTypes, 0));
}

GTFEvaluation::GTFEvaluation() {
  init();
}

GTFEvaluation::GTFEvaluation(GTF& annotation, GTF& prediction) {
  init();

  /* first make lists of all annotated and predicted exons that does
     not include duplicates (e.g., the same exon from different transcripts) */
  set<GTFFeature> predictedExons;
  set<GTFFeature> annotatedExons;
  vector<GTFFeature> correctExons;

  for (map<string, GTFTranscript>::iterator iter = prediction.transcripts.begin();
       iter != prediction.transcripts.end(); ++iter) {
    for (int i=0; i<iter->second.features.size(); i++) {
      if (iter->second.features[i].featureName == "CDS" ||
	  iter->second.features[i].featureName == "5UTR" ||
	  iter->second.features[i].featureName == "3UTR")
	predictedExons.insert(iter->second.features[i]);
    }
  }
  for (map<string, GTFTranscript>::iterator iter = annotation.transcripts.begin();
       iter != annotation.transcripts.end(); ++iter) {
    for (int i=0; i<iter->second.features.size(); i++) {
      if (iter->second.features[i].featureName == "CDS" ||
	  iter->second.features[i].featureName == "5UTR" ||
	  iter->second.features[i].featureName == "3UTR") {
	annotatedExons.insert(iter->second.features[i]);
      }
    }
  }

  insert_iterator<vector<GTFFeature> > ii(correctExons, correctExons.begin());
  set_intersection(predictedExons.begin(), predictedExons.end(),
		   annotatedExons.begin(), annotatedExons.end(),
		   ii);

  /* go through the lists and collect some stats */
  pos_t maxPos = 0;
  for (set<GTFFeature>::iterator iter = predictedExons.begin();
       iter != predictedExons.end(); ++iter) {
    if (iter->type == INITIAL || iter->type == INTERNAL ||
	iter->type == TERMINAL || iter->type == SINGLE)
      combinedExonStats[PREDICTED]++;
    individualExonStats[PREDICTED][iter->type]++;
    if (iter->end > maxPos)
      maxPos = iter->end;
  }
  for (set<GTFFeature>::iterator iter = annotatedExons.begin();
       iter != annotatedExons.end(); ++iter) {
    if (iter->type == INITIAL || iter->type == INTERNAL ||
	iter->type == TERMINAL || iter->type == SINGLE)
      combinedExonStats[ANNOTATED]++;
    individualExonStats[ANNOTATED][iter->type]++;
    if (iter->end > maxPos)
      maxPos = iter->end;
  }
  for (vector<GTFFeature>::iterator iter = correctExons.begin();
       iter != correctExons.end(); ++iter) {
    if (iter->type == INITIAL || iter->type == INTERNAL ||
	iter->type == TERMINAL || iter->type == SINGLE)
      combinedExonStats[CORRECT]++;
    individualExonStats[CORRECT][iter->type]++;
  }

  /* get transcript level stats */
  transcriptStats[ANNOTATED] = annotation.transcripts.size();
  transcriptStats[PREDICTED] = prediction.transcripts.size();
  for (map<string, GTFTranscript>::iterator annIter = annotation.transcripts.begin();
       annIter != annotation.transcripts.end(); ++annIter) {
    for (map<string, GTFTranscript>::iterator predIter = prediction.transcripts.begin();
	 predIter != prediction.transcripts.end(); ++predIter) {
      //consider a predicted transcript correct if the coding region matches an annotation
      bool match = true;

      //make sure the two transcripts have the same number of CDS features
      int annotatedCDSFeatures = 0;
      int predictedCDSFeatures = 0;
      for (int i=0; i<annIter->second.features.size(); i++) {
	if (annIter->second.features[i].featureName == "CDS")
	  annotatedCDSFeatures++;
      }
      for (int i=0; i<predIter->second.features.size(); i++) {
	if (predIter->second.features[i].featureName == "CDS")
	  predictedCDSFeatures++;
      }
      if (annotatedCDSFeatures != predictedCDSFeatures)
	match = false;

      //make sure every feature in the annotation is found in the prediction
      for (int i=0; i<annIter->second.features.size(); i++) {
	if (annIter->second.features[i].featureName == "CDS") {
	  bool featureMatch = false;
	  for (int j=0; j<predIter->second.features.size(); j++) {
	    if (annIter->second.features[i] == predIter->second.features[j])
	      featureMatch = true;
	  }
	  if (! featureMatch)
	    match = false;
	}
      }
      if (match) 
	transcriptStats[CORRECT]++;
    }
  }

  /* look for positional overlap */  
  vector<bool> cdsPosPredicted(maxPos+1, false);
  vector<vector<bool> > cdsTypePosPredicted(maxPos+1, vector<bool>(numExonTypes, false));
  vector<bool> cdsPosAnnotated(maxPos+1, false);
  vector<vector<bool> > cdsTypePosAnnotated(maxPos+1, vector<bool>(numExonTypes, false));
  
  for (set<GTFFeature>::iterator iter = predictedExons.begin();
       iter != predictedExons.end(); ++iter) {
    for (int i=iter->start; i<= iter->end; i++) {
      cdsPosPredicted[i] = true;
      cdsTypePosPredicted[i][iter->type] = true;
    }
  }
  for (set<GTFFeature>::iterator iter = annotatedExons.begin();
       iter != annotatedExons.end(); ++iter) {
    for (int i=iter->start; i<= iter->end; i++) {
      cdsPosAnnotated[i] = true;
      cdsTypePosAnnotated[i][iter->type] = true;
    }
  }

  for (int i=0; i<maxPos; i++) {
    if (cdsPosPredicted[i])
      cdsPositionalStats[PREDICTED]++;
    if (cdsPosAnnotated[i])
      cdsPositionalStats[ANNOTATED]++;
    if (cdsPosPredicted[i] && cdsPosAnnotated[i])
      cdsPositionalStats[CORRECT]++;
    for (int j=0; j<numExonTypes; j++) {
      if (cdsTypePosPredicted[i][j])
	exonTypePositionalStats[PREDICTED][j]++;
      if (cdsTypePosAnnotated[i][j])
	exonTypePositionalStats[ANNOTATED][j]++;
      if (cdsTypePosPredicted[i][j] && cdsTypePosAnnotated[i][j])
	exonTypePositionalStats[CORRECT][j]++;
    }
  }
}

void GTFEvaluation::pack(vector<int>& packed) {
  for (int i=0; i<numStatTypes; i++) {
    packed.push_back(transcriptStats[i]);
    packed.push_back(combinedExonStats[i]);
    packed.push_back(cdsPositionalStats[i]);
    for (int j=0; j<numExonTypes; j++) {
      packed.push_back(individualExonStats[i][j]);
      packed.push_back(exonTypePositionalStats[i][j]);
    }
  }
}

void GTFEvaluation::unpack(vector<int>& packed) {
  transcriptStats.clear();
  combinedExonStats.clear();
  cdsPositionalStats.clear();
  for (int i=0; i<numStatTypes; i++) {
    individualExonStats[i].clear();
    exonTypePositionalStats[i].clear();
  }			

  int index = 0;
  for (int i=0; i<numStatTypes; i++) {
    transcriptStats.push_back(packed[index++]);
    combinedExonStats.push_back(packed[index++]);
    cdsPositionalStats.push_back(packed[index++]);
    for (int j=0; j<numExonTypes; j++) {
      individualExonStats[i].push_back(packed[index++]);
      exonTypePositionalStats[i].push_back(packed[index++]);
    }
  }
}

void GTFEvaluation::operator+= (GTFEvaluation& other) {
  for (int i=0; i<numStatTypes; i++) {
    transcriptStats[i] += other.transcriptStats[i];
    combinedExonStats[i] += other.combinedExonStats[i];
    cdsPositionalStats[i] += other.cdsPositionalStats[i];
    for (int j=0; j<numExonTypes; j++) {
      individualExonStats[i][j] += other.individualExonStats[i][j];
      exonTypePositionalStats[i][j] += other.exonTypePositionalStats[i][j];
    }
  }
}

ostream& operator<<(ostream& os, GTFEvaluation& eval) {
  os << setiosflags(ios::fixed | ios::showpoint) << setprecision(1);
  
  os << "Transcript\t"
     << eval.transcriptStats[CORRECT] * 100.0 / eval.transcriptStats[ANNOTATED]
     << "\t" << eval.transcriptStats[CORRECT] * 100.0 / eval.transcriptStats[PREDICTED]
     << endl;
  os << "Exon\t\t" 
     << eval.combinedExonStats[CORRECT] * 100.0 / eval.combinedExonStats[ANNOTATED] 
     << "\t" << eval.combinedExonStats[CORRECT] * 100.0 / eval.combinedExonStats[PREDICTED] 
     << endl;
  os << "CDS Positional\t" 
     << eval.cdsPositionalStats[CORRECT] * 100.0 / eval.cdsPositionalStats[ANNOTATED] 
     << "\t" << eval.cdsPositionalStats[CORRECT] * 100.0 / eval.cdsPositionalStats[PREDICTED] 
     << endl;
  for (int i=0; i<numExonTypes; i++) {
    os << exonTypeNames[i] << "\t";
    if (exonTypeNames[i].length() < 8)
      os << "\t";
    os << eval.individualExonStats[CORRECT][i] * 100.0 / eval.individualExonStats[ANNOTATED][i]
       << "\t" << eval.individualExonStats[CORRECT][i] * 100.0 / eval.individualExonStats[PREDICTED][i]
       << "\t";
    os << eval.exonTypePositionalStats[CORRECT][i] * 100.0 / eval.exonTypePositionalStats[ANNOTATED][i]
       << "\t" << eval.exonTypePositionalStats[CORRECT][i] * 100.0 / eval.exonTypePositionalStats[PREDICTED][i]
       << endl;
  }
  os << endl;
  os << "Transcript Count\t" << eval.transcriptStats[ANNOTATED] << "\t" << eval.transcriptStats[PREDICTED] << endl;
  os << "Mean Transcript Length\t" << 1.0 * eval.cdsPositionalStats[ANNOTATED] / eval.transcriptStats[ANNOTATED]
     << "\t" << 1.0 * eval.cdsPositionalStats[PREDICTED] / eval.transcriptStats[PREDICTED] << endl;
  os << "Mean Exons Per Gene\t" << 1.0 * eval.combinedExonStats[ANNOTATED] / eval.transcriptStats[ANNOTATED]
     << "\t" << 1.0 * eval.combinedExonStats[PREDICTED] / eval.transcriptStats[PREDICTED] << endl;
  for (int i=0; i<numExonTypes; i++) {
    os << exonTypeNames[i] << " Count\t\t" << eval.individualExonStats[ANNOTATED][i] 
       << "\t" << eval.individualExonStats[PREDICTED][i] << endl;
  }
  for (int i=0; i<numExonTypes; i++) {
    os << "Mean " << exonTypeNames[i] << " Length" 
       << "\t" << 1.0 * eval.exonTypePositionalStats[ANNOTATED][i] / eval.individualExonStats[ANNOTATED][i]
       << "\t" << 1.0 * eval.exonTypePositionalStats[PREDICTED][i] / eval.individualExonStats[PREDICTED][i]
       << endl;
  }

  os << setprecision(7);
   
  return os;
}


