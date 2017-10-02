#ifndef TIMEPIX3CLUSTERING_H
#define TIMEPIX3CLUSTERING_H 1

#include "Pixel.h"
#include "Cluster.h"
#include "Algorithm.h"
#include <iostream>
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"

class Timepix3Clustering : public Algorithm {
  
public:
  // Constructors and destructors
  Timepix3Clustering(bool);
  ~Timepix3Clustering(){}

  // Functions
  void initialise(Parameters*);
  StatusCode run(Clipboard*);
  void finalise();
  void calculateClusterCentre(Cluster*);
  bool touching(Pixel*,Cluster*);
  bool closeInTime(Pixel*,Cluster*);
  
  double timingCut;
  long long int timingCutInt;

};

#endif // TIMEPIX3CLUSTERING_H