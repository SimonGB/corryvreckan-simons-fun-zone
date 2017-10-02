#ifndef ALIGNMENT_H
#define ALIGNMENT_H 1

// ROOT includes
#include "Minuit2/Minuit2Minimizer.h"
#include "Math/Functor.h"
#include "TError.h"
// Local includes
#include "Algorithm.h"
#include "Cluster.h"
#include "Track.h"

class Alignment : public Algorithm {
  
public:
  // Constructors and destructors
  Alignment(bool);
  ~Alignment(){}

  // Functions
  void initialise(Parameters*);
  StatusCode run(Clipboard*);
  void finalise();
  
  // Member variables
  Tracks m_alignmenttracks;
  int nIterations;
  int m_numberOfTracksForAlignment;

};

#endif // ALIGNMENT_H