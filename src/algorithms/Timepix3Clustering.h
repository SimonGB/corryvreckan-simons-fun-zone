#ifndef TIMEPIX3CLUSTERING_H
#define TIMEPIX3CLUSTERING_H 1

#include <iostream>
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "core/algorithm/Algorithm.h"
#include "objects/Cluster.h"
#include "objects/Pixel.h"

namespace corryvreckan {
    class Timepix3Clustering : public Algorithm {

    public:
        // Constructors and destructors
        Timepix3Clustering(Configuration config, std::vector<Detector*> detectors);
        ~Timepix3Clustering() {}

        // Functions
        void initialise();
        StatusCode run(Clipboard* clipboard);
        void finalise();
        void calculateClusterCentre(Cluster*);
        bool touching(Pixel*, Cluster*);
        bool closeInTime(Pixel*, Cluster*);

        double timingCut;
        long long int timingCutInt;
    };
}
#endif // TIMEPIX3CLUSTERING_H
