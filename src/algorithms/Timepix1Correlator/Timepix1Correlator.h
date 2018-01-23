#ifndef Timepix1Correlator_H
#define Timepix1Correlator_H 1

#include <iostream>
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "core/algorithm/Algorithm.h"

namespace corryvreckan {
    /** @ingroup Algorithms
     */
    class Timepix1Correlator : public Algorithm {

    public:
        // Constructors and destructors
        Timepix1Correlator(Configuration config, std::vector<Detector*> detectors);
        ~Timepix1Correlator() {}

        // Functions
        void initialise();
        StatusCode run(Clipboard* clipboard);
        void finalise();

        // Histograms for several devices
        std::map<std::string, TH2F*> hitmaps;
        std::map<std::string, TH2F*> hitmapsGlobal;
        std::map<std::string, TH1F*> clusterSize;
        std::map<std::string, TH1F*> clustersPerEvent;
        std::map<std::string, TH1F*> correlationPlotsX;
        std::map<std::string, TH1F*> correlationPlotsY;

        // Member variables
        int m_eventNumber;
    };
} // namespace corryvreckan
#endif // Timepix1Correlator_H
