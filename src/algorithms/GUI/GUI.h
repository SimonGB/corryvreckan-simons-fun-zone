#ifndef GUI_H
#define GUI_H 1

#include <iostream>
#include "TApplication.h"
#include "TBrowser.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TThread.h"
#include "core/algorithm/Algorithm.h"

#include "objects/Cluster.h"
#include "objects/Pixel.h"
#include "objects/Track.h"

namespace corryvreckan {
    class GUI : public Algorithm {

    public:
        // Constructors and destructors
        GUI(Configuration config, std::vector<Detector*> detectors);
        ~GUI() {}

        // Functions
        void initialise();
        StatusCode run(Clipboard* clipboard);
        void finalise();

        // Plot holders
        std::vector<TCanvas*> canvases;
        std::map<TCanvas*, std::vector<TH1*>> histograms;
        std::map<TH1*, std::string> styles;

        // Add plots and canvases
        void addPlot(TCanvas*, std::string, std::string style = "");
        void addCanvas(TCanvas*);

        // Application to allow display of canvases
        TApplication* app;

        // Misc. member objects
        int eventNumber;
        int updateNumber;
    };
} // namespace corryvreckan
#endif // GUI_H