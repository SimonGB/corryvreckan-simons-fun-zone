#include "OnlineMonitor.h"
#include <TVirtualPadEditor.h>
#include <regex>

using namespace corryvreckan;
using namespace std;

OnlineMonitor::OnlineMonitor(Configuration config, std::vector<Detector*> detectors)
    : Module(std::move(config), std::move(detectors)) {
    canvasTitle = m_config.get<std::string>("canvasTitle", "Corryvreckan Testbeam Monitor");
    updateNumber = m_config.get<int>("update", 500);

    // Set up overview plots:
    canvas_overview = m_config.getMatrix<std::string>("Overview",
                                                      {{"BasicTracking/trackChi2", ""},
                                                       {"TestAlgorithm/clusterTot_%REFERENCE%", ""},
                                                       {"TestAlgorithm/hitmap_%REFERENCE%", "colz"},
                                                       {"BasicTracking/residualsX_%REFERENCE%", ""}});

    // Set up individual plots for the DUT
    canvas_dutplots = m_config.getMatrix<std::string>("DUTPlots",
                                                      {{"Clicpix2EventLoader/hitMap", "colz"},
                                                       {"Clicpix2EventLoader/hitMapDiscarded", "colz"},
                                                       {"Clicpix2EventLoader/pixelToT", ""},
                                                       {"Clicpix2EventLoader/pixelToA", ""},
                                                       {"Clicpix2EventLoader/pixelCnt", "log"},
                                                       {"Clicpix2EventLoader/pixelsPerFrame", "log"},
                                                       {"DUTAnalysis/clusterTotAssociated", ""},
                                                       {"DUTAnalysis/associatedTracksVersusTime", ""}});
}

void OnlineMonitor::initialise() {

    // TApplication keeps the canvases persistent
    app = new TApplication("example", 0, 0);

    // Make the GUI
    gui = new GuiDisplay();

    // Make the main window object and set the attributes
    gui->m_mainFrame = new TGMainFrame(gClient->GetRoot(), 800, 600);
    gui->buttonMenu = new TGHorizontalFrame(gui->m_mainFrame, 800, 50);
    gui->canvas = new TRootEmbeddedCanvas("canvas", gui->m_mainFrame, 800, 600);
    gui->m_mainFrame->AddFrame(gui->canvas, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 10, 10, 10, 10));
    gui->m_mainFrame->SetCleanup(kDeepCleanup);
    gui->m_mainFrame->DontCallClose();

    // Add canvases and histograms

    //=== Overview canvas
    AddButton("Overview", "OverviewCanvas");
    AddPlots("OverviewCanvas", canvas_overview);

    //=== Track canvas
    AddButton("Tracking", "TrackCanvas");
    AddHisto("TrackCanvas", "/corryvreckan/BasicTracking/trackChi2");
    AddHisto("TrackCanvas", "/corryvreckan/BasicTracking/trackAngleX");

    //=== Per detector canvases
    AddButton("HitMaps", "HitmapCanvas");
    AddButton("Residuals", "ResidualCanvas");
    AddButton("EventTimes", "EventTimeCanvas");
    AddButton("CorrelationsX", "CorrelationXCanvas");
    AddButton("CorrelationsY", "CorrelationYCanvas");
    AddButton("CorrelationsX2D", "CorrelationX2DCanvas");
    AddButton("CorrelationsY2D", "CorrelationY2DCanvas");
    AddButton("ChargeDistributions", "ChargeDistributionCanvas");
    AddButton("DUTPlots", "DUTCanvas");

    // Per detector histograms
    for(auto& detector : get_detectors()) {
        string detectorID = detector->name();

        string hitmap = "/corryvreckan/TestAlgorithm/hitmap_" + detectorID;
        AddHisto("HitmapCanvas", hitmap, "colz");

        string chargeHisto = "/corryvreckan/Timepix3Clustering/clusterTot_" + detectorID;
        AddHisto("ChargeDistributionCanvas", chargeHisto);

        string eventTimeHisto = "/corryvreckan/TestAlgorithm/eventTimes_" + detectorID;
        AddHisto("EventTimeCanvas", eventTimeHisto);

        string correlationXHisto = "/corryvreckan/TestAlgorithm/correlationX_" + detectorID;
        AddHisto("CorrelationXCanvas", correlationXHisto);

        string correlationX2DHisto = "/corryvreckan/TestAlgorithm/correlationX_2Dlocal_" + detectorID;
        AddHisto("CorrelationX2DCanvas", correlationX2DHisto, "colz");

        string correlationYHisto = "/corryvreckan/TestAlgorithm/correlationY_" + detectorID;
        AddHisto("CorrelationYCanvas", correlationYHisto);

        string correlationY2DHisto = "/corryvreckan/TestAlgorithm/correlationY_2Dlocal_" + detectorID;
        AddHisto("CorrelationY2DCanvas", correlationY2DHisto, "colz");

        // Hisograms below not available for DUTs:
        if(detectorID == m_config.get<std::string>("DUT")) {
            continue;
        }

        string residualHisto = "/corryvreckan/BasicTracking/residualsX_" + detectorID;
        AddHisto("ResidualCanvas", residualHisto);
    }

    AddPlots("DUTCanvas", canvas_dutplots);

    // Set up the main frame before drawing

    // Exit button
    string exitButton = "StopMonitoring";
    gui->buttons[exitButton] = new TGTextButton(gui->buttonMenu, exitButton.c_str());
    gui->buttonMenu->AddFrame(gui->buttons[exitButton], new TGLayoutHints(kLHintsLeft, 10, 10, 10, 10));
    gui->buttons[exitButton]->Connect("Pressed()", "corryvreckan::GuiDisplay", gui, "Exit()");

    // Main frame resizing
    gui->m_mainFrame->AddFrame(gui->buttonMenu, new TGLayoutHints(kLHintsLeft, 10, 10, 10, 10));
    gui->m_mainFrame->SetWindowName(canvasTitle.c_str());
    gui->m_mainFrame->MapSubwindows();
    gui->m_mainFrame->Resize(gui->m_mainFrame->GetDefaultSize());

    // Draw the main frame
    gui->m_mainFrame->MapWindow();

    // Plot the overview tab (if it exists)
    if(gui->histograms["OverviewCanvas"].size() != 0) {
        gui->Display(const_cast<char*>(std::string("OverviewCanvas").c_str()));
    }

    gui->canvas->GetCanvas()->Paint();
    gui->canvas->GetCanvas()->Update();
    gSystem->ProcessEvents();

    // Initialise member variables
    eventNumber = 0;
}

StatusCode OnlineMonitor::run(Clipboard* clipboard) {

    // Draw all histograms
    if(eventNumber % updateNumber == 0) {
        gui->canvas->GetCanvas()->Paint();
        gui->canvas->GetCanvas()->Update();
        eventNumber++;
    }
    gSystem->ProcessEvents();

    // Get the tracks from the clipboard
    Tracks* tracks = (Tracks*)clipboard->get("tracks");
    if(tracks == NULL)
        return Success;

    // Otherwise increase the event number
    eventNumber++;
    return Success;
}

void OnlineMonitor::finalise() {

    LOG(DEBUG) << "Analysed " << eventNumber << " events";
}

void OnlineMonitor::AddPlots(std::string canvas_name, Matrix<std::string> canvas_plots) {
    for(auto plot : canvas_plots) {
        if(plot.size() != 2) {
            continue;
        }

        // Do we need to plot with a LogY scale?
        bool log_scale = (plot.back().find("log") != std::string::npos) ? true : false;

        // Replace other placeholders and add histogram
        std::string name = std::regex_replace(plot.front(), std::regex("%DUT%"), m_config.get<std::string>("DUT"));
        name = std::regex_replace(name, std::regex("%REFERENCE%"), m_config.get<std::string>("reference"));

        // Do we have a detector placeholder?
        if(name.find("%DETECTOR%") != std::string::npos) {
            LOG(DEBUG) << "Adding plot " << name << " for all detectors.";
            for(auto& detector : get_detectors()) {
                AddHisto(canvas_name,
                         "/corryvreckan/" + std::regex_replace(name, std::regex("%DETECTOR%"), detector->name()),
                         plot.back(),
                         log_scale);
            }
        } else {
            // Single histogram only.
            AddHisto(canvas_name, "/corryvreckan/" + name, plot.back(), log_scale);
        }
    }
}

void OnlineMonitor::AddHisto(string canvasName, string histoName, string style, bool logy) {

    // Add "corryvreckan" namespace:
    // histoName = "/corryvreckan/" + histoName;

    TH1* histogram = (TH1*)gDirectory->Get(histoName.c_str());
    if(histogram) {
        gui->histograms[canvasName].push_back((TH1*)gDirectory->Get(histoName.c_str()));
        gui->styles[gui->histograms[canvasName].back()] = style;
        gui->logarithmic[gui->histograms[canvasName].back()] = logy;
    } else {
        LOG(WARNING) << "Histogram " << histoName << " does not exist";
    }
}

void OnlineMonitor::AddButton(string buttonName, string canvasName) {
    gui->buttons[buttonName] = new TGTextButton(gui->buttonMenu, buttonName.c_str());
    gui->buttonMenu->AddFrame(gui->buttons[buttonName], new TGLayoutHints(kLHintsLeft, 10, 10, 10, 10));
    string command = "Display(=\"" + canvasName + "\")";
    LOG(INFO) << "Connecting button with command " << command.c_str();
    gui->buttons[buttonName]->Connect("Pressed()", "corryvreckan::GuiDisplay", gui, command.c_str());
}
