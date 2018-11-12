#include "EventLoaderEUDAQ.h"
#include "eudaq/PluginManager.hh"

#include <algorithm>

using namespace corryvreckan;
using namespace std;

EventLoaderEUDAQ::EventLoaderEUDAQ(Configuration config, std::vector<std::shared_ptr<Detector>> detectors)
    : Module(std::move(config), std::move(detectors)), m_longID(true) {

    m_filename = m_config.get<std::string>("file_name");
    m_longID = m_config.get<bool>("long_detector_id", true);
}

void EventLoaderEUDAQ::initialise() {

    // Create new file reader:
    try {
        reader = new eudaq::FileReader(m_filename, "");
    } catch(...) {
        throw ModuleError("Unable to read input file \"" + m_filename + "\"");
    }

    // Initialise member variables
    m_eventNumber = 0;
}

StatusCode EventLoaderEUDAQ::run(std::shared_ptr<Clipboard> clipboard) {

    // Read next event from EUDAQ reader:
    const eudaq::DetectorEvent& evt = reader->Event();
    LOG(TRACE) << evt;

    if(evt.IsBORE()) {
        // Process begin-of-run
        LOG(DEBUG) << "Found BORE";
        try {
            eudaq::PluginManager::Initialize(evt);
        } catch(const eudaq::Exception&) {
            throw ModuleError("Unknown event types encountered");
        }
    } else if(evt.IsEORE()) {
        LOG(INFO) << "Found EORE";
    } else {
        eudaq::StandardEvent sevt = eudaq::PluginManager::ConvertToStandard(evt);

        for(size_t iplane = 0; iplane < sevt.NumPlanes(); ++iplane) {
            // Get EUDAQ StandardPlane
            const eudaq::StandardPlane& plane = sevt.GetPlane(iplane);

            // Build Corryvreckan detector ID and check if this detector should be read:
            std::string detectorID =
                (m_longID ? (plane.Sensor() + "_") : "plane") + std::to_string(plane.ID()); // : (plane.ID()));

            if(!has_detector(detectorID)) {
                LOG(DEBUG) << "Skipping unknown detector " << detectorID;
                continue;
            }

            auto detector = get_detector(detectorID);

            // Make a new container for the data
            Pixels* deviceData = new Pixels();
            for(unsigned int ipix = 0; ipix < plane.HitPixels(); ++ipix) {
                auto col = static_cast<int>(plane.GetX(ipix));
                auto row = static_cast<int>(plane.GetY(ipix));

                // Check if this pixel is masked
                if(detector->masked(col, row)) {
                    LOG(TRACE) << "Detector " << detectorID << ": pixel " << col << "," << row << " masked";
                    continue;
                }

                Pixel* pixel = new Pixel(detectorID, row, col, static_cast<int>(plane.GetPixel(ipix)));
                pixel->setCharge(plane.GetPixel(ipix));
                pixel->timestamp(m_eventNumber);
                deviceData->push_back(pixel);
            }

            // Store on clipboard
            clipboard->put(detectorID, "pixels", reinterpret_cast<Objects*>(deviceData));
        }
    }

    // Increment event counter
    m_eventNumber++;

    // Advance to next event if possible, otherwise end this run:
    if(!reader->NextEvent()) {
        LOG(INFO) << "No more events in data stream.";
        return StatusCode::Failure;
    };

    // Return value telling analysis to keep running
    return StatusCode::Success;
}

void EventLoaderEUDAQ::finalise() {

    LOG(DEBUG) << "Read " << m_eventNumber << " events";
}