/**
 * @file
 * @brief Implementation of module Clustering4D
 *
 * @copyright Copyright (c) 2017-2020 CERN and the Corryvreckan authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 */

#include "Clustering4D.h"

using namespace corryvreckan;
using namespace std;

Clustering4D::Clustering4D(Configuration config, std::shared_ptr<Detector> detector)
    : Module(std::move(config), detector), m_detector(detector) {

    // Backwards compatibilty: also allow timing_cut to be used for time_cut_abs
    m_config.setAlias("time_cut_abs", "timing_cut", true);

    if(m_config.count({"time_cut_rel", "time_cut_abs"}) > 1) {
        throw InvalidCombinationError(
            m_config, {"time_cut_rel", "time_cut_abs"}, "Absolute and relative time cuts are mutually exclusive.");
    } else if(m_config.has("time_cut_abs")) {
        timeCut = m_config.get<double>("time_cut_abs");
    } else {
        timeCut = m_config.get<double>("time_cut_rel", 3.0) * m_detector->getTimeResolution();
    }
    neighbourRadiusRow = m_config.get<int>("neighbour_radius_row", 1);
    neighbourRadiusCol = m_config.get<int>("neighbour_radius_col", 1);
    chargeWeighting = m_config.get<bool>("charge_weighting", true);
    useEarliestPixel = m_config.get<bool>("use_earliest_pixel", false);
}

void Clustering4D::initialise() {

    // Cluster plots
    std::string title = m_detector->getName() + " Cluster size;cluster size;events";
    clusterSize = new TH1F("clusterSize", title.c_str(), 100, -0.5, 99.5);
    title = m_detector->getName() + " Cluster seed charge;cluster seed charge [e];events";
    clusterSeedCharge = new TH1F("clusterSeedCharge", title.c_str(), 256, -0.5, 255.5);
    title = m_detector->getName() + " Cluster Width - Rows;cluster width [rows];events";
    clusterWidthRow = new TH1F("clusterWidthRow", title.c_str(), 25, -0.5, 24.5);
    title = m_detector->getName() + " Cluster Width - Columns;cluster width [columns];events";
    clusterWidthColumn = new TH1F("clusterWidthColumn", title.c_str(), 100, -0.5, 99.5);
    title = m_detector->getName() + " Cluster Charge;cluster charge [e];events";
    clusterCharge = new TH1F("clusterCharge", title.c_str(), 5000, -0.5, 49999.5);
    title = m_detector->getName() + " Cluster Charge (1px clusters);cluster charge [e];events";
    clusterCharge_1px = new TH1F("clusterCharge_1px", title.c_str(), 256, -0.5, 255.5);
    title = m_detector->getName() + " Cluster Charge (2px clusters);cluster charge [e];events";
    clusterCharge_2px = new TH1F("clusterCharge_2px", title.c_str(), 256, -0.5, 255.5);
    title = m_detector->getName() + " Cluster Charge (3px clusters);cluster charge [e];events";
    clusterCharge_3px = new TH1F("clusterCharge_3px", title.c_str(), 256, -0.5, 255.5);
    title = m_detector->getName() + " Cluster Position (Global);x [mm];y [mm];events";
    clusterPositionGlobal = new TH2F("clusterPositionGlobal",
                                     title.c_str(),
                                     400,
                                     -m_detector->getSize().X() / 1.5,
                                     m_detector->getSize().X() / 1.5,
                                     400,
                                     -m_detector->getSize().Y() / 1.5,
                                     m_detector->getSize().Y() / 1.5);
    title = m_detector->getName() + " Cluster Position (Local);x [px];y [px];events";
    clusterPositionLocal = new TH2F("clusterPositionLocal",
                                    title.c_str(),
                                    m_detector->nPixels().X(),
                                    -0.5,
                                    m_detector->nPixels().X() - 0.5,
                                    m_detector->nPixels().Y(),
                                    -0.5,
                                    m_detector->nPixels().Y() - 0.5);

    title = ";cluster timestamp [ns]; # events";
    clusterTimes = new TH1F("clusterTimes", title.c_str(), 3e6, 0, 3e9);
    title = m_detector->getName() + " Cluster multiplicity;clusters;events";
    clusterMultiplicity = new TH1F("clusterMultiplicity", title.c_str(), 50, -0.5, 49.5);
    title = m_detector->getName() +
            " pixel - cluster timestamp;ts_{pixel} - ts_{cluster} [ns] (all pixels from cluster (if clusterSize>1));events";
    pixelTimeMinusClusterTime = new TH1F("pixelTimeMinusClusterTime", title.c_str(), 1000, -0.5, 999.5);

    // Get resolution in time of detector and calculate time cut to be applied
    LOG(DEBUG) << "Time cut to be applied for " << m_detector->getName() << " is "
               << Units::display(timeCut, {"ns", "us", "ms"});
}

// Sort function for pixels from low to high times
bool Clustering4D::sortByTime(Pixel* pixel1, Pixel* pixel2) {
    return (pixel1->timestamp() < pixel2->timestamp());
}

StatusCode Clustering4D::run(std::shared_ptr<Clipboard> clipboard) {

    // Get the pixels
    auto pixels = clipboard->getData<Pixel>(m_detector->getName());
    if(pixels == nullptr) {
        LOG(DEBUG) << "Detector " << m_detector->getName() << " does not have any pixels on the clipboard";
        return StatusCode::Success;
    }
    LOG(DEBUG) << "Picked up " << pixels->size() << " pixels for device " << m_detector->getName();

    // Sort the pixels from low to high timestamp
    std::sort(pixels->begin(), pixels->end(), sortByTime);
    size_t totalPixels = pixels->size();

    // Make the cluster storage
    auto deviceClusters = std::make_shared<ClusterVector>();

    // Keep track of which pixels are used
    map<Pixel*, bool> used;

    // Start to cluster
    for(size_t iP = 0; iP < pixels->size(); iP++) {
        Pixel* pixel = (*pixels)[iP];

        // Check if pixel is used
        if(used[pixel]) {
            continue;
        }

        // Make the new cluster object
        Cluster* cluster = new Cluster();
        LOG(DEBUG) << "==== New cluster";

        // Keep adding hits to the cluster until no more are found
        cluster->addPixel(pixel);
        double clusterTime = pixel->timestamp();
        used[pixel] = true;
        LOG(DEBUG) << "Adding pixel: " << pixel->column() << "," << pixel->row();
        size_t nPixels = 0;
        while(cluster->size() != nPixels) {

            nPixels = cluster->size();
            // Loop over all pixels
            for(size_t iNeighbour = (iP + 1); iNeighbour < totalPixels; iNeighbour++) {
                Pixel* neighbour = (*pixels)[iNeighbour];
                // Check if they are compatible in time with the cluster pixels
                if(abs(neighbour->timestamp() - clusterTime) > timeCut)
                    break;

                // Check if they have been used
                if(used[neighbour])
                    continue;

                // Check if they are touching cluster pixels
                if(!touching(neighbour, cluster))
                    continue;

                // Add to cluster
                cluster->addPixel(neighbour);
                clusterTime = neighbour->timestamp();
                used[neighbour] = true;
                LOG(DEBUG) << "Adding pixel: " << neighbour->column() << "," << neighbour->row() << " time "
                           << Units::display(neighbour->timestamp(), {"ns", "us", "s"});
            }
        }

        // Finalise the cluster and save it
        calculateClusterCentre(cluster);

        // Fill cluster histograms
        clusterSize->Fill(static_cast<double>(cluster->size()));
        clusterWidthRow->Fill(cluster->rowWidth());
        clusterWidthColumn->Fill(cluster->columnWidth());
        clusterCharge->Fill(cluster->charge());
        if(cluster->size() == 1) {
            clusterCharge_1px->Fill(cluster->charge());
        } else if(cluster->size() == 2) {
            clusterCharge_2px->Fill(cluster->charge());
        } else if(cluster->size() == 3) {
            clusterCharge_3px->Fill(cluster->charge());
        }
        clusterSeedCharge->Fill(cluster->getSeedPixel()->charge());
        clusterPositionGlobal->Fill(cluster->global().x(), cluster->global().y());
        clusterPositionLocal->Fill(cluster->column(), cluster->row());
        clusterTimes->Fill(static_cast<double>(Units::convert(cluster->timestamp(), "ns")));

        // to check that cluster timestamp = earliest pixel timestamp
        if(cluster->size() > 1) {
            for(auto& px : cluster->pixels()) {
                pixelTimeMinusClusterTime->Fill(
                    static_cast<double>(Units::convert(px->timestamp() - cluster->timestamp(), "ns")));
            }
        }

        deviceClusters->push_back(cluster);
    }

    clusterMultiplicity->Fill(static_cast<double>(deviceClusters->size()));

    // Put the clusters on the clipboard
    clipboard->putData(deviceClusters, m_detector->getName());
    LOG(DEBUG) << "Made " << deviceClusters->size() << " clusters for device " << m_detector->getName();

    return StatusCode::Success;
}

// Check if a pixel touches any of the pixels in a cluster
bool Clustering4D::touching(Pixel* neighbour, Cluster* cluster) {

    bool Touching = false;

    for(auto pixel : cluster->pixels()) {
        int row_distance = abs(pixel->row() - neighbour->row());
        int col_distance = abs(pixel->column() - neighbour->column());

        if(row_distance <= neighbourRadiusRow && col_distance <= neighbourRadiusCol) {
            if(row_distance > 1 || col_distance > 1) {
                cluster->setSplit(true);
            }
            Touching = true;
            break;
        }
    }
    return Touching;
}

// Check if a pixel is close in time to the pixels of a cluster
bool Clustering4D::closeInTime(Pixel* neighbour, Cluster* cluster) {

    bool CloseInTime = false;

    auto pixels = cluster->pixels();
    for(auto& px : pixels) {

        double timeDifference = abs(neighbour->timestamp() - px->timestamp());
        if(timeDifference < timeCut)
            CloseInTime = true;
    }
    return CloseInTime;
}

void Clustering4D::calculateClusterCentre(Cluster* cluster) {

    LOG(DEBUG) << "== Making cluster centre";
    // Empty variables to calculate cluster position
    double column(0), row(0), charge(0), maxcharge(0);
    double column_sum(0), column_sum_chargeweighted(0);
    double row_sum(0), row_sum_chargeweighted(0);
    bool found_charge_zero = false;

    // Get the pixels on this cluster
    auto pixels = cluster->pixels();
    string detectorID = pixels.front()->detectorID();
    double timestamp = pixels.front()->timestamp();
    LOG(DEBUG) << "- cluster has " << pixels.size() << " pixels";

    // Loop over all pixels
    for(auto& pixel : pixels) {
        // If charge == 0 (use epsilon to avoid errors in floating-point arithmetics):
        if(pixel->charge() < std::numeric_limits<double>::epsilon()) {
            // apply arithmetic mean if a pixel has zero charge
            found_charge_zero = true;
        }
        charge += pixel->charge();

        // We need both column_sum and column_sum_chargeweighted
        // as we don't know a priori if there will be a pixel with
        // charge==0 such that we have to fall back to the arithmetic mean.
        column_sum += pixel->column();
        row_sum += pixel->row();
        column_sum_chargeweighted += (pixel->column() * pixel->charge());
        row_sum_chargeweighted += (pixel->row() * pixel->charge());

        // If charge is available: cluster timestamp = pixel with largestcharge:
        if(!found_charge_zero && useEarliestPixel) {
            if(pixel->charge() > maxcharge) {
                timestamp = pixel->timestamp();
                maxcharge = pixel->charge();
            }
        } else { // use earliest pixel
            if(pixel->timestamp() < timestamp) {
                timestamp = pixel->timestamp();
            }
        }
    }

    if(chargeWeighting && !found_charge_zero) {
        // Charge-weighted centre-of-gravity for cluster centre:
        // (here it's safe to divide by the charge as it cannot be zero due to !found_charge_zero)
        column = column_sum_chargeweighted / charge;
        row = row_sum_chargeweighted / charge;
    } else {
        // Arithmetic cluster centre:
        column = column_sum / static_cast<double>(cluster->size());
        row = row_sum / static_cast<double>(cluster->size());
    }

    if(detectorID != m_detector->getName()) {
        // Should never happen...
        return;
    }

    // Calculate local cluster position
    auto positionLocal = m_detector->getLocalPosition(column, row);

    // Calculate global cluster position
    auto positionGlobal = m_detector->localToGlobal(positionLocal);

    LOG(DEBUG) << "- cluster col, row, time : " << column << "," << row << ","
               << Units::display(timestamp, {"ns", "us", "ms"});
    // Set the cluster parameters
    cluster->setColumn(column);
    cluster->setRow(row);
    cluster->setCharge(charge);

    // Set uncertainty on position from intrinstic detector spatial resolution:
    cluster->setError(m_detector->getSpatialResolution());

    cluster->setTimestamp(timestamp);
    cluster->setDetectorID(detectorID);
    cluster->setClusterCentre(positionGlobal);
    cluster->setClusterCentreLocal(positionLocal);
}
