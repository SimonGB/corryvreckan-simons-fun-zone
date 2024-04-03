/**
 * @file
 * @brief Implementation of module AlignmentMillepede
 *
 * @copyright Copyright (c) 2017-2022 CERN and the Corryvreckan authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 * SPDX-License-Identifier: MIT
 */

#include <algorithm>
#include <cmath>
#include <iomanip>

// Local
#include "AlignmentMillepede.h"
#include "objects/Cluster.hpp"

using namespace corryvreckan;
using namespace std;

//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
AlignmentMillepede::AlignmentMillepede(Configuration& config, std::vector<std::shared_ptr<Detector>> detectors)
    : Module(config, std::move(detectors)) {

    config_.setDefault<bool>("exclude_dut", false);
    config_.setDefaultArray<bool>("dofs", {true, true, false, true, true, true});
    config_.setDefault<size_t>("iterations", 5);
    config_.setDefault<double>("residual_cut", 0.05);
    config_.setDefault<double>("residual_cut_init", 0.6);
    config_.setDefault<int>("number_of_stddev", 0);
    config_.setDefault<double>("convergence", 0.00001);
    config_.setDefaultArray<double>("sigmas", {0.05, 0.05, 0.5, 0.005, 0.005, 0.005});

    m_excludeDUT = config_.get<bool>("exclude_dut");
    m_dofs = config_.getArray<bool>("dofs");
    m_nIterations = config_.get<size_t>("iterations");
    m_rescut = config_.get<double>("residual_cut");
    m_rescut_init = config_.get<double>("residual_cut_init");
    m_nstdev = config_.get<int>("number_of_stddev");
    m_convergence = config_.get<double>("convergence");
    m_sigmas = config_.getArray<double>("sigmas");

    if(m_dofs.size() != 6) {
        throw InvalidValueError(config_, "dofs", "Invalid number of degrees of freedom.");
    }

    // Check that we're not in a variable-alignment situation:
    for(auto& detector : get_regular_detectors(false)) {
        if(detector->hasVariableAlignment()) {
            throw ModuleError("Cannot perform alignment procedure with variable alignment of detector \"" +
                              detector->getName() + "\"");
        }
    }
}

//=============================================================================
// Destructor
//=============================================================================
AlignmentMillepede::~AlignmentMillepede() {}

//=============================================================================
// Initialization
//=============================================================================
void AlignmentMillepede::initialize() {

    // Renumber the planes in Millepede, ignoring masked planes.
    unsigned int index = 0;
    for(const auto& det : get_regular_detectors(!m_excludeDUT)) {
        m_millePlanes[det->getName()] = index;
        ++index;
    }

    // Print the degrees of freedom.
    LOG(INFO) << "Using the following degrees of freedom:";
    const std::vector<std::string> labels = {
        "Translation X", "Translation Y", "Translation Z", "Rotation X", "Rotation Y", "Rotation Z"};
    for(unsigned int i = 0; i < 6; ++i) {
        const std::string on = m_dofs[i] ? "ON" : "OFF";
        LOG(INFO) << labels[i] << "\t" << on;
    }
}

// During run, just pick up tracks and save them till the end
StatusCode AlignmentMillepede::run(const std::shared_ptr<Clipboard>& clipboard) {

    // Get the tracks
    auto tracks = clipboard->getData<Track>();
    TrackVector alignmenttracks;
    std::map<std::string, std::vector<Cluster*>> alignmentclusters;

    // Make a local copy and store it
    for(auto& track : tracks) {
        alignmenttracks.push_back(track);
        auto clusters = track->getClusters();
        for(auto& cluster : track->getClusters()) {
            alignmentclusters[cluster->detectorID()].push_back(cluster);
        }
    }

    // Store all tracks we want for alignment on the permanent storage:
    clipboard->putPersistentData(alignmenttracks);
    // Copy the objects of all track clusters on the clipboard to persistent storage:
    for(auto& clusters : alignmentclusters) {
        clipboard->copyToPersistentData(clusters.second, clusters.first);
    }

    return StatusCode::Success;
}

//=============================================================================
// Main alignment function
//=============================================================================
void AlignmentMillepede::finalize(const std::shared_ptr<ReadonlyClipboard>& clipboard) {

    LOG(INFO) << "Millepede alignment";
    auto alignmenttracks = clipboard->getPersistentData<Track>();

    size_t nPlanes = num_regular_detectors(!m_excludeDUT);
    LOG(INFO) << "Aligning " << nPlanes << " planes";
    const size_t nParameters = 6 * nPlanes;
    for(unsigned int iteration = 0; iteration < m_nIterations; ++iteration) {
        // Define the constraint equations.
        setConstraints(nPlanes);
        const double startfact = 100.;
        // Initialise all matrices and vectors.
        reset(nPlanes, startfact);
        LOG(INFO) << "Feeding Millepede with " << alignmenttracks.size() << " tracks...";
        // Feed Millepede with tracks.
        unsigned int nSkipped = 0;
        unsigned int nOutliers = 0;
        for(auto& track : alignmenttracks) {
            if(track->getNClusters() != nPlanes) {
                ++nSkipped;
                continue;
            }
            if(!putTrack(track.get(), nPlanes)) {
                ++nOutliers;
            }
        }
        if(nSkipped > 0) {
            LOG(INFO) << "Skipped " << nSkipped << " tracks with less than " << nPlanes << " clusters.";
        }
        if(nOutliers > 0) {
            LOG(INFO) << "Rejected " << nOutliers << " outlier tracks.";
        }
        // Do the global fit.
        LOG(INFO) << "Determining global parameters...";
        if(!fitGlobal()) {
            LOG(ERROR) << "Global fit failed.";
            break;
        }
        // Calculate the convergence metric (sum of misalignments).
        double converg = 0.;
        for(size_t i = 0; i < nParameters; ++i) {
            converg += fabs(m_dparm[i]);
        }
        converg /= static_cast<double>(nParameters);
        LOG(INFO) << "Convergence: " << converg;
        // Update the module positions and orientations.
        LOG(INFO) << "Updating geometry...";
        updateGeometry();

        // Update the cluster coordinates based on the new geometry.
        for(auto& track : alignmenttracks) {
            for(auto& cluster : track->getClusters()) {
                auto detectorID = cluster->detectorID();
                auto detector = get_detector(detectorID);
                ROOT::Math::XYZPoint pLocal(cluster->local().x(), cluster->local().y(), 0.);
                const auto pGlobal = detector->localToGlobal(pLocal);
                cluster->setClusterCentre(pGlobal);
            }
        }
        if(converg < m_convergence)
            break;
    }
}

//=============================================================================
// Setup the constraint equations.
//=============================================================================
void AlignmentMillepede::setConstraints(const size_t nPlanes) {

    // Calculate the mean z-position.
    double avgz = 0.;
    for(const auto& det : get_regular_detectors(!m_excludeDUT)) {
        avgz += det->displacement().Z();
    }
    avgz /= static_cast<double>(nPlanes);
    // Calculate the variance.
    double varz = 0.0;
    for(const auto& det : get_regular_detectors(!m_excludeDUT)) {
        const double dz = det->displacement().Z() - avgz;
        varz += dz * dz;
    }
    varz /= static_cast<double>(nPlanes);

    // Define the 9 constraints equations according to the requested geometry.
    const size_t nParameters = 6 * nPlanes;
    std::vector<double> ftx(nParameters, 0.);
    std::vector<double> fty(nParameters, 0.);
    std::vector<double> ftz(nParameters, 0.);
    std::vector<double> frx(nParameters, 0.);
    std::vector<double> fry(nParameters, 0.);
    std::vector<double> frz(nParameters, 0.);
    std::vector<double> fscaz(nParameters, 0.);
    std::vector<double> shearx(nParameters, 0.);
    std::vector<double> sheary(nParameters, 0.);

    m_constraints.clear();
    for(const auto& det : get_regular_detectors(!m_excludeDUT)) {
        const unsigned int i = m_millePlanes[det->getName()];
        const double sz = (det->displacement().Z() - avgz) / varz;
        ftx[i] = 1.0;
        fty[i + nPlanes] = 1.0;
        ftz[i + 2 * nPlanes] = 1.0;
        /*
          Do NOT change these signs, unless we figure out how to do these
          constraints properly.
         */
        frx[i + 3 * nPlanes] = i >= 4 ? 1.0 : -1.0;
        fry[i + 4 * nPlanes] = i >= 4 ? 1.0 : -1.0;
        frz[i + 5 * nPlanes] = 1.0;
        shearx[i] = sz;
        sheary[i + nPlanes] = sz;
        fscaz[i + 2 * nPlanes] = sz;
    }

    const std::vector<bool> constraints = {true, true, true, true, false, false, true, true, true};
    //  Put the constraints information in the basket.
    if(constraints[0] && m_dofs[0])
        addConstraint(ftx, 0.0);
    if(constraints[1] && m_dofs[0])
        addConstraint(shearx, 0.);
    if(constraints[2] && m_dofs[1])
        addConstraint(fty, 0.0);
    if(constraints[3] && m_dofs[1])
        addConstraint(sheary, 0.);
    if(constraints[4] && m_dofs[2])
        addConstraint(ftz, 0.0);
    // if (constraints[5] && m_dofs[2]) addConstraint(fscaz, 0.0);
    if(constraints[6] && m_dofs[3])
        addConstraint(frx, 0.0);
    if(constraints[7] && m_dofs[4])
        addConstraint(fry, 0.0);
    if(constraints[8] && m_dofs[5])
        addConstraint(frz, 0.0);
}

//=============================================================================
// Define a single constraint equation.
//=============================================================================
void AlignmentMillepede::addConstraint(const std::vector<double>& dercs, const double rhs) {

    Constraint constraint;
    // Set the right-hand side (Lagrange multiplier value, sum of equation).
    constraint.rhs = rhs;
    // Set the constraint equation coefficients.
    constraint.coefficients = dercs;
    m_constraints.push_back(constraint);
}

//=============================================================================
// Add the equations for one track to the matrix
//=============================================================================
bool AlignmentMillepede::putTrack(Track* track, const size_t nPlanes) {

    std::vector<Equation> equations;
    const size_t nParameters = 6 * nPlanes;
    // Global derivatives
    std::vector<double> dergb(nParameters, 0.);
    // Global derivatives non linearly related to residual
    std::vector<double> dernl(nParameters, 0.);
    // Local parameter indices associated with non-linear derivatives
    std::vector<int> dernli(nParameters, 0);
    // Track slopes
    std::vector<double> dernls(nParameters, 0.);

    /// Refit the track for the reference states.
    track->fit();
    const double tx = track->getState(track->getClusters().front()->detectorID()).X();
    const double ty = track->getState(track->getClusters().front()->detectorID()).Y();

    // Iterate over each cluster on the track.
    for(auto& cluster : track->getClusters()) {
        if(!has_detector(cluster->detectorID())) {
            continue;
        }
        auto detector = get_detector(cluster->detectorID());
        const auto normal = detector->normal();
        double nx = normal.x() / normal.z();
        double ny = normal.y() / normal.z();
        const double xg = cluster->global().x();
        const double yg = cluster->global().y();
        const double zg = cluster->global().z();
        // Calculate quasi-local coordinates.
        const double zl = zg - detector->displacement().Z();
        const double xl = xg - detector->displacement().X();
        const double yl = yg - detector->displacement().Y();

        std::vector<double> nonlinear = {nx, ny, 1., -yl + zl * ny, -xl + zl * nx, xl * ny - yl * nx};
        const double den = 1. + tx * nx + ty * ny;
        for(auto& a : nonlinear)
            a /= den;
        // Get the errors on the measured x and y coordinates.
        const double errx = cluster->errorX();
        const double erry = cluster->errorY();
        // Get the internal plane index in Millepede.
        const unsigned int plane = m_millePlanes[detector->getName()];
        // Set the local derivatives for the X equation.
        std::vector<double> derlc = {1., zg, 0., 0.};
        // Set the global derivatives (see LHCb-2005-101) for the X equation.
        std::fill(dergb.begin(), dergb.end(), 0.);
        std::fill(dernl.begin(), dernl.end(), 0.);
        std::fill(dernli.begin(), dernli.end(), 0);
        std::fill(dernls.begin(), dernls.end(), 0.);
        if(m_dofs[0])
            dergb[plane] = -1.;
        if(m_dofs[4])
            dergb[4 * nPlanes + plane] = -zl;
        if(m_dofs[5])
            dergb[5 * nPlanes + plane] = yl;
        for(unsigned int i = 0; i < 6; ++i) {
            if(!m_dofs[i])
                continue;
            dergb[i * nPlanes + plane] += tx * nonlinear[i];
            dernl[i * nPlanes + plane] = nonlinear[i];
            dernli[i * nPlanes + plane] = 1;
            dernls[i * nPlanes + plane] = tx;
        }
        // Store the X equation.
        addEquation(equations, derlc, dergb, dernl, dernli, dernls, xg, errx);
        // Set the local derivatives for the Y equation.
        derlc = {0., 0., 1., zg};
        // Set the global derivatives (see LHCb-2005-101) for the Y equation.
        std::fill(dergb.begin(), dergb.end(), 0.);
        std::fill(dernl.begin(), dernl.end(), 0.);
        std::fill(dernli.begin(), dernli.end(), 0);
        std::fill(dernls.begin(), dernls.end(), 0.);
        if(m_dofs[1])
            dergb[nPlanes + plane] = -1.;
        if(m_dofs[3])
            dergb[3 * nPlanes + plane] = -zl;
        if(m_dofs[5])
            dergb[5 * nPlanes + plane] = -xl;
        for(unsigned int i = 0; i < 6; ++i) {
            if(!m_dofs[i])
                continue;
            dergb[i * nPlanes + plane] += ty * nonlinear[i];
            dernl[i * nPlanes + plane] = nonlinear[i];
            dernli[i * nPlanes + plane] = 3;
            dernls[i * nPlanes + plane] = ty;
        }
        // Store the Y equation.
        addEquation(equations, derlc, dergb, dernl, dernli, dernls, yg, erry);
    }
    // Vector containing the track parameters
    std::vector<double> trackParams(2 * m_nalc + 2, 0.);
    // Fit the track.
    const unsigned int iteration = 1;
    const bool ok = fitTrack(equations, trackParams, false, iteration);
    if(ok)
        m_equations.push_back(equations);
    return ok;
}

//=============================================================================
// Store the parameters for one measurement
//=============================================================================
void AlignmentMillepede::addEquation(std::vector<Equation>& equations,
                                     const std::vector<double>& derlc,
                                     const std::vector<double>& dergb,
                                     const std::vector<double>& dernl,
                                     const std::vector<int>& dernli,
                                     const std::vector<double>& slopes,
                                     const double rmeas,
                                     const double sigma) {

    if(sigma <= 0.) {
        LOG(ERROR) << "Invalid cluster error (" << sigma << ")";
        return;
    }
    Equation equation;
    equation.rmeas = rmeas;
    equation.weight = 1. / (sigma * sigma);
    // Add non-zero local derivatives and corresponding indices.
    equation.derL.clear();
    equation.indL.clear();
    for(unsigned int i = 0; i < m_nalc; ++i) {
        if(derlc[i] == 0.)
            continue;
        equation.indL.push_back(static_cast<int>(i));
        equation.derL.push_back(derlc[i]);
    }
    // Add non-zero global derivatives and corresponding indices.
    equation.derG.clear();
    equation.indG.clear();
    equation.derNL.clear();
    equation.indNL.clear();
    equation.slopes.clear();
    const size_t nG = dergb.size();
    for(size_t i = 0; i < nG; ++i) {
        if(dergb[i] == 0.)
            continue;
        equation.indG.push_back(static_cast<int>(i));
        equation.derG.push_back(dergb[i]);
        equation.derNL.push_back(dernl[i]);
        equation.indNL.push_back(dernli[i]);
        equation.slopes.push_back(slopes[i]);
    }
    // Add the equation to the list.
    equations.push_back(std::move(equation));
}

//=============================================================================
// Track fit (local fit)
//=============================================================================
bool AlignmentMillepede::fitTrack(const std::vector<Equation>& equations,
                                  std::vector<double>& trackParams,
                                  const bool singlefit,
                                  const unsigned int iteration) {

    std::vector<double> blvec(m_nalc, 0.);
    std::vector<std::vector<double>> clmat(m_nalc, std::vector<double>(m_nalc, 0.));

    // First loop: local track fit
    for(const auto& equation : equations) {
        double rmeas = equation.rmeas;
        // Suppress the global part (only relevant with iterations).
        const size_t nG = equation.derG.size();
        for(size_t i = 0; i < nG; ++i) {
            const size_t j = static_cast<size_t>(equation.indG[i]);
            rmeas -= equation.derG[i] * m_dparm[j];
        }
        const double w = equation.weight;
        // Fill local matrix and vector.
        const size_t nL = equation.derL.size();
        for(size_t i = 0; i < nL; ++i) {
            const size_t j = static_cast<size_t>(equation.indL[i]);
            blvec[j] += w * rmeas * equation.derL[i];
            LOG(DEBUG) << "blvec[" << j << "] = " << blvec[j];

            // Symmetric matrix, don't bother k > j coefficients.
            for(unsigned int k = 0; k <= i; ++k) {
                const size_t ik = static_cast<size_t>(equation.indL[k]);
                clmat[j][ik] += w * equation.derL[i] * equation.derL[k];
                LOG(DEBUG) << "clmat[" << j << "][" << ik << "] = " << clmat[j][ik];
            }
        }
    }

    // Local parameter matrix is completed, now invert to solve.
    // Rank: number of non-zero diagonal elements.
    int rank = invertMatrixLocal(clmat, blvec, static_cast<size_t>(m_nalc));
    // Store the track parameters and errors.
    for(unsigned int i = 0; i < m_nalc; ++i) {
        trackParams[2 * i] = blvec[i];
        trackParams[2 * i + 1] = sqrt(fabs(clmat[i][i]));
    }

    // Second loop: residual calculation.
    double chi2 = 0.0;
    int ndf = 0;
    for(const auto& equation : equations) {
        double rmeas = equation.rmeas;
        // Suppress global and local terms.
        const size_t nL = equation.derL.size();
        for(size_t i = 0; i < nL; ++i) {
            const size_t j = static_cast<size_t>(equation.indL[i]);
            rmeas -= equation.derL[i] * blvec[j];
        }
        const size_t nG = equation.derG.size();
        for(size_t i = 0; i < nG; ++i) {
            const size_t j = static_cast<size_t>(equation.indG[i]);
            rmeas -= equation.derG[i] * m_dparm[j];
        }
        // Now rmeas contains the residual value.
        LOG(DEBUG) << "Residual value: " << rmeas;
        // Reject the track if the residual is too large (outlier).
        const double rescut = iteration <= 1 ? m_rescut_init : m_rescut;
        if(fabs(rmeas) >= rescut) {
            LOG(DEBUG) << "Reject track due to residual cut in iteration " << iteration;
            return false;
        }
        chi2 += equation.weight * rmeas * rmeas;
        ++ndf;
    }
    ndf -= rank;
    const bool printDetails = false;
    if(printDetails) {
        LOG(INFO);
        LOG(INFO) << "Local track fit (rank " << rank << ")";
        LOG(INFO) << " Result of local fit :      (index/parameter/error)";
        for(unsigned int i = 0; i < m_nalc; ++i) {
            LOG(INFO) << std::setprecision(6) << std::fixed << std::setw(20) << i << "   /   " << std::setw(10) << blvec[i]
                      << "   /   " << sqrt(clmat[i][i]);
        }
        LOG(INFO) << "Final chi square / degrees of freedom: " << chi2 << " / " << ndf;
    }

    // Stop here if just updating the track parameters.
    if(singlefit)
        return true;

    if(m_nstdev != 0 && ndf > 0) {
        const double chi2PerNdf = chi2 / ndf;
        const double cut = chi2Limit(m_nstdev, ndf) * m_cfactr;
        if(chi2 > cut) {
            // Reject the track.
            LOG(DEBUG) << "Rejected track because chi2 / ndof (" << chi2PerNdf << ") is larger than " << cut;
            return false;
        }
    }
    // Store the chi2 and number of degrees of freedom.
    trackParams[2 * m_nalc] = ndf;
    trackParams[2 * m_nalc + 1] = chi2;

    // Local operations are finished. Track is accepted.
    // Third loop: update the global parameters (other matrices).
    m_clcmat.assign(m_nagb, std::vector<double>(m_nalc, 0.));
    unsigned int nagbn = 0;
    std::vector<int> indnz(m_nagb, -1);
    std::vector<int> indbk(m_nagb, 0);
    for(const auto& equation : equations) {
        double rmeas = equation.rmeas;
        // Suppress the global part.
        const size_t nG = equation.derG.size();
        for(size_t i = 0; i < nG; ++i) {
            const size_t j = static_cast<size_t>(equation.indG[i]);
            rmeas -= equation.derG[i] * m_dparm[j];
        }
        const double w = equation.weight;
        // First of all, the global/global terms.
        for(size_t i = 0; i < nG; ++i) {
            const size_t j = static_cast<size_t>(equation.indG[i]);
            m_bgvec[j] += w * rmeas * equation.derG[i];
            LOG(DEBUG) << "bgvec[" << j << "] = " << m_bgvec[j];
            for(size_t k = 0; k < nG; ++k) {
                const size_t n = static_cast<size_t>(equation.indG[k]);
                m_cgmat[j][n] += w * equation.derG[i] * equation.derG[k];
                LOG(DEBUG) << "cgmat[" << j << "][" << n << "] = " << m_cgmat[j][n];
            }
        }
        // Now we have also rectangular matrices containing global/local terms.
        const size_t nL = equation.derL.size();
        for(size_t i = 0; i < nG; ++i) {
            const size_t j = static_cast<size_t>(equation.indG[i]);
            // Index of index.
            int ik = indnz[j];
            if(ik == -1) {
                // New global variable.
                indnz[j] = static_cast<int>(nagbn);
                indbk[nagbn] = static_cast<int>(j);
                ik = static_cast<int>(nagbn);
                ++nagbn;
            }
            // Now fill the rectangular matrix.
            for(size_t k = 0; k < nL; ++k) {
                const size_t ij = static_cast<size_t>(equation.indL[k]);
                m_clcmat[static_cast<size_t>(ik)][ij] += w * equation.derG[i] * equation.derL[k];
                LOG(DEBUG) << "clcmat[" << ik << "][" << ij << "] = " << m_clcmat[static_cast<size_t>(ik)][ij];
            }
        }
    }

    // Third loop is finished, now we update the correction matrices.
    multiplyAVAt(clmat, m_clcmat, m_corrm, m_nalc, nagbn);
    multiplyAX(m_clcmat, blvec, m_corrv, m_nalc, nagbn);
    for(size_t i = 0; i < nagbn; ++i) {
        const size_t j = static_cast<size_t>(indbk[i]);
        m_bgvec[j] -= m_corrv[i];
        for(size_t k = 0; k < nagbn; ++k) {
            const size_t ik = static_cast<size_t>(indbk[k]);
            m_cgmat[j][ik] -= m_corrm[i][k];
        }
    }
    return true;
}

//=============================================================================
// Update the module positions and orientations.
//=============================================================================
void AlignmentMillepede::updateGeometry() {
    auto nPlanes = num_regular_detectors(!m_excludeDUT);

    for(const auto& det : get_regular_detectors(!m_excludeDUT)) {
        auto plane = m_millePlanes[det->getName()];

        det->update(XYZPoint(det->displacement().X() + m_dparm[plane + 0 * nPlanes],
                             det->displacement().Y() + m_dparm[plane + 1 * nPlanes],
                             det->displacement().Z() + m_dparm[plane + 2 * nPlanes]),
                    det->rotation() +
                        XYZVector(m_dparm[plane + 3 * nPlanes], m_dparm[plane + 4 * nPlanes], m_dparm[plane + 5 * nPlanes]));
    }
    /*
        const unsigned int nPlanes = m_nPlanes - m_maskedPlanes.size();
        for(auto im = m_modules.cbegin(), end = m_modules.cend(); im != end; ++im) {

            Tb::SmallRotation((*im)->rotX(),
                              m_dparm[plane + 3 * nPlanes],
                              (*im)->rotY(),
                              m_dparm[plane + 4 * nPlanes],
                              m_dparm[plane + 5 * nPlanes],
                              rx,
                              ry,
                              rz);
            (*im)->setAlignment(tx, ty, tz, (*im)->rotX() - rx, (*im)->rotY() + ry, (*im)->rotZ() - rz, 0., 0., 0., 0., 0.,
       0.);
        }
        */
}

//=============================================================================
// Initialise the vectors and arrays.
//=============================================================================
bool AlignmentMillepede::reset(const size_t nPlanes, const double startfact) {

    // Reset the list of track equations.
    m_equations.clear();
    // Set the number of global derivatives.
    m_nagb = static_cast<unsigned int>(6 * nPlanes);
    // Reset matrices and vectors.
    const unsigned int nRows = m_nagb + static_cast<unsigned int>(m_constraints.size());
    m_bgvec.resize(nRows, 0.);
    m_cgmat.assign(nRows, std::vector<double>(nRows, 0.));
    m_corrv.assign(m_nagb, 0.);
    m_corrm.assign(m_nagb, std::vector<double>(m_nagb, 0.));
    m_dparm.assign(m_nagb, 0.);

    // Define the sigmas for each parameter.
    m_fixed.assign(m_nagb, true);
    m_psigm.assign(m_nagb, 0.);
    for(unsigned int i = 0; i < 6; ++i) {
        if(!m_dofs[i])
            continue;
        for(size_t j = i * nPlanes; j < (i + 1) * nPlanes; ++j) {
            m_fixed[j] = false;
            m_psigm[j] = m_sigmas[i];
        }
    }
    // Fix modules if requested.
    for(auto it = m_fixedPlanes.cbegin(); it != m_fixedPlanes.cend(); ++it) {
        LOG(INFO) << "You are fixing module " << (*it);
        // TODO: check if this is the "millePlanes" index or the regular one.
        const unsigned ndofs = m_fix_all ? 6 : 3;
        for(unsigned int i = 0; i < ndofs; ++i) {
            m_fixed[(*it) + i * nPlanes] = true;
            m_psigm[(*it) + i * nPlanes] = 0.;
        }
    }

    // Set the chi2 / ndof cut used in the local track fit.
    // Iterations are stopped when the cut factor reaches the ref. value (1).
    m_cfactref = 1.;
    m_cfactr = std::max(1., startfact);

    return true;
}

//=============================================================================
//
//=============================================================================
bool AlignmentMillepede::fitGlobal() {

    m_diag.assign(m_nagb, 0.);
    std::vector<double> bgvecPrev(m_nagb, 0.);
    std::vector<double> trackParams(2 * m_nalc + 2, 0.);
    const size_t nTracks = m_equations.size();
    std::vector<std::vector<double>> localParams(nTracks, std::vector<double>(m_nalc, 0.));

    const unsigned int nMaxIterations = 10;
    unsigned int iteration = 1;
    unsigned int nGoodTracks = static_cast<unsigned int>(nTracks);
    while(iteration <= nMaxIterations) {
        if(nGoodTracks == 0) {
            LOG(ERROR) << "No tracks to work with after outlier rejection.";
            return false;
        }

        LOG(INFO) << "Iteration " << iteration << " (using " << nGoodTracks << " tracks)";

        // Save the diagonal elements.
        for(unsigned int i = 0; i < m_nagb; ++i) {
            m_diag[i] = m_cgmat[i][i];
        }

        unsigned int nFixed = 0;
        for(unsigned int i = 0; i < m_nagb; ++i) {
            if(m_fixed[i]) {
                // Fixed global parameter.
                ++nFixed;
                for(unsigned int j = 0; j < m_nagb; ++j) {
                    // Reset row and column.
                    m_cgmat[i][j] = 0.0;
                    m_cgmat[j][i] = 0.0;
                }
            } else {
                m_cgmat[i][i] += 1. / (m_psigm[i] * m_psigm[i]);
            }
        }
        // Add the constraints equations.
        size_t nRows = static_cast<size_t>(m_nagb);
        for(const auto& constraint : m_constraints) {
            double sum = constraint.rhs;
            for(unsigned int j = 0; j < m_nagb; ++j) {
                if(m_psigm[j] == 0.) {
                    m_cgmat[nRows][j] = 0.0;
                    m_cgmat[j][nRows] = 0.0;
                } else {
                    m_cgmat[nRows][j] = nGoodTracks * constraint.coefficients[j];
                    m_cgmat[j][nRows] = m_cgmat[nRows][j];
                }
                sum -= constraint.coefficients[j] * m_dparm[j];
            }
            m_cgmat[nRows][nRows] = 0.0;
            m_bgvec[nRows] = nGoodTracks * sum;
            ++nRows;
        }

        double cor = 0.0;
        if(iteration > 1) {
            for(unsigned int j = 0; j < m_nagb; ++j) {
                for(unsigned int i = 0; i < m_nagb; ++i) {
                    if(m_fixed[i])
                        continue;
                    cor += bgvecPrev[j] * m_cgmat[j][i] * bgvecPrev[i];
                    if(i == j) {
                        cor -= bgvecPrev[i] * bgvecPrev[i] / (m_psigm[i] * m_psigm[i]);
                    }
                }
            }
        }
        LOG(DEBUG) << " Final corr. is " << cor;

        // Do the matrix inversion.
        const int rank = invertMatrix(m_cgmat, m_bgvec, nRows);
        // Update the global parameters values.
        for(unsigned int i = 0; i < m_nagb; ++i) {
            m_dparm[i] += m_bgvec[i];
            bgvecPrev[i] = m_bgvec[i];
            LOG(DEBUG) << "bgvec[" << i << "] = " << m_bgvec[i];
            LOG(DEBUG) << "dparm[" << i << "] = " << m_dparm[i];
            LOG(DEBUG) << "cgmat[" << i << "][" << i << "] = " << m_cgmat[i][i];
            LOG(DEBUG) << "err = " << sqrt(fabs(m_cgmat[i][i]));
            LOG(DEBUG) << "cgmat * diag = " << std::setprecision(5) << m_cgmat[i][i] * m_diag[i];
        }
        if(nRows - nFixed - static_cast<unsigned int>(rank) != 0) {
            LOG(WARNING) << "The rank defect of the symmetric " << nRows << " by " << nRows << " matrix is "
                         << nRows - nFixed - static_cast<unsigned int>(rank);
        }
        if(!m_iterate)
            break;
        ++iteration;
        if(iteration == nMaxIterations)
            break;
        // Update the factor in the track chi2 cut.
        const double newcfactr = sqrt(m_cfactr);
        if(newcfactr > 1.2 * m_cfactref) {
            m_cfactr = newcfactr;
        } else {
            m_cfactr = m_cfactref;
        }
        LOG(DEBUG) << "Refitting tracks with cut factor " << m_cfactr;

        // Reset global variables.
        for(unsigned int i = 0; i < nRows; ++i) {
            m_bgvec[i] = 0.0;
            for(unsigned int j = 0; j < nRows; ++j)
                m_cgmat[i][j] = 0.0;
        }
        // Refit the tracks.
        double chi2 = 0.;
        double ndof = 0.;
        nGoodTracks = 0;
        for(size_t i = 0; i < nTracks; ++i) {
            // Skip invalidated tracks.
            if(m_equations[i].empty()) {
                continue;
            }

            std::vector<Equation> equations(m_equations[i].begin(), m_equations[i].end());
            for(auto equation : equations) {
                const size_t nG = equation.derG.size();
                for(size_t j = 0; j < nG; ++j) {
                    const double t = localParams[i][static_cast<size_t>(equation.indNL[j])];
                    if(t == 0)
                        continue;
                    equation.derG[j] += equation.derNL[j] * (t - equation.slopes[j]);
                }
            }
            std::fill(trackParams.begin(), trackParams.end(), 0.);
            // Refit the track.
            bool ok = fitTrack(equations, trackParams, false, iteration);
            // Cache the track state.
            for(unsigned int j = 0; j < m_nalc; ++j) {
                localParams[i][j] = trackParams[2 * j];
            }
            if(ok) {
                // Update the total chi2.
                chi2 += trackParams[2 * m_nalc + 1];
                ndof += trackParams[2 * m_nalc];
                ++nGoodTracks;
            } else {
                // Disable the track.
                m_equations[i].clear();
            }
        }
        LOG(INFO) << "Chi2 / DOF after re-fit: " << chi2 / (ndof - static_cast<double>(nRows));
    }

    // Print the final results.
    printResults();
    return true;
}

//=============================================================================
// Obtain solution of a system of linear equations with symmetric matrix
// and the inverse (using 'singular-value friendly' GAUSS pivot).
// Solve the equation V * X = B.
// V is replaced by its inverse matrix and B by X, the solution vector
//=============================================================================
int AlignmentMillepede::invertMatrix(std::vector<std::vector<double>>& v, std::vector<double>& b, const size_t n) {
    int rank = 0;
    const double eps = 0.0000000000001;

    std::vector<double> diag(n, 0.);
    std::vector<bool> used_param(n, true);
    std::vector<bool> flag(n, true);
    for(size_t i = 0; i < n; i++) {
        for(size_t j = 0; j <= i; j++) {
            v[j][i] = v[i][j];
        }
    }

    // Find max. elements of each row and column.
    std::vector<double> r(n, 0.);
    std::vector<double> c(n, 0.);
    for(size_t i = 0; i < n; i++) {
        for(size_t j = 0; j < n; j++) {
            if(fabs(v[i][j]) >= r[i])
                r[i] = fabs(v[i][j]);
            if(fabs(v[j][i]) >= c[i])
                c[i] = fabs(v[j][i]);
        }
    }
    for(size_t i = 0; i < n; i++) {
        if(0.0 != r[i])
            r[i] = 1. / r[i];
        if(0.0 != c[i])
            c[i] = 1. / c[i];
        // Check if max elements are within requested precision.
        if(eps >= r[i])
            r[i] = 0.0;
        if(eps >= c[i])
            c[i] = 0.0;
    }

    // Equilibrate the V matrix
    for(size_t i = 0; i < n; i++) {
        for(size_t j = 0; j < n; j++) {
            v[i][j] = sqrt(r[i]) * v[i][j] * sqrt(c[j]);
        }
    }

    for(size_t i = 0; i < n; i++) {
        // Save the absolute values of the diagonal elements.
        diag[i] = fabs(v[i][i]);
        if(r[i] == 0. && c[i] == 0.) {
            // This part is empty (non-linear treatment with non constraints)
            flag[i] = false;
            used_param[i] = false;
        }
    }

    for(size_t i = 0; i < n; i++) {
        double vkk = 0.0;

        size_t k = 0;
        bool found_k = false;
        // First look for the pivot, i. e. the max unused diagonal element.
        for(size_t j = 0; j < n; j++) {
            if(flag[j] && (fabs(v[j][j]) > std::max(fabs(vkk), eps))) {
                vkk = v[j][j];
                k = j;
                found_k = true;
            }
        }

        if(found_k) {
            // Pivot found.
            rank++;
            // Use this value.
            flag[k] = false;
            // Replace pivot by its inverse.
            vkk = 1.0 / vkk;
            v[k][k] = -vkk;
            for(size_t j = 0; j < n; j++) {
                for(size_t jj = 0; jj < n; jj++) {
                    if(j != k && jj != k && used_param[j] && used_param[jj]) {
                        // Other elements (do them first as you use old v[k][j])
                        v[j][jj] = v[j][jj] - vkk * v[j][k] * v[k][jj];
                    }
                }
            }
            for(size_t j = 0; j < n; j++) {
                if(j != k && used_param[j]) {
                    v[j][k] = v[j][k] * vkk;
                    v[k][j] = v[k][j] * vkk;
                }
            }
        } else {
            // No more pivot value (clear those elements)
            for(size_t j = 0; j < n; j++) {
                if(flag[j]) {
                    b[j] = 0.0;
                    for(size_t l = 0; l < n; l++) {
                        v[j][k] = 0.0;
                        v[k][j] = 0.0;
                    }
                }
            }
            break;
        }
    }
    // Correct matrix V
    for(size_t i = 0; i < n; i++) {
        for(size_t j = 0; j < n; j++) {
            v[i][j] = sqrt(c[i]) * v[i][j] * sqrt(r[j]);
        }
    }

    std::vector<double> temp(n, 0.);
    for(size_t j = 0; j < n; j++) {
        // Reverse matrix elements
        for(size_t jj = 0; jj < n; jj++) {
            v[j][jj] = -v[j][jj];
            temp[j] += v[j][jj] * b[jj];
        }
    }

    for(size_t j = 0; j < n; j++) {
        b[j] = temp[j];
    }
    return rank;
}

//=============================================================================
// Simplified version.
//=============================================================================
int AlignmentMillepede::invertMatrixLocal(std::vector<std::vector<double>>& v, std::vector<double>& b, const size_t n) {

    int rank = 0;
    const double eps = 0.0000000000001;

    std::vector<bool> flag(n, true);
    std::vector<double> diag(n, 0.);
    for(size_t i = 0; i < n; i++) {
        // Save the absolute values of the diagonal elements.
        diag[i] = fabs(v[i][i]);
        for(size_t j = 0; j <= i; j++) {
            v[j][i] = v[i][j];
        }
    }

    for(size_t i = 0; i < n; i++) {
        double vkk = 0.0;

        size_t k = 0;
        bool found_k = false;
        // First look for the pivot, i. e. the max. unused diagonal element.
        for(size_t j = 0; j < n; j++) {
            if(flag[j] && (fabs(v[j][j]) > std::max(fabs(vkk), eps * diag[j]))) {
                vkk = v[j][j];
                k = j;
                found_k = true;
            }
        }

        if(found_k) {
            // Pivot found
            rank++;
            flag[k] = false;
            // Replace pivot by its inverse.
            vkk = 1.0 / vkk;
            v[k][k] = -vkk;
            for(size_t j = 0; j < n; j++) {
                for(size_t jj = 0; jj < n; jj++) {
                    if(j != k && jj != k) {
                        // Other elements (do them first as you use old v[k][j])
                        v[j][jj] = v[j][jj] - vkk * v[j][k] * v[k][jj];
                    }
                }
            }

            for(size_t j = 0; j < n; j++) {
                if(j != k) {
                    v[j][k] = v[j][k] * vkk;
                    v[k][j] = v[k][j] * vkk;
                }
            }
        } else {
            // No more pivot values (clear those elements).
            for(size_t j = 0; j < n; j++) {
                if(flag[j]) {
                    b[j] = 0.0;
                    for(k = 0; k < n; k++)
                        v[j][k] = 0.0;
                }
            }
            break;
        }
    }

    std::vector<double> temp(n, 0.);
    for(size_t j = 0; j < n; j++) {
        // Reverse matrix elements
        for(size_t jj = 0; jj < n; jj++) {
            v[j][jj] = -v[j][jj];
            temp[j] += v[j][jj] * b[jj];
        }
    }
    for(size_t j = 0; j < n; j++) {
        b[j] = temp[j];
    }
    return rank;
}

//=============================================================================
// Return the limit in chi^2 / nd for n sigmas stdev authorized.
// Only n=1, 2, and 3 are expected in input.
//=============================================================================
double AlignmentMillepede::chi2Limit(const int n, const int nd) const {
    constexpr double sn[3] = {0.47523, 1.690140, 2.782170};
    constexpr double table[3][30] = {{1.0000, 1.1479, 1.1753, 1.1798, 1.1775, 1.1730, 1.1680, 1.1630, 1.1581, 1.1536,
                                      1.1493, 1.1454, 1.1417, 1.1383, 1.1351, 1.1321, 1.1293, 1.1266, 1.1242, 1.1218,
                                      1.1196, 1.1175, 1.1155, 1.1136, 1.1119, 1.1101, 1.1085, 1.1070, 1.1055, 1.1040},
                                     {4.0000, 3.0900, 2.6750, 2.4290, 2.2628, 2.1415, 2.0481, 1.9736, 1.9124, 1.8610,
                                      1.8171, 1.7791, 1.7457, 1.7161, 1.6897, 1.6658, 1.6442, 1.6246, 1.6065, 1.5899,
                                      1.5745, 1.5603, 1.5470, 1.5346, 1.5230, 1.5120, 1.5017, 1.4920, 1.4829, 1.4742},
                                     {9.0000, 5.9146, 4.7184, 4.0628, 3.6410, 3.3436, 3.1209, 2.9468, 2.8063, 2.6902,
                                      2.5922, 2.5082, 2.4352, 2.3711, 2.3143, 2.2635, 2.2178, 2.1764, 2.1386, 2.1040,
                                      2.0722, 2.0428, 2.0155, 1.9901, 1.9665, 1.9443, 1.9235, 1.9040, 1.8855, 1.8681}};

    if(nd < 1)
        return 0.;
    const int m = std::max(1, std::min(n, 3));
    if(nd <= 30)
        return table[m - 1][nd - 1];
    return ((sn[m - 1] + sqrt(float(2 * nd - 3))) * (sn[m - 1] + sqrt(float(2 * nd - 3)))) / float(2 * nd - 2);
}

//=============================================================================
// Multiply general M-by-N matrix A and N-vector X
//=============================================================================
bool AlignmentMillepede::multiplyAX(const std::vector<std::vector<double>>& a,
                                    const std::vector<double>& x,
                                    std::vector<double>& y,
                                    const unsigned int n,
                                    const unsigned int m) {

    // Y = A * X, where
    //   A = general M-by-N matrix
    //   X = N vector
    //   Y = M vector
    for(unsigned int i = 0; i < m; ++i) {
        y[i] = 0.0;
        for(unsigned int j = 0; j < n; ++j) {
            y[i] += a[i][j] * x[j];
        }
    }
    return true;
}

//=============================================================================
// Multiply symmetric N-by-N matrix V from the left with general M-by-N
// matrix A and from the right with the transposed of the same general
// matrix to form a symmetric M-by-M matrix W.
//=============================================================================
bool AlignmentMillepede::multiplyAVAt(const std::vector<std::vector<double>>& v,
                                      const std::vector<std::vector<double>>& a,
                                      std::vector<std::vector<double>>& w,
                                      const unsigned int n,
                                      const unsigned int m) {

    // W = A * V * AT, where
    //   V = symmetric N-by-N matrix
    //   A = general N-by-M matrix
    //   W = symmetric M-by-M matrix
    for(unsigned int i = 0; i < m; ++i) {
        for(unsigned int j = 0; j <= i; ++j) {
            // Reset the final matrix.
            w[i][j] = 0.0;
            // Fill the upper left triangle.
            for(unsigned int k = 0; k < n; ++k) {
                for(unsigned int l = 0; l < n; ++l) {
                    w[i][j] += a[i][k] * v[k][l] * a[j][l];
                }
            }
            // Fill the rest.
            w[j][i] = w[i][j];
        }
    }

    return true;
}

//=============================================================================
// Print results
//=============================================================================
bool AlignmentMillepede::printResults() {
    const std::string line(65, '-');
    LOG(INFO) << line;
    LOG(INFO) << " Result of fit for global parameters";
    LOG(INFO) << line;
    LOG(INFO) << "   I  Difference    Last step      Error        Pull Global corr.";
    LOG(INFO) << line;
    for(unsigned int i = 0; i < m_nagb; ++i) {
        double err = sqrt(fabs(m_cgmat[i][i]));
        if(m_cgmat[i][i] < 0.0)
            err = -err;
        if(fabs(m_cgmat[i][i] * m_diag[i]) > 0) {
            // Calculate the pull.
            const double pull = m_dparm[i] / sqrt(m_psigm[i] * m_psigm[i] - m_cgmat[i][i]);
            // Calculate the global correlation coefficient
            // (correlation between the parameter and all the other variables).
            const double gcor = sqrt(fabs(1.0 - 1.0 / (m_cgmat[i][i] * m_diag[i])));
            LOG(INFO) << std::setprecision(3) << std::scientific << std::setw(3) << i << "   " << std::setw(10) << m_dparm[i]
                      << "   " << std::setw(10) << m_bgvec[i] << "   " << std::setw(8) << std::setprecision(2) << err
                      << "   " << std::setw(9) << std::setprecision(2) << pull << "   " << std::setw(9) << gcor;
        } else {
            LOG(INFO) << std::setw(3) << i << "   " << std::setw(10) << "OFF"
                      << "   " << std::setw(10) << "OFF"
                      << "   " << std::setw(8) << "OFF"
                      << "   " << std::setw(9) << "OFF"
                      << "   " << std::setw(9) << "OFF";
        }
        if((i + 1) % (m_nagb / 6) == 0)
            LOG(INFO) << line;
    }
    for(unsigned int i = 0; i < m_nagb; ++i) {
        LOG(DEBUG) << " i=" << i << "  sqrt(fabs(cgmat[i][i]))=" << sqrt(fabs(m_cgmat[i][i])) << " diag = " << m_diag[i];
    }

    return true;
}
