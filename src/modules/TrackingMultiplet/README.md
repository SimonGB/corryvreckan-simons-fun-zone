---
# SPDX-FileCopyrightText: 2017-2023 CERN and the Corryvreckan authors
# SPDX-License-Identifier: CC-BY-4.0 OR MIT
---
# TrackingMultiplet
**Maintainer**: Paul Schuetze (paul.schuetze@desy.de)
**Module Type**: *GLOBAL*  
**Status**: Functional

### Description
This module performs particle tracking based on the _Multiplet_ track model.
This track model is defined by two so-called tracklets (_upstream_ & _downstream_), which are track fits to hits in the corresponding subsets of detectors.
These tracklets are connected at a certain position along `z`, where an arbitrary kink of the track is allowed, representing a scatterer.

The upstream and downstream tracklets are determined by adding all combinations of clusters in the first and the last detector plane of the corresponding stream to track candidates.
For each plane the cluster closest to the tracklet is added, while the fit to this tracklet is continuously updated.
Time cuts are taken into account.
A tracklet is stored, when it has at least `min_hits_<up/down>stream` clusters and is not closer than `isolation_cut` to another tracklet at the position of the scatterer.

For finding `Multiplet`s the upstream and downstream tracklets are matched as follows:
For each upstream tracklet, the downstream tracklet with the lowest matching distance is chosen, where the matching distance is determined via an extrapolation of both arms to the position of the scatterer.

### Parameters
* `scatterer_position`: Position of the scatterer along `z`. Defaults to the position of the DUT, if exactly one DUT is present.
* `upstream_detectors`: Names of detectors associated to the upstream tracklet. Defaults to all detectors in front of the scatterer, except for `DUT` and `Auxiliary` detectors.
* `downstream_detectors`: Names of detectors associated to the downstream tracklet. Defaults to all detectors behind the scatterer, except for `DUT` and `Auxiliary` detectors.
* `min_hits_upstream`: Minimum number of associated clusters required to create an upstream tracklet. Default value is the number of upstream detectors.
* `min_hits_downstream`: Minimum number of associated clusters required to create an downstream tracklet. Default value is the number of downstream detectors.
* `scatterer_matching_cut`: Maximum allowed distance between the extrapolated positions of an up- and a downstream tracklet at the position of the scatterer. No default value.
* `isolation_cut`: Minimum distance for two same-side tracklets at the position of the scatterer. If closer, the two tracklets are removed. With a value of 0, the tracklets social distancing is switched off. Defaults to $2*`scatterer_matching_cut`$.
* `time_cut_rel`: Factor by which the `time_resolution` of each detector plane will be multiplied. This calculated value is then used as the maximum time difference allowed between clusters and an upstream or downstream tracklet for association to the tracklet. This allows the time cuts between different planes to be detector appropriate. By default, a relative time cut is applied. Absolute and relative time cuts are mutually exclusive. Defaults to `3.0`.
* `time_cut_abs`: Specifies an absolute value for the maximum time difference allowed between clusters and an upstream or downstream tracklet for association to the tracklet. Absolute and relative time cuts are mutually exclusive. No default value.
* `spatial_cut_rel`: Factor by which the `spatial_resolution` in x and y of each detector plane will be multiplied. These calculated value are defining an ellipse which is then used as the maximum distance in the XY plane allowed between clusters and an upstream or downstream tracklet for association to the tracklet. This allows the spatial cuts between different planes to be detector appropriate. By default, a relative spatial cut is applied. Absolute and relative spatial cuts are mutually exclusive. Defaults to `3.0`.
* `spatial_cut_abs`: Specifies a set of absolute value (x and y) which defines an ellipse for the maximum spatial distance in the XY plane between clusters and an upstream or downstream tracklet for association to the tracklet. Absolute and relative spatial cuts are mutually exclusive. No default value.
* `exclude_from_seed`: Names of detectors which should not be used as the first detector for the creation of the straigh line track seed, even if they are used in the tracking. This can be useful to build better seed tracks e.g. if the first plane in z has a much worse spatial resolution than the following tracking planes. Default is empty.
* `require_detectors`: Names of detectors which are required to have a cluster on the track. If a track does not have a cluster from all detectors listed here, it is rejected. If empty, no detector is required. Default is empty.
* `timestamp_from`: Defines the detector which provides the track timestamp. This detector is by default added to `required_detector`. If empty, the average timestamp of upstream and downstream tracklet will be used. Empty by default.
* `track_model`: Specifies the track model used for the up and downstream arms. Defaults to `straightline`
* `momentum`: Defines the beam momentum. Only required if `track_model="gbl"`
* `lorentz_beta`: Lorentz beta factor. For low momentum particles, it impacts the estimation of the scattering angle. Defaults to `1`.
* `particle_charge`: Particle charge number. Defaults to `1`.
* `refit_gbl`: Refit the multiplet tracks with GBL. Defaults to false.
* `unique_cluster_usage`: Only use a cluster for one track - in the case of multiple assignments, the track with the best chi2/ndof is kept. Defaults to `false`

### Plots produced

For both upstream and downstream tracklets (multiplet candidates and accepted tracks), the following plots are produced:

* Histogram of the tracklets per event
* Histogram of the number of clusters per tracklet
* Histograms of the tracklet angle with respect to the X/Y-axis
* Histograms of the position of the tracklet in X/Y, extrapolated to the position of the scatterer

The following plots are produced only once:

* Histograms of the matching distance of multiplet candidates in X/Y
* Histograms of the matching distance of accepted multiplets in X/Y
* Histograms of kink angles of accepted multiplets in X/Y
* Histogram of the number of accepted multiplets per event
* Histogram of the chi2(ndf) of the accepted tracks

For each detector the following plots are produced both for multiplet candidates and accepted tracks:

* Histograms of the global and local track residual in X/Y.

### Usage
```toml
[TrackingMultiplet]
spatial_cut_abs = 100um 100um
upstream_detectors = "telescope0" "telescope1" "telescope2"
downstream_detectors = "telescope3" "telescope4" "telescope5"

scatterer_position = 150mm
scatterer_matching_cut = 50um

```
