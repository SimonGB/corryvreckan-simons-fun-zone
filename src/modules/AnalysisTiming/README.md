---
# SPDX-FileCopyrightText: 2023 CERN and the Corryvreckan authors
# SPDX-License-Identifier: CC-BY-4.0 OR MIT
---
# AnalysisTiming
**Maintainer**: Stephan Lachnit (stephan.lachnit@cern.ch)
**Module Type**: *DETECTOR*
**Detector Type**: *all*
**Status**: Functional

### Description
This module allows to study the time resolution of any detector with respect to any other detector.
This can be used to study the time resolution using the time residual between two DUTs.

### Parameters
* `reference_name`: Name of the detector to use as reference. Defaults to the reference detector defined in the geometry file.
* `reference_associated_clusters`: If true, uses associated cluster for tracks. Associated clusters need to be created e.g. with the `DUTAssociation` module. If false, then cluster is taken from tracking, e.g. if the detector is in `Tracking4D`. Defaults to true, except if `reference_name` is omitted.
* `dut_associated_clusters`: Identical to `reference_associated_clusters` setting, except used for the DUT. Defaults to true.
* `chi2ndof_cut`: Acceptance criterion for the maximum telescope track Chi2/ndof, defaults to `3`.
* `time_range`: Total time window for the time residual histogram.
* `time_binning`: Bin size for the time residual histogram.
* `time_offset`: Time offset for the time residual histogram with respect to the reference detector. Defaults to `0ns`.
* `inpixel_bin_size`: The bin size for inpixel plots. Different bin sizes can be set for the x and y axis. Defaults to `1um, 1um`.

### Plots produced
For each detector the following plots are produced:

* Time residual histogram
* Time residual over time (1h)
* 2D sensor map for the mean and standard deviation of the time residual
* 2D in-pixel map for the mean and standard deviation of the time residual

### Usage
```toml
[AnalysisTiming]
name = "dSiPM_0"
reference_detector = "dSiPM_1"
reference_associated_clusters = true
time_range = 600ns
time_binning = 95ps
```
