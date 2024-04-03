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
This module allows to study the time resolution of any detector with respect to a reference, which can be any other detector or the track timestamp.
Neither the "DUT" nor the reference require to have a timestamp for every track, if one timestamp is missing the track is skipped.
This flexibility can be used for example to study the time resolution using the time residual between two DUTs, or between two telescope planes.

### Parameters
* `reference_type`: Either `dut`, `plane` or `track`. Defaults to `dut`.
  If the value is `dut`, then the reference timestamp is extracted from the closest associated cluster of the reference detector. Associated clusters need to be created e.g. with the `DUTAssociation` module.
  If the value is `plane`, then the reference timestamp is extracted from the cluster of the reference detector used in the track, e.g. if the detector provided a cluster in `Tracking4D`.
  If the value is `track`, then the reference timestamp is taken from the track timestamp. Track timestamp with a value of `0` are filtered out.
* `reference_name`: Name of the detector to use as reference if `reference_type` is not `track`.
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
reference_name = "dSiPM_1"
reference_type = "dut"
time_range = 600ns
time_binning = 95ps
```
