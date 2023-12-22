---
# SPDX-FileCopyrightText: 2017-2023 CERN and the Corryvreckan authors
# SPDX-License-Identifier: CC-BY-4.0 OR MIT
---
# ClusteringSpatial
**Maintainer**: Daniel Hynds (<daniel.hynds@cern.ch>)  
**Module Type**: *DETECTOR*  
**Detector Type**: *all*  
**Status**: Functioning

### Description
This module clusters the input data of a detector without individual hit timestamps.
The clustering method only uses positional information: either charge-weighted center-of-gravity or arithmetic mean calculation, and no timing information.
If the pixel information is binary (i.e. no valid charge-equivalent information is available), the arithmetic mean is calculated for the position.
Also, if one pixel of a cluster has charge zero, the arithmetic mean is calculated even if charge-weighting is selected because it is assumed that the zero-reading is false and does not to represent a low charge but an unknown value.
These clusters are stored on the clipboard for each device.
Split clusters can be recovered using a larger search radius for neighboring pixels.
Their width is defined as the maximum extent in column/row direction, i.e. a cluster of pixels (1,10), (1,12) would have a column width of 1 and a row width of 3.

### Parameters
* `use_trigger_timestamp`: If true, the first trigger timestamp of the Corryvreckan event is set as the cluster timestamp. Caution when using this method for very long events containing multiple triggers. If false, the last pixel added to the cluster defines the timestamp. Default value is `false`.
* `charge_weighting`: If true, calculate a charge-weighted mean for the cluster center. If false, calculate the simple arithmetic mean. Defaults to `true`.
* `reject_by_roi`: If true, clusters positioned outside the ROI set for the detector will be rejected. Defaults to `false`.
* `neighbor_radius_col`: Search radius for neighboring pixels in column direction, defaults to `1` (do not allow split clusters)
* `neighbor_radius_row`:  Search radius for neighboring pixels in row direction, defaults to `1` (do not allow split clusters)

### Plots produced
For each detector the following plots are produced:

* Histograms for cluster size, seed charge, width (columns/X and rows/Y)
* Cluster charge histogram
* 2D cluster positions in global and local coordinates
* Cluster times
* Cluster multiplicity
* Cluster uncertainty

### Usage
```toml
[ClusteringSpatial]
use_trigger_timestamp = true
```
