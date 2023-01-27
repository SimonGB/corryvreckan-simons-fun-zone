---
# SPDX-FileCopyrightText: 2017-2023 CERN and the Corryvreckan authors
# SPDX-License-Identifier: CC-BY-4.0 OR MIT
---
# ClusteringAnalog
**Maintainer**: Miljenko Suljic (<miljenko.suljic@cern.ch>), Yitao WU (<yitao.wu@cern.ch>)  
**Module Type**: *DETECTOR*  
**Detector Type**: *all*  
**Status**: Functioning

### Description
This module clusters the input data of an analog detector.
The clustering method employs thresholds (fixed, signal-to-noise, or combination) to select the seed pixels and their neighbors. This procedure is performed iteratively for neighbors to include in the cluster all pixels with the signal above a given threshold. Alternatively, neighbors can be defined as all pixels in a window centered on the seed pixel.
For analysis of signal-noise ratio, this module reads the pedestal and noise map for each pixel in the input file, which is prepared from external and specified in detector configuration as `calibration_file`. Current implementation allows only ROOT files containing `TH2F` maps for pedestal and noise.

This module also provides a cluster shape analysis. To characterize the cluster shape, the charge distribution is estimated in the entire clustering window and the neighboring pixels are ordered by 1D index number w.r.t. the local position of the seed or by decreasing charge. To understand the signal significance and consider the common shift effect, a special plot for charge ratio is accumulated with the largest N pixels in the cluster, using the sum of all positive pixels as the denominator to normalize cluster-by-cluster.

### Parameters
* `reject_by_roi`: ROI rejection with the local position of the cluster. (Default: false)
* `method`: Clustering method to reconstruct cluster position and charge. `cluster` (default): includes all adjecent pixels with signal above thresholds (see later) and calculates charge-weighted center-of-gravity as the cluster position. `seed`: equivalent to `cluster` but the cluster position and charge is given only by the seed pixel. `binary`: equivalent to `cluster` but calculates center-of-gravity of all pixels above threshold without charge weighting (as is done in binary sensors). `window`: includes all pixels in a window centered on the seed and defined by `window_size` and calculates charge-weighted center-of-gravity as the cluster position.
* `seeding_method`: Method to select seed pixels. Option `multi`  (default) selects all pixels above the seed threshold for clustering and allows multiple clusters/hits in the same event. Option `max` keeps only the single seed with the maximum charge (if above seed threshold) as seed candidate.
* `include_corners`: Consider also corner touching pixels as neighbours. (Default: `false`)
* `window_size`: Size, in pixels, of the window that defines the neighbors around the seed in `window` method. It also used to define central pixels in histograms (see later). E.g. `window_size=1` is seed plus first neighbors (aka "crown") i.e. a `3x3` matrix in case of rectanguar pixels and `window_size=2` is seed plus 2 crowns i.e. a `5x5` matrix. (Default: 1)
* `threshold_type`: `fix`, `snr`, or `mix`, use fixed threshold, signal-to-noise-ratio, or both of them to clusterize. (Default: `fix`)
* `threshold_seed`: Cut on pixel charge to find seeds. Used if `threshold_type=fix or mix`.
* `threshold_neighbor`: Cut on pixel charge to find neighbors. Used if `threshold_type=fix or mix`. (Default: `threshold_seed`)
* `threshold_iteration`: Cut on pixel charge to find neighbors of a neighbor. Used if `threshold_type=fix or mix`. (Default: `threshold_neighbor`)
* `thresholdSNR_seed`: S/N ratio cut for seed pixels. Detector `calibration_file` must be defined in geometry file and `threshold_type=snr or mix`.
* `thresholdSNR_neighbor`: S/N ratio cut to find neighbors. Detector `calibration_file` must be defined in geometry file and `threshold_type=snr or mix`.
* `thresholdSNR_iteration`: S/N ratio cut to find neighbors of neighbors. Detector `calibration_file` must be defined in geometry file and `threshold_type=snr or mix`.
* `threshold_cluster`: Cut on cluster charge, used for the optimization of seeding criteria. (Default: `threshold_seed`)
* `calibration_pedestal`: Histogram name of pedestal map in calibration file. Read as ROOT::TH2F.
* `calibration_noise`: Histogram name of noise map in calibration file. Read as ROOT::TH2F.
* `analysis_shape`: Produce more elaborate histograms for cluster shape analysis. (Default: `false`)
* `use_trigger_timestamp`: If true, the first trigger timestamp of the Corryvreckan event is set as the cluster timestamp. Caution when using this method for very long events containing multiple triggers. If false, the seed pixel defines the timestamp. Default value is `false`.

### Plots produced
For each detector, the following plots are produced:

* Histograms for the charge, S/N ratio, charge ratio of seed, and neighbors.
* Cluster size and charge distribution.
* Correlation map for seed and neighbors, seed and cluster.
* 2D cluster positions in global and local coordinates.
* Histograms for size and charge of central clusters i.e. those in which the seed is at least `window_size` away from the matrix edge.
* [Optional] Cluster shape analysis for charge sharing, including 2D histograms on the charge, charge ratio, and SNR with respect.

### Usage
The simplest configuration:

```toml
[ClusteringAnalog]
reject_by_roi=true
threshold_seed=150      # threshold_neighbor use default value as seed
```

Signal-to-noise ratio threshold configuration:

```toml
[ClusteringAnalog]
reject_by_roi=true
window_size=2
threshold_type=snr                 # Cut on SNR
thresholdSNR_seed=3
thresholdSNR_neighbor=2
method=cluster
seeding_method=multi
calibration_pedestal=hPedestalpl1  # Histogram name
calibration_noise=hnoisepl1        # Histogram name
```
