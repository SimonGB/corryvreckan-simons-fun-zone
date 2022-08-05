# ClusteringAnalog
**Maintainer**: Miljenko Suljic (<miljenko.suljic@cern.ch>), Yitao WU (<yitao.wu@cern.ch>)  
**Module Type**: *DETECTOR*  
**Detector Type**: *all*  
**Status**: Functioning

### Description
This module clusters the input data of an analog detector.
The clustering method employs thresholds (fixed, signal-to-noise, or combination) to select the seed pixels and their neighbors. This procedure is performed iteratively for neighbors to include in the cluster all pixels with the signal above a given threshold. Alternatively, neighbors can be defined as pixels in an `NxN` matrix centered on the seed pixel.
For analysis of signal-noise ratio, this module reads the pedestal and noise map for each pixel in the input file, which is prepared from external and specified in detector configuration as `calibration_file`. Current implementation allows only ROOT files containing `TH2F` maps for pedestal and noise.

This module also provides a cluster shape analysis. To characterize the cluster shape, the charge distribution is estimated in the entire clustering window and the neighboring pixels are ordered by 1D index number w.r.t. the local position of the seed or by decreasing charge. To understand the signal significance and consider the common shift effect, a special plot for charge ratio is accumulated with the largest N pixels in the cluster, using the sum of all positive pixels as the denominator to normalize cluster-by-cluster.

### Parameters
* `reject_by_roi`: ROI rejection with the local position of the cluster. (Default: false)
* `method`: Clustering method to reconstruct cluster position and charge. `cluster` (default): includes all adjecent pixels with signal above thresholds (see later) and calculates charge-weighted center-of-gravity as the cluster position. `seed`: equivalent to `cluster` but the cluster position and charge is given only by the seed pixel. `binary`: equivalent to `cluster` but calculates center-of-gravity of all pixels above threshold without charge weighting (as is done in binary sensors). `sumNxN`: includes all pixels in the `NxN` matrix around the seed, where `N=window_size` and calculates charge-weighted center-of-gravity as the cluster position.
* `seeding_method`: Method to select seed pixels. Option `multi`  (default) selects all pixels above the seed threshold for clustering and allows multiple clusters/hits in the same event. Option `max` keeps only the single seed with the maximum charge (if above seed threshold) as seed candidate.
* `window_size`: Matrix width around the seed to find neighbors in `sumNxN` method. It also used to define central pixels in histograms (see later). The module is only capable to deal with an odd number larger than 3. (Default: 3)
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

### Plots produced
For each detector, the following plots are produced:

* Histograms for the charge, S/N ratio, charge ratio of seed, and neighbors.
* Cluster size and charge distribution.
* Correlation map for seed and neighbors, seed and cluster.
* 2D cluster positions in global and local coordinates.
* Histograms for size and charge of central clusters, of which seed is not on the edge of the sensor as defined by `window_size`.
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
window_size=3
threshold_type=snr                 # Cut on SNR
thresholdSNR_seed=3
thresholdSNR_neighbor=2
method=cluster
seeding_method=multi
calibration_pedestal=hPedestalpl1  # Histogram name
calibration_noise=hnoisepl1        # Histogram name
```
