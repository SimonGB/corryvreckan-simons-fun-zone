# AlignmentTime
**Maintainer**: Finn Feindt (<finn.feindt@desy.de>)\
**Module Type**: *DETECTOR*\
**Detector Type**: *all*\
**Status**: Immature

### Description
If the time measurements from different detectors have an offset with respect to each other, it might not be possible to build correlated events based on these time stamps. The purpose of this module is to analyze the timestamps provided by a time reference and any given detector in order to find the time offset between them. To do so, a sorted array of time stamps for both detectors is taken, and a set of shifts is assumed. The shifts are applied to the time stamps of the investigated detector. For each shifted time stamp the closest time stamp in the array from the reference detector is found, their residual (difference) is calculated, and filled into a histogram. When the shift matches the offset, this distribution peaks around 0 and the corresponding shift can be added in the geometry description (`time_offset`).

To make this work, the spacing of the shifts needs to be smaller than the trigger frequency. If the offset is large, this might lead to a large number of shifts that need to be investigated. In this case it helps if the range of considered shifts can be constrained. It is recommended to use the `Metronome` module with an arbitrary spacing, so that the pixel time stamps of the considered detectors are added to the clipboard.

### Discuss: 
* At the moment the analyzed time stamps are those saved within the pixel objects for the given detectors. Could add option to take:
  * Cluster timestamps
  * Triggers
  * Event begin, end or the middle in between
  * Or a configurable combination between them
* To get the time stamps for the reference and investigated detector the metronome can be used with any spacing. This makes is tricky if one of the detectors is auxiliary, because we need to assign the timestamps to some sort of pixel hit for this to work.

### Parameters
* `reference_name`: Name of the detector providing the reference time stamps. Needs to be given for this module to work.
* `shift_start`: Sets the starting point of the shift scan.
* `shift_end`: Sets the end of the shift scan.
* `shift_n`: Sets the number of scanned shifts shifts. If this or the previous two parameters are not set the trigger period (inverse frequency) is estimated from the time stamps of the investigated detector. The scan is performed for 200 steps between -5 and 5 times the trigger period. **DoDo:** Optimize.
* `time_scale`: Sets the minimum and maximum of the residual histogram.
* `time_nbins`: Sets the number of bins for the residual histogram. If this or the previous parameter are not set the time scale defaults to 5 times the trigger period. The number of bins defaults to 200.

### Plots produced
* 1D histograms:
  * Pixel timestamps for the investigated detectors with two time scales, 3 s and 3000 s.
  * Pixel timestamps for the reference detector with two time scales, 3 s and 3000 s.

* 2D histograms:
  * Distribution of the time residuals between the investigated and the reference detector as a function of the applied shift.

### Usage
```toml
[AlignmentTime]
type = "dSiPM"
reference_name = "MIMOSA26_2"
```
