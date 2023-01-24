---
# SPDX-FileCopyrightText: 2017-2023 CERN and the Corryvreckan authors
# SPDX-License-Identifier: CC-BY-4.0 OR MIT
---
**Maintainer**: Fabian Lex (<fabianlex.fsl@gmail.com>)  
**Module Type**: *DETECTOR*  
**Detector Type**: *ALiBaVa*  
**Status**: Functional

### Description
This module allows data recorded by the ALiBaVa system and stored in either a ALiBaVa binary or a HDF5 file to be read into Corryvreckan as raw detector data. If, in addition to the run data file, a pedestal file is provided, the pedestal and noise will be calculated from it. At the moment it is not possible to use a calibration file to convert the arbitrary ADC counts into charge.

If the binary file format is chosen, the data and pedestal files need to include the run number in their name and end on `.dat` or `.ped` respectively. If the HDF5 format is chosen, the file name needs to include either `dat_run` or `ped_run`, the run number and end on `.hdf`.

It requires either another event loader of another detector type before, which defines the event start and end times by placing an Event definition on the clipboard, or an instance of the Metronome module which provides this information.

The detector needs to be defined in the geometry file with n columns and 1 row (everything else will lead to strange behavior later on).

The ROI can be set by creating a mask file for the detector and setting all channels outside the ROI to masked.

The code for reading and interpreting ALiBaVa data is based on analysis scripts originally written by [Alibava Systems S.L.](https://alibavasystems.com/alibava-system-classic/) and is published here with their consent. All related code is located in the `ALiBaVa` sub-folder of the module source.

### Parameters
* `input_directory`: Path to the directory where the input files can be found. This parameter is mandatory
* `run`: Number of the run to be analyzed. Default is 0.
* `timecut_low`: Readouts with a timestamp smaller than the lower limit are discarded. Each time a readout is triggered in the ALiBaVa system, the readout is given a timestamp by the ALiBaVa system in relation to its internal, 100ns period clock cycle. Default is 0ns
* `timecut_up`: Readouts with a timestamp larger than the upper limit are discarded. Each time a readout is triggered in the ALiBaVa system, the readout is given a timestamp by the ALiBaVa system in relation to its internal, 100ns period clock cycle. Default is 100ns
* `ignore_events`: Number of events at the start which will be ignored. This is done to ensure synchronization between ALiBaVa system and telescope. Default is 1.
* `calibration_constant`: Rudimentary way to allow for a conversion from ADC to kiloelectrons. Will change in the future. Default is 1.
* `chargecut`: If the charge of a strip is below the charge cut, the strip will not be added to the current event. Default is 0.
* `polarity`: Correction factor for the sign of the signal. Either +1 or -1. Depending on the type of detector (p-in-n or n-in-p) the signal measured by the ALiBaVa system is negative, the polarity corrects this in the analysis. Default is `-1` (needed for n-in-p sensors)

### Plots produced
For each detector the following plots are produced:

* Charge of signal
* ADC count of signal
* Signal to noise ratio
* Uncorrected pedestal
* Uncorrected noise
* Corrected pedestal
* Corrected noise
* 2D Corrected pedestal
* 2D Corrected noise
* Time profile

### Usage
```ini
[EventLoaderALiBaVa]
input_directory = "/path/to/my/alibava_data/"
run = 137
timecut_low = 2ns
timecut_up = 12ns
ignore_events = 1
chargecut = 42
polarity = -1
```
