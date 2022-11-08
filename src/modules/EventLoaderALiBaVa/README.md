**Maintainer**: Fabian Lex (“fabianlex.fsl@gmail.com”)
**Module Type**: *DETECTOR* **Detector Type**: *ALiBaVa*
**Status**: Immature

### Description
This module allows data recorded by the ALiBaVa system and stored in either a ALiBaVa binary or a HDF5 file to be read into Corryvreckan as raw detector data. If, in addition to the run data file, a pedestal file is provided, the pedestal and noise will be calculated from it. At the moment it is not possible to use a calibration file to convert the arbitrary ADC counts into charge. 

If the binary file format is chosen, the data and pedestalfile need to include the run number in their name and end on ".dat" or ".ped" respectively. If the HDF5 format is chosen, the file name needs to include either "dat_run" or "ped_run", the run number and end on ".hdf". 

It requires either another event loader of another detector type before, which defines the event start and end times by placing an Event definition on the clipboard, or an instance of the Metronome module which provides this information.

The detector needs to be defined in the geometry file with n columns and 1 row (everything else will lead to strange behaviour later on). 

The ROI can be set by creating a maskfile for the detector and setting all channels outside the ROI to masked. 

### Parameters
* `input_directory`: Path to the directory where the input files can be found. This parameter is mandatory
* `run`: Number of the run to be analysed. Default is 0.
* `timecut_low`: Readouts with a timestmap smaller than the lower limit are discarded. Each time a readout is triggered in the ALiBaVa system, the readout is given a timestamp by the ALiBaVa system in relation to its internal, 100 ns period clock cycle. Default is 0 ns
* `timecut_up`: Readouts with a timestmap larger than the upper limit are discarded. Each time a readout is triggered in the ALiBaVa system, the readout is given a timestamp by the ALiBaVa system in relation to its internal, 100 ns period clock cycle. Default is 100 ns
* `ignore_events`: Number of events at the start which will be ignored. This is done to ensure synchronisation between ALiBaVa system and telescope. Default is 1.
* `calibration_constant`: Rudimentary way to allow for a conversion from ADC to kiloelectrons. Will change in the future. Default is 1.
* `chargecut`: If the charge of a strip is below the chargecut, the strip will not be added to the current event. Default is 0. 
* `polarity`: Correction factor for the sign of the signal. Either +1 or -1. Depending on the type of detector (p-in-n or n-in-p) the signal measured by the ALiBaVa system is negative, the polarity corrects this in the analysis. Default is -1 (needed for n-in-p sensors)

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
```toml
[EventLoaderALiBaVa]
type = "ALiBaVa"
input_directory = "/home/fabian/testbeam_analysis/testdirectory/"
run = 137
timecut_low = 2ns
timecut_up = 12ns
ignore_events = 1
chargecut = 42
polarity = -1
```



