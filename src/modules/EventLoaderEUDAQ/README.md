---
# SPDX-FileCopyrightText: 2017-2023 CERN and the Corryvreckan authors
# SPDX-License-Identifier: CC-BY-4.0 OR MIT
---
# EventLoaderEUDAQ
**Maintainer**: Simon Spannagel (<simon.spannagel@cern.ch>)  
**Module Type**: *GLOBAL*  
**Status**: Functional  

### Description
This module allows data recorded by EUDAQ and stored in the EUDAQ-native raw format to be read into Corryvreckan. The EUDAQ decoder plugins are used to transform the data into the `StandardPlane` event type before storing the individual `Pixel` objects on the Corryvreckan clipboard.

The detector IDs are taken from the plane name and IDs, two possible naming options for Corryvreckan are available: When setting `long_detector_id = true`, the name of the sensor plane and the ID are used in the form `<name>_<ID>`, while only the ID is used otherwise as `plane<ID>`. Only detectors listed in the Corryvreckan geometry are decoded and stored, data from other detectors available in the same EUDAQ event are ignored.

Sometimes the user would like to give detectors different names in Corryvreckan than provided by EUDAQ, or EUDAQ stores plane IDs in the data file which are incompatible with Corryvreckan configuration section headers. In this case, the `remap_detectors` parameter can be used to assign new names to the individual EUDAQ plane IDs. The parameter is a 2D matrix, with two entries per row, the first being the EUDAQ plane ID and the second the desired Corryvreckan name of the detector. For example:

```toml
remap_detectors = [["USBPIX_GEN3_BOARD_-1_22", "FEI4"]]
```

would rename a EUDAQ plane called `USBPIX_GEN3_BOARD_-1_22`, incompatible with the Corryvreckan configuration, to a simple `FEI4`.

This module is able to process multiple runs sequentially if all the devices were integrated with EUDAQ in the data-taking. This is an important limitation preventing to add data from other devices via different event loaders, as the event building algorithm of Corryvreckan will add the data into the Event if it matches either the timestamps or trigger IDs. However, EUDAQv1.x support stopped some years ago, and this module is a bit outdated: right now is defining the Corryvreckan Event using hardcoded timestamps in the [code](https://gitlab.cern.ch/corryvreckan/corryvreckan/-/blob/master/src/modules/EventLoaderEUDAQ/EventLoaderEUDAQ.cpp#L140-143) and therefore, it is impossible to include any device from a different event loader.

### Requirements
This module requires an installation of [EUDAQ 1.x](https://github.com/eudaq/eudaq). The installation path should be set as environment variable via
```bash
export EUDAQPATH=/path/to/eudaq
```
for CMake to find the library link against and headers to include.

### Parameters
* `file_name`: File name(s) of the EUDAQ raw data file(s). At least one file is mandatory.
* `long_detector_id`: Boolean switch to configure using the long or short detector ID in Corryvreckan, defaults to `true`.
* `remap_detectors`: Matrix to map EUDAQ plane IDs to Corryvreckan detector names.

### Usage
```toml
[EventLoaderEUDAQ]
file_names = "rawdata/eudaq/run020808.raw"
long_detector_id = true
```
If all the devices were integrated in EUDAQ: you can process multiple files at once:
```toml
[EventLoaderEUDAQ]
file_names = "rawdata/eudaq/run020808.raw","rawdata/eudaq/run020810.raw"
long_detector_id = true
```

### Plots produced
For each detector the following plots are produced:

* 2D hitmaps on pixel-level
