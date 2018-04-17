## Timepix1EventLoader
**Maintainer**: Daniel Hynds (<daniel.hynds@cern.ch>)  
**Status**: Functional  

#### Description
This module loads raw data from Timepix1 devices and adds it to the clipboard as `pixel` objects. The input file must have extension `.txt`, and these files are sorted into time order using the file titles.

#### Parameters
* `inputDirectory`: Path of the directory above the data files.

#### Plots produced
No plots are produced.

#### Usage
```toml
[Timepix1EventLoader]
inputDirectory = "path/to/directory"
```
Parameters to be used in multiple modules can also be defined globally at the top of the configuration file. This is highly encouraged for parameters such as `DUT` and `reference`.