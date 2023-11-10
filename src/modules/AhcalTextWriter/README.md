# AhcalTextWriter
**Maintainer**: Jiri Kvasnicka (<jiri.kvasnicka@desy.de>)  
**Module Type**: *GLOBAL*  
**Detector Type**: *all*  
**Status**: Functional

### Description
This module extends the TextWriter. It lists the incident coordinates together with the trigger ID. This information is used for the analysis outside of the corrywreckan.

The output text file will contain
* Run number
* Event number
* Trigger number
* X coordinate
* Y coordinate
* Chi2ndof

If the event cannot be used, the text output will contain a "#skipping even" line with a explenation

### Parameters
* `file_name`:  Name of the data file to create, relative to the output directory of the framework. The file extension `.txt` will be appended if not present.
* `runnr`: EUDAQ Run number - necessary to fill when multiple runs are used
* `incident_z`: Z coordinate where the hit will be extrapolated
### Usage
```toml
[AhcalTextWriter]
file_name = "exampleFileName.txt"
```

The parameters can be changed outside of the conf file in the loop in bash:
```bash
for RUNNR in `seq 68000 70000` ; do
    EUFILENAME=`find  "${EURAWDIR}" -maxdepth 1 -iregex ".*online.Run0*${RUNNR}.*raw"`
    if [ ! -f "${EUFILENAME}" ] ; then continue ; fi
    echo "Processing $RUNNR $EUFILENAME"
    ~/git/corryvreckan/bin/corry -c batch.conf \
				 -o EventLoaderEUDAQ2.file_name=${EUFILENAME} \
				 -o AhcalTextWriter.runnr=${RUNNR} \
				 -o AhcalTextWriter.file_name="${RUNNR}_intercept" \
				 -o AhcalTextWriter.incident_z=${INTERCEPTZ} \
				 -o Tracking4D.spatial_cut_rel=100
```
