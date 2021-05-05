I set everything up on the NAF using these commands, after cloning the eudaq repo and this corryvreckan repo so that you have one folder with a subfolder for corryvreckan and one for eudaq:
```
source corryvreckan/etc/setup_lxplus.sh
export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh
lsetup cmake
lsetup hdf5
lsetup "root 6.20.02-x86_64-centos7-gcc8-opt"

cd eudaq
rm -r build
mkdir build && cd build
cmake .. -DEUDAQ_BUILD_EXECUTABLE=OFF -DEUDAQ_BUILD_GUI=OFF -DUSER_TLU_BUILD=ON 
vim ../cmake/Platform.cmake //Make the following edit: set(CMAKE_CXX_STANDARD 14)  --> set(CMAKE_CXX_STANDARD 17)
cmake ..
make install -j8

cd ../../corryvreckan
rm -r build
mkdir build && cd build
cmake .. -DBUILD_EventLoaderEUDAQ2=ON -Deudaq_DIR=**YOURPATH**/eudaq/cmake
export BUILD_EventLoaderEUDAQ2=ON
export eudaq_DIR=**YOURPATH**/eudaq/cmake
cmake ..
make install -j8

export HDF5_USE_FILE_LOCKING=FALSE
```


# EventLoaderALiBaVa
**Maintainer**: jclercx (“jclercx”)
**Module Type**: *DETECTOR* **Detector Type**: *<add types here>*  
**Status**: Immature

### Description
This is a demonstrator module only, taking data every detector on the clipboard and plots the pixel hit positions.
It serves as template to create new modules.

### Parameters
No parameters are used from the configuration file.

### Plots produced
* Histogram of event numbers

For each detector the following plots are produced:

* 2D histogram of pixel hit positions

### Usage
```toml
[EventLoaderALiBaVa]

```
