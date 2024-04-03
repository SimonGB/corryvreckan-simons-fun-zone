---
# SPDX-FileCopyrightText: 2017-2023 CERN and the Corryvreckan authors
# SPDX-License-Identifier: CC-BY-4.0 OR MIT
---
# AlignmentDUTResidual
**Maintainer**: Daniel Hynds (<daniel.hynds@cern.ch>), Simon Spannagel (<simon.spannagel@cern.ch>)  
**Module Type**: *DUT*  
**Detector Type**: *all*  
**Status**: Functional

### Description
This module performs translational and rotational DUT alignment. The alignment is performed with respect to the reference plane set in the configuration file.

This module uses tracks for alignment. The module moves the detector it is instantiated for and minimizes the unbiased residuals calculated from the track intercepts with the plane.

### Parameters
* `iterations`: Number of times the chosen alignment method is to be iterated. Default value is `3`.
* `align_position`: Boolean to select whether to align the translational displacements of the detector or not. Note that the Z displacement is never aligned. Specify the axes using `align_position_axes`. The default value is `true`.
* `align_orientation`: Boolean to select whether to align the rotations of the detector under consideration or not. Specify the axes using `align_orientation_axes`. The default value is `true`.
* `align_position_axes`: Define for which axes to perform translational alignment. The default value is `xy`, which means both X and Y displacements of the detector will be aligned.
* `align_orientation_axes`: Define for which axes to perform rotational alignment if `align_orientation = true`. The default value is `012`, which means that rotations around all axis will be aligned. `012` always represent the angles in order of the geometry file, independent of the `orientation_mode`
* `prune_tracks`: Boolean to set if tracks with a number of associated clusters > `max_associated_clusters` or with a track chi^2 > `max_track_chi2ndof` should be excluded from use in the alignment. The number of discarded tracks is written to the terminal. Default is `false`.
* `max_associated_clusters`: Maximum number of associated clusters per track allowed when `prune_tracks = true` for the track to be used in the alignment. Default value is `1`.
* `max_track_chi2ndof`: Maximum track chi^2 value allowed when `prune_tracks = true` for the track to be used in the alignment. Default value is `10.0`.
* `workers`: Specify the number of workers to use in total, should be strictly larger than zero. Defaults to the number of native threads available on the system minus one, if this can be determined, otherwise one thread is used.
* `residuals`: Array of formulas for unbiased x and y residuals. Any 2D TFormula can be used, the variables `x` and `y` represent the *track intercept* with the plane and the *cluster position* respectively. Default formulas: `x - y`. Parameters can be used (`[0]`, `[1]`, ...) and have to be separately specified (see below). It should be noted that the formula does *not* support units, values with units have to be specified as separate parameters. Both and only the functions residual_x and residual_y must be defined if used. 
* `parameters_residuals`: Array of factors, representing the parameters of the above correction function. Defaults to an empty array, i.e. by default no parameters are needed.
* `spatial_cut_sensoredge` : Define the minimal distance a track has to have from the sensors edge. Defaults to `0`

### Plots produced
For the DUT, the following plots are produced:

* Residuals in X and Y (calculated in local coordinates)
* Profile plot of residual in X vs. X, X vs. Y, Y vs. X and Y vs. Y position
* Graphs of the translational shifts along the X/Y-axis vs. the iteration number
* Graphs of the rotational shifts along the X/Y/Z-axis vs. the iteration number

### Usage
```toml
[Corryvreckan]
# The global track limit can be used to reduce the run time:
number_of_tracks = 200000

[AlignmentDUTResidual]
# example of redefinition of residuals in x and y
residuals = "(TMath::Abs(x - y) - [0])*TMath::Sign(1,(x - y))","(TMath::Abs(x - y) - [0])*TMath::Sign(1,(x - y))*[1]"
parameters_residuals = 16.9um, 17.4um, 1
```
