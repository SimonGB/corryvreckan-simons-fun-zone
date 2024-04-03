---
# SPDX-FileCopyrightText: 2017-2023 CERN and the Corryvreckan authors
# SPDX-License-Identifier: CC-BY-4.0 OR MIT
---
# AnalysisMaterialBudget
**Maintainer**: Paul Schuetze (paul.schuetze@desy.de)
**Module Type**: *GLOBAL*  
**Status**: Functional

### Description
This module analyses the kink angle distributions at the position of a scatterer. This requires the availability of `Multiplet` tracks.
A material budget image is computed, which represents the widths of the scattering angle distribution of all particles traversing a given image cell.
The widths are calculated as average absolute deviation from `0` for the selected quantile of the distribution.
At this point, no attempt is made to translate this width into the actual material budget (x/X0), given the fact that the absolute values strongly depend on the detector setup and the detector resolution and include scattering in air and neighbouring detectors. For this, a calibration would be required for every individual setup.

Further information on this technique can be found in [@material_budget_imaging].



### Parameters
* `image_size`: Two dimensional Field of view extent for all histograms. Defaults to `10mm 10 mm`.
* `cell_size`: Binning of histograms and image cell sizes for material budget images in two dimensions. Defaults to `50um 50um`.
* `angle_cut`: Maximum kink angle to evaluate. Defaults to `100 mrad`.
* `quantile`: Fraction of entries per distribution for which the width is evaluated in material budget images. Defaults to `0.9`.
* `min_cell_content`: Minimum number of registered kink angles per image cell required for the cell evaluation in the material budget image. Defaults to `20`.
* `update`: Determines whether the material budget image is updated during run time. Otherwise this process is done only once during the finalisation. Defaults to `false`.

### Plots produced
* Histogram of kink angles in x/y
* Profile of squared kink angles along x/y
* 2D profile of squared kink angles (material budget image)
* 2D histogram of angle distribution width per image cell (material budget image)
* Histogram of the number of registered kink angles per image cell

### Usage
```toml
[AnalysisMaterialBudget]
image_size = 20mm 10mm
cell_size = 100um 100um

```

[@material_budget_imaging]: https://doi.org/10.1063/1.5005503