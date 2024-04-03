---
# SPDX-FileCopyrightText: 2017-2023 CERN and the Corryvreckan authors
# SPDX-License-Identifier: CC-BY-4.0 OR MIT
---
# EtaCorrection
**Maintainer**: Daniel Hynds (<daniel.hynds@cern.ch>), Simon Spannagel (<simon.spannagel@cern.ch>), Morag Williams (<morag.williams@cern.ch>)  
**Module Type**: *DETECTOR*  
**Detector Type**: *all*  
**Status**: Functional  

### Description
This module applies previously determined $`\eta`$-corrections to cluster positions of any detector in order to correct for non-linear charge sharing. Corrections can be applied to any cluster read from the clipboard. The correction function as well as the parameters for each of the detectors can be given separately for X and Y via the configuration file.

This module does not calculate the $`\eta`$ distribution.

To use different eta correction factors for different detectors, named module instances should be used as in the example below.

### Parameters
* `eta_formula_x` / `eta_formula_y`: The formula for the $`\eta`$ correction to be applied for the X an Y coordinate, respectively. It defaults to a polynomial of fifth order, i.e. `[0] + [1]*x + [2]*x^2 + [3]*x^3 + [4]*x^4 + [5]*x^5`.
* `eta_constants_x` / `eta_constants_y`: Vector of correction factors, representing the parameters of the above correction function, in X and Y coordinates, respectively. Defaults to an empty vector, i.e. by default no correction is applied.

### Plots produced
For each detector the following plots are produced:

* Profile plot of the applied $`\eta`$-correction, for X and Y respectively

### Usage
```toml
[EtaCorrection]
name = "dut"
eta_formula_x = [0] + [1]*x + [2]*x^2 + [3]*x^3 + [4]*x^4 + [5]*x^5
eta_constants_x = 0.025 0.038 6.71 -323.58 5950.3 -34437.5
```
