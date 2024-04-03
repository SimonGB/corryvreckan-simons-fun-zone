---
# SPDX-FileCopyrightText: 2017-2023 CERN and the Corryvreckan authors
# SPDX-License-Identifier: CC-BY-4.0 OR MIT
---
# EtaCalculation
**Maintainer**: Daniel Hynds (<daniel.hynds@cern.ch>), Simon Spannagel (<simon.spannagel@cern.ch>), Morag Williams (<morag.williams@cern.ch>)  
**Module Type**: *DETECTOR*  
**Detector Type**: *all*  
**Status**: Functional  

### Description
This module performs a fit to obtain corrections for non-linear charge sharing, also know as the $`\eta`$-distribution. 
The $`\eta`$-distribution is considered separately for the X and Y axes and is only calculated for clusters with a width of 2 in the respective axis.
The position of the track intercept and the calculated centre of the associated cluster are compared relative to the 2-pixel width of the cluster, and plot in 2D and TProfile histograms.

At the end of the run, fits to the recorded profiles are performed using the provided formulas. 
A printout of the resulting fit parameters is provided in the format read by the EtaCorrection module for convenience.

In order to measure the correct $`\eta`$-distribution, no additional $`\eta`$-correction should be applied before this calculation, i.e. by using the EtaCorrection module.

### Parameters
* `chi2ndof_cut`: Track quality cut on its Chi2 over numbers of degrees of freedom. Default value is `100`.
* `calculate_x` / `calculate_y`: Calculate the correction for X or Y. Default value is `true`.
* `eta_formula_x` / `eta_formula_y`: Formula to for to the recorded $`\eta`$-distributions, defaults to a polynomial of fifth order, i.e. `[0] + [1]*x + [2]*x^2 + [3]*x^3 + [4]*x^4 + [5]*x^5`.

### Plots produced
For each detector the following plots are produced:

* 2D histogram of the calculated $`\eta`$-distribution, for X and Y respectively
* Profile plot of the calculated $`\eta`$-distribution, for X and Y respectively

### Usage
```toml
[EtaCalculation]
chi2ndof_cut = 100
eta_formula_x = [0] + [1]*x + [2]*x^2 + [3]*x^3 + [4]*x^4 + [5]*x^5
```
