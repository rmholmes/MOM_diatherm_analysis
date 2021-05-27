# MatLab MOM5 diathermal heat transport analysis routines

This repository contains MatLab analysis routines for analysing
diathermal heat transport in MOM5 simulations. These routines were
used for the analysis described in [*Holmes et al (2019) Diathermal
heat transport in a global ocean model, J. Phys. Oceanogr. 49,
141-161*](https://doi.org/10.1175/JPO-D-18-0098.1), [*Holmes et al
(2019) Atlantic ocean heat transport enabled by Indo-Pacific heat
uptake and mixing, Geophysical Research Letters, 46,
13939-13949*](https://doi.org/10.1029/2019GL085160) and [*Holmes et al
(2021) The geography of numerical mixing in a suite of global ocean
models, in review at the Journal of Advances in Modeling Earth
Systems](https://www.essoar.org/doi/10.1002/essoar.10504439.1).

## Workflow

The analysis routine `Heat_Budget_Process_Production.m` is run with
the model output .nc files as inputs (typically run on a single year
with monthly output). `Heat_Budget_Process_Production.m` contains a
number of sections for analysing different aspects of the heat
budget. The results are saved into `.mat` files which are then used to
make plots through the various plotting routines. Note that the online
heat budget diagnostics are contained in the [MOM5 code github
repository](https://github.com/mom-ocean/MOM5).

## Specific instructions for reproducing analysis in Holmes et al. (2021):

The `2020MS002333` branch contains matlab files for working with the
processed .mat data produced by `Heat_Budget_Process_Production.m` and
`Heat_Budget_Process_Variances.m`. The raw netcdf data is too much to
provide on a public server (given that 13 different model
configurations were used in the paper) but any reasonable requests for
data access will be granted. Note that general output (not including
specialised numerical mixing/heat budget diagnostics) is available
through the [Consortium for Ocean-Sea Ice Modelling in
Australia](http://cosima.org.au/), e.g. see [Research Data
Australia](https://researchdata.edu.au/cosima-model-output-collection/993052).

The processed .mat files have been provided on
(Zenodo)[https://dx.doi.org/10.5281/zenodo.4798380]. This branch of
this repository contains stripped down processing scripts for
producing most of the analysis provided in Holmes et al. (2021,
JAMES). The routines are:

- `Heat_Budget_Plot_Global.m`: This script produces the
  globally-averaged heat budget in temperature coordinates
  (e.g. Figs. 2, 7 and 15 of Holmes et al. 2021), including the
  numerical mixing term, for all 13 configurations. This script also
  outputs the globally-integrated numerical mixing and explicit mixing
  metrics Inet, Mnet and Nnet used in the summary Fig 17.

- `Heat_Budget_Plot_Spatial.m`: This script produces the
  longitude-latitude spatial structure maps (e.g. Figs. 3, 8, 10, A.1)
  of either numerical or vertical mixing. Data is provided for all 13
  configurations for the 5, 15 and 22.5C isotherms.

- `Heat_Budget_Plot_EqSlice.m`: This script produces equatorial
  Pacific longitude-depth plots of either numerical or vertical mixing
  (e.g. Figs. 6, 12 and 14).

- `Heat_Budget_Plot_ZonalAverage.m`: This script plots the global zonal
  average of vertical and numerical mixing in Fig. 4.

- `Heat_Budget_Plot_Variances.m`: This script plots the grid-scale
  variance diagnostics in Fig. 5.

- `Heat_Budget_Plot_Line.m`: This script plots the regional averages
  in the Gulf Stream and Tropical Pacific regions (Figs. 9 and 11).

- Note that the script for the spectral analysis (Fig. 16) is located
  [elsewhere](https://github.com/rmholmes/cosima-scripts/blob/master/spectra_simple.ipynb),
  but in order to use this script one needs access to the NCI system
  (if required, I can provide this access).

If you run into any problems please contact me at
`ryan.holmes@unsw.edu.au`


