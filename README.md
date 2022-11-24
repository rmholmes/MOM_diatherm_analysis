# MatLab MOM5 diathermal heat transport analysis routines

This repository contains MatLab analysis routines for analysing
diathermal heat transport in MOM5 simulations. These routines were
used for the analysis described in the following papers:

- [Holmes et al (2019) Diathermal heat transport in a global ocean model, J. Phys. Oceanogr. 49, 141-161](https://doi.org/10.1175/JPO-D-18-0098.1)

- [Holmes et al (2019) Atlantic ocean heat transport enabled by Indo-Pacific heat uptake and mixing, Geophysical Research Letters, 46, 13939-13949](https://doi.org/10.1029/2019GL085160)

- [Holmes et al (2021) The geography of numerical mixing in a suite of global ocean models, Journal of Advances in Modeling Earth Systems, 13, e2020MS002333](https://doi.org/10.1029/2020MS002333) 

- [Holmes et al. (2022) Adiabatic and diabatic signatures of ocean temperature variability, Journal of Climate, 35 (5), 1459-1477](https://doi.org/10.1175/JCLI-D-21-0695.1)

## Workflow

The analysis routine `Heat_Budget_Process_Production.m` is run with
the model output .nc files as inputs (typically run on a single year
with monthly output). `Heat_Budget_Process_Production.m` contains a
number of sections for analysing different aspects of the heat
budget. The results are saved into `.mat` files which are then used to
make plots through the various plotting routines.

## Specific instructions for reproducing the latitude-temperature analysis in Holmes et al. (2019, GRL):

The `2019GL085160` branch contains clean versions of
`Heat_Budget_Process_Production.m`, `Heat_Budget_Mask.m` and
`Heat_Budget_Plot_LatT.m`, which were used for the
latitude-temperature space analysis in Holmes et al. (2019,
GRL). Running `Heat_Budget_Process_Production.m` (which depends on
`Heat_Budget_Mask.m` to produce the Atlantic and Indo-Pacific masks)
produces .mat files which can then be used by
`Heat_Budget_Plot_LatT.m` to plot latitude-temperature slices of the
various diathermal and meridional heat fluxes.

If you already have the processed .mat files (averaged over the last
10 years of the MOM025 Control simulation) which are publically
available at [Research Data
Australia](https://doi.org/10.26190/5dc23d4b7e739), then the plotting
routine Heat_Budget_Plot_LatT.m should produce the figures used in the
GRL article.

If you run into any problems please email me at
`r.holmes@sydney.edu.au`


