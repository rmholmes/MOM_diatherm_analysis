# MatLab MOM5 diathermal heat transport analysis routines

This repository contains MatLab analysis routines for analysing
diathermal heat transport in MOM5 simulations. These routines were
used for the analysis described in [*Holmes et al (2019) Diathermal
heat transport in a global ocean model, J. Phys. Oceanogr. 49,
141-161*](https://doi.org/10.1175/JPO-D-18-0098.1), Holmes et al
(2019) Atlantic ocean heat transport enabled by Indo-Pacific heat
uptake and mixing, Geophysical Research Letters, accepted and several
other current projects.

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

If you run into any problems please email me at
`ryan.holmes@unsw.edu.au`


