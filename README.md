# MatLab MOM5 diathermal heat transport analysis routines

This repository contains MatLab analysis routines for analysing diathermal heat transport in MOM5 simulations. These routines were used for the analysis described in [*Holmes et al (2019) Diathermal heat transport in a global ocean model, J. Phys. Oceanogr. 49, 141-161*](https://doi.org/10.1175/JPO-D-18-0098.1), and for several in-prep. papers.

## Workflow

The analysis routine `Heat_Budget_Process_Production.m` is run with the model output .nc files as inputs (typically run on a single year with monthly output). `Heat_Budget_Process_Production.m` contains a number of sections for analysing different aspects of the heat budget. The results are saved into `.mat` files which are then used to make plots through the various plotting routines.


