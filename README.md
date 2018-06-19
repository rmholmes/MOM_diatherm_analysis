# MOM Diathermal Diagnostics matlab files

for_Jan branch: I've modified the Heat_Budget_Plot routine to
simplify. When you have the following .mat files in base (replace 000
with output of interest, 000 should contain 5 years of data):

* ACCESS-OM2_1deg_jra55_ryf8485_kds50_may_output000_BaseVars.mat - lon, lat, depth, temp arrays etc. 
* ACCESS-OM2_1deg_jra55_ryf8485_kds50_may_output000_GlobalHBud.mat - Globally-integrated budget
* ACCESS-OM2_1deg_jra55_ryf8485_kds50_may_output000_VertInt_T22p5C.mat - Lon-lat fields of heat fluxes at 22.5C
* ACCESS-OM2_1deg_jra55_ryf8485_kds50_may_output000_WMT_T22p5C.mat - Lon-lat fields of WMT velocities at 22.5C
* ACCESS-OM2_1deg_jra55_ryf8485_kds50_may_output000_SurfaceVars.mat - Lon-lat fields of SST, shflux etc.

(which have been processed by Heat_Budget_Process_Production.m)

Then running Heat_Budget_Plot should give you a few figures:

* The global heat budget
* temperature-time plot of surface forcing field
* Spatial Lon-lat plot of heat flux due to vertical mixing across 22.5C 
* Spatial Lon-lat plot of surface heat flux and SST contours.

There's also code there to plot SST biases from WOA13 (commented out at the moment). 

