# livneh-badger-ncc-code
This repository contains processing scripts used to evaluate future changes in the Equitable Threat Score relative to historical conditions, for snow-dominated portions of the western U.S.

Below is a description of the companion codes and data used to conduct a predictability analysis.

Grid_Point_Predictions_SWE.r: Uses publicly available downscaled hydrologic projections to apply a linear predictive equation using SVD at all grid points that meet our selection criteria (see the manuscript for details on this and for links to the raw data). 
The predictor is April 1 SWE and the predictand is AMJJ total runoff                                                        #
Each model outputs three text files that lists each point together with their geographical coordinates useful for plotting and analysis purposes. The files are:
1) General output
2) Equitable Threat Score (ETS) information
3) Goodness of fit metrics

Grid_Point_Predictions_Kitchen_Sink.r: Is an identical script to Grid_Point_Predictions_SWE.r, with the exception that it uses many predictors (hence the term "Kitchen Sink"), rather than just SWE as was the case for the aforementioned script. The specific predictors used are: April 1 SWE, April 1 water-year-to-date precipitation, April 1 water-year-to-date mean temperature, and April 1 Soil Moisture. The predictand is AMJJ total runoff. The output files and their formats are identical to Grid_Point_Predictions_SWE.r.

Jacknife_Selected_Gages_GageLoc.r: This script uses observed April 1 SWE from all combinations of 4 SNOTEL stations (e.g. x-choose-4) to predict observed streamflow for each watershed. A 'leave-one-out' cross validation is used to calibrate and validate the predictive model. The method is repeated using model data, where the SWE "observations" come from the grid cell that encloses the SNOTEL station, and where observed runoff comes from the total simulated runoff over that basin.

HUC_ID.txt: is a list of all Hydrologic Unit Codes (HUCs) in the western U.S., which is larger than, but inclusive of the 24 HUCs used in the validation

REEDS_HUC.nc: includes the polygon shapefiles of the HUCs mapped on to the model resolution.

