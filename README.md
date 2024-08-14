# Title: Prior exposure to pathogens augments host heterogeneity in susceptibility and has key epidemiological consequences
# Authors: Dana M. Hawley*, Anna M. Pérez-Umphrey, James S. Adelman, Arietta E. Fleming-Davies, Jesse Garrett-Larsen, Steven J. Geary, Lauren M. Childs$, Kate E. Langwig$,*
# *co-corresponding authors, $co-senior authors


# Files provided:
#     README.md
# Code:
#     CIs_for_fitting_doseresponse.R - makes bootstrap samples and fits dose response
#     fitting_doseresponse.R - fits basic dose response
#     fitting_eyescore.R - fits eyescore data to determine mortality
#     likelihood_functions.R - likelihood functions used in fitting
#     SIR_density.R - epidemic model trajectories from base fitting
#     SIR_dist.R - epidemic model differences from bootstrap samples
# Data:
#     day_41_removed_hofi_small_for_fitting_secondary.csv
#     deviance_params_noday41positives.csv - contains the estimated parameters for the groups
#     eye_score_primary_dose.csv - contains aggregated eye score from primary dose
#     eye_score_secondary_dose.csv - contains aggregated eye score from secondary dose
#     mortality.csv - output of fitting_eyescore.R
#     va0_gamma_distributionCIs.csv - output of CIs_for_fitting_doseresponse.R
#     va750_gamma_distributionCIs.csv - output of CIs_for_fitting_doseresponse.R
#     va30000_gamma_distributionCIs.csv - output of CIs_for_fitting_doseresponse.R



#### FIT MORTALITY
source("fitting_eyescore.R")
# requires:
#     eye_score_primary_dose.csv
#     eye_score_secondary_dose.csv
# produces:
#     mortality.csv


#### MAKE BOOTSTRAP SAMPLES AND FIT DOSE RESPONSE
source("CIs_for_fitting_doseresponse.R")
# requires:
#     day_41_removed_hofi_small_for_fitting_secondary.csv
#     deviance_params_noday41positives.csv
# 	  likelihood_functions.R
# produces:
#     all_psuedosets.csv
#     va30000_gamma_distributionCIs_new.csv
#     va750_gamma_distributionCIs_new.csv
#     va0_gamma_distributionCIs_new.csv


#### FIT BASIC DOSE RESPONSE
source("fitting_doseresponse.R")
# requires:
#     day_41_removed_hofi_small_for_fitting_secondary.csv
# 	  likelihood_functions.R
#     va30000_gamma_distributionCIs.csv
#     va750_gamma_distributionCIs.csv
#     va0_gamma_distributionCIs.csv
# produces:
#     deviance_params_noday41positives.csv
#     Fig/va_dose_response_fits_error_bars.png
#     Fig/heterogeneous_gamma_fits_dose-response.png
#     Fig/homogeneous_fits_dose-response.png


#### SIR MODEL - BAR PLOTS & TRAJECTORIES
source("SIR_dist.R")
# requires:
#       mortality.csv
#       deviance_params_noday41positives.csv
# produces:
#       SIR_Output_fitted.csv
#       SIR_Output_fixed.csv
#       Fig/Bar_fittedmortality.png
#       Fig/Bar_fixedmortality.png
#       Fig/Trajectories_fittedmortality.png
#       Fig/Trajectories_fixedmortality.png


#### SIR model - HISTOGRAMS
source("SIR_density.R")
# requires:
#       all_pseudosets.csv
#       mortality.csv
# produces:
#       Fig/Density_fittedmortality.png
#       Fig/Density_fixedmortality.png
