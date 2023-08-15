# mg_dose_response_public
data for Hawley et al. dose-response (priming and secondary) experiment 

# data files
- va_hofi_dose_response_aggregated.csv: this file contains the number of infected and total individuals used for each dose to obtain dose-response parameters. The pos and tot are used in the code. The dose_scale is the dose used. Groups are grouped according to priming dose (like in figures)
- deviance_params_21JUL2023.csv contains the estimated parameters for the groups. This file is just used for grabbing the param estimates for the kate_MG_boostrap_resids.R

# code
-kate_define_functions.R - this is all the dose-response models and the likelihood (deviance) functions for fitting
- kate_MG_bootstrap_resids.R: this file is used to obtain the 95% CI of the dose-response parameters. It needs to be better organized so the code doesn't rely on functions in a different file (or write an actual package)

