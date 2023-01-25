This folder contains the following Matlab scrips and files:


'main_1D_Monod2_rand.m' is the main program for simulations of biodegradation in heterogeneous aquifers.

'CG_BGRW_1D.m' is the function for the BGRW solver of the advection-dispersion steps.

'V_Kraichnan_Gauss_param.m' and 'V_Kraichnan_Gauss_func.m' are functions used to compute realizations of the saturated bydraulic conductivity modeled as a lognormal random function.

'initstate.m' is the function which initializes the random number generators used in the Kraichnan procedure.

'CG_reaction_Monod2.m' is the function for double Monod reactions with two mobile species and constant biomass.

'plot_CG_Monod2.m' is the function which plots the fine-scale concentration, the moving averages, and the CGST averages.

'errors_aqv.m' is the program for the computation of the maximum relative differences between volume and CGST averages.

'CG_aquifer.mat' and 'CG_aquifer_ens.mat' are data files containing results for single realization and ensembles of 100 realizations, respectively.

