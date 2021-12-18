This folder contains the following Matlab scrips and files:

#
- 'main_1D_Monod2_rand.m' is the main program for simulations of biodegradation in heterogeneous aquifers

- 'CG_BGRW_1D.m' is the function for the BGRW solver of the advection-dispersion steps.

- 'V_Kraichnan_Gauss_param.m' and 'V_Kraichnan_Gauss_func.m' are functions used to compute realizations of 
	the saturated bydraulic conductivity modeled as a lognormal random function.

- 'initstate.m' is the function which initializes the random number generators used in the Kraichnan procedure.

- 'CG_reaction_Monod2.m' is the function for double Monod reactions with two mobile species and constant biomass.

- 'plot_CG2.m' is the function which plots the fine-scale concentration, the moving averages, and the CGST averages.
