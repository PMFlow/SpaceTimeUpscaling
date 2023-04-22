This folder contains the following Matlab scrips:
#

- 'main_CG_Monod_Sat.m' is the main program for simulations of biodegradation in heterogeneous aquifers

- 'BGRW_2D_React_Monod_Sat.m' is the function for the BGRW transport solver.

- 'V_Kraichnan_Gauss_func.m' and 'V_Kraichnan_Gauss_param.m' are functions used to compute realizations of the Kraichnan velocity filed for log-bydraulic conductivity modeled as a lognormal random function.

- 'initstate.m' is the function which initializes the random number generators used in the Kraichnan procedure.

- 'CG_reaction_Monod2.m' is the function for double Monod reactions with two mobile species and constant biomass concentration.

- 'plot_solution' is the function which plots the components of the Kraichnan velocity and the fine-scale concentrations.

- 'plot_surf_c.m' is the function which plots contours of the fine-scale concentrations at the final time.

- 'plot_CG2.m' is the function which plots the CGST averages.

- 'errors.m' is the function for the computation of the maximum relative differences between volume and CGST averages.
