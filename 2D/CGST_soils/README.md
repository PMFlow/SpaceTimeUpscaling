This folder contains the following Matlab scrips and files:
#

- 'main_2D_Monod2_Richards_L.m' is the main program for simulations of biodegradation in heterogeneous soils
	with transition from unsaturated to saturated flow regime modeled by a degenerate Richards equation,
	solved with iterative GRW L-schemes.

- 'CG_BGRW_2D_L.m' is the function for the BGRW L-scheme with CGST averaging.

- 'V_Kraichnan_Gauss_param.m' and 'K_r.m' are functions used to compute realizations of the saturated bydraulic
	conductivity modeled as a lognormal random function.

- 'velocity.m' is afunction to compute velocities from pressure solutions.

- 'initstate.m' is the function which initializes the random number generators used in the Kraichnan procedure.

- 'theta_GM.m' is the function which provides the variable water content for saturated/unsaturated flow regimes
	according the the van Genuchten-Mualem model.

- 'CG_reaction_Monod2.m' is the function for double Monod reactions with two mobile species and constant biomass concentration.

- 'plot_conv_Monod' is the function which plots norms of successive approximations used to assess the convergence 
	of the flow and transport L-schemes.

- 'subplot_surf_ptht.m' is the function which plots the pressure, the water content, and the velocity at the final time.

- 'subplot_surf_c_2.m' is the function which plots the fine-scale concentrations at the final time.

- 'plot_CG2.m' is the function which plots the CGST averages.

- 'errors.m' is the program for the computation of the maximum relative differences between volume and CGST averages.

- 'convf', 'convc1', and 'convc2' are files containing convergence norms.

- '2D_T132_centr.mat' and '2D_T132_dec.mat' are data files containig the solution at the final time for centered and decentered sampling lines, respectively.
