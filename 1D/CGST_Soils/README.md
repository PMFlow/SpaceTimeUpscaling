This folder contains the following Matlab scrips and files:
#

- 'main_1D_Monod2_Richards_L.m' is the main program for simulations of biodegradation in heterogeneous soils with transition from unsaturated to saturated flow regime modeled by a degenerate Richards equation, solved with an iterative GRW L-schemes.

- 'CG_BGRW_1D_L.m' is the function for the BGRW L-scheme with CGST averaging.

- 'V_Kraichnan_Gauss_param.m' and 'K_r.m' are functions used to compute realizations of the saturated bydraulic conductivity modeled as a lognormal random function.

- 'initstate.m' is the function which initializes the random number generators used in the Kraichnan procedure.

- 'theta_GM.m' is the function which provides the variable water content for saturated/unsaturated flow regimes according the the van Genuchten-Mualem model.

- 'CG_reaction_Monod2.m' is the function for double Monod reactions with two mobile species and constant biomass concentration.

- 'plot_conv' is the function which plots norms of successive approximations used to assess the convergence of the flow and transport L-schemes.

- 'convf', 'convc1', and 'convc2' are files containing convergence norms.

- 'plot_flow.m' is the function which plots the pressure, the water content, and the velocity.

- 'plot_CG2.m' is the function which plots the fine-scale concentrations, the moving averages, and the CGST averages.

- 'errors.m' is the program for the computation of the maximum relative differences between volume and CGST averages.

- 'convf.mat', convc1.mat', and 'convc2.mat' are files containing norms of succesive corrections of the flow and transport solutions during the L-scheme iterative procedure.

- 'qp.mat', 'thtp.mat', and 'pp.mat' are files containing the flow solution (\psi,\theta,q) stored at thre CGST sampling times.

- '1D_T132.mat' is a data file containig the solution at the final time.
