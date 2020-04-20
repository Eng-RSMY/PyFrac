# -*- coding: utf-8 -*-
"""
This file is part of PyFrac.

Created by Haseeb Zia on 11.05.17.
Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, Geo-Energy Laboratory, 2016-2019.
All rights reserved. See the LICENSE.TXT file for more details.
"""

# tolerances
toleranceFractureFront = 1.0e-3         # tolerance for the fracture front position solver.
toleranceEHL = 1.0e-4                   # tolerance for the elastohydrodynamic system solver.
tol_projection = 2.5e-3                 # tolerance for the toughness iteration.

# max iterations
max_front_itrs = 25                     # maximum iterations for the fracture front.
max_solver_itrs = 140                   # maximum iterations for the elastohydrodynamic solver.
max_proj_Itrs = 10                      # maximum projection iterations.

# time and time stepping
tmStp_prefactor = 0.8                   # time step prefactor(pf) to determine the time step(dt = pf*min(dx, dy)/max(v).
req_sol_at = None                       # times at which the solution is required.
final_time = None                       # time to stop the propagation.
maximum_steps = 2000                    # maximum time steps.
timeStep_limit = None                   # limit for the time step.
fixed_time_step = None                  # constant time step.

# time step re-attempt
max_reattemps = 8                       # maximum reattempts in case of time step failure.
reattempt_factor = 0.8                  # the factor by which time step is reduced on reattempts.

# output
plot_figure = True                      # if True, figure will be plotted after the given time period.
save_to_disk = True                     # if True, fracture will be saved after the given time period.
output_folder = None                    # the address to save the output data.
plot_analytical = False                 # if True, analytical solution will also be plotted.
analytical_sol = None                   # the analytical solution to be plotted.
bck_color = None                        # the parameter according to which background is color coded (see class doc).
sim_name = None                         # name given to the simulation.
block_figure = False                    # if true, the simulation will proceed after the figure is closed.
plot_var = ['w']                        # the list of variables to be plotted during simulation.
plot_proj = '2D_clrmap'                 # projection to be plotted with.
plot_time_period = None                 # the time period after which the variables given in plot_var are plotted.
plot_TS_jump = 1                        # the number of time steps after which the given variables are plotted.
save_time_period = None                 # the time period after which the output is saved to disk.
save_TS_jump = 1                        # the number of time steps after which the output is saved to disk.

# type of solver
mech_loading = False                    # if True, the mechanical loading solver will be used.
volume_control = False                  # if True, the volume control solver will be used.
viscous_injection = True                # if True, the viscous fluid solver solver will be used.
substitute_pressure = True              # if True, the pressure will be substituted with width to make the EHL system.
solve_deltaP = True                     # if True, the change in pressure, instead of pressure will be solved.
solve_stagnant_tip = False              # if True, stagnant tip cells will also be solved for
solve_tip_corr_rib = True               # if True, the corresponding tip cells to closed ribbon cells will be solved.
solve_sparse = None                     # if True, the fluid conductivity matrix will be made with sparse matrix.
elastohydr_solver = 'anderson'          # default non-linear solver.

# miscellaneous
tip_asymptote = 'U'                     # the tip_asymptote to be used (see class documentation for details).
save_regime = False                     # if True, the the regime of the ribbon cells will also be saved.
verbosity = 2                           # the level of details about the ongoing simulation to be plotted.
enable_remeshing = True                 # if true, computation domain will be remeshed after reaching end of the domain.
remesh_factor = 2.                      # the factor by which the mesh is compressed.
front_advancing = 'predictor-corrector' # possible options include 'implicit', 'explicit' and 'predictor-corrector'.
collect_perf_data = False               # if True, performance data will be collected in the form of a tree.
param_from_tip = False                  # set the space dependant tip parameters to be taken from ribbon cell.
save_ReyNumb = False                    # if True, the Reynold's number at each edge will be saved.
save_fluid_flux = False                 # if True, the fluid flux at each edge will be saved.
save_fluid_vel = False                  # if True, the fluid vel at each edge will be saved.
gravity = False                         # if True, the effect of gravity will be taken into account.
TI_Kernel_exec_path = '../TI_Kernel/build' # the folder containing the executable to calculate TI elasticity matrix.
explicit_projection = False             # if True, direction from last time step will be used to evaluate TI parameters.
symmetric = False                       # if True, only positive quarter of the cartesian coordinates will be solved.
proj_method = 'ILSA_orig'               # set the method to evaluate projection on front to the original ILSA method.

# fracture geometry
height = None                           # fracture height to calculate the analytical solution for PKN or KGD geometry.
aspect_ratio = None                     # fracture aspect ratio to calculate the analytical solution for TI case.