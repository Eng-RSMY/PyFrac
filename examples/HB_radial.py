# -*- coding: utf-8 -*-
"""
This file is part of PyFrac.

Created by Haseeb Zia on Fri June 16 17:49:21 2017.
Copyright (c) "ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, Geo-Energy Laboratory", 2016-2019.
All rights reserved. See the LICENSE.TXT file for more details.
"""

# local imports
from mesh import CartesianMesh
from properties import MaterialProperties, FluidProperties, InjectionProperties, SimulationProperties
from fracture import Fracture
from controller import Controller
from fracture_initialization import Geometry, InitializationParameters
from visualization import *

# creating mesh
Mesh = CartesianMesh(.25, .25, 61, 61)

# solid properties
nu = 0.15                            # Poisson's ratio
youngs_mod = 3e10                 # Young's modulus
Eprime = youngs_mod / (1 - nu ** 2) # plain strain modulus
K_Ic = 2e6
sigma0 = 15e6                          # fracture toughness

# material properties
Solid = MaterialProperties(Mesh,
                           Eprime,
                           K_Ic, 
                           # Carters_coef = 1e-6
                           confining_stress=sigma0,
                           # minimum_width=1e5
                           )

# injection parameters

def sink_location(x, y):
    return x >= 1. and  x <= 1.5

def sink_vel(x, y):
    return 2e-4

Q0 = 0.001/60
Injection = InjectionProperties(Q0,
                                Mesh,
                                # sink_loc_func=sink_location, sink_vel_func=sink_vel
                                )

# fluid properties
Fluid = FluidProperties(viscosity=0.01, rheology='HBF', n=1.0, k=0.01, T0=5)

# simulation properties
simulProp = SimulationProperties()
simulProp.finalTime = 4*60                           # the time at which the simulation stops
# simulProp.frontAdvancing = 'implicit'               # to set explicit front tracking
# simulProp.saveTSJump, simulProp.plotTSJump = 5, 5   # save and plot after every five time steps
simulProp.set_outputFolder("./Data/HB") # the disk address where the files are saved
simulProp.tmStpPrefactor = 0.45
simulProp.maxSolverItrs = 250
# simulProp.plotFigure = False
simulProp.set_tipAsymptote('HBF')
simulProp.elastohydrSolver = 'anderson'
simulProp.tolFractFront = 0.01
simulProp.solveSparse = False
simulProp.set_simulation_name('Bingham_1e0_1e-1_5_61')
simulProp.plotVar = ['ev']
# simulProp.saveToDisk = False


# initialization parameters
# Fr_geometry = Geometry('radial', radius=0.2)
# from elasticity import load_isotropic_elasticity_matrix
# C = load_isotropic_elasticity_matrix(Mesh, Eprime)
# init_param = InitializationParameters(Fr_geometry, regime='static', net_pressure='3e6', elasticity_matrix=C)

Fr_geometry = Geometry('radial', radius=.2)
init_param = InitializationParameters(Fr_geometry, regime='M')

# creating fracture object
Fr = Fracture(Mesh,
              init_param,
              Solid,
              Fluid,
              Injection,
              simulProp)
# Fr = load_fractures(address="./Data/HB", sim_name='sheer_thining_75-1-7_lk')[0][-1]

# create a Controller
controller = Controller(Fr,
                        Solid,
                        Fluid,
                        Injection,
                        simulProp)

# run the simulation
# controller.run()

####################
# plotting results #
####################

from visualization import *

Fig_R = None
Fig_w = None

# loading simulation results
Fr_list, properties = load_fractures(address="./Data/HB", sim_name="Bingham_1e0_1e-1_5_61")        # load all fractures
time_srs = get_fracture_variable(Fr_list, variable='time')                      # list of times


animate_simulation_results(Fr_list, variable='w')
# plot fracture radius
plot_prop = PlotProperties(line_style='.', graph_scaling='loglog')
Fig_R = plot_fracture_list(Fr_list,
                           variable='d_mean',
                           plot_prop=plot_prop,
                           fig=Fig_R)
# plot analytical radius
Fig_R = plot_analytical_solution(regime='M',
                                 variable='d_mean',
                                 mat_prop=Solid,
                                 inj_prop=Injection,
                                 fluid_prop=Fluid,
                                 time_srs=time_srs,
                                 fig=Fig_R)

# # plot analytical radius
# plot_prop_k = PlotProperties(line_color='b')
# Fig_R = plot_analytical_solution(regime='K',
#                                  variable='d_mean',
#                                  mat_prop=Solid,
#                                  inj_prop=Injection,
#                                  fluid_prop=Fluid,
#                                  time_srs=time_srs,
#                                  fig=Fig_R,
#                                  plot_prop=plot_prop)

# # plot width at center
# Fig_w = plot_fracture_list_at_point(Fr_list,
#                                     variable='w',
#                                     plot_prop=plot_prop,
#                                     fig=Fig_w)
# # plot analytical width at center
# Fig_w = plot_analytical_solution_at_point('M',
#                                           'w',
#                                           Solid,
#                                           Injection,
#                                           fluid_prop=Fluid,
#                                           time_srs=time_srs,
#                                           fig=Fig_w)


# time_srs = np.asarray([2, 200, 5000, 30000, 100000])
# Fr_list, properties = load_fractures(address="./Data/M_radial_explicit",
#                                      time_srs=time_srs)
# time_srs = get_fracture_variable(Fr_list,
#                                  variable='time')
#
# # plot footprint
# Fig_FP = plot_fracture_list(Fr_list,
#                             variable='mesh',
#                             projection='2D')
# Fig_FP = plot_fracture_list(Fr_list,
#                             variable='footprint',
#                             projection='2D',
#                             fig=Fig_FP)
# # plot analytical footprint
# Fig_FP = plot_analytical_solution('M',
#                                   'footprint',
#                                   Solid,
#                                   Injection,
#                                   fluid_prop=Fluid,
#                                   time_srs=time_srs,
#                                   projection='2D',
#                                   fig=Fig_FP)
#
#
# # plot slice
# ext_pnts = np.empty((2, 2), dtype=np.float64)
# Fig_WS = plot_fracture_list_slice(Fr_list,
#                                   variable='w',
#                                   projection='2D',
#                                   plot_cell_center=True,
#                                   extreme_points=ext_pnts)
# # plot slice analytical
# Fig_WS = plot_analytical_solution_slice('M',
#                                         'w',
#                                         Solid,
#                                         Injection,
#                                         fluid_prop=Fluid,
#                                         fig=Fig_WS,
#                                         time_srs=time_srs,
#                                         point1=ext_pnts[0],
#                                         point2=ext_pnts[1])

# #plotting in 3D
# Fig_Fr = plot_fracture_list(Fr_list,
#                             variable='mesh',
#                             projection='3D')
# Fig_Fr = plot_fracture_list(Fr_list,
#                             variable='width',
#                             projection='3D',
#                             fig=Fig_Fr)
# Fig_Fr = plot_fracture_list(Fr_list,
#                             variable='footprint',
#                             projection='3D',
#                             fig=Fig_Fr)

plt.show(block=True)

