# -*- coding: utf-8 -*-
"""
This file is part of PyFrac.

Created by Haseeb Zia on Fri June 16 17:49:21 2017.
Copyright (c) "ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, Geo-Energy Laboratory", 2016-2019. All rights
reserved. See the LICENSE.TXT file for more details.
"""


# imports
from src.Fracture import *
from src.Controller import *
from src.anisotropy import *

# creating mesh
Mesh = CartesianMesh(6., 3., 51, 51)

# solid properties
Cij = np.zeros((6, 6), dtype=float)

# modelgam05
Cij[0, 0] = 47.5
Cij[0, 1] = 16.55
Cij[5, 5] = 0.5*(Cij[0, 0]-Cij[0, 1])
Cij[0, 2] = 20.34
Cij[2, 2] = 28.86
Cij[3, 3] = 7.32
Cij[1, 1] = Cij[0, 0]
Cij[1, 0] = Cij[0, 1]
Cij[2, 0] = Cij[0, 2]
Cij[2, 1] = Cij[0, 2]
Cij[1, 2] = Cij[0, 2]
Cij[4, 4] = Cij[3, 3]
Cij=Cij*1e9

Eprime = TI_plain_strain_modulus(np.pi/2, Cij) # plain strain modulus

# the function below will make the fracture propagate in the form of an ellipse (see Zia and Lecampion 2018)
def K1c_func(alpha):
    """ function giving the dependence of fracture toughness on propagation direction alpha"""

    K1c_3 = 2e6                     # fracture toughness along x-axis
    gamma = 2.0                     # aspect ratio

    Cij = np.zeros((6, 6), dtype=float)
    Cij[0, 0] = 47.5
    Cij[0, 1] = 16.55
    Cij[5, 5] = 0.5 * (Cij[0, 0] - Cij[0, 1])
    Cij[0, 2] = 20.34
    Cij[2, 2] = 28.86
    Cij[3, 3] = 7.32
    Cij[1, 1] = Cij[0, 0]
    Cij[1, 0] = Cij[0, 1]
    Cij[2, 0] = Cij[0, 2]
    Cij[2, 1] = Cij[0, 2]
    Cij[1, 2] = Cij[0, 2]
    Cij[4, 4] = Cij[3, 3]
    Cij = Cij * 1e9

    Eprime_ratio = TI_plain_strain_modulus(alpha,Cij) / TI_plain_strain_modulus(np.pi/2, Cij)
    beta = np.arctan(np.tan(alpha) / gamma)

    return K1c_3 * Eprime_ratio * ((np.sin(beta))**2 + (np.cos(beta)/gamma)**2)**0.25

# materila properties
Solid = MaterialProperties(Mesh,
                           Eprime=Eprime,
                           anisotropic_K1c=True,
                           Toughness=K1c_func(np.pi/2),
                           K1c_func=K1c_func,
                           TI_elasticity=True,
                           Cij=Cij)

# injection parameters
Q0 = 0.001  # injection rate
Injection = InjectionProperties(Q0, Mesh)

# fluid properties
Fluid = FluidProperties(viscosity=1.1e-5)

# aspect ratio of the elliptical fracture
gamma = 2.0

# simulation properties
simulProp = SimulationProperties()
simulProp.finalTime = 75                # the time at which the simulation stops
simulProp.set_tipAsymptote("K")         # the tip asymptote is evaluated with the toughness dominated assumption
simulProp.outputTimePeriod = 1e-3       # save after every time step
simulProp.tolFractFront = 4e-3          # increase tolerance for the anisotropic case
simulProp.maxToughnessItr = 5           # set maximum iterations to 5 for faster simulation
simulProp.remeshFactor = 1.5            # the factor by which the mesh will be compressed.
simulProp.set_volumeControl(True)       # assume inviscid fluid
simulProp.explicitProjection = True     # do not iterate on projection to find propagation direction
simulProp.aspectRatio = gamma           # aspect ratio of the fracture
simulProp.set_outputFolder("./data/TI_elasticity_ellipse")
simulProp.set_simulation_name('TI_ellasticy_benchmark')
simulProp.TI_KernelExecPath = '../src_TI_Kernel/' # path to the executable that calculates TI kernel
simulProp.verbosity = 2

# initialization parameters
Fr_geometry = Geometry('elliptical',
                       minor_axis=0.75,
                       gamma=gamma)
init_param = InitializationParameters(Fr_geometry, regime='E_E')

#  creating fracture object
Fr = Fracture(Mesh,
              init_param,
              Solid,
              Fluid,
              Injection,
              simulProp)

# create a Controller
controller = Controller(Fr,
                        Solid,
                        Fluid,
                        Injection,
                        simulProp)

#run the simulation
controller.run()

####################
# plotting results #
####################

# loading simulation results
time_srs = 2 ** np.linspace(np.log2(0.38), np.log2(74), 10)
Fr_list, properties = load_fractures(address='./data/TI_elasticity_ellipse',
                                            simulation='TI_ellasticy_benchmark',
                                            time_srs=time_srs)
time_srs = get_fracture_variable(Fr_list,
                                 variable='time')

Fig_FP = plot_fracture_list(Fr_list,
                            variable='mesh',
                            projection='2D')
Fig_FP = plot_fracture_list(Fr_list,
                            variable='footprint',
                            projection='2D',
                            fig=Fig_FP)
Fig_FP = plot_analytical_solution('E_E',
                                  'footprint',
                                  Solid,
                                  Injection,
                                  fluid_prop=Fluid,
                                  fig=Fig_FP,
                                  projection='2D',
                                  time_srs=time_srs,
                                  gamma=gamma)

time_srs = 2 ** np.linspace(np.log2(0.38), np.log2(74), 5)
Fr_list, properties = load_fractures(address='./data/TI_elasticity_ellipse',
                                            simulation='TI_ellasticy_benchmark',
                                            time_srs=time_srs)
time_srs = get_fracture_variable(Fr_list,
                                 variable='time')

ext_pnts = np.empty((2, 2), dtype=np.float64)
Fig_w_slice = plot_fracture_list_slice(Fr_list,
                                       variable='width',
                                       point1=[-Fr_list[-1].mesh.Lx, 0],
                                       plot_cell_center=True,
                                       orientation='horizontal',
                                       extreme_points=ext_pnts)
Fig_w_slice = plot_analytical_solution_slice('E_E',
                                             variable='width',
                                             mat_prop=Solid,
                                             inj_prop=Injection,
                                             fluid_prop=Fluid,
                                             fig=Fig_w_slice,
                                             point1=ext_pnts[0],
                                             point2=ext_pnts[1],
                                             time_srs=time_srs,
                                             gamma=gamma)


Fr_list, properties = load_fractures(address='./data/TI_elasticity_ellipse',
                                            simulation='TI_ellasticy_benchmark')
time_srs = get_fracture_variable(Fr_list,
                                 variable='time')
plot_prop = PlotProperties(line_style='.',
                           graph_scaling='loglog')

labels = LabelProperties('d_min', 'wm', '1D')
labels.figLabel = 'Minor axis length'
Fig_len_a = plot_fracture_list(Fr_list,
                             variable='d_min',
                             plot_prop=plot_prop,
                             labels=labels)
Fig_len_a = plot_analytical_solution('E_E',
                                   'd_min',
                                   Solid,
                                   Injection,
                                   fluid_prop=Fluid,
                                   fig=Fig_len_a,
                                   time_srs=time_srs,
                                   gamma=gamma,
                                   labels=labels)

labels.figLabel = 'Major axis length'
Fig_len_b = plot_fracture_list(Fr_list,
                             variable='d_max',
                             plot_prop=plot_prop,
                             labels=labels)
Fig_len_b = plot_analytical_solution('E_E',
                                   'd_max',
                                   Solid,
                                   Injection,
                                   fluid_prop=Fluid,
                                   fig=Fig_len_b,
                                   time_srs=time_srs,
                                   gamma=gamma,
                                   labels=labels)

plt.show(block=True)
