# -*- coding: utf-8 -*-
"""
This file is part of PyFrac.

Created by Haseeb Zia on Fri Dec 23 17:49:21 2016.
Copyright (c) "ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, Geo-Energy Laboratory", 2016-2017. All rights reserved.
See the LICENSE.TXT file for more details.
"""


# adding src folder to the path
import sys
if "win" in sys.platform:
    slash = "\\"
else:
    slash = "/"
if not '.' + slash + 'src' in sys.path:
    sys.path.append('.' + slash + 'src')

# imports
import numpy as np
from src.CartesianMesh import *
from src.Fracture import *
from src.LevelSet import *
from src.VolIntegral import *
from src.Elasticity import *
from src.Properties import *
from src.FractureFrontLoop import *
from src.Controller import *
from src.PostProcess import animate_simulation_results

# creating mesh
Mesh = CartesianMesh(3, 3, 41, 41)

# solid properties
nu = 0.4
Eprime = 3.3e10 / (1 - nu ** 2)
K_Ic = 1e6

sigma0 = np.full((Mesh.NumberOfElts,), 1e6, dtype=np.float64)
# high stressed layers
stressed_layer = np.where(abs(Mesh.CenterCoor[:,1]) > 1.5-Mesh.hx/2)[0]
sigma0[stressed_layer] = 2e6

d_grain = 1e-5
Solid = MaterialProperties(Eprime, K_Ic, 0., sigma0, d_grain, Mesh)

# injection parameters
Q0 = 0.001  # injection rate
well_location = np.array([0., 0.])
Injection = InjectionProperties(Q0, well_location, Mesh)

# fluid properties
Fluid = FluidProperties(1.1e-3, Mesh, turbulence=False)

# simulation properties
simulProp = SimulationParameters(tip_asymptote="U",
                                 output_time_period=0.02,
                                 plot_figure=False,
                                 save_to_disk=True,
                                 out_file_folder=".\\Data\\Test", # e.g. "./Data/Laminar" for linux or mac
                                 plot_analytical=True,
                                 tmStp_prefactor=0.4)


# initializing fracture
initRad = 0.8 # initial radius of fracture

# creating fracture object
Fr = Fracture(Mesh,
              initRad,
              'radius',
              'M',
              Solid,
              Fluid,
              Injection,
              simulProp)

# create a Controller
controller = Controller(Fr, Solid, Fluid, Injection, simulProp)

# run the simulation
controller.run()

# plot results
animate_simulation_results(simulProp.outFileAddress, time_period=0.2)

