# -*- coding: utf-8 -*-
"""
This file is part of PyFrac.

Created by Haseeb Zia on Tue July 10 2018.
Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, Geo-Energy Laboratory, 2016-2019.
All rights reserved. See the LICENSE.TXT file for more details.
"""

Fig_labels = {
    't': 'Time',
    'time': 'Time',
    'w': 'Fracture Width',
    'width': 'Fracture Width',
    'pf': 'Fluid Pressure',
    'fluid pressure': 'Fluid Pressure',
    'pn': 'net Pressure',
    'net pressure': 'net Pressure',
    'front velocity': 'Front Velocity',
    'v': 'Front Velocity',
    'Reynolds number': 'Reynold\'s number',
    'Re': 'Reynold\'s number',
    'fluid flux': 'Fluid Flux',
    'ff': 'Fluid Flux',
    'fluid velocity': 'Fluid Velocity',
    'fv': 'Fluid Velocity',
    'front_dist_min': 'Closest Distance to Front',
    'd_min': 'Closest Distance to Front',
    'front_dist_max': 'Farthest Distance to Front',
    'd_max': 'Farthest Distance to Front',
    'front_dist_mean': 'Mean Distance to Front',
    'd_mean': 'Mean Distance to Front',
    'V': 'Total Volume',
    'volume': 'Total Volume',
    'lk': 'Leak off',
    'leak off': 'Leak off',
    'lkt': 'Total Leaked of Volume',
    'leaked off total': 'Total Leaked of Volume',
    'ar': 'Aspect Ratio',
    'aspect ratio': 'Aspect Ratio',
    'efficiency': 'Fracture Efficiency',
    'ef': 'Fracture Efficiency',
    'mesh': 'Mesh',
    'footprint': 'Fracture Footprint',
    'surface': 'Fracture Surface',
    'regime': 'Propagation Regime (M-K)',
    'source elements': 'Source Elements',
    'se': 'Source Elements',
    'effective viscosity': 'Effective Viscosity',
    'ev': 'Effective Viscosity'
}

var_labels = {
    't': 'time',
    'time': 'time',
    'w': 'width',
    'width': 'width',
    'pf': 'pressure',
    'fluid pressure': 'pressure',
    'pn': 'pressure',
    'net pressure': 'pressure',
    'front velocity': 'front velocity',
    'v': 'Front Velocity',
    'Reynolds number': 'Reynold\'s number',
    'Re': 'Reynold\'s number',
    'fluid flux': 'fluid flux',
    'ff': 'fluid flux',
    'fluid velocity': 'fluid velocity',
    'fv': 'fluid velocity',
    'front_dist_min': '$R_{min}$',
    'd_min': '$R_{min}$',
    'front_dist_max': '$R_{max}$',
    'd_max': '$R_{max}$',
    'front_dist_mean': '$R_{mean}$',
    'd_mean': '$R_{mean}$',
    'V': 'total volume',
    'volume': 'total volume',
    'lk': 'leak off',
    'leak off': 'leak off',
    'lkt': 'total leaked off volume',
    'leaked off total': 'total leaked off volume',
    'ar': 'aspect ratio',
    'aspect ratio': 'aspect ratio',
    'efficiency': 'fracture efficiency',
    'ef': 'fracture efficiency',
    'mesh': '',
    'footprint': '',
    'surface': '',
    'regime': 'regime (M = 1, K = 0)',
    'source elements': 'source elements',
    'se': 'source elements',
    'effective viscosity': 'effective viscosity',
    'ev': 'effective viscosity'
}

units = {
    't': '($s$)',
    'time': '($s$)',
    'w': ' ($mm$)',
    'width': ' ($mm$)',
    'pf': ' ($MPa$)',
    'fluid pressure': ' ($MPa$)',
    'pn': ' ($MPa$)',
    'net pressure': ' ($MPa$)',
    'front velocity': ' ($m/s$)',
    'v': ' ($m/s$)',
    'Reynolds number': '',
    'Re': '',
    'fluid flux': ' ($m^2/s$)',
    'ff': ' ($m^2/s$)',
    'fluid velocity': ' ($m/s$)',
    'fv': ' ($m/s$)',
    'front_dist_min': ' ($meters$)',
    'd_min': ' ($meters$)',
    'front_dist_max': ' ($meters$)',
    'd_max': ' ($meters$)',
    'front_dist_mean': ' ($meters$)',
    'd_mean': ' ($meters$)',
    'V': ' $m^3$',
    'volume': ' $m^3$',
    'lk': ' $m^3$',
    'leak off': ' $m^3$',
    'lkt': ' $m^3$',
    'leaked off total': ' $m^3$',
    'ar': '',
    'aspect ratio': '',
    'efficiency': '',
    'ef': '',
    'mesh': '',
    'footprint': '',
    'surface': ' ($mm$)',
    'regime': '',
    'source elements': '',
    'se': '',
    'effective viscosity': '($Pa\cdot s$)',
    'ev': '($Pa\cdot s$)'
}

unit_conversion = {
    't': 1,
    'time': 1,
    'w': 1.e-3,
    'width': 1.e-3,
    'pf': 1.e6,
    'fluid pressure': 1.e6,
    'pn': 1.e6,
    'net pressure': 1.e6,
    'front velocity': 1.,
    'v': 1.,
    'Reynolds number': 1.,
    'Re': 1.,
    'fluid flux': 1.,
    'ff': 1.,
    'fluid velocity': 1.,
    'fv': 1.,
    'front_dist_min': 1.,
    'd_min': 1.,
    'front_dist_max': 1.,
    'd_max': 1.,
    'front_dist_mean': 1.,
    'd_mean': 1.,
    'V': 1.,
    'volume': 1.,
    'lk': 1.,
    'leak off': 1.,
    'lkt': 1.,
    'leaked off total': 1.,
    'ar': 1.,
    'aspect ratio': 1.,
    'efficiency': 100.,
    'ef': 100.,
    'mesh': None,
    'footprint': None,
    'surface': 1.e-3,
    'regime': 1.,
    'source elements': 1.,
    'se': 1.,
    'effective viscosity': 1.,
    'ev': 1.
}


supported_variables = ['w', 'width', 'pf', 'fluid pressure', 'pn', 'net pressure',
                       'front velocity', 'v', 'Reynolds number', 'Re', 'fluid flux', 'ff',
                       'fluid velocity', 'fv', 'front_dist_min', 'd_min',
                       'front_dist_max', 'd_max', 'front_dist_mean',
                       'd_mean', 'mesh', 'footprint', 't', 'time', 'volume',
                       'V', 'lk', 'leak off', 'lkt', 'leaked off total',
                       'ar', 'aspect ratio', 'efficiency', 'ef', 'surface', 'front intercepts', 'fi',
                       'regime', 'source elements', 'se', 'effective viscosity', 'ev']

unidimensional_variables = ['time', 't', 'front_dist_min', 'd_min', 'front_dist_max',
                            'd_max', 'V', 'volume', 'front_dist_mean', 'd_mean',
                            'efficiency', 'ef', 'aspect ratio', 'ar', 'lkt', 'leaked off total']
bidimensional_variables = ['w', 'width', 'pf', 'fluid pressure', 'pn', 'net pressure',
                           'front velocity', 'v', 'Reynolds number', 'Re', 'fluid flux', 'ff',
                           'fluid velocity', 'fv', 'lk', 'leak off', 'surface', 'front intercepts', 'fi',
                           'regime', 'effective viscosity', 'ev']

required_string = {
    't': '100000',
    'time': '100000',
    'w': '000100',
    'width': '000100',
    'pn': '001000',
    'net pressure': '001000',
    'front velocity': '000010',
    'v': '000010',
    'front_dist_min': '010000',
    'd_min': '010000',
    'front_dist_max': '010000',
    'd_max': '010000',
    'front_dist_mean': '010000',
    'd_mean': '010000',
    'radius': '010000',
    'r': '010000'
}

err_msg_variable = 'Given variable is not supported. Select one of the following:\n' \
                    '-- \'w\' or \'width\'\n' \
                    '-- \'pf\' or \'fluid pressure\'\n' \
                    '-- \'pn\' or \'net pressure\'\n' \
                    '-- \'Re\' or \'Reynolds number\'\n' \
                    '-- \'v\' or \'front velocity\'\n' \
                    '-- \'ff\' or \'fluid flux\'\n' \
                    '-- \'fv\' or \'fluid velocity\'\n' \
                    '-- \'d_min\' or \'front_dist_min\'\n' \
                    '-- \'d_max\' or \'front_dist_max\'\n' \
                    '-- \'d_mean\' or \'front_dist_mean\'\n' \
                    '-- \'mesh\'\n' \
                    '-- \'footprint\'\n' \
                    '-- \'V\' or \'volume\'\n' \
                    '-- \'lk\' or \'leak off\'\n' \
                    '-- \'lkt\' or \'leaked off total\'\n' \
                    '-- \'ar\' or \'aspect ratio\'\n' \
                    '-- \'ef\' or \'efficiency\'\n' \
                    '-- \'surface\'\n' \
                    '-- \'regime\'\n' \
                    '-- \'se\' or \'source elements\'' \
                    '-- \'ev\' or \'effective viscosity\''

supported_projections ={
    'w': ['2D_clrmap', '2D_contours', '3D'],
    'width': ['2D_clrmap', '2D_contours', '3D'],
    'pf': ['2D_clrmap', '2D_contours', '3D'],
    'fluid pressure': ['2D_clrmap', '2D_contours', '3D'],
    'pn': ['2D_clrmap', '2D_contours', '3D'],
    'net pressure': ['2D_clrmap', '2D_contours', '3D'],
    'front velocity': ['2D_clrmap', '2D_contours'],
    'v': ['2D_clrmap', '2D_contours'],
    'Reynolds number': ['2D_clrmap', '2D_contours', '3D'],
    'Re': ['2D_clrmap', '2D_contours', '3D'],
    'fluid flux': ['2D_clrmap', '2D_contours', '3D'],
    'ff': ['2D_clrmap', '2D_contours', '3D'],
    'fluid velocity': ['2D_clrmap', '2D_contours', '3D'],
    'fv': ['2D_clrmap', '2D_contours', '3D'],
    'front_dist_min': ['1D'],
    'd_min': ['1D'],
    'front_dist_max': ['1D'],
    'd_max': ['1D'],
    'front_dist_mean': ['1D'],
    'd_mean': ['1D'],
    'mesh': ['2D', '3D'],
    'footprint': ['2D', '3D'],
    't': ['1D'],
    'time': ['1D'],
    'volume': ['1D'],
    'V': ['1D'],
    'lk': ['2D_clrmap', '2D_contours', '3D'],
    'leak off': ['2D_clrmap', '2D_contours', '3D'],
    'lkt': ['1D'],
    'leaked off total': ['1D'],
    'ar': ['1D'],
    'aspect ratio': ['1D'],
    'efficiency': ['1D'],
    'ef': ['1D'],
    'surface': ['3D'],
    'regime': ['2D_clrmap', '2D_contours', '3D'],
    'source elements': ['2D_clrmap', '2D_contours', '3D'],
    'se': ['2D_clrmap', '2D_contours', '3D'],
    'effective viscosity': ['2D_clrmap', '2D_contours', '3D'],
    'ev': ['2D_clrmap', '2D_contours', '3D']
}

err_var_not_saved = "The required variable is not available. Probably, saving of the variable was not\n" \
                    "enabled during the simulation. Enable saving it through simulation properties."

TS_errorMessages = ["Propagation not attempted!",
                     "Time step successful!",
                     "Evaluated level set is not valid!",
                     "Front is not tracked correctly!",
                     "Evaluated tip volume is not valid!",
                     "Solution obtained from the elastohydrodynamic solver is not valid!",
                     "Did not converge after max iterations!",
                     "Tip inversion is not correct!",
                     "Ribbon element not found in the enclosure of the tip cell!",
                     "Filling fraction not correct!",
                     "Toughness iteration did not converge!",
                     "projection could not be found!",
                     "Reached end of grid!",
                     "Leak off can't be evaluated!",
                     "fracture fully closed"
                    ]

