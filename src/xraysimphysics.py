# -*- coding: utf-8 -*-
"""
Created on Mon Aug 31 14:37:10 2015

@author: Mads Thoudahl

"""
import numpy as np
from scipy.interpolate import interp1d

# Data from NIST database
# Hardcoded energies at which attenuation coefficients were measured in [MeV]

# Hardcoded measured densities and attenuation coefficients.
# Densities were measured as  rho [g/cm3], and
# Attenuation values were measured in   mu / rho [cm2/g]
# thus Attenuation coefficients returned by spline functions has unit  [1/cm]

hydrogen_energies = np.array(
    [0.001,   0.0015,  0.002,   0.003,   0.004,   0.005,   0.006,  \
     0.008,   0.01,    0.015,   0.02,    0.03,    0.04,    0.05,   \
     0.06,    0.08,    0.1,     0.15,    0.2,     0.3,     0.4,    \
     0.5,     0.6,     0.8,     1.0,     1.25,    1.5,     2.0,    \
     3.0,     4.0,     5.0,     6.0,     8.0,    10.0,    15.0,    \
    20.0    ])
hydrogen_density = 0.00008375
hydrogen_mu_per_rhos = np.array(
    [7.217,   2.148,   1.059,   0.5612,  0.4546,  0.4193,  0.4042, \
     0.3914,  0.3854,  0.3764,  0.3695,  0.357,   0.3458,  0.3355, \
     0.326,   0.3091,  0.2944,  0.2651,  0.2429,  0.2112,  0.1893, \
     0.1729,  0.1599,  0.1405,  0.1263,  0.1129,  0.1027,  0.08769,\
     0.06921, 0.05806, 0.05049, 0.04498, 0.03746, 0.03254, 0.02539,\
     0.02153 ])

titanium_energies = np.array(
    [0.001,    0.0015,  0.002,  0.003,  0.0047,  0.0049664, 0.0049664,\
     0.005,    0.006,   0.008,  0.01,   0.015,   0.02,      0.03,     \
     0.04,     0.05,    0.06,   0.08,   0.1,     0.15,      0.2,      \
     0.3,      0.4,     0.5,    0.6,    0.8,     1,         1.25,     \
     1.5,      2,       3,      4,      5,       6,         8,        \
     10,      15,      20   ])
titanium_density = 4540
titanium_mu_per_rhos = np.array(
    [5869,    2096,     986,     332.3,   151.7,    83.8,   687.8,    \
      683.8,   432.3,   202.3,   110.7,    35.87,   15.85,    4.972,  \
        2.214,   1.213,   0.7661,  0.4052,  0.2721,  0.1649,  0.1314, \
        0.1043,  0.09081, 0.08191, 0.07529, 0.06572, 0.05891, 0.05263,\
        0.04801, 0.0418,  0.03512, 0.03173, 0.02982, 0.02868, 0.02759,\
        0.02727, 0.02762, 0.02844 ])

attenuation_function = {
#    """ Hashmap holding energy attenuation functions as cubic spline
#        interpolations functions for various materials """
    0          : interp1d(hydrogen_energies, np.zeros(len(hydrogen_energies)), kind='cubic'),
    'vacuum'   : interp1d(hydrogen_energies, np.zeros(len(hydrogen_energies)), kind='cubic'),
    1          : interp1d(hydrogen_energies, hydrogen_mu_per_rhos * hydrogen_density, kind='cubic'),
    'hydrogen' : interp1d(hydrogen_energies, hydrogen_mu_per_rhos * hydrogen_density, kind='cubic'),
    22         : interp1d(titanium_energies, titanium_mu_per_rhos * titanium_density, kind='cubic'),
    'titanium' : interp1d(titanium_energies, titanium_mu_per_rhos * titanium_density, kind='cubic')
}
# example interaction with energy at 220 keV = 0.22 MeV
# attenuation_function.get('titanium')(0.22)

materials = {
  0: {'name':"vacuum",
      'Z':0,
      'mu_at_E': lambda e: 0*e },

  1: {'name':"hydrogen",
      'Z':1,
      'density':0.00008375, # (g/cm3)
      'mu_at_E': interp1d(hydrogen_energies, hydrogen_mu_per_rhos, kind='cubic') } # (1/cm)
}

def randmaterialAAscene( scenedefs ):
    """ Generates an axis aligned scene, containing random materials
        axis aligned boxes """
    xs, ys, zs = scenedefs[6:9] # resolution

    return np.random.randint(0, len(materials), xs*ys*zs).reshape((xs,ys,zs))


def emptyAAscene( scenedefs ):
    """ Generates an axis aligned scene, containing vacuum in all
        axis aligned boxes """
    xs, ys, zs = scenedefs[6:9] # resolution

    return np.zeros((xs,ys,zs))

