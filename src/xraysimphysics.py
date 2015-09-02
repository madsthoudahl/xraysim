# -*- coding: utf-8 -*-
"""
Created on Mon Aug 31 14:37:10 2015

@author: Mads Thoudahl

"""
import numpy as np


materials = {
  # material number: [name, attenuation coefficient equation]
  0:["vacuum"  , lambda e: 0 ],
  1:["hydrogen", lambda e: 1**-e ],
  2:["helium"  , lambda e: 1**-e+1 ],
  3:["iron"    , lambda e: 1**-e+2 ]  #how does this equation really look??
}

def randomaxisalignedscene( scenedefs ):
    """ Generates an axis aligned scene, containing
    axis aligned boxes """

    x0, y0, z0 = scenedefs[0:3] # position lower left
    x1, y1, z1 = scenedefs[3:6] # position upper right
    xs, ys, zs = scenedefs[6:9] # resolution

    matscene = np.random.randint(0, len(materials), xs*ys*zs).reshape((xs,ys,zs))
    cornerpositions = [(x0, y0, z0),(x1, y1, z1)]

    return (cornerpositions, matscene)
