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
