#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
Created on Mon Aug 31 14:14:15 2015

@author: Mads Thoudahl

"""
from xraysimphysics import randomscenegen, materials
from xraysimgeometry import *

## constants
power = 3
energylevel = 4

## simulation variables
scene_size = 10 #is cubed!

asource = np.array([-1,0,0, # position
              10,    # Power [Joules]
              80     # Energy level of xrays [keV]
             ])

scenedefs = np.array([-0.5,-0.5,-0.5, # position lower left
              0.5, 0.5, 0.5, # position upper right
              scene_size,    # resolution1 (#pixels)
              scene_size,    # resolution2 (#pixels)
              scene_size]    # resolution3 (#pixels)
              # voxelcount = res1 * res2 * res3
            )

detectordef = np.array([ 1,-1,-1, # corner 0 position
                 1, 1, 1, # corner h1 position
                 1, 1, 1, # corner h2 position
                 10,    # resolution1 (#pixels)
                 10     # resolution2 (#pixels)
                 # rectangular surface spanned by the 3 points, and the 2 resolutions
               ])



def xraysim_benchmark(
      sourcelist = [asource],
      scenedefinition = scenedefs,
      detectordefs = [detectordef]
      ):

    # generate a random scene from scene definitions
    scenecorners, scenematerials = randomaxisalignedscene(scenedefs)

    # generate an array of endpoints for rays (at the detector)
    detectors = []
    for ddef in detectordefs:
        detectors.append(detectorgeometry(ddef))

    for source in sourcelist:
        # preprocess the scene physics
        # building a map of attenuation coefficients
        sceneattenuates =  np.zeros(scenematerials.shape)
        for material in materials:
            sceneattenuates += (scenematerials == material) * materials[material](source[energylevel])

        #
        for detector in detectors:
            # do geometry
            rayorigin  = source[0:3]

            raynormdir = # calculate direction of every ray from source to detector and normalize


    print "xraysim execution"

    return 0


if __name__ == '__main__':
    xraysim_benchmark()