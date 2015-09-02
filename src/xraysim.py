#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
Created on Mon Aug 31 14:14:15 2015

@author: Mads Thoudahl

"""
import numpy as np
from xraysimphysics import randomaxisalignedscene, materials
from xraysimgeometry import raynormdir, detectorgeometry

## constants
pi = 3.14 #...

## indexing constants
power = 3
energylevel = 4
pixelpositions = 0
normvec = 1

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

    print "xray simulation executing"
    # generate a random scene from scene definitions
    scenecorners, scenematerials = randomaxisalignedscene(scenedefs)
    print "axisaligned scene generated"

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
        print "attenuation map for source at {0} generated".format(source[0:3])

        # prepare scene geometry
        rayorigin  = source[0:3]

        for detector in detectors:
            # do geometry
            rayunitdirections = raynormdir(rayorigin, detector[pixelpositions])
        print "raydirections to detector source at {0} generated".format(detector[pixelpositions][0,0])
            
            

    print "end of xraysim"

    return 0


if __name__ == '__main__':
    xraysim_benchmark()