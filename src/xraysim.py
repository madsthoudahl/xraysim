#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
Created on Mon Aug 31 14:14:15 2015

@author: Mads Thoudahl

"""
import numpy as np
from xraysimphysics import randomaxisalignedscene, materials
from xraysimgeometry import raygeometry, detectorgeometry

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
                 1, 1, -1, # corner h1 position
                 1, -1, 1, # corner h2 position
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

        for material_id in np.arange(len(materials)):
            materialattenuationfun = materials[material_id][1]
            sceneattenuates += (scenematerials == material_id) * materialattenuationfun(source[energylevel])

        print "attenuation map for source at {0} generated".format(source[0:3])

        # prepare scene geometry
        rayorigin  = source[0:3]

        for detector in detectors:
            # do geometry
            rayunitdirections, raydists = raygeometry(rayorigin, detector[pixelpositions])
            ray_inverse = 1.0 / rayunitdirections
            print "raydirections to detector source at {0} generated".format(detector[pixelpositions][0,0])

            ## TIME for the AABB algorithm...            
            # call function implemented in xraysimgeometry?? 
            tx1 = box.min.x - unitray.x * inv_ray.x
            tx2 = box.max.x - unitray.x * inv_ray.x
            
            # necessary?? read up and think, but explain if left out
            tmin = -inf # array shp as rays
            tmax = +inf # array shp as rays

            tmin = np.min(tx1,tx2)
            tmax = np.min(tx1,tx2)

            # do the same in y and z directions
            tmin = np.min(tmin, np.min(ty1,ty2))
            tmax = np.max(tmax, np.max(ty1,ty2))
            
    print "end of xraysim"

    return rayunitdirections, raydists , detectors, sceneattenuates


if __name__ == '__main__':
    rays, raydists, detectors, sceneatt = xraysim_benchmark()