#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
Created on Mon Aug 31 14:14:15 2015

@author: Mads Thoudahl

"""
import numpy as np
from xraysimphysics import randmaterialAAscene, materials
from xraysimgeometry import coordsAAscene, raygeometry, detectorgeometry

## constants
pi = 3.14 #...

## indexing constants
power = 3
energylevel = 4
pixelpositions = 0
normvec = 1

## simulation variables
scene_size = ss = 2 #is cubed!
res1       = 2
res2       = 2

asource = np.array([-2,0,0, # position
              10,    # Power [Joules]
              80     # Energy level of xrays [keV]
             ])

scenedefs = np.array([0, 0, 0, # position lower left
              1, 1, 1, # position upper right
              scene_size,    # resolution1 (#pixels)
              scene_size,    # resolution2 (#pixels)
              scene_size]    # resolution3 (#pixels)
              # voxelcount = res1 * res2 * res3
            )

detectordef = np.array([ 2,-2,-2, # corner 0 position
                 2, 2, -2, # corner h1 position
                 2, -2, 2, # corner h2 position
                 res1,    # resolution1 (#pixels)
                 res2     # resolution2 (#pixels)
                 # rectangular surface spanned by the 3 points, and the 2 resolutions
               ])



def xraysim_benchmark(
      sourcelist = [asource],
      scenedefinition = scenedefs,
      detectordefs = [detectordef]
      ):

    print "xray simulation executing"
    # generate a random scene from scene definitions
    scenematerials = randmaterialAAscene(scenedefs)
    scenegrid = coordsAAscene(scenedefs)
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
            rayudirs, raydists = raygeometry(rayorigin, detector[pixelpositions])
            ray_inverse = 1.0 / rayudirs
            print "raydirections to detector source at {0} generated".format(detector[pixelpositions][0,0])

            ## TIME for the AABB algorithm...
            # call function implemented in xraysimgeometry??
            ts = np.empty((6,scene_size,scene_size,scene_size,ray_inverse.shape[0]))
            txs, tys, tzs = ts[0:2], ts[2:4], ts[4:6]

            infs = np.array([-np.inf, +np.inf]) # array shp as rays
            tmin, tmax = infs[0], infs[1]

            
            txs[0] = ( scenegrid[0][:-1,:-1,:-1].reshape(ss,ss,ss,1) - rayorigin[0] ) * ray_inverse[:,0]
            txs[1] = ( scenegrid[0][1:,1:,1:].reshape(ss,ss,ss,1) - rayorigin[0] ) * ray_inverse[:,0] 
            #txs.shape should be (scene_size, scene_size, scene_size, res1*res2)??

            tys[0] = ( scenegrid[1][:-1,:-1,:-1].reshape(ss,ss,ss,1) - rayorigin[1] ) * ray_inverse[:,1]
            tys[1] = ( scenegrid[1][1:,1:,1:].reshape(ss,ss,ss,1) - rayorigin[1] ) * ray_inverse[:,1]            

            tzs[0] = ( scenegrid[2][:-1,:-1,:-1].reshape(ss,ss,ss,1) - rayorigin[2] ) * ray_inverse[:,2]
            tzs[1] = ( scenegrid[2][1:,1:,1:].reshape(ss,ss,ss,1) - rayorigin[2] ) * ray_inverse[:,2]            
            
#            #tmin seems to be ok
#            tmin = np.max([np.min(txs, axis=0),np.min(tys, axis=0),np.min(tzs, axis=0)],axis=0)
#            
#            #tmax seems to be flawed somehow            
#            tmax = np.min([np.max(txs, axis=0),np.max(tys, axis=0),np.max(tzs, axis=0)],axis=0)
#            
#            raydst = (tmax-tmin)*(tmax>=tmin) # wrong result! all rays should hit! 
            #raydst.shape should be (scene_size, scene_size, scene_size, res1*res2)?
            
            #raydst is now to be correlated with material/attenuation grid


    print "end of xraysim"

    return rayudirs, raydists , detectors, sceneattenuates, scenegrid, ray_inverse, rayorigin, ts, tmin, tmax


if __name__ == '__main__':
    rayudirs, raydists, detectors, sceneattenuates, scenegrid, ray_inverse, rayorigin, ts, tmin, tmax = xraysim_benchmark()
    
    