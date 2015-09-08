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
eps = 1e-6 #fault tolerance epsilon

## indexing constants
power = 3
energylevel = 4
pixelpositions = 0
normvec = 1

## simulation variables
scene_size = ss = 60 #is cubed!
res1       = 20 # y-dir?
res2       = 20 # z-dir?
#(60**3)*20*20 = memory error :-) @ 8GB

asource = np.array([
              0.0, 0.0, 0.0, # position
              10,    # Power [Joules]
              80     # Energy level of xrays [keV]
             ])

scenedefs = np.array([
              2.0, 0.0, 0.0, # position lower left
              3.0, 1.0, 1.0, # position upper right
              scene_size,    # resolution1 (#pixels)
              scene_size,    # resolution2 (#pixels)
              scene_size]    # resolution3 (#pixels)
              # voxelcount = res1 * res2 * res3
            )

detectordef = np.array([
                 4, 0, 0, # corner 0 position
                 4, 1.4, 0, # corner h1 position
                 4, 0, 1.4, # corner h2 position
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

#            ## AABB algorithm... (AxisAlignedBoundingBox)
#            # call function implemented in xraysimgeometry??
#            ts = np.empty((6,ss,ss,ss,ray_inverse.shape[0]))
#            txs, tys, tzs = ts[0:2], ts[2:4], ts[4:6]
#
#            txs[0] = ( scenegrid[0][:-1,:-1,:-1].reshape(ss,ss,ss,1) - rayorigin[0] ) * ray_inverse[:,0]
#            txs[1] = ( scenegrid[0][1:,1:,1:].reshape(ss,ss,ss,1) - rayorigin[0] ) * ray_inverse[:,0]
#
#            tys[0] = ( scenegrid[1][:-1,:-1,:-1].reshape(ss,ss,ss,1) - rayorigin[1] ) * ray_inverse[:,1]
#            tys[1] = ( scenegrid[1][1:,1:,1:].reshape(ss,ss,ss,1) - rayorigin[1] ) * ray_inverse[:,1]
#
#            tzs[0] = ( scenegrid[2][:-1,:-1,:-1].reshape(ss,ss,ss,1) - rayorigin[2] ) * ray_inverse[:,2]
#            tzs[1] = ( scenegrid[2][1:,1:,1:].reshape(ss,ss,ss,1) - rayorigin[2] ) * ray_inverse[:,2]
#
#
#            t_in = np.max(np.array([np.min(txs, axis=0),np.min(tys, axis=0),np.min(tzs, axis=0)]),axis=0)
#            t_out = np.min(np.array([np.max(txs, axis=0),np.max(tys, axis=0),np.max(tzs, axis=0)]),axis=0)
#
#            raydst1 = (t_out-t_in)*(t_out>=t_in)
#            #raydst.shape should be (ss, ss, ss, res1*res2) :-)


            #COMPACT AABB ALGORITHM
            sgs = scenegrid.shape
            if sgs[1]==sgs[2]==sgs[3]: ##compact is safe

                ## Initialize & name aliases for calculation array
                tss = np.empty((6,ss,ss,ss,ray_inverse.shape[0]))
                txs, tys, tzs = tss[0:2], tss[2:4], tss[4:6]

                ## Ray Geometry
                tss[::2] = ( scenegrid[:,:-1,:-1,:-1] - rayorigin.reshape(3,1,1,1) ).reshape(3,ss,ss,ss,1) * ray_inverse[:,:].T.reshape(3,1,1,1,res1*res2)
                tss[1::2] = ( scenegrid[:,1:,1:,1:] - rayorigin.reshape(3,1,1,1) ).reshape(3,ss,ss,ss,1) * ray_inverse[:,:].T.reshape(3,1,1,1,res1*res2)

                ## Figure incoming and outgoing t_values (t is travelled distance) at each voxel
                t_in = np.max([np.min(txs, axis=0),np.min(tys, axis=0),np.min(tzs, axis=0)],axis=0)
                t_out = np.min([np.max(txs, axis=0),np.max(tys, axis=0),np.max(tzs, axis=0)],axis=0)

                ## Figure distances travelled in each voxel
                raydst = (t_out-t_in)*(t_out>=t_in) ## shape (ss,ss,ss,res1*res2)

#                ## Determine if compact calculation is the same as
#                print "Compact is valid: {}".format(np.max(raydst1-raydst)<eps)



            #raydst is now to be correlated with material/attenuation grid


    print "end of xraysim"

    return rayudirs, raydists , detectors, sceneattenuates, scenegrid, ray_inverse, rayorigin, ts, t_out, t_in, raydst


if __name__ == '__main__':
    rayudirs, raydists, detectors, sceneattenuates, scenegrid, ray_inverse, rayorigin, ts, t_out, t_in, raydst = xraysim_benchmark()
    txs, tys, tzs = ts[0:2], ts[2:4], ts[4:6]


