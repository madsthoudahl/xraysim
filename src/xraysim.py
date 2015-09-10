#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
Created on Mon Aug 31 14:14:15 2015

@author: Mads Thoudahl

"""
import numpy as np
from xraysimphysics import randmaterialAAscene, materials, attenuation_function
from xraysimgeometry import coordsAAscene, raygeometry, \
                            detectorgeometry, runAABB, runAABBcompact

## constants
eps = 1e-6 #fault tolerance epsilon
invsqr = 1.0 / 4 * np.pi

## indexing constants
power = 3
energylevel = 4
pixelpositions = 0
normvec = 1
pixelarea = 2
result = 3

## simulation variables
scene_size = ss = 10 #is cubed!
res1       = 10 # y-dir?
res2       = 10 # z-dir?
#(60**3)*20*20 = memory error :-) @ 8GB

asource = np.array([
              0.0, 0.0, 0.0, # position
              1,     # Power [relative to other sources]
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
    sgs = np.array(scenegrid.shape)[1:4] -1 #just xyz dims and voxel# not edge#
    print "axisaligned scene generated"

    # generate an array of endpoints for rays (at the detector)
    detectors = []
    for ddef in detectordefs:
        detectors.append(detectorgeometry(ddef))


    for source in sourcelist:
        _sox, _soy, _soz, spower, senergy = source
        # preprocess the scene physics

        # building a map of attenuation coefficients
        sceneattenuates =  np.zeros(scenematerials.shape)

        for material_id in np.arange(len(materials)):
            materialattenuationfun = materials[material_id][1]
            sceneattenuates += (scenematerials == material_id) * materialattenuationfun(senergy)

        print "attenuation map for source at {0} generated".format(source[0:3])

        # prepare scene geometry
        rayorigin  = source[0:3]

        for pixelpositions, pixelareavector, dshape, result in detectors:
            # do geometry
            rayudirs, raylengths, rayinverse = raygeometry(rayorigin, pixelpositions)
            print "raydirections to detector source at {0} generated".format( pixelpositions[0,0] )

            if sgs[0]==sgs[1]==sgs[2]: ##compact is safe
                #?COMPACT? AABB ALGORITHM
                raydst = runAABBcompact(scenegrid, rayudirs, rayorigin, rayinverse)
            else:
                #Normal AABB ALGORITHM
                raydst = runAABB(scenegrid, rayudirs, rayorigin, rayinverse)

            #raydst is now to be correlated with material/attenuation grid
#            print "shapes:"
#            print "pixelareavector:\t{}".format(pixelareavector)

            dtectattenuates = np.sum((sceneattenuates.reshape(sgs[0],sgs[1],sgs[2],1)* raydst).reshape(sgs[0]*sgs[1]*sgs[2],dshape[0],dshape[1]), axis = 0) #shape?
#            print "dtectattenuates:\t{}".format(dtectattenuates.shape)

            pixelintensity = (source[power] * invsqr / raylengths).reshape(dshape) #shape?
#            print "pixelintensity:\t{}".format(pixelintensity.shape)

            area = np.dot( rayudirs, pixelareavector.reshape(3,1) ).reshape(dshape) #correct? #shapes
#            print "area:\t\t\t{}".format(area.shape)

            result += pixelintensity * area * np.exp(- dtectattenuates) #shapes


    print "end of xraysim"

    return rayudirs, raydists , detectors, sceneattenuates, scenegrid, ray_inverse, rayorigin, tss, t_out, t_in, raydst


if __name__ == '__main__':
    rayudirs, raydists, detectors, sceneattenuates, scenegrid, ray_inverse, rayorigin, tss, t_out, t_in, raydst = xraysim_benchmark()
    txs, tys, tzs = tss[0:2], tss[2:4], tss[4:6]


