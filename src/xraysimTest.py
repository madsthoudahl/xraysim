#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
Created on Mon Aug 31 14:14:15 2015

@author: Mads Thoudahl

"""
import time
#import bohrium
import numpy as np
from xraysim import xraysim, buildscene
from xraysimphysics import Material, visualize #emptyAAscene, randomAAscene, matname, attenuation
from xraysimgeometry import Shape, Reference   #coordsAAscene, raygeometry, \
#                            detectorgeometry, runAABB, runAABBcompact



SCENE_SIZE = 20
RES1 = RES2 = 20


asource = np.array([
              0.0, 0.0, 0.0,    # position
              1,                # Power [relative to other sources]
              0.080             # Energy level of xrays [MeV]
             ])

scenedefs = np.array([
              1.0, 0.0, 0.0,    # position lower left
              2.0, 1.0, 1.0,    # position upper right
              SCENE_SIZE,       # resolution1 (#pixels)
              SCENE_SIZE,       # resolution2 (#pixels)
              SCENE_SIZE]       # resolution3 (#pixels)
                                # voxelcount = res1 * res2 * res3
            )

adetector = np.array([
                 2.4, 0, 0,       # corner 0 position
                 2.4, 1.4, 0,     # corner h1 position
                 2.4, 0, 1.4,     # corner h2 position
                 RES1,          # resolution1 (#pixels)
                 RES2           # resolution2 (#pixels)
                                # rectangular surface spanned by the 3 points, and the 2 resolutions
               ])

firstobject = [  Shape.cube,                    # Shape
                 Material.hydrogen,             # Material from which it consists
                 Reference.absolute,            # Reference point of view
                 [[2.2,0.1,0.1],[0.4,0.4,0.4]]  # Shape specific geometric characteristics
              ]

secondobject = [ Shape.cube,                    # Shape
                 Material.titanium,             # Material from which it consists
                 Reference.relative,            # Reference point of view
                 [np.array([1,1,1]),np.array([5,5,5])]  # Shape specific geometric characteristics
               ]



def xraysim_benchmark(
      sourcelist      = [asource],
      detectorlist    = [adetector],
      scenedefinition = scenedefs,
      objectlist      = [firstobject, secondobject]
      ):

    # time building the scene
    start1 = time.time()
    scenegrid, scenematerials = buildscene(scenedefs, objectlist)
    buildtiming = time.time() - start1 

    # time calculating the rays
    start2 = time.time()
    sim = xraysim( sourcelist, detectorlist, scenegrid, scenematerials )
    calctiming = time.time() - start2 # TODO

    # display timings
    print "time taken to build scene of {0:0.0f} voxels containing {1:0.0f} objects:\n\t{2:f} seconds"\
            .format(np.product(scenedefs[6:9]), len(objectlist), buildtiming)
    print "time taken to calculate effect of {0:0.0f} rays:\n\t{1:f} seconds"\
            .format(np.product(scenedefs[6:9])*np.product(detectorlist[0][9:11]), calctiming)
    

    # visualize detectors
    visualize(sim[2][0])

    return sim


if __name__ == '__main__':
    rayudirs, raylengths, detectors, sceneattenuates, scenegrid, rayinverse, rayorigin, raydst, scenemat = xraysim_benchmark()

