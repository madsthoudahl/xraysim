#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
Created on Mon Aug 31 14:14:15 2015

@author: Mads Thoudahl

"""
#import bohrium as np
import numpy as np
from xraysimphysics import emptyAAscene, addobjtoscene, matname, attenuation#, randomAAscene
from xraysimgeometry import coordsAAscene, raygeometry, \
                            detectorgeometry, runAABB#, runAABBcompact


class Const:
    EPS = 1e-5
    invsqr = 1.0 / (np.pi * 4)

def buildscene(
    scenedefs,
    objlist
    ):
    # generate a random scene from scene definitions
    scenegrid = coordsAAscene(scenedefs)
    scenematerials = emptyAAscene(scenedefs)
    print "axisaligned scene generated"
    for obj in objlist:
        if not addobjtoscene(scenegrid, scenematerials, obj):
            print "object NOT included"
    print "axisaligned scene inhabited with objects"
    return scenegrid, scenematerials



def xraysim(sourcelist,
            detectordeflist,
            scenegrid,
            scenematerials
              ):
    """ performs the calculations figuring out what is detected
        INPUT:
        sourcelist: list of np.array([
                        px,py,pz,       position
                        relative_power, relative to other sources simulated
                        energy])        MeV

        detectorlist: list of np.array([
                        px0,py0,pz0,    position lower left corner
                        px1,py1,pz1,    position upper right corner
                        res1,res2 ])    resolution of detector

        scenegrid:      a numpy array 'meshgrid' shaped (xs+1,ys+1,zs+1)
                        containing absolute coordinates of grid at 'intersections'
                        as returned by buildscene

        scenematerials: a numpy array shaped (xs,ys,zs)
                        containing information of which MATERIAL inhibits this voxel
                        ## for a better model this should also contain
                        ## information of how much of the voxel is filled..
    """
    ## indexing constants
    power, pixelpositions, result =  [3,0,3]

    print "xray simulation executing"
#    sgs = np.array(scenegrid.shape)[1:4] -1 #just xyz dims and voxel# not edge#
    xres, yres, zres = (len(scenegrid[0])-1, len(scenegrid[1])-1, len(scenegrid[2])-1)
    sgs = np.array([xres,yres,zres])
    # generate an array of endpoints for rays (at the detector)
    detectors = []
    for ddef in detectordeflist:
        detectors.append(detectorgeometry(ddef))

    for source in sourcelist:
        _sox, _soy, _soz, spower, senergy = source
        # preprocess the scene physics

        # building a map of attenuation coefficients
        sceneattenuates =  np.zeros(scenematerials.shape)

        for material_id in matname.keys():
            sceneattenuates += (scenematerials == material_id) \
                    * attenuation.get(material_id)(senergy)

        print "attenuation map for source at {0} generated".format(source[0:3])

        # prepare scene geometry
        rayorigin  = source[0:3]

        for pixelpositions, pixelareavector, dshape, result in detectors:
            # do geometry
            rayudirs, raylengths, rayinverse = raygeometry(rayorigin, pixelpositions)
            print "raydirections to detector source at {0} generated"\
                    .format( pixelpositions[0,0] )

            print "runtime warning in runAABB:"
            raydst = runAABB(scenegrid, rayudirs, rayorigin, rayinverse)

            #raydst is now to be correlated with material/attenuation grid
#            print "shapes:"
#            print "pixelareavector:\t{}".format(pixelareavector)
            print sceneattenuates.shape
            print (sgs[0],sgs[1],sgs[2],1)
            print "shape {}".format((sceneattenuates.reshape(sgs[0],sgs[1],sgs[2],1) ).shape)
            dtectattenuates = np.add.reduce( (sceneattenuates.reshape(sgs[0],sgs[1],sgs[2],1)\
                    * raydst).reshape(sgs[0]*sgs[1]*sgs[2],dshape[0],dshape[1]), axis = 0) #shape?
            dtectattenuates = np.add.reduce( (sceneattenuates.reshape(sgs[0],sgs[1],sgs[2],1)\
                    * raydst).reshape(sgs[0]*sgs[1]*sgs[2],dshape[0],dshape[1]), axis = 0) #shape?
#            print "dtectattenuates:\t{}".format(dtectattenuates.shape)

            pixelintensity = (source[power] * Const.invsqr / raylengths).reshape(dshape) #shape?
#            print "pixelintensity:\t{}".format(pixelintensity.shape)

            area = np.dot( rayudirs, pixelareavector.reshape(3,1) ).reshape(dshape) #correct? #shapes
#            print "area:\t\t\t{}".format(area.shape)

            result += pixelintensity * area * np.exp(- dtectattenuates) #shapes


    print "end of xraysim"

    return (rayudirs, raylengths , detectors, sceneattenuates, scenegrid, rayinverse, rayorigin, raydst, scenematerials)

