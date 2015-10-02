# -*- coding: utf-8 -*-
"""
Created on Mon Aug 31 14:38:12 2015

@author: Mads Thoudahl

"""
#import bohrium as np
import numpy as np
import numpy


## for enumeration purposes
class Shape:
    cube   = 0
    sphere = 1

class Reference:
    absolute = 0
    relative = 1


def coordsAAscene(scenedefs):
    """ returns a meshgrid of the scene coordinates from its definitions """
    x0, y0, z0 = scenedefs[0:3] # position lower left
    x1, y1, z1 = scenedefs[3:6] # position upper right
    xres, yres, zres = scenedefs[6:9] # resolution

    xgrid = np.array(np.linspace(x0, x1, (xres+1)))
    ygrid = np.array(np.linspace(y0, y1, (yres+1)))
    zgrid = np.array(np.linspace(z0, z1, (zres+1)))
    print len(xgrid), len(ygrid), len(zgrid)
    print xgrid
    print ygrid
    print zgrid
    return np.meshgrid(xgrid,ygrid,zgrid)


def raygeometry(src, detpixpos):
    """ calculates an array of normalized (length=1) vectors representing
        geometric ray directions, from the source and detector
        pixel positions

        parameters:
            src          is source position
            detpixpos    is a 3D array, 2 pixel dimensions,
                         one spatial xyz positions
        returns:
            unitrays     is the ray directions in unit lengths,
                         flattened to 1 pixel dimension, and 1 spatial dimension
            distances    is the distances belonging to the ray-vector
            rayinverses  serves to minimize division operations
        """
    rayshp         = detpixpos.shape
    raycount       = rayshp[0]*rayshp[1]
    origin         = np.array([src[0],src[1],src[2]])
    rays           = detpixpos.reshape(raycount,3)[:,0:3] - origin
    print "rays.shape {}".format(rays.shape)
    print type(rays)
    np.add.reduce(np.array(rays,copy=True))
    raydistances   = np.sqrt( np.add.reduce(rays*rays, axis=1) ).reshape(raycount,1)

    unitrays       = rays*1.0/raydistances
    rayinverses    = 1.0 / unitrays
    return  unitrays, raydistances, rayinverses


def detectorgeometry(ddef, npdtype='float32'):
    """ calculates the geometric properties of the detector
        from its definitions
        parameters:
            ddef is array detector definitions

        returns:
            pixelpositions    The endpoint of all the vectors
            unitnormalvector  of the detector  """
    c0 = np.array([ddef[0],ddef[1],ddef[2]]) # corner c1-c0-c2
    c1 = np.array([ddef[3],ddef[4],ddef[5]]) # corner c0-c1-c3 or 1st axis endposition
    c2 = np.array([ddef[6],ddef[7],ddef[8]]) # corner c0-c2-c3 or 2nd axis endposition
    r1 = ddef[9]   # resolution in 1st dimension
    r2 = ddef[10]  # resolution in 2nd dimension

    dshape = (r2,r1)  # CONTROL of SEQUENCE
    # unit direction vectors of detector sides
    di = (c1 - c0) * (1.0/r1)
    dj = (c2 - c0) * (1.0/r2)
    def pcfun(j,i,k): return  c0[k] + (i+0.5) * di[k] + (j+0.5) * dj[k]
    def pcfunx(j,i):  return pcfun(j,i,0)
    def pcfuny(j,i):  return pcfun(j,i,1)
    def pcfunz(j,i):  return pcfun(j,i,2)
    pxs  = numpy.fromfunction(pcfunx, shape=dshape, dtype=npdtype )
    pys  = numpy.fromfunction(pcfuny, shape=dshape, dtype=npdtype )
    pzs  = numpy.fromfunction(pcfunz, shape=dshape, dtype=npdtype )
    pixelpositions = np.array(numpy.dstack((pxs,pys,pzs))) # shape = (r2,r1,3)

    pixelareavector = np.array(numpy.cross(di, dj))
    result = np.zeros(dshape)

    return  pixelpositions, pixelareavector, dshape, result


def runAABB(scenegrid, rayudirs, rayorigin, rayinverse):
    """ Returns the distances that every ray travels in every voxel """
    ssx, ssy, ssz = (len(scenegrid[0])-1, len(scenegrid[1])-1, len(scenegrid[2])-1)

    ts = np.empty((6, ssx, ssy, ssz, rayinverse.shape[0]))
    txs, tys, tzs = ts[0:2], ts[2:4], ts[4:6]

    txs[0] = ( scenegrid[0][:-1,:-1,:-1].reshape(ssx, ssy, ssz, 1) - rayorigin[0] ) * rayinverse[:,0]
    txs[1] = ( scenegrid[0][1:,1:,1:].reshape(ssx, ssy, ssz, 1) - rayorigin[0] ) * rayinverse[:,0]

    tys[0] = ( scenegrid[1][:-1,:-1,:-1].reshape(ssx, ssy, ssz, 1) - rayorigin[1] ) * rayinverse[:,1]
    tys[1] = ( scenegrid[1][1:,1:,1:].reshape(ssx, ssy, ssz, 1) - rayorigin[1] ) * rayinverse[:,1]

    tzs[0] = ( scenegrid[2][:-1,:-1,:-1].reshape(ssx, ssy, ssz, 1) - rayorigin[2] ) * rayinverse[:,2]
    tzs[1] = ( scenegrid[2][1:,1:,1:].reshape(ssx, ssy, ssz, 1) - rayorigin[2] ) * rayinverse[:,2]

    print "RUNTIME WARNING in t_in = ... when using bohrium"
 #   t_in = numpy.max(numpy.array([numpy.min(numpy.array(txs), axis=0),numpy.min(numpy.array(tys), axis=0),numpy.min(numpy.array(tzs), axis=0)]),axis=0)
    t_in = np.max(np.array([np.min(txs, axis=0),np.min(tys, axis=0),np.min(tzs, axis=0)]),axis=0)
    t_out = np.min(np.array([np.max(txs, axis=0),np.max(tys, axis=0),np.max(tzs, axis=0)]),axis=0)

    voxelraydistances = (t_out-t_in)*(t_out>=t_in)
    return voxelraydistances


def runAABBcompact(scenegrid, rayudirs, rayorigin, rayinverse):
    """ Returns the distances that every ray travels in every voxel"""
    ss, _ssy, _ssz = np.array(scenegrid.shape[1:4]) - 1

    raycount = rayinverse.shape[0]
    ## Initialize & name aliases for calculation array
    tss = np.empty((6, ss, ss, ss, raycount))
    txs, tys, tzs = tss[0:2], tss[2:4], tss[4:6]
#    print "aabbcompactsg = {}, rayorigin: {}".format(scenegrid.shape, rayorigin.shape)
    ## Ray Geometry
    tss[::2] = ( scenegrid[:,:-1,:-1,:-1] - rayorigin.reshape(3,1,1,1) )\
            .reshape(3,ss,ss,ss,1) * rayinverse.T.reshape(3,1,1,1,raycount)

    tss[1::2] = ( scenegrid[:,1:,1:,1:] - rayorigin.reshape(3,1,1,1) )\
            .reshape(3,ss,ss,ss,1) * rayinverse[:,:].T.reshape(3,1,1,1,raycount)

    ## Figure incoming and outgoing t_values (t is travelled distance) at each voxel
    t_in = np.max([np.min(txs, axis=0),np.min(tys, axis=0),np.min(tzs, axis=0)],axis=0)
    t_out = np.min([np.max(txs, axis=0),np.max(tys, axis=0),np.max(tzs, axis=0)],axis=0)

    ## Figure distances travelled in each voxel
    voxelraydistances = (t_out-t_in)*(t_out>=t_in)
    return voxelraydistances



