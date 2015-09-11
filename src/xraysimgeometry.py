# -*- coding: utf-8 -*-
"""
Created on Mon Aug 31 14:38:12 2015

@author: Mads Thoudahl

"""
import numpy as np



def coordsAAscene(scenedefs, samexyzres=True):
    """ returns a meshgrid of the scene coordinates from its definitions """
    x0, y0, z0 = scenedefs[0:3] # position lower left
    x1, y1, z1 = scenedefs[3:6] # position upper right
    xres, yres, zres = scenedefs[6:9] # resolution

    samexyzres = (xres == yres == zres)

    if not samexyzres:
        xgrid = np.linspace(x0, x1, xres+1)
        ygrid = np.linspace(y0, y1, yres+1)
        zgrid = np.linspace(z0, z1, zres+1)
        return np.meshgrid(xgrid,ygrid,zgrid)

    res = xres + 1
    grid = np.empty((3,res,res,res))
    grid[0] = np.repeat(np.repeat(np.linspace(x0, x1, res).reshape(res,1),res,axis=1).reshape(res,res,1),res,axis=2)
    grid[1] = np.repeat(np.repeat(np.linspace(y0, y1, res).reshape(res,1),res,axis=1).reshape(res,res,1),res,axis=2)
    grid[2] = np.repeat(np.repeat(np.linspace(z0, z1, res).reshape(res,1),res,axis=1).reshape(res,res,1),res,axis=2)
    return grid


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
    rays           = detpixpos.reshape(raycount,3)[:,0:3] - src[0:3]
    raydistances   = np.sqrt(np.sum(rays*rays,1)).reshape(raycount,1)
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
    c0 = ddef[0:3] # corner c1-c0-c2
    c1 = ddef[3:6] # corner c0-c1-c3 or 1st axis endposition
    c2 = ddef[6:9] # corner c0-c2-c3 or 2nd axis endposition
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
    pxs  = np.fromfunction(pcfunx, shape=dshape, dtype=npdtype )
    pys  = np.fromfunction(pcfuny, shape=dshape, dtype=npdtype )
    pzs  = np.fromfunction(pcfunz, shape=dshape, dtype=npdtype )
    pixelpositions = np.dstack((pxs,pys,pzs)) # shape = (r2,r1,3)

#    normalvector = np.cross(di,dj)
#    unitnormalvector = normalizevector(normalvector)
    pixelareavector = np.cross(di, dj)
    result = np.zeros(dshape)

    return  pixelpositions, pixelareavector, dshape, result


def runAABB(scenegrid, rayudirs, rayorigin, rayinverse):
    """ Returns the distances that every ray travels in every voxel """
    ssx, ssy, ssz = np.array(scenegrid.shape[1:4]) - 1
    ts = np.empty((6, ssx, ssy, ssz, rayinverse.shape[0]))
    txs, tys, tzs = ts[0:2], ts[2:4], ts[4:6]


    txs[0] = ( scenegrid[0][:-1,:-1,:-1].reshape(ssx, ssy, ssz, 1) - rayorigin[0] ) * rayinverse[:,0]
    txs[1] = ( scenegrid[0][1:,1:,1:].reshape(ssx, ssy, ssz, 1) - rayorigin[0] ) * rayinverse[:,0]

    tys[0] = ( scenegrid[1][:-1,:-1,:-1].reshape(ssx, ssy, ssz, 1) - rayorigin[1] ) * rayinverse[:,1]
    tys[1] = ( scenegrid[1][1:,1:,1:].reshape(ssx, ssy, ssz, 1) - rayorigin[1] ) * rayinverse[:,1]

    tzs[0] = ( scenegrid[2][:-1,:-1,:-1].reshape(ssx, ssy, ssz, 1) - rayorigin[2] ) * rayinverse[:,2]
    tzs[1] = ( scenegrid[2][1:,1:,1:].reshape(ssx, ssy, ssz, 1) - rayorigin[2] ) * rayinverse[:,2]


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
    tss[::2] = ( scenegrid[:,:-1,:-1,:-1] - rayorigin.reshape(3,1,1,1) ).reshape(3,ss,ss,ss,1) * rayinverse.T.reshape(3,1,1,1,raycount)
    tss[1::2] = ( scenegrid[:,1:,1:,1:] - rayorigin.reshape(3,1,1,1) ).reshape(3,ss,ss,ss,1) * rayinverse[:,:].T.reshape(3,1,1,1,raycount)

    ## Figure incoming and outgoing t_values (t is travelled distance) at each voxel
    t_in = np.max([np.min(txs, axis=0),np.min(tys, axis=0),np.min(tzs, axis=0)],axis=0)
    t_out = np.min([np.max(txs, axis=0),np.max(tys, axis=0),np.max(tzs, axis=0)],axis=0)

    ## Figure distances travelled in each voxel
    voxelraydistances = (t_out-t_in)*(t_out>=t_in)
    return voxelraydistances
