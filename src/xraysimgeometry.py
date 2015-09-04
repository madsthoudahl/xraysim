# -*- coding: utf-8 -*-
"""
Created on Mon Aug 31 14:38:12 2015

@author: Mads Thoudahl

"""
import numpy as np

def normalizevector(vector):
    """ returns a vector with frombenius norm = 1 """
    return vector * (1.0 / np.linalg.norm(vector))


def coordsAAscene(scenedefs):
    """ returns a meshgrid of the scene coordinates from its definitions """
    x0, y0, z0 = scenedefs[0:3] # position lower left
    x1, y1, z1 = scenedefs[3:6] # position upper right
    xres, yres, zres = scenedefs[6:9] # resolution

    xgrid = np.linspace(x0, x1, xres+1)
    ygrid = np.linspace(y0, y1, yres+1)
    zgrid = np.linspace(z0, z1, zres+1)
    return np.meshgrid(xgrid,ygrid,zgrid)#, sparse=True, indexing='xy')
#    return [xgrid,ygrid,zgrid]


def raygeometry(src, detpixpos):
    """ calculates an array of normalized (length=1) vectors representing
        geometric ray directions, from the source and detector
        pixel positions

        parameters:
            src        is source position
            detpixpos  is a 3D array, 2 pixel dimensions,
                       one spatial xyz positions
        returns:
            unitrays   is the ray directions in unit lengths,
                       flattened to 1 pixel dimension, and 1 spatial dimension
            distances  is the distances belonging to the ray-vector
        """
    rayshp = detpixpos.shape
    raycount = rayshp[0]*rayshp[1]
    rays = detpixpos.reshape(raycount,3)[:,0:3] - src[0:3]
    dists = np.sqrt(np.sum(rays*rays,1)).reshape(raycount,1)
    return rays * 1.0/dists, dists


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

    normalvector = np.cross(di,dj)
    unitnormalvector = normalizevector(normalvector)

    return (pixelpositions, unitnormalvector)