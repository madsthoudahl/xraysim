# -*- coding: utf-8 -*-
"""
Created on Mon Aug 31 14:38:12 2015

@author: Mads Thoudahl

"""
import numpy as np

def normalizevector(vector):
    """ returns a vector with frombenius norm = 1 """
    return vector * (1.0 / np.linalg.norm(vector))

def detectorpixelpos(detectordefs):
    """ calculates the pixel positions of the detector, 
        from the detector definitions"""

def raynormdir(src, detpixpos):
    """ calculates an array of normalized (length=1) vectors representing 
        geometric ray directions, from the source and detector 
        pixel positions"""
    return np.norm(detpixpos - src[0:3])
    
    
def detectorgeometry(ddef, npdtype='float32'):
    """ calculates the geometric properties of the detector 
        from its definitions """
    c0 = ddef[0:3] # corner c1-c0-c2
    c1 = ddef[3:6] # corner c0-c1-c3 or 1st axis endposition
    c2 = ddef[6:9] # corner c0-c2-c3 or 2nd axis endposition
    r1 = ddef[9]   # resolution in 1st dimension
    r2 = ddef[10]  # resolution in 2nd dimension
    
    
    ## TO BE CHECKED AND TESTED!
    dshape = (r2,r1)  # CONTROL SEQUENCE
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