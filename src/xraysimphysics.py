# -*- coding: utf-8 -*-
"""
Created on Mon Aug 31 14:37:10 2015

@author: Mads Thoudahl

"""
import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from xraysimgeometry import Shape, Reference

class Material:
    vacuum   = 0
    hydrogen = 1
    titanium = 22
    

matname = {
    Material.vacuum   : 'vacuum',
    Material.hydrogen : 'hydrogen',
    Material.titanium : 'titanium'
}

# Data from NIST database
# Hardcoded energies at which attenuation coefficients were measured in [MeV]

# Hardcoded measured densities and attenuation coefficients.
# Densities were measured as  rho [g/cm3], and
# Attenuation values were measured in   mu / rho [cm2/g]
# thus Attenuation coefficients returned by spline functions has unit  [1/cm]

hydrogen_energies = np.array(
    [0.001,   0.0015,  0.002,   0.003,   0.004,   0.005,   0.006,  \
     0.008,   0.01,    0.015,   0.02,    0.03,    0.04,    0.05,   \
     0.06,    0.08,    0.1,     0.15,    0.2,     0.3,     0.4,    \
     0.5,     0.6,     0.8,     1.0,     1.25,    1.5,     2.0,    \
     3.0,     4.0,     5.0,     6.0,     8.0,    10.0,    15.0,    \
    20.0    ])
hydrogen_density = 0.00008375
hydrogen_mu_per_rhos = np.array(
    [7.217,   2.148,   1.059,   0.5612,  0.4546,  0.4193,  0.4042, \
     0.3914,  0.3854,  0.3764,  0.3695,  0.357,   0.3458,  0.3355, \
     0.326,   0.3091,  0.2944,  0.2651,  0.2429,  0.2112,  0.1893, \
     0.1729,  0.1599,  0.1405,  0.1263,  0.1129,  0.1027,  0.08769,\
     0.06921, 0.05806, 0.05049, 0.04498, 0.03746, 0.03254, 0.02539,\
     0.02153 ])

titanium_energies = np.array(
    [0.001,    0.0015,  0.002,  0.003,  0.0047,  0.0049664, 0.0049664,\
     0.005,    0.006,   0.008,  0.01,   0.015,   0.02,      0.03,     \
     0.04,     0.05,    0.06,   0.08,   0.1,     0.15,      0.2,      \
     0.3,      0.4,     0.5,    0.6,    0.8,     1,         1.25,     \
     1.5,      2,       3,      4,      5,       6,         8,        \
     10,      15,      20   ])
titanium_density = 4540
titanium_mu_per_rhos = np.array(
    [5869,    2096,     986,     332.3,   151.7,    83.8,   687.8,    \
      683.8,   432.3,   202.3,   110.7,    35.87,   15.85,    4.972,  \
        2.214,   1.213,   0.7661,  0.4052,  0.2721,  0.1649,  0.1314, \
        0.1043,  0.09081, 0.08191, 0.07529, 0.06572, 0.05891, 0.05263,\
        0.04801, 0.0418,  0.03512, 0.03173, 0.02982, 0.02868, 0.02759,\
        0.02727, 0.02762, 0.02844 ])


#    """ Hashmap holding energy attenuation functions as cubic spline
#        interpolations functions for various materials """
attenuation = {
    Material.vacuum   : interp1d(hydrogen_energies, \
            np.zeros(len(hydrogen_energies)), kind='cubic'),
    Material.hydrogen : interp1d(hydrogen_energies, \
            hydrogen_mu_per_rhos * hydrogen_density, kind='cubic'),
    Material.titanium : interp1d(titanium_energies, \
            titanium_mu_per_rhos * titanium_density, kind='cubic'),
}

# example interaction with energy at 220 keV = 0.22 MeV
# attenuation_function.get(22)(0.232)



def randomAAscene( scenedefs ):
    """ Generates an axis aligned scene, containing random materials
        axis aligned boxes """
    xs, ys, zs = scenedefs[6:9] # resolution
    scene = emptyAAscene(scenedefs)
    for i in np.arange(0, np.random.randint(1,5)):
        corner = cx, cy, cz = np.array([
                            np.random.randint(0,xs),
                            np.random.randint(0,ys),
                            np.random.randint(0,zs)])

        size   = np.array([ np.random.randint(0,xs-cx),
                            np.random.randint(0,ys-cy),
                            np.random.randint(0,zs-cz)])

        material = np.random.choice(attenuation.keys())
        
        cubedescribers = (corner, size)

        if not addAAcube(scenedefs, scene, cubedescribers, material, ref=Reference.relative):
            print "randmaterialAAscene: adding a aa-cube failed"
    return scene



def emptyAAscene( scenedefs ):
    """ Generates an axis aligned scene, containing vacuum in all
        axis aligned boxes
        Model: all boxes are *fully* occupied by exactly one (1)
               material, vacuum (Material.vacuum=0) represents zero attenuation"""
    shp = xs, ys, zs = scenedefs[6:9] # resolution
    return np.zeros(shp)


def addobjtoscene(scene, matscene, obj):
    shape           = obj[0]
    material        = obj[1]
    reference       = obj[2]
    shapedescribers = obj[3]
    if shape == Shape.cube:
        return addAAcube(scene, matscene, shapedescribers, material, reference)
    else:
        print "Cube shaped objects are the only ones supported at this time"
        print "Object dropped"
    


def addAAcube(scene, matscene, describers, material, ref=Reference.relative):
#        addAAcube(scene, matscene, shapedescribers, material, reference)
    """ adds a cubic shape to scene, replacing any material in same location

        scene:         scene metadata p0, p1, res
        matscene:      initiated (xsize, y_size, z_size) scene in resolution
        cubedescribers:
            mincorner: reference dependant starting point
            size:      reference dependant side length (3,) np array
        material:      Integer from Material Class (enumerator)
        ref:           input arguments refer to the following
                       Reference.relative    resolution relative references
                       Reference.absolute    real world coordinate references

       returns:        Bool - success
                       changes materialscene as sideeffect
       """
    mincorner, size = describers
    if ref == Reference.relative:
        x0,y0,z0 = mincorner
        x1,y1,z1 = maxcorner = mincorner + size
        sshp = matscene.shape
        if np.any(sshp < maxcorner):
            print "WARNING, PART OF ADDED CUBE OUTSIDE SCENE"
            # reduce size to match scene representation
            difs = maxcorner - sshp
            size -= (difs>0) * difs
        matscene[x0:x1,y0:y1,z0:z1] = np.ones(size) * material
        return True

    if ref == Reference.absolute:
        return False

    return False



### old visualization function to be modified
def visualize(detector):
    """ Visualize a detector output """
    _pixelpositions, _pixelareavector, _dshape, result = detector
    title = "Detector image"
    plt.figure()
    plt.title(title)
    plt.imshow(result, cmap='gray', origin='lower')
    plt.colorbar(shrink=.92)
    plt.xticks([])
    plt.yticks([])



def visualize2(scene, outputdir, scenefile, outputarg, conf):
    """ Visualizes a simulation """
    fig  = 0
    afig = 100
    energies = [str(source.energy) for source in scene.sources]
    srcstr   = [str(source) for source in scene.sources]
    for det in scene.detectors:
        fig  += 1
        title = "Detector image: {}\n{}\n{}"
        path =  outputdir + scenefile + '/' + outputarg # savepath
        res = '@R' + str(det.a_res) +'x'+ str(det.b_res) # resolution
        nrg = '-E' + '-E'.join(energies)
        mat = conf.MATTER if conf.MATTER else ''
        fmt = "." + conf.OUTPUTFORMAT                   # format
        fname =  path + res + nrg + mat + fmt

        plt.figure(fig)
        plt.title(title .format(scenefile, str(det), '\n'.join(srcstr) ) )
        plt.imshow(det.result, cmap='gray', origin='lower')
        plt.colorbar(shrink=.92)
        plt.xticks([])
        plt.yticks([])
        plt.show()

        if conf.ATTENUATION:
            afig  += 1
            atitle = "X-ray attenuation image: {}\n{}\n{}"
            path = outputdir + scenefile + '/' + outputarg
            res = '@R' + str(det.a_res) +'x'+ str(det.b_res)
            nrg = "-E" + '-E'.join(energies)
            mat = conf.MATTER if conf.MATTER else ''
            fmt = "-A." + conf.OUTPUTFORMAT
            afname = path + res + nrg + mat + fmt
            plt.figure(afig)
            plt.title(atitle .format(scenefile, str(det), '\n'.join(srcstr) ) )
            plt.imshow(det.attenuation, cmap='gray', origin='lower')
            plt.colorbar(shrink=.92)
            plt.xticks([])
            plt.yticks([])

        try: # Store it on disk
            plt.figure(fig)
            plt.savefig( fname, format=conf.OUTPUTFORMAT )
            if conf.ATTENUATION:
                plt.figure(afig)
                plt.savefig( afname, format=conf.OUTPUTFORMAT )
        except IOError:
            try:
                os.makedirs( outputdir + scenefile )
                plt.savefig( fname, format=conf.OUTPUTFORMAT )
            except OSError:
                if not os.path.isdir(outputdir + scenefile):
                    raise
        if conf.VISUALIZE:
            plt.show()
