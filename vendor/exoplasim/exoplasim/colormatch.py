import numpy as np
from scipy import interpolate


cie_wvl,cie_xx,cie_yy,cie_zz = np.loadtxt("/".join(__file__.split("/")[:-1])+"/cmf.csv",unpack=True,delimiter=',')


def _loadcie():
    wvl,xx,yy,zz = np.loadtxt("/".join(__file__.split("/")[:-1])+"/cmf.csv",unpack=True,delimiter=',')
    return wvl,xx,yy,zz
    
def interpolant(wvl):
    '''Construct interpolated color-matching functions for a series of wavelengths.
    
    Parameters
    ----------
    wvl : array-like
        Array of N wavelengths in nanometers
    
    Returns
    -------
    (array-like,array-like,array-like)
        (fx(lambda),fy(lambda),fz(lambda))
    '''
    #w0,xx,yy,zz = _loadcie()
    w0 = cie_wvl
    xx = cie_xx
    yy = cie_yy
    zz = cie_zz
    fx = interpolate.interp1d(w0,xx)
    fy = interpolate.interp1d(w0,yy)
    fz = interpolate.interp1d(w0,zz)
    imin = np.where(wvl>np.amin(w0))[0][0]
    imax = np.where(wvl<np.amax(w0))[0][-1]
    wn = wvl[imin:imax+1]
    xn = fx(wn)
    yn = fy(wn)
    zn = fz(wn)
    return (xn,yn,zn)
    
def makexyz(wvl,spec,interpolant=None):
    '''Convert a spectrum to XYZ colour coordinates.
    
    The XYZ colour coordinate system is related to how the human eye's colour-receptive cells respond to
    light. The XYZ coordinates are computed by convolving the spectrum with three different empirically-derived
    response functions. These coordinates can then be transformed into RGB colour coordinates. Note that in this 
    system, (x,y,z) are brightness-normalized, while (X,Y,Z) are not. Additionally, z=1-x-y. Therefore, the
    three coordinates that are needed to produce an RGB tuple are not (x,y,z), but (x,y,Y). 
    
    Parameters
    ----------
    wvl : array-like
        Wavelengths in nanometers, of shape (N,)
    spec : array-like
        Spectrum, of shape (N,); units are arbitrary, but it should be given in flux, not flux density.
    interpolant : array-like, optional
        2D array of shape (3,N), corresponding to interpolated color-matching functions
        
    Returns
    -------
    numpy.ndarray, numpy.ndarray, numpy.ndarray
        x,y,Y--z can be inferred from x and y (z = 1-x-y), but Y preserves intensity information.
    '''
    
    if np.amin(wvl)<1.0e-3: #probably meters not nanometers
        wvl*=1.0e9
    #w0,xx,yy,zz = _loadcie()
    w0 = cie_wvl
    xx = cie_xx
    yy = cie_yy
    zz = cie_zz
    imin = np.where(wvl>np.amin(w0))[0][0]
    imax = np.where(wvl<np.amax(w0))[0][-1]
    wn = wvl[imin:imax+1]
    specn = spec[imin:imax+1]
    if interpolant is None:
        fx = interpolate.interp1d(w0,xx)
        fy = interpolate.interp1d(w0,yy)
        fz = interpolate.interp1d(w0,zz)
        
        xn = fx(wn)
        yn = fy(wn)
        zn = fz(wn)
    else:
        xn = interpolant[0]
        yn = interpolant[1]
        zn = interpolant[2]
    
    XI = np.trapz(xn[~np.isnan(specn)]*specn[~np.isnan(specn)],x=wn[~np.isnan(specn)])
    YI = np.trapz(yn[~np.isnan(specn)]*specn[~np.isnan(specn)],x=wn[~np.isnan(specn)])
    ZI = np.trapz(zn[~np.isnan(specn)]*specn[~np.isnan(specn)],x=wn[~np.isnan(specn)])
    xyzmin = np.amin((XI,YI,ZI))
    if xyzmin<0:
        XI -=xyzmin
        YI -=xyzmin
        ZI -=xyzmin
    if (XI+YI+ZI)>0:
        xnu = XI/(XI+YI+ZI)
        ynu = YI/(XI+YI+ZI)
        znu = 1.0-(xnu+ynu)
    else:
        xnu=0
        ynu=0
        znu=0
    
    return xnu,ynu,YI

def xyz2rgb(x,y,normalization,gamut="sRGB"):
    '''Convert (x,y) coordinates to RGB tuples, normalized to a given value.
    
    Note that z=1-x-y. This routine uses a wide gamut colourspace found at http://www.brucelindbloom.com/index.html?Eqn_RGB_XYZ_Matrix.html.
    
    Parameters
    ----------
    x : array-like
        x colour-coordinate
    y : array-like
        y colour-coordinate.
    normalization : array-like or float
        Normalization factor for scaling RGB values
    gamut : str or np.ndarray(3,3), optional
        Color gamut to be used. For available built-in color gamuts, see colormatch.colorgamuts.
    
    Returns
    -------
    array-like, array-like, array-like
        R,G,B colour values.
    '''
    
    if gamut in colorgamuts:
        colorgamut = colorgamuts[gamut]
        if "%s_norm"%gamut in colorgamuts:
            extranorm = colorgamuts["%s_norm"%gamut]
        else:
            extranorm = 1.0
    else:
        if gamut.shape==np.array((3,3)).shape:
            colorgamut = gamut
            extranorm = 1.0
        else:
            raise Exception("Error: must specify valid RGB color gamut")
    
    z = 1-(x+y)
    
    r = np.sum(colorgamut[0,:]*x)*normalization
    g = np.sum(colorgamut[1,:]*y)*normalization
    b = np.sum(colorgamut[2,:]*z)*normalization
    cmax = np.amax((r,g,b))
    
    #return r/cmax,g/cmax,b/cmax
    return r,g,b


#From http://www.brucelindbloom.com/index.html?Eqn_RGB_XYZ_Matrix.html
#Assume whitepoint is D50 unless indicated otherwise
colorgamuts = {"wide": np.array([[ 1.4628067, -0.1840623, -0.2743606],
                                 [-0.5217933,  1.4472381,  0.0677227],
                                 [ 0.0349342, -0.0968930,  1.2884099]]),
               "sRGB_D65": np.array([[ 3.2404542, -1.5371385, -0.4985314],
                                     [-0.9692660,  1.8760108,  0.0415560],
                                     [ 0.0556434, -0.2040259,  1.0572252]]),
               "sRGB": np.array([[ 3.1338561, -1.6168667, -0.4906146],
                                 [-0.9787684,  1.9161415,  0.0334540],
                                 [ 0.0719453, -0.2289914,  1.4052427]]),
               "CIE_E": np.array([[ 2.3706743, -0.9000405, -0.4706338],
                                  [-0.5138850,  1.4253036,  0.0885814],
                                  [ 0.0052982, -0.0146949,  1.0093968]]),
               "NTSC_C": np.array([[ 1.9099961, -0.5324542, -0.2882091],
                                   [-0.9846663,  1.9991710, -0.0283082],
                                   [ 0.0583056, -0.1183781,  0.8975535]]),
               "ProPhoto": np.array([[ 1.3459433, -0.2556075, -0.0511118],
                                     [-0.5445989,  1.5081673,  0.0205351],
                                     [ 0.0000000,  0.0000000,  1.2118128]]),
               "CIE": np.array([[ 2.3638081, -0.8676030, -0.4988161],
                                [-0.5005940,  1.3962369,  0.1047562],
                                [ 0.0141712, -0.0306400,  1.2323842]]),
               "Apple": np.array([[ 2.8510695, -1.3605261, -0.4708281],
                                  [-1.0927680,  2.0348871,  0.0227598],
                                  [ 0.1027403, -0.2964984,  1.4510659]]),
               "Adobe_D65": np.array([[ 2.0413690, -0.5649464, -0.3446944],
                                      [-0.9692660,  1.8760108,  0.0415560],
                                      [ 0.0134474, -0.1183897,  1.0154096]]),
               "Apple_D65": np.array([[ 2.9515373, -1.2894116, -0.4738445],
                                      [-1.0851093,  1.9908566,  0.0372026],
                                      [ 0.0854934, -0.2694964,  1.0912975]]),
               "Best": np.array([[ 1.6832270, -0.4282363, -0.2360185],
                                 [-0.7710229,  1.7065571,  0.0446900],
                                 [ 0.0400013, -0.0885376,  1.2723640]]),
               "Bruce_D65": np.array([[ 2.7454669, -1.1358136, -0.4350269],
                                      [-0.9692660,  1.8760108,  0.0415560],
                                      [ 0.0112723, -0.1139754,  1.0132541]]),
               "ColorMatch": np.array([[ 2.6422874, -1.2234270, -0.3930143],
                                       [-1.1119763,  2.0590183,  0.0159614],
                                       [ 0.0821699, -0.2807254,  1.4559877]]),
               "Don": np.array([[ 1.7603902, -0.4881198, -0.2536126],
                                [-0.7126288,  1.6527432,  0.0416715],
                                [ 0.0078207, -0.0347411,  1.2447743]]),
               "ECI": np.array([[ 1.7827618, -0.4969847, -0.2690101],
                                [-0.9593623,  1.9477962, -0.0275807],
                                [ 0.0859317, -0.1744674,  1.3228273]]),
               "Ekta-Space-PS5": np.array([[ 2.0043819, -0.7304844, -0.2450052],
                                           [-0.7110285,  1.6202126,  0.0792227],
                                           [ 0.0381263, -0.0868780,  1.2725438]]),
               "PAL/SECAM_D65": np.array([[ 3.0628971, -1.3931791, -0.4757517],
                                          [-0.9692660,  1.8760108,  0.0415560],
                                          [ 0.0678775, -0.2288548,  1.0693490]]),
               "SMPTE-C_D65": np.array([[ 3.5053960, -1.7394894, -0.5439640],
                                        [-1.0690722,  1.9778245,  0.0351722],
                                        [ 0.0563200, -0.1970226,  1.0502026]]),
               "Adobe": np.array([[ 1.9624274, -0.6105343, -0.3413404],
                                  [-0.9787684,  1.9161415,  0.0334540],
                                  [ 0.0286869, -0.1406752,  1.3487655]]),
               "Bruce": np.array([[ 2.6502856, -1.2014485, -0.4289936],
                                  [-0.9787684,  1.9161415,  0.0334540],
                                  [ 0.0264570, -0.1361227,  1.3458542]]),
               "NTSC": np.array([[ 1.8464881, -0.5521299, -0.2766458],
                                 [-0.9826630,  2.0044755, -0.0690396],
                                 [ 0.0736477, -0.1453020,  1.3018376]]),
               "PAL/SECAM": np.array([[ 2.9603944, -1.4678519, -0.4685105],
                                      [-0.9787684,  1.9161415,  0.0334540],
                                      [ 0.0844874, -0.2545973,  1.4216174]]),
               "SMPTE-C": np.array([[ 3.3921940, -1.8264027, -0.5385522],
                                    [-1.0770996,  2.0213975,  0.0207989],
                                    [ 0.0723073, -0.2217902,  1.3960932]])
               }
               
#XYZ components for each whitepoint illuminant, again from Bruce Lindbloom
illuminants = {"A"   :np.array([1.09850,1.00000,0.35585]),
               "B"   :np.array([0.99072,1.00000,0.85223]),
               "C"   :np.array([0.98074,1.00000,1.18232]),
               "D50" :np.array([0.96422,1.00000,0.82521]),
               "D55" :np.array([0.95682,1.00000,0.92149]),
               "D65" :np.array([0.95047,1.00000,1.08883]),
               "D75" :np.array([0.94972,1.00000,1.22638]),
               "E"   :np.array([1.00000,1.00000,1.00000]),
               "F2"  :np.array([0.99186,1.00000,0.67393]),
               "F7"  :np.array([0.95041,1.00000,1.08747]),
               "F11" :np.array([1.00962,1.00000,0.64350])}

illuminantsxy = {"A"   : np.array([illuminants["A"  ][0]/np.sum(illuminants["A"  ]),
                                   illuminants["A"  ][1]/np.sum(illuminants["A"  ])]),    
                 "B"   : np.array([illuminants["B"  ][0]/np.sum(illuminants["B"  ]),
                                   illuminants["B"  ][1]/np.sum(illuminants["B"  ])]),    
                 "C"   : np.array([illuminants["C"  ][0]/np.sum(illuminants["C"  ]),
                                   illuminants["C"  ][1]/np.sum(illuminants["C"  ])]),    
                 "D50" : np.array([illuminants["D50"][0]/np.sum(illuminants["D50"]),
                                   illuminants["D50"][1]/np.sum(illuminants["D50"])]),    
                 "D55" : np.array([illuminants["D55"][0]/np.sum(illuminants["D55"]),
                                   illuminants["D55"][1]/np.sum(illuminants["D55"])]),    
                 "D65" : np.array([illuminants["D65"][0]/np.sum(illuminants["D65"]),
                                   illuminants["D65"][1]/np.sum(illuminants["D65"])]),    
                 "D75" : np.array([illuminants["D75"][0]/np.sum(illuminants["D75"]),
                                   illuminants["D75"][1]/np.sum(illuminants["D75"])]),    
                 "E"   : np.array([illuminants["E"  ][0]/np.sum(illuminants["E"  ]),
                                   illuminants["E"  ][1]/np.sum(illuminants["E"  ])]),    
                 "F2"  : np.array([illuminants["F2" ][0]/np.sum(illuminants["F2" ]),
                                   illuminants["F2" ][1]/np.sum(illuminants["F2" ])]),    
                 "F7"  : np.array([illuminants["F7" ][0]/np.sum(illuminants["F7" ]),
                                   illuminants["F7" ][1]/np.sum(illuminants["F7" ])]),    
                 "F11" : np.array([illuminants["F11"][0]/np.sum(illuminants["F11"]),
                                   illuminants["F11"][1]/np.sum(illuminants["F11"])])}    
    
    
#Compute the internal normlization factor for each colorspace    
_gamuts = list(colorgamuts.keys())               
for gamut in _gamuts:
    if "_" in gamut:
        il = gamut.split("_")[-1]
    else:
        il = "D50"
    white = np.array(xyz2rgb(illuminantsxy[il][0],illuminantsxy[il][1],1.0,gamut=gamut))
    extranorm = 1.0/white.max() #So the equal-power colour has a max RGB value of 1.0
    colorgamuts["%s_norm"%gamut] = extranorm
#Compute the purity of each colorspace--i.e. given the white illuminant, how white is it? Pure white is
#(1,1,1).
for gamut in _gamuts:
    if "_" in gamut:
        il = gamut.split("_")[-1]
    else:
        il = "D50"
    white = np.array(xyz2rgb(illuminantsxy[il][0],illuminantsxy[il][1],1.0,gamut=gamut))
    purity = 1.0 - np.mean(abs(1.0-white))
    colorgamuts["%s_purity"%gamut] = purity
    
def spec2rgb(wvl,spec,normalization=None,gamma=True,gamut="sRGB"):
    '''Convert a spectrum to (R,G,B) tuple, with optional normalization
    
    Parameters
    ----------
    wvl : array-like
        Array of length N containing wavelengths in nanometers
    spec : array-like
        Array of length N containing fluxes
    normalization : float, optional
        Maximum value. If not specified, defaults to 1.0.
    gamma : bool or float, optional
        If True, use the piecewise gamma-function defined for sRGB; otherwise if a float, use rgb^(1/gamma).
        If None, gamma=1.0 is used.
    gamut : str or np.ndarray(3,3)
        Color gamut to be used. For available built-in color gamuts, see colormatch.colorgamuts.
        
    Returns
    -------
    (float,float,float)
        (R,G,B).
    '''
    x,y,I = makexyz(wvl,spec)
    if normalization:
        norm = normalization
    else:
        norm = 1.0
    
    r,g,b = xyz2rgb(x,y,norm,gamut=gamut)
    colors = np.array((r,g,b))
    
    if gamma is True:
        colors[colors<0.0031308] = 12.92*colors[colors<0.0031308]
        colors[colors>=0.0031308] = 1.055*colors[colors>=0.0031308]**(1./2.4)-0.055
    elif gamma is not None:
        colors = colors**(1./gamma)
    
    r,g,b = colors
    
    return r,g,b


def specs2rgb(wvl,specs,gamma=True,gamut='sRGB'):
    '''Convert a set of spectra into RGB colour values, so that relative intensity is preserved.
    
    Parameters
    ----------
    wvl : array-like
        Wavelengths in nanometers. Must have shape (N,)
    specs : array-like
        Spectra (fluxes), of shape (M,N) where M is the number of spectra.
    gamma : bool or float, optional
        If True, use the piecewise gamma-function defined for sRGB; otherwise if a float, use rgb^(1/gamma).
        If None, gamma=1.0 is used.
    gamut : str or np.ndarray(3,3)
        Color gamut to be used. For available built-in color gamuts, see colormatch.colorgamuts.
        
    Returns
    -------
    array-like
        (M,3)-shape numpy array of R/G/B values
    '''
    interpol = interpolant(wvl)
    intensities = np.zeros((len(specs),3))
    for n in range(0,len(specs)):
        intensities[n,:] = makexyz(wvl,specs[n,:],interpolant=interpol)
    norms = intensities[:,2]/np.amax(intensities[:,2])
    colors = np.zeros((len(specs),3))
    for n in range(0,len(specs)):
        colors[n,:] = xyz2rgb(intensities[n,0],intensities[n,1],norms[n],gamut=gamut)
    if gamma is True:
        colors[colors<0.0031308] = 12.92*colors[colors<0.0031308]
        colors[colors>=0.0031308] = 1.055*colors[colors>=0.0031308]**(1./2.4)-0.055
    elif gamma is not None:
        colors = colors**(1./gamma)
    return colors
    
    