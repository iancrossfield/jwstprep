from pylab import *
import pandas as pd
import os
import pdb

import analysis as an

_tessdata = os.path.expanduser('~/docs/papers/contributing/k2_catalog_c0-4/paper/')
_k2data   = os.path.expanduser('~/proj/mdwarfs/kepler2/catalogs/')

def k2_table(_k2fn=_k2data + 'output28.csv'):
    import k2fits
    """Load my K2 planet candidate tables.
    """
    k2 = pd.read_csv(_k2fn, skipinitialspace=True)
    try:
        candidates = []
        for starname in k2.starname:
            candidates.append(k2fits.getCandidates(starname, k2fits._starcat, k2fits._planetcat)[0][0])
        k2['Teff'] = np.array([can.teff for can in candidates])
        k2['Bmag'] = np.array([can.Bmag for can in candidates])
        k2['Vmag'] = np.array([can.Vmag for can in candidates])
        k2['Jmag'] = np.array([can.Jmag for can in candidates])
        k2['Hmag'] = np.array([can.Hmag for can in candidates])
        k2['Kmag'] = np.array([can.Kmag for can in candidates])
    except:
        print "Tried to load K2 candidate parameters, but couldn't..."
        
    colnames = [col.lower() for col in k2.columns]
    if not 'rplanet' in colnames: k2['Rplanet'] = k2['$R_P$']
    if not 'rstar' in colnames:   k2['Rstar']   = k2['$R_*$']
    if not 'mstar' in colnames:   k2['Mstar']   = k2['$M_*$']
    if not 'teq' in colnames:     k2['Teq'] = 250 * k2['$S_{inc}$']**0.25
    k2['T14'] = k2['$T_{14}$'] / 24.
    if not 'a' in colnames:
        k2['a'] = ((k2['$P$']/365.24)**2 * k2.Mstar)**0.333333        
    if not 'teff' in colnames:
        k2['Teff'] = 5777. * k2['$S_{inc}$'] * k2.a**2 / k2.Rstar**2
    if not 'mplanet' in colnames:
        Mplanet = np.zeros(len(k2), float)
        #large = (k2.Rplanet > 4).nonzero()[0]
        #small = (k2.Rplanet < 4).nonzero()[0]
        #tiny  = (k2.Rplanet < 1.5).nonzero()[0]
        #Mplanet[small]  = 3.53 * k2.Rplanet[small]**0.6
        #Mplanet[tiny] = 1.333*np.pi*(k2.Rplanet**3 * (1.65 + 3.69*k2.Rplanet))[tiny]
        #Mplanet[large]  = k2.Rplanet[large]**2.06
        Mplanet = k2.Rplanet**2.06
        Mplanet[(Mplanet > 6000).nonzero()[0]] = 6000.
        k2['Mplanet'] = Mplanet
        #pdb.set_trace()
    if not 'transitdepth' in colnames:
        k2['transitdepth'] = (k2['$R_P/R_*$']/100)**2

    if not 'k' in colnames:
        jupmass = k2.Mplanet/317.8
        k2['K'] = 203 * k2['$P$']**-0.3333 * (jupmass) * (k2.Mstar + 9.55e-4 * jupmass)**-0.66667
        
    return k2

def tess_table(_tessfn=_tessdata + 'Sullivan+2015.csv'):
    """Load the Sullivan et al. TESS simulation and add extra columns.
    """
    tess = pd.read_csv(_tessfn, skipinitialspace=True, delimiter=',')
    colnames = [col.lower() for col in tess.columns]
    if not 'mstar' in colnames:
        tess['Mstar'] = tess.Rstar
        tess.Mstar[(tess.Mstar > 1.5) * (tess.Rstar > 1.5)] = 1.5
    if not 'a' in colnames:
        tess['a'] = ((tess.P/365.24)**2 * tess.Mstar)**0.333333
    if not 'mplanet' in colnames:
        tess['Mplanet'] = (an.mjup/an.mearth) * (tess.K / 28.43) * (tess.P/365.24)**0.3333 * tess.Mstar**0.66667
        #tess['Mplanet'] = tess.Rplanet**2.06
    if not 'teq' in colnames:
        tess['Teq'] = 250 * tess['S/SEarth']**0.25
    tess['transitdepth'] = ((tess.Rplanet * an.rearth) / (tess.Rstar * an.rsun))**2
    tess['T14'] = (tess.P / np.pi) * np.arcsin((tess.Rstar * an.rsun) / (tess.a * an.AU) * np.sqrt((1 - tess.transitdepth**0.5)**2 - 0.5**2))

    return tess

def atmoparams(planetTable=None, lam=None, MMW=None, mh=None, **k2):
    """Compute atmospheric parameters for a table of planets.

    Necessary columns in planetTable:
       Rplanet - in R_Earth
       Mplanet - in M_Earth
       Teq  
       Rstar   - in R_Sun
       Teff

    If it lacks a column "MMW" [Mean Molecular Weight, in amu] we
    naively try to interpolate based on planet mass and/or radius.
    """
    
    if planetTable is None:
        planetTable = table(**k2)
        
    r_transition = 1.5, 2.2
    mmw_transition = 18, 3
    if MMW is None:
        MMW = np.zeros(len(planetTable))
        noH2 = planetTable.Rplanet < r_transition[0]
        mostlyH2 = planetTable.Rplanet > r_transition[1]
        MMW = mmw_transition[0] + np.diff(mmw_transition) / np.diff(r_transition) * (planetTable.Rplanet - r_transition[0])
        MMW[noH2] = mmw_transition[0]
        MMW[mostlyH2] = mmw_transition[1]

    if mh is None:
        mh = 3 * (planetTable.Mplanet / 318.)**-1.4  # SOlar System mass-metallicity relation
        noH2 = planetTable.Rplanet < r_transition[0]
        mh[noH2] = 1e4
        
    planetTable['mh'] = mh
    planetTable['MMW'] = MMW
    planetTable['g_mks'] = an.G * an.mearth/an.rearth**2 * (planetTable.Mplanet / planetTable.Rplanet**2)
    planetTable['H_mks'] = (an.k * planetTable.Teq) / (planetTable.MMW / 6e26 * planetTable.g_mks)
    planetTable['scaleheightarea'] = planetTable.H_mks * (planetTable.Rplanet * an.rearth) / (planetTable.Rstar * an.rsun)**2
    planetTable['SHA_multiplier'] = 5.
    planetTable['transmission_amplitude'] = planetTable.SHA_multiplier * planetTable.scaleheightarea
    
    if lam is not None:
        eclipsedepth = []
        if hasattr(planetTable, 'index'):
            ind = planetTable.index
        else:
            ind = range(len(planetTable))
        for ii in ind:
            eclipsedepth.append(blam(planetTable.Teq[ii], lam) / blam(planetTable.Teff[ii], lam) * planetTable.transitdepth[ii])
        planetTable['eclipsedepth'] = eclipsedepth
        
    
    return planetTable


def MOST(table):
    """A crude model of MOST observing for a table of planets.
    """
    # Compute a 2D map of the CVZ:
    ra0, dec0 = np.linspace(0, 24, 721), np.linspace(-90, 90, 361.)
    raa, decc = np.meshgrid(ra0, dec0)
    equatorial = (decc <= 35) * (decc >= -17.)
    ecliptic = np.sin(ra0*np.pi*15./180.) * 23.
    sunsensor = (decc <= (ecliptic+32.)) * (decc >= (ecliptic-32.))
    cvz = np.logical_and(sunsensor, equatorial)
    cvzs = np.array([dec0[cvzcol].min() for cvzcol in cvz.T])
    cvzn = np.array([dec0[cvzcol].max() for cvzcol in cvz.T])
    #viz = cvz.astype(float)
    #viz[decc>cvzn] = 1.0-(decc-cvzn)[decc>cvzn]*0.7/45.
    #viz[decc<cvzs] = 1.0-(cvzs-decc)[decc<cvzs]*0.7/45.
    #viz[viz<0.3] = 0.

    # Intepolate to get the CVZ's North and South edge at each Dec
    cvzs2 = np.interp(table.RA/15., ra0, cvzs)
    cvzn2 = np.interp(table.RA/15., ra0, cvzn)

    # Now compute the equatorial, sunsensor, and CVZ regions:
    equatorial = (table.Dec <= 35) * (table.Dec >= -17.)
    ecliptic = np.sin(table.RA*np.pi/180.) * 23.
    viz_sunsensor = (table.Dec <= (ecliptic+32.)) * (table.Dec >= (ecliptic-32.))
    vizz = np.logical_and(viz_sunsensor, equatorial)
    sep_from_cvz = 1. - vizz.astype(float)
    sep_from_cvz[table.Dec>cvzn2] = (table.Dec-cvzn2)[table.Dec>cvzn2]
    sep_from_cvz[table.Dec<cvzs2] = (cvzs2-table.Dec)[table.Dec<cvzs2]

    vizz[table.Dec>cvzn2] = 1.0-(table.Dec-cvzn2)[table.Dec>cvzn2]*0.7/45.
    vizz[table.Dec<cvzs2] = 1.0-(cvzs2-table.Dec)[table.Dec<cvzs2]*0.7/45.
    vizz[vizz<0.3] = 0.
    
    table['MOST_vis'] = vizz
    table['MOST_sunsensor'] = viz_sunsensor
    table['MOST_sep_from_cvz'] = sep_from_cvz
    
    return table

def estimatePlanetParams(teff, rstar, per, a, rp, redist=0.3, ab=0.2, mrmode='wolfgang2016'):
    """
    Given a set of planet parameters, estimate/calculate additional parameters.

        planetparams = estimatePlanetParams(validtic.iloc[ii].Teff, validtic.iloc[ii].rad, \
                                        thisres.p, thisres.a, thisres.r)
        ty.atmoparams(planetparams);

    wolfgang2016
    weiss2016

    TO-DO: use composition models to set minimum allowable mass
    """
    teq = teff/np.sqrt(215.1*a/rstar) * (redist * (1. - ab))**0.25
    if mrmode.lower()=='wolfgang2016':
        c, gamma, sm1, beta = 2.7, 1.3, 1.9, 0.
        c, gamma, sm1, beta = 1.6, 1.8, 2.9, 0.
        mp0 = c*rp**gamma
        u_mp0 = np.sqrt(sm1**2 + beta*(rp-1))
        mp = np.random.normal(mp0, u_mp0)
    if mrmode.lower()=='weiss2016':
        ind = rp > 1.5
        mp1 = 3.53* rp**0.60
        mp2 = (1.65 + 3.69 * rp) * (1.65 + 3.69 * 1.5) * (1.333*3.14*(1.5*6.378e6)**3)/5.97e21  # radius & mass of earth
        mp = 0. + mp2
        mp[ind] = mp1[ind]
        
    mp[mp<0] = 0.001  # minimum allowable mass
    res = pd.DataFrame(dict(per=per, Rplanet=rp))
    res['a']   = a
    res['Teq'] = teq
    res['Mplanet']  = mp
    res['gplanet']  = 9.8*mp/res.Rplanet**2
    res['Teff'] = teff
    res['Rstar'] = rstar
    return res
