import pandas as pd
import numpy as np
import os
import pickle as pk

import pandexo.engine.justdoit as jdi # THIS IS THE HOLY GRAIL OF PANDEXO
import pandexo.engine.justplotit as jpi

import analysis as an  # just for constants
import tess_yield as ty

_modulepath =  os.path.split(__file__)[0]
if _modulepath=='': _modulepath = '.'

#_nexsciPlanetTable = os.path.expanduser('~/proj/jwst/go1/playground/nexsci_planets_20180123_ed.csv')
_nexsciPlanetTable = _modulepath + '/nexsci_planets_20180123_ed.csv'



"""A suite of useful tools for planning JWST observations.

 :FUNCTIONS:

  :func:`run_pandexo` -- Run a simulation of PandExo and avoid messing
                         with scripts.

  :func:`loadplanets` -- Load the NASA Exoplanet Archive's list of
                         exoplanets & their properties.


 :DEPENDENCIES:

  Python : pandas, numpy, pickle, os

  Other : PandExo:
             https://natashabatalha.github.io/PandExo/
             https://github.com/natashabatalha/PandExo/

  IJMC's : analysis (just for constants)
           tess_yield (for :func:`tess_yield.estimatePlanetParams` function)

2018-01-26 14:19 IJMC: Created.

"""

def run_pandexo(planetname='WASP-12 b', mode='transit', instrument = 'MIRI LRS', subarray=None, _modelspectrum=None, modelwave='um', ntransits=1, noise_floor=0., refband='k', plotres=30, mrmode='weiss2016', _outputpath='.', retpandexo=True, reddict=False):
    """Run a simulation of PandExo -- and avoid messing with scripts.
    
    INPUTS:
 
     planetname : str
       Name of a valid planet: 'WASP-12 b', 'HD 209458 b', 'GJ 9827 d',
       etc.  Loaded from :func:`loadplanets`
 
     mode : str
       'transit' or 'eclipse'.  
 
     instrument : str
       Choose from: 'NIRCam F444W', 'NIRSpec Prism', 'NIRSpec G395M',
       'NIRCam F322W2', 'NIRSpec G395H', 'NIRSpec G235H', 
       'NIRSpec G235M', 'NIRSpec G140M', 'NIRSpec G140H', 'MIRI LRS',
       'NIRISS SOSS' (or 'WFC3 G141' for HST).

     subarray : str or None
        Detector subarray to use with your choice of instrument. Leave
        as None unless you know what you're doing.
     
     _modelspectrum : str or None
        Name of model spectrum to use. This should be a two-column,
        space-delimited ASCII file with wavelength in the first column
        and depth (transit or eclipse) in the second column.
 
       OR:
 
        If None, then assumes a simple, flat spectrum in transit (using
        size of planet) or eclipse (using equilibrium temperature of
        planet)

     _outputpath : str
        location where output pickle will be saved.
 
     modelwave : str
        Wavelength units (readable by AstroPy) for input
        _modelspectrum: "um","nm" ,"Angs", "secs" (for phase curves)
 
     ntransits : int
        Number of transits (or eclipses) to observe.
 
     noise_floor : scalar
        Assumed noise floor, as accepted for input by PandExo JWST
        simulator.  This can be a fixed level or it can be a filepath
 
     refband : str
        'j' or 'h' or 'k' for stellar flux.  Magnitude is read from
        table of planet objects.
 
     plotres : int
        Final spectral resolution for returned simulated spectrum.
 
     mrmode : str
        Mass-Radius mode for estimating planet masses when they aren't
        known. See :func:`loadplanets` for more details.

     retpandexo : bool
        Whether to also return pandexo 'results' dictionary.

     retdict : bool
        Whether to also return the exo_dict and inst_dict inputs used
        to run PandExo (mainly for debugging/troubleshooting)

    :EXAMPLE:
      ::

       import jwst
       import pylab as py

       out = jwst.run_pandexo(planetname='WASP-12 b', mode='eclipse', instrument='NIRISS SOSS', retpandexo=True)

       py.figure()
       py.plot(out['result']['OriginalInput']['model_wave'], 1e6*out['result']['OriginalInput']['model_spec'], '-r')
       py.errorbar(out['wavelength'][0], 1e6*out['spectrum'][0], 1e6*out['uspectrum'][0], fmt='.k', elinewidth=2, mew=2)
       py.ylabel('Depth [ppm]', fontsize=20)
       py.xlabel('Wavelength' , fontsize=20)
       py.minorticks_on()
       py.xlim(out['wavelength'][0].min()-0.1, out['wavelength'][0].max()+0.1)


    :TO-DO:
      Phase curve mode.

      User-specified planet parameters

      Modern error handling (instead of just crashing)

    """
    # 2018-01-26 11:59 IJMC: Created
    
    knownmass = False
    planets = loadplanets(knownmass=knownmass, mrmode=mrmode)
    planetindex = planets.pl_name==planetname
    if not planetindex.any():
        print "Planet %s not found -- bombing out" % planetname
        stopp
    elif planetindex.sum()>1:
        print "More than one copy of planet %s found in table -- check your file. Bombing out." % planetname
        stopp
    else:
        planet = planets[planetindex]

    # Create filename:
    iter = 0
    outputfilename = ('PE_%s_%s_%s_%s_n%i_%04i.pickle' % (planetname, mode, instrument, subarray, ntransits, iter)).replace(' ', '_')
    if _modelspectrum is None or _modelspectrum=='constant':
        outputfilename = outputfilename.replace(mode, mode+'-constant')
    while os.path.isfile(outputfilename):
        iter += 1
        outputfilename = outputfilename.replace('%04i.pickle' % (iter-1), '%04i.pickle' % iter)

    # Set other options (under the hood):
    if refband.lower()=='k':
        ref_wave = 2.2
    elif refband.lower()=='h':
        ref_wave = 1.6
    elif refband.lower()=='j':
        ref_wave = 1.25

    mag = float(getattr(planet, 'st_' + refband))

    if hasattr(planet, 'st_metallicity'): # what is it really called?
        metal = float(planet.st_metallicity)
    else:
        metal = 0.

    if mode=='transit':
        f_unit = 'rp^2/r*^2'
    elif mode=='eclipse':
        f_unit = 'fp/f*'
    else:
        print "I don't know mode '%s' -- bombing out" % mode
        stopp


    # Pandexo: Begin!
    exo_dict = jdi.load_exo_dict()

    exo_dict['observation']['sat_level'] = 80                  #saturation level in percent of full well
    exo_dict['observation']['sat_unit'] = '%'                  # other option = 'e' for electrons
    exo_dict['observation']['noccultations'] = ntransits  
    # 2018-02-11 16:03 IJMC: Changed from 4.0*60.0*60.0/'total':
    exo_dict['observation']['baseline'] = 1.0
    exo_dict['observation']['baseline_unit'] = 'frac'         #total obersving time, other option 'frac' = in/out
    exo_dict['observation']['noise_floor'] = noise_floor  

    exo_dict['star']['type'] = 'phoenix'     
    exo_dict['star']['mag'] = mag            
    exo_dict['star']['ref_wave'] = ref_wave  
    exo_dict['star']['temp'] = int(planet.Teff)
    exo_dict['star']['metal'] = float(metal)   
    exo_dict['star']['logg'] = 2+np.log10(float(an.G*an.msun*planet.st_mass / \
                                                (an.rsun*planet.st_rad)**2))
    exo_dict['star']['radius'] = float(planet.st_rad)
    exo_dict['star']['r_unit'] = 'R_sun'




    if _modelspectrum is None or _modelspectrum=='constant': 
        exo_dict['planet']['type'] = 'constant'
        if mode=='transit':
            exo_dict['planet']['f_unit'] = 'rp^2/r*^2'      
        elif mode=='eclipse':
            exo_dict['planet']['f_unit'] = 'fp/f*'        
    else:
        exo_dict['planet']['type'] ='user'
        exo_dict['planet']['exopath'] = _modelspectrum
        exo_dict['planet']['w_unit'] = modelwave         
        exo_dict['planet']['f_unit'] = f_unit         

    exo_dict['planet']['transit_duration']      = float(planet.pl_trandur*24*3600.)
    exo_dict['planet']['td_unit'] = 's'      

    exo_dict['planet']['temp'] = float(planet.Teq)
    exo_dict['planet']['radius'] = float(planet.pl_radj)
    exo_dict['planet']['r_unit'] = 'R_jup'            
    exo_dict['planet']['i']       = float(planet.pl_orbincl)             
    #exo_dict['planet']['ars']     = float((planet.pl_orbsmax*an.AU) / (planet.st_rad*an.rsun))
    exo_dict['planet']['period']  = float(planet.pl_orbper)                  


    # Now load in the instrument:
    inst_dict = jdi.load_mode_dict(instrument)
    if subarray is not None: inst_dict['configuration']['detector']['subarray']     = subarray   
    ## 2018-02-11 15:51 IJMC: Commented out:
    #inst_dict['configuration']['detector']['nsamp']        = None  
    #inst_dict['configuration']['detector']['samp_seq']     = None  


    result = jdi.run_pandexo(exo_dict, inst_dict, output_file=_outputpath + outputfilename)
    if not np.isfinite([result['timing'][key] for key in result['timing'].keys()]).all():
        print "Something went wrong with simulation. Maybe star is too bright? Bombing out."
        stopp

    x,y, e = jpi.jwst_1d_spec(result, R=plotres, model=True, plot=False) #, x_range=[.8,1.28]) # num_tran=10

    # Prepare dictionary for returning to user:
    
    ret = dict(wavelength=x, spectrum=y, uspectrum=e, planet=planet)
    if retpandexo:
        ret['result'] = result
        ret['outputfile'] = _outputpath + outputfilename
    
    if reddict:
        ret['exo_dict'] = exo_dict
        ret['inst_dict'] = inst_dict

    return ret

    

def loadplanets(lam=4.5, _table=_nexsciPlanetTable, knownmass=False, mrmode='weiss2016'):
    """
    Load the NASA Exoplanet Archive's list of exoplanets & their properties.

    lam -- wavelength in microns for BLACKBODY eclipse depths. Can be
           a vector or scalar.

    _table -- path to IPAC planet table. Be sure to get all the right 
              fields

    knownmass (bool):
      whether to only use planets with KNOWN, measured masses.
      Otherwise, infer masses using 'mrmode'

    mrmod (str);
      'weiss2016' or 'wolfgang2016' (for tess.estimatePlanetParams())

    """
    # 2018-01-25 19:51 IJMC: Created
    
    planets = pd.read_csv(_table, comment='#')
    planets = planets.iloc[(planets.pl_tranflag<>0).nonzero()]
    planets['Rplanet'] = planets.pl_radj * an.rjup/an.rearth

    if knownmass:
        planets = planets.iloc[np.isfinite(planets.pl_bmassj).nonzero()]

    # Estimate semimajor axes (often missing from Nexsci):
    noa = np.logical_not(np.isfinite(planets.pl_orbsmax))
    planets.pl_orbsmax[noa] = (((planets.pl_orbper/365.24)**2 * planets.st_mass)**0.3333)[noa]

    pparams = ty.estimatePlanetParams(planets.st_teff, planets.st_rad, planets.pl_orbper, planets.pl_orbsmax, planets.Rplanet, mrmode=mrmode)

    planets['Mplanet'] = planets.pl_bmassj * an.mjup/an.mearth

    if not knownmass:
        unknownmass = np.logical_not(np.isfinite(planets.Mplanet))
        planets.Mplanet[unknownmass] = pparams.Mplanet[unknownmass]

    for key in ['a', 'Teq', 'Teff', 'Rstar']:  planets[key] = pparams[key]
    planets['transitdepth'] = ((planets.Rplanet * an.rearth)/ (planets.Rstar*an.rsun))**2

    return planets

