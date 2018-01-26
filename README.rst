# jwstprep

"""A suite of useful tools for planning JWST observations.

 :FUNCTIONS:

  :func:`run_pandexo` -- Run a simulation of PandExo and avoid messing
                         with scripts.

  :func:`loadplanets` -- Load the NASA Exoplanet Archive's list of
                         exoplanets & their properties.

 :DEPENDENCIES:

  - Python : pandas, numpy, pickle, os

  - Other : PandExo:
            - https://natashabatalha.github.io/PandExo/
            - https://github.com/natashabatalha/PandExo/

  - IJMC's :
           - analysis (just for constants)
           - tess_yield (for :func:`tess_yield.estimatePlanetParams` function)

:INSTALLATION:
  - Download all the files in this repository 
  - Put them somewhere in your Python path
  - Follow the example in :func:`run_pandexo`

:TO-DO:
  - Anything with phase curves
  - Make this a pythonic 'module'/'package' instead of some files
  
  

2018-01-26 14:19 IJMC: Created.

"""
