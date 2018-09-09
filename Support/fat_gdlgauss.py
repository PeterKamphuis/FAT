from numarray import *

def fat_gdlgauss(axis, profile):
   """
    NAME:
          FAT_GDLGAUSS
   
    PURPOSE:
          Program to simulate GAussfit in the GDL environment
   
    CATEGORY:
          Support
   
    CALLING SEQUENCE:
          FAT_GDLGAUSS, axis, profile, coeff
   
   
    INPUTS:
          axis = xaxis of the profile to be fitted
       profile = yaxis of the profile to be fitted
    OPTIONAL INPUTS:
          -
   
    KEYWORD PARAMETERS:
          -
   
    OUTPUTS:
          -
   
    OPTIONAL OUTPUTS:
          -
   
    PROCEDURES CALLED:
          EXP(), MAX(), MEAN(), MPFITFUN(), REPLICATE()
   
    EXAMPLE:
   
   
    MODIFICATION HISTORY:
          Written 02-06-2016 P.Kamphuis v1.0
   
    NOTE:
   
   
   """

   n_params = 2
   
   error = (0.5)*ones([array(axis, copy=0).nelements()])
   start = concatenate([max(profile), mean(axis), max(axis) / 2.])
   coeff = mpfitfun('FATGAUSS', axis, profile, error, start, quiet=True)
   f = coeff[0] * exp(-0.5 * ((axis[:] - coeff[1]) / coeff[2]) ** 2)
   
   return f


def fatgauss(x, p):
   n_params = 2
   
   return p[0] * exp(-0.5 * ((x - p[1]) / p[2]) ** 2)

