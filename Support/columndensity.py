import numpy as np

def columndensity(levels, beam, vsys=1., vwidth=1., ncolumn=False, densities=None, arcsquare=False, msolar=False):
   """
    NAME:
          COLUMNDENSITY
   
    PURPOSE:
          Program to calculate the column density based on the given flux level
   
    CATEGORY:
          Support
   
    CALLING SEQUENCE:
          COLUMNDENSITY,levels,vsys,beam,VWIDTH=vwidth,/NCOLUMN,DENSITIES=densities,/ARCSQUARE
   
   
    INPUTS:
          levels = Is an array or scalar with the flux levels (mJy/beam) that are
          to be transformed to column densities
          vsys = is scalar with the systemic velocity of the system in km/s
          beam = the beam of the observations in arcsec if it is 1 parameter then a
          circular beam is assumed
   
    OPTIONAL INPUTS:
          -
   
    KEYWORD PARAMETERS:
          vwidth = the channel width of the observations must be given when
          flux is in mJy/beam if not given it is assumed that the flux level
          is mJy/beam*km/s
          /NCOLUMN - if set it is assumed that level is column densities
          /ARCSQUARE - if set it is assumed that levels is in
                       mJy/arcsec^2 this is specifically useful
                       for converting tirific's SBR to
                       columndensities program to calculate the column
                       density
           /MSOLAR - Give output in M_solar pc^2,  if ncolumn is set it
                     is assumed the input is in M_solar pc^2
   
    OUTPUTS:
          levels = the transformed units if densities is unset
   
    OPTIONAL OUTPUTS:
          DENSITIES = if given the result will be published there
          otherwise levels will be converted. IF /NCOLUM is given and
          densities is given it will be assumed that densities is the
          input and the levels output that have to be transforemed to flux levels in mJy/beam
   
    PROCEDURES CALLED:
          -
   
    EXAMPLE:
   
   
    MODIFICATION HISTORY:
          16-12-2016 P.Kamphuis; Added the option to have the levels
          converted to Msolar pc^2
          Written by P.Kamphuis 01-01-2015
   
    NOTE:
         A more rough shod way would be n(HI)=1e+21 F(mJy/beam km/s)/bmaj(arcsec)/bmin(arcsec)
   """
   
   # COMPILE_OPT IDL2
   levels= np.array(levels)
   f0 = 1.420405751786e9           #Hz rest freq
   c = 299792.458                  #light speed in km/s 
   pc = 3.086e+18                  #parsec in cm
   solarmass = 1.98855e30          #Solar mass in kg
   mhi = 1.6737236e-27 #neutral hydrogen mass in kg
   if vsys > 10000:   
      vsys = vsys / 1000.
   f = f0 * (1 - (vsys / c))             #Systemic frequency
   if arcsquare:   
      #we do not want to use the beam but
      #need to correct for that with a
      #factor the difference is a factor
      #1.1330900354567984..., which comes
      #from the normalisation to the unit
      #"beam". It is the Integral of a 2D
      #Gaussian with unity HPBW. From
      #Josh's e-mail.
      hiconv = 605.7383 * 1.823e18 * (2. * np.pi / (np.log(256.)))
      if ncolumn:   
         #I can't find the redshift dependence
         #so let's leave it out in this
         #one note that this one gives you
         #column densities in mJy/arcsec^2
         
         if densities != None:   
            levels[:] = densities[:]
         #if msolar is set we want to convert
         #back to column densities
         if msolar:   
            levels[:] = levels[:] * solarmass / (mhi * pc ** 2)         
         levels[:] = levels[:] / (hiconv * vwidth)         
      else:   
         
         nh = hiconv * levels * vwidth
                  
         if densities != None:   
            densities[:] = nh[:]
         else:   
            levels[:] = nh[:]
      
   else:   
      if len(beam) == 1:   
         beam = [beam,beam]
      b = beam[0] * beam[1]
      if ncolumn:   
         if densities != None:   
            levels[:] = densities[:]
         #if msolar is set we want to convert
         #back to column densities
         if msolar:   
            levels[:] = levels[:] * solarmass / (mhi * pc ** 2)
         
       
         
         tk = levels[:] / (1.823e18 * vwidth)
         levels[:] = tk[:] / (((605.7383) / (b)) * (f0 / f) ** 2)
      else:   
         tk = ((605.7383) / (b)) * (f0 / f) ** 2 * levels[:]
         nh = 1.823e18 * tk * vwidth
         
         if densities != None:
            densities[:] = nh[:]
         else:   
            levels[:] = nh[:]
   
   if ncolumn and msolar:   
      levels[:] = levels[:] * mhi * pc ** 2 / solarmass
   
   return levels

