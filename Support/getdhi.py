from numarray import *

def getdhi(momentmap, header, pa, center, dhi):
   """
    NAME:
          GETDHI
   
    PURPOSE:
          Routine to obtain the diameter in of the moment 0 map at the 1e20 columndensity
   
    CATEGORY:
          Support
   
    CALLING SEQUENCE:
          GETDHI,momentmap,header,PA,center,DHI
   
   
    INPUTS:
          momentmap = 2D array with the values of the momentmap pixels
          header = fits header of the moment map
          PA = the position angle of the galaxy in degree
          center = the center of the galaxies in degrees or [hh:mm:ss,km/s]
   
    OPTIONAL INPUTS:
          -
   
    KEYWORD PARAMETERS:
          -
   
    OUTPUTS:
          DHI =  the diameter of the HI disk at 1e20 columndensity in arcsec
   
    OPTIONAL OUTPUTS:
          -
   
    PROCEDURES CALLED:
          COLUMNDENSITY, INT_PROFILEV2, INTERPOLATE, STRTRIM() STRCOMPRESS(),
          STR_SEP(), SXPAR(), STRUPCASE()
   
    EXAMPLE:
   
   
    MODIFICATION HISTORY:
          Written by P.Kamphuis 01-01-2015
   
    NOTE:
   
   """

   n_params = 5
   def _ret():  return (momentmap, header, pa, center, dhi)
   
   # COMPILE_OPT IDL2
   tmp = str_sep(strtrim(strcompress(string(center[0])), 2), ':')
   if array(tmp, copy=0).nelements() > 1:   
      x = string(center[0])
      y = string(center[1])
      x, y, True = convertradec(x, y, invert=True)
      center[0:2] = concatenate([x, y])
      center[2] = array(center[2], copy=0).astype(Float64)
   center[0] = sxpar(header, 'CRPIX1') + ((center[0] - sxpar(header, 'CRVAL1')) / sxpar(header, 'CDELT1'))
   center[1] = sxpar(header, 'CRPIX2') + ((center[1] - sxpar(header, 'CRVAL2')) / sxpar(header, 'CDELT2'))
   yrange = concatenate([-0.5 * sxpar(header, 'BMAJ'), 0.5 * sxpar(header, 'BMAJ')])
   test = 1
   int_profilev2(momentmap, majprofile, header=header, xcenter=center[0], ycenter=center[1], pa=pa, range=yrange, axis=xaxis, rotimage=test)
   clevels = 1.24773e20
   clevels, center[2], concatenate([sxpar(header, 'BMAJ') * 3600., sxpar(header, 'BMIN') * 3600.]), True = columndensity(clevels, center[2], concatenate([sxpar(header, 'BMAJ') * 3600., sxpar(header, 'BMIN') * 3600.]), ncolumn=True)
   clevels = clevels / 1000.
   
   #First we need to cut the major axis
   #profile so that it does not include
   #noise
   #From the max if it drops below 0 we
   #assume that is the end of the galaxy
   #This could go wrong with two galaxies
   #in the cube
   
   maxprof = max(majprofile)
   maxisat = where(ravel(maxprof == majprofile))[0]
   #determine lower bound
   stop = 0
   lowerbound = maxisat[0]
   while bitwise_and(stop == 0, lowerbound > 0):
      lowerbound = lowerbound - 1
      if majprofile[lowerbound] <= 0.:   
         stop = 1.
   stop = 0
   upperbound = maxisat[0]
   while bitwise_and(stop == 0, upperbound < array(majprofile, copy=0).nelements() - 2):
      upperbound = upperbound + 1
      if majprofile[upperbound] <= 0.:   
         stop = 1.
   majprofile[0:(lowerbound)+1] = 0.
   majprofile[upperbound:(array(majprofile, copy=0).nelements() - 1)+1] = 0.
   tmp = where(ravel(majprofile >= clevels))[0]
   lowrad = 1.
   if tmp[0] > -1:   
      interpolate(concatenate([xaxis[tmp[0]], xaxis[tmp[0] - 1]]), concatenate([majprofile[tmp[0]], majprofile[tmp[0] - 1]]), newradii=clevels, output=lowrad)
      highrad = 1.
      interpolate(concatenate([xaxis[tmp[array(tmp, copy=0).nelements() - 1]], xaxis[tmp[array(tmp, copy=0).nelements() - 1] + 1]]), concatenate([majprofile[tmp[array(tmp, copy=0).nelements() - 1]], majprofile[tmp[array(tmp, copy=0).nelements() - 1] + 1]]), newradii=clevels, output=highrad)
      tmp = str_sep(strtrim(strcompress(sxpar(header, 'CTYPE1')), 2), '---')
      ctype1 = tmp[0]
      if strupcase(ctype1) == 'RA':   
         dhi = (abs(highrad - lowrad)) * cos(sxpar(header, 'CRVAL2') * _sys_pi / 180.) * 3600.
      else:   
         dhi = abs(highrad - lowrad)
   else:   
      dhi = 0.
   
   return _ret()




