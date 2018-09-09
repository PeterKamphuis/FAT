
from math import floor
def convertradec(ra, dec, Invert=False, Colon=False):
   """
    NAME:
          CONVERTRADEC
   
    PURPOSE:
          Program to convert DEG to RA and DEC or vice versa
   
    CATEGORY:
          Support
   
    CALLING SEQUENCE:
          CONVERTRADEC,RA,DEC,/INVERT,/COLON
   
   
    INPUTS:
          RA = Right Ascension in Degrees
          DEC = Declination in degrees
   
    OPTIONAL INPUTS:
          -
   
    KEYWORD PARAMETERS:
          /INVERT - if set it is assumed that RA and DEC are in hour angle
          /COLON - if set it the output will be hh:mm:ss instead of
                   correctly formatted
   
    OUTPUTS:
          RA = the transformed units
          DEC = the transformed units
   
    OPTIONAL OUTPUTS:
   
    PROCEDURES CALLED:
          STRTRIM(),STRCOMPRESS(),STRING(),STR_SEP(),FLOOR()
   
    EXAMPLE:
   
   
    MODIFICATION HISTORY:
          15-12-2016 P. Kamphuis; Improved handling of -00d
          Written by P.Kamphuis 01-01-2015
   
    NOTE:
   
   """

   n_params = 2
 
   
   # COMPILE_OPT IDL2
   xpos = []
   ypos =[]
   #make sure ra and dec are lists
   try:
      len(ra)
   except TypeError:
      ra = [ra]
      dec= [dec]
   if not Invert:
      for i in range(len(ra)):
         xposh = floor((ra[i] / 360.) * 24.)
         xposm = floor((((ra[i] / 360.) * 24.) - xposh) * 60.)
         xposs = (((((ra[i] / 360.) * 24.) - xposh) * 60.) - xposm) * 60
         yposh = floor(abs(dec[i] * 1.))
         yposm = floor((((abs(dec[i] * 1.)) - yposh) * 60.))
         yposs = (((((abs(dec[i] * 1.)) - yposh) * 60.) - yposm) * 60)
         sign = (round(dec[i] / abs(dec[i])))
         if Colon:   
            xpos.append("{}:{:02d}:{:.1f}".format(xposh,xposm,xposs))
            ypos.append("{}:{:02d}:{:.1f}".format(yposh,yposm,yposs))
            if sign < 0.:   
               ypos[i] = '-{}'.format(ypos[i])
         else:   
            xpos.append("{}h{:02d}m{:.1f}s".format(xposh,xposm,xposs))
            ypos.append("{}d{:02d}'{:.1f}\"".format(yposh,yposm,yposs))
            if sign < 0.:   
              ypos[i] = '-{}'.format(ypos[i])
   else:
      import re
      for i in range(len(ra)):
         tmp = [float(item) for item in list(filter(None,re.split("[:hdms\'\"]+",ra[i])))]
         xpos.append((tmp[0]+((tmp[1]+(tmp[2]/60.))/60.))*15.)
         tmp = [float(item) for item in list(filter(None,re.split("[:hdms\'\"]+",dec[i])))]
         ypos.append(abs(tmp[0])+((tmp[1]+(tmp[2]/60.))/60.))
         if dec[i][0] == '-':   
            ypos[i] = -ypos[i]
   
   return xpos,ypos

