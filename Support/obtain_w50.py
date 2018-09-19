from numarray import *

def obtain_w50(cube, mask, hed, w50):
   """
    NAME:
          OBTAIN_W50
   
    PURPOSE:
          Routine to obtain W50 from a 3D fits cube array
   
    CATEGORY:
          Support
   
    CALLING SEQUENCE:
          OBTAIN_W50,Cube,mask,hed,W50
   
   
    INPUTS:
          Cube =  Array of the cube to be analysed
          Mask = Mask to identify the emission in the cube
          hed = the header of the cube
   
    OPTIONAL INPUTS:
          -
   
    KEYWORD PARAMETERS:
          -
   
    OUTPUTS:
          W50 = the velocity width at 50% from the peak value in the
          spectrum in km/s
   
    OPTIONAL OUTPUTS:
          -
   
    PROCEDURES CALLED:
          BUILDAXII, STRUPCASE(), STRTRIM(), SXPAR()
   
    EXAMPLE:
   
   
    MODIFICATION HISTORY:
          Written 01-01-2015 P.Kamphuis v1.0
   
    NOTE:
   
   """

   n_params = 4
   def _ret():  return (cube, mask, hed, w50)
   
   # COMPILE_OPT IDL2
   buildaxii(hed, xaxis, yaxis, zaxis=zaxis)
   if strupcase(strtrim(sxpar(hed, 'CUNIT3'), 2)) == 'M/S':   
      divide = 1000.
   else:   
      divide = 1.
   zaxis = zaxis / divide
   profile = dblarr(array(zaxis, copy=0).nelements())
   maskedcube = dblarr(array(cube[0,0,:], copy=0).nelements(), array(cube[0,:,0], copy=0).nelements(), array(cube[:,0,0], copy=0).nelements())
   maskedcube[where(ravel(mask == 1))[0]] = cube[where(ravel(mask == 1))[0]]
   for i in arange(0, (array(zaxis, copy=0).nelements() - 1)+(1)):
      profile[i] = total(maskedcube[i,:,:])
   
   
   maxprof = max(profile)
   
   above50 = where(ravel(profile > maxprof / 2.))[0]
   if array(above50, copy=0).nelements() > 2:   
      w50 = zaxis[above50[array(above50, copy=0).nelements() - 1]] - zaxis[above50[0]]
   else:   
      w50 = sxpar(hed, 'CDELT3') / divide
   
   
   
   return _ret()

