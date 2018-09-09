from numarray import *

def fat_smooth(map, sigma, gaussian=None, box=None, mask=None):
   """
    NAME:
          FAT_SMOOTH
   
    PURPOSE:
          Routine to apply a gaussian or box averaging smooth kernel to an image
   
    CATEGORY:
          Support
   
    CALLING SEQUENCE:
          FAT_SMOOTH,map,sigma
   
   
    INPUTS:
          map = moment 1 map of the galaxy
        sigma = the standard deviation of the gaussian to be used, the
                gaussian kernel will have the size of 10 times the FWHM defined
                by this sigma.
    OPTIONAL INPUTS:
          -
   
    KEYWORD PARAMETERS:
          /GAUSSIAN - Gaussian Smoothing
          /BOX      - Box averaging
          /MASK     - Input is a mask cube instead of image.
   
    OUTPUTS:
          Result = the smoothed map
   
    OPTIONAL OUTPUTS:
          -
   
    PROCEDURES CALLED:
   
   
    EXAMPLE:
   
   
    MODIFICATION HISTORY:
                  10-09-2016: Added a mask expansion routine
          Written 29-12-2016 by P.Kamphuis
   
    NOTE:
   
   """

   n_params = 2
   
   mapor = map
   if bitwise_and(bitwise_and(bitwise_not((gaussian is not None)), bitwise_not((box is not None))), bitwise_not((mask is not None))):   
      gaussian = 1
   
   if (gaussian is not None):   
      xgauss = zeros([array(3 * 2.3548 * sigma, copy=0).astype(Int32)], Float32)
      ygauss = zeros([array(3 * 2.3548 * sigma, copy=0).astype(Int32)], Float32)
      xgauss = findgen(array(3 * 2.3548 * sigma, copy=0).astype(Int32)) - array(3 * 2.3548 * sigma / 2., copy=0).astype(Int32)
      ygauss = findgen(array(3 * 2.3548 * sigma, copy=0).astype(Int32)) - array(3 * 2.3548 * sigma / 2., copy=0).astype(Int32)
      gauss = zeros([array(3 * 2.3548 * sigma, copy=0).astype(Int32), array(3 * 2.3548 * sigma, copy=0).astype(Int32)], Float32)
      for i in arange(0, (array(3 * 2.3548 * sigma, copy=0).astype(Int32) - 1)+(1)):
      #From http://mathworld.wolfram.com/GaussianFunction.html
         gauss[:,i] = (1. / (2. * _sys_pi * sigma ** 2)) * exp(-((0.5 * xgauss[i] ** 2) / (sigma ** 2)) - ((0.5 * ygauss[:] ** 2) / (sigma ** 2)))
      tmp_smoothed_map = dblarr(array(3 * 2.3548 * sigma, copy=0).astype(Int32), array(3 * 2.3548 * sigma, copy=0).astype(Int32))
      startpixx = array(3 * 2.3548 * sigma / 2., copy=0).astype(Int32)
      endpixx = array(map[0,:], copy=0).nelements() - startpixx - 1
      startpixy = array(3 * 2.3548 * sigma / 2., copy=0).astype(Int32)
      endpixy = array(map[:,0], copy=0).nelements() - startpixy - 1
      for xax in arange(startpixx, (endpixx)+(1)):
         startx = xax - startpixx
         endx = xax + startpixx
         for yax in arange(startpixy, (endpixy)+(1)):
            starty = yax - startpixy
            endy = yax + startpixy
            map[yax,xax] = total(mapor[starty:(endy)+1,startx:(endx)+1] * gauss[:,:])
      return map
   if (box is not None):   
      inf = where(ravel(finite(mapor) == 0))[0]
      if inf[0] != -1:   
         mapor[inf] = 0.
      startpixx = array(sigma / 2., copy=0).astype(Int32)
      endpixx = array(map[0,:], copy=0).nelements() - startpixx - 1
      startpixy = array(sigma / 2., copy=0).astype(Int32)
      endpixy = array(map[:,0], copy=0).nelements() - startpixy - 1
      for xax in arange(startpixx, (endpixx)+(1)):
         startx = xax - startpixx
         endx = xax + startpixx
         for yax in arange(startpixy, (endpixy)+(1)):
            starty = yax - startpixy
            endy = yax + startpixy
            zero = where(ravel(mapor[starty:(endy)+1,startx:(endx)+1] == 0))[0]
            if zero[0] == -1:   
               map[yax,xax] = total(mapor[starty:(endy)+1,startx:(endx)+1]) / array(mapor[starty:(endy)+1,startx:(endx)+1], copy=0).nelements()
            else:   
               if total(mapor[starty:(endy)+1,startx:(endx)+1]) == 0.:   
                  map[yax,xax] = _sys_values.f_nan
               else:   
                  map[yax,xax] = total(mapor[starty:(endy)+1,startx:(endx)+1]) / (array(mapor[starty:(endy)+1,startx:(endx)+1], copy=0).nelements() - array(zero, copy=0).nelements())
      return map
   if (mask is not None):   
      
      for z in arange(0, (array(map[:,0,0], copy=0).nelements() - 1)+(1)):
         for x in arange(3, (array(map[0,0,:], copy=0).nelements() - 4)+(1)):
            for y in arange(3, (array(map[0,:,0], copy=0).nelements() - 4)+(1)):
               if mapor[z,y,x] == 1:   
                  map[z,y - sigma:(y + sigma)+1,x - sigma:(x + sigma)+1] = 1
      return map
   

