from numarray import *

def obtain_pav2(map, pa, inclination=None, center=None, noise=None, iterations=None):
   """
    NAME:
          OBTAIN_PAV2
   
    PURPOSE:
          Routine to fit five ellipses to the regions defined by
          the ten area's between the maximum and the noise/minimum
   
    CATEGORY:
          Support
   
    CALLING SEQUENCE:
          OBTAIN_PAV2,map,PA,INCLINATION=inclination,CENTER=center,NOISE=noise,ITERATIONS=norings
   
   
    INPUTS:
          map = moment 0 map of the galaxy
   
    OPTIONAL INPUTS:
          -
   
    KEYWORD PARAMETERS:
          INCLINATION = first guess at the inclination from the
          ellipse's axis ratios.
          CENTER = center of the galaxy in pixels
          NOISE = estimate of the noise of the moment 0 map
          ITERATIONS = The amount of rings to be fitted, the default is 5.
   
    OUTPUTS:
          PA = the positional angle
   
    OPTIONAL OUTPUTS:
          -
   
    PROCEDURES CALLED:
          FIT_ELLIPSE()
   
    EXAMPLE:
   
   
    MODIFICATION HISTORY:
          16-03-2017 P.Kamphuis; Fit ellipse crashes when only one index
                                 is presented to it hence the condition has been set for more
                                 than one index to be present
          18-02-2016 P.Kamphuis; Replaced sigma with STDDEV
          Written 01-01-2015 P.Kamphuis v1.0
   
    NOTE:
   
   """

   n_params = 2
   norings = iterations
   _opt = (inclination, center, noise, norings)
   def _ret():
      _optrv = zip(_opt, [inclination, center, noise, norings])
      _rv = [map, pa]
      _rv += [_o[1] for _o in _optrv if _o[0] is not None]
      return tuple(_rv)
   
   # COMPILE_OPT IDL2
   tmp0 = where(ravel(map > 0.))[0]
   mapmax = max(map[tmp0], min=mapmin)
   if array(noise, copy=0).nelements() < 1:   
      noise = mapmin
   if array(norings, copy=0).nelements() == 0:   
      norings = 5.
   maxint = mapmax - mapmax / 2.
   step = (maxint - (3. * noise)) / norings
   pafound = dblarr(norings)
   inclfound = dblarr(norings)
   for i in arange(0, (norings - 1)+(1)):
      minint = maxint - step
      tmp = where(ravel(bitwise_and(map > minint, map < maxint)))[0]
      #fit_ellipse does not work when we have only one element
      if array(tmp, copy=0).nelements() > 1:   
         if array(center, copy=0).nelements() == 0:   
            ell = fit_ellipse(tmp, xsize=array(map[0,:], copy=0).nelements(), ysize=array(map[:,0], copy=0).nelements(), orientation=paring, axes=inclring)
            pafound[i] = paring
            ratio = array(inclring[1] / inclring[0], copy=0).astype(Float64)
            if ratio > 1.:   
               ratio = 1.
            if ratio < 0.21:   
               ratio = 0.21
            inclfound[i] = acos(sqrt((ratio ** 2 - 0.2 ** 2) / 0.96)) * 180. / _sys_pi + 2.
         else:   
            ell = fit_ellipse(tmp, center=center, xsize=array(map[0,:], copy=0).nelements(), ysize=array(map[:,0], copy=0).nelements(), orientation=paring, axes=inclring)
            pafound[i] = paring
            ratio = array(inclring[1] / inclring[0], copy=0).astype(Float64)
            if ratio > 1.:   
               ratio = 1.
            if ratio < 0.21:   
               ratio = 0.21
            inclfound[i] = acos(sqrt((ratio ** 2 - 0.2 ** 2) / 0.96)) * 180. / _sys_pi + 2.
      maxint = minint
   tmp = where(ravel(pafound < 0))[0]
   if tmp[0] != -1:   
      pafound[tmp] = pafound[tmp] + 360
   for i in arange(1, (array(pafound, copy=0).nelements() - 1)+(1)):
      diff = pafound[i - 1] - pafound[i]
      if diff > 120.:   
         pafound[i] = pafound[i] + 180
      if diff < -120:   
         pafound[i] = pafound[i] - 180
   tmp = where(ravel(pafound != 0.))[0]
   pa = dblarr(2)
   inclination = dblarr(2)
   if tmp[0] != -1:   
      
      pa[0] = total(pafound[tmp]) / array(pafound[tmp], copy=0).nelements() + 90.
      if pa[0] < 0:   
         pa[0] = pa[0] + 360
      if pa[0] > 360:   
         pa[0] = pa[0] - 360
      pa[1] = stddev(pafound[tmp])
      
      tmp2 = inclfound[2:(array(inclfound, copy=0).nelements() - 1)+1]
      tmp = where(ravel(tmp2 != 0))[0]
      if tmp[0] != -1:   
         inclination[0] = total(tmp2[tmp]) / (array(tmp, copy=0).nelements())
         inclination[1] = stddev(tmp2[tmp])
      
   if inclination[0] == 0.:   
      inclination[1] = 45.
      inclination[0] = 45.
   
   if pa[0] == 0.:   
      pa[1] = 180.
      pa[0] = 90.
   
   
   return _ret()

