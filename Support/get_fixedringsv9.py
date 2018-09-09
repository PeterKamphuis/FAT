from numarray import *

def get_fixedringsv9(parametersin, rings, smooth=None, debug=None, warp_output=None, radii=None, sbr=None):
   """
    NAME:
          GET_FIXEDRINGSV9
   
    PURPOSE:
          Little routine to determine how many rings of the parameters
          should be considered as the flat inner part of the model based
          on the tiltogram of the current variations
   
    CATEGORY:
          Support
   
    CALLING SEQUENCE:
          GET_FIXEDRINGSV9,Parametersin,rings
   
   
    INPUTS:
          Parametersin = The INCL and PA parameters
   
    OPTIONAL INPUTS:
          -
   
    KEYWORD PARAMETERS:
          -
   
    OUTPUTS:
          rings = the amount of rings that should remain fixed
   
    OPTIONAL OUTPUTS:
          -
   
    PROCEDURES CALLED:
          ACOS(),ATAN(),MAX(),FINITE(),MAX() STDDEV(), ROBUST_SIGMA(), FLOOR()
   
    EXAMPLE:
   
   
    MODIFICATION HISTORY:
          06-05-2017 P.Kamphuis; Added the posibility to have the warp
                                 parameters written into a directory
                                 called Warp_Info
          18-02-2016 P.Kamphuis; Replaced sigma with STDDEV
          Written by P.Kamphuis 01-01-2015
   
    NOTE:
   
   """

   n_params = 2
   _opt = (smooth, debug, warp_output, radii, sbr)
   def _ret():
      _optrv = zip(_opt, [smooth, debug, warp_output, radii, sbr])
      _rv = [parametersin, rings]
      _rv += [_o[1] for _o in _optrv if _o[0] is not None]
      return tuple(_rv)
   
   # COMPILE_OPT IDL2
   parameters = parametersin
   if (warp_output is not None):   
      angle = dblarr(3, 2)
      tmp = where(ravel(sbr[0,:] > 1e-8))[0]
      if tmp[0] != -1:   
         outer1 = tmp[array(tmp, copy=0).nelements() - 1]
      else:   
         outer1 = array(tirresult[0,:], copy=0).nelements() - 1
      tmp = where(ravel(sbr[1,:] > 1e-8))[0]
      if tmp[0] != -1:   
         outer2 = tmp[array(tmp, copy=0).nelements() - 1]
      else:   
         outer2 = array(tirresult[0,:], copy=0).nelements() - 1
      max_index = concatenate([outer1, outer2])
   
   #First we need to calculate the normal
   #vectors for all rings on both sides
   
   #[PA,PA_2,INCL,INCL_2]
   pa1 = parameters[0,:]
   pa2 = parameters[1,:]
   incl1 = parameters[2,:]
   incl2 = parameters[3,:]
   x1 = dblarr(array(pa1, copy=0).nelements())
   y1 = dblarr(array(pa1, copy=0).nelements())
   z1 = dblarr(array(pa1, copy=0).nelements())
   x2 = dblarr(array(pa2, copy=0).nelements())
   y2 = dblarr(array(pa2, copy=0).nelements())
   z2 = dblarr(array(pa2, copy=0).nelements())
   add1 = dblarr(array(pa1, copy=0).nelements())
   add2 = dblarr(array(pa2, copy=0).nelements())
   for i in arange(0, (array(pa1, copy=0).nelements() - 1)+(1)):
      add1[i] = 0.
      while pa1[i] >= 90:
         pa1[i] = pa1[i] - 90
         add1[i] = add1[i] + 90.
      if pa1[i] == 0.:   
         pa1[i] = 0.00001
      theta = atan(tan(incl1[i] * (_sys_pi / 180.)) * tan(pa1[i] * (_sys_pi / 180.)))
      phi = atan(tan(pa1[i] * (_sys_pi / 180.)) / sin(theta))
      #So in cartesian coordinates
      x1[i] = sin(theta) * cos(phi)
      y1[i] = sin(theta) * sin(phi)
      z1[i] = cos(theta)
   for i in arange(0, (array(pa2, copy=0).nelements() - 1)+(1)):
      add2[i] = 0.
      while pa2[i] >= 90:
         pa2[i] = pa2[i] - 90
         add2[i] = add2[i] + 90.
      if pa2[i] == 0.:   
         pa2[i] = 0.00001
      theta = atan(tan(incl2[i] * (_sys_pi / 180.)) * tan(pa2[i] * (_sys_pi / 180.)))
      phi = atan(tan(pa2[i] * (_sys_pi / 180.)) / sin(theta))
      #So in cartesian coordinates
      x2[i] = sin(theta) * cos(phi)
      y2[i] = sin(theta) * sin(phi)
      z2[i] = cos(theta)
   #Next we calculate the tiltogram
   tiltogram1 = dblarr(array(pa1, copy=0).nelements(), array(pa1, copy=0).nelements())
   tiltogram2 = dblarr(array(pa2, copy=0).nelements(), array(pa2, copy=0).nelements())
   orings = findgen(array(pa1, copy=0).nelements()) + 0.5
   rrings = concatenate([findgen(array(pa1, copy=0).nelements() * 10.) / 10., 11])
   for i in arange(0, (array(pa1, copy=0).nelements() - 1)+(1)):
      for j in arange(0, (array(pa1, copy=0).nelements() - 1)+(1)):
         tiltogram1[j,i] = acos(x1[i] * x1[j] + y1[i] * y1[j] + z1[i] * z1[j]) / _sys_dtor
         if bitwise_not(finite(tiltogram1[j,i])):   
            tiltogram1[j,i] = 0.
   for i in arange(0, (array(pa2, copy=0).nelements() - 1)+(1)):
      for j in arange(0, (array(pa2, copy=0).nelements() - 1)+(1)):
         tiltogram2[j,i] = acos(x2[i] * x2[j] + y2[i] * y2[j] + z2[i] * z2[j]) / _sys_dtor
         if bitwise_not(finite(tiltogram2[j,i])):   
            tiltogram2[j,i] = 0.
   #If we run in debug mode we write fits files for the tiltograms
   if (warp_output is not None):   
      #  CD,workingdir,CURRENT=old_dir
      exist = file_test('Warp_Info', directory=True)
      if exist == 0:   
         spawn('mkdir Warp_Info')
      
      cd('Warp_Info')
      mkhdr(hed, tiltogram1)
      writefits('Receding_Tiltogram.fits', tiltogram1, hed)
      mkhdr(hed, tiltogram2)
      writefits('Approaching_Tiltogram.fits', tiltogram2, hed)
      cd('../')
   if (debug is not None):   
      mkhdr(hed, tiltogram1)
      writefits('test1.fits', tiltogram1, hed)
      mkhdr(hed, tiltogram2)
      writefits('test2.fits', tiltogram2, hed)
   tmp = dblarr(array(pa1, copy=0).nelements(), array(pa1, copy=0).nelements() * 10 + 1)
   rtiltogram1 = dblarr(array(pa1, copy=0).nelements() * 10. + 1, array(pa1, copy=0).nelements() * 10 + 1)
   rtiltogram2 = dblarr(array(pa2, copy=0).nelements() * 10. + 1, array(pa2, copy=0).nelements() * 10. + 1)
   for i in arange(0, (array(pa1, copy=0).nelements() - 1)+(1)):
      val = dblarr(array(tiltogram1[:,i], copy=0).nelements())
      val[:] = tiltogram1[:,i]
      output = 1
      interpolate(val, orings, output=output, newradii=rrings)
      tmp[:,i] = output[:]
   for i in arange(0, (array(tmp[:,0], copy=0).nelements() - 1)+(1)):
      val = tmp[i,:]
      interpolate(val, orings, output=output, newradii=rrings)
      rtiltogram1[i,:] = output[:]
   tmp = dblarr(array(pa2, copy=0).nelements(), array(pa2, copy=0).nelements() * 10 + 1)
   for i in arange(0, (array(pa2, copy=0).nelements() - 1)+(1)):
      val = dblarr(array(tiltogram2[:,i], copy=0).nelements())
      val[:] = tiltogram2[:,i]
      output = 1
      interpolate(val, orings, output=output, newradii=rrings)
      tmp[:,i] = output[:]
   for i in arange(0, (array(tmp[:,0], copy=0).nelements() - 1)+(1)):
      val = tmp[i,:]
      interpolate(val, orings, output=output, newradii=rrings)
      rtiltogram2[i,:] = output[:]
   if (debug is not None):   
      mkhdr(hed, rtiltogram1)
      writefits('tests1.fits', rtiltogram1, hed)
      mkhdr(hed, rtiltogram2)
      writefits('tests2.fits', rtiltogram2, hed)
      help(rtiltogram1)
   #now that we have the tiltograms we
   #need to calculate the average at each
   #radius. we want at 5 15 25 35 45
   thetain = dblarr(array(tiltogram1[0,:], copy=0).nelements())
   thetamut = dblarr(array(tiltogram1[0,:], copy=0).nelements())
   thetaout = dblarr(array(tiltogram1[0,:], copy=0).nelements())
   for i in arange(0, (array(tiltogram1[0,:], copy=0).nelements() - 1)+(1)):
      thetain[i] = mean(rtiltogram1[0:((i + 1) * 10)+1,0:((i + 1) * 10)+1])
      thetamut[i] = mean(rtiltogram1[0:((i + 1) * 10)+1,(i + 1) * 10:(array(rtiltogram1[0,:], copy=0).nelements() - 1)+1])
      thetaout[i] = mean(rtiltogram1[(i + 1) * 10:(array(rtiltogram1[0,:], copy=0).nelements() - 1)+1,(i + 1) * 10:(array(rtiltogram1[0,:], copy=0).nelements() - 1)+1])
   
   #And then we want to apply the rules
   #(i) the difference between thetain
   #and thetamut is larger than the differences observed at other radii
   #(ii) thetain < 5 deg
   # (iii) thetamut > 15 deg
   diff = abs(thetain[:] - thetamut[:])
   found = 0
   mxdif = 1
   if (debug is not None):   
      print 'This is the difference in disk1'
      print diff
   while bitwise_and(found == 0, mxdif != 0):
      mxdif = max(diff)
      rings = where(ravel(mxdif == diff))[0]
      #     rings=WHERE(diff GT 10)
      #     IF keyword_set(debug) then print,rings
      #We want the smallest ring that satifies our conditions
      if array(rings, copy=0).nelements() > 1:   
         rings = rings[0]
      if bitwise_and(thetain[rings] < 5, thetamut[rings] > 15):   
         found = 1
      else:   
         diff[rings] = 0
   if total(diff) == 0:   
      rings = array(diff, copy=0).nelements() - 1
   rings1 = rings
   if (debug is not None):   
      print 'This is the inner theta'
      print thetain
      print 'This is the edge theta'
      print thetamut
      print 'This is the outer theta'
      print thetaout
      print 'this is what we found in disk 1', rings1
   
   if (debug is not None):   
      mkhdr(hed, tiltogram1)
      writefits('test1.fits', tiltogram1, hed)
      mkhdr(hed, tiltogram2)
      writefits('test2.fits', tiltogram2, hed)
   if (warp_output is not None):   
      angle[0,:] = concatenate([tiltogram1[max_index[0],0], tiltogram2[max_index[1],0], (tiltogram1[max_index[0],0] + tiltogram2[max_index[1],0]) / 2.])
      tmp1 = max(tiltogram1[0:(max_index[0])+1,0])
      tmp2 = max(tiltogram2[0:(max_index[1])+1,0])
      tmp3 = max((tiltogram1[0:(max_index[0])+1,0] + tiltogram2[0:(max_index[1])+1,0]) / 2.)
      
      angle[1,:] = concatenate([tmp1[0], tmp2[0], tmp3[0]])
   
   tmp = dblarr(array(pa1, copy=0).nelements(), array(pa1, copy=0).nelements() * 10 + 1)
   rtiltogram1 = dblarr(array(pa1, copy=0).nelements() * 10. + 1, array(pa1, copy=0).nelements() * 10 + 1)
   rtiltogram2 = dblarr(array(pa2, copy=0).nelements() * 10. + 1, array(pa2, copy=0).nelements() * 10. + 1)
   for i in arange(0, (array(pa1, copy=0).nelements() - 1)+(1)):
      val = dblarr(array(tiltogram1[:,i], copy=0).nelements())
      val[:] = tiltogram1[:,i]
      output = 1
      interpolate(val, orings, output=output, newradii=rrings)
      tmp[:,i] = output[:]
   for i in arange(0, (array(tmp[:,0], copy=0).nelements() - 1)+(1)):
      val = tmp[i,:]
      interpolate(val, orings, output=output, newradii=rrings)
      rtiltogram1[i,:] = output[:]
   tmp = dblarr(array(pa2, copy=0).nelements(), array(pa2, copy=0).nelements() * 10 + 1)
   for i in arange(0, (array(pa2, copy=0).nelements() - 1)+(1)):
      val = dblarr(array(tiltogram2[:,i], copy=0).nelements())
      val[:] = tiltogram2[:,i]
      output = 1
      interpolate(val, orings, output=output, newradii=rrings)
      tmp[:,i] = output[:]
   for i in arange(0, (array(tmp[:,0], copy=0).nelements() - 1)+(1)):
      val = tmp[i,:]
      interpolate(val, orings, output=output, newradii=rrings)
      rtiltogram2[i,:] = output[:]
   if (debug is not None):   
      mkhdr(hed, rtiltogram1)
      writefits('tests1.fits', rtiltogram1, hed)
      mkhdr(hed, rtiltogram2)
      writefits('tests2.fits', rtiltogram2, hed)
      help(rtiltogram1)
   #now that we have the tiltograms we
   #need to calculate the average at each
   #radius. we want at 5 15 25 35 45
   thetain1 = dblarr(array(tiltogram1[0,:], copy=0).nelements())
   thetamut1 = dblarr(array(tiltogram1[0,:], copy=0).nelements())
   thetaout1 = dblarr(array(tiltogram1[0,:], copy=0).nelements())
   for i in arange(0, (array(tiltogram1[0,:], copy=0).nelements() - 1)+(1)):
      thetain1[i] = mean(rtiltogram1[0:((i + 1) * 10)+1,0:((i + 1) * 10)+1])
      thetamut1[i] = mean(rtiltogram1[0:((i + 1) * 10)+1,(i + 1) * 10:(array(rtiltogram1[0,:], copy=0).nelements() - 1)+1])
      thetaout1[i] = mean(rtiltogram1[(i + 1) * 10:(array(rtiltogram1[0,:], copy=0).nelements() - 1)+1,(i + 1) * 10:(array(rtiltogram1[0,:], copy=0).nelements() - 1)+1])
   
   #And then we want to apply the rules
   #(i) the difference between thetain
   #and thetamut is larger than the differences observed at other radii
   #(ii) thetain < 5 deg
   # (iii) thetamut > 15 deg
   diff = abs(thetain1[:] - thetamut1[:])
   found = 0
   mxdif = 1
   if (debug is not None):   
      print 'This is the difference in disk1'
      print diff
   while bitwise_and(found == 0, mxdif != 0):
      mxdif = max(diff)
      rings = where(ravel(mxdif == diff))[0]
      #     rings=WHERE(diff GT 10)
      #     IF keyword_set(debug) then print,rings
      #We want the smallest ring that satifies our conditions
      if array(rings, copy=0).nelements() > 1:   
         rings = rings[0]
      if bitwise_and(thetain1[rings] < 5, thetamut1[rings] > 15):   
         found = 1
      else:   
         diff[rings] = 0
   if total(diff) == 0:   
      rings = array(diff, copy=0).nelements() - 1
   rings1 = rings
   if (debug is not None):   
      print 'This is the inner theta'
      print thetain1
      print 'This is the edge theta'
      print thetamut1
      print 'This is the outer theta'
      print thetaout1
      print 'this is what we found in disk 1', rings1
   if (warp_output is not None):   
      if rings1 > max_index[0]:   
         rings1 = max_index[0]
      cd('Warp_Info', current=old_dir)
      openw(1, 'Warp_Info.txt')
      printf(1, 'Side', 'Inner Theta', 'Outer Theta', 'Warp Radius', 'Max Angle', 'Outer Angle', format='(6A20)')
      printf(1, 'Receding', thetain1[rings1], thetaout1[rings1], radii[rings1], angle[1,0], angle[0,0], format='(A20,5F20.3)')
      close(1)
      cd(old_dir)
   #now that we have the tiltograms we
   #need to calculate the average at each
   #radius. for the other side
   thetain = dblarr(array(tiltogram2[0,:], copy=0).nelements())
   thetamut = dblarr(array(tiltogram2[0,:], copy=0).nelements())
   thetaout = dblarr(array(tiltogram2[0,:], copy=0).nelements())
   for i in arange(0, (array(tiltogram2[0,:], copy=0).nelements() - 1)+(1)):
      thetain[i] = mean(rtiltogram2[0:((i + 1) * 10)+1,0:((i + 1) * 10)+1])
      thetamut[i] = mean(rtiltogram2[0:((i + 1) * 10)+1,(i + 1) * 10:(array(rtiltogram2[0,:], copy=0).nelements() - 1)+1])
      thetaout[i] = mean(rtiltogram2[(i + 1) * 10:(array(rtiltogram2[0,:], copy=0).nelements() - 1)+1,(i + 1) * 10:(array(rtiltogram2[0,:], copy=0).nelements() - 1)+1])
   
   #And then we want to apply the rules
   #to the otherside as well
   #(i) the difference between thetain
   #and thetamut is larger than the differences observed at other radii
   #(ii) thetain < 5 deg
   # (iii) thetamut > 15 deg
   diff = abs(thetain[:] - thetamut[:])
   found = 0
   mxdif = 1
   if (debug is not None):   
      print 'This is the difference in disk2'
      print diff
   while bitwise_and(found == 0, mxdif != 0):
      mxdif = max(diff)
      rings = where(ravel(mxdif == diff))[0]
      #    rings = WHERE(diff GT 10)
      if (debug is not None):   
         print rings
      if array(rings, copy=0).nelements() > 1:   
         rings = rings[0]
      if bitwise_and(thetain[rings] < 5, thetamut[rings] > 15):   
         found = 1
      else:   
         diff[rings] = 0
   if total(diff) == 0:   
      rings = array(diff, copy=0).nelements() - 1
   rings2 = rings
   
   if rings1 < rings2:   
      rings = rings1
   else:   
      rings = rings2
   if (debug is not None):   
      print 'This is the inner theta'
      print thetain
      print 'This is the edge theta'
      print thetamut
      print 'This is the outer theta'
      print thetaout
      print 'this is what we found in disk 2', rings2
   if (warp_output is not None):   
      if rings2 > max_index[1]:   
         rings2 = max_index[1]
      cd('Warp_Info', current=old_dir)
      openu(1, 'Warp_Info.txt', append=True)
      printf(1, 'Approaching', thetain[rings2], thetaout[rings2], radii[rings2], angle[1,1], angle[0,1], format='(A20,5F20.3)')
      printf(1, 'Average', (thetain1[rings1] + thetain[rings2]) / 2., (thetaout1[rings1] + thetaout[rings2]) / 2., (radii[rings1] + radii[rings2]) / 2., angle[1,2], angle[0,2], format='(A20,5F20.3)')
      close(1)
      cd(old_dir)
   if rings < 3:   
      rings = 3
   
   if bitwise_not((smooth is not None)):   
      if rings > array(array(parameters[0,:], copy=0).nelements() / 2., copy=0).astype(Int32):   
         if array(array(parameters[0,:], copy=0).nelements() / 2., copy=0).astype(Int32) >= 3:   
            rings = array(array(parameters[0,:], copy=0).nelements() / 2., copy=0).astype(Int32)
         else:   
            rings = 3
      if array(parameters[0,:], copy=0).nelements() > 12:   
         if rings < array(parameters[0,:], copy=0).nelements() / 4.:   
            ring = ceil(array(parameters[0,:], copy=0).nelements() / 4.)
   else:   
      if rings > array(array(parameters[0,:], copy=0).nelements() / 2., copy=0).astype(Int32):   
         if array(parameters[0,:], copy=0).nelements() > 10:   
            correction = 0.3
            while rings - array(correction * rings, copy=0).astype(Int32) < array(array(parameters[0,:], copy=0).nelements() / 2., copy=0).astype(Int32):
               correction = correction / 1.1
            rings = rings - array(correction * rings, copy=0).astype(Int32)
         else:   
            rings = rings - 1
   
   
   return _ret()

