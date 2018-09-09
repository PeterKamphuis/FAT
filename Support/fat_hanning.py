from numarray import *

def fat_hanning(sbrin, rings=None):
   """
    NAME:
          FAT_HANNING
   
    PURPOSE:
          Routine to apply a hanning smoothing to the SBR Profile
   
    CATEGORY:
          Support
   
    CALLING SEQUENCE:
          Result = FAT_HANNING(Profile)
   
   
    INPUTS:
         SBRin = the profile to be smoothed
         rings = the amount of rings that are valid
    OPTIONAL INPUTS:
          -
   
    KEYWORD PARAMETERS:
   
    OUTPUTS:
          Result = the smoothed profile
   
    OPTIONAL OUTPUTS:
          -
   
    PROCEDURES CALLED:
   
   
    EXAMPLE:
   
   
    MODIFICATION HISTORY:
          Written 16-06-2017 by P.Kamphuis, S. Kurapati
   
    NOTE:
   
   """

   n_params = 1
   
   # COMPILE_OPT IDL2
   defsysv('!GDL', exists=gdlidl) #is 1 when running GDL
   #  goto,skipcatch
   catch(error_status)
   if error_status != 0.:   
      
      print ' '
      print 'Oops the following went wrong in FAT_HANNING:'
      print '-------------------------------------------'
      help(last_message=True, output=theerrmess)
      #This gives trace back information in
      #IDL but unfortunately not in GDL
      for j in arange(0, (array(theerrmess, copy=0).nelements() - 1)+(1)):
         print theerrmess[j]
      print '-------------------------------------------'
      if gdlidl:   
         print ' '
         print 'Unfortunately GDL does not provide trace back information in help'
         print 'In order to fix this error please open FAT.pro and '
         print 'uncomment line 51 (goto,skipcatch).'
         print ' '
         print 'Then reproduce the error with trace back information '
         print 'and submit an issue at:'
         print ' '
         print '       https://github.com/PeterKamphuis/FAT/issues'
         print ' '
         print 'If the error occured while fitting a galaxy, please attach you fitting'
         print 'log as well.'
      else:   
         print ' '
         print 'If this message completely baffles you then please submit the last 20 lines of output as a bug report to:'
         print ' '
         print '       https://github.com/PeterKamphuis/FAT/issues'
         print ' '
         print 'If the error occured while fitting a galaxy, please attach you fitting'
         print 'log as well.'
      stop()
   # skipcatch:
   sbr = sbrin
   #If we provide a number of valid rings
   #then we want to set all rings outside
   #that to 0.
   #  IF n_elements(rings) GT 0 then begin
   #     IF rings[0] LT n_elements(SBR) then SBR[rings[0]-1:n_elements(SBR)-1]=0.
   #  ENDIF
   
   sbrout = dblarr(array(sbr, copy=0).nelements())
   _expr = 1
   if _expr == (array(sbr, copy=0).nelements() <= 3):   
      return sbr
   elif _expr == (array(sbr, copy=0).nelements() < 15):   
      points = 5
   else:   
      points = 7
   
   window = dblarr(points)
   xaxis = findgen(points)
   window = 0.5 * (1 - cos(_sys_pi * 2. * xaxis[:] / (points - 1)))
   window = window[:] / total(window)
   for i in arange(array((points - 1) / 2., copy=0).astype(Int32), (array(sbr, copy=0).nelements() - (array((points - 1) / 2., copy=0).astype(Int32) + 1))+(1)):
      sbrout[i] = total(sbr[i - array((points - 1) / 2., copy=0).astype(Int32):(i + array((points - 1) / 2., copy=0).astype(Int32))+1] * window)
   #We'll just use a truncated edge in the center
   for i in arange(0, (array((points - 1) / 2., copy=0).astype(Int32) - 1)+(1)):
      sbrout[i] = total(sbr[0:(i + array((points - 1) / 2., copy=0).astype(Int32))+1] * window[array((points - 1) / 2., copy=0).astype(Int32) - i:(points - 1)+1])
   #at the outer edge we can use
   #zer0-padding
   for i in arange(array(sbr, copy=0).nelements() - array((points - 1) / 2., copy=0).astype(Int32), (array(sbr, copy=0).nelements() - 1)+(1)):
      sbrout[i] = total(concatenate([sbr[i - array((points - 1) / 2., copy=0).astype(Int32):(array(sbr, copy=0).nelements() - 1)+1], (0.)*ones([i - (array(sbr, copy=0).nelements() - (array((points - 1) / 2., copy=0).astype(Int32) + 1))])]) * window)
   if array(rings, copy=0).nelements() > 0:   
      if rings[0] + 1 < array(sbr, copy=0).nelements():   
         sbrout[rings[0] + 1:(array(sbrout, copy=0).nelements() - 1)+1] = 1e-16
   tmp = where(ravel(sbrout < 1e-8))[0]
   if tmp[0] != -1:   
      sbrout[tmp] = 1e-16
   
   return sbrout

