from numarray import *

def ra_names(values, tick_value=None, tick_name=None, n_ticks=None, increment=None, force=None):
#+
# NAME:
#   RA_NAMES
# PURPOSE:
#   To generate the locations and names of appropriately places RA
#   ticks on an axis of declination in a graph.
#
# CALLING SEQUENCE:
#   RA_NAMES, values, TICK_VAL = TICK_VAL, TICK_NAME = TICK_NAME
#
# INPUTS:
#   VALUES -- The values of DEC along the axis.
#
# KEYWORD PARAMETERS:
#   N_TICKS -- Forces labelling to this number of ticks.
#   INCREMENT -- Spacing between tickmarks
#   FORCE -- force the given  number of ticks instead of
#            converting to nice values
#
# OUTPUTS:
#   TICK_VAL -- The value (in decimal degrees) of the location of the ticks.
#   TICK_NAME -- An array of strings naming the ticks in the
#                appropriate format.
# MODIFICATION HISTORY:
#       28-02-2016 P.Kamphuis; Replaced !U with !E for GDL
#       compatibility. Which didn't work.
#       august 2009 Decimals command unknown changed to string with decimals
#       Written in a deeply fatigued state.
#       Wed Jul 24 02:13:37 2002, Erik Rosolowsky <eros@cosmic>
#
# NOTE: P. Kamphuis, GDL does not accept superscripts in the tickmarks
#       of plot. This leads to bad output in the GDL plot.
#-

# Set up some information about the coordinate axis.
   n_params = 1
   tv = tick_value
   tn = tick_name
   _opt = (tv, tn, n_ticks, increment, force)
   def _ret():
      _optrv = zip(_opt, [tv, tn, n_ticks, increment, force])
      _rv = [values]
      _rv += [_o[1] for _o in _optrv if _o[0] is not None]
      return tuple(_rv)
   
   values = array(values, copy=0).astype(Float64)
   nelts = array(values, copy=0).nelements()
   start = values[0]
   finish = values[nelts - 1]
   range = finish - start
   
   # Establish the array of legitimate values for an increment.
   legit = concatenate([6, 4, 2, 1, 1 / 60e0 * concatenate([30, 15, 10, 5, 2, 1]), 1 / 3.6e3 * concatenate([30, 15, 10, 5, 2, 1]), 1 / 3.6e4 * concatenate([5, 2, 1]), 1 / 3.6e5 * concatenate([5, 2, 1]), 1 / 3.6e6 * concatenate([5, 2, 1])]) * 15
   
   # Determine the spacing.  If number of ticks not specified, assume
   # between 3 and 6 ticks will do.
   if array(n_ticks, copy=0).nelements() > 0:   
      incr = range / n_ticks
      nticks = n_ticks + 2
      if bitwise_not((force is not None)):   
         # Choose ideal increment as largest increment smaller than required
         # increment.  Measure spacing logarithmically.
         diff = alog10(abs(incr)) - alog10(legit)
         min_diff = array(diff[where(ravel(diff > 0))[0]], copy=0).min()
         ind = where(ravel(diff == min_diff))[0]
         increment = legit[ind] * incr / abs(incr)
      else:   
         increment = incr
   else:   
      # As above but choose number of ticks as well as the ideal increment.
      incr = range / (dindgen(4) + 3)
      legit_array = transpose(matrixmultiply(transpose((zeros([4], Float32) + 1)), transpose(legit)))
      incr_array = transpose(matrixmultiply(transpose(incr), transpose((zeros([array(legit, copy=0).nelements()], Float32) + 1))))
      diff = alog10(abs(incr_array)) - alog10(legit_array)
      
      ind = where(ravel(diff > 0))[0]
      
      if ct > 0:   
         min_diff = array(diff[ind], copy=0).min()
         target = where(ravel(diff == min_diff))[0]
         nticks = max((target % 4) + 3) + 2
         increment = array(legit[target / 4] * incr / abs(incr), copy=0).min()
      else:   
         message('This condition not yet implemented. I Suck.', con=True)
         return _ret()
   increment = increment[0]
   if increment >= 0:   
      tv = start + increment * (start > 0) - (start % increment)
      while tv[array(tv, copy=0).nelements() - 1] < finish:
         tv = concatenate([tv, tv[array(tv, copy=0).nelements() - 1] + increment])
   else:   
      tv = start - (start % increment)
      while tv[array(tv, copy=0).nelements() - 1] > finish:
         tv = concatenate([tv, tv[array(tv, copy=0).nelements() - 1] + increment])
   
   # Separate into nice Hour, Min, Sec.  This includes some
   # black magic to  avoid irritating floor results with the DOUBLE values.
   
   tv_hold = tv / 15
   ind = where(ravel(tv_hold < 0))[0]
   if ct > 0:   
      tv_hold[ind] = tv_hold[ind] + 24
   hour = floor(tv_hold)
   minute = floor(array((tv_hold - hour) * 60e0, copy=0).astype(Float64))
   sec = array(((array((tv_hold - hour) * 60e0, copy=0).astype(Float64) - minute) * 60), copy=0).astype(Float32)
   
   #Carrying
   c_sec = sec >= 60
   sec = sec - 60 * c_sec
   minute = minute + c_sec
   c_min = minute >= 60
   minute = minute - 60 * c_min
   hour = hour + c_min
   
   if abs(increment / 15 * 3600) < 1:   
      decs = ceil(-alog10(abs(increment / 15) * 3600))
   else:   
      decs = 0
   a = string(0140)
   b = string(042)
   # Now, the tedium of checking all the
   # formatting and only marking
   # differences.
   
   hour_name = strarr(array(tv, copy=0).nelements())
   min_name = hour_name
   sec_name = hour_name
   hour_name[0] = string(hour[0], format='(i2.2)') + "!Eh!N"
   min_name[0] = string(minute[0], format='(i2.2)') + "!Em!N"
   sec_name[0] = string(sec[0], format='(i2.2)') + "." + string(decs, format='(i1.1)') + "!Es!N"
   for i in arange(1, (array(tv, copy=0).nelements() - 1)+(1)):
      if hour[i] - hour[i - 1] != 0:   
         hour_name[i] = string(hour[i], format='(i2.2)') + "!Eh!N"
      if minute[i] - minute[i - 1] != 0:   
         min_name[i] = string(minute[i], format='(i2.2)') + "!Em!N"
      if sec[i] - sec[i - 1] != 0:   
         sec_name[i] = string(sec[i], format='(i2.2)') + "." + string(decs, format='(i1.1)') + "!Es!N"
      
   if total(sec) == 0:   
      sec_name = strarr(array(tv, copy=0).nelements())
   if total(minute) == 0:   
      min_name = strarr(array(tv, copy=0).nelements())
   if total(sec) == 0:   
      tn = hour_name + min_name
   else:   
      tn = hour_name + min_name + sec_name
   return _ret()


