from numarray import *

def dec_names(values, tick_value=None, tick_name=None, n_ticks=None, increment=None, force=None):
#+
# NAME:
#   DEC_NAMES
# PURPOSE:
#   To generate the locations and names of appropriately places DEC
#   ticks on an axis of declination in a graph.
#
# CALLING SEQUENCE:
#   DEC_NAMES, values, TICK_VAL = TICK_VAL, TICK_NAME = TICK_NAME
#
# INPUTS:
#   VALUES -- The values of DEC along the axis.
#   FORCE -- force the given tick marks instead of converting to nice values
#
# KEYWORD PARAMETERS:
#   N_TICKS -- Forces labelling to this number of ticks.
#   INCREMENT -- Spacing between tick marks
#
# OUTPUTS:
#   TICK_VAL -- The value (in decimal degrees) of the location of the ticks.
#   TICK_NAME -- An array of strings naming the ticks in the
#                appropriate format.
# MODIFICATION HISTORY:
#       18-07-2011
#       Modified to correctly produce the degree symbol when using
#       True Type and postscript fonts
#       P.Kamphuis
#
#       15-09-2009
#       Modified to correctly handle southern declinations.
#       P.Kamphuis
#
#       Modified to correctly handle southern declinations.????
#       Mon Mar 3, 2003, Nate McCrady <nate@astro>
#
#       Written in a deeply fatigued state.
#       Wed Jul 24 02:13:37 2002, Erik Rosolowsky <eros@cosmic>
#
# NOTE: P. Kamphuis, GDL does not recognize the degree symbol not does it accept
#       superscripts in the tickmarks of plot. This leads to bad output in
#       the GDL plot.
#  -

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
   legit = concatenate([concatenate([60, 30, 15, 10, 5, 2, 1]), 1 / 60e0 * concatenate([30, 15, 10, 5, 2, 1]), 1 / 3.6e3 * concatenate([30, 15, 10, 5, 2, 1]), 1 / 3.6e4 * concatenate([5, 2, 1]), 1 / 3.6e5 * concatenate([5, 2, 1]), 1 / 3.6e6 * concatenate([5, 2, 1])])
   
   # Determine the spacing.  If number of ticks not specified, assume
   # between 3 and 6 ticks will do.
   if array(n_ticks, copy=0).nelements() > 0:   
      incr = range / n_ticks
      nticks = n_ticks + 2
      if bitwise_not((force is not None)):   
         # Choose ideal increment as largest increment smaller than required
         # increment.  Measure spacing logarithmically.
         diff = alog10(incr[0]) - alog10(legit)
         min_diff = array(diff[where(ravel(diff > 0))[0]], copy=0).min()
         ind = where(ravel(diff == min_diff))[0]
         increment = legit[ind]
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
         increment = array(legit[target / 4], copy=0).min() * incr / abs(incr)
      else:   
         message('This condition not yet implemented. I Suck.', con=True)
         return _ret()
   increment = increment[0]
   
   if increment >= 0:   
      tv = start + increment * (bitwise_and(start > 0, increment > 0)) - (start % increment)
      while tv[array(tv, copy=0).nelements() - 1] < finish:
         tv = concatenate([tv, tv[array(tv, copy=0).nelements() - 1] + increment])
   else:   
      tv = start + increment * (bitwise_and(start > 0, increment > 0)) - (start % increment)
      while tv[array(tv, copy=0).nelements() - 1] > finish:
         tv = concatenate([tv, tv[array(tv, copy=0).nelements() - 1] + increment])
   #  tv = dindgen(nticks)*increment+start+increment*(start mod increment gt 0)-(start mod increment)
   
   # Separate into nice Deg, Min, Sec.  This includes some
   # black magic to avoid irritating floor results with the DOUBLE values.
   
   deg = floor(tv) * (tv >= 0) + ceil(tv) * (tv < 0)
   minute = floor(array((tv - deg) * 60, copy=0).astype(Float32)) * (tv >= 0) + (tv < 0) * floor(array(abs((tv - deg)) * 60, copy=0).astype(Float32))
   sec = (array((array((tv - deg) * 60, copy=0).astype(Float64) - minute) * 60, copy=0).astype(Float64)) * (tv >= 0) + abs((array((array((tv - deg) * 60, copy=0).astype(Float64) + minute) * 60, copy=0).astype(Float64))) * (tv < 0)
   
   sec = array(sec, copy=0).astype(Float64)
   # Carrying
   
   c_sec = sec >= 60
   sec = sec - 60. * c_sec
   minute = minute + c_sec
   c_min = minute >= 60
   minute = minute - 60 * c_min
   deg = deg + c_min
   if abs(increment * 3600) < 1:   
      decs = ceil(-alog10(abs(increment) * 3600))
   else:   
      decs = 0
   
   a = string(047)
   b = string(042)
   # Now, the tedium of checking all the
   # formatting and only marking
   # differences.
   
   deg_name = strarr(array(tv, copy=0).nelements())
   min_name = deg_name
   sec_name = deg_name
   if _sys_p.font == -1:   
      degreesymbol = '!9%!X'
   else:   
      degreesymbol = '!9' + string(0260) + '!X'
   if tv[0] < 0.:   
      if array(tv[0], copy=0).astype(Int32) == 0:   
         deg_name[0] = '-' + strcompress(string(deg[0]), rem=True) + degreesymbol
      else:   
         deg_name[0] = strcompress(string(deg[0]), rem=True) + degreesymbol
   else:   
      deg_name[0] = strcompress(string(deg[0]), rem=True) + degreesymbol
   #  deg_name[0] = strcompress(string(deg[0]), /rem)+"!9%!X"
   min_name[0] = string(minute[0], format='(i2.2)') + a
   #  sec_name[0] = decimals(sec[0], decs)+b
   sec_name[0] = string(sec[0], format='(i2.2)') + "." + string(decs, format='(i1.1)') + b
   tracksec = 0.
   for i in arange(1, (array(tv, copy=0).nelements() - 1)+(1)):
      if deg[i] - deg[i - 1] != 0:   
         if tv[i] < 0.:   
            if array(tv[i], copy=0).astype(Int32) == 0:   
               deg_name[i] = '-' + strcompress(string(deg[i]), rem=True) + degreesymbol
            else:   
               deg_name[i] = strcompress(string(deg[i]), rem=True) + degreesymbol
         else:   
            deg_name[i] = strcompress(string(deg[i]), rem=True) + degreesymbol
      
      if minute[i] - minute[i - 1] != 0:   
         min_name[i] = string(minute[i], format='(i2.2)') + a
      if sec[i] - sec[i - 1] != 0:   
         sec_name[i] = string(sec[i], format='(i2.2)') + "." + string(decs, format='(i1.1)') + b
      tracksec = tracksec + array(string(sec[i], format='(i2.2)') + "." + string(decs, format='(i1.1)'), copy=0).astype(Float64)
   
   if total(sec) == 0:   
      sec_name = strarr(array(tv, copy=0).nelements())
   if total(minute) == 0:   
      min_name = strarr(array(tv, copy=0).nelements())
   if tracksec == 0.:   
      tn = deg_name + min_name
   else:   
      tn = deg_name + min_name + sec_name
   return _ret()

