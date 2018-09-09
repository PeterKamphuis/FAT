from numarray import *

def check_cflux(nopoints, norings, tirificfirst, tirificfirstvars, cfluxadjusted, log=None):
   """
    NAME:
          CHECK_CFLUX
   
    PURPOSE:
          Routine to adjust the cflux such that the number of points in
          the models is between 0.5E6 and 2.2E6
   
    CATEGORY:
          Support
   
    CALLING SEQUENCE:
          CHECK_CFLUX,nopoints,tirificfirst,tirificfirstvars,cfluxadjusted,log=log
   
   
    INPUTS:
          nopoints = The number of points in the model
          tirificfirst = String array with the template of a tirific def file
          tirificfirstvars = String array with the variables in the previous template
   
    OPTIONAL INPUTS:
          -
   
    KEYWORD PARAMETERS:
          log= general fitting log to document the action the routine has taken
   
    OUTPUTS:
          cfluxadjusted = a trigger that indicates whether the cflux is
          adjusted 0 no, 1 yes
   
    OPTIONAL OUTPUTS:
          -
   
    PROCEDURES CALLED:
          LINENUMBER(),STRTRIM(),STRCOMPRESS(),STRING(),STR_SEP()
   
    EXAMPLE:
   
   
    MODIFICATION HISTORY:
        06-03-2017 P. Kamphuis: Added a condition that allows for higher
        number of point sources for larger models. Additionally also
        check that nopoints is finite.
        23-07-2015 added a upper limit of 0.005 and a check that CFLUX
        is never set to 0.
        Written by P.Kamphuis 01-01-2015
   
    NOTE:
   """

   n_params = 5
   _opt = (log,)
   def _ret():
      _optrv = zip(_opt, [log])
      _rv = [nopoints, norings, tirificfirst, tirificfirstvars, cfluxadjusted]
      _rv += [_o[1] for _o in _optrv if _o[0] is not None]
      return tuple(_rv)
   
   # COMPILE_OPT IDL2
   cfluxadjusted = 0.
   
   for j in arange(0, (array(nopoints, copy=0).nelements() - 1)+(1)):
   
      if finite(nopoints[j]) == 0:   
         if size(log, type=True) == 7:   
            openu(66, log, append=True)
            printf(66, linenumber() + "CHECK_CFLUX: We had an infinite number of points somehow. We stop the code")
            close(66)
         stop()
      else:   
         if norings[0] < 15:   
            _expr = (1)
            if _expr == (nopoints[j] < 0.5e6):   
               if j == 0:   
                  tmppos = where(ravel('CFLUX' == tirificfirstvars))[0]
                  currentcflux = str_sep(strtrim(strcompress(tirificfirst[tmppos]), 2), '=')
                  if array(string(array(currentcflux[1], copy=0).astype(Float64) * nopoints[j] / 1e6), copy=0).astype(Float64) != 0.:   
                     tirificfirst[tmppos] = 'CFLUX= ' + string(array(currentcflux[1], copy=0).astype(Float64) * nopoints[j] / 1e6)
                  else:   
                     tirificfirst[tmppos] = 'CFLUX= 1e-5'
               else:   
                  tmppos = where(ravel('CFLUX_' + strtrim(strcompress(array(j + 1, copy=0).astype(Int32)), 2) == tirificfirstvars))[0]
                  currentcflux1 = str_sep(strtrim(strcompress(tirificfirst[tmppos]), 2), '=')
                  if array(string(array(currentcflux1[1], copy=0).astype(Float64) * nopoints[j] / 1e6), copy=0).astype(Float64) != 0.:   
                     tirificfirst[tmppos] = currentcflux1[0] + '= ' + string(array(currentcflux1[1], copy=0).astype(Float64) * nopoints[j] / 1e6)
                  else:   
                     tirificfirst[tmppos] = currentcflux1[0] + '= 1e-5'
               if j == array(nopoints, copy=0).nelements() - 1:   
                  if array(currentcflux, copy=0).nelements() == 0:   
                     tmppos = where(ravel('CFLUX' == tirificfirstvars))[0]
                     currentcflux = str_sep(strtrim(strcompress(tirificfirst[tmppos]), 2), '=')
                  
                  if size(log, type=True) == 7:   
                     openu(66, log, append=True)
                     printf(66, linenumber() + "CHECK_CFLUX: CFLUX is scaled down to" + strjoin(concatenate([string(array(currentcflux[1], copy=0).astype(Float64) * nopoints[0] / 1e6), string(array(currentcflux1[1], copy=0).astype(Float64) * nopoints[1] / 1e6)]), ' and '))
                     close(66)
               cfluxadjusted = 1.
            elif _expr == (nopoints[j] > 2.2e6):   
               if j == 0:   
                  tmppos = where(ravel('CFLUX' == tirificfirstvars))[0]
                  currentcflux = str_sep(strtrim(strcompress(tirificfirst[tmppos]), 2), '=')
                  if array(string(array(currentcflux[1], copy=0).astype(Float64) * nopoints[j] / 1e6), copy=0).astype(Float64) < 0.005:   
                     tirificfirst[tmppos] = 'CFLUX= ' + string(array(currentcflux[1], copy=0).astype(Float64) * nopoints[j] / 1e6)
                  else:   
                     tirificfirst[tmppos] = 'CFLUX= 0.005'
               else:   
                  tmppos = where(ravel('CFLUX_' + strtrim(strcompress(array(j + 1, copy=0).astype(Int32)), 2) == tirificfirstvars))[0]
                  currentcflux1 = str_sep(strtrim(strcompress(tirificfirst[tmppos]), 2), '=')
                  if array(string(array(currentcflux1[1], copy=0).astype(Float64) * nopoints[j] / 1e6), copy=0).astype(Float64) < 0.005:   
                     tirificfirst[tmppos] = currentcflux1[0] + '= ' + string(array(currentcflux1[1], copy=0).astype(Float64) * nopoints[j] / 1e6)
                  else:   
                     tirificfirst[tmppos] = currentcflux1[0] + '= 0.005'
               if j == array(nopoints, copy=0).nelements() - 1:   
                  if array(currentcflux, copy=0).nelements() == 0:   
                     tmppos = where(ravel('CFLUX' == tirificfirstvars))[0]
                     currentcflux = str_sep(strtrim(strcompress(tirificfirst[tmppos]), 2), '=')
                  if size(log, type=True) == 7:   
                     openu(66, log, append=True)
                     printf(66, linenumber() + "CHECK_CFLUX: CFLUX is scaled up to" + strjoin(concatenate([string(array(currentcflux[1], copy=0).astype(Float64) * nopoints[0] / 1e6), string(array(currentcflux1[1], copy=0).astype(Float64) * nopoints[1] / 1e6)]), ' and '))
                     close(66)
            else:   
               if j == array(nopoints, copy=0).nelements() - 1:   
                  if size(log, type=True) == 7:   
                     openu(66, log, append=True)
                     printf(66, linenumber() + "CHECK_CFLUX: CFLUX is ok.")
                     close(66)
            
         else:   
            factor = (norings[0] / 15.) ** 1.5
            _expr = (1)
            if _expr == (nopoints[j] < 0.5e6 * factor):   
               if j == 0:   
                  tmppos = where(ravel('CFLUX' == tirificfirstvars))[0]
                  currentcflux = str_sep(strtrim(strcompress(tirificfirst[tmppos]), 2), '=')
                  if array(string(array(currentcflux[1], copy=0).astype(Float64) * nopoints[j] / (1e6 * factor)), copy=0).astype(Float64) != 0.:   
                     tirificfirst[tmppos] = 'CFLUX= ' + string(array(currentcflux[1], copy=0).astype(Float64) * nopoints[j] / (1e6 * factor))
                  else:   
                     tirificfirst[tmppos] = 'CFLUX= 1e-5'
               else:   
                  tmppos = where(ravel('CFLUX_' + strtrim(strcompress(array(j + 1, copy=0).astype(Int32)), 2) == tirificfirstvars))[0]
                  currentcflux1 = str_sep(strtrim(strcompress(tirificfirst[tmppos]), 2), '=')
                  if array(string(array(currentcflux1[1], copy=0).astype(Float64) * nopoints[j] / (1e6 * factor)), copy=0).astype(Float64) != 0.:   
                     tirificfirst[tmppos] = currentcflux1[0] + '= ' + string(array(currentcflux1[1], copy=0).astype(Float64) * nopoints[j] / (1e6 * factor))
                  else:   
                     tirificfirst[tmppos] = currentcflux1[0] + '= 1e-5'
               if j == array(nopoints, copy=0).nelements() - 1:   
                  if array(currentcflux, copy=0).nelements() == 0:   
                     tmppos = where(ravel('CFLUX' == tirificfirstvars))[0]
                     currentcflux = str_sep(strtrim(strcompress(tirificfirst[tmppos]), 2), '=')
                  if size(log, type=True) == 7:   
                     openu(66, log, append=True)
                     printf(66, linenumber() + "CHECK_CFLUX: CFLUX is scaled down to" + strjoin(concatenate([string(array(currentcflux[1], copy=0).astype(Float64) * nopoints[0] / (1e6 * factor)), string(array(currentcflux1[1], copy=0).astype(Float64) * nopoints[1] / (1e6 * factor))]), ' and '))
                     close(66)
               cfluxadjusted = 1.
            elif _expr == (nopoints[j] > 2.2e6 * factor):   
               if j == 0:   
                  tmppos = where(ravel('CFLUX' == tirificfirstvars))[0]
                  currentcflux = str_sep(strtrim(strcompress(tirificfirst[tmppos]), 2), '=')
                  if array(string(array(currentcflux[1], copy=0).astype(Float64) * nopoints[j] / (1e6 * factor)), copy=0).astype(Float64) < 0.005:   
                     tirificfirst[tmppos] = 'CFLUX= ' + string(array(currentcflux[1], copy=0).astype(Float64) * nopoints[j] / (1e6 * factor))
                  else:   
                     tirificfirst[tmppos] = 'CFLUX= 0.005'
               else:   
                  tmppos = where(ravel('CFLUX_' + strtrim(strcompress(array(j + 1, copy=0).astype(Int32)), 2) == tirificfirstvars))[0]
                  currentcflux1 = str_sep(strtrim(strcompress(tirificfirst[tmppos]), 2), '=')
                  if array(string(array(currentcflux1[1], copy=0).astype(Float64) * nopoints[j] / (1e6 * factor)), copy=0).astype(Float64) < 0.005:   
                     tirificfirst[tmppos] = currentcflux1[0] + '= ' + string(array(currentcflux1[1], copy=0).astype(Float64) * nopoints[j] / (1e6 * factor))
                  else:   
                     tirificfirst[tmppos] = currentcflux1[0] + '= 0.005'
               if j == array(nopoints, copy=0).nelements() - 1:   
                  if array(currentcflux, copy=0).nelements() == 0:   
                     tmppos = where(ravel('CFLUX' == tirificfirstvars))[0]
                     currentcflux = str_sep(strtrim(strcompress(tirificfirst[tmppos]), 2), '=')
                  if size(log, type=True) == 7:   
                     openu(66, log, append=True)
                     printf(66, linenumber() + "CHECK_CFLUX: CFLUX is scaled up to" + strjoin(concatenate([string(array(currentcflux[1], copy=0).astype(Float64) * nopoints[0] / (1e6 * factor)), string(array(currentcflux1[1], copy=0).astype(Float64) * nopoints[1] / (1e6 * factor))]), ' and '))
                     close(66)
            else:   
               if j == array(nopoints, copy=0).nelements() - 1:   
                  if size(log, type=True) == 7:   
                     openu(66, log, append=True)
                     printf(66, linenumber() + "CHECK_CFLUX: CFLUX is ok.")
                     close(66)
            
   
   return _ret()

