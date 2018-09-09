from numarray import *

def get_progress(name, ac, nopoints, loops, models):
   """
    NAME:
          GET_PROGRESS
   
    PURPOSE:
          This routine analyses the progress log from tirific and
          extracts the acceppted parameter, number of points, loops and
          total amount of models
   
    CATEGORY:
          Support
   
    CALLING SEQUENCE:
          GET_PROGRESS,name,AC,nopoints,loops,models
   
   
    INPUTS:
          name = name of the progress file
   
    OPTIONAL INPUTS:
          -
   
    KEYWORD PARAMETERS:
          -
   
    OUTPUTS:
          AC = value of the acceptance parameter
          nopoints = number of point sources in the model. 2D array for
          both sides
          loops = amount of completed big loops
          models = number of models produced
   
    OPTIONAL OUTPUTS:
          -
   
    PROCEDURES CALLED:
          STRUPCASE(), STRTRIM() STRCOMPRESS(), STR_SEP()
   
    EXAMPLE:
   
   
    MODIFICATION HISTORY:
          Written by P.Kamphuis 01-01-2015
   
    NOTE:
   
   """

   n_params = 5
   def _ret():  return (name, ac, nopoints, loops, models)
   
   # COMPILE_OPT IDL2
   ac = dblarr(1)
   nopoints = dblarr(1)
   openr(1, name)
   h = ' '
   #Just need the first line
   readf(1, h)
   tmp = str_sep(strtrim(strcompress(h), 2), ' ')
   close(1)
   #AC is the ninth item (In my modified Tirific)
   for i in arange(0, (array(tmp, copy=0).nelements() - 1)+(1)):
      tmp1 = str_sep(strtrim(strcompress(tmp[i]), 2), ':')
      if strupcase(tmp1[0]) == 'NP':   
         tmp2 = str_sep(strtrim(strcompress(tmp1[1]), 2), '/')
         nopoints = dblarr(array(tmp2, copy=0).nelements())
         nopoints = array(tmp2, copy=0).astype(Float64)
      if strupcase(tmp1[0]) == 'AC':   
         ac = array(tmp1[1], copy=0).astype(Float64)
      if strupcase(tmp1[0]) == 'BL':   
         loops = array(tmp1[1], copy=0).astype(Float64)
      if strupcase(tmp1[0]) == 'TM':   
         models = array(tmp1[1], copy=0).astype(Float64)
   
   return _ret()

