from numarray import *

def sbr_check(template, templatevars, sbrarr, sbrarr2, cutoff):
   """
    NAME:
          SBR_CHECK
   
    PURPOSE:
          Routine to check that the sbr arrays are reasonable for the fit
   
    CATEGORY:
          Support
   
    CALLING SEQUENCE:
          SBR_CHECK,template,templatevars,sbrarr,sbrarr2,cutoff
   
   
    INPUTS:
          template = tirific template
          templatevars = template variable
          SBRarr = approaching sbrarr
          SBRarr2 = receding sbrarr
          cutoff = sbr cutoff values
   
    OPTIONAL INPUTS:
          -
   
    KEYWORD PARAMETERS:
          -
   
    OUTPUTS:
          template = modified tirific template
   
    OPTIONAL OUTPUTS:
          -
   
    PROCEDURES CALLED:
          -
   
    EXAMPLE:
   
   
    MODIFICATION HISTORY:
          Written 01-01-2015 P.Kamphuis v1.0
   
    NOTE:
   
   """

   n_params = 5
   def _ret():  return (template, templatevars, sbrarr, sbrarr2, cutoff)
   
   # COMPILE_OPT IDL2
   trigger = 0.
   if bitwise_and(bitwise_and(bitwise_and(sbrarr[0] > 2. * sbrarr[2], sbrarr2[0] > 2. * sbrarr2[2]), sbrarr[0] > 2. * sbrarr[1]), sbrarr2[0] > 2. * sbrarr2[1]):   
      if (sbrarr[2] + sbrarr2[2]) / 2. > cutoff[2]:   
         sbrarr[0] = (sbrarr[2] + sbrarr2[2]) / 2.
         sbrarr2[0] = (sbrarr[2] + sbrarr2[2]) / 2.
         sbrarr[1] = (sbrarr[1] + sbrarr2[1]) / 4.
         sbrarr2[1] = (sbrarr[1] + sbrarr2[1]) / 4.
      else:   
         sbrarr[0] = (sbrarr[2] + sbrarr2[2]) / 2.
         sbrarr[1] = 1.5 * cutoff[2]
         sbrarr[2] = 1.5 * cutoff[2]
         sbrarr2[0] = (sbrarr[2] + sbrarr2[2]) / 2.
         sbrarr2[1] = 1.5 * cutoff[2]
         sbrarr2[2] = 1.5 * cutoff[2]
      trigger = 1
   tmp = where(ravel(sbrarr < 0))[0]
   if tmp[0] != -1:   
      if tmp[0] == 0.:   
         sbrarr[0] = sbrarr[1]
      for i in arange(0, (array(tmp, copy=0).nelements() - 1)+(1)):
         if bitwise_and(tmp[i] != 0, tmp[i] != array(sbrarr, copy=0).nelements() - 1):   
            sbrarr[tmp[i]] = (sbrarr[tmp[i] - 1] + sbrarr[tmp[i] + 1]) / 2.
         if tmp[i] == 0:   
            sbrarr[tmp[i]] = sbrarr[tmp[i] + 1]
         if tmp[i] == array(sbrarr, copy=0).nelements() - 1:   
            sbrarr[tmp[i]] = sbrarr[tmp[i] - 1]
      trigger = 1
      
   tmp = where(ravel(sbrarr2 < 0))[0]
   if tmp[0] != -1:   
      if tmp[0] == 0.:   
         sbrarr2[0] = sbrarr2[1]
      for i in arange(0, (array(tmp, copy=0).nelements() - 1)+(1)):
         if bitwise_and(tmp[i] != 0, tmp[i] != array(sbrarr2, copy=0).nelements() - 1):   
            sbrarr2[tmp[i]] = (sbrarr2[tmp[i] - 1] + sbrarr2[tmp[i] + 1]) / 2.
         if tmp[i] == 0:   
            sbrarr2[tmp[i]] = sbrarr2[tmp[i] + 1]
         if tmp[i] == array(sbrarr2, copy=0).nelements() - 1:   
            sbrarr2[tmp[i]] = sbrarr2[tmp[i] - 1]
      trigger = 1
      
   if sbrarr[0] > sbrarr[1]:   
      sbrarr[0] = sbrarr[1]
      sbrarr2[0] = sbrarr[1]
      trigger = 1
      
   if trigger:   
      tmppos = where(ravel('SBR' == templatevars))[0]
      template[tmppos] = 'SBR=' + strjoin(strtrim(strcompress(string(sbrarr))), ' ')
      tmppos = where(ravel('SBR_2' == templatevars))[0]
      template[tmppos] = 'SBR_2=' + strjoin(strtrim(strcompress(string(sbrarr2))), ' ')
   
   return _ret()

