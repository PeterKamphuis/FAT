from numarray import *

def set_sbr(sbrinput1, sbrinput2, sbrinput3, sbrinput4, sbrinput5, sbrinput6, sbrarr, cutoff, norings, finishafter, log=None, initial=None, doubled=None):
   """
    NAME:
          SET_SBR
   
    PURPOSE:
          Routine to update the sbr fitting parameters based on the cutoff and the ring numbers
   
    CATEGORY:
          Support
   
    CALLING SEQUENCE:
          SET_SBR,SBRinput1,SBRinput2,SBRinput3,SBRinput4,SBRinput5,SBRinput6,SBRarr,cutoff,norings,finishafter,LOG=log,/INITIAL
   
   
    INPUTS:
          SBRinput1 = template fitting setting
          cutoff = cutoff values for the nodes
          norings = number of rings in model
          finishafter = indicator of type of fit
   
   
    OPTIONAL INPUTS:
          LOG = name of the tracing log
   
    KEYWORD PARAMETERS:
          /INITIAL - indicates it is the first time parameters are being set
   
    OUTPUTS:
          SBRinput1-6= Fitting parameters for different ranges of the model
   
    OPTIONAL OUTPUTS:
          -
   
    PROCEDURES CALLED:
          STR_SEP(), STRTRIM(), STRCOMPRESS(), STRING()
   
    EXAMPLE:
   
   
    MODIFICATION HISTORY:
          07-03-2017 P.Kamphuis; Added a condition that when noring[0]
                                 gt 15 half the minpar
          Written 01-01-2015 P.Kamphuis v1.0
   
    NOTE:
   
   """

   n_params = 10
   _opt = (log, initial, doubled)
   def _ret():
      _optrv = zip(_opt, [log, initial, doubled])
      _rv = [sbrinput1, sbrinput2, sbrinput3, sbrinput4, sbrinput5, sbrinput6, sbrarr, cutoff, norings, finishafter]
      _rv += [_o[1] for _o in _optrv if _o[0] is not None]
      return tuple(_rv)
   
   # COMPILE_OPT IDL2
   
   
   #first define some stepsizes based on
   #the cutoff values
   minstep = strtrim(strcompress(string(cutoff[norings[0]] / 2., format='(E8.1)')), 1)
   startstep = strtrim(strcompress(string(5. * cutoff[norings[0] - 1], format='(E8.1)')), 1)
   if (initial is not None):   
      satlevel = strtrim(strcompress(string(5. * array(startstep, copy=0).astype(Float64), format='(E8.1)')), 1)
   else:   
      satlevel = strtrim(strcompress(string(2. * array(startstep, copy=0).astype(Float64), format='(E8.1)')), 1)
   sbrinput1 = concatenate(['!SBR ' + strtrim(strcompress(string(norings[0], format='(F7.4)')), 1) + ':3', '1', strtrim(strcompress(string(cutoff[array(cutoff, copy=0).nelements() - 1], format='(E12.5)')), 1), startstep, minstep, satlevel, minstep, '3', '70', '70'])
   if norings[0] > 15.:   
      sbrinput1[2] = strtrim(strcompress(string(cutoff[norings[0]] / 2., format='(E8.1)')), 1)
   if doubled:   
      sbrinput1[2] = strtrim(strcompress(string(cutoff[norings[0]] / 1.2, format='(E8.1)')), 1)
   #Then copy the array into multiple arrays
   sbrinput2 = sbrinput1
   sbrinput3 = sbrinput1
   sbrinput4 = sbrinput1
   #and adjust them properly
   sbrinput2[0] = '!SBR_2 ' + strtrim(strcompress(string(norings[0], format='(F7.4)')), 1) + ':3'
   sbrinput3[0] = 'SBR  1 SBR_2 1 '
   sbrinput3[1] = strtrim(strcompress(string(sbrarr[1], format='(E12.5)')), 1)
   sbrinput4[0] = 'SBR 2 SBR_2 2 '
   #If the disk is too small
   if bitwise_or(norings[0] <= 4, finishafter == 1.1):   
      sbrinput1[2] = strtrim(strcompress(string(cutoff[array(cutoff, copy=0).nelements() - 1] / 4., format='(E12.5)')), 1)
      sbrinput2[2] = strtrim(strcompress(string(cutoff[array(cutoff, copy=0).nelements() - 1] / 4., format='(E12.5)')), 1)
      sbrinput3[2] = strtrim(strcompress(string(cutoff[array(cutoff, copy=0).nelements() - 1], format='(E12.5)')), 1)
      sbrinput4[2] = strtrim(strcompress(string(cutoff[array(cutoff, copy=0).nelements() - 1] * 2., format='(E12.5)')), 1)
   #IF the disk is large than treat inner and outer parts differently
   if bitwise_and(norings[0] > 6, bitwise_not((initial is not None))):   
      sbrinput5 = sbrinput1
      sbrinput6 = sbrinput2
      halfrings = floor((norings[0] - 2) / 2.)
      if halfrings <= 4:   
         sbrinput1[0] = '!SBR ' + strtrim(strcompress(string(norings[0] - 2, format='(F7.4)')), 1) + ':3'
         sbrinput1[2] = strtrim(strcompress(string(cutoff[norings[0] - 1], format='(E12.5)')), 1)
         sbrinput5[0] = '!SBR ' + strtrim(strcompress(string(norings[0], format='(F7.4)')), 1) + ':' + strtrim(strcompress(string(norings[0] - 1, format='(F7.4)')), 1)
         sbrinput5[2] = strtrim(strcompress(string(cutoff[norings[0] - 1] / 4., format='(E12.5)')), 1)
         sbrinput2[0] = '!SBR_2 ' + strtrim(strcompress(string(norings[0] - 2, format='(F7.4)')), 1) + ':3'
         sbrinput2[2] = strtrim(strcompress(string(cutoff[norings[0] - 1], format='(E12.5)')), 1)
         sbrinput6[0] = '!SBR_2 ' + strtrim(strcompress(string(norings[0], format='(F7.4)')), 1) + ':' + strtrim(strcompress(string(norings[0] - 1, format='(F7.4)')), 1)
         sbrinput6[2] = strtrim(strcompress(string(cutoff[norings[0] - 1] / 4., format='(E12.5)')), 1)
      else:   
         sbrinput1[0] = '!SBR ' + strtrim(strcompress(string(halfrings + 2, format='(F7.4)')), 1) + ':3'
         sbrinput1[2] = strtrim(strcompress(string(cutoff[norings[0] - 1], format='(E12.5)')), 1)
         sbrinput5[0] = '!SBR ' + strtrim(strcompress(string(norings[0], format='(F7.4)')), 1) + ':' + strtrim(strcompress(string(halfrings + 3, format='(F7.4)')), 1)
         sbrinput5[2] = strtrim(strcompress(string(cutoff[norings[0] - 1] / 4., format='(E12.5)')), 1)
         sbrinput2[0] = '!SBR_2 ' + strtrim(strcompress(string(halfrings + 2, format='(F7.4)')), 1) + ':3'
         sbrinput2[2] = strtrim(strcompress(string(cutoff[norings[0] - 1], format='(E12.5)')), 1)
         sbrinput6[0] = '!SBR_2 ' + strtrim(strcompress(string(norings[0], format='(F7.4)')), 1) + ':' + strtrim(strcompress(string(halfrings + 3, format='(F7.4)')), 1)
         sbrinput6[2] = strtrim(strcompress(string(cutoff[norings[0] - 1] / 4., format='(E12.5)')), 1)
   
   return _ret()

