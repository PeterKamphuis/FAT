from numarray import *

def set_warp_slopev3(sbrarr, sbrarr2, cutoff, inclinput2, painput2, inclinput3, painput3, norings, log=None, innerfix=None, sloped=None):
   """
    NAME:
          SET_WARP_SLOPEV3
   
    PURPOSE:
          This little routine takes the Surface brightness profiles and
          determines where they are below the cutoff values and sets
          those rings  to be only fitted with a slope in the Inclination and PA parameters
   
    CATEGORY:
          Support
   
    CALLING SEQUENCE:
          SET_WARP_SLOPEV3,SBRarr,SBRarr2,cutoff,INCLinput2,PAinput2,INCLinput3,PAinput3,norings,LOG=log,INNERFIX=innerfix
   
   
    INPUTS:
          SBRarr = approaching sbr array
          SBRarr2 = receding sbr array
          cutoff = cutoff values for the nodes
          norings = number of rings in model
          INCLinput2-3 = Template for fitting parameters for inclination
          PAinput2-3 = Template for fitting parameters for pa
   
    OPTIONAL INPUTS:
          LOG = name of the tracing log
          INNERFIX = number of rings that should be fitted as one in the
          central parameters
   
    KEYWORD PARAMETERS:
          -
   
    OUTPUTS:
          INCLinput2-3 = Updated fitting parameters for inclination
          PAinput2-3 = Updated fitting parameters for pa
   
   
    OPTIONAL OUTPUTS:
          -
   
    PROCEDURES CALLED:
          STR_SEP(), STRTRIM(), STRCOMPRESS(), STRING()
   
    EXAMPLE:
   
   
    MODIFICATION HISTORY:
          14-06-2016 P.Kamphuis; Added the rings that are sloped as output.
          Written 01-01-2015 P.Kamphuis v1.0
   
    NOTE:
   
   """

   n_params = 8
   warpconstused = sloped
   _opt = (log, innerfix, warpconstused)
   def _ret():
      _optrv = zip(_opt, [log, innerfix, warpconstused])
      _rv = [sbrarr, sbrarr2, cutoff, inclinput2, painput2, inclinput3, painput3, norings]
      _rv += [_o[1] for _o in _optrv if _o[0] is not None]
      return tuple(_rv)
   
   # COMPILE_OPT IDL2
   
   warpconstused = array(sbrarr, copy=0).nelements()
   if array(innerfix, copy=0).nelements() != 0:   
      inclinput2[0] = '!INCL ' + strtrim(strcompress(string(norings[0], format='(F7.4)')), 1) + ':' + strtrim(strcompress(string(innerfix + 1, format='(F7.4)')), 1)
      inclinput3[0] = '!INCL_2 ' + strtrim(strcompress(string(norings[0], format='(F7.4)')), 1) + ':' + strtrim(strcompress(string(innerfix + 1, format='(F7.4)')), 1)
      painput2[0] = '!PA ' + strtrim(strcompress(string(norings[0], format='(F7.4)')), 1) + ':' + strtrim(strcompress(string(innerfix + 1, format='(F7.4)')), 1)
      painput3[0] = '!PA_2 ' + strtrim(strcompress(string(norings[0], format='(F7.4)')), 1) + ':' + strtrim(strcompress(string(innerfix + 1, format='(F7.4)')), 1)
   else:   
      inclinput2[0] = '!INCL ' + strtrim(strcompress(string(norings[0], format='(F7.4)')), 1) + ':5'
      inclinput3[0] = '!INCL_2 ' + strtrim(strcompress(string(norings[0], format='(F7.4)')), 1) + ':5'
      painput2[0] = '!PA ' + strtrim(strcompress(string(norings[0], format='(F7.4)')), 1) + ':5'
      painput3[0] = '!PA_2 ' + strtrim(strcompress(string(norings[0], format='(F7.4)')), 1) + ':5'
   inclinput2[10] = ' '
   inclinput3[10] = ' '
   painput2[10] = ' '
   painput3[10] = ' '
   #Set warp slope of this fit
   #And we determine which part of the
   #warp should be a slope. i.e. when SBR
   #is lower than cutoff.
   sbrarr, sbrarr2, cutoff, warpconstused, True = get_newringsv9(sbrarr, sbrarr2, cutoff, warpconstused, individual=True)
   if warpconstused[0] < 6:   
      warpconstused[0] = 6
   if warpconstused[1] < 6:   
      warpconstused[1] = 6
   if bitwise_and(array(innerfix, copy=0).nelements() != 0, warpconstused[0] < innerfix + 2):   
      warpconstused[0] = innerfix + 2
   if bitwise_and(array(innerfix, copy=0).nelements() != 0, warpconstused[1] < innerfix + 2):   
      warpconstused[1] = innerfix + 2
   if warpconstused[0] < norings[0] - 1:   
      inclinput2[10] = 'INCL ' + strtrim(strcompress(string(norings[0] - 1, format='(I3)')), 1) + ':' + strtrim(strcompress(string(warpconstused[0] + 1, format='(I3)')), 1)
      painput2[10] = 'PA ' + strtrim(strcompress(string(norings[0] - 1, format='(I3)')), 1) + ':' + strtrim(strcompress(string(warpconstused[0] + 1, format='(I3)')), 1)
   if warpconstused[1] < norings[0] - 1:   
      inclinput3[10] = 'INCL_2 ' + strtrim(strcompress(string(norings[0] - 1, format='(I3)')), 1) + ':' + strtrim(strcompress(string(warpconstused[1] + 1, format='(I3)')), 1)
      painput3[10] = 'PA_2 ' + strtrim(strcompress(string(norings[0] - 1, format='(I3)')), 1) + ':' + strtrim(strcompress(string(warpconstused[1] + 1, format='(I3)')), 1)
   
   if size(log, type=True) == 7:   
      openu(66, log, append=True)
      printf(66, linenumber() + "SET_WARP_SLOPEV3: We fix the warp on side 1 from ring " + strtrim(string(warpconstused[0]), 2) + " on and side 2 from ring " + strtrim(string(warpconstused[1]), 2) + " on.")
      close(66)
   
   return _ret()

