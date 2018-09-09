from numarray import *

def linenumber():
   """
    NAME:
          LINENUMBER
   
    PURPOSE:
          A function to print the linenumber in the routine directly below main
   
    CATEGORY:
          Support
   
    CALLING SEQUENCE:
          result=LINENUMBER()
   
   
    INPUTS:
          -
   
    OPTIONAL INPUTS:
          -
   
    KEYWORD PARAMETERS:
          -
   
    OUTPUTS:
          result = the line number of the program that is being run. If
          in main it will return (You are in the main)
   
    OPTIONAL OUTPUTS:
          -
   
    PROCEDURES CALLED:
          STR_SEP(), STRTRIM(), STRCOMPRESS()
   
    EXAMPLE:
   
   
    MODIFICATION HISTORY:
          Written by P.Kamphuis 01-01-2015
   
    NOTE:
   
   """

   n_params = 0
   
   help(calls=cc)		# Get list of all calls.
   tmp = str_sep(strtrim(strcompress(cc[array(cc, copy=0).nelements() - 2]), 2), '(')
   if array(tmp, copy=0).nelements() > 1:   
      tmp3 = str_sep(strtrim(strcompress(tmp[0]), 2), ' ')
      if tmp3[0] != 'LINENUMBER':   
         tmp2 = str_sep(strtrim(strcompress(tmp[1]), 2), ')')
         linenumber = tmp2[0]
      else:   
         linenumber = 'You are in the main'
   else:   
      linenumber = 'You are in the main'
   return ' (' + linenumber + ') '
   

