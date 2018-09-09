from numarray import *

def isnumeric(input):
   """
    NAME:
          ISNUMERIC
   
    PURPOSE:
          Routine to determine wether a scalar is numeric
   
    CATEGORY:
          Support
   
    CALLING SEQUENCE:
          Result = isnumeric(input)
   
    INPUTS:
          input = the value to be checked
   
    OPTIONAL INPUTS:
          -
   
    KEYWORD PARAMETERS:
          -
   
    OUTPUTS:
          result= binary to check wether the value is numeric 0 is no 1
          is yes.
   
    OPTIONAL OUTPUTS:
          -
   
    PROCEDURES CALLED:
          -
   
    EXAMPLE:
   
   
    MODIFICATION HISTORY:
          15-11-2016 P.Kamphuis; Added full string length capability to
                                 check against numeric excepting '+-.ED'
                                 as otherwise a single digit is enough
                                 to make the string numeric.
          Written by ??????
   
    NOTE:
   
   """

   n_params = 1
   
   # COMPILE_OPT IDL2
   #if double returns an error then the string is not numeric
   # ON_IOERROR, FALSE
   #Get the length of the string
   tmp = strlen(input)
   #Keep a set of counters for numeric signs
   pluscount = 0
   minuscount = 0
   ecount = 0
   dcount = 0
   dotcount = 0
   #Loop through each element of the
   #string to check it is numeric
   for i in arange(0, (tmp - 1)+(1)):
      tmpcheck = strmid(input, i, 1)
      _expr = 1
      if _expr == (bitwise_and(strupcase(tmpcheck) == '+', pluscount < 2)):   
         pluscount += 1
      elif _expr == (bitwise_and(strupcase(tmpcheck) == '-', minuscount < 2)):   
         minuscount += 1
      elif _expr == (bitwise_and(strupcase(tmpcheck) == '.', dotcount < 2)):   
         dotcount += 1
      elif _expr == (bitwise_and(strupcase(tmpcheck) == 'E', ecount < 1)):   
         ecount += 1
      elif _expr == (bitwise_and(strupcase(tmpcheck) == 'D', dcount < 1)):   
         dcount += 1
      else:   
         test = array(tmpcheck, copy=0).astype(Float64)
      
   #Make sure the total is numeric
   test = array(input, copy=0).astype(Float64)
   #Return 1 on succes 0 on error
   return 1
   # false:
   return 0

