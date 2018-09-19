def writefittingvariables(inputarray, v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12, v13, v14, v15, v16, v17, v18, v19, v20, v21, v22, v23, v24, v25, v26, v27, v28, v29, v30, v31, v32, v33, v34, v35, v36, v37, v38, v39, v40):
   """
    NAME:
          WRITEFITTINGVARIABLES
   
    PURPOSE:
          Routine to write the fitting variables to the tirific template
          array
   
    CATEGORY:
          Support
   
    CALLING SEQUENCE:
          WRITEFITTINGVARIABLES,inputarray,v1[,v2,v3,...v40]
   
   
    INPUTS:
          inputarray = tirific template
          v1-v40 = Arrays with the various fitting parameters that
          should be written to the template
   
    OPTIONAL INPUTS:
          -
   
    KEYWORD PARAMETERS:
          -
   
    OUTPUTS:
          inputarray = updated tirific template
   
    OPTIONAL OUTPUTS:
          -
   
    PROCEDURES CALLED:
          STR_SEP(), STRTRIM(), STRCOMPRESS(), STRING()
   
    EXAMPLE:
   
   
    MODIFICATION HISTORY:
          Written 01-01-2015 P.Kamphuis v1.0
   
    NOTE:
         Fitting VARIABLE input order= [Vary, Parmax,parmin, delstart,delend,satdelt,mindelt,moderate,ITESTART,ITEEND,VARINDX]
   """

   n_params = 41
   def _ret():  return (inputarray)
   
   # COMPILE_OPT IDL2
   # ON_ERROR, 2                    #Return to caller
   #We need at least three input parameters
   if n_params < 2:   
      print 'WRITEFITTINGVARIABLES: You need to change at least one variable'
      return _ret()
   
   nvar = n_params - 1         #Number of variables to deal with
   #Let's check wether the amount of limits and variables match
   
   #let's build  some arrays from the input
   strings = strarr(11)
   strings[0] = 'VARY= '
   strings[1] = 'PARMAX= '
   strings[2] = 'PARMIN= '
   strings[3] = 'MODERATE= '
   strings[4] = 'DELSTART= '
   strings[5] = 'DELEND= '
   strings[6] = 'ITESTART= '
   strings[7] = 'ITEEND= '
   strings[8] = 'SATDELT= '
   strings[9] = 'MINDELTA= '
   strings[10] = 'VARINDX= '
   vv = 'v' + strtrim(indgen(nvar) + 1, 2)
   
   
   for i in arange(0, (nvar - 1)+(1)):
      limits = scope_varfetch(vv[i], level=0)
      if array(limits, copy=0).nelements() == 10:   
         strings[0] = strings[0] + strtrim(strcompress(string(limits[0]))) + ','
         strings[1] = strings[1] + strtrim(strcompress(string(limits[1]))) + ' '
         strings[2] = strings[2] + strtrim(strcompress(string(limits[2]))) + ' '
         strings[3] = strings[3] + strtrim(strcompress(string(limits[7]))) + ' '
         strings[4] = strings[4] + strtrim(strcompress(string(limits[3]))) + ' '
         strings[5] = strings[5] + strtrim(strcompress(string(limits[4]))) + ' '
         strings[6] = strings[6] + strtrim(strcompress(string(limits[8]))) + ' '
         strings[7] = strings[7] + strtrim(strcompress(string(limits[9]))) + ' '
         strings[8] = strings[8] + strtrim(strcompress(string(limits[5]))) + ' '
         strings[9] = strings[9] + strtrim(strcompress(string(limits[6]))) + ' '
      if array(limits, copy=0).nelements() == 11:   
         strings[0] = strings[0] + strtrim(strcompress(string(limits[0]))) + ','
         strings[1] = strings[1] + strtrim(strcompress(string(limits[1]))) + ' '
         strings[2] = strings[2] + strtrim(strcompress(string(limits[2]))) + ' '
         strings[3] = strings[3] + strtrim(strcompress(string(limits[7]))) + ' '
         strings[4] = strings[4] + strtrim(strcompress(string(limits[3]))) + ' '
         strings[5] = strings[5] + strtrim(strcompress(string(limits[4]))) + ' '
         strings[6] = strings[6] + strtrim(strcompress(string(limits[8]))) + ' '
         strings[7] = strings[7] + strtrim(strcompress(string(limits[9]))) + ' '
         strings[8] = strings[8] + strtrim(strcompress(string(limits[5]))) + ' '
         strings[9] = strings[9] + strtrim(strcompress(string(limits[6]))) + ' '
         strings[10] = strings[10] + strtrim(strcompress(string(limits[10]))) + ' '
      
   strings[0] = strmid(strings[0], 0, strlen(strings[0]) - 1)
   
   
   
   for i in arange(0, (array(inputarray, copy=0).nelements() - 1)+(1)):
      tmp = str_sep(strtrim(strcompress(inputarray[i]), 2), '=')
      _expr = tmp[0]
      if _expr == VARY:   
         inputarray[i] = strings[0]
      elif _expr == PARMAX:   
         inputarray[i] = strings[1]
      elif _expr == PARMIN:   
         inputarray[i] = strings[2]
      elif _expr == MODERATE:   
         inputarray[i] = strings[3]
      elif _expr == DELSTART:   
         inputarray[i] = strings[4]
      elif _expr == DELEND:   
         inputarray[i] = strings[5]
      elif _expr == ITESTART:   
         inputarray[i] = strings[6]
      elif _expr == ITEEND:   
         inputarray[i] = strings[7]
      elif _expr == SATDELT:   
         inputarray[i] = strings[8]
      elif _expr == MINDELTA:   
         inputarray[i] = strings[9]
      elif _expr == VARINDX:   
         inputarray[i] = strings[10]
         
      else:   
         print 'empty'
         
      
   
   
   return _ret()


