from numarray import *

def changeradii(tirifictemplate, numberofrings):
   """
    NAME:
          CHANGERADII
   
    PURPOSE:
          Program to reduce the number of rings in a tirific template by 1 ring
   
    CATEGORY:
          Support
   
    CALLING SEQUENCE:
          CHANGERADII,tirifictemplate,numberofrings
   
   
    INPUTS:
          tirifictemplate =  a string array with all the input of a
          tirific def file.
          numberofrings = the new numberof ring requested
   
    OPTIONAL INPUTS:
          -
   
    KEYWORD PARAMETERS:
          -
   
    OUTPUTS:
          tirifictemplate =  the modified string array with all the
          input of a tirific def file
   
    OPTIONAL OUTPUTS:
          -
   
    PROCEDURES CALLED:
          STRTRIM(),STRCOMPRESS(),STRING(),STR_SEP(),WHERE(),DOUBLE(),FLOOR()
   
    EXAMPLE:
   
   
    MODIFICATION HISTORY:
        Written by P.Kamphuis 01-01-2015
   
    NOTE:
   
   """

   n_params = 2
   def _ret():  return (tirifictemplate, numberofrings)
   
   # COMPILE_OPT IDL2
   
   possiblevars = concatenate(['VROT', 'Z0', 'SBR', 'INCL', 'PA', 'XPOS', 'YPOS', 'VSYS', 'AZ1W', 'AZ1P', 'SDIS'])
   
   for i in arange(0, (array(tirifictemplate, copy=0).nelements() - 1)+(1)):
      tmp = str_sep(strtrim(strcompress(tirifictemplate[i]), 2), '=')
      if tmp[0] == 'NUR':   
         oldrings = tmp[1]
         tirifictemplate[i] = tmp[0] + '=' + strtrim(strcompress(string(numberofrings)), 1)
         break
   
   if oldrings > numberofrings:   
      oldrings = numberofrings
   
   for i in arange(0, (array(tirifictemplate, copy=0).nelements() - 1)+(1)):
   #first we split at the is sign
      tmp = str_sep(strtrim(strcompress(tirifictemplate[i]), 2), '=')
      if array(tmp, copy=0).nelements() > 1:   
         #then we check wether we have a disk variable
         checkdisk = str_sep(strtrim(strcompress(tmp[0]), 2), '_')
         if array(checkdisk, copy=0).nelements() > 1:   
            currentvar = checkdisk[0]
         else:   
            currentvar = tmp[0]
         
         found = where(ravel(currentvar == possiblevars))[0]
         if found[0] != -1:   
            tmp2 = str_sep(strtrim(strcompress(tmp[1]), 2), ' ')
            outputstring = strarr(1)
            outputstring = tmp[0] + '='
            if array(tmp2, copy=0).nelements() > 1:   
               if array(tmp2, copy=0).nelements() < numberofrings:   
                  amount = array(tmp2, copy=0).nelements() - 1
               else:   
                  amount = numberofrings - 1
               for j in arange(0, (amount)+(1)):
                  outputstring = outputstring + tmp2[j] + ' '
               tirifictemplate[i] = outputstring
            else:   
               tirifictemplate[i] = tmp[0] + '=' + tmp[1]
            
         if tmp[0] == 'RADI':   
            tmp2 = str_sep(strtrim(strcompress(tmp[1]), 2), ' ')
            outputstring = strarr(1)
            outputstring = tmp[0] + '='
            if array(tmp2, copy=0).nelements() < numberofrings:   
               for j in arange(0, (array(tmp2, copy=0).nelements() - 1)+(1)):
                  outputstring = outputstring + tmp2[j] + ' '
               for j in arange(array(tmp2, copy=0).nelements(), (numberofrings - 1)+(1)):
                  outputstring = outputstring + strtrim(strcompress(string(array(tmp2[array(tmp2, copy=0).nelements() - 1], copy=0).astype(Float64) + (j - array(array(tmp2, copy=0).nelements(), copy=0).astype(Float64) + 1) * (array(tmp2[array(tmp2, copy=0).nelements() - 1], copy=0).astype(Float64) - array(tmp2[array(tmp2, copy=0).nelements() - 2], copy=0).astype(Float64)))), 1) + ' '
            else:   
               for j in arange(0, (numberofrings - 1)+(1)):
                  outputstring = outputstring + tmp2[j] + ' '
            tirifictemplate[i] = outputstring
         if tmp[0] == 'VARY':   
            tirifictemplate[i] = tmp[0] + '= '
            tmpsep1 = str_sep(strtrim(strcompress(tmp[1]), 2), ',')
            for j in arange(0, (array(tmpsep1, copy=0).nelements() - 1)+(1)):
               tmpsep2 = str_sep(strtrim(strcompress(tmpsep1[j]), 2), ' ')
               for x in arange(0, (array(tmpsep2, copy=0).nelements() - 1)+(1)):
                  tmpsep3 = str_sep(strtrim(strcompress(tmpsep2[x]), 2), ':')
                  if array(tmpsep3, copy=0).nelements() == 1:   
                     if isnumeric(tmpsep2[x]):   
                        if array(tmpsep2[x], copy=0).astype(Float64) >= array(oldrings, copy=0).astype(Float64):   
                           tmpsep2[x] = strtrim(strcompress(string(numberofrings, format='(F7.4)')), 1)
                  else:   
                     for y in arange(0, (array(tmpsep3, copy=0).nelements() - 1)+(1)):
                        if isnumeric(tmpsep3[y]):   
                           if array(tmpsep3[y], copy=0).astype(Float64) >= array(oldrings, copy=0).astype(Float64):   
                              tmpsep3[y] = strtrim(strcompress(string(numberofrings, format='(F7.4)')), 1)
                        if y == 0:   
                           tmpsep2[x] = tmpsep3[y]
                        else:   
                           tmpsep2[x] = tmpsep2[x] + ':' + tmpsep3[y]
                  
                  if x == 0:   
                     tmpsep1[j] = tmpsep2[x]
                  else:   
                     tmpsep1[j] = tmpsep1[j] + ' ' + tmpsep2[x]
               if j == 0:   
                  tirifictemplate[i] = tirifictemplate[i] + tmpsep1[j]
               else:   
                  tirifictemplate[i] = tirifictemplate[i] + ',' + tmpsep1[j]
               
            
         #this part should be adjusted to accomodate any number of
         #regularization parameters
         if tmp[0] == 'REGNUME':   
            tirifictemplate[i] = 'REGNUME= ' + string(floor((numberofrings) / 2.)) + ',' + string(floor((numberofrings) / 2.))
         if tmp[0] == 'REGDENO':   
            #first we need to get the low number
            
            tmpnow = str_sep(strtrim(strcompress(tmp[1]), 2), ',')
            tmplow = str_sep(strtrim(strcompress(tmpnow[0]), 2), ':')
            if array(tmplow[0], copy=0).astype(Float64) == floor(numberofrings / 2.):   
               tmplow[0] = string(tmplow - 1)
            tirifictemplate[i] = 'REGDENO= ' + tmplow[0] + ':' + strtrim(strcompress(string(floor((numberofrings) / 2.))), 2) + ',' + tmplow[0] + ':' + strtrim(strcompress(string(floor((numberofrings) / 2.))), 2)
   
   
   
   return _ret()

