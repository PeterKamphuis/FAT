from numarray import *

def read_template(name, array, variables, sofia=None):
   """
    NAME:
          READ_TEMPLATE
   
    PURPOSE:
          Routine to read in a specific template file with = defined commands.
   
    CATEGORY:
          Support
   
    CALLING SEQUENCE:
          READ_TEMPLATE, name, array, variables, /SOFIA
   
   
    INPUTS:
          Name = name of the file to be read.
   
    OPTIONAL INPUTS:
          -
   
    KEYWORD PARAMETERS:
          /SOFIA - trigger to indicate it is a sofia input template not tirific
   
    OUTPUTS:
          array = An array with each line of the input file
          variable = an array with all the variables (i.e. the words
          appaearing before =) of the file
    OPTIONAL OUTPUTS:
          -
   
    PROCEDURES CALLED:
          FILE_LINES(), STR_SEP(), STRTRIM(), STRCOMPRESS()
   
    EXAMPLE:
   
   
    MODIFICATION HISTORY:
          Written 01-01-2015 P.Kamphuis v1.0
   
    NOTE:
   
   """

   n_params = 3
   _opt = (sofia,)
   def _ret():
      _optrv = zip(_opt, [sofia])
      _rv = [name, array, variables]
      _rv += [_o[1] for _o in _optrv if _o[0] is not None]
      return tuple(_rv)
   
   # COMPILE_OPT IDL2
   if bitwise_not((sofia is not None)):   
      openr(1, name)
      filelength = file_lines(name)
      h = ' '
      array = strarr(filelength)
      variables = strarr(filelength)
      for j in arange(0, (filelength - 1)+(1)):
         readf(1, h)
         tmp = str_sep(strtrim(strcompress(h), 2), '=')
         if array(tmp, copy=0).nelements() > 1:   
            variables[j] = tmp[0]
         array[j] = h
      close(1)
   else:   
      openr(1, name)
      filelength = file_lines(name)
      h = ' '
      array = strarr(filelength)
      variables = intarr(4)
      for j in arange(0, (filelength - 1)+(1)):
         readf(1, h)
         tmp = strtrim(strcompress(str_sep(h, '=')), 2)
         if array(tmp, copy=0).nelements() > 1:   
            _expr = tmp[0]
            if _expr == import.inFile:   
               variables[0] = j
            elif _expr == steps.doReliability:   
               variables[1] = j
            elif _expr == parameters.dilatePixMax:   
               variables[2] = j
            elif _expr == SCfind.threshold:   
               variables[3] = j
            else:   
               pass
            
         array[j] = h
      close(1)
   
   return _ret()

