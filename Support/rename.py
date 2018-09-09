from numarray import *

def rename(stringin, stringreplace):
   """
    NAME:
          RENAME
   
    PURPOSE:
         Rename Tirific Output that is produced while the pipeline runs
   
    CATEGORY:
          Support
   
    CALLING SEQUENCE:
          RENAME,stringin,stringreplace,files
   
   
    INPUTS:
         stringin = the string to be replaced
         stringreplace = The string that will be put in
   
    OPTIONAL INPUTS:
   
    KEYWORD PARAMETERS:
          -
    OUTPUTS:
   
    OPTIONAL OUTPUTS:
          -
   
    PROCEDURES CALLED:
          SPAWN
   
    EXAMPLE:
   
   
    MODIFICATION HISTORY:
          Written 04-01-2016 P.Kamphuis v1.0
   
    NOTE: This only works for the set of files produced by tirific this
    is not a generic rename
   
   """

   n_params = 2
   def _ret():  return (stringin, stringreplace)
   
   extensions = concatenate(['def', 'log', 'ps', 'fits'])
   
   for i in arange(0, (array(extensions, copy=0).nelements() - 1)+(1)):
      if file_test(stringreplace + extensions[i]):   
         spawn('rm -f ' + stringreplace + extensions[i])
      if file_test(stringin + extensions[i]):   
         spawn('mv ' + stringin + extensions[i] + ' ' + stringreplace + extensions[i])
   
   
   return _ret()



