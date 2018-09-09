import numpy as np
def writenewtotemplate(Template=None, Newfilename= None, Variables=None, Arrays=[], Variablechange=['BMIN', 'BMAJ', 'BPA', 'RMS', 'DISTANCE', 'NUR', 'RADI', 'VROT', 'Z0', 'SBR', 'INCL', 'PA', 'XPOS', 'YPOS', 'VSYS', 'SDIS', 'VROT_2', 'Z0_2', 'SBR_2', 'INCL_2', 'PA_2', 'XPOS_2', 'YPOS_2', 'VSYS_2', 'SDIS_2', 'CONDISP', 'CFLUX', 'CFLUX_2'], Extract=False):
   """
    NAME:
          WRITENEWTOTEMPLATE
   
    PURPOSE:
          Routine to write a tirific def file into the template used
   
    CATEGORY:
          Support
   
    CALLING SEQUENCE:
          WRITENEWTOTEMPLATE,Template,NewFileName,VARIABLES=TemplateVariables,ARRAYS=Arrays,VARIABLECHANGE=VariableChange,/EXTRACT
   
   
    INPUTS:
          Template = Array with the Template That the New file needs to be written to
          NewFileName = Is the File Name Of the Tirific file containing the new parameters
   
    OPTIONAL INPUTS:
          Variables = gives an array indicating the variables in the
          Template it is not necessary to provide but saves time
          Arrays = is a 2D array which will contain the arrays
          [*,#changed variables] the ordering is the same as the order
          of the variables to be changed. However string values will be
          excluded even if they are switched. In that case they provide
          an empty array there to get the default ordering of
          variableChange and their names
          VARIABLECHANGE = a string array with the names of the variable
          that should be changed or extracted
   
    KEYWORD PARAMETERS:
          /EXTRACT - Set this keyword to not write the new values to the
                     template but merely extract them from the file
   
    OUTPUTS:
          inputarray = updated tirific template
   
    OPTIONAL OUTPUTS:
          -
   
    PROCEDURES CALLED:
          ISNUMERIC(), STR_SEP(), STRTRIM(), STRCOMPRESS()
   
    EXAMPLE:
   
   
    MODIFICATION HISTORY:
          01-08-2018 P.Kamphuis; started translation to python
          07-03-2017 P.Kamphuis; added the cflux parameters to the
                                 default variables
          25-02-2016 P.Kamphuis; Made an adjustment to always look for
                                 the presence of error parameters
          Written 01-01-2015 P.Kamphuis v1.0
   
    NOTE:
   
   """
   if Template == None and not Extract:
      print("You did not provide a Template and do not want to extract.\nThere is no point in this exercise\n")
      return
   #First we need to know where which variable is in the template when not given. If the template is not given and we're only extracting this is not necessary
   if Variables is None and Template != None:
      for line in Template:
         try:
            Variables.append(line.split('=')[0].strip())
         except AttributeError:
            Variables = [line.split('=')[0].strip()]
   # To be flexible we want to first create the Arrays to write values into
   # this array should be of the shape (NUR,Variables to change)
   f = open(Newfilename, 'r')
   for line in f:
      tmp = line.split('=')
      if tmp[0].strip() == 'NUR':
         Arrays=np.zeros([int(tmp[1]),len(Variablechange)])
         break
   f.close()
   #then open the previous fit 
   f = open(Newfilename, 'r')
   for line in f:
      #check which variable and if we want it changed
      tmp = line.split('=')
      try:
         pos_in_change = Variablechange.index(tmp[0].strip())
      except ValueError:
         pos_in_change = -1
      # we want it changed   
      if pos_in_change != -1:
         # if we have a template
         if Template != None and not Extract:
            try:
               pos_in_template = Variables.index(tmp[0].strip())
            except ValueError:
               pos_in_template = -1
            if pos_in_template != -1:
               Template[pos_in_template] = line
         try:
            Arrays[0:len(tmp[1].strip().split()),pos_in_change] =  tmp[1].strip().split()
         except ValueError:
            pass
      else:
         if tmp[0].split(' ')[0].strip() == '#':
            try:
               err_pos_in_change =  Variablechange.index(tmp[0].split(' ')[1].strip())
            except ValueError:
               err_pos_in_change = -1
            if err_pos_in_change != -1:
               try:
                  Arrays[0:len(tmp[1].strip().split()),err_pos_in_change] =  tmp[1].strip().split()
               except ValueError:
                  pass
   return Arrays
