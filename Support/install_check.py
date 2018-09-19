from numarray import *

def install_check(gdlidl):
   """
    NAME:
          INSTALL_CHECK
   
    PURPOSE:
          Function to check the fit of a user of N2903 against the ones
          run on IDL 7.1 and GDL 0.9.6
   
    CATEGORY:
          Support
   
    CALLING SEQUENCE:
          result = install_check(gdlidl)
   
    INPUTS:
            gdlidl = identifier of running IDL or GDL. is one when
            running GDL
   
    OPTIONAL INPUTS:
   
    KEYWORD PARAMETERS:
   
    OUTPUTS:
   
    OPTIONAL OUTPUTS:
          -
   
    PROCEDURES CALLED:
          WRITENEWTOTEMPLATE
   
    MODIFICATION HISTORY:
          Written 17-05-2017 P.Kamphuis v1.0
   
    NOTE:
   
   """

   n_params = 1
   
   # COMPILE_OPT IDL2
   #IF anything goes wrong in this routine then the installation failed
   catch(error_status)
   if error_status != 0.:   
      print ' '
      print 'Oops the following went wrong:'
      print '-------------------------------------------'
      help(last_message=True, output=theerrmess)
      #This gives trace back information in
      #IDL but unfortunately not in GDL
      for j in arange(0, (array(theerrmess, copy=0).nelements() - 1)+(1)):
         print theerrmess[j]
      print '-------------------------------------------'
      catch(cancel=True)
      return 1
   
   
   
   #Let's first read in the newly fitted parameters
   checkpara = concatenate(['RADI', 'SBR', 'SBR_2', 'VROT', 'VROT_2', 'PA', 'PA_2', 'INCL', 'INCL_2', 'BMAJ', 'SDIS', 'XPOS', 'XPOS_2', 'YPOS', 'VSYS'])
   limits = concatenate([0.05, 5e-4, 5e-4, 1, 1, 1, 1, 1, 1, 1e-6, 1, 1e-4, 1e-4, 1e-4, 0.1])
   templatefit = 1.
   writenewtotemplate(templatefit, 'Installation_Check/Def_Files/Finalmodel.def', arrays=arraysfit, variablechange=checkpara, extract=True)
   
   #First we check that the rotation curve has indeed been held symmetrical
   tmppos = where(ravel(checkpara == 'VROT'))[0]
   vrot1 = arraysfit[tmppos[0],:]
   tmppos = where(ravel(checkpara == 'VROT_2'))[0]
   vrot2 = arraysfit[tmppos[0],:]
   diffrot = total(vrot1 - vrot2)
   if diffrot != 0:   
      return 2
   #And the same for the RA position
   tmppos = where(ravel(checkpara == 'XPOS'))[0]
   xpos1 = arraysfit[tmppos[0],:]
   tmppos = where(ravel(checkpara == 'XPOS_2'))[0]
   xpos2 = arraysfit[tmppos[0],:]
   diffrot = total(xpos1 - xpos2)
   if diffrot != 0:   
      return 2
   #Now let's read in the provided
   #def files
   if gdlidl:   
      writenewtotemplate(templatefit, 'Installation_Check/Finalmodel_GDL.def', arrays=arraysprovided, variablechange=checkpara, extract=True)
   else:   
      writenewtotemplate(templatefit, 'Installation_Check/Finalmodel_IDL.def', arrays=arraysprovided, variablechange=checkpara, extract=True)
   if array(arraysfit[0,:], copy=0).nelements() != array(arraysprovided[0,:], copy=0).nelements():   
      return 2
   for i in arange(0, (array(arraysfit[:,0], copy=0).nelements() - 1)+(1)):
      diff = total(arraysfit[i,:] - arraysprovided[i,:]) / array(arraysfit[0,:], copy=0).nelements()
      if diff > limits[i]:   
         print '    __________________INSTALLATION_CHECK___________________ '
         print '    ------------------------------------------------------- '
         print '    ' + checkpara[i] + ' differs too much in the fit from the provided input'
         print '    The average difference is ' + strtrim(string(diff), 2) + '.'
         print '    This is outside the provided limits.'
         print '    -------------------------------------------------------'
         print ' '
         return 2
   return 0



