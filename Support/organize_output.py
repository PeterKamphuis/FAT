from numarray import *

def organize_output(names, version, directories):
   """
    NAME:
          ORGANIZE_OUTPUT
   
    PURPOSE:
          ordering of FAT output
   
    CATEGORY:
          Support
   
    CALLING SEQUENCE:
          ORGANIZE_OUTPUT,names,version,directories
   
   
    INPUTS:
         names = the file names
         version = the requested amount of files. 0 just organize the
         output and keep all (This will also happen when a fit is
         unsuccesful); 1 remove optimized files, log files and input
         files; 2  remove optimized files, log files, input
         files, ps files and unsmoothed files; 3 (Default) remove optimized files, log files, input
         files, ps files, unsmoothed files and all model fits files
         except the final model; 4 keep only the def files and remove
         all other output. 5 indicates a failed fit clean up. >6 is the same as 0. Residuals are created for
         all cases where the fits files are maintained. If 0.5 is added
         the final fit was the first fit.
         directories = the directories to organize things into.
   
    OPTIONAL INPUTS:
          -
   
    KEYWORD PARAMETERS:
          -
   
    OUTPUTS:
   
    OPTIONAL OUTPUTS:
          -
   
    PROCEDURES CALLED:
   
   
    EXAMPLE:
   
   
    MODIFICATION HISTORY:
          Written 24-07-2015 P.Kamphuis v1.0
   
    NOTE:
   
   """

   n_params = 3
   def _ret():  return (names, version, directories)
   
   # COMPILE_OPT IDL2
   for index in arange(0, (array(directories, copy=0).nelements() - 1)+(1)):
   
      exist = file_test(directories[index], directory=True)
      if exist == 0:   
         spawn('mkdir ' + directories[index], isthere)
      _expr = directories[index]
      if _expr == Optimized:   
         spawn('ls -l *_opt*', list)
         if array(list, copy=0).nelements() > 1:   
            if array(version, copy=0).astype(Int32) == version:   
               if file_test('1stfit_opt.log'):   
                  spawn('mv 1stfit_opt.log Optimized/No_Warp_opt.log', isthere)
               if file_test('1stfit_opt.def'):   
                  spawn('mv 1stfit_opt.def Optimized/No_Warp_opt.def', isthere)
               if file_test('1stfit_opt.fits'):   
                  spawn('mv 1stfit_opt.fits Optimized/No_Warp_opt.fits', isthere)
               if file_test('1stfit_opt.ps'):   
                  spawn('mv 1stfit_opt.ps Optimized/No_Warp_opt.ps', isthere)
               if file_test('2ndfit_opt.log'):   
                  spawn('mv 2ndfit_opt.log Optimized/Finalmodel_opt.log', isthere)
               if file_test('2ndfit_opt.def'):   
                  spawn('mv 2ndfit_opt.def Optimized/Finalmodel_opt.def', isthere)
               if file_test('2ndfit_opt.fits'):   
                  spawn('mv 2ndfit_opt.fits Optimized/Finalmodel_opt.fits', isthere)
               if file_test('2ndfit_opt.ps'):   
                  spawn('mv 2ndfit_opt.ps Optimized/Finalmodel_opt.ps', isthere)
            else:   
               if file_test('1stfit_opt.log'):   
                  spawn('mv 1stfit_opt.log Optimized/Finalmodel_opt.log', isthere)
               if file_test('1stfit_opt.def'):   
                  spawn('mv 1stfit_opt.def Optimized/Finalmodel_opt.def', isthere)
               if file_test('1stfit_opt.fits'):   
                  spawn('mv 1stfit_opt.fits Optimized/Finalmodel_opt.fits', isthere)
               if file_test('1stfit_opt.ps'):   
                  spawn('mv 1stfit_opt.ps Optimized/Finalmodel_opt.ps', isthere)
            if file_test('1stfit_opt.log'):   
               spawn('mv 1stfit_opt.ps Optimized/Finalmodel_opt.ps', isthere)
         else:   
            if exist == 0:   
               spawn('rm -Rf Optimized')
         if file_test(names[0] + '_opt.fits'):   
            spawn('mv ' + names[0] + '_opt.fits Optimized/', isthere)
         spawn('ls Optimized', filled)
         if filled[0] == '':   
            spawn('rm -Rf Optimized')
      elif _expr == Finalmodel:   
         if array(version, copy=0).astype(Int32) == version:   
            if file_test('2ndfit.log'):   
               spawn('mv 2ndfit.log Finalmodel/Finalmodel.log', isthere)
            if file_test('2ndfit.def'):   
               spawn('mv 2ndfit.def Finalmodel/Finalmodel.def', isthere)
            if file_test('2ndfit.fits'):   
               spawn('mv 2ndfit.fits Finalmodel/Finalmodel.fits', isthere)
            if file_test('2ndfit.ps'):   
               spawn('mv 2ndfit.ps Finalmodel/Finalmodel.ps', isthere)
         else:   
            if file_test('1stfit.log'):   
               spawn('mv 1stfit.log Finalmodel/Finalmodel.log', isthere)
            if file_test('1stfit.def'):   
               spawn('mv 1stfit.def Finalmodel/Finalmodel.def', isthere)
            if file_test('1stfit.fits'):   
               spawn('mv 1stfit.fits Finalmodel/Finalmodel.fits', isthere)
            if file_test('1stfit.ps'):   
               spawn('mv 1stfit.ps Finalmodel/Finalmodel.ps', isthere)
         spawn('ls Finalmodel', filled)
         if filled[0] == '':   
            spawn('rm -Rf Finalmodel')
      elif _expr == No_Warp:   
         if array(version, copy=0).astype(Int32) == version:   
            if file_test('1stfit.log'):   
               spawn('mv 1stfit.log No_Warp/No_Warp.log', isthere)
            if file_test('1stfit.def'):   
               spawn('mv 1stfit.def No_Warp/No_Warp.def', isthere)
            if file_test('1stfit.fits'):   
               spawn('mv 1stfit.fits No_Warp/No_Warp.fits', isthere)
            if file_test('1stfit.ps'):   
               spawn('mv 1stfit.ps No_Warp/No_Warp.ps', isthere)
         else:   
            if exist == 0:   
               spawn('rm -Rf No_Warp')
         spawn('ls No_Warp', filled)
         if filled[0] == '':   
            spawn('rm -Rf No_Warp')
      elif _expr == Sofia_Output:   
         check = str_sep(names[3], '/')
         if bitwise_and(file_test(names[3] + '.fits'), check[0] != 'Sofia_Output'):   
            spawn('mv ' + names[3] + '.fits Sofia_Output/')
         check = str_sep(names[5], '/')
         if bitwise_and(file_test(names[5]), check[0] != 'Sofia_Output'):   
            spawn('mv ' + names[5] + ' Sofia_Output/')
         spawn('ls Sofia_Output', filled)
         if filled[0] == '':   
            spawn('rm -Rf Sofia_Output')
      elif _expr == Moments:   
         check = str_sep(names[1], '/')
         if bitwise_and(file_test(names[1] + '.fits'), check[0] != 'Moments'):   
            spawn('mv ' + names[1] + '.fits' + ' Moments/')
         check = str_sep(names[2], '/')
         if bitwise_and(file_test(names[2] + '.fits'), check[0] != 'Moments'):   
            spawn('mv ' + names[2] + '.fits' + ' Moments/')
         check = str_sep(names[4], '/')
         if bitwise_and(file_test(names[4] + '.fits'), check[0] != 'Moments'):   
            spawn('mv ' + names[4] + '.fits' + ' Moments/')
         if array(version, copy=0).astype(Int32) == version:   
            if file_test('1stfit_mom0.fits'):   
               spawn('mv 1stfit_mom0.fits Moments/No_Warp_mom0.fits')
            if file_test('2ndfit_mom0.fits'):   
               spawn('mv 2ndfit_mom0.fits Moments/Finalmodel_mom0.fits')
            if file_test('1stfit_mom1.fits'):   
               spawn('mv 1stfit_mom1.fits Moments/No_Warp_mom1.fits')
            if file_test('2ndfit_mom1.fits'):   
               spawn('mv 2ndfit_mom1.fits Moments/Finalmodel_mom1.fits')
         else:   
            if file_test('1stfit_mom0.fits'):   
               spawn('mv 1stfit_mom0.fits Moments/Finalmodel_mom0.fits')
            if file_test('1stfit_mom1.fits'):   
               spawn('mv 1stfit_mom1.fits Moments/Finalmodel_mom1.fits')
         spawn('ls Moments', filled)
         if filled[0] == '':   
            spawn('rm -Rf Moments')
      elif _expr == PV-Diagrams:   
         if array(version, copy=0).astype(Int32) == version:   
            if file_test('1stfit_xv.fits'):   
               spawn('mv 1stfit_xv.fits PV-Diagrams/No_Warp_xv.fits', isthere)
            if file_test('2ndfit_xv.fits'):   
               spawn('mv 2ndfit_xv.fits PV-Diagrams/Finalmodel_xv.fits', isthere)
         else:   
            if file_test('1stfit_xv.fits'):   
               spawn('mv 1stfit_xv.fits PV-Diagrams/Finalmodel_xv.fits', isthere)
         for i in arange(0, 3):
            if file_test(names[0] + '_' + strtrim(string(i, format='(i1)'), 2) + '_xv.fits'):   
               spawn('mv ' + names[0] + '_' + strtrim(string(i, format='(i1)'), 2) + '_xv.fits PV-Diagrams/', isthere)
         spawn('ls PV-Diagrams', filled)
         if filled[0] == '':   
            spawn('rm -Rf PV-Diagrams')
      elif _expr == Intermediate:   
         ext = concatenate(['log', 'def', 'ps', 'fits'])
         infiles1 = concatenate(['_opt', 'old', 'all'])
         infiles2 = concatenate(['_opt', 'old', 'unsmooth', 'uncor', 'slop'])
         outfiles1 = concatenate(['_opt', 'prev', 'first_correct_center'])
         outfiles2 = concatenate(['opt', 'prev', 'unsmoothed', 'uncorrected', 'sloped'])
         if array(version, copy=0).astype(Int32) == version:   
            for i in arange(0, (array(ext, copy=0).nelements() - 1)+(1)):
               for j in arange(0, (array(infiles1, copy=0).nelements() - 1)+(1)):
                  if file_test('1stfit' + infiles1[j] + '.' + ext[i]):   
                     spawn('mv 1stfit' + infiles1[j] + '.' + ext[i] + ' Intermediate/No_Warp_' + outfiles1[j] + '.' + ext[i], isthere)
            for i in arange(0, (array(ext, copy=0).nelements() - 1)+(1)):
               for j in arange(0, (array(infiles2, copy=0).nelements() - 1)+(1)):
                  if file_test('2ndfit' + infiles2[j] + '.' + ext[i]):   
                     spawn('mv 2ndfit' + infiles2[j] + '.' + ext[i] + ' Intermediate/Finalmodel_' + outfiles2[j] + '.' + ext[i], isthere)
         else:   
            for i in arange(0, (array(ext, copy=0).nelements() - 1)+(1)):
               for j in arange(0, (array(infiles1, copy=0).nelements() - 1)+(1)):
                  if file_test('1stfit' + infiles1[j] + '.' + ext[i]):   
                     spawn('mv 1stfit' + infiles1[j] + '.' + ext[i] + ' Intermediate/Finalmodel_' + outfiles1[j] + '.' + ext[i], isthere)
         
         if file_test('progress1.txt'):   
            spawn('mv progress1.txt Intermediate/', isthere)
         if file_test('progress2.txt'):   
            spawn('mv progress2.txt Intermediate/', isthere)
         if file_test('sofia_input.txt'):   
            spawn('mv sofia_input.txt Intermediate/', isthere)
         if file_test('tirific.def'):   
            spawn('mv tirific.def Intermediate/', isthere)
         spawn('ls Intermediate', filled)
         if filled[0] == '':   
            spawn('rm -Rf Intermediate')
      elif _expr == Def_Files:   
         #     spawn,'mv *.def Def_Files/',isthere
         if array(version, copy=0).astype(Int32) == version:   
            if file_test('1stfit.def'):   
               spawn('mv 1stfit.def Def_Files/No_Warp.def', isthere)
            if file_test('1stfit_opt.def'):   
               spawn('mv 1stfit_opt.def Def_Files/No_Warp_opt.def', isthere)
            if file_test('1stfitold.def'):   
               spawn('mv 1stfitold.def Def_Files/No_Warp_prev.def', isthere)
            if file_test('1stfitall.def'):   
               spawn('mv 1stfitall.def Def_Files/No_Warp_first_correct_center.def', isthere)
            if file_test('2ndfit.def'):   
               spawn('mv 2ndfit.def Def_Files/Finalmodel.def', isthere)
            if file_test('2ndfitold.def'):   
               spawn('mv 2ndfitold.def Def_Files/Finalmodel_prev.def', isthere)
            if file_test('2ndfitunsmooth.def'):   
               spawn('mv 2ndfitunsmooth.def Def_Files/Finalmodel_unsmoothed.def', isthere)
            if file_test('2ndfituncor.def'):   
               spawn('mv 2ndfituncor.def Def_Files/Finalmodel_uncorrected.def', isthere)
            if file_test('2ndfitslop.def'):   
               spawn('mv 2ndfitslop.def Def_Files/Finalmodel_sloped.def', isthere)
         else:   
            if file_test('1stfitall.def'):   
               spawn('mv 1stfitall.def Def_Files/Finalmodel_first_correct_center.def', isthere)
            if file_test('1stfit.def'):   
               spawn('mv 1stfit.def Def_Files/Finalmodel.def', isthere)
            if file_test('1stfit_opt.def'):   
               spawn('mv 1stfit_opt.def Def_Files/Finalmodel_opt.def', isthere)
            if file_test('1stfitold.def'):   
               spawn('mv 1stfitold.def Def_Files/Finalmodel_prev.def', isthere)
         if file_test('tirific.def'):   
            spawn('mv tirific.def Def_Files/', isthere)
         spawn('ls Def_Files', filled)
         if filled[0] == '':   
            spawn('rm -Rf Def_Files')
      else:   
         print 'ORGANIZE_OUTPUT: That is not a proper directory for output'
         if exist == 0:   
            spawn('rm -Rf ' + directories[index])
      
   
   return _ret()

