from numarray import *

def book_keeping(filenames, version, distance, gdlidl, log=None, noise=None, finishafter=None, debug=None):
   """
    NAME:
          BOOK_KEEPING
   
    PURPOSE:
          Clean up, ordering and the creation of residuals if requested.
   
    CATEGORY:
          Support
   
    CALLING SEQUENCE:
          BOOK_KEEPING,filenames,version,distance,gdlidl,log=log,noise=noise,finishafter=finishafter,debug=debug
   
   
   
    INPUTS:
         filenames = the names of the files that are variable the order
                     is [Cube,moment0,moment1,mask,noisemap,sofia_catalog,basicinfofilename].
           version = the requested amount of files. 0 just organize the
                     output and keep all (This will also happen when a
                     fit is unsuccesful); 1 remove optimized files, log
                     files and input files; 2  remove optimized files,
                     log files, input files, ps files and unsmoothed
                     files; 3 (Default) remove optimized files, log
                     files, input files, ps files, unsmoothed files and
                     all model fits files except the final model; 4 keep
                     only the def files and remove all other output. 5
                     indicates a failed fit clean up. >6 is the same as
                     0. Residuals are created for all cases where the
                     fits files are maintained. If 0.5 is added the
                     final fit was the first fit.
          distance = The distance to the galaxy for overview_plot.pro.
            gdlidl = Identifier whether running idl or gdl.
   
    OPTIONAL INPUTS:
               LOG = name of the tracing log.
             noise = The noise used for the fitting.
       finishafter = key for what kind of fitting was done.
   
    KEYWORD PARAMETERS:
            /DEBUG = Option for making the routine verbose.
   
    OUTPUTS:
   
    OPTIONAL OUTPUTS:
          -
   
    PROCEDURES CALLED:
          CREATE_RESIDUALS,ORGANIZE_OUTPUT, OVERVIEW_PLOT
   
    EXAMPLE:
   
   
    MODIFICATION HISTORY:
          03-12-2017 P.Kamphuis; Added a debug option for printing
                                 output.
          17-11-2017 P.Kamphuis; Split out the 4 and 4.5 condition to
                                 ensure they always work correctly.
          27-02-2016 P.Kamphuis; Modified to accomodate the creation of
                                 an overview plot. This does mean
                                 creating the full output first and
                                 deleting subsequently.
          Written 24-07-2015 P.Kamphuis v1.0
   
    NOTE:
   
   """

   n_params = 4
   _opt = (log, noise, finishafter, debug)
   def _ret():
      _optrv = zip(_opt, [log, noise, finishafter, debug])
      _rv = [filenames, version, distance, gdlidl]
      _rv += [_o[1] for _o in _optrv if _o[0] is not None]
      return tuple(_rv)
   
   # COMPILE_OPT IDL2
   spawn('pwd', currentdir)
   if (debug is not None):   
      print "The input names are"
      print filenames
   if bitwise_and(version != 5, finishafter != 0):   
      create_residuals(filenames, version)
   organize_output(filenames, version, concatenate(['Optimized', 'Intermediate', 'Finalmodel', 'No_Warp', 'Moments', 'PV-Diagrams', 'Sofia_Output']))
   if bitwise_and(version != 5, finishafter != 0):   
      overview_plot(distance, gdlidl, noise=noise, finishafter=finishafter, filenames=filenames, version=version)
   if size(log, type=True) == 7:   
      openu(66, log, append=True)
      printf(66, linenumber() + "BOOK_KEEPING: Removing the following files from " + currentdir)
      close(66)
   if finishafter == 0:   
      spawn('rm -Rf Optimized Intermediate Finalmodel No_Warp Residuals', isthere)
      spawn('rm -Rf PV-Diagrams/Finalmodel_xv.fits PV-Diagrams/No_Warp_xv.fits PV-Diagrams/Cube_preprocessed_small_1_xv.fits  PV-Diagrams/Cube_preprocessed_small_2_xv.fits')
      if size(log, type=True) == 7:   
         openu(66, log, append=True)
         printf(66, linenumber() + 'BOOK_KEEPING:rm -Rf Optimized Intermediate Finalmodel No_Warp Residuals')
         printf(66, linenumber() + 'BOOK_KEEPING:rm -Rf PV-Diagrams/Finalmodel_xv.fits PV-Diagrams/No_Warp_xv.fits PV-Diagrams/Cube_preprocessed_small_1_xv.fits  PV-Diagrams/Cube_preprocessed_small_2_xv.fits')
         close(66)
   else:   
      _expr = version
      if _expr == 1:   
         spawn('rm -Rf Optimized')
         if size(log, type=True) == 7:   
            openu(66, log, append=True)
            printf(66, linenumber() + 'BOOK_KEEPING: rm -Rf Optimized')
            close(66)
         spawn('rm -f No_Warp/No_Warp.log Intermediate/No_Warp_first_correct_center.log Intermediate/No_Warp_prev.log Finalmodel/Finalmodel.log Intermediate/Finalmodel_sloped.log Intermediate/Finalmodel_prev.log Intermediate/Finalmodel_uncorrected.log Intermediate/Finalmodel_unsmoothed.log Intermediate/sofia_input.txt Intermediate/tirific.def', isthere)
         if size(log, type=True) == 7:   
            openu(66, log, append=True)
            printf(66, linenumber() + 'BOOK_KEEPING: rm -f No_Warp/No_Warp.log Intermediate/No_Warp_first_correct_center.log Intermediate/No_Warp_prev.log Finalmodel/Finalmodel.log Intermediate/Finalmodel_sloped.log Intermediate/Finalmodel_prev.log Intermediate/Finalmodel_uncorrected.log Intermediate/Finalmodel_unsmoothed.log Intermediate/sofia_input.txt Intermediate/tirific.def')
            close(66)
      elif _expr == 1.5:   
         spawn('rm -Rf Optimized')
         if size(log, type=True) == 7:   
            openu(66, log, append=True)
            printf(66, linenumber() + 'BOOK_KEEPING: rm -Rf Optimized')
            close(66)
         spawn('rm -f Finalmodel/Finalmodel.log Intermediate/Finalmodel_first_correct_center.log Intermediate/Finalmodel_prev.log  Intermediate/sofia_input.txt Intermediate/tirific.def', isthere)
         if size(log, type=True) == 7:   
            openu(66, log, append=True)
            printf(66, linenumber() + 'BOOK_KEEPING: rm -f Finalmodel/Finalmodel.log Intermediate/Finalmodel_first_correct_center.log Intermediate/Finalmodel_prev.log  Intermediate/sofia_input.txt Intermediate/tirific.def')
            close(66)
      elif _expr == 2:   
         spawn('rm -Rf Optimized')
         if size(log, type=True) == 7:   
            openu(66, log, append=True)
            printf(66, linenumber() + 'BOOK_KEEPING: rm -Rf Optimized')
            close(66)
         spawn('rm -f No_Warp/No_Warp.log Intermediate/No_Warp_first_correct_center.log Intermediate/No_Warp_prev.log Finalmodel/Finalmodel.log Intermediate/Finalmodel_sloped.log Intermediate/Finalmodel_prev.log Intermediate/Finalmodel_uncorrected.log Intermediate/Finalmodel_unsmoothed.log Intermediate/sofia_input.txt Intermediate/tirific.def No_Warp/No_Warp.ps Intermediate/No_Warp_prev.ps Finalmodel/Finalmodel.ps  Intermediate/Finalmodel_prev.ps Intermediate/Finalmodel_sloped.ps  Intermediate/Finalmodel_uncorrected.ps Intermediate/Finalmodel_unsmoothed.ps  Intermediate/Finalmodel_unsmoothed.fits Intermediate/Finalmodel_uncorrected.fits Intermediate/progress1.txt Intermediate/progress2.txt', isthere)
         if size(log, type=True) == 7:   
            openu(66, log, append=True)
            printf(66, linenumber() + 'BOOK_KEEPING: rm -f No_Warp/No_Warp.log Intermediate/No_Warp_first_correct_center.log Intermediate/No_Warp_prev.log Finalmodel/Finalmodel.log Intermediate/Finalmodel_sloped.log Intermediate/Finalmodel_prev.log Intermediate/Finalmodel_uncorrected.log Intermediate/Finalmodel_unsmoothed.log Intermediate/sofia_input.txt Intermediate/tirific.def No_Warp/No_Warp.ps Intermediate/No_Warp_prev.ps Finalmodel/Finalmodel.ps  Intermediate/Finalmodel_prev.ps Intermediate/Finalmodel_sloped.ps  Intermediate/Finalmodel_uncorrected.ps Intermediate/Finalmodel_unsmoothed.ps  Intermediate/Finalmodel_unsmoothed.fits Intermediate/Finalmodel_uncorrected.fits Intermediate/progress1.txt Intermediate/progress2.txt')
            close(66)
      elif _expr == 2.5:   
         spawn('rm -Rf Optimized')
         if size(log, type=True) == 7:   
            openu(66, log, append=True)
            printf(66, linenumber() + 'BOOK_KEEPING: rm -Rf Optimized')
            close(66)
         spawn('rm -f Finalmodel/Finalmodel.log Intermediate/Finalmodel_first_correct_center.log Intermediate/Finalmodel_prev.log  Intermediate/sofia_input.txt Intermediate/tirific.def Finalmodel/Finalmodel.ps Intermediate/Finalmodel_prev.ps Intermediate/progress1.txt Intermediate/progress2.txt', isthere)
         if size(log, type=True) == 7:   
            openu(66, log, append=True)
            printf(66, linenumber() + 'BOOK_KEEPING:rm -f Finalmodel/Finalmodel.log Intermediate/Finalmodel_first_correct_center.log Intermediate/Finalmodel_prev.log  Intermediate/sofia_input.txt Intermediate/tirific.def Finalmodel/Finalmodel.ps Intermediate/Finalmodel_prev.ps Intermediate/progress1.txt Intermediate/progress2.txt')
            close(66)
      elif _expr == 3:   
         spawn('rm -Rf Optimized')
         if size(log, type=True) == 7:   
            openu(66, log, append=True)
            printf(66, linenumber() + 'BOOK_KEEPING: rm -Rf Optimized')
            close(66)
         spawn('rm -f No_Warp/No_Warp.log Intermediate/No_Warp_first_correct_center.log Intermediate/No_Warp_prev.log Finalmodel/Finalmodel.log Intermediate/Finalmodel_prev.log   Intermediate/Finalmodel_uncorrected.log Intermediate/Finalmodel_unsmoothed.log  Intermediate/Finalmodel_sloped.log Intermediate/No_Warp_first_correct_center.ps Intermediate/No_Warp_first_correct_center.fits Intermediate/sofia_input.txt Intermediate/tirific.def No_Warp/No_Warp.ps Intermediate/No_Warp_prev.ps Finalmodel/Finalmodel.ps Intermediate/Finalmodel_prev.ps Intermediate/Finalmodel_sloped.ps  Intermediate/Finalmodel_uncorrected.ps Intermediate/Finalmodel_unsmoothed.ps Intermediate/Finalmodel_sloped.fits  Intermediate/Finalmodel_unsmoothed.fits Intermediate/Finalmodel_uncorrected.fits  Intermediate/No_Warp_prev.def  Intermediate/No_Warp_prev.fits Intermediate/Finalmodel_prev.def Intermediate/Finalmodel_prev.fits Intermediate/progress1.txt Intermediate/progress2.txt', isthere)
         if size(log, type=True) == 7:   
            openu(66, log, append=True)
            printf(66, linenumber() + 'BOOK_KEEPING:rm -f No_Warp/No_Warp.log Intermediate/No_Warp_first_correct_center.log Intermediate/No_Warp_prev.log Finalmodel/Finalmodel.log Intermediate/Finalmodel_prev.log   Intermediate/Finalmodel_uncorrected.log Intermediate/Finalmodel_unsmoothed.log  Intermediate/Finalmodel_sloped.log Intermediate/No_Warp_first_correct_center.ps Intermediate/No_Warp_first_correct_center.fits Intermediate/sofia_input.txt Intermediate/tirific.def No_Warp/No_Warp.ps Intermediate/No_Warp_prev.ps Finalmodel/Finalmodel.ps Intermediate/Finalmodel_prev.ps Intermediate/Finalmodel_sloped.ps  Intermediate/Finalmodel_uncorrected.ps Intermediate/Finalmodel_unsmoothed.ps Intermediate/Finalmodel_sloped.fits  Intermediate/Finalmodel_unsmoothed.fits Intermediate/Finalmodel_uncorrected.fits  Intermediate/No_Warp_prev.def  Intermediate/No_Warp_prev.fits Intermediate/Finalmodel_prev.def Intermediate/Finalmodel_prev.fits Intermediate/progress1.txt Intermediate/progress2.txt')
            close(66)
      elif _expr == 3.5:   
         spawn('rm -Rf Optimized')
         if size(log, type=True) == 7:   
            openu(66, log, append=True)
            printf(66, linenumber() + 'BOOK_KEEPING: rm -Rf Optimized')
            close(66)
         spawn('rm -f Finalmodel/Finalmodel.log Intermediate/Finalmodel_first_correct_center.log Intermediate/Finalmodel_prev.log Intermediate/Finalmodel_first_correct_center.ps Intermediate/Finalmodel_first_correct_center.fits Intermediate/sofia_input.txt Intermediate/tirific.def Finalmodel/Finalmodel.ps Intermediate/Finalmodel_prev.ps Intermediate/Finalmodel_prev.def  Intermediate/Finalmodel_prev.fits Intermediate/progress1.txt Intermediate/progress2.txt', isthere)
         if size(log, type=True) == 7:   
            openu(66, log, append=True)
            printf(66, linenumber() + 'BOOK_KEEPING:rm -f Finalmodel/Finalmodel.log Intermediate/Finalmodel_first_correct_center.log Intermediate/Finalmodel_prev.log Intermediate/Finalmodel_first_correct_center.ps Intermediate/Finalmodel_first_correct_center.fits Intermediate/sofia_input.txt Intermediate/tirific.def Finalmodel/Finalmodel.ps Intermediate/Finalmodel_prev.ps Intermediate/Finalmodel_prev.def  Intermediate/Finalmodel_prev.fits Intermediate/progress1.txt Intermediate/progress2.txt')
            close(66)
      elif _expr == 4.0:   
         spawn('mkdir Def_Files', isthere)
         spawn('mv No_Warp/*.def Def_Files', isthere)
         spawn('mv Finalmodel/*.def Def_Files', isthere)
         spawn('mv Intermediate/*.def Def_Files', isthere)
         spawn('mv Optimized/*.def Def_Files', isthere)
         spawn('rm -Rf Optimized Intermediate Finalmodel No_Warp Moments PV-Diagrams Sofia_Output Residuals', isthere)
         if size(log, type=True) == 7:   
            openu(66, log, append=True)
            printf(66, linenumber() + 'BOOK_KEEPING:rm -Rf Optimized Intermediate Finalmodel No_Warp Moments PV-Diagrams Sofia_Output Residuals')
            close(66)
      elif _expr == 4.5:   
         spawn('mkdir Def_Files', isthere)
         spawn('mv No_Warp/*.def Def_Files', isthere)
         spawn('mv Finalmodel/*.def Def_Files', isthere)
         spawn('mv Intermediate/*.def Def_Files', isthere)
         spawn('mv Optimized/*.def Def_Files', isthere)
         spawn('rm -Rf Optimized Intermediate Finalmodel No_Warp Moments PV-Diagrams Sofia_Output Residuals', isthere)
         if size(log, type=True) == 7:   
            openu(66, log, append=True)
            printf(66, linenumber() + 'BOOK_KEEPING:rm -Rf Optimized Intermediate Finalmodel No_Warp Moments PV-Diagrams Sofia_Output Residuals')
            close(66)
      else:   
         if size(log, type=True) == 7:   
            openu(66, log, append=True)
            printf(66, linenumber() + "BOOK_KEEPING: none " + currentdir)
            close(66)
         
         
         
         
      
   
   
   return _ret()

