from numarray import *

def create_residual_file(outname, data, model, header):
   n_params = 4
   def _ret():  return (outname, data, model, header)
   
   # COMPILE_OPT IDL2
   cmax = max(data - model, min=cmin)
   if finite(cmax):   
      sxaddpar(header, 'DATAMAX', cmax)
   if finite(cmin):   
      sxaddpar(header, 'DATAMIN', cmin)
   writefits(outname, data - model, header)
   
   return _ret()


def create_residuals(filenames, version):
   """
    NAME:
          CREATE_RESIDUALS,
   
    PURPOSE:
          Creation of residuals if requested.
   
    CATEGORY:
          Support
   
    CALLING SEQUENCE:
          CREATE_RESIDUALS,names,version
   
   
    INPUTS:
         filenames = the names of the files that are variable the order
         is [Cube,moment0,moment1,mask,noisemap,sofia_catalog,basicinfofilename]
         version = the requested amount of files. 0 just organize the
         output and keep all (This will also happen when a fit is
         unsuccesful); 1 remove optimized files, log files and input
         files; 2  remove optimized files, log files, input
         files, ps files and unsmoothed files; 3 (Default) remove optimized files, log files, input
         files, ps files, unsmoothed files and all model fits files
         except the final model; 4 keep only the def files and remove
         all other output. >5 is the same as 0. Residuals are created for
         all cases where the fits files are maintained. If 0.5 is added
         the final fit was the first fit.
   
    OPTIONAL INPUTS:
          LOG = name of the tracing log
   
    KEYWORD PARAMETERS:
          -
   
    OUTPUTS:
   
    OPTIONAL OUTPUTS:
          -
   
    PROCEDURES CALLED:
          CREATE_RESIDUALS,ORGANIZE_OUTPUT
   
    EXAMPLE:
   
   
    MODIFICATION HISTORY:
          08-06-2016 P.Kamphuis; Added a finite check for datamax and
                                 datamin in the residual to avoid crashes.
          Written 24-07-2015 P.Kamphuis v1.0
   
    NOTE:
   
   """

   n_params = 2
   def _ret():  return (filenames, version)
   
   # COMPILE_OPT IDL2
   exist = file_test('Residuals', directory=True)
   if exist == 0:   
      spawn('mkdir Residuals', isthere)
   cube = readfits(filenames[0] + '.fits', hed, silent=True)
   xv = readfits(filenames[0] + '_1_xv.fits', hed1xv, silent=True)
   mom0 = readfits(filenames[1] + '.fits', mom0hed, silent=True)
   mom1 = readfits(filenames[2] + '.fits', mom1hed, silent=True)
   if version == array(version, copy=0).astype(Int32):   
      xv2 = readfits(filenames[0] + '_2_xv.fits', hed2xv, silent=True)
      cubem1 = readfits('1stfit.fits', dummy, silent=True)
      cubem2 = readfits('2ndfit.fits', dummy, silent=True)
      xvm1 = readfits('1stfit_xv.fits', dummy, silent=True)
      xvm2 = readfits('2ndfit_xv.fits', dummy, silent=True)
      mom0m1 = readfits('1stfit_mom0.fits', dummy, silent=True)
      mom1m1 = readfits('1stfit_mom1.fits', dummy, silent=True)
      mom0m2 = readfits('2ndfit_mom0.fits', dummy, silent=True)
      mom1m2 = readfits('2ndfit_mom1.fits', dummy, silent=True)
      'Residuals/Cube_No_Warp.fits', cube, cubem1, hed = create_residual_file('Residuals/Cube_No_Warp.fits', cube, cubem1, hed)
      'Residuals/Cube_No_Warp_xv.fits', xv, xvm1, hed1xv = create_residual_file('Residuals/Cube_No_Warp_xv.fits', xv, xvm1, hed1xv)
      'Residuals/Cube_No_Warp_mom0.fits', mom0, mom0m1, mom0hed = create_residual_file('Residuals/Cube_No_Warp_mom0.fits', mom0, mom0m1, mom0hed)
      'Residuals/Cube_No_Warp_mom1.fits', mom1, mom1m1, mom1hed = create_residual_file('Residuals/Cube_No_Warp_mom1.fits', mom1, mom1m1, mom1hed)
      'Residuals/Cube_Finalmodel.fits', cube, cubem2, hed = create_residual_file('Residuals/Cube_Finalmodel.fits', cube, cubem2, hed)
      'Residuals/Cube_Finalmodel_xv.fits', xv2, xvm2, hed2xv = create_residual_file('Residuals/Cube_Finalmodel_xv.fits', xv2, xvm2, hed2xv)
      'Residuals/Cube_Finalmodel_mom0.fits', mom0, mom0m2, mom0hed = create_residual_file('Residuals/Cube_Finalmodel_mom0.fits', mom0, mom0m2, mom0hed)
      'Residuals/Cube_Finalmodel_mom1.fits', mom1, mom1m2, mom1hed = create_residual_file('Residuals/Cube_Finalmodel_mom1.fits', mom1, mom1m2, mom1hed)
   else:   
      cubem1 = readfits('1stfit.fits', dummy, silent=True)
      xvm1 = readfits('1stfit_xv.fits', dummy, silent=True)
      mom0m1 = readfits('1stfit_mom0.fits', dummy, silent=True)
      mom1m1 = readfits('1stfit_mom1.fits', dummy, silent=True)
      'Residuals/Cube_Finalmodel.fits', cube, cubem1, hed = create_residual_file('Residuals/Cube_Finalmodel.fits', cube, cubem1, hed)
      'Residuals/Cube_Finalmodel_xv.fits', xv, xvm1, hed1xv = create_residual_file('Residuals/Cube_Finalmodel_xv.fits', xv, xvm1, hed1xv)
      'Residuals/Cube_Finalmodel_mom0.fits', mom0, mom0m1, mom0hed = create_residual_file('Residuals/Cube_Finalmodel_mom0.fits', mom0, mom0m1, mom0hed)
      'Residuals/Cube_Finalmodel_mom1.fits', mom1, mom1m1, mom1hed = create_residual_file('Residuals/Cube_Finalmodel_mom1.fits', mom1, mom1m1, mom1hed)
   
   return _ret()



