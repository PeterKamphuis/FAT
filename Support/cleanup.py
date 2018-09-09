from numarray import *

def cleanup(name):
   """
    NAME:
          CLEANUP
   
    PURPOSE:
          Clean up any existing files, it will only remove specific
          files not directories or non-FAT files
   
    CATEGORY:
          Support
   
    CALLING SEQUENCE:
          CLEANUP
   
   
    INPUTS:
   
   
    OPTIONAL INPUTS:
   
   
    KEYWORD PARAMETERS:
          -
   
    OUTPUTS:
   
    OPTIONAL OUTPUTS:
          -
   
    PROCEDURES CALLED:
          FILE_TEST()
   
    EXAMPLE:
   
   
    MODIFICATION HISTORY:
          28-04-2017 P. Kamphuis; Added the removal of the intermediate
                                  files. This will cause problems when
                                  using testing NE 0.
          Written 12-08-2015 P.Kamphuis v1.0
   
    NOTE:
   
   """

   n_params = 1
   def _ret():  return name
   
   dirs = concatenate(['Optimized', 'Intermediate', 'Finalmodel', 'No_Warp', 'Moments', 'PV-Diagrams', 'Sofia_Output', 'Def_Files'])
   ext = concatenate(['.fits', '.log', '.ps', '.def'])
   for i in arange(0, (array(dirs, copy=0).nelements() - 1)+(1)):
      if file_test(dirs[i], directory=True):   
         _expr = 1
         if _expr == (bitwise_or(dirs[i] == 'Optimized', dirs[i] == 'Intermediate')):   
            cd(dirs[i])
            spawn('rm -f Finalmodel* No_Warp*')
            cd('../')
         elif _expr == (bitwise_or(dirs[i] == 'Finalmodel', dirs[i] == 'No_Warp')):   
            for j in arange(0, (array(ext, copy=0).nelements() - 1)+(1)):
               spawn('rm -f ' + dirs[i] + '/' + dirs[i] + ext[j])
         elif _expr == (dirs[i] == 'Moments'):   
            spawn('rm -f Moments/No_Warp_mom0.fits')
            spawn('rm -f Moments/Finalmodel_mom0.fits')
            spawn('rm -f Moments/No_Warp_mom1.fits')
            spawn('rm -f Moments/Finalmodel_mom1.fits')
            spawn('rm -f Moments/' + name + '*_mom*.fits')
         elif _expr == (dirs[i] == 'PV-Diagrams'):   
            spawn('rm -f PV-Diagrams/No_Warp_xv.fits', isthere)
            spawn('rm -f PV-Diagrams/Finalmodel_xv.fits', isthere)
         elif _expr == (dirs[i] == 'Sofia_Output'):   
            spawn('rm -f Sofia_Output/' + name + '*_binmask*', isthere)
            spawn('rm -f Sofia_Output/' + name + '*.ascii', isthere)
         elif _expr == (dirs[i] == 'Def_Files'):   
            spawn('rm -f Def_Files/*.def')
         else:
            raise RuntimeError('no match found for expression')
   spawn('rm -f ' + name + '*_small*', isthere)
   spawn('rm -f 1stfit.* 1stfitold.* 1stfit_mom1.fits 1stfit_mom0.fits 1stfit_xv.fits 1stfit_opt.*', isthere)
   spawn('rm -f 2ndfit.* 2ndfitold.* 2ndfit_mom1.fits 2ndfit_mom0.fits 2ndfit_xv.fits 2ndfit_opt.*', isthere)
   spawn('rm -f tirific.def', isthere)
   spawn('rm -f ' + name + '*_cut*', isthere)
   spawn('rm -f ' + name + '*_binmask*', isthere)
   spawn('rm -f ' + name + '*_mom*', isthere)
   spawn('rm -f ' + name + 'BasicInfo*txt', isthere)
   
   return _ret()

