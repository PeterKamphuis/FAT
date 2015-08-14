Pro book_keeping,filenames,version,log=log
  
;+
; NAME:
;       BOOK_KEEPING
;
; PURPOSE:
;       Clean up, ordering and the creation of residuals if requested.
;
; CATEGORY:
;       Support
; 
; CALLING SEQUENCE:
;       BOOK_KEEPING,filenames,version,log=log
;
;
; INPUTS:
;      filenames = the names of the files that are variable the order
;      is [Cube,moment0,moment1,mask,noisemap,sofia_catalog,basicinfofilename]   
;      version = the requested amount of files. 0 just organize the
;      output and keep all (This will also happen when a fit is
;      unsuccesful); 1 remove optimized files, log files and input
;      files; 2  remove optimized files, log files, input
;      files, ps files and unsmoothed files; 3 (Default) remove optimized files, log files, input
;      files, ps files, unsmoothed files and all model fits files
;      except the final model; 4 keep only the def files and remove
;      all other output. 5 indicates a failed fit clean up. >6 is the same as 0. Residuals are created for
;      all cases where the fits files are maintained. If 0.5 is added
;      the final fit was the first fit.
;
; OPTIONAL INPUTS:
;       LOG = name of the tracing log 
;
; KEYWORD PARAMETERS:
;       -
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;       -
; 
; PROCEDURES CALLED:
;       CREATE_RESIDUALS,ORGANIZE_OUTPUT
;
; EXAMPLE:
;      
;
; MODIFICATION HISTORY:
;       Written 24-07-2015 P.Kamphuis v1.0
;
; NOTE:
;     
;-
  COMPILE_OPT IDL2
  spawn,'pwd',currentdir
  IF size(log,/TYPE) EQ 7 then begin
     openu,66,log,/APPEND
     printf,66,linenumber()+"Removing the following files from "+currentdir
     close,66
  ENDIF 
  case version of
     1 OR 1.5: begin
        spawn,'rm -f 1stfit_opt.log 1stfit.log 1stfitall.log 1stfitold.log 2ndfit.log 2ndfitold.log 2ndfituncor.log 2ndfitunsmooth.log 2ndfitslop.log 1stfit_opt.def 1stfit_opt.fits 1stfit_opt.ps 2ndfit_opt.def 2ndfit_opt.fits 2ndfit_opt.ps '+filenames[0]+'_opt.fits sofia_input.txt tirific.def',isthere 
        create_residuals,filenames,version
        organize_output,filenames,version, ['Intermediate','Finalmodel','No_Warp','Moments','PV-Diagrams','Sofia_Output']
        IF size(log,/TYPE) EQ 7 then begin
           openu,66,log,/APPEND
           printf,66,linenumber()+'rm -f 1stfit_opt.log 1stfit.log 1stfitall.log 1stfitold.log 2ndfit.log 2ndfitold.log 2ndfituncor.log 2ndfitunsmooth.log 2ndfitslop.log 1stfit_opt.def 1stfit_opt.fits 1stfit_opt.ps 2ndfit_opt.def 2ndfit_opt.fits 2ndfit_opt.ps '+filenames[0]+'_opt.fits sofia_input.txt tirific.def' 
           close,66
        ENDIF
       
     end    
     2 OR 2.5: begin
        spawn,'rm -f 1stfit_opt.log 1stfit.log 1stfitall.log 1stfitold.log 2ndfit.log 2ndfituncor.log 2ndfitold.log 2ndfitunsmooth.log 2ndfitslop.log 1stfit_opt.def 1stfit_opt.fits 1stfit_opt.ps 2ndfit_opt.def 2ndfit_opt.fits 2ndfit_opt.ps '+filenames[0]+'_opt.fits sofia_input.txt tirific.def 1stfit_opt.ps 1stfit.ps 1stfitold.ps 2ndfit.ps 2ndfitold.ps 2ndfitslop.ps 2ndfitunsmooth.ps  2ndfitunsmooth.fits  2ndfituncor.fits  progress1.txt progress2.txt',isthere 
        create_residuals,filenames,version
        organize_output,filenames,version, ['Intermediate','Finalmodel','No_Warp','Moments','PV-Diagrams','Sofia_Output']
        IF size(log,/TYPE) EQ 7 then begin
           openu,66,log,/APPEND
           printf,66,linenumber()+'rm -f 1stfit_opt.log 1stfit.log 1stfitall.log 1stfitold.log 2ndfitold.log 2ndfit.log 2ndfituncor.log 2ndfitunsmooth.log 2ndfitslop.log 1stfit_opt.def 1stfit_opt.fits 1stfit_opt.ps 2ndfit_opt.def 2ndfit_opt.fits 2ndfit_opt.ps '+filenames[0]+'_opt.fits sofia_input.txt tirific.def 1stfit_opt.ps 1stfit.ps 1stfitold.ps 2ndfit.ps 2ndfitold.ps 2ndfitslop.ps 2ndfitunsmooth.ps  2ndfitunsmooth.fits  2ndfituncor.fits  progress1.txt progress2.txt'
           close,66
        ENDIF
     end
     3 OR 3.5: begin
           spawn,'rm -f 1stfit_opt.log 1stfit.log 1stfitall.log 1stfitold.log 2ndfit.log 2ndfitold.log 2ndfituncor.log 2ndfitunsmooth.log 2ndfitslop.log 1stfit_opt.def 1stfit_opt.fits 1stfit_opt.ps 1stfitall.ps 1stfitall.fits 2ndfit_opt.def 2ndfit_opt.fits 2ndfit_opt.ps '+filenames[0]+'_opt.fits sofia_input.txt tirific.def 1stfit_opt.ps 1stfit.ps 1stfitold.ps 2ndfit.ps 2ndfitold.ps 2ndfitslop.ps 2ndfitunsmooth.ps 2ndfitslop.fits 2ndfitunsmooth.fits  2ndfituncor.fits 1stfitold.def 1stfitold.fits 2ndfitold.def 2ndfitold.fits  progress1.txt progress2.txt'+filenames[0]+'_0_xv.fits '+filenames[0]+'_1_xv.fits ',isthere 
        create_residuals,filenames,version
        organize_output,filenames,version, ['Intermediate','Finalmodel','No_Warp','Moments','PV-Diagrams','Sofia_Output']
        IF size(log,/TYPE) EQ 7 then begin
           openu,66,log,/APPEND
           printf,66,linenumber()+ spawn,'rm -f 1stfit_opt.log 1stfit.log 1stfitall.log 1stfitold.log 2ndfit.log 2ndfitold.log 2ndfituncor.log 2ndfitunsmooth.log 2ndfitslop.log 1stfit_opt.def 1stfit_opt.fits 1stfit_opt.ps 1stfitall.ps 1stfitall.fits 2ndfit_opt.def 2ndfit_opt.fits 2ndfit_opt.ps '+filenames[0]+'_opt.fits sofia_input.txt tirific.def 1stfit_opt.ps 1stfit.ps 1stfitold.ps 2ndfit.ps 2ndfitold.ps 2ndfitslop.ps 2ndfitunsmooth.ps 2ndfitslop.fits 2ndfitunsmooth.fits  2ndfituncor.fits 1stfitold.def 1stfitold.fits 2ndfitold.def 2ndfitold.fits  progress1.txt progress2.txt'+filenames[0]+'_0_xv.fits '+filenames[0]+'_1_xv.fits ' 
           close,66
        ENDIF
     end
     4 OR 4.5: begin        
             spawn,'rm -f 1stfit_opt.log 1stfit.log 1stfitold.log 1stfitall.log 2ndfitold.log 2ndfit.log 2ndfituncor.log 2ndfitunsmooth.log 2ndfitslop.log 1stfit_opt.def 1stfit_opt.fits 1stfit_opt.ps 1stfitall.ps 1stfitall.fits 2ndfit_opt.def 2ndfit_opt.fits 2ndfit_opt.ps '+filenames[0]+'_opt.fits sofia_input.txt tirific.def 1stfit_opt.ps 1stfit.ps 1stfitold.ps 2ndfit.ps 2ndfitold.ps 2ndfitslop.ps 2ndfitunsmooth.ps  2ndfitunsmooth.fits  2ndfituncor.fits 1stfitold.def 1stfitold.fits 2ndfitold.def 2ndfitold.fits'+filenames[1]+'.fits '+filenames[2]+'.fits '+filenames[0]+'_0_xv.fits '+filenames[0]+'_1_xv.fits '+filenames[3]+'.fits '+filenames[4]+'.fits '+filenames[5]+'.fits '+filenames[6]+'.fits 2ndfit.fits 1stfit.fits  progress1.txt progress2.txt',isthere 
        organize_output,filenames,version, ['Def_Files']
        IF size(log,/TYPE) EQ 7 then begin
           openu,66,log,/APPEND
           printf,66,linenumber()+'rm -f 1stfit_opt.log 1stfit.log 1stfitold.log 1stfitall.log 2ndfitold.log 2ndfit.log 2ndfituncor.log 2ndfitunsmooth.log 2ndfitslop.log 1stfit_opt.def 1stfit_opt.fits 1stfit_opt.ps 1stfitall.ps 1stfitall.fits 2ndfit_opt.def 2ndfit_opt.fits 2ndfit_opt.ps '+filenames[0]+'_opt.fits sofia_input.txt tirific.def 1stfit_opt.ps 1stfit.ps 1stfitold.ps 2ndfit.ps 2ndfitold.ps 2ndfitslop.ps 2ndfitunsmooth.ps  2ndfitunsmooth.fits  2ndfituncor.fits 1stfitold.def 1stfitold.fits 2ndfitold.def 2ndfitold.fits'+filenames[1]+'.fits '+filenames[2]+'.fits '+filenames[0]+'_0_xv.fits '+filenames[0]+'_1_xv.fits '+filenames[3]+'.fits '+filenames[4]+'.fits '+filenames[5]+'.fits '+filenames[6]+'.fits 2ndfit.fits 1stfit.fits  progress1.txt progress2.txt'
           close,66
        ENDIF
     end
     5:begin
        organize_output,filenames,version, ['Intermediate','Moments','Sofia_Output']
     end
     else:begin       
        create_residuals,filenames,version
        organize_output,filenames,version, ['Optimized','Intermediate','Finalmodel','No_Warp','Moments','PV-Diagrams','Sofia_Output']
        IF size(log,/TYPE) EQ 7 then begin
           openu,66,log,/APPEND
           printf,66,linenumber()+"none "+currentdir
           close,66
        ENDIF
     END


     
  endcase


end
