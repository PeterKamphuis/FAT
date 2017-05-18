Pro cleanup,name

;+
; NAME:
;       CLEANUP
;
; PURPOSE:
;       Clean up any existing files, it will only remove specific
;       files not directories or non-FAT files
;
; CATEGORY:
;       Support
; 
; CALLING SEQUENCE:
;       CLEANUP
;
;
; INPUTS:
;     
;
; OPTIONAL INPUTS:
;       
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
;       FILE_TEST()
;
; EXAMPLE:
;      
;
; MODIFICATION HISTORY:
;       28-04-2017 P. Kamphuis; Added the removal of the intermediate
;                               files. This will cause problems when
;                               using testing NE 0.   
;       Written 12-08-2015 P.Kamphuis v1.0
;
; NOTE:
;     
;-

dirs=['Optimized','Intermediate','Finalmodel','No_Warp','Moments','PV-Diagrams','Sofia_Output','Def_Files']
ext=['.fits','.log','.ps','.def']
for i=0,n_elements(dirs)-1 do begin
   IF FILE_TEST(dirs[i],/DIRECTORY) THEN begin
      CASE 1 of
        dirs[i] EQ 'Optimized' OR dirs[i] EQ 'Intermediate': begin
            CD,dirs[i]
            spawn,'rm -f Finalmodel* No_Warp*'
            CD,'../'
         end
         dirs[i] EQ 'Finalmodel' or dirs[i] EQ 'No_Warp':begin
            for j=0,n_elements(ext)-1 do begin
               spawn,'rm -f '+dirs[i]+'/'+dirs[i]+ext[j]
            endfor
         end
         dirs[i] EQ 'Moments':begin
            spawn,'rm -f Moments/No_Warp_mom0.fits'
            spawn,'rm -f Moments/Finalmodel_mom0.fits'
            spawn,'rm -f Moments/No_Warp_mom1.fits'
            spawn,'rm -f Moments/Finalmodel_mom1.fits'
            spawn,'rm -f Moments/'+name+'*_mom*.fits'
         end
         dirs[i] EQ 'PV-Diagrams':begin
            spawn,'rm -f PV-Diagrams/No_Warp_xv.fits',isthere
            spawn,'rm -f PV-Diagrams/Finalmodel_xv.fits',isthere
         end
         dirs[i] EQ 'Sofia_Output':begin
            spawn,'rm -f Sofia_Output/'+name+'*_binmask*',isthere
            spawn,'rm -f Sofia_Output/'+name+'*.ascii',isthere
         end
         dirs[i] EQ 'Def_Files':spawn,'rm -f Def_Files/*.def'
      endcase
   endif
endfor
spawn,'rm -f '+name+'*_small*',isthere
spawn,'rm -f 1stfit.* 1stfitold.* 1stfit_mom1.fits 1stfit_mom0.fits 1stfit_xv.fits 1stfit_opt.*',isthere
spawn,'rm -f 2ndfit.* 2ndfitold.* 2ndfit_mom1.fits 2ndfit_mom0.fits 2ndfit_xv.fits 2ndfit_opt.*',isthere
spawn,'rm -f tirific.def',isthere
spawn,'rm -f '+name+'*_cut*',isthere	
spawn,'rm -f '+name+'*_binmask*',isthere
spawn,'rm -f '+name+'*_mom*',isthere
spawn,'rm -f '+name+'BasicInfo*txt',isthere
end
