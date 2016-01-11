Pro organize_output,names,version,directories

    
;+
; NAME:
;       ORGANIZE_OUTPUT
;
; PURPOSE:
;       ordering of FAT output
;
; CATEGORY:
;       Support
; 
; CALLING SEQUENCE:
;       ORGANIZE_OUTPUT,names,version,directories
;
;
; INPUTS:
;      names = the file names
;      version = the requested amount of files.
;      directories = the directories to organize things into.   
;
; OPTIONAL INPUTS:
;       -
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
;       
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
  for index=0,n_elements(directories)-1 do begin
     
     exist=FILE_TEST(directories[index],/DIRECTORY)
     IF exist EQ 0 then spawn,'mkdir '+directories[index],isthere
     case directories[index] of
        'Optimized':begin
           spawn,'ls -l *_opt*',list
           if n_elements(list) GT 1 then begin
              if fix(version) EQ version then begin
                 IF FILE_TEST('1stfit_opt.log') then spawn,'mv 1stfit_opt.log Optimized/No_Warp_opt.log',isthere
                 IF FILE_TEST('1stfit_opt.def') then spawn,'mv 1stfit_opt.def Optimized/No_Warp_opt.def',isthere
                 IF FILE_TEST('1stfit_opt.fits') then spawn,'mv 1stfit_opt.fits Optimized/No_Warp_opt.fits',isthere
                 IF FILE_TEST('1stfit_opt.ps') then spawn,'mv 1stfit_opt.ps Optimized/No_Warp_opt.ps',isthere
                 IF FILE_TEST('2ndfit_opt.log') then spawn,'mv 2ndfit_opt.log Optimized/Finalmodel_opt.log',isthere
                 IF FILE_TEST('2ndfit_opt.def') then spawn,'mv 2ndfit_opt.def Optimized/Finalmodel_opt.def',isthere
                 IF FILE_TEST('2ndfit_opt.fits') then spawn,'mv 2ndfit_opt.fits Optimized/Finalmodel_opt.fits',isthere
                 IF FILE_TEST('2ndfit_opt.ps') then spawn,'mv 2ndfit_opt.ps Optimized/Finalmodel_opt.ps',isthere
              ENDIF ELSE BEGIN
                 IF FILE_TEST('1stfit_opt.log') then spawn,'mv 1stfit_opt.log Optimized/Finalmodel_opt.log',isthere
                 IF FILE_TEST('1stfit_opt.def') then spawn,'mv 1stfit_opt.def Optimized/Finalmodel_opt.def',isthere
                 IF FILE_TEST('1stfit_opt.fits') then spawn,'mv 1stfit_opt.fits Optimized/Finalmodel_opt.fits',isthere
                 IF FILE_TEST('1stfit_opt.ps') then spawn,'mv 1stfit_opt.ps Optimized/Finalmodel_opt.ps',isthere
              ENDELSE
              IF FILE_TEST('1stfit_opt.log') then spawn,'mv 1stfit_opt.ps Optimized/Finalmodel_opt.ps',isthere
           endif else begin
              if exist EQ 0 then spawn,'rm -Rf Optimized'
           endelse
           IF FILE_TEST(names[0]+'_opt.fits') then spawn,'mv '+names[0]+'_opt.fits Optimized/',isthere
           spawn,'ls Optimized',filled
           IF filled[0] EQ '' then spawn,'rm -Rf Optimized'
        end
        'Finalmodel':begin
           if fix(version) EQ version then begin
              IF FILE_TEST('2ndfit.log') then spawn,'mv 2ndfit.log Finalmodel/Finalmodel.log',isthere
              IF FILE_TEST('2ndfit.def') then spawn,'mv 2ndfit.def Finalmodel/Finalmodel.def',isthere
              IF FILE_TEST('2ndfit.fits') then spawn,'mv 2ndfit.fits Finalmodel/Finalmodel.fits',isthere
              IF FILE_TEST('2ndfit.ps') then spawn,'mv 2ndfit.ps Finalmodel/Finalmodel.ps',isthere
           ENDIF ELSE BEGIN
              IF FILE_TEST('1stfit.log') then spawn,'mv 1stfit.log Finalmodel/Finalmodel.log',isthere
              IF FILE_TEST('1stfit.def') then spawn,'mv 1stfit.def Finalmodel/Finalmodel.def',isthere
              IF FILE_TEST('1stfit.fits') then spawn,'mv 1stfit.fits Finalmodel/Finalmodel.fits',isthere
              IF FILE_TEST('1stfit.ps') then spawn,'mv 1stfit.ps Finalmodel/Finalmodel.ps',isthere
           ENDELSE
           spawn,'ls Finalmodel',filled
           IF filled[0] EQ '' then spawn,'rm -Rf Finalmodel'
        end
        'No_Warp':begin
           if fix(version) EQ version then begin
              IF FILE_TEST('1stfit.log') then spawn,'mv 1stfit.log No_Warp/No_Warp.log',isthere
                 IF FILE_TEST('1stfit.def') then spawn,'mv 1stfit.def No_Warp/No_Warp.def',isthere
                 IF FILE_TEST('1stfit.fits') then spawn,'mv 1stfit.fits No_Warp/No_Warp.fits',isthere
                 IF FILE_TEST('1stfit.ps') then spawn,'mv 1stfit.ps No_Warp/No_Warp.ps',isthere
           endif else begin
              if exist EQ 0 then spawn,'rm -Rf No_Warp'
           endelse
           spawn,'ls No_Warp',filled
           IF filled[0] EQ '' then spawn,'rm -Rf No_Warp'
        end
        'Sofia_Output':begin
           check=str_sep(names[3],'/')
           IF FILE_TEST(names[3]+'.fits') AND check[0] NE 'Sofia_Output' then spawn,'mv '+names[3]+'.fits Sofia_Output/'
           check=str_sep(names[5],'/')
           IF FILE_TEST(names[5]) AND check[0] NE 'Sofia_Output' then spawn,'mv '+names[5]+' Sofia_Output/'
           spawn,'ls Sofia_Output',filled
           IF filled[0] EQ '' then spawn,'rm -Rf Sofia_Output'
        end
        'Moments':begin
           check=str_sep(names[1],'/')
           IF FILE_TEST(names[1]+'.fits') AND check[0] NE 'Moments' then spawn,'mv '+names[1]+'.fits'+' Moments/'
           check=str_sep(names[2],'/')
           IF FILE_TEST(names[2]+'.fits') AND check[0] NE 'Moments' then spawn,'mv '+names[2]+'.fits'+' Moments/'
           check=str_sep(names[4],'/')
           IF FILE_TEST(names[4]+'.fits') AND check[0] NE 'Moments' then spawn,'mv '+names[4]+'.fits'+' Moments/'
           if fix(version) EQ version then begin
              IF FILE_TEST('1stfit_mom0.fits') then spawn,'mv 1stfit_mom0.fits Moments/No_Warp_mom0.fits'
              IF FILE_TEST('2ndfit_mom0.fits') then spawn,'mv 2ndfit_mom0.fits Moments/Finalmodel_mom0.fits'
              IF FILE_TEST('1stfit_mom1.fits') then spawn,'mv 1stfit_mom1.fits Moments/No_Warp_mom1.fits'
              IF FILE_TEST('2ndfit_mom1.fits') then spawn,'mv 2ndfit_mom1.fits Moments/Finalmodel_mom1.fits'
           ENDIF ELSE begin
              IF FILE_TEST('1stfit_mom0.fits') then spawn,'mv 1stfit_mom0.fits Moments/Finalmodel_mom0.fits'
              IF FILE_TEST('1stfit_mom1.fits') then spawn,'mv 1stfit_mom1.fits Moments/Finalmodel_mom1.fits'
           ENDELSE
           spawn,'ls Moments',filled
           IF filled[0] EQ '' then spawn,'rm -Rf Moments'
        end
        'PV-Diagrams':begin
           if fix(version) EQ version then begin
              IF FILE_TEST('1stfit_xv.fits') then spawn,'mv 1stfit_xv.fits PV-Diagrams/No_Warp_xv.fits',isthere
              IF FILE_TEST('2ndfit_xv.fits') then spawn,'mv 2ndfit_xv.fits PV-Diagrams/Finalmodel_xv.fits',isthere
           ENDIF ELSE begin
              IF FILE_TEST('1stfit_xv.fits') then spawn,'mv 1stfit_xv.fits PV-Diagrams/Finalmodel_xv.fits',isthere
           ENDELSE
           for i=0,2 do begin
              IF FILE_TEST(names[0]+'_'+strtrim(string(i,format='(i1)'),2)+'_xv.fits') then spawn,'mv '+names[0]+'_'+strtrim(string(i,format='(i1)'),2)+'_xv.fits PV-Diagrams/',isthere
           endfor
           spawn,'ls PV-Diagrams',filled
           IF filled[0] EQ '' then spawn,'rm -Rf PV-Diagrams'
        end
        'Intermediate':begin
           ext=['log','def','ps','fits']
           infiles1=['_opt','old','all']
           infiles2=['_opt','old','unsmooth','uncor','slop']
           outfiles1=['_opt','prev','first_correct_center']
           outfiles2=['opt','prev','unsmoothed','uncorrected','sloped']
           IF fix(version) EQ version then begin
              for i=0,n_elements(ext)-1 do begin
                 for j=0,n_elements(infiles1)-1 do begin
                    IF FILE_TEST('1stfit'+infiles1[j]+'.'+ext[i]) then spawn,'mv 1stfit'+infiles1[j]+'.'+ext[i]+' Intermediate/No_Warp_'+outfiles1[j]+'.'+ext[i],isthere
                 endfor
              endfor
              for i=0,n_elements(ext)-1 do begin
                 for j=0,n_elements(infiles2)-1 do begin
                    IF FILE_TEST('2ndfit'+infiles2[j]+'.'+ext[i]) then spawn,'mv 2ndfit'+infiles2[j]+'.'+ext[i]+' Intermediate/Finalmodel_'+outfiles2[j]+'.'+ext[i],isthere
                 endfor
              endfor
           ENDIF ELSE BEGIN
              for i=0,n_elements(ext)-1 do begin
                 for j=0,n_elements(infiles1)-1 do begin
                    IF FILE_TEST('1stfit'+infiles1[j]+'.'+ext[i]) then spawn,'mv 1stfit'+infiles1[j]+'.'+ext[i]+' Intermediate/Finalmodel_'+outfiles1[j]+'.'+ext[i],isthere
                 endfor
              endfor
           ENDELSE
          
           IF FILE_TEST('progress1.txt') then spawn,'mv progress1.txt Intermediate/',isthere
           IF FILE_TEST('progress2.txt') then spawn,'mv progress2.txt Intermediate/',isthere
           IF FILE_TEST('sofia_input.txt') then spawn,'mv sofia_input.txt Intermediate/',isthere
           IF FILE_TEST('tirific.def') then spawn,'mv tirific.def Intermediate/',isthere
           spawn,'ls Intermediate',filled
           IF filled[0] EQ '' then spawn,'rm -Rf Intermediate'
        end
        'Def_Files':begin
           spawn,'mv *.def Def_Files/',isthere
           IF fix(version) EQ version then begin 
              IF FILE_TEST('1stfit.def') then spawn,'mv 1stfit.def Def_Files/No_Warp.def',isthere
              IF FILE_TEST('1stfit_opt.def') then spawn,'mv 1stfit_opt.def Def_Files/No_Warp_opt.def',isthere
              IF FILE_TEST('1stfitold.def') then spawn,'mv 1stfitold.def Def_Files/No_Warp_prev.def',isthere
              IF FILE_TEST('1stfitall.def') then spawn,'mv 1stfitall.def Def_Files/No_Warp_first_correct_center.def',isthere
              IF FILE_TEST('2ndfit.def') then spawn,'mv 2ndfit.def Def_Files/Finalmodel.def',isthere
              IF FILE_TEST('2ndfitold.def') then spawn,'mv 2ndfitold.def Def_Files/Finalmodel_prev.def',isthere
              IF FILE_TEST('2ndfitunsmooth.def') then spawn,'mv 2ndfitunsmooth.def Def_Files/Finalmodel_unsmoothed.def',isthere
              IF FILE_TEST('2ndfituncor.def') then spawn,'mv 2ndfituncor.def Def_Files/Finalmodel_uncorrected.def',isthere
              IF FILE_TEST('2ndfitslop.def') then spawn,'mv 2ndfitslop.def Def_Files/Finalmodel_sloped.def',isthere
           ENDIF ELSE BEGIN
              IF FILE_TEST('1stfitall.def') then spawn,'mv 1stfitall.def Def_Files/Finalmodel_first_correct_center.def',isthere
              IF FILE_TEST('1stfit.def') then spawn,'mv 1stfit.def Def_Files/Finalmodel.def',isthere
              IF FILE_TEST('1stfit_opt.def') then spawn,'mv 1stfit_opt.def Def_Files/Finalmodel_opt.def',isthere
              IF FILE_TEST('1stfitold.def') then spawn,'mv 1stfitold.def Def_Files/Finalmodel_prev.def',isthere
           ENDELSE
           IF FILE_TEST('tirific.def') then spawn,'mv tirific.def Def_Files/',isthere
           spawn,'ls Def_Files',filled
           IF filled[0] EQ '' then spawn,'rm -Rf Def_Files'
        end
        else:begin
           print,'ORGANIZE_OUTPUT: That is not a proper directory for output'
           if exist eq 0 then spawn,'rm -Rf '+directories[index]
        end
     endcase
  endfor
end
