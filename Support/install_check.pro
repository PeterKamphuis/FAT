FUNCTION install_check,gdlidl

  
;+
; NAME:
;       INSTALL_CHECK
;
; PURPOSE:
;       Function to check the fit of a user of N2903 against the ones
;       run on IDL 7.1 and GDL 0.9.6
;
; CATEGORY:
;       Support
; 
; CALLING SEQUENCE:	
;       result = install_check(gdlidl)
;
; INPUTS:
;         gdlidl = identifier of running IDL or GDL. is one when
;         running GDL
;
; OPTIONAL INPUTS:
;  
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;       -
; 
; PROCEDURES CALLED:
;       WRITENEWTOTEMPLATE
;
; MODIFICATION HISTORY:
;       Written 17-05-2017 P.Kamphuis v1.0
;
; NOTE:
;     
;-

  COMPILE_OPT IDL2
                                ;IF anything goes wrong in this routine then the installation failed
  CATCH,Error_status  
  IF  Error_status NE 0. THEN BEGIN
     print, ' '
     print, 'Oops the following went wrong:'
     print,'-------------------------------------------'
     help,/Last_Message, Output = theerrmess
                                ;This gives trace back information in
                                ;IDL but unfortunately not in GDL
     for j=0,n_elements(theerrmess)-1 do print,theerrmess[j]
     print,'-------------------------------------------'
     CATCH,/cancel
     return,1
  endif


  
                                ;Let's first read in the newly fitted parameters
  checkpara=['RADI','SBR','SBR_2','VROT','VROT_2','PA','PA_2','INCL','INCL_2','BMAJ','SDIS','XPOS','XPOS_2','YPOS','VSYS']
  limits=[0.05,5e-4,5e-4,1,1,1,1,1,1,1e-6,1,1e-4,1e-4,1e-4,0.1]
  TemplateFit=1.
  WriteNewToTemplate,TemplateFit,'Installation_Check/Finalmodel/Finalmodel.def',ARRAYS=ArraysFit,VARIABLECHANGE=checkpara,/EXTRACT

                                ;First we check that the rotation curve has indeed been held symmetrical
  tmppos=WHERE(checkpara EQ 'VROT')
  vrot1=ArraysFit[*,tmppos[0]]
  tmppos=WHERE(checkpara EQ 'VROT_2')
  vrot2=ArraysFit[*,tmppos[0]]
  diffrot=TOTAL(vrot1-vrot2)
  IF diffrot NE 0 then return,2
                                ;And the same for the RA position
  tmppos=WHERE(checkpara EQ 'XPOS')
  xpos1=ArraysFit[*,tmppos[0]]
  tmppos=WHERE(checkpara EQ 'XPOS_2')
  xpos2=ArraysFit[*,tmppos[0]]
  diffrot=TOTAL(xpos1-xpos2)
  IF diffrot NE 0 then return,2
                                ;Now let's read in the provided
                                ;def files
  if gdlidl then begin
     WriteNewToTemplate,TemplateFit,'Installation_Check/Finalmodel_GDL.def',ARRAYS=ArraysProvided,VARIABLECHANGE=checkpara,/EXTRACT
  endif else begin
     WriteNewToTemplate,TemplateFit,'Installation_Check/Finalmodel_IDL.def',ARRAYS=ArraysProvided,VARIABLECHANGE=checkpara,/EXTRACT
  ENDELSE
  IF n_elements(ArraysFIT[*,0]) NE n_elements(ArraysProvided[*,0]) then return,2 
  for i=0,n_elements(ArraysFit[0,*])-1 do begin
     diff=TOTAL(ArraysFit[*,i]-ArraysProvided[*,i])/n_elements(ArraysFit[*,0])
     IF diff GT limits[i] then begin
        print,'    __________________INSTALLATION_CHECK___________________ '
        print,'    ------------------------------------------------------- '
        print,'    '+checkpara[i]+ ' differs too much in the fit from the provided input'
        print,'    The average difference is '+strtrim(string(diff),2)+'.'
        print,'    This is outside the provided limits.'
        print,'    -------------------------------------------------------'
        print,' '
        return,2
     ENDIF
  endfor  
  return,0
end

  
