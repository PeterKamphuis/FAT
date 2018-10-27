Pro WriteNewToTemplate,Template,NewFileName,VARIABLES=TemplateVariables,ARRAYS=Arrays,VARIABLECHANGE=VariableChange,EXTRACT=extract
;+
; NAME:
;       WRITENEWTOTEMPLATE
;
; PURPOSE:
;       Routine to write a tirific def file into the template used
;       !!!!!!!!!!!!!!!!!Be aware it will only write the variables
;       indicated in VariableChange!!!!!!!!!!!!!!!!!!!!!!!!!!!
; CATEGORY:
;       Support
; 
; CALLING SEQUENCE:
;       WRITENEWTOTEMPLATE,Template,NewFileName,VARIABLES=TemplateVariables,ARRAYS=Arrays,VARIABLECHANGE=VariableChange,/EXTRACT
;
;
; INPUTS:
;       Template = Array with the Template That the New file needs to be written to
;       NewFileName = Is the File Name Of the Tirific file containing the new parameters
;
; OPTIONAL INPUTS:
;       Variables = gives an array indicating the variables in the
;       Template it is not necessary to provide but saves time 
;       Arrays = is a 2D array which will contain the arrays
;       [*,#changed variables] the ordering is the same as the order
;       of the variables to be changed. However string values will be
;       excluded even if they are switched. In that case they provide
;       an empty array there to get the default ordering of
;       variableChange and their names
;       VARIABLECHANGE = a string array with the names of the variable
;       that should be changed or extracted
;       
; KEYWORD PARAMETERS:
;       /EXTRACT - Set this keyword to not write the new values to the
;                  template but merely extract them from the file
;
; OUTPUTS:
;       inputarray = updated tirific template
;
; OPTIONAL OUTPUTS:
;       -
; 
; PROCEDURES CALLED:
;       ISNUMERIC(), STR_SEP(), STRTRIM(), STRCOMPRESS()
;
; EXAMPLE:
;      
;
; MODIFICATION HISTORY:
;       07-03-2017 P.Kamphuis; added the cflux parameters to the
;                              default variables 
;       25-02-2016 P.Kamphuis; Made an adjustment to always look for
;                              the presence of error parameters 
;       Written 01-01-2015 P.Kamphuis v1.0
;
; NOTE:
;      
;-
  COMPILE_OPT IDL2                                  
                                ;First we need to know where which variable is when not give
  IF NOT keyword_set(TemplateVariables) then begin
     TemplateVariables=strarr(n_elements(Template))
     for i=0,n_elements(Template)-1 do begin
        tmp=strtrim(str_sep(strtrim(strcompress(Template[i]),2),'='),2)     
        TemplateVariables[i]=tmp[0]
     ENDFOR
  ENDIF
  

;Lets make an default array with variables we want to transfer if it is not given
  IF NOT keyword_set(VariableChange) then begin
     VariableChange=['BMIN','BMAJ','BPA','RMS','DISTANCE','NUR','RADI','VROT',  'Z0',$
                     'SBR', 'INCL','PA','XPOS','YPOS','VSYS','SDIS','VROT_2',  'Z0_2','SBR_2', 'INCL_2','PA_2','XPOS_2','YPOS_2','VSYS_2','SDIS_2','CONDISP','CFLUX','CFLUX_2']
  ENDIF
  
                                ;then open the previous fit ;or when looping back open the previous fit

  close,1

  h=' '
  openr,1,NewFileName
  rings=0
  WHILE rings EQ 0 Do begin
     readf,1,h
     tmp=strtrim(str_sep(strtrim(strcompress(h),2),'='),2)
     IF tmp[0] EQ 'NUR' then begin
        Arrays=dblarr(fix(tmp[1]),n_elements(VariableChange))
        rings++
     ENDIF
  ENDWHILE
  close,1
  h=' '
  openr,1,NewFileName

  WHILE (NOT EOF(1)) DO BEGIN
     readf,1,h
     tmp=strtrim(str_sep(strtrim(strcompress(h),2),'='),2)
     varpos = WHERE(tmp[0] EQ VariableChange)
     IF varpos[0] NE -1 then begin
        tmppos=where(string(VariableChange[varpos[0]]) EQ string(TemplateVariables))
        IF NOT keyword_set(extract) then Template[tmppos]=h
        arr=str_sep(strtrim(strcompress(tmp[1]),2),' ')
        IF isnumeric(arr[0]) then begin
           Arrays[0:n_elements(arr)-1,varpos]=double(arr)
        ENDIF
     endif ELSE BEGIN
        tmp1=str_sep(strtrim(strcompress(tmp[0]),2),' ')
        IF tmp1[0] EQ '#' then begin
           tmp3=WHERE(tmp1[1] EQ VariableChange)
           IF tmp3[0] NE -1 then begin
              arr=str_sep(strtrim(strcompress(tmp[1]),2),' ')
              IF isnumeric(arr[0]) then begin
                 IF n_elements(arr) GT n_elements(Arrays[*,0]) then arr=[0,0]
                 Arrays[0:n_elements(arr)-1,tmp3]=double(arr)
              ENDIF
           ENDIF
        ENDIF
     ENDELSE  
  ENDWHILE
  close,1  
END
