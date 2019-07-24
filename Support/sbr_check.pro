Pro sbr_check,template,templatevars,sbrarr,sbrarr2,cutoff,debug=debug

;+
; NAME:
;       SBR_CHECK
;
; PURPOSE:
;       Routine to check that the sbr arrays are reasonable for the fit
;
; CATEGORY:
;       Support
; 
; CALLING SEQUENCE:
;       SBR_CHECK,template,templatevars,sbrarr,sbrarr2,cutoff
;
;
; INPUTS:
;       template = tirific template
;       templatevars = template variable
;       SBRarr = approaching sbrarr
;       SBRarr2 = receding sbrarr
;       cutoff = sbr cutoff values
;
; OPTIONAL INPUTS:
;       - 
;
; KEYWORD PARAMETERS:
;       -
;
; OUTPUTS:
;       template = modified tirific template
;
; OPTIONAL OUTPUTS:
;       -
; 
; PROCEDURES CALLED:
;       -
;
; EXAMPLE:
;      
;
; MODIFICATION HISTORY:
;       15-03-2019 P. Kamphuis; Always write arrays back to template.            
;       Written 01-01-2015 P.Kamphuis v1.0
;
; NOTE:
;     
;-
  COMPILE_OPT IDL2 
  trigger=0.
  IF SBRarr[0] GT 2.*SBRarr[2] AND  SBRarr2[0] GT 2.*SBRarr2[2] AND SBRarr[0] GT 2.*SBRarr[1] AND  SBRarr2[0] GT 2.*SBRarr2[1] then begin
     IF (SBRarr[2]+SBRarr2[2])/2. GT cutoff[2] then begin
        SBRarr[0]=(SBRarr[2]+SBRarr2[2])/2. 
        SBRarr2[0]=(SBRarr[2]+SBRarr2[2])/2.
        SBRarr[1]=(SBRarr[1]+SBRarr2[1])/4.
        SBRarr2[1]=(SBRarr[1]+SBRarr2[1])/4.
     endif else begin
        SBRarr[0]=(SBRarr[2]+SBRarr2[2])/2.
        SBRarr[1]=1.5*cutoff[2]
        SBRarr[2]=1.5*cutoff[2]
        SBRarr2[0]=(SBRarr[2]+SBRarr2[2])/2.
        SBRarr2[1]=1.5*cutoff[2]
        SBRarr2[2]=1.5*cutoff[2]
     ENDELSE  
     trigger=1
  ENDIF
  tmp=WHERE(SBRarr LT 0)
  IF tmp[0] NE -1 then begin
     IF tmp[0] EQ 0. then SBRarr[0]=SBRarr[1] 
     for i=0,n_elements(tmp)-1 do begin
        IF tmp[i] NE 0 AND tmp[i] NE n_elements(SBRarr)-1 then SBRarr[tmp[i]]=(SBrarr[tmp[i]-1]+SBRarr[tmp[i]+1])/2.
        IF tmp[i] EQ 0 then SBRarr[tmp[i]]=SBRarr[tmp[i]+1]
        IF tmp[i] EQ n_elements(SBRarr)-1 then SBRarr[tmp[i]]=SBRarr[tmp[i]-1]
     endfor
     trigger=1
    
  ENDIF
  tmp=WHERE(SBRarr2 LT 0)
  IF tmp[0] NE -1 then begin
     IF tmp[0] EQ 0. then SBRarr2[0]=SBRarr2[1] 
     for i=0,n_elements(tmp)-1 do begin
        IF tmp[i] NE 0 AND tmp[i] NE n_elements(SBRarr2)-1 then SBRarr2[tmp[i]]=(SBrarr2[tmp[i]-1]+SBRarr2[tmp[i]+1])/2.
        IF tmp[i] EQ 0 then SBRarr2[tmp[i]]=SBRarr2[tmp[i]+1]
        IF tmp[i] EQ n_elements(SBRarr2)-1 then SBRarr2[tmp[i]]=SBRarr2[tmp[i]-1]
     endfor
     trigger=1
  
  ENDIF
  IF SBRarr[0] GT SBRarr[1] then begin
     SBRarr[0]= SBRarr[1]
     SBRarr2[0]= SBRarr[1]
     trigger=1    
  ENDIF
  
  SBRarr[0:1]=(SBRarr[0:1]+SBRarr2[0:1])/2.  
  SBRarr2[0:1]=SBRarr[0:1]
  if keyword_set(debug) then begin
     print,'!!!!!!!!!! How Does This Happen !!!!!!!!!!!!!!!!!!!!!'
     print,SBRarr,SBRarr2
     print,'!!!!!!!!!! How Does This Happen After !!!!!!!!!!!!!!!!!!!!!'
  endif
  tmppos=where('SBR' EQ templatevars)
  template[tmppos]='SBR= '+STRJOIN(strtrim(strcompress(string(SBRarr))),' ')
  tmppos=where('SBR_2' EQ templatevars)
  template[tmppos]='SBR_2= '+STRJOIN(strtrim(strcompress(string(SBRarr2))),' ')
end
