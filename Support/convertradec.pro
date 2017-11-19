Pro convertradec,RA,DEC,INVERT=invert,COLON=colon

;+
; NAME:
;       CONVERTRADEC
;
; PURPOSE:
;       Program to convert DEG to RA and DEC or vice versa
;
; CATEGORY:
;       Support
; 
; CALLING SEQUENCE:
;       CONVERTRADEC,RA,DEC,/INVERT,/COLON
;
;
; INPUTS:
;       RA = Right Ascension in Degrees
;       DEC = Declination in degrees
;
; OPTIONAL INPUTS:
;       - 
;
; KEYWORD PARAMETERS:
;       /INVERT - if set it is assumed that RA and DEC are in hour angle
;       /COLON - if set it the output will be hh:mm:ss instead of
;                correctly formatted    
;
; OUTPUTS:
;       RA = the transformed units
;       DEC = the transformed units
;
; OPTIONAL OUTPUTS:
; 
; PROCEDURES CALLED:
;       STRTRIM(),STRCOMPRESS(),STRING(),STR_SEP(),FLOOR()
;
; EXAMPLE:
;      
;
; MODIFICATION HISTORY:
;       15-12-2016 P. Kamphuis; Improved handling of -00d  
;       Written by P.Kamphuis 01-01-2015 
;
; NOTE:
;     
;-
  COMPILE_OPT IDL2 
  IF NOT keyword_set(invert) then begin
     xpos=RA
     ypos=DEC
     RA=strarr(n_elements(xpos))
     DEC=strarr(n_elements(xpos))
     for i=0,n_elements(RA)-1 do begin      
        xposh=floor((xpos[i]/360.)*24.)
        xposm=floor((((xpos[i]/360.)*24.)-xposh)*60.)
        xposs=(((((xpos[i]/360.)*24.)-xposh)*60.)-xposm)*60
        yposh=floor(ABS(ypos[i]*1.))
        yposm=floor((((ABS(ypos[i]*1.))-yposh)*60.))
        yposs=(((((ABS(ypos[i]*1.))-yposh)*60.)-yposm)*60)
        sign=string(round(ypos[i]/ABS(ypos[i])))
        IF keyword_set(colon) then begin
           RA[i]=strtrim(string(xposh,format='(i02)'),2)+':'+strtrim(string(xposm,format='(i02)'),2)+':'+strtrim(string(xposs,format='(f05.2)'),2)
           DEC[i]=strtrim(string(yposh,format='(i02)'),2)+':'+strtrim(string(yposm,format='(i02)'),2)+':'+strtrim(string(yposs,format='(f05.2)'),2)
           IF sign LT 0. then DEC[i]='-'+DEC[i]
        ENDIF ELSE BEGIN
           RA[i]=string(xposh,format='(i02)')+'h'+string(xposm,format='(i02)')+'m'+string(xposs,format='(f05.2)')+'s'
           DEC[i]=string(yposh,format='(i02)')+'d'+string(yposm,format='(i02)')+"'"+string(yposs,format='(f05.2)')+'"'
           IF sign LT 0. then DEC[i]='-'+DEC[i]
        ENDELSE
     endfor
  ENDIF ELSE BEGIN
     xpos=RA
     ypos=DEC
     RA=dblarr(n_elements(xpos))
     DEC=dblarr(n_elements(xpos))
     for i=0,n_elements(RA)-1 do begin
        tmp=str_sep(strtrim(strcompress(xpos[i]),2),':')
        RA[i]=(tmp[0]+((tmp[1]+(tmp[2]/60.))/60.))*15.
        tmp=str_sep(strtrim(strcompress(ypos[i]),2),':')
        DEC[i]=ABS(tmp[0])+((tmp[1]+(tmp[2]/60.))/60.)
        sign=STRMID(ypos[i],0,1)
        IF sign EQ '-' then DEC[i]=-DEC[i]
     endfor
  ENDELSE
end
