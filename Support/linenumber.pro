Function linenumber

;+
; NAME:
;       LINENUMBER
;
; PURPOSE:
;       A function to print the linenumber in the routine directly below main
;
; CATEGORY:
;       Support
; 
; CALLING SEQUENCE:
;       result=LINENUMBER()
;
;
; INPUTS:
;       -
;
; OPTIONAL INPUTS:
;       - 
;
; KEYWORD PARAMETERS:
;       -
;
; OUTPUTS:
;       result = the line number of the program that is being run. If
;       in main it will return (You are in the main)
;
; OPTIONAL OUTPUTS:
;       -
; 
; PROCEDURES CALLED:
;       STR_SEP(), STRTRIM(), STRCOMPRESS()
;
; EXAMPLE:
;      
;
; MODIFICATION HISTORY:
;       Written by P.Kamphuis 01-01-2015 
;
; NOTE:
;     
;-

help,calls=cc		; Get list of all calls.
tmp= str_sep(strtrim(strcompress(cc[n_elements(cc)-2]),2),'(')
if n_elements(tmp) GT 1 then begin
   tmp3=str_sep(strtrim(strcompress(tmp[0]),2),' ')
   IF tmp3[0] NE 'LINENUMBER' then begin
      tmp2= str_sep(strtrim(strcompress(tmp[1]),2),')')
      Linenumber=tmp2[0]
   ENDIF ELSE BEGIN
      Linenumber='You are in the main'
   ENDELSE
endif else begin
   Linenumber='You are in the main'
Endelse
return,' ('+Linenumber+') '

end
