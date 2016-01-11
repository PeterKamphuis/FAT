Pro set_warp_slopev2,SBRarr,SBRarr2,cutoff,INCLinput2,PAinput2,INCLinput3,PAinput3,norings,log=log,innerfix=innerfix

;+
; NAME:
;       SET_WARP_SLOPEV2
;
; PURPOSE:
;       This little routine takes the Surface brightness profiles and
;       determines where they are below the cutoff values and sets
;       those rings  to be only fitted with a slope in the Inclination and PA parameters
;
; CATEGORY:
;       Support
; 
; CALLING SEQUENCE:
;       SET_WARP_SLOPEV2,SBRarr,SBRarr2,cutoff,INCLinput2,PAinput2,INCLinput3,PAinput3,norings,LOG=log,INNERFIX=innerfix
;
;
; INPUTS:
;       SBRarr = approaching sbr array
;       SBRarr2 = receding sbr array
;       cutoff = cutoff values for the nodes
;       norings = number of rings in model
;       INCLinput2-3 = Template for fitting parameters for inclination
;       PAinput2-3 = Template for fitting parameters for pa
;
; OPTIONAL INPUTS:
;       LOG = name of the tracing log 
;       INNERFIX = number of rings that should be fitted as one in the
;       central parameters
;
; KEYWORD PARAMETERS:
;       -
;
; OUTPUTS:
;       INCLinput2-3 = Updated fitting parameters for inclination
;       PAinput2-3 = Updated fitting parameters for pa

;
; OPTIONAL OUTPUTS:
;       -
; 
; PROCEDURES CALLED:
;       STR_SEP(), STRTRIM(), STRCOMPRESS(), STRING()
;
; EXAMPLE:
;      
;
; MODIFICATION HISTORY:
;       Written 01-01-2015 P.Kamphuis v1.0
;
; NOTE:
;     
;-
  COMPILE_OPT IDL2 
  
  warpconstused=n_elements(SBRarr)	   
  if n_elements(innerfix) NE 0 then begin
     INCLinput2[0]='!INCL '+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+':'+strtrim(strcompress(string(innerfix+1,format='(F7.4)')),1)
     INCLinput3[0]='!INCL_2 '+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+':'+strtrim(strcompress(string(innerfix+1,format='(F7.4)')),1)
     PAinput2[0]='!PA '+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+':'+strtrim(strcompress(string(innerfix+1,format='(F7.4)')),1)
     PAinput3[0]='!PA_2 '+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+':'+strtrim(strcompress(string(innerfix+1,format='(F7.4)')),1)
  ENDIF ELSE BEGIN
     INCLinput2[0]='!INCL '+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+':5'
     INCLinput3[0]='!INCL_2 '+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+':5'
     PAinput2[0]='!PA '+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+':5'
     PAinput3[0]='!PA_2 '+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+':5'
  ENDELSE
  INCLinput2[10]=' '
  INCLinput3[10]=' '
  PAinput2[10]=' '
  PAinput3[10]=' '
                                ;Set warp slope of this fit
                                ;And we determine which part of the
                                ;warp should be a slope. i.e. when SBR
                                ;is lower than cutoff.
  get_newringsv8,SBRarr,SBRarr2,cutoff,warpconstused,/individual
  IF warpconstused[0] LT 6 then warpconstused[0]=6
  IF warpconstused[1] LT 6 then warpconstused[1]=6
  IF n_elements(innerfix) NE 0 AND warpconstused[0] LT innerfix+2 then warpconstused[0]=innerfix+2
  IF n_elements(innerfix) NE 0 AND warpconstused[1] LT innerfix+2 then warpconstused[1]=innerfix+2
  IF warpconstused[0] LT norings[0]-1 then begin
     INCLinput2[10]='INCL '+strtrim(strcompress(string(norings[0]-1,format='(I3)')),1)+':'+strtrim(strcompress(string(warpconstused[0]+1,format='(I3)')),1)
     PAinput2[10]='PA '+strtrim(strcompress(string(norings[0]-1,format='(I3)')),1)+':'+strtrim(strcompress(string(warpconstused[0]+1,format='(I3)')),1)
  ENDIF
  IF warpconstused[1] LT norings[0]-1 then begin
     INCLinput3[10]='INCL_2 '+strtrim(strcompress(string(norings[0]-1,format='(I3)')),1)+':'+strtrim(strcompress(string(warpconstused[1]+1,format='(I3)')),1)
     PAinput3[10]='PA_2 '+strtrim(strcompress(string(norings[0]-1,format='(I3)')),1)+':'+strtrim(strcompress(string(warpconstused[1]+1,format='(I3)')),1)
  ENDIF
  IF size(log,/TYPE) EQ 7 then begin
     openu,66,log,/APPEND
     printf,66,linenumber()+"SET_WARP_SLOPEV2: We fix the warp on side 1 from ring "+strtrim(string(warpconstused[0]),2)+" on and side 2 from ring "+strtrim(string(warpconstused[1]),2)+" on."
     close,66
  ENDIF 
end
