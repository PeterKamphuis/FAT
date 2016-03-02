Pro get_fixedringsv8,Parametersin,rings

;+
; NAME:
;       GET_FIXEDRINGSV8
;
; PURPOSE:
;       Little routine to determine how many rings of the parameters
;       should be considered as the flat inner part of the model based
;       on the variation in INCL and PA on both sides
;
; CATEGORY:
;       Support
; 
; CALLING SEQUENCE:
;       GET_FIXEDRINGSV8,Parametersin,rings
;
;
; INPUTS:
;       Parametersin = The INCL and PA parameters
;
; OPTIONAL INPUTS:
;       - 
;
; KEYWORD PARAMETERS:
;       -
;
; OUTPUTS:
;       rings = the amount of rings that should remain fixed
;
; OPTIONAL OUTPUTS:
;       -
; 
; PROCEDURES CALLED:
;       MAX(), STDDEV(), ROBUST_SIGMA(), FLOOR()
;
; EXAMPLE:
;      
;
; MODIFICATION HISTORY:
;       18-02-2016 P.Kamphuis; Replaced sigma with STDDEV   
;       Written by P.Kamphuis 01-01-2015 
;
; NOTE:
;     
;-
  COMPILE_OPT IDL2
  Parameters=Parametersin
  indrings=dblarr(n_elements(Parameters[0,*]))
  half=fix(n_elements(Parameters[*,0])/2)
  for j=0,n_elements(Parameters[0,*])-1 do begin
     newPA=dblarr(n_elements(Parameters[*,0]))
     newPA[0]=(Parameters[0,j]+Parameters[1,j])/2.
     for i=1,n_elements(Parameters[*,j])-2 do begin
        newPA[i]=(Parameters[i-1,j]+Parameters[i,j]+Parameters[i+1,j])/3.
     endfor
     newPA[n_elements(Parameters[*,0])-1]=(Parameters[n_elements(Parameters[*,0])-2,j]+Parameters[n_elements(Parameters[*,0])-1,j])/2.
     MAXdiff=MAX(Parameters[*,j]-newPA)
     IF  n_elements(Parameters[0:half,j]) LT 8. then begin
        rms=ROBUST_SIGMA(Parameters[*,j]-newPA)
     ENDIF ELSE BEGIN
        rms=(STDDEV(Parameters[*,j]-newPA)+6*ROBUST_SIGMA(Parameters[*,j]-newPA))/7.
     ENDELSE
     IF rms GT MAXdiff/4. and rms GT 10. then begin
        newPA2=dblarr(n_elements(newPA[*]))
        newPA2[0]=(newPA[0]+newPA[1])/2.
        for i=1,n_elements(newPA[*])-2 do begin
           newPA2[i]=(newPA[i-1]+newPA[i]+newPA[i+1])/3.
        endfor
        newPA2[n_elements(newPA[*])-1]=(newPA[n_elements(newPA[*])-2]+newPA[n_elements(newPA[*])-1])/2.
        IF n_elements(Parameters[0:half,j]) LT 8. then begin
           rms=ROBUST_SIGMA(newPA-newPA2)
        ENDIF ELSE BEGIN
           rms=(STDDEV(newPA-newPA2)+6*ROBUST_SIGMA(newPA-newPA2))/7.
        ENDELSE
     ENDIF
     av1=Parameters[0,j]
     av2=MEAN(Parameters[0:half,j])
     IF ABS(av1-av2) GT rms AND n_elements(Parameters[0:half,j]) LT 8. then av=av1 else av=av2
     tmp=WHERE(ABS(Parameters[*,j]-av) GT rms)
     IF n_elements(tmp) GT 1 then begin
        for i=0,n_elements(tmp)-2 do begin
           IF tmp[i+1]-tmp[i] EQ 1 then begin
              indrings[j]=tmp[i]-1
              break
           ENDIF ELSE BEGIN
              IF i LT n_elements(tmp)-3 then begin
                 IF tmp[i+1]-tmp[i] EQ 2 AND tmp[i+2]-tmp[i+1] EQ 1 then begin
                    indrings[j]=tmp[i]-1
                    break
                 ENDIF
              ENDIF
           ENDELSE
        endfor
     ENDIF ELSE begin
        indrings[j]=n_elements(Parameters[*,j])-1
     ENDELSE
  endfor
  rings=floor(TOTAL(double(indrings))/(n_elements(indrings)))-1
  if rings LT 3 then rings=3
  IF rings GT fix(n_elements(Parameters[*,0])/2.) then begin
     IF fix(n_elements(Parameters[*,0])/2.) GE 3 then rings=fix(n_elements(Parameters[*,0])/2.) else rings=3
  ENDIF
end
