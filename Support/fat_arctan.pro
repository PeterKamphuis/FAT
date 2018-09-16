FUNCTION FAT_ARCTAN, axis, profile,error=error,gdlidl=gdlidl


;+
; NAME:
;       FAT_ARCTAN
;
; PURPOSE:
;       Program to fit an arctan
;
; CATEGORY:
;       Support
; 
; CALLING SEQUENCE:
;       FAT_GDLGAUSS, axis, profile, coeff
;
;
; INPUTS:
;       axis = xaxis of the profile to be fitted
;    profile = yaxis of the profile to be fitted
; OPTIONAL INPUTS:
;       - 
;
; KEYWORD PARAMETERS:
;       -
; 
; OUTPUTS:
;       -
;
; OPTIONAL OUTPUTS:
;       -
; 
; PROCEDURES CALLED:
;       EXP(), MAX(), MEAN(), MPFITFUN(), REPLICATE()
;
; EXAMPLE:
;      
;
; MODIFICATION HISTORY:
;       Written 02-06-2016 P.Kamphuis v1.0
;
; NOTE:
;
;     
;-
  
if n_elements(error) EQ 0. then error=REPLICATE(0.5,n_elements(axis))
start=[MEAN(axis)/3.,0.,4.,MEAN(Profile)]
if gdlidl EQ 1 then begin
   coeff = MPFITFUN('FATARCTAN', axis, profile, error, start,/QUIET)
endif else begin
   coeff = CURVEFIT(axis, profile, error, start,FUNCTION_NAME='FATARCTAN')
ENDELSE
if coeff[1] LT 0. then coeff[1]=abs(coeff[1])
C=abs(axis[fix(n_elements(axis)/3.*2)]-axis[fix(n_elements(axis)/3.)])

profile=-1.*ATAN((axis-coeff[0])/(C+coeff[1]))/!pi*abs(coeff[2])+coeff[3]

return,coeff
end

FUNCTION FATARCTAN, X, P
  C=abs(X[fix(n_elements(x)/3.*2)]-X[fix(n_elements(x)/3.)])
  RETURN, -1.*ATAN((X-P[0])/(C+abs(P[1])))/!pi*abs(P[2])+P[3] 
END
