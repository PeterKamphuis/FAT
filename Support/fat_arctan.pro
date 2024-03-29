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
;       Written 10-10-2018 P.Kamphuis v1.0
;
; NOTE:
;
;     
;-
  
if n_elements(error) EQ 0. then error=REPLICATE(0.5,n_elements(axis))
start=[MEAN(axis)/3.,0.,4.,MEAN(Profile)]
coeff = MPFITFUN('FATARCTAN', axis, profile, error, start,/QUIET)
if coeff[1] LT 0. then coeff[1]=abs(coeff[1])
;C=abs(axis[fix(n_elements(axis)/3.*1.5)]-axis[fix(n_elements(xaxis)/3.)])
C=axis[n_elements(axis)-1]*0.1
profile=-1.*ATAN((axis-coeff[0])/(C+coeff[1]))/!pi*abs(coeff[2])+coeff[3]

return,coeff
end

FUNCTION FATARCTAN, X, P
  ;C=abs(X[fix(n_elements(x)/3.*1.5)]-X[fix(n_elements(x)/3.)])
  C=X[n_elements(X)-1]*0.1
  RETURN, -1.*ATAN((X-P[0])/(C+abs(P[1])))/!pi*abs(P[2])+P[3] 
END
