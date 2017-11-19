FUNCTION FAT_GDLGAUSS, axis, profile


;+
; NAME:
;       FAT_GDLGAUSS
;
; PURPOSE:
;       Program to simulate GAussfit in the GDL environment
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
  
error=REPLICATE(0.5,n_elements(axis))
start=[MAX(profile),MEAN(axis),MAX(axis)/2.]
coeff = MPFITFUN('FATGAUSS', axis, profile, error, start,/QUIET)
f=coeff[0]*EXP(-0.5*((axis[*]-coeff[1])/coeff[2])^2)

return,f
end

FUNCTION FATGAUSS, X, P
  RETURN, P[0] *EXP(-0.5*((X-P[1])/P[2])^2)
END
