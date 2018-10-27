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
;       17-10-2018 P.Kamphuis ; Made a better estimate of sigma and center  
;       Written 02-06-2016 P.Kamphuis v1.0
;
; NOTE:
;
;     
;-
  
  error=REPLICATE(0.5,n_elements(axis))
  tmp=WHERE(profile GT MAX(profile)/2.)
  sig=(axis[tmp[n_elements(tmp)-1]]-axis[tmp[0]])/(2*SQRT(2*ALOG(2.)))
  start=[MAX(profile),axis[fix(n_elements(axis)/2.)],sig]
  coeff = MPFITFUN('FATGAUSS', axis, profile, error, start,/quiet)
  f=coeff[0]*EXP(-0.5*((axis[*]-coeff[1])/coeff[2])^2)

  return,f
end

FUNCTION FATGAUSS, X, P
  RETURN, P[0] *EXP(-0.5*((X-P[1])/P[2])^2)
END
