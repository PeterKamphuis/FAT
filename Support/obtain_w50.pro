Pro obtain_w50,Cube,mask,hed,W50

;+
; NAME:
;       OBTAIN_W50
;
; PURPOSE:
;       Routine to obtain W50 from a 3D fits cube array
;
; CATEGORY:
;       Support
; 
; CALLING SEQUENCE:
;       OBTAIN_W50,Cube,mask,hed,W50
;
;
; INPUTS:
;       Cube =  Array of the cube to be analysed
;       Mask = Mask to identify the emission in the cube
;       hed = the header of the cube
;       
; OPTIONAL INPUTS:
;       - 
;
; KEYWORD PARAMETERS:
;       -
;
; OUTPUTS:
;       W50 = the velocity width at 50% from the peak value in the
;       spectrum in km/s
;
; OPTIONAL OUTPUTS:
;       -
; 
; PROCEDURES CALLED:
;       BUILDAXII, STRUPCASE(), STRTRIM(), SXPAR() 
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
  buildaxii,hed,xaxis,yaxis,zaxis=zaxis
  IF STRUPCASE(strtrim(sxpar(hed,'CUNIT3'),2)) EQ 'M/S' then divide=1000. else divide=1. 
  zaxis=zaxis/divide
  profile=dblarr(n_elements(zaxis))
  maskedCube=dblarr(n_elements(Cube[*,0,0]),n_elements(Cube[0,*,0]),n_elements(Cube[0,0,*]))
  maskedCube[WHERE(mask EQ 1)]=Cube[WHERE(mask EQ 1)]
  for i=0,n_elements(zaxis)-1 do begin
     profile[i]=TOTAL(maskedCube[*,*,i])
  endfor


  maxprof=MAX(profile)

  above50=WHERE(profile GT maxprof/2.)
  if n_elements(above50) GT 2 then W50=zaxis[above50[n_elements(above50)-1]]-zaxis[above50[0]] else W50=sxpar(hed,'CDELT3')/divide


end
