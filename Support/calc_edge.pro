Pro calc_edge,noise,radii,bm,cutoffrings

;+
; NAME:
;       CALC_EDGE
;
; PURPOSE:
;       subroutine to calculate whether a ring can be reliably fitted or not
;
; CATEGORY:
;       Support
; 
; CALLING SEQUENCE:
;       CALC_EDGE,noise,radii,bm,cutoffringss
;
;
; INPUTS:
;       noise = is a scalar with the rms of the cube to be fitted 
;       radii = an array with the nodes of the model in arcsec
;       bm    = a two dimensional array with the FWHM of the major and
;       minor axis of the beam of the observation
;
; OPTIONAL INPUTS:
;       -
;
; KEYWORD PARAMETERS:
;       -
;
; OUTPUTS:
;       cutoffrings = the values at which nodes with an sbr that is
;       lower is no longer reliable. 
;
; OPTIONAL OUTPUTS:
;       -  
;
; PROCEDURES CALLED:
;       ALOG()
;
;
; EXAMPLE:
;
;
; MODIFICATION HISTORY:
;     Written by P.Kamphuis 01-01-2015 
; NOTE:
;     We take the sensitivity limit provided in equation 4 of Josza et
;     al. However this limit is independent of the actual noise of the
;     cube which cannot be right. Therefore we scale the limit by the
;     ratio noise/3.6 mJy/beam as 3.6 mJy/beam is the noise used in
;     their test cubes. Specifics of Josza et al. 2007 test cubes
;     vwidth=4.12 km/s FHWM Beam=12*14 noise=3.6 mJy/beam most likely
;     the actual maj beam was 14.2 which gives a 202 arcsec^2 solid
;     angle as stated in the paper check wether we need the area correction factor but I don't think so
;-
  COMPILE_OPT IDL2 

ratio=noise/0.0036
beamsolid=(!pi*bm[0]*bm[1])/(4.*ALOG(2.))
ringarea=dblarr(n_elements(radii))
IF RADII[0] GT 0 then begin
   ringarea[0]=!pi*((radii[0]+radii[1])/2.)^2
endif else ringarea[0]=0.
for i=1,n_elements(radii)-2 do begin
   ringarea[i]=!pi*(((radii[i]+radii[i+1])/2.)^2-((radii[i]+radii[i-1])/2.)^2)  
endfor
ringarea[n_elements(radii)-1]=!pi*((radii[n_elements(radii)-1]+0.5*(radii[n_elements(radii)-1]-radii[n_elements(radii)-2]))^2)-!pi*(((radii[n_elements(radii)-1]+radii[n_elements(radii)-2])/2.)^2)
cutoffrings=9E-4*(ringarea/beamsolid)^(-0.82)*ratio
if ringarea[0] EQ 0. then cutoffrings[0]=0.


end
