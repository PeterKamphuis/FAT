Pro obtain_pav2,map,PA,INCLINATION=inclination,CENTER=center,NOISE=noise,ITERATIONS=norings

;+
; NAME:
;       OBTAIN_PAV2
;
; PURPOSE:
;       Routine to fit five ellipses to the regions defined by
;       the ten area's between the maximum and the noise/minimum
;
; CATEGORY:
;       Support
; 
; CALLING SEQUENCE:
;       OBTAIN_PAV2,map,PA,INCLINATION=inclination,CENTER=center,NOISE=noise,ITERATIONS=norings
;
;
; INPUTS:
;       map = moment 0 map of the galaxy
;
; OPTIONAL INPUTS:
;       - 
;
; KEYWORD PARAMETERS:
;       INCLINATION = first guess at the inclination from the
;       ellipse's axis ratios.
;       CENTER = center of the galaxy in pixels
;       NOISE = estimate of the noise of the moment 0 map
;       ITERATIONS = The amount of rings to be fitted, the default is 5.
;
; OUTPUTS:
;       PA = the positional angle
;
; OPTIONAL OUTPUTS:
;       -
; 
; PROCEDURES CALLED:
;       FIT_ELLIPSE()
;
; EXAMPLE:
;      
;
; MODIFICATION HISTORY:
;       16-03-2017 P.Kamphuis; Fit ellipse crashes when only one index
;                              is presented to it hence the condition has been set for more
;                              than one index to be present    
;       18-02-2016 P.Kamphuis; Replaced sigma with STDDEV   
;       Written 01-01-2015 P.Kamphuis v1.0
;
; NOTE:
;     
;- 
  COMPILE_OPT IDL2 
  tmp0=WHERE(map GT 0.)
  mapmax=MAX(map[tmp0],MIN=mapmin)
  IF n_elements(noise) LT 1 then noise=mapmin 
  IF n_elements(norings) EQ 0 then norings=5.
  maxint=mapmax-mapmax/2.
  step=(maxint-(3.*noise))/norings
  pafound=dblarr(norings)
  inclfound=dblarr(norings)
  for i=0,norings-1 do begin
     minint=maxint-step
     tmp=WHERE(map GT minint AND map LT maxint)
     ;fit_ellipse does not work when we have only one element
     IF n_elements(tmp) GT 1 then begin
        IF n_elements(center) EQ 0 then begin
           ell=fit_ellipse(tmp,xsize=n_elements(map[*,0]),ysize=n_elements(map[0,*]),orientation=paring,axes=inclring)
           pafound[i]=paring
           ratio=double(inclring[1]/inclring[0])
           IF ratio GT 1. then ratio=1.
           IF ratio LT 0.21 then ratio=0.21
           inclfound[i]=acos(SQRT((ratio^2-0.2^2)/0.96))*180./!pi+2.
        ENDIF else begin
           ell=fit_ellipse(tmp,center=center,xsize=n_elements(map[*,0]),ysize=n_elements(map[0,*]),orientation=paring,axes=inclring)
           pafound[i]=paring
           ratio=double(inclring[1]/inclring[0])
           IF ratio GT 1. then ratio=1.
           IF ratio LT 0.21 then ratio=0.21
           inclfound[i]=acos(SQRT((ratio^2-0.2^2)/0.96))*180./!pi+2.
        ENDELSE
     ENDIF
     maxint=minint
  endfor
  tmp=WHERE(pafound LT 0)
  IF tmp[0] NE -1 then pafound[tmp]=pafound[tmp]+360
  for i=1,n_elements(pafound)-1 do begin
     diff=pafound[i-1]-pafound[i]
     IF diff GT 120. then pafound[i]=pafound[i]+180
     IF diff LT -120 then pafound[i]=pafound[i]-180
  endfor
  tmp=WHERE(pafound NE 0.) 
  PA=dblarr(2)
  inclination=dblarr(2)
  IF tmp[0] NE -1 then begin
     
     PA[0]=TOTAL(pafound[tmp])/n_elements(pafound[tmp])+90.
     IF PA[0] LT 0 then PA[0]=PA[0]+360
     IF PA[0] GT 360 then PA[0]=PA[0]-360
     PA[1]=STDDEV(pafound[tmp])

     tmp2=inclfound[2:n_elements(inclfound)-1]
     tmp=WHERE(tmp2 NE 0)
     IF tmp[0] NE -1 then begin 
        inclination[0]=TOTAL(tmp2[tmp])/(n_elements(tmp))
        inclination[1]=STDDEV(tmp2[tmp])
     ENDIF

  ENDIF
  IF inclination[0] EQ 0. then begin
     inclination[1]=45.
     inclination[0]=45.
  ENDIF

  IF PA[0] EQ 0. then begin
     PA[1]=180.
     PA[0]=90.
  ENDIF                         

END
