FUNCTION obtain_ratios,angles,map,center=center,MAJ_AXIS=tmpwidth,gdlidl=gdlidl,noise=noise,beam=beam,debug=debug

;+
; NAME:
;       OBTAIN_RATIOS
;
; PURPOSE:
;       A function to get the minor and major axis ratio under a set
;       of angles
;
; CATEGORY:
;       Support
; 
; CALLING SEQUENCE:
;       OBTAIN_RATIOS,angles ,map,center=center,MAJ_AXIS=tmpwidth,gdlidl=gdlidl,noise=noise,beam=beam,/debug

;
;
; INPUTS:
;       angles: The angles to check
;       map: the moment0 map  
;
; OPTIONAL INPUTS:
;       - 
;
; KEYWORD PARAMETERS:
;      center = central coordinates in pixels  
;      NOISE = noise in the momentmap (This is required)
;      BEAM = Beam of the observations in pixels. If unset it is thought to be
;       1 pixel beam, if a single value a circular beam is assumed
;      GDLIDL = Ar we doing idl or gdl  
;      /DEBUG - flag for getting additional output to screen
;
; OUTPUTS:
;      MAJ_AXIS = the major axis width in pixels
;
; OPTIONAL OUTPUTS:
;       -
; 
; PROCEDURES CALLED:
;       INT_PROFILEV2, INTERPOLATE, GAUSSFIT(), FAT_GAUSS
;
; EXAMPLE:
;      
;
; MODIFICATION HISTORY: 
;       Written 01-01-2019 P.Kamphuis v1.0
;
; NOTE:   
;-

  IF n_elements(beam) EQ 0. then beam = [1.,1]
  IF n_elements(gdlidl) EQ 0. then gdlidl=0.
  IF n_elements(center) EQ 0 then center=[fix(n_elements(map[*,0])/2.),fix(n_elements(map[0,*])/2.)]
  if n_elements(noise) EQ 0. then begin
     tmp=WHERE(map GT 0.)
     noise= 2.* MIN(map[tmp])
  ENDIF
  if keyword_set(debug) then begin
     print,'This is the input'
     print,'beam'
     print,beam
     print,'center'
     print,center
     print,'noise'
     print,noise
     print,'map size and total'
     print,n_elements(map[*,0]),n_elements(map[0,*]),TOTAL(map)
     print,'Angles to check'
     print,angles
  endif

  
  ratios=dblarr(n_elements(angles))
  tmpwidth=dblarr(n_elements(angles))
                                ;Experiment
  axis1=findgen(beam[0]*40.)/10.-beam[0]*2.
  psf1=(1./(2.*!pi*(beam[0]*10.)^2))*EXP(-((0.5*axis1^2)/((beam[0]*10.)^2)))
  psf2=(1./(2.*!pi*(beam[1]*10.)^2))*EXP(-((0.5*axis1^2)/((beam[1]*10.)^2)))

  for i=0,n_elements(angles)-1 do begin
                                ;first get the profiles from the image
     int_profilev2,map,xprofile,pa=angles[i],xcenter=center[0],ycenter=center[1]
     int_profilev2,map,yprofile,pa=angles[i],xcenter=center[0],ycenter=center[1],/minor
    
                                ;Fit a Gaussian to reduce noise effects
                                ;get the maximum
     axis=findgen(n_elements(xprofile))
     IF gdlidl then begin
        xfit = FAT_GDLGAUSS(axis, xprofile)
        yfit = FAT_GDLGAUSS(axis, yprofile)
     ENDIF ELSE BEGIN
        xfit = GAUSSFIT(axis, xprofile, coeff, NTERMS=3)
        yfit = GAUSSFIT(axis, yprofile, coeff, NTERMS=3)
     ENDELSE
     newaxis=findgen(n_elements(xprofile)*10.)/10.
     xprof=1
     yprof=1

     interpolate,xfit,axis,newradii=newaxis,output=xprofile
     interpolate,yfit,axis,newradii=newaxis,output=yprofile 
     
     maxxprof=MAX(xprofile)
     maxyprof=MAX(yprofile)
                                ;and the limits
     limit=2.*noise
     xFWHM=0.
     if limit LT maxxprof and limit LT maxyprof  then begin
        tmp=WHERE(xprofile GT limit)
        
        IF tmp[n_elements(tmp)-1]-tmp[0] LT beam[1]*10. then xFWHM=beam[0]*10. else $
           xFWHM=SQRT(double(tmp[n_elements(tmp)-1]-tmp[0])^2-2.*(beam[0]*10.)^2)
       
        tmp=WHERE(yprofile GT limit)
        IF tmp[n_elements(tmp)-1]-tmp[0] LT beam[1]*10. then yFWHM=beam[1]*10. else $
           yFWHM=SQRT(double(tmp[n_elements(tmp)-1]-tmp[0])^2-2.*(beam[1]*10.)^2)
                                ;IF the sizes are very close to
                                ;the beam sizes let's but a bit
                                ;more effort  into deconvolution
        ;print,xFWHM,yFWHM,beam[0]*10., beam[1]*10.
        if xFWHM/(beam[0]*10.) LT 5. then begin
           ;print,'Old xFWHM',xFWHM
           xFWHM = xFWHM- (beam[0]*10)^2/(2.*xFWHM)
           ;print,'New xFWHM',xFWHM
        endif
        if yFWHM/(beam[1]*10.) LT 5. then begin
           ;print,'Old yFWHM',yFWHM
           yFWHM = yFWHM- (beam[1]*10.)^2/(2.*yFWHM)
           ;print,'New yFWHM',yFWHM
        endif
        if keyword_set(debug) then begin
           print,'ANGLE is '+string(angles[i])
           print,'Found FWHM'
           print,xFWHM,yFWHM
        endif
        ratios[i]=double(yFWHM/xFWHM)
     endif else begin
        if keyword_set(debug) then print,'Limits are too high'
        ratios[i]=0.
     endelse
     ;This is supposedly the FWHM so let's take 3sigma
     tmpwidth[i]=(xFWHM/(2.*SQRT(2*ALOG(2))))*3.
     
     
  endfor
 
  return,ratios
end
