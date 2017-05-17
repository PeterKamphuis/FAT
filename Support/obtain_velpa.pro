Pro obtain_velpa,map,velpa,CENTER=center,debug=debug,intensity=intensity

;+
; NAME:
;       OBTAIN_VELPA
;
; PURPOSE:
;       Routine to obtain the PA from a moment 1 map of a galaxy
;
; CATEGORY:
;       Support
; 
; CALLING SEQUENCE:
;       OBTAIN_VELPA,map,velpa,CENTER=center
;
;
; INPUTS:
;       map = moment 1 map of the galaxy
;
; OPTIONAL INPUTS:
;       - 
;
; KEYWORD PARAMETERS:
;       CENTER = center of the galaxy in pixels
;    INTENSITY = Moment0 map to identify edges.   
;        DEBUG = option to have verbose running
;  
; OUTPUTS:
;       PA = the positional angle
;
; OPTIONAL OUTPUTS:
;       -
; 
; PROCEDURES CALLED:
;       FAT_GAUSS_SMOOTH(), WRITEFITS
;
; EXAMPLE:
;      
;
; MODIFICATION HISTORY:
;       15-05-2017 P.Kamphuis; Added the intensity keyword which
;                              allows for a moment 0 map to be
;                              given. If given the velocity field will
;                              be cut at 5% of the maximum intensity.   
;       21-03-2017 P.Kamphuis; Added debug keyword for debugging  
;       18-02-2016 P.Kamphuis; Replaced sigma with STDDEV   
;       06-01-2016 P.Kamphuis; Added NE 0. Condition to detecting the
;                              min and max values   
;       29-12-2016 P.Kamphuis: Replaced GAUSS_SMOOTH and SMOOTH with FAT_SMOOTH
;                              for increased compatibility with IDL 7.0 and GDL 
;       Written 01-01-2015 P.Kamphuis v1.0
;
; NOTE:
;     
;-
  COMPILE_OPT IDL2 
  velpa=dblarr(2)
  if n_elements(intensity) NE 0 then begin
     tmp=WHERE(intensity LT 0.05*MAX(intensity))
     map[tmp]=!values.f_nan
  endif

  If n_elements(center) EQ 0 then center=[0,0]
  tmpor=WHERE(FINITE(map))
  ormax=MAX(map[tmpor],min=ormin)
  if keyword_set(debug) then begin
     print,'These are the maximum and minimum velocities in the original field'
     print,ormax,ormin
  ENDIF
  smoothfield=map
  pixsmooth=3
smoothagain:
  smoothfield=FAT_SMOOTH(smoothfield,pixsmooth,/GAUSSIAN)
  tmp=WHERE(FINITE(smoothfield) AND smoothfield NE 0.)
  IF tmp[0] NE -1 AND n_elements(tmp) GT 10. then begin
     maxvel=MAX(smoothfield[tmp],min=minvel)
     if keyword_set(debug) then begin
        print,'These are the maximum and minimum velocities in the smoothed field'
        print,maxvel,minvel
     ENDIF
     POS1=WHERE(maxvel EQ smoothfield)
     POS2=WHERE(minvel EQ smoothfield)
   ;  print,n_elements(POS1),n_elements(POS2)
     IF n_elements(POS1) GT 1 or n_elements(POS2) GT 1 then begin
        POS1=POS1[0]
        POS2=POS2[0]
     ENDIF
     ;goto,smoothagain else begin
        s = SIZE(smoothfield)
        ncol = s[1]
        x1 = POS1 MOD ncol
        y1 = POS1 / ncol
        x2 = POS2 MOD ncol
        y2 = POS2 / ncol
        if keyword_set(debug) then begin
           print,'These are the location of the maximum and minimum'
           print,'max',x1,y1
           print,'min',x2,y2
        ENDIF
        velpamax=ATAN((center[0]-x1)/(center[1]-y1))
        velpamin=ATAN((center[0]-x2)/(center[1]-y2))
        velnocen=ATAN(double(x1-x2)/double(y1-y2))
        case (1) of
           center[0]-x1 LT 0 AND velpamax LT 0:velpamax=ABS(velpamax)+180.*!DtoR
           center[0]-x1 GT 0 AND velpamax LT 0:velpamax=ABS(velpamax)
           center[0]-x1 LT 0 AND velpamax GT 0:velpamax=360.*!DtoR-(velpamax)
           center[0]-x1 GT 0 AND velpamax GT 0:velpamax=180.*!DtoR-(velpamax)
           center[0]-x1 EQ 0 and center[1]-y1 LT 0:velpamax=0.
           center[0]-x1 EQ 0 and center[1]-y1 GT 0:velpamax=180.*!DtoR
           else:velpamax=!values.f_nan
        endcase
        case (1) of
           center[0]-x2 LT 0 AND velpamin LT 0:velpamin=ABS(velpamin)
           center[0]-x2 GT 0 AND velpamin LT 0:velpamin=ABS(velpamin)+180.*!DtoR
           center[0]-x2 LT 0 AND velpamin GT 0:velpamin=180.*!DtoR-(velpamin)
           center[0]-x2 GT 0 AND velpamin GT 0:velpamin=360.*!DtoR-(velpamin)
           center[0]-x2 EQ 0 and center[1]-y2 LT 0:velpamin=180.*!DtoR
           center[0]-x2 EQ 0 and center[1]-y2 GT 0:velpamin=0.
           else:velpamin=!values.f_nan
        endcase
        case (1) of
           velnocen EQ 0 AND center[1]-y2 GT 0: velnocen=0.
           velnocen EQ 0 AND center[1]-y2 LT 0: velnocen=180.*!DtoR
           velnocen LT 0 AND x2-x1 GT 0: velnocen=ABS(velnocen)
           velnocen LT 0 AND x2-x1 LT 0: velnocen=ABS(velnocen)+180.*!DtoR
           velnocen GT 0 AND x2-x1 GT 0: velnocen=180.*!DtoR-velnocen
           velnocen GT 0 AND x2-x1 LT 0: velnocen=ABS(velnocen)+180.*!DtoR
           else:velnocen=!values.f_nan
        endcase
        
     

  ENDIF else begin
     IF pixsmooth GT 2 then begin
        pixsmooth=2
        smoothfield=map
        goto,smoothagain
     ENDIF ELSE begin
        smoothfield=map
        smoothfield=FAT_SMOOTH(smoothfield,3,/BOX)
        tmp=WHERE(FINITE(smoothfield) AND smoothfield NE 0.)      
        maxvel=MAX(smoothfield[tmp],min=minvel)
        if keyword_set(debug) then begin
           print,'These are the maximum and minimum velocities in the Box field'
           print,maxvel,minval
        ENDIF
        POS1=WHERE(maxvel EQ smoothfield)
        POS2=WHERE(minvel EQ smoothfield)
        IF n_elements(POS1) GT 1 or  n_elements(POS2) GT 1 then begin
           velpa=[!values.f_nan,!values.f_nan]
           goto,endthis
        ENDIF else begin
           s = SIZE(smoothfield)
           ncol = s[1]
           x1 = POS1 MOD ncol
           y1 = POS1 / ncol
           x2 = POS2 MOD ncol
           y2 = POS2 / ncol
           
           velpamax=ATAN((center[0]-x1)/(center[1]-y1))
           velpamin=ATAN((center[0]-x2)/(center[1]-y2))
           velnocen=ATAN(double(x1-x2)/double(y1-y2))
           case (1) of
              center[0]-x1 LT 0 AND velpamax LT 0:velpamax=ABS(velpamax)+180.*!DtoR
              center[0]-x1 GT 0 AND velpamax LT 0:velpamax=ABS(velpamax)
              center[0]-x1 LT 0 AND velpamax GT 0:velpamax=360.*!DtoR-(velpamax)
              center[0]-x1 GT 0 AND velpamax GT 0:velpamax=180.*!DtoR-(velpamax)
              center[0]-x1 EQ 0 and center[1]-y1 LT 0:velpamax=0.
              center[0]-x1 EQ 0 and center[1]-y1 GT 0:velpamax=180.*!DtoR
              else:velpamax=!values.f_nan
           endcase
           case (1) of
              center[0]-x2 LT 0 AND velpamin LT 0:velpamin=ABS(velpamin)
              center[0]-x2 GT 0 AND velpamin LT 0:velpamin=ABS(velpamin)+180.*!DtoR
              center[0]-x2 LT 0 AND velpamin GT 0:velpamin=180.*!DtoR-(velpamin)
              center[0]-x2 GT 0 AND velpamin GT 0:velpamin=360.*!DtoR-(velpamin)
              center[0]-x2 EQ 0 and center[1]-y2 LT 0:velpamin=180.*!DtoR
              center[0]-x2 EQ 0 and center[1]-y2 GT 0:velpamin=0.
              else:velpamin=!values.f_nan
           endcase
           case (1) of
              velnocen EQ 0 AND center[1]-y2 GT 0: velnocen=0.
              velnocen EQ 0 AND center[1]-y2 LT 0: velnocen=180.*!DtoR
              velnocen LT 0 AND x2-x1 GT 0: velnocen=ABS(velnocen)
              velnocen LT 0 AND x2-x1 LT 0: velnocen=ABS(velnocen)+180.*!DtoR
              velnocen GT 0 AND x2-x1 GT 0: velnocen=180.*!DtoR-velnocen
              velnocen GT 0 AND x2-x1 LT 0: velnocen=ABS(velnocen)+180.*!DtoR
              else:velnocen=!values.f_nan
           endcase
           arr=[velpamax,velpamin,velnocen]
           tmp=WHERE(FINITE(arr))
           IF tmp [0] NE -1 then velpam=MEAN(arr[tmp])*!RADEG else begin
              velpa=[!values.f_nan,!values.f_nan]
              goto,endthis
           ENDELSE
           paerr=STDDEV(arr[tmp]*!RADEG)*3.
           IF paerr GT 10. then begin
              velpa=[!values.f_nan,!values.f_nan]
              goto,endthis
           ENDIF ELSE BEGIN
              velpa=[velpam,paerr]
              goto,endthis
           ENDELSE
        ENDELSE
     ENDELSE
  ENDELSE
  arr=[velpamax,velpamin,velnocen]
  maxang=MAX(arr*!RADEG,min=minang)
  if keyword_set(debug) then begin
     print,'These are  max angle and the min angle'
     print,maxang,minang
  ENDIF
  IF maxang GT 270 AND minang LT 90 then begin
     for i=0,n_elements(arr)-1 do begin
        if arr[i]*!RADEG GT 270 then arr[i]=arr[i]-2.*!pi
     endfor
  endif
  if keyword_set(debug) then begin
     print,'These are the angles to the max, the min and max to min'
     print,arr*!RADEG
  ENDIF
  tmp=WHERE(FINITE(arr))
  IF tmp [0] NE -1 then velpam=MEAN(arr[tmp])*!RADEG else begin
     velpa=[!values.f_nan,!values.f_nan]
     goto,endthis
  ENDELSE
  IF n_elements(tmp) GT 1 then paerr=STDDEV(arr[tmp]*!RADEG)*3. ELSE paerr=20.
  if keyword_set(debug) then begin
     print,'This is the final angle and error'
     print,velpam,paerr
  ENDIF
  
  velpa=[velpam,paerr]

endthis:
end
