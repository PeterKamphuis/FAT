Pro obtain_inclinationv8,map,inPA,inclination,center,EXTEND=extend,NOISE=noise,BEAM=beam,DEBUG=debug

;+
; NAME:
;       OBTAIN_INCLINATIONV8
;
; PURPOSE:
;       Program to get the inclination of a galaxy based on it's axis ratio at FWHM FW1.5M FW2.5M the median of that. 
;
; CATEGORY:
;       Support
; 
; CALLING SEQUENCE:
;       OBTAIN_INCLINATIONV8,map,inPA,inclination,center,extend=extend,noise=noise,beam=beam,/DEBUG
;
;
; INPUTS:
;       map  = a moment 0 map of the galaxy for which to determine the inclination 
;       inPA = The PA. It is taken from north counter clockwise in
;       degrees. Either a 1D to 2D array. In case of the latter the
;       second value should be the error on the PA if omitted an error
;       of 2 deg is assumed
;       center = center of the galaxy in pixels
;
; OPTIONAL INPUTS:
;       - 
;
; KEYWORD PARAMETERS:
;       EXTEND = If set this will return the mean with of the profiles
;       in pixels (Yes it is misspelled)
;       NOISE = noise in the momentmap (This is required)
;       BEAM = Beam of the observations in pixels. If unset it is thought to be
;       1 pixel beam, if a single value a circular beam is assumed
;       /DEBUG - flag for getting additional output to screen
;
; OUTPUTS:
;       inclination = the determined inclination. A 2D array with
;       inclination and the error on the inclination
;
; OPTIONAL OUTPUTS:
;       -
; 
; PROCEDURES CALLED:
;       INT_PROFILEV2, INTERPOLATE, GAUSSFIT()
;
; EXAMPLE:
;      
;
; MODIFICATION HISTORY:
;       Written 01-01-2015 P.Kamphuis v1.0
;
; NOTE:
;       please note that at low inclinations e.g < 10 deg there is a
;       random element to the determination due to the axis width
;       falling into the same pixels.
;     
;-

IF keyword_set(debug) then begin
   print,'The PA'
   print,inPA
ENDIF
IF n_elements(beam) EQ 0 then beam=1
IF n_elements(beam) EQ 1 then beam=[beam,beam] 
xpos=dblarr(n_elements(map[0,*]))
ypos=dblarr(n_elements(map[0,*]))
mapmax=MAX(map)
tmp=WHERE(map GT mapmax/2.)
maphigh=dblarr(n_elements(map[*,0]),n_elements(map[0,*]))
maphigh[tmp]=map[tmp]
;determine the PAs
IF n_elements(inPA) EQ 2 then PA=[inPA[0]-inPA[1],inPA[0]-inPA[1]/2.,inPA[0],inPA[0]+inPA[1]/2.,inPA[0]+inPA[1]] $
                                  else  PA=[inPA[0]-2,inPA[0]-1,inPA[0],inPA[0]+1,inPA[0]+2]
tmpinc=dblarr(n_elements(PA))
tmpincerr=dblarr(3,n_elements(PA))
tmpwidth=dblarr(n_elements(PA))
for i=0,n_elements(PA)-1 do begin
                                ;first get the profiles from the image
   int_profilev2,map,xprofile,pa=PA[i],xcenter=center[0],ycenter=center[1]
   int_profilev2,map,yprofile,pa=PA[i],xcenter=center[0],ycenter=center[1],/minor

   ;Fit a Gaussian to reduce noise effects
   

                                ;get the maximum
   axis=findgen(n_elements(xprofile))
   xfit = GAUSSFIT(axis, xprofile, coeff, NTERMS=3)
   yfit = GAUSSFIT(axis, yprofile, coeff, NTERMS=3)
 
   newaxis=findgen(n_elements(xprofile)*10.)/10.
   xprof=1
   yprof=1
   interpolate,xfit,axis,newradii=newaxis,output=xprofile
   interpolate,yfit,axis,newradii=newaxis,output=yprofile 

 
   maxxprof=MAX(xprofile)
   maxyprof=MAX(yprofile)
                                ;and the limits
   limit=maxxprof/3.
   IF limit LT 2.*noise then limit=2.*noise 
   IF limit GT maxxprof then begin
      inclin3=0.
      inclin2=0.
      inclin15=0.
      goto,skip
   ENDIF
   tmp=WHERE(xprofile GT limit)
   IF tmp[n_elements(tmp)-1]-tmp[0] LT beam[0] then xFWHM=1. else $
      xFWHM=SQRT(double(tmp[n_elements(tmp)-1]-tmp[0])^2-beam[0]^2)
   tmpwidth[i]=xFWHM*2.5/(2.*SQRT(2*ALOG(2)))
   limit=maxyprof/3.
   IF limit LT 2*noise then limit=2.*noise 
   IF limit GT maxyprof then begin
      inclin3=0.
      inclin2=0.
      inclin15=0.
      goto,skip
   ENDIF
   tmp=WHERE(yprofile GT limit)
   IF tmp[n_elements(tmp)-1]-tmp[0] LT beam[0] then yFWHM=1. else $
      yFWHM=SQRT(double(tmp[n_elements(tmp)-1]-tmp[0])^2-beam[0]^2)
   ratio=double(yFWHM/xFWHM)
   WHILE yFWHM EQ xFWHM do begin
      yFWHM=yFWHM+(randomu(seed)-0.5)
      xFWHM=xFWHM+(randomu(seed)-0.5)
      ratio=double(yFWHM/xFWHM)
   ENDWHILE
   IF ratio GT 1.0 then begin
                                ;at very low inclination x and y
                                ;easily get confused so if so we need
                                ;to flip as 0.2 corresponds to 40 deg
                                ;if it higher than that we leave the
                                ;measurement out
      IF ratio LT 1.2 then begin
         ratio=2-ratio
         inclin3=double(acos(SQRT((ratio^2-0.2^2)/0.96))*180./!pi+2.)
      ENDIF ELSE BEGIN
         inclin3 = 0.
      ENDELSE 
   ENDIF else begin
      IF ratio LT 0.204 then ratio=0.204
      inclin3=double(acos(SQRT((ratio^2-0.2^2)/0.96))*180./!pi+2.)
   ENDELSE
   IF keyword_set(debug) then begin
      print,'Ratio final 3 is',ratio,inclin3,xFWHM,yFWHM
   ENDIF
                                ;new limit
   limit=maxxprof/2.
   IF limit LT 3.*noise then limit=3.*noise 
   IF limit GT maxxprof then begin
      inclin2=0.
      inclin15=0.
      goto,skip
   ENDIF
   tmp=WHERE(xprofile GT limit)
   IF tmp[n_elements(tmp)-1]-tmp[0] LT beam[0] then xFWHM=1. else $
      xFWHM=SQRT(double(tmp[n_elements(tmp)-1]-tmp[0])^2-beam[0]^2)
   tmpwidth[i]=(tmpwidth[i]+xFWHM*2.5/(2.*SQRT(2*ALOG(2.))))/2.
   limit=maxyprof/2.
   IF limit LT 3*noise then limit=3.*noise 
   IF limit GT maxyprof then begin
      inclin2=0.
      inclin15=0.
      goto,skip
   ENDIF
   tmp=WHERE(yprofile GT limit)
   IF tmp[n_elements(tmp)-1]-tmp[0] LT beam[1] then yFWHM=1. else $
      yFWHM=SQRT(double(tmp[n_elements(tmp)-1]-tmp[0])^2-beam[1]^2)
   ratio=double(yFWHM/xFWHM)
   WHILE yFWHM EQ xFWHM do begin
      yFWHM=yFWHM+(randomu(seed)-0.5)
      xFWHM=xFWHM+(randomu(seed)-0.5)
      ratio=double(yFWHM/xFWHM)
   ENDWHILE
   IF ratio GT 1.0 then begin
                                ;at very low inclination x and y
                                ;easily get confused so if so we need
                                ;to flip as 0.2 corresponds to 40 deg
                                ;if it higher than that we leave the
                                ;measurement out
      IF ratio LT 1.2 then begin
         ratio=2-ratio
         inclin2=double(acos(SQRT((ratio^2-0.2^2)/0.96))*180./!pi+2.)
      ENDIF ELSE BEGIN
         inclin2= 0.
      ENDELSE
    
   ENDIF else begin
      IF ratio LT 0.204 then ratio=0.204
      inclin2=double(acos(SQRT((ratio^2-0.2^2)/0.96))*180./!pi+2.)
   ENDELSE
   IF keyword_set(debug) then begin
      print,'Ratio final 2 is',ratio,inclin2,xFWHM,yFWHM
   ENDIF
   limit=maxxprof/1.5
   IF limit LT 4*noise then limit=4.*noise 
   IF limit GT maxxprof then begin
      inclin15=0.
      goto,skip
   ENDIF
   tmp=WHERE(xprofile GT limit)
   IF tmp[n_elements(tmp)-1]-tmp[0] LT beam[0] then xFWHM=1. else $
      xFWHM=SQRT(double(tmp[n_elements(tmp)-1]-tmp[0])^2-beam[0]^2)
   tmpwidth[i]=(tmpwidth[i]+xFWHM*2.5/(2.*SQRT(2.*ALOG(2.))))/2.
   limit=maxyprof/1.5
   IF limit LT 4*noise then limit=4.*noise 
   IF limit GT maxyprof then begin
      inclin15=0.
      goto,skip
   ENDIF
   tmp=WHERE(yprofile GT limit)
   IF tmp[n_elements(tmp)-1]-tmp[0] LT beam[1] then yFWHM=1. else $
      yFWHM=SQRT(double(tmp[n_elements(tmp)-1]-tmp[0])^2-beam[1]^2)
   ratio=double(yFWHM/xFWHM)
  
   WHILE yFWHM EQ xFWHM do begin
      yFWHM=yFWHM+(randomu(seed)-0.5)
      xFWHM=xFWHM+(randomu(seed)-0.5)
      ratio=double(yFWHM/xFWHM)
   ENDWHILE
   IF ratio GT 0.99999  then begin
                                ;at very low inclination x and y
                                ;easily get confused so if so we need
                                ;to flip as 0.2 corresponds to 40 deg
                                ;if it higher than that we leave the
                                ;measurement out
      IF ratio LT 1.2 then begin
         ratio=2-ratio
         inclin15=double(acos(SQRT((ratio^2-0.2^2)/0.96))*180./!pi+2.)
      ENDIF ELSE BEGIN
         inclin15 = 0.
      ENDELSE
   ENDIF else begin
                                ;IF ratio is smaller than 0.2 we get none coherent answers
      IF ratio LT 0.204 then  ratio=0.204 
      inclin15=double(acos(SQRT((ratio^2-0.2^2)/0.96))*180./!pi+2.)
   ENDELSE
   IF keyword_set(debug) then begin
      print,'Ratio final 1 is',ratio,inclin15,xFWHM,yFWHM
   ENDIF
   skip:
   incarr=[inclin15,inclin2,inclin3]
   IF keyword_set(debug) then begin
      print,'This is incarr',incarr
   ENDIF
   tmpzero=WHERE(incarr GT 0.1)
   IF tmpzero[0] NE -1 then begin
      incarrnozer=dblarr(n_elements(tmpzero))
      incarrnozer=incarr[tmpzero]
      tmp=incarrnozer[SORT(incarrnozer)]
      IF n_elements(tmpzero) GT 1 then begin
         tmpinc[i]=tmp[0]
      endif else begin
         tmpinc[i]=tmp[0]
      endelse
      tmpincerr[0:2,i]=incarr
   ENDIF ELSE BEGIN
      tmpinc[i]=0.
      tmpincerr[0:2,i]=[randomu(seed)*90.,randomu(seed)*90.,randomu(seed)*90.]
   ENDELSE

endfor
inclination=dblarr(2)
nozer=WHERE(tmpinc GT 1.)
If nozer[0] EQ -1 then begin
   inclination[0]=20.
   inclination[1]=40.
endif else begin
   inclination[0]=TOTAL(tmpinc[nozer])/n_elements(nozer)
   inclination[1]=SIGMA(tmpincerr[WHERE(FINITE(tmpincerr) EQ 1)])
   IF keyword_set(debug) then begin
      print,tmpincerr,'this is tmpincer'
      print,inclination[1],MEDIAN(tmpincerr),MEAN(tmpincerr)
   ENDIF
endelse
;let's not do inclinations lower than 1
IF inclination[1] LT 2 AND  inclination [0] GE 20 then inclination[1]=2. 
IF inclination[1] LT 2 AND  inclination [0] LT 20 then inclination[1]=5.
IF inclination[0] LT 80 AND inclination [0] GT 20 then inclination[0]=inclination[0]-2
extend=TOTAL(tmpwidth)/n_elements(tmpwidth)
end


