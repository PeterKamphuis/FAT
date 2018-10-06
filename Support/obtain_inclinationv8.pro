Pro obtain_inclinationv8,map,inPA,inclination,center,EXTEND=extend,NOISE=noise,BEAM=beam,DEBUG=debug,gdlidl=gdlidl

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
;       09-05-2017 P.Kamphuis; If no values could be found the code
;                              would add random values to the errors.
;                              This introduced a randomness on the
;                              error leading to wildly varying cutoff
;                              values. Have changed this to set the
;                              errors in this case with a wide spread
;                              thus increasing the final error towards
;                              90.  
;       02-05-2017 P.Kamphuis; There was a bug where the profile was
;                              determined every tenth of a pixel but
;                              the beam correction was still the beam
;                              in normal pixels. Hence the estimates
;                              for small galaxies were far off. Additionally,
;                              as the beam smearing leads to lower
;                              inclinations already we do not apply the
;                              -2 correction to small galaxies.      
;       28-04-2017 P.Kamphuis; In case of a failed fit we now take a
;                              inclination from the ratio of the shape
;                              of the moment 0 map.  
;       02-06-2016 P.Kamphuis; Added GDL compatibility by replacing
;                              GAUSSFIT with MPFITFUN   
;       18-02-2016 P.Kamphuis; Replaced sigma with STDDEV   
;       Written 01-01-2015 P.Kamphuis v1.0
;
; NOTE:
;       please note that at low inclinations e.g < 10 deg there is a
;       random element to the determination due to the axis width
;       falling into the same pixels.
;     
;-

IF keyword_set(debug) then begin
   print,'OBTAIN_INCLINATIONV8: The PA'
   print,inPA
   print,'OBTAIN_INCLINATIONV8: The Noise'
   print,noise
ENDIF
IF n_elements(gdlidl) EQ 0 then gdlidl=0
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
beam=beam*10.
for i=0,n_elements(PA)-1 do begin
                                ;first get the profiles from the image
   int_profilev2,map,xprofile,pa=PA[i],xcenter=center[0],ycenter=center[1]
   int_profilev2,map,yprofile,pa=PA[i],xcenter=center[0],ycenter=center[1],/minor

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
   IF keyword_set(debug) then begin
      print,'OBTAIN_INCLINATIONV8: Max X Profile'
      print,maxxprof
      print,'OBTAIN_INCLINATIONV8: Max. Y Profile'
      print,maxyprof
   ENDIF
 
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
 
   IF tmp[n_elements(tmp)-1]-tmp[0] LT beam[1] then xFWHM=beam[1] else $
      xFWHM=SQRT(double(tmp[n_elements(tmp)-1]-tmp[0])^2-2.*beam[1]^2)
   tmpwidth[i]=xFWHM*2.5/(2.*SQRT(2*ALOG(2)))
   limit=maxyprof/3.
   IF keyword_set(debug) then begin
      print,'OBTAIN_INCLINATIONV8: X Width, Beam and Beam corrected'
      print,double(tmp[n_elements(tmp)-1]-tmp[0]),beam[1],xFWHM
   ENDIF
   IF limit LT 2*noise then limit=2.*noise 
   IF limit GT maxyprof then begin
      inclin3=0.
      inclin2=0.
      inclin15=0.
      goto,skip
   ENDIF
   tmp=WHERE(yprofile GT limit)
   IF tmp[n_elements(tmp)-1]-tmp[0] LT beam[1] then yFWHM=beam[1] else $
      yFWHM=SQRT(double(tmp[n_elements(tmp)-1]-tmp[0])^2-2.*beam[1]^2)
   IF keyword_set(debug) then begin
      print,'OBTAIN_INCLINATIONV8: Y Width, Beam and Beam corrected'
      print,double(tmp[n_elements(tmp)-1]-tmp[0]),beam[1],yFWHM
   ENDIF
   ratio=double(yFWHM/xFWHM)
   WHILE yFWHM EQ xFWHM do begin
      yFWHM=yFWHM+(randomu(seed)-0.5)
      xFWHM=xFWHM+(randomu(seed)-0.5)
      ratio=double(yFWHM/xFWHM)
   ENDWHILE
   IF keyword_set(debug) then begin
      print,'OBTAIN_INCLINATIONV8: Modified random'
      print,double(tmp[n_elements(tmp)-1]-tmp[0]),beam[1],yFWHM
   ENDIF
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
      print,'OBTAIN_INCLINATIONV8: Ratio final 3 is',ratio,inclin3,xFWHM,yFWHM
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
   IF tmp[n_elements(tmp)-1]-tmp[0] LT beam[1] then xFWHM=beam[1] else $
      xFWHM=SQRT(double(tmp[n_elements(tmp)-1]-tmp[0])^2-2*beam[1]^2)
   tmpwidth[i]=(tmpwidth[i]+xFWHM*2.5/(2.*SQRT(2*ALOG(2.))))/2.
   limit=maxyprof/2.
   IF limit LT 3*noise then limit=3.*noise 
   IF limit GT maxyprof then begin
      inclin2=0.
      inclin15=0.
      goto,skip
   ENDIF
   tmp=WHERE(yprofile GT limit)
   IF tmp[n_elements(tmp)-1]-tmp[0] LT beam[1] then yFWHM=beam[1] else $
      yFWHM=SQRT(double(tmp[n_elements(tmp)-1]-tmp[0])^2-2*beam[1]^2)
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
      print,'OBTAIN_INCLINATIONV8: Ratio final 2 is',ratio,inclin2,xFWHM,yFWHM
   ENDIF
   limit=maxxprof/1.5
   IF limit LT 4*noise then limit=4.*noise 
   IF limit GT maxxprof then begin
      inclin15=0.
      goto,skip
   ENDIF
   tmp=WHERE(xprofile GT limit)
   IF tmp[n_elements(tmp)-1]-tmp[0] LT beam[1] then xFWHM=beam[1] else $
      xFWHM=SQRT(double(tmp[n_elements(tmp)-1]-tmp[0])^2-2*beam[1]^2)
   tmpwidth[i]=(tmpwidth[i]+xFWHM*2.5/(2.*SQRT(2.*ALOG(2.))))/2.
   limit=maxyprof/1.5
   IF limit LT 4*noise then limit=4.*noise 
   IF limit GT maxyprof then begin
      inclin15=0.
      goto,skip
   ENDIF
   tmp=WHERE(yprofile GT limit)
   IF tmp[n_elements(tmp)-1]-tmp[0] LT beam[1] then yFWHM=beam[1] else $
      yFWHM=SQRT(double(tmp[n_elements(tmp)-1]-tmp[0])^2-2*beam[1]^2)
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
      print,'OBTAIN_INCLINATIONV8: Ratio final 1 is',ratio,inclin15,xFWHM,yFWHM
   ENDIF
   skip:
   incarr=[inclin15,inclin2,inclin3]
   IF keyword_set(debug) then begin
      print,'OBTAIN_INCLINATIONV8: This is incarr',incarr
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
;      tmpincerr[0:2,i]=!values.f_nan
      tmpincerr[0:2,i]=[85.,45,20.]
   ENDELSE

endfor  
IF keyword_set(debug) then begin
   print,'OBTAIN_INCLINATIONV8: This is tmpinc',tmpinc
ENDIF
inclination=dblarr(2)
nozer=WHERE(tmpinc GT 1.)
If nozer[0] EQ -1 then begin
   xstart=-1
   xend=-1
   ystart=-1
   yend=-1
   newmap=dblarr(n_elements(map[*,0]),n_elements(map[0,*]))
   newmap[*,*]=ROT(map[*,*],inpa[0]-90.,1.0,center[0],center[1],missing=!values.f_nan,cubic=-1,/PIVOT) ;For HI
   for i=0,n_elements(newmap[*,0])-1 do begin
      IF xstart EQ -1  then begin
         tmp=Where(newmap[i,*] NE 0 and FINITE(newmap[i,*]) NE 0)
         IF n_elements(tmp) GT 5 then xstart=i
      ENDIF
      IF xend EQ -1  then begin
         tmp=Where(newmap[n_elements(newmap[*,0])-1-i,*] NE 0 and FINITE(newmap[n_elements(newmap[*,0])-1-i,*]) NE 0)
         IF n_elements(tmp) GT 5 then xend=n_elements(newmap[*,0])-1-i
      ENDIF
      IF xstart NE -1 and xend NE -1 then break
   endfor
   for i=0,n_elements(newmap[0,*])-1 do begin
      IF ystart EQ -1  then begin
         tmp=Where(newmap[*,i] NE 0 and FINITE(newmap[*,i]) NE 0)
         IF n_elements(tmp) GT 5 then ystart=i
      ENDIF
      IF yend EQ -1  then begin
         tmp=Where(newmap[*,n_elements(newmap[*,0])-1-i] NE 0 and FINITE(newmap[*,n_elements(newmap[*,0])-1-i]) NE 0)
         IF n_elements(tmp) GT 5 then yend=n_elements(newmap[0,*])-1-i
      ENDIF
      IF ystart NE -1 and yend NE -1 then break
   endfor
    IF keyword_set(debug) then begin
      print,'OBTAIN_INCLINATIONV8: We failed but found the maps size as ',xstart,xend,ystart,yend
   ENDIF
   IF yend EQ -1 OR ystart EQ -1 or xend EQ -1 or xstart EQ -1 or xend EQ xstart then begin
      inclination[0]=20.
      inclination[1]=70.
   ENDIF ELSE BEGIN
      ratio=(double(yend-ystart))/(double(xend-xstart))
     
      IF keyword_set(debug) then begin
         print,'OBTAIN_INCLINATIONV8: We failed but found a ratio ',ratio
      ENDIF
      IF ratio GT 0.99999  then begin
                                ;at very low inclination x and y
                                ;easily get confused so if so we need
                                ;to flip as 0.2 corresponds to 40 deg
                                ;if it higher than that we leave the
                                ;measurement out
         IF ratio LT 1.2 then begin
            ratio=2-ratio
            inclination[0]=double(acos(SQRT((ratio^2-0.2^2)/0.96))*180./!pi+2.)
            inclination[1]=30.
         ENDIF ELSE BEGIN
            inclination[0]=20.
            inclination[1]=70.
            ;Else the position angle is very of 
          
         ENDELSE
      ENDIF ELSE BEGIN
         IF ratio LT 0.204 then  ratio=0.204  
         inclination[0]=double(acos(SQRT((ratio^2-0.2^2)/0.96))*180./!pi+2.)
         inclination[1]=30.
      ENDELSE
   ENDELSE
endif else begin
   IF keyword_set(debug) then begin
      print,'OBTAIN_INCLINATIONV8: tmpin[nozer]',tmpinc[nozer]
   ENDIF
   inclination[0]=TOTAL(tmpinc[nozer])/(n_elements(nozer))
   inclination[1]=STDDEV(tmpincerr[WHERE(FINITE(tmpincerr) EQ 1)])
   ;*(n_elements(tmpincerr)-TOTAL(FINITE(tmpincerr)))/3.
   IF keyword_set(debug) then begin
      print,'OBTAIN_INCLINATIONV8:',tmpincerr,'this is tmpincer'
      print,'OBTAIN_INCLINATIONV8:',inclination[0],inclination[1],MEDIAN(tmpincerr),MEAN(tmpincerr)
   ENDIF
endelse
;let's not do inclinations lower than 1
IF inclination[1] LT 2 AND  inclination [0] GE 20 then inclination[1]=2. 
IF inclination[1] LT 2 AND  inclination [0] LT 20 then inclination[1]=5.
extend=TOTAL(tmpwidth)/n_elements(tmpwidth)/10.
IF keyword_set(debug) then begin
   print,'OBTAIN_INCLINATIONV8: EXTENT',extend,3*beam[0]/10.
ENDIF

If inclination[0] LT 40. AND inclination[1] LT abs(40.-inclination[0])/2. then inclination[1]=abs(40.-inclination[0])/2.
If inclination[0] LT 5. then inclination[0]=5.
IF inclination[0] LT 80 AND inclination [0] GT 20 AND extend GT 3.*beam[0]/10. then inclination[0]=inclination[0]-2
end


