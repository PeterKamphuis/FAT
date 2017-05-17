Pro get_fixedringsv9,Parametersin,rings,smooth=smooth,debug=debug,warp_output=warp_output,radii=radii,SBR=SBR

;+
; NAME:
;       GET_FIXEDRINGSV9
;
; PURPOSE:
;       Little routine to determine how many rings of the parameters
;       should be considered as the flat inner part of the model based
;       on the tiltogram of the current variations
;
; CATEGORY:
;       Support
; 
; CALLING SEQUENCE:
;       GET_FIXEDRINGSV9,Parametersin,rings
;
;
; INPUTS:
;       Parametersin = The INCL and PA parameters
;
; OPTIONAL INPUTS:
;       - 
;
; KEYWORD PARAMETERS:
;       -
;
; OUTPUTS:
;       rings = the amount of rings that should remain fixed
;
; OPTIONAL OUTPUTS:
;       -
; 
; PROCEDURES CALLED:
;       ACOS(),ATAN(),MAX(),FINITE(),MAX() STDDEV(), ROBUST_SIGMA(), FLOOR()
;
; EXAMPLE:
;      
;
; MODIFICATION HISTORY:
;       06-05-2017 P.Kamphuis; Added the posibility to have the warp
;                              parameters written into a directory
;                              called Warp_Info   
;       18-02-2016 P.Kamphuis; Replaced sigma with STDDEV   
;       Written by P.Kamphuis 01-01-2015 
;
; NOTE:
;     
;-
  COMPILE_OPT IDL2
  Parameters=Parametersin
  IF keyword_set(warp_output) then begin
     angle=dblarr(3,2)
     tmp=WHERE(SBR[*,0] GT 1e-8)
     IF tmp[0] NE -1 then outer1=tmp[n_elements(tmp)-1] else outer1=n_elements(Tirresult[*,0])-1
     tmp=WHERE(SBR[*,1] GT 1e-8)
     IF tmp[0] NE -1 then outer2=tmp[n_elements(tmp)-1] else outer2=n_elements(Tirresult[*,0])-1
     max_index=[outer1,outer2] 
  ENDIF
  
                                ;First we need to calculate the normal
                                ;vectors for all rings on both sides

  ;[PA,PA_2,INCL,INCL_2]
  PA1=Parameters[*,0]
  PA2=Parameters[*,1]
  INCL1=Parameters[*,2]
  INCL2=Parameters[*,3]
  x1=dblarr(n_elements(PA1))
  y1=dblarr(n_elements(PA1))
  z1=dblarr(n_elements(PA1))
  x2=dblarr(n_elements(PA2))
  y2=dblarr(n_elements(PA2))
  z2=dblarr(n_elements(PA2))
  add1=dblarr(n_elements(PA1))
  add2=dblarr(n_elements(PA2))
  for i=0,n_elements(PA1)-1 do begin
     add1[i]=0.
     while PA1[i] GE 90 do begin
        PA1[i]=PA1[i]-90
        add1[i]=add1[i]+90.
     endwhile
     if PA1[i] EQ 0. then PA1[i]=0.00001
     Theta=ATAN(tan(INCL1[i]*(!pi/180.))*tan(PA1[i]*(!pi/180.)))
     phi=ATAN(tan(PA1[i]*(!pi/180.))/sin(theta))
                                ;So in cartesian coordinates
     x1[i]=sin(theta)*cos(phi)
     y1[i]=sin(theta)*sin(phi)
     z1[i]=cos(theta)
  endfor
  for i=0,n_elements(PA2)-1 do begin
     add2[i]=0.
     while PA2[i] GE 90 do begin
        PA2[i]=PA2[i]-90
        add2[i]=add2[i]+90.
     endwhile
     if PA2[i] EQ 0. then PA2[i]=0.00001
     Theta=ATAN(tan(INCL2[i]*(!pi/180.))*tan(PA2[i]*(!pi/180.)))
     phi=ATAN(tan(PA2[i]*(!pi/180.))/sin(theta))
                                ;So in cartesian coordinates
     x2[i]=sin(theta)*cos(phi)
     y2[i]=sin(theta)*sin(phi)
     z2[i]=cos(theta)
  endfor
                                ;Next we calculate the tiltogram
  tiltogram1=dblarr(n_elements(PA1),n_elements(PA1))
  tiltogram2=dblarr(n_elements(PA2),n_elements(PA2))
  orings=findgen(n_elements(PA1))+0.5
  rrings=[findgen(n_elements(PA1)*10.)/10.,11]
  for i=0,n_elements(PA1)-1 do begin
     for j=0,n_elements(PA1)-1 do begin
        tiltogram1[i,j]=ACOS(x1[i]*x1[j]+y1[i]*y1[j]+z1[i]*z1[j])/!DtoR
         IF NOT FINITE(tiltogram1[i,j]) then tiltogram1[i,j]=0.
     endfor
  endfor
  for i=0,n_elements(PA2)-1 do begin
     for j=0,n_elements(PA2)-1 do begin
        tiltogram2[i,j]=ACOS(x2[i]*x2[j]+y2[i]*y2[j]+z2[i]*z2[j])/!DtoR
        IF NOT FINITE(tiltogram2[i,j]) then tiltogram2[i,j]=0.
     endfor
  endfor
                                ;If we run in debug mode we write fits files for the tiltograms
  if keyword_set(warp_output) then begin
   ;  CD,workingdir,CURRENT=old_dir
     exist=FILE_TEST('Warp_Info',/DIRECTORY)
     if exist EQ 0 then spawn,'mkdir Warp_Info'
     
     CD,'Warp_Info'
     mkhdr,hed,tiltogram1
     writefits,'Receding_Tiltogram.fits',tiltogram1,hed
     mkhdr,hed,tiltogram2
     writefits,'Approaching_Tiltogram.fits',tiltogram2,hed
     CD,'../'
  endif  
  if keyword_set(debug) then begin
     mkhdr,hed,tiltogram1
     writefits,'test1.fits',tiltogram1,hed
     mkhdr,hed,tiltogram2
     writefits,'test2.fits',tiltogram2,hed
  ENDIF
  tmp=dblarr(n_elements(PA1),n_elements(PA1)*10+1)
  rtiltogram1=dblarr(n_elements(PA1)*10.+1,n_elements(PA1)*10+1)
  rtiltogram2=dblarr(n_elements(PA2)*10.+1,n_elements(PA2)*10.+1)
  for i=0,n_elements(PA1)-1 do begin
     val=dblarr(n_elements(tiltogram1[i,*]))
     val[*]=tiltogram1[i,*]
     output=1
     Interpolate,Val,orings,output=output,newradii=rrings
     tmp[i,*]=output[*]
  endfor
  for i=0,n_elements(tmp[0,*])-1 do begin
     val=tmp[*,i]
     Interpolate,Val,orings,output=output,newradii=rrings
     rtiltogram1[*,i]=output[*]
  endfor
  tmp=dblarr(n_elements(PA2),n_elements(PA2)*10+1)
  for i=0,n_elements(PA2)-1 do begin
     val=dblarr(n_elements(tiltogram2[i,*]))
     val[*]=tiltogram2[i,*]
     output=1
     Interpolate,Val,orings,output=output,newradii=rrings
     tmp[i,*]=output[*]
  endfor
  for i=0,n_elements(tmp[0,*])-1 do begin
     val=tmp[*,i]
     Interpolate,Val,orings,output=output,newradii=rrings
     rtiltogram2[*,i]=output[*]
  endfor
  IF keyword_set(debug) then begin
     mkhdr,hed,rtiltogram1
     writefits,'tests1.fits',rtiltogram1,hed
     mkhdr,hed,rtiltogram2
     writefits,'tests2.fits',rtiltogram2,hed
     help,rtiltogram1
  endif
                                ;now that we have the tiltograms we
                                ;need to calculate the average at each
                                ;radius. we want at 5 15 25 35 45
  thetain=dblarr(n_elements(tiltogram1[*,0]))
  thetamut=dblarr(n_elements(tiltogram1[*,0]))
  thetaout=dblarr(n_elements(tiltogram1[*,0]))
  for i=0,n_elements(tiltogram1[*,0])-1 do begin
     thetain[i]=MEAN(rtiltogram1[0:(i+1)*10,0:(i+1)*10])
     thetamut[i]=MEAN(rtiltogram1[(i+1)*10:n_elements(rtiltogram1[*,0])-1,0:(i+1)*10])
     thetaout[i]=MEAN(rtiltogram1[(i+1)*10:n_elements(rtiltogram1[*,0])-1,(i+1)*10:n_elements(rtiltogram1[*,0])-1])
  endfor

                                ;And then we want to apply the rules
                                ;(i) the difference between thetain
                                ;and thetamut is larger than the differences observed at other radii
                                ;(ii) thetain < 5 deg
                                ; (iii) thetamut > 15 deg
  diff=ABS(thetain[*]-thetamut[*])
  found=0
  mxdif=1
  IF  keyword_set(debug) then begin
     print,'This is the difference in disk1'
     print,diff
  ENDIF
  While found EQ 0 AND mxdif NE 0 DO BEGIN
     mxdif=MAX(diff)
     rings=WHERE(mxdif EQ diff)
;     rings=WHERE(diff GT 10)
;     IF keyword_set(debug) then print,rings
     ;We want the smallest ring that satifies our conditions
     IF n_elements(rings) GT 1 then rings=rings[0]
     IF thetain[rings] LT 5 AND thetamut[rings] GT 15 then begin
        found=1
     ENDIF ELSE diff[rings]=0
  ENDWHILE
  IF TOTAL(diff) EQ 0 then rings=n_elements(diff)-1
  rings1=rings
  IF keyword_set(debug) then begin
     print,'This is the inner theta'
     print,thetain
     print,'This is the edge theta'
     print,thetamut
     print,'This is the outer theta'
     print,thetaout
     print,'this is what we found in disk 1',rings1
  ENDIF
 
  if keyword_set(debug) then begin
     mkhdr,hed,tiltogram1
     writefits,'test1.fits',tiltogram1,hed
     mkhdr,hed,tiltogram2
     writefits,'test2.fits',tiltogram2,hed
  ENDIF
  if keyword_set(warp_output) then begin
     angle[*,0]=[tiltogram1[0,max_index[0]],tiltogram2[0,max_index[1]],(tiltogram1[0,max_index[0]]+tiltogram2[0,max_index[1]])/2.]
     tmp1=MAX(tiltogram1[0,0:max_index[0]])
     tmp2=MAX(tiltogram2[0,0:max_index[1]])
     tmp3=MAX((tiltogram1[0,0:max_index[0]]+tiltogram2[0,0:max_index[1]])/2.)
  
     angle[*,1]=[tmp1[0],tmp2[0],tmp3[0]]
  ENDIF
  
  tmp=dblarr(n_elements(PA1),n_elements(PA1)*10+1)
  rtiltogram1=dblarr(n_elements(PA1)*10.+1,n_elements(PA1)*10+1)
  rtiltogram2=dblarr(n_elements(PA2)*10.+1,n_elements(PA2)*10.+1)
  for i=0,n_elements(PA1)-1 do begin
     val=dblarr(n_elements(tiltogram1[i,*]))
     val[*]=tiltogram1[i,*]
     output=1
     Interpolate,Val,orings,output=output,newradii=rrings
     tmp[i,*]=output[*]
  endfor
  for i=0,n_elements(tmp[0,*])-1 do begin
     val=tmp[*,i]
     Interpolate,Val,orings,output=output,newradii=rrings
     rtiltogram1[*,i]=output[*]
  endfor
  tmp=dblarr(n_elements(PA2),n_elements(PA2)*10+1)
  for i=0,n_elements(PA2)-1 do begin
     val=dblarr(n_elements(tiltogram2[i,*]))
     val[*]=tiltogram2[i,*]
     output=1
     Interpolate,Val,orings,output=output,newradii=rrings
     tmp[i,*]=output[*]
  endfor
  for i=0,n_elements(tmp[0,*])-1 do begin
     val=tmp[*,i]
     Interpolate,Val,orings,output=output,newradii=rrings
     rtiltogram2[*,i]=output[*]
  endfor
  IF keyword_set(debug) then begin
     mkhdr,hed,rtiltogram1
     writefits,'tests1.fits',rtiltogram1,hed
     mkhdr,hed,rtiltogram2
     writefits,'tests2.fits',rtiltogram2,hed
     help,rtiltogram1
  endif
                                ;now that we have the tiltograms we
                                ;need to calculate the average at each
                                ;radius. we want at 5 15 25 35 45
  thetain1=dblarr(n_elements(tiltogram1[*,0]))
  thetamut1=dblarr(n_elements(tiltogram1[*,0]))
  thetaout1=dblarr(n_elements(tiltogram1[*,0]))
  for i=0,n_elements(tiltogram1[*,0])-1 do begin
     thetain1[i]=MEAN(rtiltogram1[0:(i+1)*10,0:(i+1)*10])
     thetamut1[i]=MEAN(rtiltogram1[(i+1)*10:n_elements(rtiltogram1[*,0])-1,0:(i+1)*10])
     thetaout1[i]=MEAN(rtiltogram1[(i+1)*10:n_elements(rtiltogram1[*,0])-1,(i+1)*10:n_elements(rtiltogram1[*,0])-1])
  endfor

                                ;And then we want to apply the rules
                                ;(i) the difference between thetain
                                ;and thetamut is larger than the differences observed at other radii
                                ;(ii) thetain < 5 deg
                                ; (iii) thetamut > 15 deg
  diff=ABS(thetain1[*]-thetamut1[*])
  found=0
  mxdif=1
  IF  keyword_set(debug) then begin
     print,'This is the difference in disk1'
     print,diff
  ENDIF
  While found EQ 0 AND mxdif NE 0 DO BEGIN
     mxdif=MAX(diff)
     rings=WHERE(mxdif EQ diff)
;     rings=WHERE(diff GT 10)
;     IF keyword_set(debug) then print,rings
     ;We want the smallest ring that satifies our conditions
     IF n_elements(rings) GT 1 then rings=rings[0]
     IF thetain1[rings] LT 5 AND thetamut1[rings] GT 15 then begin
        found=1
     ENDIF ELSE diff[rings]=0
  ENDWHILE
  IF TOTAL(diff) EQ 0 then rings=n_elements(diff)-1
  rings1=rings
  IF keyword_set(debug) then begin
     print,'This is the inner theta'
     print,thetain1
     print,'This is the edge theta'
     print,thetamut1
     print,'This is the outer theta'
     print,thetaout1
     print,'this is what we found in disk 1',rings1
  ENDIF
  if keyword_set(warp_output) then begin
     IF rings1 GT max_index[0] then rings1=max_index[0]
     CD,'Warp_Info',CURRENT=old_dir    
     openw,1,'Warp_Info.txt'
     printf,1,format='(6A20)','Side','Inner Theta','Outer Theta','Warp Radius','Max Angle','Outer Angle'
     printf,1,format='(A20,5F20.3)','Receding',thetain1[rings1],thetaout1[rings1],radii[rings1],angle[0,1], angle[0,0]
     close,1
     CD,old_dir
  ENDIF
                               ;now that we have the tiltograms we
                                ;need to calculate the average at each
                                ;radius. for the other side
  thetain=dblarr(n_elements(tiltogram2[*,0]))
  thetamut=dblarr(n_elements(tiltogram2[*,0]))
  thetaout=dblarr(n_elements(tiltogram2[*,0]))
  for i=0,n_elements(tiltogram2[*,0])-1 do begin
     thetain[i]=MEAN(rtiltogram2[0:(i+1)*10,0:(i+1)*10])
     thetamut[i]=MEAN(rtiltogram2[(i+1)*10:n_elements(rtiltogram2[*,0])-1,0:(i+1)*10])
     thetaout[i]=MEAN(rtiltogram2[(i+1)*10:n_elements(rtiltogram2[*,0])-1,(i+1)*10:n_elements(rtiltogram2[*,0])-1])
  endfor

                                ;And then we want to apply the rules
                                ;to the otherside as well
                                ;(i) the difference between thetain
                                ;and thetamut is larger than the differences observed at other radii
                                ;(ii) thetain < 5 deg
                                ; (iii) thetamut > 15 deg
  diff=ABS(thetain[*]-thetamut[*])
  found=0
  mxdif=1
   IF  keyword_set(debug) then begin
     print,'This is the difference in disk2'
     print,diff
  ENDIF
  While found EQ 0 AND mxdif NE 0 DO BEGIN
     mxdif=MAX(diff)
     rings=WHERE(mxdif EQ diff)
 ;    rings = WHERE(diff GT 10)
     IF keyword_set(debug) then print,rings
     IF n_elements(rings) GT 1 then rings=rings[0]
     IF thetain[rings] LT 5 AND thetamut[rings] GT 15 then begin
        found=1
     ENDIF ELSE diff[rings]=0
  ENDWHILE
  IF TOTAL(diff) EQ 0 then rings=n_elements(diff)-1
  rings2=rings

  IF rings1 LT rings2 then rings=rings1 else rings=rings2
  IF keyword_set(debug) then begin
     print,'This is the inner theta'
     print,thetain
     print,'This is the edge theta'
     print,thetamut
     print,'This is the outer theta'
     print,thetaout
     print,'this is what we found in disk 2',rings2
  ENDIF
  if keyword_set(warp_output) then begin
     IF rings2 GT max_index[1] then rings2=max_index[1]
     CD,'Warp_Info',CURRENT=old_dir    
     openu,1,'Warp_Info.txt',/append
     printf,1,format='(A20,5F20.3)','Approaching',thetain[rings2],thetaout[rings2],radii[rings2],angle[1,1], angle[1,0]
     printf,1,format='(A20,5F20.3)','Average',(thetain1[rings1]+thetain[rings2])/2.,(thetaout1[rings1]+thetaout[rings2])/2.,(radii[rings1]+radii[rings2])/2.,angle[2,1], angle[2,0]
     close,1
     CD,old_dir
  ENDIF
  if rings LT 3 then rings=3
 
  IF NOT keyword_set(smooth) then begin
     IF rings GT fix(n_elements(Parameters[*,0])/2.) then begin
        IF fix(n_elements(Parameters[*,0])/2.) GE 3 then rings=fix(n_elements(Parameters[*,0])/2.) else rings=3
     ENDIF
     IF n_elements(Parameters[*,0]) GT 12 then begin
        IF rings LT n_elements(Parameters[*,0])/4. then ring=ceil(n_elements(Parameters[*,0])/4.)
     ENDIF
  ENDIF ELSE BEGIN
     IF rings GT fix(n_elements(Parameters[*,0])/2.) then begin
        IF n_elements(Parameters[*,0]) GT 10 then begin
           correction=0.3
           WHILE rings-fix(correction*rings) LT fix(n_elements(Parameters[*,0])/2.) DO BEGIN
              correction=correction/1.1
           ENDWHILE
           rings=rings-fix(correction*rings)
        ENDIF ELSE rings=rings-1
     ENDIF
  ENDELSE
 
end
