FUNCTION obtain_ratios,angles,map,center=center,MAJ_AXIS=tmpwidth,gdlidl=gdlidl,noise=noise,beam=beam
  IF n_elements(beam) EQ 0. then beam = [1.,1]
  IF n_elements(gdlidl) EQ 0. then gdlidl=0.
  IF n_elements(center) EQ 0 then center=[fix(n_elements(map[*,0])/2.),fix(n_elements(map[0,*])/2.)]
  if n_elements(noise) EQ 0. then begin
     tmp=WHERE(map GT 0.)
     noise= 2.* MIN(map[tmp])
  ENDIF
  ratios=dblarr(n_elements(angles))
  tmpwidth=dblarr(n_elements(angles))
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
        
        IF tmp[n_elements(tmp)-1]-tmp[0] LT beam[1] then xFWHM=beam[1] else $
           xFWHM=SQRT(double(tmp[n_elements(tmp)-1]-tmp[0])^2-2.*beam[1]^2)
        
        tmp=WHERE(yprofile GT limit)
        IF tmp[n_elements(tmp)-1]-tmp[0] LT beam[1] then yFWHM=beam[1] else $
           yFWHM=SQRT(double(tmp[n_elements(tmp)-1]-tmp[0])^2-2.*beam[1]^2)
        ratios[i]=double(yFWHM/xFWHM)
     endif else ratios[i]=0.
     tmpwidth[i]=xFWHM*2.5/(2.*SQRT(2*ALOG(2)))/2.
     
     
  endfor
  return,ratios
end
