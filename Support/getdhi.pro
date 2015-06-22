Pro getdhi,momentmap,header,PA,center,DHI

;+
; NAME:
;       GETDHI
;
; PURPOSE:
;       Routine to obtain the diameter in of the moment 0 map at the 1e20 columndensity
;
; CATEGORY:
;       Support
; 
; CALLING SEQUENCE:
;       GETDHI,momentmap,header,PA,center,DHI
;
;
; INPUTS:
;       momentmap = 2D array with the values of the momentmap pixels
;       header = fits header of the moment map
;       PA = the position angle of the galaxy in degree
;       center = the center of the galaxies in degrees or [hh:mm:ss,km/s]
;
; OPTIONAL INPUTS:
;       - 
;
; KEYWORD PARAMETERS:
;       -
;
; OUTPUTS:
;       DHI =  the diameter of the HI disk at 1e20 columndensity in arcsec
;
; OPTIONAL OUTPUTS:
;       -
; 
; PROCEDURES CALLED:
;       COLUMNDENSITY, INT_PROFILEV2, INTERPOLATE, STRTRIM() STRCOMPRESS(),
;       STR_SEP(), SXPAR(), STRUPCASE() 
;
; EXAMPLE:
;      
;
; MODIFICATION HISTORY:
;       Written by P.Kamphuis 01-01-2015 
;
; NOTE:
;     
;-

  COMPILE_OPT IDL2 
  tmp=str_sep(strtrim(strcompress(string(center[0])),2),':')
  IF n_elements(tmp) GT 1 then begin
     x=string(center[0])
     y=string(center[1])
     convertRADEC,x,y,/INVERT
     center[0:1]=[x,y]
     center[2]=double(center[2])
  ENDIF
  center[0]=sxpar(header,'CRPIX1')+((center[0]-sxpar(header,'CRVAL1'))/sxpar(header,'CDELT1'))
  center[1]=sxpar(header,'CRPIX2')+((center[1]-sxpar(header,'CRVAL2'))/sxpar(header,'CDELT2'))
  yrange=[-0.5*sxpar(header,'BMAJ'),0.5*sxpar(header,'BMAJ')]
  test=1
  int_profilev2,momentmap,majprofile,header=header,xcenter=center[0],ycenter=center[1],pa=pa,range=yrange,axis=xaxis,rotimage=test
  clevels=1E20
  Columndensity,clevels,center[2],[sxpar(header,'BMAJ')*3600.,sxpar(header,'BMIN')*3600.],/NCOLUMN
  clevels=clevels/1000.

                                ;First we need to cut the major axis
                                ;profile so that it does not include
                                ;noise
                                ;From the max if it drops below 0 we
                                ;assume that is the end of the galaxy
                                ;This could go wrong with two galaxies
                                ;in the cube

  maxprof=MAX(majprofile)
  maxisat=WHERE(maxprof EQ majprofile)
                                ;determine lower bound
  stop=0
  lowerbound=maxisat[0]
  WHILE stop EQ 0 and lowerbound GT 0 do begin
     lowerbound=lowerbound-1
     IF majprofile[lowerbound] LE 0. then stop=1.
  ENDWHILE
  stop=0
  upperbound=maxisat[0]
  WHILE stop EQ 0 and upperbound LT n_elements(majprofile)-2 do begin
     upperbound=upperbound+1
     IF majprofile[upperbound] LE 0. then stop=1.
  ENDWHILE
  majprofile[0:lowerbound]=0.
  majprofile[upperbound:n_elements(majprofile)-1]=0.
  tmp=WHERE(majprofile GE clevels)
  lowrad=1.
  IF tmp[0] GT -1 then begin
     interpolate,[xaxis[tmp[0]],xaxis[tmp[0]-1]],[majprofile[tmp[0]],majprofile[tmp[0]-1]],newradii=clevels,output=lowrad  
     highrad=1.
     interpolate,[xaxis[tmp[n_elements(tmp)-1]],xaxis[tmp[n_elements(tmp)-1]+1]],[majprofile[tmp[n_elements(tmp)-1]],majprofile[tmp[n_elements(tmp)-1]+1]],newradii=clevels,output=highrad
     tmp=str_sep(strtrim(strcompress(sxpar(header,'CTYPE1')),2),'---')
     ctype1=tmp[0]
     IF strupcase(ctype1) EQ 'RA' then DHI=(ABS(highrad-lowrad))*cos(sxpar(header,'CRVAL2')*!pi/180.)*3600. else $
        DHI=ABS(highrad-lowrad)
  ENDIF ELSE DHI=0. 
end



