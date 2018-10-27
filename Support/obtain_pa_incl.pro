Pro obtain_pa_incl,map,PA,incl,center,EXTENT=extent,NOISE=noise,BEAM=beam,DEBUG=debug,gdlidl=gdlidl

;+
; NAME:
;       OBTAIN_PA_INCL
;
; PURPOSE:
;       Program to get the Pa and inclination of a galaxy based
;       on it's axis ratios in a full circle and their minimum
;       and maximums at FWHM FW1.5M FW2.5M the median of that. 
;
; CATEGORY:
;       Support
; 
; CALLING SEQUENCE:
;       OBTAIN_PA_INCL,map,inPA,inclination,center,extend=extend,noise=noise,beam=beam,/DEBUG
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
;       05-10-2018 P.Kamphuis; Starts from obtain_inclinationv8  
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


  IF n_elements(gdlidl) EQ 0 then gdlidl=0
  IF n_elements(beam) EQ 0 then beam=1
  IF n_elements(beam) EQ 1 then beam=[beam,beam] 

  corrected_map = 0.
  space2=1.
use_corrected:
  xpos=dblarr(n_elements(map[0,*]))
  ypos=dblarr(n_elements(map[0,*]))
  mapmax=MAX(map)
  tmp=WHERE(map GT mapmax/2.)
  maphigh=dblarr(n_elements(map[*,0]),n_elements(map[0,*]))
  maphigh[tmp]=map[tmp]
;determine the PAs
  space=5
  angles=findgen(fix(180/space))*space
  ratios = obtain_ratios(angles,map,center=center,gdlidl=gdlidl,noise=noise,beam=beam)


  tmp=WHERE(ratios GT 0.)
  maxrat=MAX(ratios[tmp],min=minrat)

  indexmax=WHERE(ratios EQ maxrat)
  indexmin=where(ratios EQ minrat)
  IF keyword_set(debug) then begin
     print,'indexmin'
     print,indexmin
  ENDIF
  pamax=angles[indexmax]-90
  tmp=WHERE(pamax LT 0.)
  IF tmp[0] NE -1 then pamax[tmp]=pamax[tmp]+180.
  IF keyword_set(debug) then print,pamax,angles[indexmin]
  pa=MEAN([pamax,angles[indexmin]])
  IF keyword_set(debug) then begin
     print,'Initial PA'
     print,pa
  ENDIF

wrongpa:


  angles=[findgen(fix(4*space/space2))*space2+pa-2.*space,findgen(fix(4*space/space2))*space2+pa-2*space+90]


  IF keyword_set(debug) then begin
     print,'We check these angles'
     print,angles
  ENDIF
  ratios=obtain_ratios(angles,map,center=center,MAJ_AXIS=tmpwidth,gdlidl=gdlidl,noise=noise,beam=beam)

  tmp=WHERE(ratios GT 0.)
  maxrat=MAX(ratios[tmp],min=minrat)
  indexmax=WHERE(ratios EQ maxrat)
  indexmin=where(ratios EQ minrat)

  IF keyword_set(debug) then print,ratios
  if indexmin[0] LT 2 then begin
     IF keyword_set(debug) then begin
        print,indexmin[0]
        print,pa
     ENDIF
     pa=pa[0]-space
     IF keyword_set(debug) then begin
        print,'weird'
        print,pa
     ENDIF
     goto, wrongpa
  ENDIF
  if n_elements(indexmin) GT 1 then indexmin=indexmin[fix(n_elements(indexmin)/2.)]
  pa=dblarr(2)
  incl=dblarr(2)

  pa[0]=MEAN([angles[indexmax]-90,angles[indexmin]])
  tmp=WHERE(ratios LE minrat AND ratios NE 0.)
  IF N_elements(tmp) GT 3 then begin
     pa[1]=ABS(angles[tmp[n_elements(tmp)-1]]-angles[tmp[0]])/2.
  ENDIF ELSE BEGIN
     pa[1]=ABS(angles[indexmin-1]-angles[indexmin+1])/2.
  ENDELSE
  IF keyword_set(debug) then print,'This is the diff max min,'+string(maxrat-minrat )
  IF keyword_set(debug) then begin
     print,'indexmin 2'
     print,indexmin
     print,angles[indexmax]-90,angles[indexmin]
     for i=0,n_elements(ratios)-1 do begin
        print, angles[i],ratios[i]
     endfor
  endif

  IF keyword_set(debug) then begin
     print,minrat,1./maxrat
     plot,angles,ratios,charsize=2.,charthick=1
  endif
  IF minrat LT 0.204 then  minrat=0.204
  IF 1./maxrat LT 0.204 then  maxrat=1./0.204
  
  incl[0]=MEAN([double(acos(SQRT((minrat^2-0.2^2)/0.96))*180./!pi+2.),double(acos(SQRT(((1./maxrat)^2-0.2^2)/0.96))*180./!pi+2.)])
  IF keyword_set(debug) then begin
     print,'hjkg understyand no?>'
     print,double(acos(SQRT((minrat^2-0.2^2)/0.96))*180./!pi+2.),double(acos(SQRT(((1./maxrat)^2-0.2^2)/0.96))*180./!pi+2.)
  endif
  IF N_elements(tmp) GT 3 then begin
     ratioerr=ratios[tmp[n_elements(tmp)-1]]
     IF ratioerr LT 0.225 then ratioerr=0.225
     errlow=double(acos(SQRT((ratioerr^2-0.2^2)/0.96))*180./!pi+2.)
     ratioerr=ratios[tmp[0]]
     IF ratioerr LT 0.225 then ratioerr=0.225
     errhigh=double(acos(SQRT((ratioerr^2-0.2^2)/0.96))*180./!pi+2.)
  ENDIF ELSE BEGIN
     ratioerr=ratios[indexmin-1]
     IF ratioerr LT 0.225 then ratioerr=0.225
     errlow=double(acos(SQRT((ratioerr^2-0.2^2)/0.96))*180./!pi+2.)
     ratioerr=ratios[indexmin+1]
     IF ratioerr LT 0.225 then ratioerr=0.225
     errhigh=double(acos(SQRT((ratioerr^2-0.2^2)/0.96))*180./!pi+2.)
  ENDELSE
  incl[1]=(abs(errlow-incl[0])+abs(errhigh-incl[0]))/2.
;ENDELSE

;If the inclination is low we become rather sensitive to
;inhomogeneities so we should try to correct for them
;let's not do inclinations lower than 1
  IF incl[1] LT 2 AND  incl[0] GE 20 then incl[1]=2. 
  IF incl[1] LT 2 AND  incl[0] LT 20 then incl[1]=5.
  extent=TOTAL(tmpwidth[indexmin-1:indexmin+1])/30.
  If incl[0] LT 40. AND incl[1] LT abs(40.-incl[0])/2. then incl[1]=abs(40.-incl[0])/2.
  If incl[0] LT 5. then incl[0]=5.
  IF incl[0] LT 80 AND incl[0] GT 20 AND extent GT 3.*beam[0]/10. then incl[0]=incl[0]-2
;If the inclination is low we become rather sensitive to
;inhomogeneities so we should try to correct for them

  If incl[0] LT 50 and corrected_map LT 1 then begin
                                ;we assume that the initial values are
                                ;accurate enough to use them to
                                ;deproject and rotate the map
     
     result=ROT(map, pa[0]-90,1.0,center[0],center[1],CUBIC=-0.5,/PIVOT)
     IF keyword_set(debug) then begin
        mkhdr,hed,result
        writefits,'or_map.fits',map,hed
        writefits,'rot_map.fits',result,hed
     ENDIF
                                ;we then need do deproject the yprofiles
                                ;our original axis =
     axis=findgen(n_elements(result[0,*]))-center[1]
     newaxis=axis*cos(incl[0]*!DtoR)
     dep_map=result
     dep_map[*]=0.
     rot_map=dep_map
     for i=0,n_elements(result[*,0])-1 do begin
        prof=dblarr(n_elements(result[0,*]))
        prof[*]=result[i,*]
        newprof=1
        interpolate,prof,axis,newradii=newaxis,output=newprof
        dep_map[i,*]=newprof     
     endfor
     

     
                                ;Now let's rotate our extract
                                ;line profile s at all angles
                                ;and create a median profile deprojected map and average out any inhomgeneities

     
     
     writefits,'dep_map.fits',dep_map,hed

     angles=(findgen(fix(355./5.)))*5.
     
     maps=dblarr(n_elements(angles)+1,n_elements(dep_map[*,0]),n_elements(dep_map[0,*]))
     maps[0,*,*]=dep_map[*,*]
     for i=0,n_elements(angles)-1 do begin
        result=ROT(dep_map, angles[i],1.0,center[0],center[1],CUBIC=-0.5,/PIVOT)
        maps[i,*,*]=result[*,*]
     endfor
     
     rot_map=dblarr(n_elements(dep_map[*,0]),n_elements(dep_map[0,*]))
     for i=0,n_elements(rot_map[0,*])-1 do begin
        for j=0,n_elements(rot_map[*,0])-1 do begin
           tmp=WHERE(maps[*,i,j] GT noise)
           rot_map[i,j]=MIN(maps[tmp,i,j])
        endfor
     endfor
                                ;rot_map=fat_smooth(rot_map,beam[0],/GAUSSIAN)
     writefits,'rot_map.fits',rot_map,hed
                                ;Let's subtract this minimum
                                ;map from the dep map
     
     res_map=dep_map-rot_map
     
     tmp=WHERE(res_map LT 2.*MEDIAN(res_map[WHERE(res_map GT 0.)]))
     res_map[tmp]=0.
     res_map=fat_smooth(res_map,3,/GAUSSIAN)
     
     writefits,'res_map.fits',res_map,hed

     
                                ;then we deproject the circular map
     inc_map=rot_map
     inc_map[*]=0.
     newaxis=axis/cos(incl[0]*!DtoR)
     for i=0,n_elements(inc_map[*,0])-1 do begin
        prof=dblarr(n_elements(rot_map[0,*]))
        prof[*]=res_map[i,*]
        newprof=1
        interpolate,prof,axis,newradii=newaxis,output=newprof
        inc_map[i,*]=newprof     
     endfor
     
     writefits,'inc_map.fits',inc_map,hed
                                ;Let's get a correction factor
                                ;by measuring the inclination
                                ;of the features

     



     
                                ;and we rotate it back to the estimated pa
     
     in_map=ROT(inc_map, -1*pa[0]-90,1.0,center[0],center[1],CUBIC=-0.5,/PIVOT)
     writefits,'corr_res_map.fits',in_map,hed
                                ;We need to check whether we have
                                ;significan inhomogeneities
     tmp=WHERE(in_map NE 0.)
     tmpinmap=WHERE(map NE 0.)
     if keyword_set(debug) then begin
        print, MEAN(in_map[tmp]),2.*noise,n_elements(tmp),n_elements(tmpinmap)*0.1
     ENDIF
     IF MEAN(in_map[tmp]) GT noise and n_elements(tmp) GT n_elements(tmpinmap)*0.075 then begin

     
                                ;Let's get the inclination and
                                ;pa from the inhomogeneities.
        ratios=obtain_ratios(angles,in_map,center=center,gdlidl=gdlidl,noise=MAX(in_map)/20.,beam=beam)
        tmp=WHERE(ratios GT 0.)
        maxrat=MAX(ratios[tmp],min=minrat)
        indexmax=WHERE(ratios EQ maxrat)
        indexmin=where(ratios EQ minrat)
        pa_in=MEAN([angles[indexmax]-90,angles[indexmin]])
     
        incl_in=MEAN([double(acos(SQRT((minrat^2-0.2^2)/0.96))*180./!pi+2.),double(acos(SQRT(((1./maxrat)^2-0.2^2)/0.96))*180./!pi+2.)])
        if keyword_set(debug) then begin
           print,'These are the pa and incl of the inhomogenieties'
           print,pa_in, incl_in
        endif
        tmp=WHERE(map NE 0.)
        bright_in=MEAN(in_map[tmp])
        map=map-in_map
        bright=MEAN(map[tmp])
                                ;map=fat_smooth(map,beam[0],/GAUSSIAN)
        writefits,'corr_map.fits',map,hed
     
        space2=space2/2.      
        corrected_map++
        goto,use_corrected
     ENDIF 
  endif
;If we looked for inhomogeneties we want to correct
  if corrected_map EQ 1 then begin
     IF keyword_set(debug) then begin
        print,'Before correcting for the inhomogeneities'
        print,pa,incl,pa_in,incl_in,bright,bright_in
     ENDIF
     incl[0]=incl[0]+(incl[0]-incl_in)*bright_in/bright*2.
     incl[1]=incl[1]+abs(incl[0]-incl_in)/2.
     IF abs(pa[0]-pa_in) GT 90 then pa_in=pa_in-180.
     pa[0]=pa[0]+(pa[0]-pa_in)*bright_in/bright*2.
     pa[1]=pa[1]+abs(pa[0]-pa_in)/2.
     IF keyword_set(debug) then begin
        print,pa,incl,pa_in,incl_in,bright,bright_in
     ENDIF
     
  endif

end






