Pro revised_regularisation_com,PAin,SBRin,RADIIin,error=errorin,fixedrings=fixedringsin,REVERSE=rev,NOWEIGHT=noweight,Difference=DDivin,cutoff=cutoffin,arctan=arctanin,debug=debug,order=order,max_deviation=maxdevin,max_par=PAmaxin,min_par=PAminin,accuracy=accuracy,extending=extending,gdlidl=gdlidl,sloped=slopedrings,log=log

;+
; NAME:
;       REVISED_REGULARISTION_COM
;
; PURPOSE:
;       Routine to regularize the parameters of a tirific fit
;
; CATEGORY:
;       Support
; 
; CALLING SEQUENCE:
;        REVISED_REGULARISTION, PAin, SBRin, RADIIin, ERROR=errorin,
;        FIXEDRINGS=fixedringsin, /REVERSE, NOWEIGHT=noweight,
;        DIFFERENCE=DDivin, CUTOFF=cutoffin, ARCTAN=arctanin, /DEBUG,
;        /NOCENTRAL, ORDER=order, MAX_DEVIATION=maxdevin, MAX_PAR =
;        PAmaxin, MIN_PAR = PAminin, ACCURACY = accuracy, /EXTENDING
;
;
; INPUTS:
;       PAin =  the array of parameters as function of radius
;       SBRin = The Surface brightness profile going with the
;       parameters
;       RADIIin = The location of the nodes in arcsec
;
; OPTIONAL INPUTS:
;        FIXEDRINGS = the amount of rings in the center that should
;        have the same value
;        DIFFERENCE = The basic error on the parameters. If variations
;        between nodes are less than this they are assumed to have the
;        same value.
;        CUTOFF = cutoff values for each node
;        ARCTAN = Type of regularisation requested. 0 = polynomial
;        fit, 1 is simple averaging, 2 = avreaging with a check that
;        the outer part does not decline strongly 
;        MAX_DEVIATION = The maximum a parameter is allowed to vary
;        from 1 node to the next
;        MAX_PAR = The maximum allowed value for the parameters
;        MIN_PAR =  The minimum allowed value for the parameters
;        ACCURACY = parameter to indicate the accuracy of the profile.
;        The smaller this number the closer the regularistion will
;        follow the input
;        gdlidl = 1 if running gdl 0 if idl  
;
; KEYWORD PARAMETERS:
;        /REVERSE - Set this keyword to reverse the parameter before
;                   regularistion. Do this when the larger radii are
;                   expected to be flat.
;        /NOWEIGHT - Set this keyword in order cancel the weighing
;                    based on the SBR profile and treat all point as
;                    equally accurate.
;        /DEBUG - Keyword to get debugging information.
;        /NOCENTRAL - Set this keyword to exclude the inner most point
;                     from the regularisation 
;        /EXTENDING - Keyword for indicating that we are use the
;                     regularisation to get an estimate for an
;                     additional outer ring
; OUTPUTS:
;       PAin = the regularised parameter
;
; OPTIONAL OUTPUTS:
;       ERROR = the errors on the regularised parameter
;       ORDER = the order of the applied polynomial;     
; 
; PROCEDURES CALLED:
;       MEAN(), FAT_FIT(), STDDEV(), ROBUST_SIGMA()
;
; EXAMPLE:
;      
;
; MODIFICATION HISTORY:
;       13-09-2018 P.Kamphuis; Added increased errors for outer rings
;                              when galaxies >25 rings. This increase
;                              is based on their deviation from the
;                              mean value.
;       08-03-2017 P.Kamphuis; Added the condition that when the
;                              initial order is 0 when refitting do
;                              not take the lowest order. Aditionally
;                              when minchi is too high check that the
;                              errors do not increase beyond 3* maxdev
;                              and removed the condition that when the
;                              profile approaches minpar the error is
;                              the distance to minpar as it actually
;                              reduces the error.
;       01-03-2017 P.Kamphuis; Added a final condition where the
;                              polynomial is re-fitted with the
;                              estimated errors from an initial fit.  
;       13-06-2016 P.Kamphuis; Increased the option to check for
;                              sloped warp parts and increase the error.  
;       30-03-2016 P.Kamphuis; A complete overhaul of this routine
;                              build from scratch based on
;                              parrameterreguv87. This new routine
;                              will evaluate the PA, INCL at the
;                              same time as unreliable rings will
;                              appear in both hence increase
;                              stability.
;                              As the evalution of the rotation curve
;                              is signifcantly different these will
;                              now be evaluated in a different routine  
;       08-03-2016 P.Kamphuis; Replaced poly_fit with FAT_FIT, this
;                              routine is now GDL compatible   
;       18-02-2016 P.Kamphuis; Replaced sigma with STDDEV   
;       11-08-2015 P.Kamphuis; Modified the handling of single value
;       rotation curves
;       Written 01-01-2015 P.Kamphuis v1.0
;
; NOTE:
;This routine has an extensive debugging mode that can be turned on
;with the keyword, /debug
;     
;-
  COMPILE_OPT IDL2 

; First thing is to check all the input whetehr it is reasonable
  SBR=SBRin
  RADII=RADIIin
  PA=PAin
 
  If n_elements(errorin) GT 0 then error=errorin
  If n_elements(fixedringsin) GT 0 then fixedrings=fixedringsin else fixedrings=0
  If n_elements(DDivin) GT 0 then DDiv=DDivin
  If n_elements(cutoffin) GT 0 then cutoff=cutoffin
                                ;arctan can have the values of 0: for a polynomial fit
                                ;1: for a smoothing of PA or INCL  (default)
                                ;2: for a smoothing of VROT
  If n_elements(arctanin) GT 0 then arctan=arctanin else arctan=1
  If n_elements(maxdevin) GT 0 then  maxdev=maxdevin
  If n_elements(PAmaxin) GT 0 then PAmax=double(PAmaxin)
  If n_elements(PAminin) GT 0 then PAmin=double(PAminin)
  IF n_elements(accuracy) EQ 0 then accuracy=[1.,1.]
  IF n_elements(order) EQ 0 then order=dblarr(1)
  IF n_elements(slopedrings) EQ 0 then slopedrings=0.
  IF slopedrings EQ n_elements(PAin[*,0]) then slopedrings=0.
  IF slopedrings NE 0. then slopedrings=slopedrings+1
  for i=0,n_elements(accuracy)-1 do $
     IF accuracy[i] LT 0.1 then accuracy[i]=0.1                          
                                ;If the maxdeviation is not set  we
                                ;want to create one based on the
                                ;accuracy
  IF n_elements(maxdev) LT n_elements(PAmin) then begin
     tmp=dblarr(n_elements(PAmin))
     i=0
     WHILE i LE n_elements(maxdev)-1 do begin
        tmp[i]=maxdev[i]
        i++
     ENDWHILE
     WHILE i LE n_elements(PAmin)-1 do begin        
        tmp[i]=(PAmax[i]-PAmin[i])/10.*SQRT(accuracy[i])
        IF tmp[i] LT ddiv[i]*5*SQRT(accuracy[i]) then tmp[i]=ddiv[i]*5.*SQRT(accuracy[i])
        i++
     ENDWHILE
     maxdev=tmp
     IF keyword_set(debug) then begin
           print,'This is the maximum deviation we have created'
           print,maxdev,ddiv[*]*5*SQRT(accuracy[*]),ddiv[*]*5./SQRT(accuracy[*]),(PAmax[*]-PAmin[*])/5.       
     ENDIF
  ENDIF
                                ;some tracking variables to use
  IF arctan EQ 0 then attempts=0. else attempts=1
                             
  errincfact=1.5
  checkfit=0
  largerad=0.
  fitstat=0
  prevreduce=0.
  errfact=0.75
  finorder=dblarr(2)
                                ;If the radius has more elements than
                                ;the parameters to be smoothed we want
                                ;to cut back the radius

  IF n_elements(RADII) GT n_elements(PA[*,0]) then begin
     tmp=RADII
     addrad=RADII[n_elements(RADII)-1]
     RADII=dblarr(n_elements(PA))
     RADII=tmp[0:n_elements(PA)-1]
     largerad=1.
  ENDIF
  IF n_elements(SBR) GT n_elements(PA[*,0]) then begin
     tmp=SBR
     SBR=dblarr(n_elements(PA))
     SBR=tmp[0:n_elements(PA)-1]
  ENDIF
  if n_elements(cutoff) LT n_elements(PA[*,0]) then begin
     IF n_elements(cutoff) EQ 0 then cutoff=1E-9
     tmp=cutoff
     cutoff=dblarr(n_elements(PA[*,0]))
     last=n_elements(tmp)-1
     cutoff[0:last]=tmp[0:last]
     cutoff[last+1:n_elements(PA)-1]=cutoff[last]
  endif
  IF cutoff[0] EQ 0 then cutoff[0]=cutoff[1]
  IF keyword_set(rev) and fixedrings GE n_elements(RADII)-1 then fixedrings=n_elements(RADII)-2
  if keyword_set(debug) then begin
     print,'A new regularization starts here'
     print,'PA'
     print,PA
     help,PA
     print,'SBR'
     print,SBR
     help,SBR
     print,'RADII'
     print,RADII
     help,RADII
     if n_elements(error) NE 0 then begin
        print,'error'
        print,error
        help,error
     ENDIF
     if n_elements(fixedrings) NE 0 then begin
        print,'fixedrings'
        print,fixedrings
        help,fixedrings
     ENDIF
     if n_elements(DDiv) NE 0 then begin
        print,'Difference'
        print,DDiv
        help,DDiv
     ENDIF
     if n_elements(maxdev) NE 0 then begin
        print,'Maxdev'
        print,maxdev
        help,maxdev
     ENDIF
     if n_elements(cutoff) NE 0 then begin
        print,'Cutoff'
        print,cutoff
        help,cutoff
     ENDIF
     if n_elements(arctan) NE 0 then begin
        print,'Function to fit'
        print,arctan
        help,arctan
     ENDIF
     if keyword_set(nocentral) then begin
        print,'We remove the central value'
     ENDIF
  ENDIF
                                ;IF we want to exclude the central
                                ;value we strip the PA Not for PA and INCL
;  calculate the smoothed profiles
  PAsmooth=dblarr(n_elements(PA[*,0]),N_elements(PA[0,*]))
  for i=0,n_elements(PA[0,*])-1 do begin
     IF PA[0,i] EQ 0. then PAsmooth[0,i]=PA[0,i] else PAsmooth[0,i]=(PA[0,i]+PA[1,i])/2.                               
     for j=1,n_elements(PA[*,0])-2 do begin         
        PAsmooth[j,i]=(PA[j-1,i]+PA[j,i]+PA[j+1,i])/3           
     endfor
     PAsmooth[n_elements(PA[*,0])-1,i]=(PA[n_elements(PA[*,0])-2,i]+PA[n_elements(PA[*,0])-1,i])/2.
  endfor
                                ;If arctan is 1 then we only want smooth and return the new values
  IF arctan EQ 1 OR n_elements(PA[*,0]) LE 5 then begin
     newPA=dblarr(n_elements(PAin[*,0]),n_elements(PAin[0,*]))
     errors=dblarr(n_elements(PAin[*,0]),n_elements(PAin[0,*]))
                                ;We do not smooth small profiles as
                                ;this supresseses warps too much
     IF n_elements(PAin[*,0]) LT 15 then begin
        for i=0,n_elements(PAin[0,*])-1 do begin
           newPA[*,i]=PAin[*,i]
           errors[*,i]=ABS(PA[*,i]-PAsmooth[*,i])
        endfor
        
     endif else begin
        for i=0,n_elements(PAin[0,*])-1 do begin
           newPA[*,i]=PAsmooth[*,i]
           errors[*,i]=ABS(PA[*,i]-PAsmooth[*,i])
        endfor
     endelse
     fitstat=-1.
     arctan=1
     finorder[*]=!values.f_nan
     fixedrings=[fixedrings,fixedrings]
     goto,cleanup
  ENDIF        
                               
restartall:

;IF we have maxdev then first we are going to do a small check on
;ridiculous outliers.
  
  IF n_elements(maxdev) GT 0. then begin
     for j=0,n_elements(PA[0,*])-1 do begin
        lastmod=1
        start=n_elements(PA[*,0])-1
        for i=1,n_elements(PA[*,0])-2 do begin
           diff1=ABS(PAin[i,j]-PAin[i-1,j])
           diff2=ABS(PAin[i,j]-PAin[i+1,j])
           diff3=ABS(PAin[i+1,j]-PAin[i-1,j])
           IF diff1 GT maxdev[j] AND  diff2 GT maxdev[j] AND diff3 LT (diff1+diff2)/2. then begin
              if keyword_set(debug) then begin
                 print,"These are the differences from ring to ring +the value +maxdev"
                 print,diff1,diff2,diff3,i,PA[i],maxdev
              ENDIF
              PA[i,j]=(PAin[i-1,j]+PAin[i+1,j])/2
              lastmod=i
           ENDIF ELSE start=i
        endfor
        IF ABS(PA[n_elements(PA[*,0])-1,j]-PA[n_elements(PA[*,0])-2,j]) GT maxdev[j] then begin
           diff1=(PAin[n_elements(PA[*,0])-2,j]-PAin[n_elements(PA[*,0])-3,j])
           diff2=(PAin[n_elements(PA[*,0])-3,j]-PAin[n_elements(PA[*,0])-4,j])
           IF diff1/diff2 LT 0 then PA[n_elements(PA[*,0])-1,j]=(PAin[n_elements(PA[*,0])-2,j]+PAin[n_elements(PA[*,0])-3,j])/2. else $
              PA[n_elements(PA[*,0])-1,j]=PAin[n_elements(PA[*,0])-2,j]+(PAin[n_elements(PA[*,0])-2,j]-PAin[n_elements(PA[*,0])-3,j])
           resetlasterror=1
           if keyword_set(debug) then print,'last value is out of bounds'
        ENDIF
        IF keyword_set(debug) then begin
           print,start,lastmod,n_elements(PA),PA,fixedrings
        ENDIF
                                ;if we have these extremes untill the last point just take the average
        if lastmod EQ n_elements(PA[*,0])-2 AND (start LT n_elements(PA[*,0])-3 OR start LE n_elements(PA[*,0])-fixedrings AND  start+1 LE n_elements(PA[*,0])-1) then begin
           IF keyword_set(debug) then begin
              print,'big jigs till the end adjusting before'
              print, PA[start+1:n_elements(PA)-1]
           ENDIF
           PA[start+1:n_elements(PA[*,0])-1,j]=TOTAL(PAin[start+1:n_elements(PA[*,0])-1,j])/n_elements(PAin[start+1:n_elements(PA[*,0])-1,j])
           IF keyword_set(debug) then begin
              print,'big jigs till the end adjusting'
              print, PA[start+1:n_elements(PA)-1]
           ENDIF
        ENDIF
     ENDFOR
  ENDIF

  IF keyword_set(debug) then begin
     print,'The Profiles after adjusting for big fluctuation'
     print,PA
  ENDIF






                                ;Then we will calculate the
                                ;errors. The base for this is the
                                ;ratio between SBR and the cutoff and
                                ;the accuracy
  
  errors=dblarr(n_elements(radii),n_elements(ddiv))
  IF keyword_set(noweight) then begin
     errors[*,0]=ddiv[0]
     errors[*,1]=ddiv[1]
  ENDIF ELSE BEGIN
     IF n_elements(RADII) LT 15 then ratio=cutoff[0:n_elements(SBR)-1]/(SBR[*]*RADII[*]) else ratio=cutoff[0:n_elements(SBR)-1]/SBR
 
     tmp=MAX(ratio,min=norm)
     tmp=WHERE(ratio EQ norm)
     ratio=ratio/norm
     ratio[0:tmp]=1
     
     for i=0,n_elements(ddiv)-1 do begin
       
        errors[*,i]=ddiv[i]*ratio[*]*accuracy[i]
        tmp=WHERE(FINITE(errors[*,i]) EQ 0.)
        maxer=MAX(errors[*,i])
        IF tmp[0] NE -1 then errors[tmp,i]=maxer
     endfor
     IF keyword_set(debug) then begin
        print,'These are the initial errors'
        print,errors
     ENDIF
                                ;if the profile is larger than fifteen
                                ;rings we will modify the errors based
                                ;on their distance from the smoothed profile 
     IF n_elements(PA[*,0]) GT 15 then begin
        for i=0,n_elements(PA[0,*])-1 do begin
           rms=ROBUST_SIGMA(PA[*,i]-PAsmooth[*,i])
           tmperrors=ABS(SQRT(ABS(PA[*,i]-PAsmooth[*,i])))
           IF keyword_set(debug) then begin
              print,'These are the errors based on the smoothed profile',rms
              print,tmperrors
           ENDIF
           if n_elements(PA[*,0]) GT 25 then begin
                 mean=MEAN(PA[*,i])
                 tmperrors[fix(n_elements(PA[*,i])/2.):n_elements(PA[*,i])-1]=tmperrors[fix(n_elements(PA[*,i])/2.):n_elements(PA[*,i])-1]+SQRT(ABS(PA[fix(n_elements(PA[*,i])/2.):n_elements(PA[*,i])-1,i]-mean))
           endif
           errors[*,i]=(errors[*,i]+tmperrors[*])/2.
           tmp=WHERE(FINITE(errors[*,i]) EQ 0)
           IF tmp[0] NE -1 then errors[tmp,i]=ddiv[i]
        endfor
        IF keyword_set(debug) then begin
           print,'Are the errors modified on the smoothed profile errors'
           print,errors
        ENDIF
     ENDIF 
                                ;If there are rings that are fixed we
                                ;will set them to ddiv adittionally if
                                ;we have any errors less than ddiv
                                ;they will also be set to ddiv
     if n_elements(fixedrings) GT 0 then begin
        for i=0,n_elements(PA[0,*])-1 do begin
           errors[0:fixedrings-1,i]=ddiv[i]
                                ;if the errors are less than ddiv in
                                ;the inner half set them to ddiv in
                                ;the outer half we average between the
                                ;surrounding errors
           tmp=WHERE(errors[0:fix(n_elements(errors[*,i])/2.),i] LT DDiv[i])
           IF tmp[0] NE -1 then errors[tmp,i]=DDiv[i]
           tmp=WHERE(errors[fix(n_elements(errors[*,i])/2.+1):n_elements(errors[*,i])-1,i] LT DDiv[i])
           IF tmp[0] NE -1 then begin
              tmp=tmp+fix(n_elements(errors[*,i])/2.+1)
              for j=0,n_elements(tmp)-1 do begin
                 IF tmp[j] LT n_elements(errors[*,i])-1 then errors[tmp[j],i]=(errors[fix(tmp[j]-1),i]+errors[fix(tmp[j]+1),i])/2. else errors[tmp[j],i]=errors[fix(tmp[j]-1),i]
              endfor
           ENDIF
        endfor
        IF keyword_set(debug) then begin
           print,'The errors modified for the rings that are fixed'
           print,errors
        ENDIF
     endif
                                ;This part concludes constructing the
                                ;basic error for the fitting. The next
                                ;section will modify the the profile
                                ;based on the saw pattern or ring to
                                ;ring variations. Remember that we
                                ;always have defined a maxdev
     deltaslope=[0.,0.]
     deltasign=[0.,0.]
     prevsign=[0.,0.]
     
                                ;we will check this from the outside
                                ;in and ignore the inner parts
     counter=0
     errratio=errors[0,0]/errors[0,1]
     for j=fix(n_elements(PA[*,0])/2.),n_elements(PA[*,0])-2 do begin
                                ;first check if we satisfy the
                                ;condition
        adjust=0
        sawtooth=0
        maxdevfac=[0.,0.]
        trend=[0,0.]
        for i=0,n_elements(PA[0,*])-1 do begin
           deltaslope[i]=(PA[j+1,i]-PA[j,i])-(PA[j,i]-PA[j-1,i])
           deltasign[i]=(PA[j+1,i]-PA[j,i])/(PA[j,i]-PA[j-1,i])
                                ;Here we will determine the current
                                ;trend of the change
           sign=dblarr(3)
           IF j GT 2 then begin
              for x=0,2 do begin
                 sign[x]=PA[j-x,i]-PA[j-x-1,i]
              endfor
           ENDIF
           tmp=WHERE(sign LT 0)
           case 1 of
              ;if all elements are less than 0 the curve is growing
              n_elements(tmp) EQ 3:begin
                 trend[i]=1
                 IF (PA[j+1,i]-PA[j,i]) GT 0 AND ABS(deltaslope[i]) GT maxdev[i] then adjust=1
                 
              end
              ;IF all elements are larger than 0 the trend is declining
              tmp[0] EQ -1:begin
                 trend[i]=-1
                 IF (PA[j+1,i]-PA[j,i]) LT 0 AND ABS(deltaslope[i]) GT maxdev[i] then adjust=1
                 
              end
              ;IF there are positive and negatives we have a sawtooth pattern
              else:begin
                 trend[i]=0
                 tmp=WHERE(ABS(sign) GT ddiv[i])
                 IF tmp[0] NE -1 then sawtooth=1
                 IF ABS(deltaslope[i]) GT maxdev[i] and deltasign[i] LT 0 then adjust=1
              end
           endcase
           if Keyword_set(debug) then begin
              print,'We obtained these values'
              print,adjust,sawtooth,deltasign[i],trend[i],deltaslope[i],maxdev[i],ddiv[i],sign,j
              help,adjust,sawtooth,deltasign[i],trend[i],deltaslope[i],maxdev[i],ddiv[i],sign,j
           endif
        endfor
        
        IF adjust then begin
           IF j EQ n_elements(PA[*,0])-2 then begin
              IF sawtooth then begin
                 tmp=ABS((PA[j,*]-PA[j+1,*])/2.)
                 IF tmp[0]/2. LT Ddiv[0] then errors[j+1,0]=errors[j,0] else  errors[j+1,0]=tmp[0]
                 IF tmp[1]/2. LT Ddiv[1] then errors[j+1,1]=errors[j,1] else  errors[j+1,1]=tmp[1]
           ;      errors[j+1,*]=ABS((PA[j,*]-PA[j+1,*])/2.)
                 IF ABS((PA[j,0]-PA[j+1,0])) GT maxdev[0]*2 then begin
                    errors[j+1,0]=ABS((PA[j,0]-PA[j+1,0]))*(ABS((PA[j,0]-PA[j+1,0])/maxdev[0]))
                    ;*ratio[j+1]
                    errors[j+1,1]=ABS((PA[j,1]-PA[j+1,1]))*(ABS((PA[j,0]-PA[j+1,0])/maxdev[0]))
                    ;*ratio[j+1]
                 ENDIF
                 IF ABS((PA[j,1]-PA[j+1,1])) GT maxdev[1]*2 then begin
                    errors[j+1,0]=ABS((PA[j,0]-PA[j+1,0]))*(ABS((PA[j,1]-PA[j+1,1])/maxdev[1]))
                    ;*ratio[j+1]
                    errors[j+1,1]=ABS((PA[j,1]-PA[j+1,1]))*(ABS((PA[j,1]-PA[j+1,1])/maxdev[1]))
                    ;*ratio[j+1]
                 ENDIF
                 IF ABS((PA[j,0]-PA[j+1,0])) GT maxdev[0]  OR ABS((PA[j,1]-PA[j+1,1])) GT maxdev[1] then begin
                    errors[j+1,0]=errors[j+1,0]^2
                    errors[j+1,1]=errors[j+1,1]^2
                 ENDIF
                     
              ENDIF ELSE BEGIN                 
                 IF ABS((PA[j,0]-PA[j+1,0])) GT maxdev[0] then begin
                    errors[j+1,0]=errors[j+1,0]*((ABS((PA[j+1,0]-PA[j,0])-(PA[j,0]-PA[j-1,0])))/(maxdev[0]))^2.5
                    errors[j+1,1]=errors[j+1,1]*((ABS((PA[j+1,0]-PA[j,0])-(PA[j,0]-PA[j-1,0])))/(maxdev[0]))^2.5
                 ENDIF
                 IF ABS((PA[j,1]-PA[j+1,1])) GT maxdev[1] then begin
                    errors[j+1,0]=errors[j+1,0]*((ABS((PA[j+1,1]-PA[j,1])-(PA[j,1]-PA[j-1,1])))/(maxdev[1]))^2.5
                    errors[j+1,1]=errors[j+1,1]*((ABS((PA[j+1,1]-PA[j,1])-(PA[j,1]-PA[j-1,1])))/(maxdev[1]))^2.5
                 ENDIF              
              ENDELSE
              errors[j+1,0]=errors[j+1,0]*maxdev[0]
              errors[j+1,1]=errors[j+1,1]*maxdev[1]
           ENDIF
         
           IF sawtooth then begin
              IF counter GT 0 then errors[j,*]=ABS(ABS(PA[j-1,*]+PA[j+1,*])/2.-PA[j,*])/2.
              maxdevfac[*]=ABS(deltaslope[*]/(maxdev[*]/accuracy[*]))
              tmp=WHERE(maxdevfac LT 1.2)
              IF tmp[0] NE -1 then maxdevfac[tmp]=1.2
              penalty=MAX(maxdevfac)
              IF MAX([errors[j,0],errors[j,1]*errratio]) EQ errors[j,0] then begin
                 IF keyword_set(debug) then begin
                    print,'We are applying the 0 point errors'
                    print,j,errors[j,0]*(1+penalty*(counter/(n_elements(PA[*,0])/4.))),errratio,counter
                 endif 
                 errors[j,0]=errors[j,0]*(1+penalty*(counter/(n_elements(PA[*,0])/4.)))
                 errors[j,1]=errors[j,0]*(1+penalty*(counter/(n_elements(PA[*,0])/4.)))/errratio
              ENDIF ELSE BEGIN
                 IF keyword_set(debug) then begin
                    print,'We are applying the 1 point errors'
                    print,j,errors[j,1]*(1+penalty*(counter/(n_elements(PA[*,0])/4.))),errratio,counter
                 ENDIF
                 errors[j,0]=errors[j,1]*(1+penalty*(counter/(n_elements(PA[*,0])/4.)))*errratio
                 errors[j,1]=errors[j,1]*(1+penalty*(counter/(n_elements(PA[*,0])/4.)))
              ENDELSE
             ; print,maxdevfac,errors[j,*]
              counter++
           ENDIF ELSE BEGIN                    
              IF ABS(deltaslope[0]) GT maxdev[0]/accuracy[0] then begin
                 errors[j,0]=errors[j,0]*(ABS(deltaslope[0])/(maxdev[0]))
                 errors[j,1]=errors[j,1]*(ABS(deltaslope[0])/(maxdev[0]))
              ENDIF
              IF ABS(deltaslope[1]) GT maxdev[1]/accuracy[1] then begin
                 errors[j,0]=errors[j,0]*(ABS(deltaslope[1])/(maxdev[1]))
                 errors[j,1]=errors[j,1]*(ABS(deltaslope[1])/(maxdev[1]))
              ENDIF
          
                                ; ENDIF
              prevsign=0
           ENDELSE
        ENDIF ELSE BEGIN
           IF sawtooth then begin
              errors[j,*]=ABS(ABS(PA[j-1,*]+PA[j+1,*])/2.-PA[j,*])*ratio[j]
              IF errors[j,0] LT ddiv[0] then errors[j,0]=ddiv[0]
              IF errors[j,1] LT ddiv[1] then errors[j,1]=ddiv[1]
           endif
        ENDELSE
     ENDFOR
  ENDELSE
   IF keyword_set(debug) then begin
     print,'The errors after correcting sawtooth are:'
     print,errors
  ENDIF
                                ; We will also penalize the sloped rings and beyond
  IF slopedrings NE 0 AND slopedrings-1 LE n_elements(errors[*,0])-1 then begin
     for i=slopedrings-1,n_elements(errors[*,0])-1 do begin
        errors[i,*]=MAX([errors[i,*]*1.75,errors[i-1,*]*1.75])
     endfor
  ENDIF
                                ; the maximum error allowed is maxdev
                                ; except when SBR is less than cutoff
                                 ;No errors can be larger then maxdev
  IF n_elements(maxdev) GT 0 then begin
     for i=0,1 do begin
        tmp=WHERE(errors[*,i] GT 2.*maxdev[i])
        if tmp[0] NE -1 then begin
           for j=0,n_elements(tmp)-1 do begin
              IF cutoff[tmp[j]] LT SBR[tmp[j]] then errors[tmp[j],i]=2*maxdev[i]
           endfor
           tmp=WHERE(errors[*,i] GT 3.*maxdev[i])
           IF tmp[0] NE -1 then errors[tmp,i]=3*maxdev[i]
        ENDIF
     endfor
  ENDIF
                                ;Where the ratios SBR/cutoff is less then 1.5 we want the 3*maxdev
  for i=0,1 do begin
     j=n_elements(errors[*,i])-1
     
     WHILE SBR[j]/cutoff[j] LT 1.5 do begin
        IF keyword_set(debug) then begin
           print,'Here we are penalising low SBR rings'
           print,j,SBR[j],cutoff[j],SBR[j]/cutoff[j]
        ENDIF
        errors[j,i]=3*maxdev[i]/(SBR[j]/cutoff[j])
        j--
     ENDWHILE
  endfor
                                ;errors should kind off be increasing unless there is a drop in SBR
;  for j=0,1 do begin
;     for i=1,n_elements(SBR)-1 do begin
;        if errors[i,j] LT errors[i-1,j] AND SBR[i]/cutoff[i] LT SBR[i-1]/cutoff[i-1] then begin
;           errors[i,j]=errors[i-1,j]
;        ENDIF
;        if n_elements(ddiv) GT 0 then begin
;           if errors[i,j] GT 2*ddiv[j]*errors[i-1,j] then begin
;              errors[i,j]=2*ddiv[j]*errors[i-1,j]
;           ENDIF
;        ENDIF
;        
 ;    endfor
;  endfor
  IF keyword_set(extending) then errors[n_elements(errors[*,0])-1,*]=errors[n_elements(errors[*,0])-2,*]
; There can't be zeros in the errors let's check
  for i=0,1 do begin
     tmp=WHERE(errors[*,i] LT ddiv[i])
     IF tmp[0] NE -1 THEN begin
        for j=0,n_elements(tmp)-1 do begin
           IF tmp[j] NE n_elements(errors[*,i])-1 then errors[tmp[j],i]=ddiv[i] else errors[tmp[j],i]=errors[tmp[j]-1,i]
        endfor
     endif
  endfor
  
  IF keyword_set(debug) then begin
     print,'The final errors are:'
     print,errors
  ENDIF
  

                                ;we need to determine where the profile falls below the cutoff
  IF n_elements(cutoff) GT 0. then begin
     for i=n_elements(SBR)-1,0,-1 do begin
        IF SBR[i] GT cutoff[i] then break
     endfor
     cutoffring=i 
     IF cutoffring+1 LT n_elements(SBR)-1 then cutoffring2=cutoffring+2 else cutoffring2=n_elements(SBR)-1
  endif else begin
     cutoffring=n_elements(SBR)-1
     cutoffring2=n_elements(SBR)-1
  ENDELSE
  IF cutoffring GT n_elements(PA[*,0])-1 then cutoffring=n_elements(PA[*,0])
  IF cutoffring2 GT n_elements(PA[*,0])-1 then cutoffring2=n_elements(PA[*,0])
  IF cutoffring LT 1 then cutoffring=1
  IF cutoffring2 LT 1 then cutoffring2=1

  if keyword_set(debug) then begin
     print,'This is the cutoff ring'
     print,cutoffring
     print,'This is the number of rings'
     print,n_elements(SBR)-1
  ENDIF
  tmp=fixedrings
  fixedrings=dblarr(2)
  fixedrings[*]=tmp
  IF n_elements(DDiv) NE 0 then begin
     for j=0,n_elements(PA[0,*])-1 do begin
        diff=dblarr(n_elements(PA[*,0]))
        diff[*]=PA[*,j]-PA[0,j]
                                ;If there are only a few rings in the
                                ;model we only want to apply ddiv
                                ;check in the inner half of the rings
        IF n_elements(PA[*,j]) LT 15 then ddivend=fix(n_elements(PA[*,j])/2.) else ddivend=fix(n_elements(PA[*,j]))-1
        tmpchange=WHERE(ABS(diff[0:ddivend]) LT DDiv[j])
        count=0.
        IF n_elements(tmpchange) GT 1 then begin
           for i=1,n_elements(tmpchange)-1 do begin
              IF tmpchange[i] - 1 NE tmpchange[i-1] then break else count++
           endfor
           PA[0:tmpchange[count],j]=PA[0,j]
           ;Should fixedrings always be the same for PA and INCL
           IF fixedrings[j] LT tmpchange[count] then fixedrings[j]=tmpchange[count]
        ENDIF
     endfor
  ENDIF
  IF keyword_set(debug) then begin
     print,'PA after flattening out minor changes'
     print,PA
  ENDIF
  
  
                                ;We average the fixed rings to a mean
                                ;value, However for INCL and PA the
                                ;inner four rings are always fixed if
                                ;this is a small part of the whole
                                ;model this is not necesarrily
                                ;accurate and might have to be
                                ;corrected for.
  for i=0,n_elements(PA[0,*])-1 do begin
     IF fixedrings[i] NE 0. then begin
        IF RADIIin[3]/RADIIin[n_elements(PA[*,0])-1] GT 0.2 then begin
           PA[0:fixedrings[i],i]=TOTAL(PA[0:fixedrings[i],i],1)/(fixedrings[i]+1)
        ENDIF else begin
           checkrms=0.
           checkmean=0.
           checkrms=STDDEV(PA[4:9,i])
           print,checkmean
           print,'impossible'
           checkmean=MEAN(PA[4:9,i])
           print,checkmean
           print,'impossible'
           print,i
           print,MEAN(PA[4:9,0])
           print,'WTF'
           
           if keyword_set(debug) then begin
              print,'This is the rms, mean and 0 value'
              print,checkrms,checkmean,PA[0,i]
           ENDIF
           IF checkrms LT DDIV[i] then begin
              print,'This is odd'
              print,checkmean,checkrms
              print,'what'
              print,PA[0:9,i]
              print,'This is not possible'
              print,MEAN(PA[4:9,i])
              IF ABS(PA[0,i]-checkmean) GE checkrms then PA[0:3,i]=checkmean 
           ENDIF
           PA[0:fixedrings[i],i]=TOTAL(PA[0:fixedrings[i],i])/(fixedrings[i]+1)
           
        ENDELSE
     ENDIF
  endfor
  
  IF keyword_set(debug) then begin
     print,'PA After correcting for the fixedrings'
     print,PA
  ENDIF
   
refit:

;  errors[*,0]=cutoff/SBR*Ddiv[0]/(cutoff[1]/SBR[1])
;  errors[*,1]=cutoff/SBR*Ddiv[1]/(cutoff[1]/SBR[1])
  shifterrors=errors
  if keyword_set(debug) then begin
     print,'Input to the polynomial fit'
     print,'PA'
     print,PA
     help,PA
     print,'RADII'
     print,RADII
     help,RADII
     
     print,'error'
     print,errors
     help,errors
     print,'SBR/cutoff'
     print,SBR/cutoff
     
     
  ENDIF
                                ;If all reliable rings are fixed we
                                ;can skip the fitting of polynomials
                                ;and just want to use a straight line
                                ;the same if at this point all values are the same
 
  newPA=PA
  newPA[*]=0.
  coefffound=dblarr(6,n_Elements(PA[0,*]))
  triesorder=0
  fixed=[0,0]
  regendorder=[0,0]
  changedorder:
  for par=0,n_Elements(PA[0,*])-1 do begin
     addminchi=0.
     erroradjusted=0
     newPAcoeff=0.
     tmp=WHERE(PA[0,par] EQ PA[*,par])
     IF fixedrings[par] GE cutoffring+1 OR n_elements(tmp) EQ n_elements(PA[*,par]) then begin           
        IF keyword_set(debug) then begin
           print,'As all rings above the cutoff are fixed we fit a straight line'             
        ENDIF
        newPA[*,par]=PA[*,par]
        newPAcoeff=PA[0,par]
        coefffound[0,par]=newPA[0,par]
        finorder[par]=0
        fixed[par]=1
        goto,allfixed
     ENDIF
     endorder=ceil((n_elements(PA[*,par])-3.)/2.)
        
     IF keyword_set(debug) then begin
        print,'this is the fixedrings',fixedrings[par]
        print,'This is our endorder',endorder
     ENDIF
     IF endorder LT 2 then endorder=2
     IF endorder GT 5 then endorder=5
     maxendorder=endorder
     
     fitPAoriginal=PA[*,par]
     fiterrors=errors[*,par]
     shifterrorsor=errors[*,par]
     newendorder:
     IF endorder GT maxendorder then endorder=maxendorder
     fitPA=fitPAoriginal
     devflat=dblarr(endorder-1)
     Chi=dblarr(endorder-1)
     ReturnPA=1
     mcerrors=dblarr(n_elements(PA[*,0]),endorder-1)
     finalcoeff=dblarr(6,endorder-1)
 
     regendorder[par]=endorder
     for order=endorder,2,-1 do begin
        shifterrors=shifterrorsor
        tmp=dblarr(1)
        IF keyword_set(debug) then begin
           print,'starting loop for order',order
        ENDIF
    
        fitfailed=0.
       
        newPAcoeff=FAT_FIT(RADII,fitPA,order,RCHI_SQR=tmp,errors=fiterrors,STATUS=fitstat,/Monte_Carlo,mc_iters=150000.,$
                           mc_errors=shifterrors,maximum_deviation=maxdev[par],accuracy=accuracy[par],fixedrings=fixedrings[par],$
                           mc_y_min=pamin[par],mc_y_max=pamax[par],mc_maxfails=100.,log=log,newy=ReturnPA,/debug)
      
        newPA[*,par]=fitPA[*]
                                ;let's calculate how much our
                                ;solution deviates from the central
                                ;value        
        devflat[order-2]=TOTAL(ABS(ReturnPA[*]-ReturnPA[0]))
        if fitstat EQ 1 then Chi[order-2]=!values.f_nan else begin
           IF tmp GT 1e9 then Chi[order-2]=1e9 else Chi[order-2]=tmp
           for j=0,order do begin
              finalcoeff[j,order-2]=newPAcoeff[j]
           endfor
            
           mcerrors[*,order-2]=shifterrors[*]
        endelse
        IF keyword_set(debug) then begin
           print,'Printing reduced chi, order'
           print,tmp,order
           print,'Finished this order'
        ENDIF
     endfor
     doneorder:
     IF keyword_set(debug) then begin
        print,'Checking the order fitted'
     ENDIF
     checkorderfitted=FINITE(Chi,/NAN)
     IF TOTAL(checkorderfitted) EQ n_elements(Chi) then begin
        IF attempts EQ 0 then begin 
                                ; order=order-1
           IF keyword_set(debug) then begin
              print,'As all Chis are infinite we try to fit the smoothed profile'
           ENDIF
           attempts++
           fitPAoriginal=PAsmooth[*,par]
           shifterrorsor=errors[*,par]
           goto,newendorder
        ENDIF ELSE BEGIN
           IF keyword_set(debug) then begin
              print,'We return the smoothed profile as all Chis are infinite'
           ENDIF
           IF n_elements(newPA) EQ 0 then newPA=dblarr(n_elements(PAin[*,0]),n_elements(PAin[0,*]))
           IF n_elements(errors) EQ 0 then errors=dblarr(n_elements(PAin[*,0]),n_elements(PAin[0,*]))
                                ;We do not smooth small profiles as
                                ;this supresseses warps too much
           IF n_elements(PAin[*,0])-fixedrings[0] LT 5 then begin
             
                 newPA[*,par]=PAin[*,par]
                 errors[*,par]=ABS(PA[*,par]-PAsmooth[*,par])
                
           endif else begin
                 newPA[*,par]=PAsmooth[*,par]
                 errors[*,par]=ABS(PA[*,par]-PAsmooth[*,par])
           endelse
           fitstat=-1.
           arctan=1
           finorder[par]=!values.f_nan
           fixedrings[par]=fixedrings[par]
           goto,cleanuppar
        ENDELSE
     ENDIF
       
     maxchi=MAX(Chi,MIN=minchi)
     IF keyword_set(debug) then begin
        print,'Max min',maxchi,minchi
     ENDIF

     case 1 of 
        minchi EQ 1E9: begin
           IF erroradjusted NE 1 then begin
              IF attempts EQ 0 then begin 
                                ; order=order-1
                 
                 IF keyword_set(debug) then begin
                    print,'We will try to fit the smoothed profile because the error is not adjusted and minchi is GT 1e9'
                 ENDIF
                 attempts++
                 fitPAoriginal=PAsmooth[*,par]
                 shifterrorsor=errors[*,par]
                 goto,newendorder
              ENDIF ELSE BEGIN
                 IF keyword_set(debug) then begin
                    print,'We will return the smoothed profile because the error is not adjusted and minchi is GT 1e9'
                 ENDIF
                 newPA=dblarr(n_elements(PAin[*,0]),n_elements(PAin[0,*]))
                 errors=dblarr(n_elements(PAin[*,0]),n_elements(PAin[0,*]))
                                ;We do not smooth small profiles as
                                ;this supresseses warps too much
                 IF n_elements(PAin[*,0]) LT 15 then begin
                    for i=0,n_elements(PAin[0,*])-1 do begin
                       newPA[*,i]=PAin[*,i]
                       errors[*,i]=ABS(PA[*,i]-PAsmooth[*,i])
                    endfor
                 endif else begin
                    for i=0,n_elements(PAin[0,*])-1 do begin
                       newPA[*,i]=PAsmooth[*,i]
                       errors[*,i]=ABS(PA[*,i]-PAsmooth[*,i])
                    endfor
                 endelse
                 fitstat=-1.
                 arctan=1
                 finorder[*]=!values.f_nan
                 fixedrings=[fixedrings,fixedrings]
                 goto,cleanup
              ENDELSE
           ENDIF ELSE BEGIN
              fiterrors=fiterrors*1.2
              IF keyword_set(debug) then begin
                 print,'We go to newendorder because the error is adjusted and minchi is GT 1e9'
              ENDIF
              goto,newendorder
           ENDELSE
        end
        maxchi LT 1:begin 
           IF addminchi EQ 1 then begin
              fiterrors=shifterrors/errfact
              addminchi=2
              goto,newendorder
           ENDIF
          
          
           
           
           adjustzeroth=0.
           erroradjusted=0
           xxx=MAX([minchi,0.8],min=reduceerr)
           reduceerr=round(reduceerr*10.)/10.
           IF prevreduce then reduceerr=reduceerr*2.
           If reduceerr LT 0.1 then reduceerr=0.1
           If reduceerr GT 0.9 then reduceerr=0.9
           
           minChi=1e9
           
           fiterrors=fiterrors*reduceerr
           IF keyword_set(debug) then begin
              print,'Is it too small tmp',maxchi,n_elements(fitPA),reduceerr 
           ENDIF

                                ;for the rotation curve we always want
                                ;the higher polynomials to be considered
           
           
           
           erroradjusted=1
           IF keyword_set(debug) then begin
              print,'We go to refit cause maxchi is too small'
           ENDIF
           goto,newendorder
           
              
        end
        minchi LT 1:begin
           IF addminchi EQ 1 then begin
              IF keyword_set(debug) then begin
                 print,'We go to newendorder cause addminchi EQ 1'
              ENDIF
              fiterrors=shifterrors/errfact
              addminchi=2.
              
              goto,newendorder
           ENDIF
                                ;for the rotation curve we always want
                                ;the higher polynomials to be considered
                              
        ;   xxx=MAX([minchi,0.8],min=reduceerr)
        ;   reduceerr=round(reduceerr*10.)/10.
        ;   IF prevreduce then reduceerr=reduceerr*2.
        ;   If reduceerr LT 0.1 then reduceerr=0.1
        ;   If reduceerr GT 0.9 then reduceerr=0.9
        ;   fiterrors=fiterrors*reduceerr
        ;   erroradjusted=1
        ;   IF keyword_set(debug) then begin
        ;      print,'We go to newendorder cause in minchi is too small'
        ;   ENDIF
        ;   prevreduce=1
                                ;   goto,newendorder

                                ;IF we are here it means that at least
                                ;a low order does have a reasonable
                                ;chi square since we are now in PA and
                                ;INCL only it is ok to not use the
                                ;higher orders.
                                ; ENDIF
           IF keyword_set(debug) then begin
              print,'We gobbl',minchi
           ENDIF
           WHILE minchi LT 1 do begin  
              IF keyword_set(debug) then begin
                 print,'adjusting dunno'
              ENDIF
              order=WHERE(Chi EQ minchi)
              Chi[order]=maxchi+10.
              xxx=MAX(Chi,MIN=minchi)
              IF keyword_set(debug) then begin
                 print,'adjusting chi'
              ENDIF
           ENDWHILE
        END
        minchi GT 20 and addminchi LT 1:begin
           incerr=minchi/20
           If incerr GT 10 then incerr=10.
           fiterrors=fiterrors*incerr
           tmp=WHERE(fiterrors GT 3.*maxdev[par])
           IF tmp[0] NE -1 then fiterrors[tmp]=3.*maxdev[par]
           IF keyword_set(debug) then begin
              print,'We go to newendorder cause minchi is too big'
           ENDIF
           addminchi=1
           goto,newendorder
        end
        else:begin
           IF keyword_set(debug) then begin
              print,'Minimum and maximum Chi all good'
           ENDIF
        end
     endcase
     adjustzeroth=0.
                                ;If we are doing the INCLination we
                                ;want the variations to be as small as
                                ;possible and hence we penalize the
                                ;order furthest removed from flat
                                ;However we only want to do this if
                                ;the original input is fairly flat
                                ;else we just penalize the higher orders
   
     if par EQ 1  then begin
        devflatin=TOTAL(ABS(PA[*,par]-PA[0,par]))/(n_elements(PA[*,par])-fixedrings[par])
        IF devflatin LT 4.*ddiv[par]and n_elements(PAin[*,0]) LT 15 then begin
                                ;first we normalize the deviations from flat
           if keyword_set(debug) then begin
              print,'this is devflat before normalising'
              print,devflat
           endif
           flatmax=MAX(devflat,min=flatmin)
           devflat[*]=devflat[*]/(n_elements(ReturnPA)-fixedrings[par])
           if keyword_set(debug) then begin
              print,'This is the normalised deviation from flat before taking out flats'
              print,devflat
           ENDIF
           
           tmp=WHERE(devflat LT ddiv[pa]/2.)
           tmp2=WHERE(devflat LT ddiv[pa]/4.)
           
           if tmp[0] NE -1 then devflat[tmp]=devflat[tmp]*1.25
           if tmp2[0] NE -1 then devflat[tmp2]=devflat[tmp2]*5.
           if keyword_set(debug) then begin
              print,'This is the normalised deviation from flat for the orders'
              print,devflat
              print,'This is the Chi before'
              print,Chi
           ENDIF
                                ;make sure there is no devflat less than 1
;           tmp=WHERE(devflat LT ddiv[par])
;           ordercor=SQRT(findgen(endorder-1)+1)
;           if tmp[0] NE -1 then devflat[tmp]=devflat[tmp]*ordercor[tmp]
                                ; and then penalize the ones deviating
                                ; most as well as the higher orders if
                                ; their devflat is less than ddiv
;However we will reduce the Chi's if the PA is not 0.
                                ;      if finorder[0] NE 0 then
           chi[*]=chi[*]*devflat[*]
           ;*ordercor[*]
           if keyword_set(debug) then begin
              print,'this is the chi after correcting for flatness'
              print,chi
           endif
        
        endif else begin
         ;  IF Chi[*]=Chi[*]*SQRT(findgen(endorder-1)+1.)
        endelse
        ;However we will reduce the Chi's if the PA is not 0.
        
        if finorder[0] EQ 0  then Chi[*]=Chi[*]*SQRT(findgen(endorder-1)+1.)*2.
        ;SQRT(findgen(endorder-1)+2)
        maxchi=MAX(Chi,min=minchi)
        
     endif
     
                                ;Let's check a 0th order
                                ;polynomial if it is not the
                                ;rotation curve
     
     adjustederr=0.
     straighterrors=fiterrors
     againzero:
     IF keyword_set(debug) then begin
     
        print,'These errors were used for straight'
        print,straighterrors
     ENDIF
     newPAcoeff=FAT_FIT(RADII,fitPA,0.,RCHI_SQR=redChi,errors= straighterrors,STATUS=fitstat,log=log)      
     IF redChi LT 1 then begin
        tmpeq=WHERE(fitPA NE PA[0])
        IF tmpeq[0] EQ -1 then begin
           redchi=1.0001
           goto,toogood
        ENDIF
        xxx=MAX([redChi,0.85],min=reduceerr)
        IF keyword_set(debug) then begin
           print,'Is it too small tmp',tmp,n_elements(fitPA) ,fiterrors,redChi
        ENDIF
        adjustederr=1
        straighterrors=straighterrors*reduceerr
        goto,againzero
     ENDIF else begin
        tmp=WHERE(fitPA NE PA[0,par])
        IF not FINITE(TOTAL(newPAcoeff)) and tmp[0] EQ -1 then begin
           redchi=1.0001
           newPAcoeff[0]=PA[0,par]
        ENDIF
     ENDELSE
     toogood:
     IF keyword_set(debug) then begin
        print,'The 0th order reduced ChisQ is ',redChi,fitstat
        print,'with this constant',newPAcoeff
        print,'The other CHISq are'
        help,Chi
        for hj=0,n_elements(Chi)-1 do print,hj+2,Chi[hj]
     ENDIF
     IF redChi LT minchi AND redChi GT 1. then begin
           
        order=0
        fiterrors= straighterrors
        newPAcoeff=PAin[0,par]
     ENDIF ELSE begin 
           
        order=WHERE(Chi EQ minchi)+2
        tmp=order
        order=intarr(1)
        order=tmp[0]
        newPAcoeff=dblarr(order)
        newPAcoeff=finalcoeff[0:order,order-2]
     ENDELSE
     
     IF fitstat LT 1 then begin
        finorder[par]=order
        newPA[*,par]=newPAcoeff[0]
        for i=1,n_elements(newPAcoeff)-1 do begin
           newPA[*,par]= newPA[*,par]+newPAcoeff[i]*RADII[*]^i 
        endfor
        IF order GT 0 then errors[*,par]=mcerrors[*,order-2] else errors[*,par]=straighterrors[*]
     ENDIF Else begin
        IF attempts GE 1 then begin
           finorder[par]=!values.f_nan
           newPA[*,par]=PAsmooth[*,par]
           errors[*,par]=fiterrors[*]
           arctan=1
        ENDIF else begin
           attempts++
           fitPAoriginal=PAsmooth[*,par]
           shifterrorsor=errors[*,par]
           goto,newendorder
        ENDELSE
     ENDELSE

  
     if keyword_set(debug) then begin
        print,'We find this as the new profile without checking the variations'
        print,newPA[*,par]
        print,'With the order '+string(finorder[par])
     endif

     
                                ;we will check the newPA against ddiv
                                ;again if changes are continuously
                                ;smaller then the order =0
     if order GT 0 then begin
        IF n_elements(DDiv) NE 0 then begin
           
           diff=dblarr(n_elements(PA[*,0]))
           diff[*]=newPA[*,par]-newPA[0,par]
           tmpchange=WHERE(ABS(diff) LT DDiv[par])
           IF n_elements(tmpchange) EQ n_elements(newPA[*,par]) then begin
              order=0
              finorder[par]=order
              newPA[*,par]=PAin[0,par]
              coefffound[0,par]=PAin[0,par]
              errors[*,par]=straighterrors[*]
           ENDIF ELSE coefffound[0:order,par]=finalcoeff[0:order,order-2]
        ENDIF ELSE coefffound[0:order,par]=finalcoeff[0:order,order-2]
     ENDIF else  coefffound[0,par]=newPA[0,par]
     allfixed:
                                ;we will check the newPA against ddiv
                                ;again if changes are continuously
                                ;smaller then the order =0
     
     if keyword_set(debug) then begin
        print,'We find this as the new profile',ddiv[par]
        print,newPA
        print,'With the order '+string(finorder[par])
     endif
     erroradjusted=0
     cleanuppar:
  endfor   

 
  cleanup:
  for i=0,1 do begin
     IF fixedrings[i] GT n_elements(newPA[*,i])-1 then fixedrings[i]=n_elements(newPA[*,i])-1
     IF fixedrings[i] GT 0. AND finorder[i] NE 0. then begin
        IF arctan NE 0 then begin
           newPA[0:fixedrings[i],i]=TOTAL(newPA[0:fixedrings[i],i])/(fixedrings[i]+1.)
        ENDIF ELSE BEGIN
           IF fixedrings[i] GT fix(n_elements(newPA[*,i])/2.) then begin
              if finorder[i] NE 0 then begin
                 newPA[0:fix(n_elements(newPA[*,i])/2.)-1,i]=PAin[0,i]
                 newPA[fix(n_elements(newPA[*,i])/2.):fixedrings[i],i]=(PAin[0,i]+ newPA[fix(n_elements(newPA[*,i])/2.):fixedrings[i],i])/2.
              endif
           
           ENDIF else begin
              If finorder[i] NE 0. then newPA[0:fixedrings[i],i]=PAin[0,i]
           ENDELSE
        ENDELSE
     ENDIF
  endfor
  if keyword_set(debug) then begin
     print,'What up?'
     print,errors
  endif
  for i=0,n_elements(PA[0,*])-1 do begin
     for j=0,n_elements(PA[*,i])-1 do begin
        errors[j,i]=MAX([errors[j,i]*errincfact,ABS(PAin[j,i]-newPA[j,i]),ddiv[i]])        
     endfor
  endfor
  if keyword_set(debug) then begin
     print,'The errors before min max stuff'
     print,errors
  ENDIF
 
  tmp=WHERE(FINITE(newPA) EQ 0.)
  IF tmp[0] NE -1 then newPA[WHERE(FINITE(newPA) EQ 0.)]=PAin[WHERE(FINITE(newPA) EQ 0.)]
  for i=0,n_elements(PA[0,*])-1 do begin 
   ;  IF n_elements(pamin[i]) gt 0 then begin
   ;     tmp=WHERE(newPA[*,i]-errors[*,i] LT pamin[i])
   ;  ENDIF else tmp=-1
   ;  IF tmp[0] NE -1  then begin
   ;     for j=0,n_elements(tmp)-1 do begin
   ;        errors[tmp[j],i]=MAX([newPA[tmp[j],i]-pamin[i],ddiv[i]])
   ;     endfor
   ;  endif
     IF SBRin[n_elements(PAin[*,i])-1] LT cutoff[n_elements(PAin[*,i])-1] then errors[n_elements(PAin[*,i])-1,i]=MAX(errors[*,i])
     If errors[n_elements(errors[*,0])-1,i] LT errors[n_elements(errors[*,0])-2,i] then errors[n_elements(errors[*,0])-1,i]=errors[n_elements(errors[*,0])-2,i]
  endfor
  if arctan eq 0 and checkfit LT 1 then begin
  ;   neworderfind=dblarr(endorder+1)
  ;   newchi=dblarr(endorder+1)
     IF keyword_set(debug) then print,'maybe the sbr is the same of course'
     ratio=SBRin[*]/cutoff
     IF keyword_set(debug) then print,ratio
     tmp=WHERE(ratio LT 2.)
     IF n_elements(tmp) EQ n_elements(ratio) then tmp=WHERE(ratio LT 1.5)
     endrings=n_elements(ratio)-1
     IF keyword_set(debug) then print,'the endrings',endrings
     IF tmp[0] NE -1 then begin
        WHILE tmp[n_elements(tmp)-1] EQ endrings AND n_elements(tmp) GT 3 DO BEGIN
           tmp2=tmp
           tmp=tmp2[0:n_elements(tmp)-2]
           endrings--
        ENDWHILE
     ENDIF
     IF keyword_set(debug) then print,endrings
     
     for j=0,1 do begin
        endorder=regendorder[j]
        newchi=dblarr(endorder+1)
        IF fixed[j] NE 1 then begin
           fitPA=PAin[0:endrings,j]
           IF fixedrings[j] LE endrings then fitPA[0:fixedrings[j]]=PAin[0,j] else begin
              fitPA[0:endrings]=PAin[0,j]
              fixedrings[j]=endrings
           ENDELSE
           if keyword_set(debug) then begin
              print,'This is the PA we will fit'
              print,fitPA
              print,'With these errors'
              print,errors[0:endrings,j]
           ENDIF
           fact=1
           for i=0,endorder do begin
              IF i EQ 1 then begin
                 newchi[i]=5e9
                 i++
              ENDIF
              newcoeff=FAT_FIT(RADII[0:endrings],fitPA,i,RCHI_SQR=redChi,errors= errors[0:endrings,j],STATUS=fitstat,log=log,fixedrings=fixedrings[j])
              IF keyword_set(debug) then print,i,redchi
              newchi[i]=redChi
           endfor
          
           checkerr=MAX(newchi,min=minchi)
           IF minchi LT 1. then begin
              fact=fact/minchi
              errors[*,j]=errors[*,j]/fact
           ENDIF
           tmp=WHERE(newchi EQ minchi)
           refitcounter=0
           refitad:
           refitcounter++
           if keyword_set(debug) then begin
              print,'This is the PA we will fit'
              print,fitPA
              print,'With these errors'
              print,errors[0:endrings,j]
              print,'We selected this order'
              print,tmp[0]
           ENDIF
           newcoeff=FAT_FIT(RADII[0:endrings],fitPA,tmp[0],RCHI_SQR=redChi,errors= errors[0:endrings,j],STATUS=fitstat,log=log,newy=newPAfit,fixedrings=fixedrings[j])
           IF refitcounter LT 100. then begin
              IF tmp[0] EQ 0 AND FINITE(newcoeff[0]) EQ 1 then begin
                 redchi=1.001
                 newcoeff[0]=PAin[0,j]
              ENDIF
              IF keyword_set(debug) then print,redChi,tmp[0],finorder[j]
              IF redChi LT 1. then begin
                 IF redChi GT 0.9 then fact=1.1 else fact=1./(redChi+0.1)
                 errors[*,j]=errors[*,j]/fact
                 goto,refitad
              ENDIF
              IF redCHi GT 3 then begin
                 fact=redChi/2.
                 errors[*,j]=errors[*,j]*fact
                 goto,refitad
              ENDIF
              IF tmp[0] NE finorder[j] then begin
                 errors[*,j]=errors[*,j]*1.5
                 fg=MAX([tmp[0],finorder[j]],min=neworder)
                 if finorder[j] NE 0 then begin
                    tmp[0]=neworder
                    finorder[j]=neworder
                    goto,refitad
                 ENDIF
              ENDIF
           ENDIF
           IF keyword_set(debug) then print,'worrisome'
           IF endrings EQ n_elements(newPA[*,j])-1 then newPA[*,j]=newPAfit else begin
              IF keyword_set(debug) then print,'what ishjhjh'
              newy=dblarr(n_elements(newPA[*,j]))
              newy[*]=newcoeff[0]
              for i=1,n_elements(newcoeff)-1 do newy[*]=newy[*]+newcoeff[i]*RADII[*]^i
              IF fixedrings[j] GT 0. and finorder[j] NE 0 then begin
                 IF fixedrings[j] GT fix(n_elements(newy)/2.)-1 then begin
                    newy[0:fix(n_elements(newy)/2.)-1]=PAin[0,j]
                    newy[fix(n_elements(newy)/2.):fixedrings[j]]=(PAin[0,j]+newy[fix(n_elements(newy)/2.):fixedrings[j]])/2.                    
                 ENDIF else begin
                    newy[0:fixedrings[j]]=PAin[0,j]
                    newy[fixedrings[j]+1]=(newy[fixedrings[j]+1]+newy[fixedrings[j]])/2.
                    newy[fixedrings[j]+2]=(newy[fixedrings[j]+2]*2.+newy[fixedrings[j]+1])/3.
                    newy[fixedrings[j]+3]=(newy[fixedrings[j]+3]*3.+newy[fixedrings[j]+2])/4.
                 ENDELSE
              ENDIF
              newPA[*,j]=newy
           ENDELSE
           finorder[j]=tmp[0]
           IF j EQ 0 then begin
              IF errors[0,j] LT 1 then begin
                 errors[*,j]=errors[*,j]/errors[0,j]
              ENDIF
           ENDIF ELSE BEGIN
              IF errors[0,j] LT 2. then begin
                 errors[*,j]=errors[*,j]/errors[0,j]
              ENDIF
           ENDELSE
           IF keyword_set(debug) then print,'is it finished?',finorder[j]
        ENDIF
     endfor
     checkfit=1
     errincfact=1.
     goto,cleanup
  endif
  

  errorin=errors
  order=finorder
  if keyword_set(debug) then begin
     print,'The estimated  errors'
     print,errorin
  ENDIF
  PAin=newPA
  if keyword_set(debug) then begin
     print,'The final fitted output'
     print,PAin
  ENDIF
  tmp=WHERE(FINITE(Pain) EQ 0)
  IF tmp[0] NE -1 then stop
  IF size(log,/TYPE) EQ 7 then begin
     openu,66,log,/APPEND
     IF arctan EQ 0 then begin
        tmp=where(coefffound[*,0] NE 0.)
        a=findgen(n_elements(tmp))
        printf,66,linenumber()+'REVISED_REGULARISATION_COM: We have fitted the PA curve with a polynomial of order '+strtrim(string(order[0]),2)
        printf,66,linenumber()+'REVISED_REGULARISATION_COM: With the following coefficients '+STRJOIN('c'+strtrim(string(a,format='(I1)'),2)+'='+strtrim(string(coefffound[tmp,0]),2),', ')
        
        printf,66,linenumber()+'REVISED_REGULARISATION_COM: We have fixed the inner '+strtrim(string(fixedrings[0]),2)+' rings'
        tmp=where(coefffound[*,1] NE 0.)
        a=findgen(n_elements(tmp))
        printf,66,linenumber()+'REVISED_REGULARISATION_COM: We have fitted the INCL curve with a polynomial of order '+strtrim(string(order[1]),2)
        printf,66,linenumber()+'REVISED_REGULARISATION_COM: With the following coefficients '+STRJOIN('c'+strtrim(string(a,format='(I1)'),2)+'='+strtrim(string(coefffound[tmp,1]),2),', ')
        
        printf,66,linenumber()+'REVISED_REGULARISATION_COM: We have fixed the inner '+strtrim(string(fixedrings[1]),2)+' rings'
        
     ENDIF ELSE BEGIN
        printf,66,linenumber()+'REVISED_REGULARISATION_COM: We have only smoothed the PA and INCL'
     ENDELSE
     close,66
  ENDIF

end
