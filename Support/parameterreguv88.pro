Pro parameterreguv88,PAin,SBRin,RADIIin,error=errorin,fixedrings=fixedringsin,REVERSE=rev,NOWEIGHT=noweight,Difference=DDivin,cutoff=cutoffin,arctan=arctanin,debug=debug,nocentral=nocentral,order=order,max_deviation=maxdevin,max_par=PAmaxin,min_par=PAminin,accuracy=accuracy,extending=extending

;+
; NAME:
;       PARAMETERREGUV87
;
; PURPOSE:
;       Routine to regularize the parameters of a tirific fit
;
; CATEGORY:
;       Support
; 
; CALLING SEQUENCE:
;        PARAMETERREGUV87, PAin, SBRin, RADIIin, ERROR=errorin,
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
;       15-03-2016 P.Kamphuis; Replaced POLY_FIT with FAT_FIT  
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

;let's make sure no things are changed except PA


  SBR=SBRin
  RADII=RADIIin
  PA=PAin
  If n_elements(errorin) GT 0 then error=errorin
  If n_elements(fixedringsin) GT 0 then fixedrings=fixedringsin else fixedrings=1
  If n_elements(DDivin) GT 0 then DDiv=DDivin
  If n_elements(cutoffin) GT 0 then cutoff=cutoffin
  If n_elements(arctanin) GT 0 then arctan=arctanin
  If n_elements(maxdevin) GT 0 then  maxdev=maxdevin
  If n_elements(PAmaxin) GT 0 then PAmax=PAmaxin
  If n_elements(PAminin) GT 0 then PAmin=PAminin
  IF n_elements(accuracy) EQ 0 then accuracy=1.
  IF n_elements(order) EQ 0 then order=dblarr(1)
  IF accuracy LT 0.1 then accuracy=0.1

                                ;arctan can have the values of 0: for a polynomial fit
                                ;1: for a smoothing of PA or INCL
                                ;2: for a smoothing of VROT
                                ;IF the initial fit fails it will try
                                ;to smooth and fit again
  largerad=0.
  IF n_elements(PAmax) GT 0 then PAmax=double(PAmax)
  IF n_elements(PAmin) GT 0 then PAmin=double(PAmin)
  IF n_elements(PAmin) GT 0 AND  n_elements(PAmax) GT 0 AND n_elements(maxdev) EQ 0 then begin
     maxdev=(PAmax-PAmin)/10.*SQRT(accuracy)
     IF maxdev LT ddiv*5*SQRT(accuracy) then maxdev=ddiv*5.*SQRT(accuracy)
     IF keyword_set(debug) then begin
        print,maxdev,ddiv*5*SQRT(accuracy),ddiv*5./SQRT(accuracy),(PAmax-PAmin)/10.
     ENDIF
  ENDIF
  IF n_elements(RADII) GT n_elements(PA) then begin
     tmp=RADII
     addrad=RADII[n_elements(RADII)-1]
     RADII=dblarr(n_elements(PA))
     RADII=tmp[0:n_elements(PA)-1]
     largerad=1.
  ENDIF
  IF n_elements(SBR) GT n_elements(PA) then begin
     tmp=SBR
     SBR=dblarr(n_elements(PA))
     SBR=tmp[0:n_elements(PA)-1]
  ENDIF
  if n_elements(cutoff) LT n_elements(PA) then begin
     IF n_elements(cutoff) EQ 0 then cutoff=1E-9
     tmp=cutoff
     cutoff=dblarr(n_elements(PA))
     IF n_elements(tmp) GT n_elements(PA) then last=n_elements(PA)-1 else last=n_elements(tmp)-1
     cutoff[0:last]=tmp[0:last]
     if last LT n_elements(PA)-1 then begin
        cutoff[last+1:n_elements(PA)-1]=cutoff[last]
     endif
  endif
  IF cutoff[0] EQ 0 then cutoff[0]=cutoff[1]
  IF n_elements(fixedrings) EQ 0. then fixedrings=0.
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


;We are going to estimate the error based on their distance from the
;smoothed profile If we have more than 15 rings
  restartall:
  IF n_elements(PA) GT 15 then decline=0.9 else decline=0.8
  IF n_elements(PA) GT 15 then begin
     smoothPA=dblarr(n_elements(PA))
     smoothPA[0]=(PA[0]+PA[1])/2.
     for i=1,n_elements(PA)-2 do begin
        smoothPA[i]=(PA[i-1]+PA[i]+PA[i+1])/3.
     endfor
     smoothPA[n_elements(PA)-1]=(PA[n_elements(PA)-2]+PA[n_elements(PA)-1])/2.
     rms=ROBUST_SIGMA(PA-smoothPA)
     IF rms LT 0.5 then rms = 0.5
     IF keyword_set(debug) then begin
        print,'why this?'
        print,rms,ABS((PA-smoothPA))
     ENDIF
     MAXRAD=MAX(RADII)
     errors=ABS(SQRT(ABS(PA-smoothPA))/rms) 
     IF keyword_set(debug) then begin
        print,'Initital errors'
        print,errors
     ENDIF
     IF keyword_set(rev) then begin
        IF keyword_set(nocentral) then begin
           smoothPA[1]=PA[1]
           smoothPA[0]=PA[0]
        ENDIF else BEGIN
           IF PA[0] EQ 0 then smoothPA[0]=0
        ENDELSE
        tmp=MEAN(errors[2:n_elements(PA)-1-fixedrings])
        errors[0:1]=tmp
        errors[n_elements(PA)-fixedrings:n_elements(PA)-1]=tmp
     endIF
     IF keyword_set(extending) then errors[n_elements(errors)-1]=errors[n_elements(errors)-2]
     IF keyword_set(debug) then begin
        print,errors
     ENDIF
     tmp=WHERE(FINITE(errors) EQ 0.)
     IF tmp[0] NE -1 then errors[tmp]=1.                              
     IF keyword_set(rev) then begin
        maxerr=MAX(errors[0:n_elements(errors)-2-fixedrings],min=minerr)
        IF minerr EQ 0. then begin
           minerr=0.1
           errors[WHERE(errors EQ 0)]=0.1
        ENDIF
        errors[n_elements(errors)-1-fixedrings:n_elements(errors)-1]=minerr
     ENDIF ELSE begin
        maxerr=MAX(errors[fixedrings:n_elements(errors)-1],min=minerr)
        IF minerr EQ 0. then begin
           minerr=0.1
           errors[WHERE(errors EQ 0)]=0.1
        ENDIF
        errors[0:fixedrings-1]=minerr/2.
     ENDELSE
     IF keyword_set(debug) then begin
        print,'before 0 and nonfinite'
        print,errors
     ENDIF
     tmp=WHERE(FINITE(errors) EQ 0.)
     IF tmp[0] NE -1 then errors[WHERE(FINITE(errors) EQ 0.)]=maxerr*2
     tmp=WHERE(errors EQ 0.)
     IF tmp[0] NE -1 then errors[WHERE(errors EQ 0.)]=minerr/2.
     errfact=1.
     SBRfact=cutoff[0:n_elements(SBRin)-1]/SBRin
     IF not keyword_set(rev) then sbrmax=MAX(SBRfact[fixedrings:n_elements(SBRin)-1],min=sbrmin) else  sbrmax=MAX(SBRfact[0:n_elements(SBRin)-1-fixedrings],min=sbrmin)
     SBRfact=SBRfact/(sbrmin)
     IF not keyword_set(rev) then begin
        SBRfact[0:fixedrings]=1. 
        maxsbr=MAX(SBRfact[fixedrings+1:n_elements(SBRfact)-1],min=minsbr)
        IF keyword_set(maxdev) then SBRfact[fixedrings+1:n_elements(SBRfact)-1]=(SBRfact[fixedrings+1:n_elements(SBRfact)-1]-minsbr)/(maxsbr/(maxdev/2.))+1. else $
           SBRfact[fixedrings+1:n_elements(SBRfact)-1]=(SBRfact[fixedrings+1:n_elements(SBRfact)-1]-minsbr)/(maxsbr/(5.))+1.
     ENDIF else begin
        maxsbr=MAX(SBRfact[2:n_elements(SBRin)-fixedrings-1],min=minsbr)
        SBRfact[2:n_elements(SBRin)-fixedrings-1]=(SBRfact[2:n_elements(SBRin)-fixedrings-1]-minsbr)/(maxsbr/(5.))+1.
        SBRfact[n_elements(SBRin)-fixedrings:n_elements(SBRin)-1]=2
        SBRFact[0:1]=1

     ENDELSE
     tmp=WHERE(SBRfact LT 1)
     IF tmp[0] NE -1 then SBRfact[tmp]=1

     IF keyword_set(debug) then begin
        print,'before SBRfact'
        print,errors,SBRfact,SBRin
     endif
     errors=errors*SBRfact
     IF keyword_set(rev) then begin
        errmax=MAX(errors,min=errmin)
        errors=((errors-errmin)/((errmax-errmin)/5)+1.)*ddiv
     ENDIF
     IF keyword_set(debug) then begin
        print,'before DDiv'
        print,errors
     ENDIF
     tmp=WHERE(errors[0:fix(n_elements(errors)/2.)] LT DDiv)
     IF tmp[0] NE -1 then errors[tmp]=DDiv
     IF keyword_set(debug) then begin
        print,tmp
     ENDIF
     tmp=WHERE(errors[fix(n_elements(errors)/2.+1):n_elements(errors)-1] LT DDiv)
     IF keyword_set(debug) then begin
        print,tmp
     ENDIF
 

    IF tmp[0] NE -1 then begin
        tmp=tmp+fix(n_elements(errors)/2.+1)
        for j=0,n_elements(tmp)-1 do begin
           IF tmp[j] LT n_elements(errors)-1 then errors[tmp[j]]=(errors[fix(tmp[j]-1)]+errors[fix(tmp[j]+1)])/2. else errors[tmp[j]]=errors[fix(tmp[j]-1)]
        endfor
     ENDIF
     tmp=WHERE(errors LT DDiv)
     IF tmp[0] NE -1 then errors[tmp]=DDiv
     IF keyword_set(debug) then begin
        print,'before maxdev'
        print,errors
     ENDIF
                                ;if the change in slope from 1 ring to
                                ;the next is more tha maxdev we want
                                ;to increase the error
     IF keyword_set(maxdev) then begin
        deltaslope=0
        deltasign=0
        prevsign=0
        
                                ;we will check this from the outside in and ignore the inner parts
        for j=n_elements(PA)-2,fix(n_elements(PA)/2.),-1 do begin
                                ;first check if we satisfy the condition
           deltaslope=(PA[j+1]-PA[j])-(PA[j]-PA[j-1])
           deltasign=(PA[j+1]-PA[j])/(PA[j]-PA[j-1])
         
           IF ABS(deltaslope) GT maxdev then begin
                                ;Are we doing a rotation curve then
                                ;compare it to the total flatness
              if keyword_set(rev) then begin
                 IF keyword_set(debug) then begin
                    print,'We print the increase in error for the reversedfit', (ABS(deltaslope)/(maxdev))^5 
                 ENDIF
                 IF j EQ n_elements(PA)-1 then begin
                    IF ABS(PA[j+1]-PA[j]) GT maxdev then errors[j+1]=errors[j+1]*((ABS(deltaslope))/(maxdev))^5
                 ENDIF
                 errors[j]=errors[j]*((ABS(deltaslope))/(maxdev))^5
                                ;If this is a point in the middle of
                                ;the curve we flatten it
                 IF j LT n_elements(PA)-2 AND j GT 6 then begin
                    IF keyword_set(debug) then begin
                       print,'do we get here normally'
                    ENDIF
                                ;If these points are flattish then we
                                ;just smooth everything to a single value
                    IF (PA[j]-PA[j+1] LT 0.5 AND PA[j]-TOTAL(PA[3:j-1])/n_elements(PA[3:j-2]) LT 10.) OR ABS(PA[j]-TOTAL(PA[3:j-1])/n_elements(PA[3:j-1])) GT maxdev then begin
                       IF keyword_set(debug) then begin
                          print,'do we do this'
                       ENDIF
                       PA[j:n_elements(PA)-1]=TOTAL(PA[3:j-1])/n_elements(PA[3:j-1])
                    ENDIF
                 ENDIF
                 IF keyword_set(debug) then begin
                    print,'this is after'
                    print,      ((ABS(deltaslope))/(maxdev))^5
                 ENDIF
              ENDIF ELSE begin
                                ;for not reversed we  apply this stronger
                                ;if it concerns the last two points or
                                ;if it is a jigsaw
                 ;We need to take care of the last ring
                 IF j EQ n_elements(PA)-2 then begin
                    IF deltasign LT 0 then begin
                       errors[j+1]=ABS((PA[j]-PA[j+1])/2.)
                       ;If it is way out then penalize
                       IF keyword_set(debug) then begin
                          print,prevsign,'this is the diff of the last ring',ABS((PA[j]-PA[j+1])), maxdev*2
                       ENDIF
                       IF ABS((PA[j]-PA[j+1])) GT maxdev*2 then errors[j+1]=ABS((PA[j]-PA[j+1]))*(ABS((PA[j]-PA[j+1])/maxdev))
                    ENDIF ELSE BEGIN
                       IF ABS((PA[j]-PA[j+1])) GT maxdev then errors[j+1]=errors[j+1]*((ABS((PA[j+1]-PA[j])-(PA[j]-PA[j-1])))/(maxdev))^2.5
                    ENDELSE
                 ENDIF
                 IF deltasign LT 0 then begin
                    maxdevfac=ABS(deltaslope/(maxdev/accuracy))
                    IF maxdevfac LT 1 then maxdevfac=1.
                    errors[j]=ABS((deltaslope)/4.)*maxdevfac
                    IF keyword_set(debug) then begin
                       print,'jigsaw',PA[j-1:j+1],deltaslope,errors[j],maxdev,maxdevfac ;   wait,5
                    ENDIF
                 ENDIF ELSE BEGIN
                    IF keyword_set(debug) then begin
                       print,prevsign,'this is the sign before'
                    ENDIF
                    IF ABS(deltaslope) GT maxdev/accuracy then begin
                       errors[j]=errors[j]*(ABS(deltaslope)/(maxdev))
                    ENDIF
                    IF j GE n_elements(PA)-2 then errors[j]=errors[j]*((ABS(deltaslope))/(maxdev))^2.5
                   ; ENDIF
                    prevsign=0
                 ENDELSE
                 IF keyword_set(debug) then begin
                    print,'we penalize',PA[j-1:j+1],deltaslope,errors[j],maxdev,maxdev/accuracy,deltasign ;   wait,5
                 ENDIF
                 IF keyword_set(debug) then begin
                    print,PA[j-1:j],((ABS(deltaslope))),errors[j],maxdev ;   wait,5
                 ENDIF
              ENDELSE
           ENDIF
        endfor
        
        IF keyword_set(rev) then errors[n_elements(PA)-1]=errors[n_elements(PA)-2]
        IF keyword_set(debug) then begin
           print,'final errors'
           print,errors
        ENDIF
     ENDIF
  ENDIF
;IF we have maxdev then first we are going to do a small check on
;ridiculous outliers.
 IF n_elements(maxdev) GT 0. then begin 
    if keyword_set(debug) then print,'we do the maxdev stuff',ABS(PA[n_elements(PA)-1]-PA[n_elements(PA)-2])
    lastmod=1
    start=n_elements(PA)-1
    for i=1,n_elements(PA)-2 do begin
       if keyword_set(nocentral) and i EQ 1 then begin
          diff2=ABS(PAin[i]-PAin[i+1])
          IF diff2 GT maxdev then PA[i:i+1]=(PA[i]+PA[i+1])/2.
       ENDIF ELSE BEGIN
          diff1=ABS(PAin[i]-PAin[i-1])
          diff2=ABS(PAin[i]-PAin[i+1])
          diff3=ABS(PAin[i+1]-PAin[i-1])
          
          if keyword_set(debug) then print,diff1,diff2,diff3,i,PA[i],maxdev
          IF diff1 GT maxdev AND  diff2 GT maxdev AND diff3 LT (diff1+diff2)/2. then begin
             PA[i]=(PAin[i-1]+PAin[i+1])/2
             lastmod=i
             IF keyword_set(debug) then begin
                print,diff1,diff2,diff3,i,PA[i]
             ENDIF
          ENDIF ELSE start=i
       ENDELSE
    endfor
    IF ABS(PA[n_elements(PA)-1]-PA[n_elements(PA)-2]) GT maxdev AND not keyword_set(rev) then begin
       diff1=(PAin[n_elements(PA)-2]-PAin[n_elements(PA)-3])
       diff2=(PAin[n_elements(PA)-3]-PAin[n_elements(PA)-4])
       IF diff1/diff2 LT 0 then PA[n_elements(PA)-1]=(PAin[n_elements(PA)-2]+PAin[n_elements(PA)-3])/2. else $
          PA[n_elements(PA)-1]=PAin[n_elements(PA)-2]+(PAin[n_elements(PA)-2]-PAin[n_elements(PA)-3])
       if keyword_set(debug) then print,'last value is out of bounds'
    ENDIF
    IF keyword_set(debug) then begin
       print,start,lastmod,n_elements(PA),PA,fixedrings
    ENDIF
    ;if we have these extremes untill the last point just take the average
    if lastmod EQ n_elements(PA)-2 AND (start LT n_elements(PA)-3 OR start LE n_elements(PA)-fixedrings AND  start+1 LE n_elements(PA)-1) then begin
       IF keyword_set(debug) then begin
          print,'big jigs till the end adjusting before'
          print, PA[start+1:n_elements(PA)-1]
       ENDIF
       PA[start+1:n_elements(PA)-1]=TOTAL(PAin[start+1:n_elements(PA)-1])/n_elements(PAin[start+1:n_elements(PA)-1])
       IF keyword_set(debug) then begin
          print,'big jigs till the end adjusting'
          print, PA[start+1:n_elements(PA)-1]
       ENDIF
    ENDIF
 ENDIF

  IF keyword_set(debug) then begin
     print,'What on earth'
     print,PA
  ENDIF
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
  IF cutoffring GT n_elements(PA)-1 then cutoffring=n_elements(PA)
  IF cutoffring2 GT n_elements(PA)-1 then cutoffring2=n_elements(PA)
  IF cutoffring LT 1 then cutoffring=1
  IF cutoffring2 LT 1 then cutoffring2=1

  if keyword_set(debug) then begin
     print,'This is the cutoff ring'
     print,cutoffring
     print,'This is the number of rings'
     print,n_elements(SBR)-1
  ENDIF
  IF n_elements(arctan) LT 1 then arctan=0
                                ;We avereage the fixed rings to a mean
                                ;value, However for INCL and PA the
                                ;inner four rings are always fixed if
                                ;this is a small part of the whole
                                ;model this is not necesarrily
                                ;accurate and might have to be
                                ;corrected for.
  IF NOT keyword_set(rev) then begin
     IF fixedrings NE 0. then begin
        IF RADIIin[3]/RADIIin[n_elements(PA)-1] GT 0.2 then PA[0:fixedrings]=TOTAL(PA[0:fixedrings])/(fixedrings+1) else begin
           checkrms=STDDEV(PA[4:9])
           checkmean=MEAN(PA[4:9])
           if keyword_set(debug) then begin
              print,'This is the rms, mean and 0 value'
              print,checkrms,checkmean,PA[0]
           ENDIF
           IF checkrms LT DDIV then begin
              IF ABS(PA[0]-checkmean) GE checkrms then PA[0:3]=checkmean 
           ENDIF
           PA[0:fixedrings]=TOTAL(PA[0:fixedrings])/(fixedrings+1)
        ENDELSE
     ENDIF
  ENDIF else begin
     IF fixedrings NE 0. then begin
        IF keyword_set(debug) then begin
           print,'doing this'
           print,PA
           print,n_elements(PA)-2,n_elements(PA)-fixedrings
        ENDIF
                                ;This procedure is to make sure that the slope of the fixed rings is
                                ;not steeper than ddiv. if it is we average into a flat slope.
        for i=n_elements(PA)-2,n_elements(PA)-fixedrings,-1 do begin
           if ABS(PA[i]-PA[i+1]) GT ddiv then begin
              IF keyword_set(debug) then begin
                 print,ABS(PA[i]-PA[i+1])
                 print,'ye here'
              ENDIF
           ENDIF ELSE goto,skipav
        endfor
        IF keyword_set(debug) then begin
           print,'We average these values'
           Print,PA[i:n_elements(PA)-1]
        ENDIF

        PA[i:n_elements(PA)-1]=TOTAL(PA[i:n_elements(PA)-1])/(n_elements(PA)-i)
     ENDIF
     skipav:
                                ; The while routine here will become
                                ; infinite if there are negative
                                ; values in the array so if these
                                ; things are present we will replace
                                ; them with closest value that is not
                                ; negative
     IF keyword_set(debug) then begin
        print,'This is so weird'
        Print,PA
     ENDIF

     tmp=WHERE(PA LT 0)
     IF tmp[0] NE -1 AND (n_elements(tmp) GT 1 OR tmp[0] NE 0 ) then begin
        IF tmp[0] EQ 0 then PA[tmp]=PA[tmp[1]-1] ELSE  PA[tmp]=PA[tmp[0]-1]
     ENDIF
     for i=fix(n_elements(PA)/2.),n_elements(PA)-2 do begin          
        WHILE PA[i+1] LT PA[i]*decline do begin
           PA[i:i+1]=ABS(PA[i+1]+PA[i])/2.
           
        ENDWHILE
        
     endfor
  ENDELSE
  IF keyword_set(debug) then begin
     print,'PA before  DDiv'
     print,PA
  ENDIF
  IF n_elements(DDiv) NE 0 then begin
     IF NOT keyword_set(rev) then begin
        diff=PA-PA[0]
        tmpchange=WHERE(ABS(diff) LT DDiv)
        count=0.
        IF n_elements(tmpchange) GT 1 then begin
           for i=1,n_elements(tmpchange)-1 do begin
              IF tmpchange[i] - 1 NE tmpchange[i-1] then break else count++
           endfor
           PA[0:tmpchange[count]]=PA[0]
           IF fixedrings LT tmpchange[count] then fixedrings=tmpchange[count]
        ENDIF
     ENDIF ELSE begin
        If fixedrings EQ 0 then begin
           diff=PA-PA[n_elements(PA)-1]
           tmpchange=WHERE(ABS(diff) LT DDiv)
           count=tmpchange[n_elements(tmpchange)-1]
           IF keyword_set(debug) then begin
              print,'Count in DDiv with rev'
              print,count
           ENDIF
           IF  n_elements(tmpchange) GT 1 AND count GE n_elements(PA)-1 then begin 
              for i=n_elements(tmpchange)-1,1,-1 do begin
                 IF tmpchange[i] - 1 NE tmpchange[i-1] then break else count=tmpchange[i-1]
              endfor
              PA[count:n_elements(PA)-1]=TOTAL(PA[count:n_elements(PA)-1])/(n_elements(PA)-count)
              fixedrings=n_elements(PA)-count
           ENDIF
        ENDIF ELSE BEGIN
           diff=PA-PA[n_elements(PA)-fixedrings]
           tmpchange=WHERE(ABS(diff) LT 0.25*DDiv)
           count=tmpchange[n_elements(tmpchange)-1]
           IF keyword_set(debug) then begin
              print,'Count in DDiv'
              print,count,tmpchange
           ENDIF
           IF  n_elements(tmpchange) GT 1 AND count GE n_elements(PA)-fixedrings-1 then begin 
              for i=n_elements(tmpchange)-1,1,-1 do begin
                 IF tmpchange[i] LT n_elements(PA)-fixedrings-1 then begin
                    
                    IF tmpchange[i] - 1 NE tmpchange[i-1] then break else  count=tmpchange[i-1]
                    IF keyword_set(debug) then begin
                       print,count,i
                    ENDIF
                 ENDIF
              endfor
              
              IF count LT n_elements(PA)-1-fixedrings then begin
                 PA[count:n_elements(PA)-1]=TOTAL(PA[count:n_elements(PA)-1])/(n_elements(PA)-count)
                 IF count EQ 0 then fixedrings=n_elements(PA)-1 else fixedrings=n_elements(PA)-count
              ENDIF
              
           ENDIF
        ENDELSE
     ENDELSE

  ENDIF
  IF keyword_set(debug) then begin
     print,'PA after DDiv'
     print,PA
     print,'And count'
     print,count
  ENDIF
  IF n_elements(PA) LE 15 then begin
     errors=dblarr(n_elements(RADII))
  ENDIF
  exptry=0.
  polytry=0.
  arctry=0.
  fitstat=0
  numofit=200.
  fitstat=0.
  prevreduce=0.
  IF arctan EQ 0 then attempts=0. else attempts=1
  RAD=MAX(RADII)
  IF keyword_set(debug) then begin
     print,'What about here',PA
  ENDIF
refit:
  IF keyword_set(noweight) then begin
     errors[*]=ddiv
  ENDIF ELSE BEGIN
     rings=n_elements(RADII)
     
     fact=1.+0.5*accuracy
     IF  n_elements(PA) LE 15 then begin
        IF fixedrings EQ 0. then begin
                                ;The errors should be based on the
                                ;brightness of the ring compared to
                                ;its cutoff however they should always
                                ;increase towards larger radii.

           errors=cutoff/(SBR*RADII) 
           maxerr=MAX(errors[WHERE(FINITE(errors) EQ 1)],MIN=minerr)
           tmp=WHERE(errors GT 2.)
           IF tmp[0] NE -1 then errors[tmp]=errors[tmp]/2.
           tmp=WHERE(errors LT 0.)
           IF tmp[0] NE -1 then errors[tmp]=minerr
           tmp=WHERE(FINITE(errors) EQ 0.)
           IF tmp[0] NE -1 then errors[tmp]=minerr
           maxerr=MAX(errors[1:n_elements(errors)-1],MIN=minerr)
           errors=errors/minerr
           
         ;  IF (n_elements(PA)-fixedrings GT 10 AND fixedrings GT 5) OR n_elements(PA)-fixedrings GT 15 then begin
         ;     for j=1,n_elements(errors)-1 do begin
         ;        IF keyword_set(extending) then begin
         ;           errors[j]=errors[j-1]*fact
         ;        ENDIF ELSE BEGIN
         ;           errors[j]=(errors[j-1]+errors[j])/2.*fact
         ;        ENDELSE
         ;     endfor
         ;  ENDIF else begin
              for j=1,n_elements(errors)-1 do begin
                 WHILE errors[j-1] GT errors[j] do errors[j]=errors[j]*fact
                 
              endfor
         ;  ENDELSE
           IF keyword_set(rev) then begin
              IF keyword_set(nocentral) then begin
                 errors[0:1]=errors[1]*10.
              ENDIF ELSE errors[0]=errors[1]
           ENDIF
        endif else begin
           SBRsmooth=dblarr(n_elements(SBR))
           sbrsmooth[0]=(SBR[0]+SBR[1])/2.
           for i=1,fix(n_elements(SBR)/2.) do begin
              sbrsmooth[i]=(SBR[i-1]+SBR[i]+SBR[i+1])/3.
           endfor
           sbrsmooth[i:n_elements(SBR)-1]=SBR[i:n_elements(SBR)-1]
;Decreasing surface brightness reduce the accuracy but increasing ring
;size increases the accuracy. This should be captured in Cutoff so the
;ratio or difference between cutoff and SBR should be indicative of
;the error.
           errors=cutoff/(SBRsmooth*RADII) 
         ;  errors=cutoff/(SBRsmooth)
           IF keyword_set(debug) then begin
              print,'1'
              print,errors
           ENDIF
           maxerr=MAX(errors[WHERE(FINITE(errors) EQ 1)],MIN=minerr)
           tmp=WHERE(errors GT 2.)    
           IF tmp[0] NE -1 then errors[tmp]=errors[tmp]/2.
           tmp=WHERE(errors LT 0.)
           IF tmp[0] NE -1 then errors[tmp]=maxerr
           tmp=WHERE(FINITE(errors) EQ 0.)
           IF tmp[0] NE -1 then errors[tmp]=maxerr
           IF not keyword_set(rev) then begin
              maxerr=MAX(errors[fixedrings:n_elements(RADII)-1],MIN=minerr)
              errors=errors/minerr
              maxerr=TOTAL(errors[fixedrings:n_elements(RADII)-1])/(n_elements(RADII)-fixedrings-1)^2
;              IF cutoffring LT n_elements(SBR)-1 then begin
 ;                errors[cutoffring:n_elements(SBR)-1]= errors[cutoffring:n_elements(SBR)-1]/3.
  ;            ENDIF
              IF keyword_set(debug) then begin
                 print,'This is the cutoff ring',cutoffring,n_elements(SBR)-1
                 print,errors
              ENDIF
              errors[0:fixedrings-1]=5.
              errors[fixedrings]=1.
              IF (n_elements(PA)-fixedrings GT 10 AND fixedrings GT 5) OR n_elements(PA)-fixedrings GT 15 then begin
                 for j=fixedrings+2,n_elements(errors)-1 do begin
                    IF keyword_set(extending) then begin
                       errors[j]=errors[j-1]*fact
                    ENDIF ELSE BEGIN
                       errors[j]=(errors[j-1]+errors[j])/2.*fact
                    ENDELSE
                 endfor
              ENDIF else begin
                 for j=fixedrings+1,n_elements(errors)-1 do begin
                    WHILE errors[j-1] GT errors[j] do errors[j]=errors[j]*fact
                 endfor
              ENDELSE
              IF keyword_set(debug) then begin
                 print,'This the errors after doing something which shouldnt matter'
                 print,errors
              ENDIF
              IF fixedrings LT n_elements(errors)-2 then maxerr=MAX(errors[fixedrings:n_elements(errors)-2]) else maxerr=MAX(errors)
              If errors[n_elements(errors)-1] GT 5*maxerr then errors[n_elements(errors)-1]=5*maxerr
              IF keyword_set(debug) then begin
                 print,'checking maxerr',maxerr,errors[n_elements(errors)-1]
                 print,errors
              ENDIF
;We also need to penalize the rings that are below the cutoff as we do
;not want to get rid of them completely
              IF cutoffring+1 LT n_elements(errors)-1 then begin
                 errors[cutoffring+1:n_elements(errors)-1]=errors[cutoffring+1:n_elements(errors)-1]*SQRT((cutoff[cutoffring+1:n_elements(errors)-1]/SBR[cutoffring+1:n_elements(errors)-1]))
              ENDIF
              IF keyword_set(debug) then begin
                 print,'checking cutoff',cutoffring+1, n_elements(errors)-1
                 print,errors
              ENDIF
           ENDIF ELSE  BEGIN
              maxerr=MAX(errors[0:n_elements(RADII)-fixedrings-1],MIN=minerr)
              errors=errors/minerr
              IF keyword_set(debug) then begin
                 print,'2'
                 print,errors
              ENDIF
              maxerr=TOTAL(errors[0:n_elements(RADII)-fixedrings-1])/(n_elements(RADII)-fixedrings-1)^2
              IF n_elements(PA)-cutoffring GT fixedrings then startrings=n_elements(PA)-cutoffring else startrings=fixedrings
              IF (n_elements(PA)-fixedrings GT 10 AND fixedrings GT 5) OR n_elements(PA)-fixedrings GT 15 then begin
                 for j=startrings+2,n_elements(errors)-1 do begin
                    errors[j]=errors[j-1]*fact
                 endfor
              ENDIF else begin
                 for j=n_elements(errors)-1-startrings,1,-1 do begin
                    WHILE errors[j-1] LT errors[j] do errors[j-1]=errors[j-1]*fact
                 endfor
              ENDELSE
              IF keyword_set(debug) then begin
                 print,'before no cent' 
                 print,errors   ;       endfor
              ENDIF
              IF not keyword_set(nocentral) then begin
                 maxerr2=TOTAL(errors[0:n_elements(RADII)-fixedrings-1])/(n_elements(RADII)-fixedrings-1)^2
                 errors[0:1]=maxerr2
              ENDIF ELSE BEGIN
                 IF n_elements(RADII)-fixedrings-1 GT 2 then maxerr2=TOTAL(errors[2:n_elements(RADII)-fixedrings-1])/(n_elements(RADII)-fixedrings-3) else  maxerr2=TOTAL(errors[0:n_elements(RADII)-fixedrings-1])/(n_elements(RADII)-fixedrings-1)
                 

                 errors[0:1]=maxerr2*2.
              ENDELSE
              IF keyword_set(debug) then begin
                 print,'after no cent' 
                 print,errors   
              ENDIF
           ENDELSE
        endelse
                                ;also for smaller galaxies we want to penalize large changes in slope
        IF keyword_set(maxdev) then begin
           deltaslope=0.
           deltasign=0
           prevsign=0
                                ;only the outer parts
           for j=n_elements(PA)-2,fixedrings,-1 do begin
              deltaslope=(PA[j+1]-PA[j])-(PA[j]-PA[j-1])
              deltasign=(PA[j+1]-PA[j])/(PA[j]-PA[j-1])
                                ;  prevsign=0.
              IF ABS(deltaslope) GT maxdev  then begin
                                ;reversed or not
                 if keyword_set(rev) then begin
                    IF j EQ n_elements(PA)-2 then begin
                       IF ABS(PA[j+1]-PA[j]) GT maxdev then errors[j+1]=errors[j+1]*((ABS(deltaslope))/(maxdev))^5
                    ENDIF
                    IF j LT n_elements(PA)-1 AND j GT 6 then begin
                       errors[j]=errors[j]*(ABS(deltaslope)/(maxdev))^5
                                ;if it is a jigsaw we want to flatten
                                ;and penalize the previous point
                     ;  IF deltasign LT 0 then begin
                      ;    errors[j-1]=errors[j-1]*(ABS(deltaslope)/(maxdev))^2.5
                      ; ENDIF
                       IF (PA[j+1]-PA[j] LT 0.5 AND PA[j]-TOTAL(PA[3:j-1])/n_elements(PA[3:j-1]) LT 10.) OR ABS(PA[j+1]-TOTAL(PA[3:j-1])/n_elements(PA[3:j-1])) GT maxdev then begin
                          PA[j+1:n_elements(PA)-1]=TOTAL(PA[3:j])/n_elements(PA[3:j])
                       ENDIF
                    ENDIF
                 ENDIF ELSE begin
                                ;in this case we want to penalize all equally
                    divfac=3.
                                ;however if it is a jigsaw we also
                                ;want to penalize the previous point
                    IF j EQ n_elements(PA)-2 then begin
                       IF deltasign LT 0 then begin
                          errors[j+1]=ABS((PA[j]-PA[j+1])/1.5)
                                ;If it is way out then penalize
                          IF keyword_set(debug) then begin
                             print,prevsign,'this is the diff of the last ring',ABS((PA[j]-PA[j+1])), maxdev*2
                          ENDIF
                          IF ABS((PA[j]-PA[j+1])) GT maxdev then begin
                             errors[j+1]=errors[j+1]*(ABS((PA[j]-PA[j+1])/maxdev))
                             divfac=2.
                          ENDIF
                       ENDIF ELSE BEGIN
                          IF ABS((PA[j]-PA[j+1])) GT maxdev then errors[j+1]=errors[j+1]*((ABS((PA[j+1]-PA[j])-(PA[j]-PA[j-1])))/(maxdev))^2.5
                       ENDELSE
                    ENDIF
                    IF deltasign LT 0 then begin
                       maxdevfac=ABS(deltaslope/(maxdev/accuracy))
                       IF maxdevfac LT 1 then maxdevfac=1.
                       errors[j]=ABS((deltaslope)/divfac)*maxdevfac
                       IF keyword_set(debug) then begin
                          print,'jigsaw small',PA[j-1:j+1],deltaslope,errors[j],maxdev,maxdevfac ;   wait,5
                       ENDIF
  ;                     errors[j]=ABS((PA[j]-PA[j-1])/2.)
                      ; errors[j-1]=ABS((PA[j]-PA[j-1])/2.)
                       prevsign=1
                    ENDIF ELSE BEGIN
                       IF keyword_set(debug) then begin
                          print,prevsign,'this is the sign before'
                       ENDIF
                       IF not prevsign then begin
                          errors[j]=errors[j]*(ABS(deltaslope)/(maxdev))^2
                          IF j EQ n_elements(PA)-2 then errors[j+1]=errors[j+1]*((ABS(deltaslope))/(maxdev))
                       ENDIF
                       prevsign=0
                    ENDELSE
                    IF keyword_set(debug) then begin
                       print,'we penalize',PA[j-1:j],deltaslope,errors[j],maxdev,diff[j],deltasign ;   wait,5
                    ENDIF
                 ENDELSE
              ENDIF
           endfor
           IF keyword_set(rev) then errors[n_elements(PA)-1]=errors[n_elements(PA)-2]
        ENDIF
        IF keyword_set(debug) then begin
           print,errors
        ENDIF
        
        errfact=0.5
        tmp=WHERE(errors LT DDiv)
        if tmp[0] NE -1 then begin
           max=MAX(errors,min=Minerr)
           errors=errors*(ddiv/Minerr)
        ENDIF
        if keyword_set(fixedrings) then begin
           if keyword_set(rev) then begin
              errors[n_elements(PA)-fixedrings:n_elements(PA)-1]=ddiv
              IF NOT keyword_set(nocentral) then errors[0]=ddiv
           endif else begin
              errors[0:fixedrings]=ddiv
           ENDELSE
        ENDIF
     ENDIF
     
  ENDELSE
   
   
  IF keyword_set(debug) then begin
     print,'PA after Errors'
     print,PA
     print,'the final errors'
     print,errors
  ENDIF
 
  if cutoffring LT n_elements(RADII)-1 then begin
     IF keyword_set(rev) then begin
        
        IF n_elements(PA) LE 15 then fitPA=PA else fitPA=smoothPA
        
        fitRADII=RADII
        fitSBR=SBR
        fiterrors=errors
        IF fixedrings  GT  n_elements(RADII)-1-cutoffring then begin
           fitPA[cutoffring+1:n_elements(RADII)-1]=TOTAL(fitPA[n_elements(RADII)-fixedrings-1:n_elements(RADII)-1])/(fixedrings+1)
        ENDIF ELSE begin
           fitPA[cutoffring2:n_elements(RADII)-1]=TOTAL(fitPA[cutoffring+1:n_elements(RADII)-1])/(n_elements(RADII)-cutoffring-1)
           fitPA[cutoffring+1]=(TOTAL(fitPA[cutoffring+1:n_elements(RADII)-1])/(n_elements(RADII)-cutoffring-1)+fitPA[cutoffring+1])/2.
        ENDELSE


     ENDIF ELSE BEGIN
        fitPA=PA
        fitRADII=RADII
        fitSBR=SBR
        fiterrors=errors
     ENDELSE
  ENDIF ELSE BEGIN
     fitRADII=RADII
     IF keyword_set(debug) then begin
        print,'we modify PA',PA
     ENDIF
     IF n_elements(PA) LE 15 then fitPA = PA else begin
        IF not keyword_set(rev) then fitPA=PA else fitPA=smoothPA 
     ENDELSE
     fitSBR=SBR
     fiterrors=errors
  ENDELSE
  IF keyword_set(nocentral) then begin
     tmp=fitPA
     PA0=[fitPA[0]]
     fitPA=dblarr(n_elements(tmp)-1)
     fitPA=tmp[1:n_elements(tmp)-1]
     tmp=fitSBR
     SBR0=[fitSBR[0]]
     fitSBR=dblarr(n_elements(tmp)-1)
     fitSBR=tmp[1:n_elements(tmp)-1]
     tmp=fitRADII
     rad0=[fitRADII[0]]
     fitRADII=dblarr(n_elements(tmp)-1)
     fitRADII=tmp[1:n_elements(tmp)-1]
     tmp=fiterrors
     errors0=[fiterrors[0]]
     fiterrors=dblarr(n_elements(tmp)-1)
     fiterrors=tmp[1:n_elements(tmp)-1]
  ENDIF
  IF keyword_set(rev) then begin
     tmp=WHERE(fitPA NE fitPA[0])
     IF tmp[0] EQ -1 then fitPA[0]=fitPA[0]/2.
     fitPA=reverse(fitPA) 
     fiterrors=reverse(fiterrors)
  ENDIF 
  IF largerad EQ 1 then begin
     tmp=RADII
     RADII=dblarr(n_elements(tmp)+1)
     RADII=[tmp,addrad]
  ENDIF
  tmp=dblarr(n_elements(fitPA)-1)
  MAXPA=MAX(fitPA,min=minPA)
  IF n_elements(fitPA) LE 5 then begin
     IF arctan LE 0. then begin
        IF keyword_set(rev) then arctan = 2 else arctan = 1
     endif
   
     attempts = 1
     IF keyword_set(debug) then begin
        print,'did we do this then?'
     ENDIF
  ENDIF
  IF keyword_set(debug) then begin
     print,fitpa,'this here'
  ENDIF
  addminchi=0.
  erroradjusted=0
shifterrors=fiterrors
  case arctan of
     0:begin 
        if keyword_set(debug) then begin
           print,'Input to the polynomial fit'
           print,'PA'
           print,fitPA
           help,fitPA
           print,'RADII'
           print,fitRADII
           help,fitRADII
           print,'error'
           print,fiterrors
           help,fiterrors
        ENDIF
                                ;If all reliable rings are fixed we
                                ;can skip the fitting of polynomials
                                ;and just want to use a straight line
      
       
        endorder=ceil((n_elements(fitPA)-(fixedrings-3))/2)
        
        IF keyword_set(debug) then begin
           print,'this is the fixedrings',fixedrings
           print,'This is our endorder',endorder
        ENDIF
        IF endorder LT 2 then endorder=2
        IF endorder GT 5 then endorder=5
        maxendorder=endorder
        fitPAoriginal=fitPA
        checkPA=dblarr(n_elements(fitPA),10,endorder-1)
        IF keyword_set(debug) then begin
           print,'This is check PA'
           help,checkPA
           print,n_elements(PA),n_elements(fitPA)
           print,'IT should fit', n_elements(fitRADII)
        ENDIF
        shifterrors=fiterrors*errfact
        IF fixedrings GE cutoffring+1 AND  not keyword_set(reverse) then begin           
           order=0
           minchi=1e9
           IF keyword_set(debug) then begin
              print,'As all rings above the cutoff are fixed we fit a straight line'
              
           ENDIF
           goto,allfixed
        ENDIF
        newendorder:
        IF endorder GT maxendorder then endorder=maxendorder
        fitPA=fitPAoriginal
        Chi=dblarr(endorder-1)
        finalcoeff=dblarr(6,endorder-1)
        
        for order=endorder,2,-1 do begin
           tmp=dblarr(1)
           IF keyword_set(debug) then begin
              print,'starting loop for order',order
           ENDIF
           again:
           its=10000.
           maxfailedfits=1000.
           tmpcoeff=dblarr(its*10,order+1)
           tmptmp=dblarr(its*10)
           IF keyword_set(debug) then begin 
              ortmp=dblarr(its*10) 
           ENDIF
           fitfailed=0.
           For x=0,9 do begin
              IF keyword_set(debug) then begin
                 print,'starting x',order
              ENDIF
              mccount=0.
              failedfit=0.
              WHILE mccount LT its do begin
                 
                 newPAcoeff=FAT_FIT(fitRADII,fitPA,order,CHI_SQR=tmp,errors=fiterrors,STATUS=fitstat)
                 IF keyword_set(debug) then begin 
                    ortmp[mccount+x*its]=tmp
                 ENDIF
                 IF fitstat GE 1 then begin
                    IF keyword_set(debug) then begin
                       print,'fitfailed'
                    ENDIF
                    failedfit++
                    IF failedfit GT maxfailedfits then begin
                       Chi[order-2]=!values.f_nan
                       order++
                       IF keyword_set(debug) then begin
                          print,'is it refit'
                       ENDIF
                       IF order LT endorder then goto,again else goto,doneorder
                    ENDIF ELSE goto,skipthisfit
                 ENDIF
                 
                                ; As we ignore underdetermined fits and
                                ; penalization actualy might lift the
                                ; reduced Chi^2 above 1 we skip
                                ; penalizing when the fit is underdetermined
                 
                 newPA=dblarr(n_elements(fitRADII))
                 newPA[*]=newPAcoeff[0]
                 for i=1,n_elements(newPAcoeff)-1 do begin
                    newPA[*]=newPA[*]+newPAcoeff[i]*fitRADII[*]^i
                 endfor
                 
                 IF tmp LT (n_elements(fitPA)-order) then begin
                    goto,skippenalize
                 ENDIF
                 IF keyword_set(maxdev) then begin
                    deltaslope=0.
                    IF keyword_set(rev) then begin
                       for j=0,fix(n_elements(newPA)/2.) do begin
                          deltaslope=(newPA[j]-newPA[j+1])-(newPA[j+1]-newPA[j+2])
                          IF ABS(deltaslope) GT maxdev then begin
                             tmp=tmp*((ABS(deltaslope))/(maxdev))^5
                          ENDIF
                       endfor
                    ENDIF ELSE BEGIN
                       for j=n_elements(newPA)-2,fix(n_elements(newPA)/2.),-1 do begin
                          deltaslope=(newPA[j+1]-newPA[j])-(newPA[j]-newPA[j-1])
                          deltaslope2=(fitPA[j+1]-fitPA[j])-(fitPA[j]-fitPA[j-1])
                          IF (ABS(deltaslope/deltaslope2) GT 2.5 OR ABS(deltaslope/deltaslope2) LT 1./2.5) AND deltaslope GT maxdev/accuracy then begin
                             IF j EQ n_elements(newPA)-2 then begin
                                tmp=tmp*((ABS(deltaslope))/(maxdev/accuracy))^2*2
                             ENDIF ELSE BEGIN
                                tmp=tmp*((ABS(deltaslope))/(maxdev/accuracy))^2
                             ENDELSE
                          ENDIF
                       endfor
                    ENDELSE
                 ENDIF
                 
                 if n_elements(pamin) GT 0 and n_elements(pamax) GT 0 then begin
                    toolarge=WHERE(newPA GT pamax)
                    toosmall=WHERE(newPA LT pamin)
                    IF keyword_set(rev) and not keyword_set(nocentral) then begin
                       IF toosmall[0] NE -1 AND toosmall[n_elements(toosmall)-1] EQ n_elements(newPA)-1 then begin
                          endsmall=n_elements(toosmall)-2
                          IF newPA[n_elements(newPA)-2] LT 0. then begin
                             IF keyword_set(debug) then begin
                                print,'at order',order
                                print,'we get neg values for ',newPA[n_elements(newPA)-2],n_elements(newPA)-2
                             ENDIF
                             tmp=tmp^2
                          ENDIF
                       ENDIF else endsmall=n_elements(toosmall)-1
                    ENDIF ELSE endsmall=n_elements(toosmall)-1
                    IF toolarge[0] NE -1 then begin
                       for i=0,n_elements(toolarge)-1 do tmp=tmp*SQRT(ABS((newPA[toolarge[i]]-pamax)))+1.
                    ENDIF

                    IF toosmall[0] NE -1 then begin
                       
                       for i=0,endsmall do  tmp=tmp*SQRT(ABS((pamin-newPA[toosmall[i]])))+1.
                       
                    ENDIF
                 ENDIF
                 
;we dislike declining rotation curves so let's check the last 1/4.
                 IF keyword_set(rev) then begin
                    for i=1,fix(n_elements(newPA)/4.) do begin
                       IF newPA[i-1]-newPA[i] LT 0. then begin
                          figgy=(newPA[i-1]-newPA[i])
                          IF figgy GT ddiv then begin
                             tmp=tmp*(figgy)
                          ENDIF 
                       ENDIF
                    endfor
                 ENDIF
                 
                 skippenalize:
                 tmptmp[mccount+x*its]=tmp
                 tmpcoeff[mccount+x*its,*]=newPAcoeff
                 mccount++
                 skipthisfit:
                 for j=0,n_elements(fitPA)-1 do begin
                    fitPA[j]=fitPAoriginal[j]+(randomu(seed,/Double)-0.5)*shifterrors[j]
                 ENDFOR
              ENDWHILE
              checkPA[*,x,order-2]=newPA[*]
           endfor
           newtmp=TOTAL(tmptmp)/n_elements(tmptmp)
           for j=0,order do begin
              finalcoeff[j,order-2]=TOTAL(tmpcoeff[*,j]*tmptmp[*])/TOTAL(tmptmp)
           endfor
           IF keyword_set(debug) then begin
              print,'Printing chi,nopen chi, order and reduce chi'
              print,newtmp,Total(ortmp)/n_elements(ortmp),order,newtmp/(n_elements(fitPA)-order)
           ENDIF
           IF keyword_set(debug) then begin
              help,Chi         
           ENDIF
           IF newtmp GT 1e9 then Chi[order-2]=1e9 else Chi[order-2]=newtmp/(n_elements(fitPA)-order)
        
           IF keyword_set(debug) then begin
              print,'Finished this order'
              
           ENDIF
        endfor
        doneorder:
        IF keyword_set(debug) then begin
           print,'Checking the order fitted'
        ENDIF
        checkorderfitted=FINITE(Chi,/NAN)
        IF TOTAL(checkorderfitted) EQ n_elements(Chi) then begin
           order=-1
           If keyword_set(rev) then arctan=2 else arctan=1
           IF keyword_set(debug) then begin
              print,'We go to refit because all chis are infinite'
           ENDIF
           goto,refit
        ENDIF
       
        maxchi=MAX(Chi,MIN=minchi)
        IF keyword_set(debug) then begin
           print,'Max min',maxchi,minchi
        ENDIF

        case 1 of 
            minchi EQ 1E9: begin
                IF keyword_set(debug) then begin
                   print,'We go toQ 1'
                ENDIF
               IF erroradjusted NE 1 then begin
                  order=0.
                  If keyword_set(rev) then arctan=2 else arctan=1
                  IF keyword_set(debug) then begin
                     print,'We go to refit because the error is not adjusted and minchi is GT 1e9'
                  ENDIF
                  goto,refit
               ENDIF ELSE BEGIN
                  fiterrors=fiterrors*1.2
                  IF keyword_set(debug) then begin
                     print,'We go to newendoreder because the error is adjusted and minchi is GT 1e9'
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
               xxx=MAX([maxchi,0.8],min=reduceerr)
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
               IF keyword_set(rev) then begin
                  erroradjusted=1
                  IF keyword_set(debug) then begin
                     print,'We go to newendorder cause in rev and maxchi is too small'
                  ENDIF
                  goto,newendorder
               ENDIF

               tmp=WHERE(maxchi EQ Chi)
               IF n_elements(PA) LE 15 then begin
                  tmp=tmp[0]+3
               
               endif else begin
                  IF not keyword_set(rev) then endorder=MAX([tmp[0]+2,endorder-1])
               ENDELSE
             ;  endorder=MAX([tmp[0]+2,endorder-1])
               erroradjusted=1
               IF keyword_set(debug) then begin
                  print,'We go to newendorder cause and maxchi is too small'
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
             ;  IF keyword_set(rev) then begin
                  xxx=MAX([minchi,0.8],min=reduceerr)
                  reduceerr=round(reduceerr*10.)/10.
                  IF prevreduce then reduceerr=reduceerr*2.
                  If reduceerr LT 0.1 then reduceerr=0.1
                  If reduceerr GT 0.9 then reduceerr=0.9
                  fiterrors=fiterrors*reduceerr
                  erroradjusted=1
                  IF keyword_set(debug) then begin
                     print,'We go to newendorder cause in minchi is too small'
                  ENDIF
                  prevreduce=1
                  goto,newendorder
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
        erroradjusted=0

                                ;Let's check a 0th order
                                ;polynomial if it is not the
                                ;rotation curve
        IF not keyword_set(rev) then begin
           adjustederr=0.
           allfixed:
           straighterrors=fiterrors
           againzero:
           newPAcoeff=FAT_FIT(fitRADII,fitPA,0.,CHI_SQR=tmp,errors=straighterrors,STATUS=fitstat)
           redChi=tmp/n_elements(fitPA)
           continue0:
           IF tmp LT n_elements(fitPA) then begin
              tmpeq=WHERE(fitPA NE PA[0])
              IF tmpeq[0] EQ -1 then begin
                 redchi=1.0001
                 goto,toogood
              ENDIF
              xxx=MAX([(tmp/(n_elements(fitPA))),0.85],min=reduceerr)
              IF keyword_set(debug) then begin
                 print,'Is it too small tmp',tmp,n_elements(fitPA) ,fiterrors,redChi
              ENDIF
              adjustederr=1
              straighterrors=straighterrors*reduceerr
              goto,againzero
           ENDIF else begin
              IF not FINITE(newPAcoeff) and WHERE(fitPA NE PA[0]) EQ -1 then begin
                 redchi=1.0001
                 newPAcoeff[0]=PA[0]
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
        ENDIF ELSE redChi=1e9
   ;     REDChi=0.5
        IF redChi LT minchi AND redChi GT 1. then begin
          
           order=0
           fiterrors= straighterrors
        ENDIF ELSE begin 
           
           order=WHERE(Chi EQ minchi)+2
           tmp=order
           order=intarr(1)
           order=tmp[0]
    ;       order=4
           newPAcoeff=dblarr(order)
           newPAcoeff=finalcoeff[0:order,order-2]
        ENDELSE
        if keyword_set(debug) then begin
           print,'Output of the polynomial fit'
           print,'PA Coefficients'
           print,newPAcoeff
           help,newPAcoeff
        ENDIF
        IF fitstat LT 1 then begin
           
           newPA=dblarr(n_elements(fitRADII))
           newPA[*]=newPAcoeff[0]
           for i=1,n_elements(newPAcoeff)-1 do begin
              newPA[*]=newPA[*]+newPAcoeff[i]*fitRADII[*]^i 
           endfor
           
           if keyword_set(nocentral) then begin
              tmp=newPA
              newPA=dblarr(n_elements(tmp)+1)
              if keyword_set(rev) then newPA=[tmp,PA0] else newPA=[PA0,tmp]
           endif
        ENDIF Else begin
           If keyword_set(rev) then arctan=2 else arctan=1
           goto,refit
        ENDELSE
     END    
     1: begin 
        attempts++
        IF keyword_set(debug) then begin
           print,'We arrive at the smoothing as'
           print,PA,fitPA,attempts
        ENDIF
        newPA=dblarr(n_elements(fitPA))
        IF keyword_set(rev) then fitPA=REVERSE(fitPA)
        IF fitPA[0] EQ 0. then newPA[0]=fitPA[0] else newPA[0]=(fitPA[0]+fitPA[1])/2.
                                ;first a check that there are not a
                                ;total outlier more than 5*the max
                                ;deviation
        
       
;and then the normal averaging
        for i=1,n_elements(fitPA)-2 do begin
           
           newPA[i]=(fitPA[i-1]+fitPA[i]+fitPA[i+1])/3
           
        endfor
        newPA[n_elements(fitPA)-1]=(fitPA[n_elements(fitPA)-2]+fitPA[n_elements(fitPA)-1])/2.
        if keyword_set(nocentral) then begin
           tmp=newPA
           newPA=dblarr(n_elements(tmp)+1)
           newPA=[PA0,tmp]
        endif
        IF keyword_set(rev) then newPA=REVERSE(newPA)
        IF attempts EQ 1 then begin
           arctan=0
           PA=newPA
           smoothPA=newPA
           IF keyword_set(debug) then begin
              print,'Is it refit 3'
              print,PA,newPA
           ENDIF
           goto,restartall
        ENDIF
        
        fitstat=-1.
        order=!values.f_nan
     end
     2: begin 
        attempts++
        IF keyword_set(debug) then begin
           print,fitPA,'this is the input PA'
        ENDIF
        newPA=dblarr(n_elements(fitPA))
        IF keyword_set(rev) then  fitPA=REVERSE(fitPA)
        IF keyword_set(debug) then begin
           print,fitPA,'this is the input PA for section 2'
        ENDIF
        IF fitPA[0] EQ 0. then newPA[0]=fitPA[0] else newPA[0]=(fitPA[0]+fitPA[1])/2.
        for i=1,n_elements(fitPA)-2 do begin
           newPA[i]=(fitPA[i-1]+fitPA[i]+fitPA[i+1])/3
           
           WHILE newPA[i] LT newPA[i-1]*decline do begin
              IF keyword_set(debug) then begin
                 print,newPA[i],newPA[i-1],'increasing'
              ENDIF
              IF newPA[i] LT 0. then begin
                 IF newPA[i-1] GT 0. then newPA[i]=newPA[i-1] else begin
                    IF n_elements(PAmin) GT 0 then newPA[i]=PAmin*2. else begin
                       tmpind=WHERE(newPA GT 0.)
                       if tmpind[0] NE -1 then newPA[i]=MEAN(newPA[tmpind]) else newPA[i]=1.
                    ENDELSE
                 ENDELSE
                       
              endif else newPA[i]=newPA[i]*1.05
           ENDWHILE
           
        endfor
        newPA[n_elements(fitPA)-1]=(fitPA[n_elements(fitPA)-2]+fitPA[n_elements(fitPA)-1])/2.
        WHILE newPA[n_elements(fitPA)-1] LT newPA[n_elements(fitPA)-2]*decline do begin
           IF newPA[n_elements(fitPA)-1] LT 0 then newPA[n_elements(fitPA)-1]=newPA[n_elements(fitPA)-2] else newPA[n_elements(fitPA)-1]=newPA[n_elements(fitPA)-1]*1.05
        ENDWHILE
 
        if keyword_set(nocentral) then begin
           tmp=newPA
           newPA=dblarr(n_elements(tmp)+1)
           newPA=[PA0,tmp]
        endif
        IF keyword_set(rev) then newPA=REVERSE(newPA)
        IF attempts EQ 1 then begin
           arctan=0
           PA=newPA
           smoothPA=newPA
           IF keyword_set(debug) then begin
              print,'Is it refit 3'
              print,PA,newPA
           ENDIF
           goto,restartall
        ENDIF
        fitstat=-2.
        order=!values.f_nan
     end
     ELSE:begin
        print,'PARAMETERREGUV87: this should never happen'
        stop
     ENDELSE
  ENDCASE
  IF fitstat GT 0. then begin
     IF keyword_set(rev) then newPA=REVERSE(PA) ELSE newPA=PA
  ENDIF
  IF fixedrings GT n_elements(newPA)-1 then fixedrings=n_elements(newPA)-1
  IF fixedrings GT 0. and order NE 0 then begin
     IF not keyword_set(rev) then begin
        IF arctan NE 0 then newPA[0:fixedrings]=TOTAL(newPA[0:fixedrings])/(fixedrings+1.) ELSE BEGIN
           if keyword_set(debug) then print,'we doing this?', fix(n_elements(newPA)/2.),fixedrings
           IF fixedrings GT fix(n_elements(newPA)/2.) then begin
              newPA[0:fix(n_elements(newPA)/2.)-1]=PA[0]
              newPA[fix(n_elements(newPA)/2.):fixedrings]=(PA[0]+ newPA[fix(n_elements(newPA)/2.):fixedrings])/2.
              if keyword_set(debug) then print,'we doing this?', fix(n_elements(newPA)/2.)
           ENDIF else begin
              newPA[0:fixedrings]=PA[0]
           ENDELSE
        ENDELSE
     ENDIF

  ENDIF

  IF n_elements(errorin) NE 0. AND n_elements(checkPA) GT 0. AND order GT 0. AND FINITE(order) then begin     
     errorin=dblarr(n_elements(PA)) 
     tmperr=dblarr(n_elements(PA)) 
     tmperr2=dblarr(n_elements(PA)) 
     IF keyword_set(rev) then begin
        for j=0,n_elements(checkPA[*,0,0])-1 do begin  
           tmperr[j]=STDDEV(checkPA[j,*,order-2])*3.
        endfor
        tmperr=REVERSE(tmperr)
        tmperr2=ABS(PAin-REVERSE(newPA))
        for j=0,n_elements(PA)-1 do errorin[j]=MAX([tmperr[j],tmperr2[j]])
     ENDIF ELSE BEGIN
        for j=0,n_elements(checkPA[*,0,0])-1 do begin
           if keyword_set(nocentral) then begin
              tmperr[j+1]=STDDEV(checkPA[j,*,order-2])*3.
           ENDIF ELSE BEGIN
              tmperr[j]=STDDEV(checkPA[j,*,order-2])*3.
           ENDELSE
        endfor
        tmperr2=ABS(PAin-newPA)
        for j=0,n_elements(PA)-1 do errorin[j]=MAX([tmperr[j],tmperr2[j],ddiv])
     ENDELSE
     IF keyword_set(nocentral) then errorin[0]=DDiv
  ENDIF  else begin
     errorin=dblarr(n_elements(PA)) 
     tmperr=dblarr(n_elements(PA)) 
     tmperr2=dblarr(n_elements(PA))
     IF keyword_set(rev) then begin
        tmperr[*]=DDiv
     ENDIF ELSE BEGIN
        if keyword_set(nocentral) then begin
           tmperr[0:n_elements(fitPA)-1]=shifterrors*1.5
           tmperr[n_elements(fitPA)]=DDiv*1.5
        ENDIF ELSE BEGIN
           tmperr[*]=shifterrors*1.5
        ENDELSE
     ENDELSE

     tmperr2=ABS(PAin-newPA)

     for j=0,n_elements(PA)-1 do begin
         if keyword_set(debug) then begin
            print,'The final errors'
            print,tmperr[j],tmperr2[j]
         ENDIF
        errorin[j]=MAX([tmperr[j],tmperr2[j],ddiv])
        if keyword_set(debug) then begin
            print,'The final errors'
            print,errorin[j]
         ENDIF
     endfor
  ENDELSE

  tmp=WHERE(FINITE(newPA) EQ 0.)
  IF tmp[0] NE -1 then newPA[WHERE(FINITE(newPA) EQ 0.)]=PA[WHERE(FINITE(newPA) EQ 0.)]
  PAin=newPA
  IF keyword_set(rev) then begin
     IF n_elements(PAin) GT 4 then PAin[0]=PAin[1]
     PAin=REVERSE(PAin)
  ENDIF
  if keyword_set(debug) then begin
     print,'The final fitted output'
     print,PAin
  ENDIF
  IF n_elements(pamin) gt 0 then tmp=WHERE(PAin-errorin LT pamin) else tmp=-1
  IF tmp[0] NE -1 then errorin[tmp]=PAin[tmp]-pamin
  IF n_elements(pamax) gt 0 then tmp=WHERE(PAin+errorin GT pamax) else tmp=-1
  IF tmp[0] NE -1 then errorin[tmp]=pamax-PAin[tmp]
  IF n_elements(pamin) gt 0 and n_elements(pamax) gt 0 then tmp=WHERE(errorin LT 0)else tmp=-1
  IF tmp[0] NE -1 then errorin[tmp]=pamax-pamin


end
