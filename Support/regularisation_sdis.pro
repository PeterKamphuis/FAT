Pro regularisation_sdis,PAin,SBRin,RADIIin,error=errorin,fixedrings=fixedringsin,NOWEIGHT=noweight,Difference=DDivin,cutoff=cutoffin,arctan=arctanin,debug=debug,order=order,max_deviation=maxdevin,max_par=PAmaxin,min_par=PAminin,accuracy=accuracy,extending=extending,gdlidl=gdlidl,log=log

;+
; NAME:
;       REGULARISTION_SDIS
;
; PURPOSE:
;       Routine to regularize the dispersion of a tirific fit
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
;       Written 14-09-2018 P.Kamphuis v1.0
;
; NOTE:
;This routine has an extensive debugging mode that can be turned on
;with the keyword, /debug
;     
;-
  COMPILE_OPT IDL2 

; First thing is to check all the input whether it is reasonable

  SBR=SBRin
  RADII=RADIIin
  PA=PAin
  IF n_elements(cutoffin) GT 0 then cutoff=cutoffin[0:n_elements(SBRin)-1]
  If n_elements(fixedringsin) GT 0 then fixedrings=fixedringsin else fixedrings=0
  If n_elements(DDivin) GT 0 then DDiv=DDivin
  If n_elements(errorin) GT 0 then error=errorin
                                ;arctan can have the values of 0: for
                                ;a arctan fit
                                ;1: for a smoothing

  If n_elements(arctanin) GT 0 then arctan=arctanin else arctan=1
  If n_elements(maxdevin) GT 0 then  maxdev=maxdevin
  If n_elements(PAmaxin) GT 0 then PAmax=double(PAmaxin) else PAmax=MAX(PA)
  If n_elements(PAminin) GT 0 then PAmin=double(PAminin) else PAmin=MIN(PA)
  IF PAmax LE PAmin then begin
     PAmax=MAX(PAin,min=PAmin)
  ENDIF
  IF n_elements(accuracy) EQ 0 then accuracy=1.
  IF n_elements(order) EQ 0 then order=dblarr(1)  
  IF accuracy LT 0.1 then accuracy=0.1
  if n_elements(fixedrings) GT 1 then fixedrings=fixedrings[0]
  fact=1.+0.5*accuracy
  order=!values.f_nan
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
           print,maxdev,ddiv[*]*5*SQRT(accuracy[*]),ddiv[*]*5./SQRT(accuracy[*]),(PAmax[*]-PAmin[*])/10.       
     ENDIF
  ENDIF
                                ;some tracking variables to use
  IF arctan EQ 0 then attempts=0. else attempts=1
  largerad=0.
  fitstat=0
  prevreduce=0.
  errfact=0.75
  
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
     cutoff[last:n_elements(PA)-1]=tmp[0:last]
     cutoff[0:last-1]=cutoff[last]
  endif
  IF cutoff[n_elements(PA[*])-1] EQ 0 then cutoff[n_elements(PA[*])-1]=cutoff[n_elements(PA[*])-2]
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
  IF n_elements(PA) GT 5 then begin
     IF PA[0] LT MEDIAN(PA) then PA[0:2] = MEAN(PA[0:2])
     PA[n_elements(PA)-3:n_elements(PA)-1] = MEAN(PA[n_elements(PA)-3:n_elements(PA)-1])
  ENDIF
  PAsmooth=dblarr(n_elements(PA[*]))
  PAsmooth[0]=PA[0]                               
  for j=1,n_elements(PA[*])-2 do begin         
     PAsmooth[j]=(PA[j-1]+PA[j]+PA[j+1])/3     
  endfor
  PAsmooth[n_elements(PA[*])-1]=(PA[n_elements(PA[*])-2]+PA[n_elements(PA[*])-1])/2.
  
                                ;If arctan is 1 then we only want smooth and return the new values
  IF arctan EQ 1 OR n_elements(PA) LE 5 then begin
     newPA=PAsmooth
     errors=ABS(PA[*]-PAsmooth[*])
     fitstat=-1.
     arctan=1
     finorder=!values.f_nan
     goto,cleanup
  ENDIF        
                                ;Then we will calculate the
                                ;errors. The base for this is the
                                ;ratio between SBR and the cutoff
restartall:
  errors=dblarr(n_elements(radii))
  IF keyword_set(noweight) then begin
     errors[*]=ddiv[0]
  ENDIF ELSE BEGIN
     IF n_elements(RADII) LT 15 then ratio=cutoff[0:n_elements(SBR)-1]/(SBR[*]*RADII[*]) else ratio=cutoff[0:n_elements(SBR)-1]/SBR
     tmp=MAX(ratio,min=norm)
     tmp=WHERE(ratio EQ norm)
     ratio=ratio/norm
;     ratio[0:tmp]=1
     IF tmp[0] LT n_elements(ratio)-1 then ratio[tmp[0]+1:n_elements(ratio)-1]=1.
     for i=0,n_elements(ddiv)-1 do errors[*]=ddiv*ratio[*]
     tmp=WHERE(FINITE(errors) EQ 0.)
     maxer=MAX(errors)
     IF tmp[0] NE -1 then errors[tmp]=maxer
     IF keyword_set(debug) then begin
        print,'These are the initial errors'
        print,errors
     ENDIF
                                ;if the profile is larger than fifteen
                                ;rings we will modify the errors based
                                ;on their distance from the smoothed profile 
     IF n_elements(PA) GT 15 then begin
        rms=ROBUST_SIGMA(PA[*]-PAsmooth[*])
        tmperrors=ABS(SQRT(ABS(PA[*]-PAsmooth[*]))/rms)
        IF keyword_set(debug) then begin
           print,'These are the errors based on the smoothed profile'
           print,tmperrors
        ENDIF
        errors[*]=(errors[*]+tmperrors[*])/2.
        tmp=WHERE(FINITE(errors[*]) EQ 0)
        IF tmp[0] NE -1 then errors[tmp]=ddiv
   
        IF keyword_set(debug) then begin
           print,'Are the errors modified on the smoothed profile errors'
           print,errors
        ENDIF
     ENDIF
                                ;If there are rings that are fixed we
                                ;will set the error to that of the
                                ;last unfixed ring. adittionally if
                                ;we have any errors less than ddiv
                                ;they will also be set to ddiv
     if n_elements(fixedrings) GT 0 then begin
       
           errors[0:fixedrings-1]=errors[fixedrings]
                                ;if the errors are less than ddiv in
                                ;the inner half set them to ddiv in
                                ;the outer half we average between the
                                ;surrounding errors
           tmp=WHERE(errors[0:fix(n_elements(errors[*])/2.)] LT DDiv)
           IF tmp[0] NE -1 then errors[tmp]=DDiv
           tmp=WHERE(errors[fix(n_elements(errors[*])/2.+1):n_elements(errors[*])-1] LT DDiv)
           IF tmp[0] NE -1 then begin
              tmp=tmp+fix(n_elements(errors[*])/2.+1)
              for j=0,n_elements(tmp)-1 do begin
                 IF tmp[j] LT n_elements(errors[*])-1 then errors[tmp[j]]=(errors[fix(tmp[j]-1)]+errors[fix(tmp[j]+1)])/2. else errors[tmp[j]]=errors[fix(tmp[j]-1)]
              endfor
           ENDIF
       
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
     deltaslope=[0.]
     deltasign=[0.]
     prevsign=[0.]
     
                                ;we will check this from the outside
                                ;in and ignore the inner parts
     counter=0
     
     for j=fix(n_elements(PA[*])/2.),n_elements(PA[*])-2 do begin
                                ;first check if we satisfy the
                                ;condition
        adjust=0
        sawtooth=0
        maxdevfac=[0.]
        
        deltaslope=(PA[j+1]-PA[j])-(PA[j]-PA[j-1])
        deltasign=(PA[j+1]-PA[j])/(PA[j]-PA[j-1])
                                ;Here we will determine the current
                                ;trend of the change
        sign=dblarr(3)
        for x=0,2 do begin
           sign[x]=PA[j-x]-PA[j-x-1]
        endfor
        tmp=WHERE(sign LT 0)
        if keyword_set(debug) then begin
           print,'This is where the sig is negative'
           print,tmp
        endif
        case 1 of
                                ;if all elements are less than 0 the curve is growing
           n_elements(tmp) EQ 3:begin
              trend=1
              IF (PA[j+1]-PA[j]) GT 0 AND ABS(deltaslope) GT maxdev then adjust=1
           end
                                ;IF all elements are larger than 0 the trend is declining
           tmp[0] EQ -1:begin
              trend=-1
              IF (PA[j+1]-PA[j]) GT 0 AND ABS(deltaslope) GT maxdev then adjust=1
           end
                                ;IF there are positive and negatives we have a sawtooth pattern
           else:begin
              trend=0
              tmp=WHERE(sign GT ddiv)
              IF tmp[0] NE -1 then sawtooth=1
              IF ABS(deltaslope) GT maxdev and deltasign LT 0 then adjust=1
           end
        endcase
           
          
        
        
        IF adjust then begin  
           IF sawtooth then begin
              IF counter GT 0 then errors[j]=ABS(ABS(PA[j-1]+PA[j+1])/2.-PA[j])
              maxdevfac=ABS(deltaslope/(maxdev/accuracy))
              tmp=WHERE(maxdevfac LT 1.2)
              IF tmp[0] NE -1 then maxdevfac[tmp]=1.2
              IF keyword_set(debug) then begin
                 print,'do we dothis'
                 print,j,errors[j],n_elements(errors)-1
              ENDIF
              errors[j]=errors[j]*(1+maxdevfac*(counter/(n_elements(PA[*])/4.)))
              counter++
           ENDIF ELSE BEGIN                    
              IF ABS(deltaslope) GT maxdev/accuracy then begin
                 errors[j]=errors[j]*(ABS(deltaslope)/(maxdev))
              ENDIF
          
                                ; ENDIF
              prevsign=0
           ENDELSE
        ENDIF ELSE BEGIN
           IF sawtooth then begin
              errors[j]=ABS(ABS(PA[j-1]+PA[j+1])/2.-PA[j])
             
           endif
        ENDELSE
     ENDFOR
                                ;If the model has more then 15 rings
                                ;we also want to check the last 4
                                ;elements for the sawtooth pattern
   
   

     
   
 
  ENDELSE
                                ;No errors can be larger then maxdev
  IF n_elements(maxdev) GT 0 then begin
     tmp=WHERE(errors GT 2*maxdev[0])
     IF keyword_set(debug) then begin
        print,tmp,maxdev
        print,'here here'
     ENDIF
     if tmp[0] NE -1 then begin
        for j=0,n_elements(tmp)-1 do begin
           IF cutoff[tmp[j]] LT SBR[tmp[j]] then errors[tmp[j]]=2*maxdev[0]
        endfor
        tmp=WHERE(errors GT 3.*maxdev[0])
        IF tmp[0] NE -1 then errors[tmp]=3*maxdev[0]
     ENDIF
  ENDIF
  tmp=WHERE(errors LT ddiv)
  IF tmp[0] NE -1 then errors[tmp]=ddiv

  IF keyword_set(debug) then begin
     print,'The final errors are:'
     print,errors
  ENDIF
  
;IF we have maxdev then first we are going to do a small check on
;ridiculous outliers.
  IF keyword_set(debug) then begin
     print,'Before Adjusting Big Jigs'
     print, PA
  ENDIF
  PAtmp=PA
  IF n_elements(maxdev) GT 0. then begin
     lastmod=n_elements(PA[*])-2
     start=1
     for i=1,n_elements(PA[*])-2 do begin
        diff1=ABS(PAtmp[i]-PAtmp[i-1])
        diff2=ABS(PAtmp[i]-PAtmp[i+1])
        diff3=ABS(PAtmp[i+1]-PAtmp[i-1])
        IF diff1 GT maxdev AND  diff2 GT maxdev AND diff3 LT (diff1+diff2)/2. then begin
           if keyword_set(debug) then begin
              print,"These are the differences from ring to ring +the value +maxdev"
              print,diff1,diff2,diff3,i,PA[i],maxdev
           ENDIF
           PA[i]=((PAtmp[i-1]+PAtmp[i+1])/2+PAtmp[i])/2.
           lastmod=i
        ENDIF ELSE start=i
     endfor
    
     IF keyword_set(debug) then begin
        print,start,lastmod,n_elements(PA),PA,fixedrings
     ENDIF
                                ;if we have these extremes untill the last point just take the average
     if lastmod EQ n_elements(PA[*])-2  AND (start GT 3 OR start GE fixedrings AND  start+1 GE 1) then begin
        IF keyword_set(debug) then begin
           print,'big jigs till the end adjusting before'
           print, PA[start+1:n_elements(PA)-1]
        ENDIF
        PA[start+1:n_elements(PA[*])-1]=TOTAL(PAin[start+1:n_elements(PA[*])-1])/n_elements(PAin[start+1:n_elements(PA[*])-1])
        IF keyword_set(debug) then begin
           print,'big jigs till the end adjusting'
           print, PA[start+1:n_elements(PA[*])-1]
        ENDIF
     ENDIF
     
  ENDIF

  IF keyword_set(debug) then begin
     print,'The Profiles after adjusting for big fluctuation'
     print,PA
  ENDIF
                                ;we need to determine where the profile falls below the cutoff
  IF n_elements(cutoff) GT 0. then begin
     for i=0,n_elements(SBR)-1 do begin
        IF SBR[i] GT cutoff[i] then break
     endfor
     cutoffring=i 
     IF cutoffring+1 GT 1 then cutoffring2=cutoffring+2 else cutoffring2=0
  endif else begin
     cutoffring=0
     cutoffring2=0
  ENDELSE
  IF cutoffring GT n_elements(PA[*,0])-1 then cutoffring=n_elements(PA[*,0])-1
  IF cutoffring2 GT n_elements(PA[*,0])-1 then cutoffring2=n_elements(PA[*,0])-1
  IF cutoffring LT 1 then cutoffring=0
  IF cutoffring2 LT 1 then cutoffring2=0

  if keyword_set(debug) then begin
     print,'This is the cutoff ring'
     print,cutoffring
     print,'This is the number of rings'
     print,n_elements(SBR)-1
  ENDIF
 
                                ;We average the fixed rings to a mean
                                ;value, However for INCL and PA the
                                ;inner four rings are always fixed if
                                ;this is a small part of the whole
                                ;model this is not necesarrily
                                ;accurate and might have to be
                                ;corrected for.
 
   
refit:

  MAXPA=MAX(PA,min=minPA)
  addminchi=0.
  erroradjusted=0
  shifterrors=errors
  if keyword_set(debug) then begin
     print,'Input to the arctan fit'
     print,'PA'
     print,PA
     help,PA
     print,'RADII'
     print,RADII
     help,RADII
     print,'error'
     print,errors
     help,errors
  ENDIF
                                ;If all reliable rings are fixed we
                                ;can skip the fitting of polynomials
                                ;and just want to use a straight line
 
  newPA=PA
  newPA[*]=0.
  fitPAoriginal=PA[*]
  fiterrors=errors[*]
  shifterrorsor=errors[*]
                                ;IF the galaxy is large then we do not
                                ;want the lower orders as the are
                                ;often pushing the rotation curve in
                                ;the wrong direction
 
  fitPA=PA
                                ;We will always fit a arctangent to
  attempt=0                              ;increase stability
  tryagain:
  para=fat_arctan(RADII,fitPA,error=fiterrors,gdlidl=gdlidl)
  if keyword_set(debug) then begin
     print,'These are the fitted coefficients'
     print,para
  endif
  if para[0] LT RADII[2] OR  para[0] GT RADII[n_elements(RADII)-2] OR $
     para[3] LT PAmin OR  para[3] GT PAmax OR abs(para[2]) GT 20. then begin
     IF size(log,/TYPE) EQ 7 then begin
        openu,66,log,/APPEND
        printf,66,linenumber()+'REGULARISATION_SDIS: The arctan fit has failed because one of the following was true:'
        printf,66,linenumber()+'REGULARISATION_SDIS: '+string(para[0])+' LT '+string(RADII[2])
        printf,66,linenumber()+'REGULARISATION_SDIS: '+string(para[0])+' GT '+string(RADII[n_elements(RADII)-2])
        printf,66,linenumber()+'REGULARISATION_SDIS: '+string(para[2])+' LT '+string(PAmin)
        printf,66,linenumber()+'REGULARISATION_SDIS: '+string(para[2])+' GT '+string(PAmax)
        printf,66,linenumber()+'REGULARISATION_SDIS: '+string(para[1])+' GT '+string(20.)
        printf,66,fitPA
        printf,66,linenumber()+'REGULARISATION_SDIS: With the following coefficients '+STRJOIN('c'+strtrim(string(a,format='(I1)'),2)+'='+strtrim(string(para),2),', ')
    
        close,66
     endif
     if attempt LT 3 then begin
        attempt++
        IF size(log,/TYPE) EQ 7 then begin
           openu,66,log,/APPEND
           printf,66,linenumber()+'REGULARISATION_SDIS: We will try again with a smoother profile'
           close,66
        endif
        case attempt of
           1:begin
              fitPA=PAsmooth
              IF fitPA[0] LT maxpa*0.8 then fitpa[0:1]=maxpa*0.8
           end
           2:begin
              fitPA=PAsmooth
              IF fitPA[0] LT maxpa*0.8 then fitpa[0:1]=maxpa*0.8
              fitPA[0:fix(n_elements(fitPA)/2.)]=MEAN(fitPA[0:fix(n_elements(fitPA)/2.)])
              fitPA[fix(n_elements(fitPA)/2.)+1:n_elements(fitPA)-1]=MEDIAN(fitPA[fix(n_elements(fitPA)/2.)+1:n_elements(fitPA)-1])
           end
           3:begin
              fitPA=PAsmooth
              fitPA[*]=MEAN(fitPA)
              
           end
           else:
        endcase
        goto,tryagain             
     endif
     fitPA[n_elements(fitPA)-3:n_elements(fitPA)-1]=MEAN(fitPA[n_elements(fitPA)-3:n_elements(fitPA)-1])
     newPA=PAsmooth
     IF newPA[0] LT maxpa*0.8 then newpa[0:1]=maxpa*0.8
     arctan=1 
  endif else begin
     newPA=fitPA
  ENDELSE


  newPA[0:2]=MEAN(newPA[0:2])
  if n_elements(PA) GT 6 then begin
     for i=n_elements(fitPA)-3,n_elements(fitPA)-1 do newPA[i]=(fitPA[i]+MEAN(fitPA[n_elements(fitPA)-3:n_elements(fitPA)-1]))/2.
  endif else newPA[n_elements(fitPA)-2:n_elements(fitPA)-1]=MEAN(newPA[n_elements(fitPA)-2:n_elements(fitPA)-1])

  
  cleanup:
   
  IF fixedrings GT n_elements(newPA)-1 then fixedrings=n_elements(newPA)-1

  for j=0,n_elements(PA[*])-1 do errors[j]=MAX([errors[j],ABS(PA[j]-newPA[j]),ddiv])
  
  if keyword_set(debug) then begin
     print,'The estimated  errors before the end'
     print,errors
  ENDIF
 
  tmp=WHERE(FINITE(newPA) EQ 0.)
  IF tmp[0] NE -1 then newPA[WHERE(FINITE(newPA) EQ 0.)]=PAin[WHERE(FINITE(newPA) EQ 0.)]
  IF n_elements(pamin) gt 0 then tmp=WHERE(PA[*]-errors[*] LT pamin) else tmp=-1
  IF tmp[0] NE -1  then begin
     for j=0,n_elements(tmp)-1 do begin
        IF tmp[j] NE n_elements(errors)-1 then begin
           IF keyword_set(debug) then print,errors[tmp[j]]
           errors[tmp[j]]=PA[tmp[j]]-pamin
           IF keyword_set(debug) then print,errors[tmp[j]],PAin[tmp[j]],pamin
        ENDIF
     ENDFOR
  endif
  
  If errors[0] LT errors[1] then errors[0]=errors[1]
  
  if n_elements(errorin) NE 0 then errorin=errors
  PAin=newPA
  if keyword_set(debug) AND n_elements(errorin) NE 0  then begin
     print,'The estimated  errors'
     print,errorin
  ENDIF
  if keyword_set(debug) then begin
     print,'The final fitted output'
     print,PAin
  ENDIF
  IF size(log,/TYPE) EQ 7 then begin
     openu,66,log,/APPEND
     IF arctan EQ 0 then begin
        a=findgen(n_elements(para))
        printf,66,linenumber()+'REGULARISATION_SDIS: We have fitted the dispersion with an arctan profile'
        printf,66,linenumber()+'REGULARISATION_SDIS: With the following coefficients '+STRJOIN('c'+strtrim(string(a,format='(I1)'),2)+'='+strtrim(string(para),2),', ')
       
     ENDIF ELSE BEGIN
        printf,66,linenumber()+'REGULARISATION_SDIS: We have only smoothed the dispersion profile'
     ENDELSE
     close,66
  ENDIF


  
end
