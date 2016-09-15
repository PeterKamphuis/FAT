Pro revised_regularisation_rot,PAin,SBRin,RADIIin,error=errorin,fixedrings=fixedringsin,REVERSE=rev,NOWEIGHT=noweight,Difference=DDivin,cutoff=cutoffin,arctan=arctanin,debug=debug,nocentral=nocentral,order=order,max_deviation=maxdevin,max_par=PAmaxin,min_par=PAminin,accuracy=accuracy,extending=extending,gdlidl=gdlidl,log=log

;+
; NAME:
;       REVISED_REGULARISTION_ROT
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
  IF KEYWORD_SET(NOCENTRAL) then begin 
     SBR=SBRin[1:n_elements(SBRin)-1]
     RADII=RADIIin[1:n_elements(SBRin)-1]
     PA0=PAin[0]
     PA=PAin[1:n_elements(SBRin)-1]
     IF n_elements(cutoffin) GT 0 then cutoff=cutoffin[1:n_elements(SBRin)-1]
  ENDIF ELSE BEGIN
     SBR=SBRin
     RADII=RADIIin
     PA=PAin
     IF n_elements(cutoffin) GT 0 then cutoff=cutoffin[0:n_elements(SBRin)-1]
  ENDELSE
  PA=REVERSE(PA)
  SBR=REVERSE(SBR)
  RADII=REVERSE(RADII)
  cutoff=REVERSE(cutoff)

  IF n_elements(cutoffin) GT 0 then cutoff=Reverse(cutoff)
  If n_elements(fixedringsin) GT 0 then fixedrings=fixedringsin else fixedrings=0
  If n_elements(DDivin) GT 0 then DDiv=DDivin
  If n_elements(errorin) GT 0 then error=errorin
  If n_elements(cutoffin) GT 0 then cutoff=cutoffin
                                ;arctan can have the values of 0: for a polynomial fit
                                ;1: for a smoothing of PA or INCL  (default)
                                ;2: for a smoothing of VROT
  If n_elements(arctanin) GT 0 then arctan=arctanin else arctan=1
  If n_elements(maxdevin) GT 0 then  maxdev=maxdevin
  If n_elements(PAmaxin) GT 0 then PAmax=double(PAmaxin)
  If n_elements(PAminin) GT 0 then PAmin=double(PAminin)
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
  IF n_elements(PA) GT 15 then decline=0.9 else decline=0.8
  PAsmooth=dblarr(n_elements(PA[*]))
  IF PA[0] EQ 0. then PAsmooth[0]=PA[0] else PAsmooth[0]=(PA[0]+PA[1])/2.                               
  for j=1,n_elements(PA[*])-2 do begin         
     PAsmooth[j]=(PA[j-1]+PA[j]+PA[j+1])/3
     WHILE PAsmooth[j] LT PAsmooth[j-1]*decline do begin
        IF keyword_set(debug) then begin
           print,PA[j],PA[j-1],'increasing'
        ENDIF
        IF PAsmooth[j] LT 0. then begin
           IF PAsmooth[j-1] GT 0. then PAsmooth[j]=PAsmooth[j-1] else begin
              IF n_elements(PAmin) GT 0 then PAsmooth[j]=PAmin*2. else begin
                 tmpind=WHERE(PAsmooth GT 0.)
                 if tmpind[0] NE -1 then PAsmooth[j]=MEAN(PAsmooth[tmpind]) else PAsmooth[j]=1.
              ENDELSE
           ENDELSE
           
        endif else PAsmooth[j]=PAsmooth[j]*1.05
     ENDWHILE
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
     ratio[0:tmp]=1
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
                                ;will set them to ddivadittionally if
                                ;we have any errors less than ddiv
                                ;they will also be set to ddiv
     if n_elements(fixedrings) GT 0 then begin
       
           errors[0:fixedrings-1]=ddiv
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
    
     for j=1,fix(n_elements(PA[*])/2.) do begin
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
              sign[x]=PA[j+x]-PA[j+x-1]
           endfor
           tmp=WHERE(sign LT 0)
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
     IF PA[n_elements(PA)-1] EQ 0 then errors[n_elements(PA)-1]=ddiv/4.
     for j=fixedrings+1,0,-1 do begin
         IF errors[j+1] EQ 0 then errors[j+1] = errors[j]/2.
        WHILE errors[j] GT errors[j+1] do begin
           errors[j+1]=errors[j+1]*fact
        ENDWHILE
     endfor
     
  ENDELSE
  IF keyword_set(debug) then begin
     print,'The final errors are:'
     print,errors
  ENDIF
  
;IF we have maxdev then first we are going to do a small check on
;ridiculous outliers.
  
  IF n_elements(maxdev) GT 0. then begin
     lastmod=1
     start=n_elements(PA[*])-2
     for i=n_elements(PA[*])-4,1,-1 do begin
        diff1=ABS(PAin[i]-PAin[i-1])
        diff2=ABS(PAin[i]-PAin[i+1])
        diff3=ABS(PAin[i+1]-PAin[i-1])
        IF diff1 GT maxdev AND  diff2 GT maxdev AND diff3 LT (diff1+diff2)/2. then begin
           if keyword_set(debug) then begin
              print,"These are the differences from ring to ring +the value +maxdev"
              print,diff1,diff2,diff3,i,PA[i],maxdev
           ENDIF
           PA[i]=(PAin[i-1]+PAin[i+1])/2
           lastmod=i
        ENDIF ELSE start=i
     endfor
    
     IF keyword_set(debug) then begin
        print,start,lastmod,n_elements(PA),PA,fixedrings
     ENDIF
                                ;if we have these extremes untill the last point just take the average
     if lastmod EQ 1 AND (start GT 3 OR start GE fixedrings AND  start+1 GE 1) then begin
        IF keyword_set(debug) then begin
           print,'big jigs till the end adjusting before'
           print, PA[start+1:n_elements(PA)-1]
        ENDIF
        PA[0:start+1]=TOTAL(PAin[0:start+1])/n_elements(PAin[0:start+1])
        IF keyword_set(debug) then begin
           print,'big jigs till the end adjusting'
           print, PA[0:start+1]
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
  ENDIF
                                ;If all reliable rings are fixed we
                                ;can skip the fitting of polynomials
                                ;and just want to use a straight line
 
  newPA=PA
  newPA[*]=0.
  endorder=n_elements(PA[*])
        
  IF keyword_set(debug) then begin
     print,'this is the fixedrings',fixedrings
     print,'This is our endorder',endorder
  ENDIF
  IF endorder LT 2 then endorder=2
  IF endorder GT 5 then endorder=5
  maxendorder=endorder
  fitPAoriginal=PA[*]
  fiterrors=errors[*]
  shifterrorsor=errors[*]
  newendorder:
  IF endorder GT maxendorder then endorder=maxendorder
  fitPA=fitPAoriginal
     
  Chi=dblarr(endorder-1)
  mcerrors=dblarr(n_elements(PA[*]),endorder-1)
  finalcoeff=dblarr(6,endorder-1)
  for order=endorder,2,-1 do begin
     shifterrors=shifterrorsor
     tmp=dblarr(1)
     IF keyword_set(debug) then begin
        print,'starting loop for order',order
     ENDIF
    
     fitfailed=0.
     newPAcoeff=FAT_FIT(RADII,fitPA,order,RCHI_SQR=tmp,errors=fiterrors,STATUS=fitstat,/Monte_Carlo,mc_iters=150000.,$
                           mc_errors=shifterrors,maximum_deviation=maxdev,accuracy=accuracy,fixedrings=fixedrings,$
                           mc_y_min=pamin,mc_y_max=pamax,mc_ddiv=ddiv,nocentral=nocentral,mc_maxfails=100.,log=log)
           
     newPA[*]=fitPA[*]   
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
     order=-1
     IF keyword_set(debug) then begin
        print,'We return the smoothed profile as all Chis are infinite'
     ENDIF
     attempts++
     fitPAoriginal=PAsmooth[*]
     shifterrorsor=errors[*]
     IF attempts LT 2 then goto,newendorder
  ENDIF
       
  maxchi=MAX(Chi,MIN=minchi)
  IF keyword_set(debug) then begin
     print,'Max min',maxchi,minchi
  ENDIF

  case 1 of 
     minchi EQ 1E9: begin
        IF erroradjusted NE 1 then begin
           order=-1
           IF keyword_set(debug) then begin
              print,'We will return the smoothed profile because the error is not adjusted and minchi is GT 1e9'
           ENDIF
           attempts++
           fitPAoriginal=PAsmooth[*]
           shifterrorsor=errors[*]
           goto,newendorder
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
  
     
                                ;Let's check a 0th order
                                ;polynomial if it is not the
                                ;rotation curve
     
    
           
  order=WHERE(Chi EQ minchi)+2
  tmp=order
  order=intarr(1)
  order=tmp[0]
  newPAcoeff=dblarr(order)
  newPAcoeff=finalcoeff[0:order,order-2]
    
  IF fitstat LT 1 then begin
     finorder=order
     newPA[*]=newPAcoeff[0]
     for i=1,n_elements(newPAcoeff)-1 do begin
        newPA[*]= newPA[*]+newPAcoeff[i]*RADII[*]^i 
     endfor
     locmax=WHERE(MAX(PA) EQ PA)
     if locmax[n_elements(locmax)-1] GT n_elements(PA)/2. then newPA[0:fixedrings]=TOTAL(PA[0:fixedrings])/n_elements(PA[0:fixedrings])
     endch=floor(n_elements(PA)/3.)
     IF endch GT 3 then endch=3
     decline=0
     for ch=1,endch do begin
        if newPA[ch-1] lt newPA[ch]*0.9 then decline=ch else break
     endfor
     if decline gt 0. then newPA[0:ch-1]=newPA[ch]*0.9    
     errors[*]=mcerrors[*,order-2] 
  ENDIF Else begin
     IF attempts GE 1 then begin
        finorder=!values.f_nan
        newPA[*]=PAsmooth[*]
        errors[*]=fiterrors[*]
        arctan=1
     ENDIF else begin
        attempts++
        fitPAoriginal=PAsmooth[*]
        shifterrorsor=errors[*]
        goto,newendorder
     ENDELSE
  ENDELSE
  if keyword_set(debug) then begin
     print,'We find this as the new profile'
     print,newPA[*]
     print,'With the order '+string(finorder)
  endif
  erroradjusted=0
    
  
  cleanup:
  IF fixedrings GT n_elements(newPA)-1 then fixedrings=n_elements(newPA)-1
  IF fixedrings GT 0. and order NE 0 then begin
     locmax=WHERE(MAX(REVERSE(PAin)) EQ REVERSE(PAin))
     if locmax[n_elements(locmax)-1] GT fixedrings then newPA[0:fixedrings]=TOTAL(newPA[0:fixedrings])/n_elements(newPA[0:fixedrings])
  ENDIF
 
  
 
  for j=0,n_elements(PA[*])-1 do errors[j]=MAX([errors[j],ABS(PA[j]-newPA[j]),ddiv])
  
  if keyword_set(debug) then begin
     print,'The estimated  errors before the end'
     print,errors
  ENDIF
 
  tmp=WHERE(FINITE(newPA) EQ 0.)
  IF tmp[0] NE -1 then newPA[WHERE(FINITE(newPA) EQ 0.)]=PAin[WHERE(FINITE(newPA) EQ 0.)]
  IF n_elements(pamin) gt 0 then tmp=WHERE(PA[*]-errors[*] LT pamin) else tmp=-1
  IF tmp[0] NE -1  then begin
     errors[tmp]=PAin[tmp]-pamin
  endif
  IF SBRin[n_elements(PAin[*])-1] LT cutoff[n_elements(PAin[*])-1] then errors[0]=MAX(errors[*])
  If errors[0] LT errors[1] then errors[0]=errors[1]
  IF n_elements(pamax) gt 0 then tmp=WHERE(newPA[*]+errors[*] GT pamax) else tmp=-1
  IF tmp[0] NE -1 then errors[tmp]=pamax-PA[tmp]
  IF n_elements(pamin) gt 0 and n_elements(pamax) gt 0 then tmp=WHERE(errors[*] LT 0) else tmp=-1
  IF tmp[0] NE -1 AND n_elements(errorin) NE 0 then errorin[tmp]=pamax-pamin
  IF n_elements(errorin) NE 0 then errorin=REVERSE(errors)
  
  order=finorder
 
  
  PAin=REVERSE(newPA)
  IF KEYWORD_SET(NOCENTRAL) then begin
     PAin=[PA0,PAin]
     IF n_elements(errorin) NE 0 then errorin=[errorin[0],errorin]
  ENDIF
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
        a=findgen(n_elements(newPAcoeff))
        printf,66,linenumber()+'REVISED_REGULARISATION_ROT: We have fitted the Rotation curve with a polynomial of order '+strtrim(string(order),2)
        printf,66,linenumber()+'REVISED_REGULARISATION_ROT: With the following coefficients '+STRJOIN('c'+strtrim(string(a,format='(I1)'),2)+'='+strtrim(string(newPAcoeff),2),', ')
        printf,66,linenumber()+'REVISED_REGULARISATION_ROT: This is fitted in Reverse'
        if locmax[n_elements(locmax)-1] GT fixedrings then printf,66,linenumber()+'REVISED_REGULARISATION_ROT: We have flattened the last '+strtrim(string(fixedrings),2)+' rings'
     ENDIF ELSE BEGIN
        printf,66,linenumber()+'REVISED_REGULARISATION_ROT: We have only smoothed the rotation curve'
     ENDELSE
     close,66
  ENDIF


  
end
