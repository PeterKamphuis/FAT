Function FAT_FIT,xin,yin,order,RCHI_SQR=rchisqr,newy=newy,CHI_SQR=chisqr,errors=errors,STATUS=status,gdlidl=gdlidl,$
                 fixedrings=fixedrings,Monte_Carlo=mc,mc_iters=mc_iters,mc_errors=mc_errors,mc_maxfails=mc_maxfailedfit,$
                 maximum_deviation=mc_maxdev,nocentral=nocentral,accuracy=accuracy,mc_y_min=ymin,$
                 mc_y_max=ymax,mc_ddiv=ddiv,debug=debug,log=log
  ;+
; NAME:
;       FAT_FIT
;
; PURPOSE:
;       Function to fit polynomials through inverted matrix linear
;       regression method.
;
; CATEGORY:
;       Support
; 
; CALLING SEQUENCE:
;       FAT_FIT,x,y,order,RCHI_SQR=rchisqr,newy=newy,CHI_SQR=chisqr,errors=errors
;
;
; INPUTS:
;       x = array of x-axis data points
;       y = array of dependent y -points
;   order = order of the requested polynomial to fit
;  
; OPTIONAL INPUTS:
;  RCHI_SQR = Reduced Chi Square
;   CHI_SQR = Chi Square  
;    errors = errors on the y data points. Points will be weighted with
;             the inverse squared of the errors.
;      newy = array that contains the calculated y points
;    gdlidl = 1 when running gdl 0 when running idl
;fixedrings = the number of fixedrings that will be reset. If set this will
;             be taken into account into the chi square evaluation as that is the  
;             function returned
;  mc_iters = the iterations to be performed in the evaluation
;  mc_errors= the sigma errors to be used in the evaluation if unset
;             it is assumed to be errors
;
; KEYWORD PARAMETERS:
; /NO_CENTRAL -
; /MONTE_CARLO - Set this to berform a set of iterations with slightly
;                different values in order to level of the variations
;                due to errors. IF set all mc_ inputs are required.
  
; OUTPUTS:
;       p = coefficients of the fitted polynomial 
;
; OPTIONAL OUTPUTS:
;       -
;
; PROCEDURES CALLED:
;       INVERSE(),TOTAL(),TRANSPOSE()
;
; EXAMPLE:
;      
;
; MODIFICATION HISTORY:
;       20-05-2015 P.Kamphuis; Added a condition for real massive
;                              galaxies to always be declining.  
;       Written 29-02-2016 by P.Kamphuis, N. Giese
;
; NOTE:
;     
;-

 
  status=0
  chisqr=1
  rchisqr=1
  y=yin
  x=xin
  IF n_elements(errors) NE 0 then sqerrors=double(errors^2)
  IF n_elements(fixedrings) EQ 0 then fixedrings=0
  p=dblarr(order+1)
  IF n_elements(x) NE n_elements(y) then begin
     print,'FAT_FIT: x should be equal to y'
     help,x,y
     status=1
     stop
     return,!values.f_nan
  ENDIF
  ;build a matrix
  matr=dblarr(n_elements(x),order+1)
  matr[*,0]=1.
  for i=1,order do begin
     matr[*,i]=double(x[*]^i)
  endfor
  ;transpose it
  tmatr=transpose(matr)
  yor=double(y)
  if n_elements(ddiv) GT 0 then begin
     locmax=WHERE(MAX(yor) EQ yor)
     if locmax[n_elements(locmax)-1] GT n_elements(yor)/2. then yor[0:fixedrings]=TOTAL(yor[0:fixedrings])/n_elements(yor[0:fixedrings])
  endif
  y=double(yor)
  newy=dblarr(n_elements(x))
  testy=dblarr(n_elements(x))
  if n_elements(mc_maxfailedfits) EQ 0 then mc_maxfailedfits=1000.
                                ;Since in our mc chain we only need to
                                ;change y we are much better off doing
                                ;all iterations inside here instead of
                                ;building all matrices over and over again.
                                ;calculate coefficients
  IF keyword_set(mc) AND n_elements(mc_iters) GT 0 then begin
     IF n_elements(errors) EQ 0 then begin
        print,'FAT_FIT: you have omitted the errors. MC is not possible.' 
        p=!values.f_nan
        return,p
     ENDIF ELSE BEGIN
        IF  n_elements(mc_errors) EQ 0 then mc_errors=errors
        weight=dblarr(n_elements(sqerrors),n_elements(sqerrors))
        ind=findgen(n_elements(sqerrors))
        weight[ind,ind]=double(1./sqerrors)
        invmatr=la_invert(tmatr#weight#matr)
        iterno=0L
        failedfit=0
        chiarr=dblarr(mc_iters)
        totp=dblarr(n_elements(p))
        prevp=dblarr(n_elements(p))
        totpun=dblarr(n_elements(p))
        prevpun=dblarr(n_elements(p))
        prevp[*]=0.
        prevpun[*]=0.
        coeff=dblarr(mc_iters,order+1)
        checky=dblarr(n_elements(y),mc_iters)
        satisf=0.
        gftim=0.
        check = 0
        WHILE iterno LT mc_iters-1 AND satisf NE 1 DO BEGIN
           IF iterno/10000. EQ fix(iterno/10000.) AND iterno GT 50000 then begin
               IF size(log,/TYPE) EQ 7 then begin
                  openu,66,log,/APPEND
                   IF check EQ 0 then begin
                     printf,66,linenumber()+'FAT_FIT: Using a large number of iterations.'
                     check++
                  ENDIF
                  printf,66,linenumber()+'FAT_FIT: We are currently at iteration number '+strtrim(string(iterno),2)
                  IF iterno GT 75000 then begin
                     printf,66,linenumber()+'FAT_FIT: The current difference is '+strtrim(STRJOIN(div[*]/totp[*]*100.,' '),2)
                     printf,66,linenumber()+'FAT_FIT: And we have satisfied '+strtrim(string(gftim),2)+' times.'
                  ENDIF
                  close,66
               ENDIF ELSE BEGIN
                  IF check EQ 0 then begin
                     print,'FAT_FIT: Using a large number of iterations.'
                     check++
                  ENDIF
                  print,'FAT_FIT: We are currently at iteration number '+strtrim(string(iterno),2)
                  IF iterno GT 75000 then begin
                     print,'FAT_FIT: The current difference is '+strtrim(STRJOIN(div[*]/totp[*]*100.,' '),2)
                     print,'FAT_FIT: And we have satisfied '+strtrim(string(gftim),2)+' times.'
                  ENDIF
               ENDELSE
           ENDIF
           p=tmatr#weight#y#invmatr
           chp=WHERE(FINITE(p))
           IF chp[0] EQ -1 then begin
              failedfit++
              IF failedfit GT mc_maxfailedfit then begin
                 p=!values.f_nan
                 rchisqr=!values.f_nan
                 status=1
                 return,p
              ENDIF
              
           ENDIF ELSE BeGIN
              newy[*]=p[0]
              for i=1,n_elements(p)-1 do newy[*]=newy[*]+p[i]*x[*]^i
 
              IF fixedrings GT 0. and order NE 0 then begin
                 IF n_elements(ddiv) EQ 0 then begin
                    IF fixedrings GT fix(n_elements(newy)/2.)-1 then begin
                       newy[0:fix(n_elements(newy)/2.)-1]=yor[0]
                       newy[fix(n_elements(newy)/2.):fixedrings]=(yor[0]+newy[fix(n_elements(newy)/2.):fixedrings])/2.                    
                    ENDIF else begin
                       newy[0:fixedrings]=yor[0]
                       newy[fixedrings+1]=(newy[fixedrings+1]+newy[fixedrings])/2.
                       newy[fixedrings+2]=(newy[fixedrings+2]*2.+newy[fixedrings+1])/3.
                       newy[fixedrings+3]=(newy[fixedrings+3]*3.+newy[fixedrings+2])/4.
                    ENDELSE
                 ENDIF  ELSE BEGIN
                                ; IF the roation curve is supermassive
                                ; we want it to decline
                    IF (TOTAL(newy[n_elements(newy)-2:n_elements(newy)-1])/2. GT 150 AND newy[n_elements(newy)-1] GT newy[n_elements(newy)-2] AND newy[n_elements(newy)-1] GT newy[n_elements(newy)-3]) OR $
                       (MEAN(newy[0:n_elements(newy)-1]) GT 250.) then begin
                       x=0
                       WHILE newy[x] GT  newy[x+1] AND x LT fix(n_elements(newy)/2) DO x++
                       IF x GT 0 then begin
                          newy[0:x]=newy[x]       
                       ENDIF ELSE BEGIN
                          IF MEAN(newy[1:n_elements(newy)-1]) GT 250. then begin
                             min=MAX(newy)
                             xind=0
                             for x=0,fix(n_elements(newy)/2) do begin
                                IF newy[x] LT min then begin
                                   min=newy[x]
                                   xind=x
                                ENDIF
                             ENDFOR
                             IF xind GT 0 then begin          
                                newy[0:xind]=newy[xind]
                                
                             ENDIF
                          endif
                       ENDELSE
                    ENDIF
                
   

                    
                                ;If we have a declining rotation curve
                                ;we want the outer part to be reset to
                                ;flat.
         
                    if locmax[n_elements(locmax)-1] GT fixedrings then begin
                       IF locmax[n_elements(locmax)-1] LT fix(n_elements(yor)/2.) then newy[0:locmax[n_elements(locmax)-1]]=newy[locmax[n_elements(locmax)-1]]
                    ENDIF
                    IF keyword_set(nocentral) AND n_elements(xin) GT 15. then newy[0:fixedrings]=newy[fixedrings+1]
                                ; Also if it is supermassive we want
                                ; it flat 
                    
                 ENDELSE
              ENDIF
              
              IF n_elements(errors) GT 0 then chisqr=TOTAL((newy - y)^2/sqerrors) else chisqr=TOTAL((newy - y)^2)
                                ;If underdetermined do not penalize
              
              IF chisqr LT (n_elements(fitPA)-order) then goto,skippenalize
             
                                ;Let's see if any of the
                                ;changes per ring are outside
                                ;the maximum.
                                ;If so penalize the chi square
              IF n_elements(maxdev) then begin
                 deltaslope=0.
                 IF n_elements(ddiv) GT 0 then begin
                    for j=0,fix(n_elements(newy)/2.) do begin
                       deltaslope=(newy[j]-newy[j+1])-(newy[j+1]-newy[j+2])
                       IF ABS(deltaslope) GT maxdev then begin
                          chisqr=chisqr*((ABS(deltaslope))/(maxdev))^5
                       ENDIF
                    endfor
                 ENDIF ELSE BEGIN
                    for j=n_elements(newy)-2,fix(n_elements(newy)/2.),-1 do begin
                       deltaslope=(newy[j+1]-newy[j])-(newy[j]-newy[j-1])
                       deltaslope2=(y[j+1]-y[j])-(y[j]-y[j-1])
                       IF (ABS(deltaslope/deltaslope2) GT 2.5 OR ABS(deltaslope/deltaslope2) LT 1./2.5) AND deltaslope GT maxdev/accuracy then begin
                          IF j EQ n_elements(newy)-2 then begin
                             chisqr=chisqr*((ABS(deltaslope))/(maxdev/accuracy))^2*2
                          ENDIF ELSE BEGIN
                             chisqr=chisqr*((ABS(deltaslope))/(maxdev/accuracy))^2
                          ENDELSE
                       ENDIF
                    endfor
                 ENDELSE
              ENDIF
                                ;If the fitted polynomial exceeds the
                                ;set minimum and maximum we penalize it
              if n_elements(ymin) GT 0 and n_elements(ymax) GT 0 then begin
                 toolarge=WHERE(newy GT ymax)
                 toosmall=WHERE(newy LT ymin)
                 IF n_elements(ddiv) NE 0 and not keyword_set(nocentral) then begin
                    IF toosmall[0] NE -1 AND toosmall[n_elements(toosmall)-1] EQ n_elements(newy)-1 then begin
                       endsmall=n_elements(toosmall)-2
                       IF newy[n_elements(newy)-2] LT 0. then begin
                          IF keyword_set(debug) then begin
                             print,'at order',order
                             print,'we get neg values for ',newy[n_elements(newy)-2],n_elements(newy)-2
                          ENDIF
                          chisqr=chisqr^2
                       ENDIF
                    ENDIF else endsmall=n_elements(toosmall)-1
                 ENDIF ELSE endsmall=n_elements(toosmall)-1
                 IF toolarge[0] NE -1 then begin
                    for i=0,n_elements(toolarge)-1 do chisqr=chisqr*SQRT(ABS((newy[toolarge[i]]-ymax)))+1.
                 ENDIF

                 IF toosmall[0] NE -1 then begin
                    
                    for i=0,endsmall do  chisqr=chisqr*SQRT(ABS((ymin-newy[toosmall[i]])))+1.
                    
                 ENDIF
              ENDIF
;we dislike declining rotation curves so let's check the last 1/4.
              IF n_elements(ddiv) NE 0 then begin
                 for i=1,fix(n_elements(newy)/4.) do begin
                    change=newY[i-1]-newY[i]
                    IF change LT 0. AND ABS(change) GT ddiv then chisqr=chisqr*ABS(change)
                 endfor
              ENDIF

skippenalize:
     
              chiarr[iterno]=chisqr
              coeff[iterno,*]=p
              totp[*]=TOTAL(chiarr[0:iterno]#coeff[0:iterno,*],1)/TOTAL(chiarr[0:iterno])
            
              div=ABS(totp[*]-prevp[*])
              
              tmp=WHERE(ABS(div[*]/totp[*]*100.) LT 1e-2)
              
              
              IF n_elements(tmp) GE n_elements(totp)/2. AND TOTAL(ABS(div[*]/totp[*]*100.))/n_elements(totp) LT 1e-2 then begin
                 gftim++
               
                 IF gftim gt 1000 AND iterno GT 10000 then satisf=1
              ENDIF
             ; if iterno GT 10000 then stop
              prevp[*]=totp[*]
            
              iterno++
              checky[*,iterno]=newy[*]
           ENDELSE
           ini=randomu(seed,/Double)
           y[*]=yor[*]+(randomu(seed,/Double,n_elements(y))-0.5)*mc_errors[*]
        ENDWHILE
     ;   print,'How many iterations?',iterno
        p[*]=totp[*]        
      
        for j=0,n_elements(y)-1 do mc_errors[j]=STDDEV(checky[j,0:iterno-1])
        newy[*]=p[0]
        for i=1,n_elements(p)-1 do newy[*]=newy[*]+p[i]*x[*]^i
           
        IF fixedrings GT 0. and order NE 0 then begin
           IF n_elements(ddiv) EQ 0 then begin
              IF fixedrings GT fix(n_elements(newy)/2.)-1 then begin
                 newy[0:fix(n_elements(newy)/2.)-1]=yor[0]
                 newy[fix(n_elements(newy)/2.):fixedrings]=(yor[0]+newy[fix(n_elements(newy)/2.):fixedrings])/2.                    
              ENDIF else begin
                 newy[0:fixedrings]=yor[0]
                 newy[fixedrings+1]=(newy[fixedrings+1]+newy[fixedrings])/2.
                 newy[fixedrings+2]=(newy[fixedrings+2]*2.+newy[fixedrings+1])/3.
                 newy[fixedrings+3]=(newy[fixedrings+3]*3.+newy[fixedrings+2])/4.
              ENDELSE
           ENDIF  ELSE BEGIN
                                ; IF the roation curve is supermassive
                                ; we want it to decline
              IF (TOTAL(newy[n_elements(newy)-2:n_elements(newy)-1])/2. GT 150 AND newy[n_elements(newy)-1] GT newy[n_elements(newy)-2] AND newy[n_elements(newy)-1] GT newy[n_elements(newy)-3]) OR $
                 (MEAN(newy[0:n_elements(newy)-1]) GT 250.) then begin
                 x=0
                 WHILE newy[x] GT  newy[x+1] AND x LT fix(n_elements(newy)/2) DO x++
                 IF x GT 0 then begin
                    newy[0:x]=newy[x]       
                 ENDIF ELSE BEGIN
                    IF MEAN(newy[1:n_elements(newy)-1]) GT 250. then begin
                       min=MAX(newy)
                       xind=0
                       for x=0,fix(n_elements(newy)/2)-1 do begin
                          IF newy[x] LT min then begin
                             min=newy[x]
                             xind=x
                          ENDIF
                       ENDFOR
                       IF xind GT 0 then begin          
                          newy[0:xind]=newy[xind]
                          
                       ENDIF
                    ENDIF
                 ENDELSE
              ENDIF
              
                                ;If we have a declining rotation curve
                                ;we want the outer part to be reset to
                                ;flat.
              if locmax[n_elements(locmax)-1] GT fixedrings then begin
                 IF locmax[n_elements(locmax)-1] LT fix(n_elements(yor)/2.) then newy[0:locmax[n_elements(locmax)-1]]=newy[locmax[n_elements(locmax)-1]]
              ENDIF
              IF keyword_set(nocentral) AND n_elements(xin) GT 15. then newy[0:fixedrings]=newy[fixedrings+1]
;              if locmax[n_elements(locmax)-1] GT n_elements(yor)/2. then newy[0:fixedrings]=TOTAL(y[0:fixedrings])/n_elements(y[0:fixedrings])
 
           ENDELSE
        ENDIF
          
           ;IF n_elements(errors) GT 0 then chisqr=TOTAL((newy - y)^2/(sqerrors)) else chisqr=TOTAL((newy - y)^2)
       
        ;rchisqr=chisqr/(n_elements(y)-order)
        chisqr=TOTAL(chiarr[0:iterno-1])/n_elements(chiarr[0:iterno-1])
        rchisqr=chisqr/(n_elements(y)-order)
     ENDELSE
     chiarr=0.
     coeff=0.
     checky=0.
  ENDIF ELSE begin  
     if n_elements(errors) GT 0 then begin
        weight=dblarr(n_elements(sqerrors),n_elements(sqerrors))
        ind=findgen(n_elements(sqerrors))
        weight[ind,ind]=double(1./sqerrors)
        IF order GT 0 then p=tmatr#weight#y#la_invert(tmatr#weight#matr) else begin
           inv=tmatr#weight#matr
           tmpmatr=tmatr#weight#y
           p=tmpmatr/inv
        ENDELSE
     endif else begin
        IF order GT 0 then p=tmatr#y#la_invert(tmatr#matr) else begin
           inv=tmatr#matr
           tmpmatr=tmatr#y
           p=tmpmatr/inv
        ENDELSE
     ENDELSE
     
                                ;calculate new fit
     newy=dblarr(n_elements(x))
     newy[*]=p[0]
     for i=1,n_elements(p)-1 do newy[*]=newy[*]+p[i]*x[*]^i
     IF fixedrings GT 0. and order NE 0 then begin
        IF fixedrings GT fix(n_elements(newy)/2.)-1 then begin
           newy[0:fix(n_elements(newy)/2.)-1]=yor[0]
           newy[fix(n_elements(newy)/2.):fixedrings]=(yor[0]+newy[fix(n_elements(newy)/2.):fixedrings])/2.                    
        ENDIF else begin
           newy[0:fixedrings]=yor[0]
           newy[fixedrings+1]=(newy[fixedrings+1]+newy[fixedrings])/2.
           newy[fixedrings+2]=(newy[fixedrings+2]*2.+newy[fixedrings+1])/3.
           newy[fixedrings+3]=(newy[fixedrings+3]*3.+newy[fixedrings+2])/4.
        ENDELSE
     ENDIF
                                ;calculate chi-squares    
     IF n_elements(errors) GT 0 then chisqr=TOTAL((newy - y)^2/sqerrors) else chisqr=TOTAL((newy - y)^2)
     rchisqr=chisqr/(n_elements(y)-order)
     check=WHERE(FINITE(p))
     IF check[0] EQ -1 then status=1
  ENDELSE
  finishearly:
  return,p
end
