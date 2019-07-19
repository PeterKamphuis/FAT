Function fat_savgol,SBRin,Radin,rings=rings,step=step

;+
; NAME:
;       FAT_SAVGOL
;
; PURPOSE:
;       Routine to apply a SAVITSKY-GOLAY smoothing to the SBR Profile
;
; CATEGORY:
;       Support
; 
; CALLING SEQUENCE:
;       Result = FAT_SAVGOL(Profile)
;
;
; INPUTS:
;      SBRin = the profile to be smoothed 
;      rings = the amount of rings that are valid
; OPTIONAL INPUTS:
;       - 
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;       Result = the smoothed profile 
;
; OPTIONAL OUTPUTS:
;       -
;
; PROCEDURES CALLED:
;       
;
; EXAMPLE:
;      
;
; MODIFICATION HISTORY:
;       Written 19-07-2017 by P.Kamphuis
;
; NOTE:
;     
;-
  COMPILE_OPT IDL2
  DEFSYSV, '!GDL', EXISTS = gdlidl ;is 1 when running GDL
  goto,skipcatch
  CATCH,Error_status
  IF  Error_status NE 0. THEN BEGIN
    
     print, ' '
     print, 'Oops the following went wrong in FAT_SAVGOL:'
     print,'-------------------------------------------'
     help,/Last_Message, Output = theerrmess
                                ;This gives trace back information in
                                ;IDL but unfortunately not in GDL
     for j=0,n_elements(theerrmess)-1 do print,theerrmess[j]
     print,'-------------------------------------------'
     if gdlidl then begin
        print, ' '
        print,'Unfortunately GDL does not provide trace back information in help'
        print,'In order to fix this error please open FAT.pro and '
        print,'uncomment line 51 (goto,skipcatch).'
        print, ' '
        print,'Then reproduce the error with trace back information '
        print, 'and submit an issue at:'
        print, ' '
        print,'       https://github.com/PeterKamphuis/FAT/issues'
        print, ' '
        print,'If the error occured while fitting a galaxy, please attach you fitting'
        print,'log as well.' 
     endif else begin
        print, ' '
        print,'If this message completely baffles you then please submit the last 20 lines of output as a bug report to:'
        print, ' '
        print,'       https://github.com/PeterKamphuis/FAT/issues'
        print, ' '
        print,'If the error occured while fitting a galaxy, please attach you fitting'
        print,'log as well.' 
     endelse
     stop
  ENDIF
  skipcatch:
  if n_elements(step) EQ 0 then step=0
  ;We will initially ignore the central point.
 
                                ;If we provide a number of valid rings
                                ;then we want to set all rings outside
                                ;that to 0.
;  IF n_elements(rings) GT 0 then begin
;     IF rings[0] LT n_elements(SBR) then SBR[rings[0]-1:n_elements(SBR)-1]=0.
;  ENDIF
  

  case 1 of
     n_elements(SBRin) LE 4:begin
        SBR = SBRin
        SBR[0]=(SBRin[0]+SBRin[1]*2.)/3.
        
        case  n_elements(SBRin) of
           2:begin
              SBR[1]=SBRin[1]
           end
           3:begin
              SBR[1]=(SBRin[0]+SBRin[1]*2+SBRin[2]*2.)/5.
              SBR[2]=(SBRin[1]+SBRin[2]*2)/3.
           end
           4:begin
              SBR[1]=(SBRin[0]+SBRin[1]*2+SBRin[2]*2.)/5.
              SBR[2]=(SBRin[1]+SBRin[2]+SBRin[3])/3.
              SBR[3]=(SBRin[2]+SBRin[3]*2.)/3.
           end
           else:print,'This should never happen'
        endcase
        return,SBR
     end
     n_elements(SBRin) LT 10: begin
        points=3
        order= 1
     end
     n_elements(SBRin) LT 15: begin
        points=7
        order=2
     end
     n_elements(SBRin) LT 20: begin
        points=9
        order=3
     end
     else: begin
        points=11
        order=4
     end
  endcase
  
  SBR=SBRin[0:n_elements(SBRin)-1]
  RAD=radin[0:n_elements(SBRin)-1]
  extRAD=dblarr(n_elements(SBR)+fix(points/2.+1))
  extSBR=dblarr(n_elements(SBR)+fix(points/2.+1))
  extRAD[0:n_elements(SBR)-1]=RAD[*]
;As we want to make sure the SBR
                                ;profile tapers of nicely we pad we
                                ;1E-16 the outer rings
  extSBR[*]=1E-8
  extSBR[0:n_elements(SBR)-1]=SBR[*]
;If the inner to points are more then 1.5* point 3 than make them half
  IF extSBR[0] GT extSBR[2]*1.5 then extSBR[0]=extSBR[0]/2.
  IF extSBR[1] GT extSBR[2]*1.5 then extSBR[1]=extSBR[1]/2.

  
  for i=0,fix(points/2.) do begin
     extRAD[n_elements(SBR)+i]=RAD[n_elements(SBR)-1]+(i+1)*(RAD[n_elements(SBR)-1]-RAD[n_elements(SBR)-2])
  endfor
  IF n_elements(rings) GT 0 then begin
     IF rings[0]+1 LT n_elements(SBR) then SBR[rings[0]+1:n_elements(SBRout)-1]=1e-16
  ENDIF
  SBRout=dblarr(n_elements(SBR))
  for i=0,n_elements(SBR)-1 do begin
     case 1 of
        i LT points/2.-1: begin
           x=extRAD[0:points-1]
           y=extSBR[0:points-1]
           retr_index=i
          
        end
   
        else: begin
           x=extRAD[i-(points-1)/2:i+(points-1)/2]
           y=extSBR[i-(points-1)/2:i+(points-1)/2]
           retr_index=fix(points/2.+0.5)-1
        end
     endcase
     dummy = fat_fit(x,y,order,newy=ourpoint)
   
     SBRout[i]=ourpoint[retr_index]
  endfor
  
  IF n_elements(rings) GT 0 then begin
     IF rings[0]+1 LT n_elements(SBR) then SBRout[rings[0]+1:n_elements(SBRout)-1]=1e-16
  ENDIF
  tmp=WHERE(SBRout LT 1e-8)
  if tmp[0] NE -1 then SBRout[tmp]=1e-16

                                ;Finally we need to add the central
                                ;point again which we will merely
                                ;extrapolate from the profile
 
  neg_ind=WHERE(SBRout LT 0.)
  if neg_ind[0] NE -1 then SBRout[neg_ind]=0.
  WHILE n_elements(SBRout) NE n_elements(SBRin) do begin
     IF n_elements(SBRout) LT n_elements(SBRin) then SBROut=[SBRout[0],SBRout]
     IF n_elements(SBRout) GT n_elements(SBRin) then SBROut=SBRout[0:n_elements(SBRin)-2]
  ENDWHILE
  return,SBRout
end
