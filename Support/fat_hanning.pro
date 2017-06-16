Function fat_hanning,SBRin

;+
; NAME:
;       FAT_HANNING
;
; PURPOSE:
;       Routine to apply a hanning smoothing to the SBR Profile
;
; CATEGORY:
;       Support
; 
; CALLING SEQUENCE:
;       Result = FAT_HANNING(Profile)
;
;
; INPUTS:
;      SBRin = the profile to be smoothed 
;
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
;       Written 16-06-2017 by P.Kamphuis, S. Kurapati
;
; NOTE:
;     
;-
  COMPILE_OPT IDL2
  DEFSYSV, '!GDL', EXISTS = gdlidl ;is 1 when running GDL
 ;  goto,skipcatch
  CATCH,Error_status  
  IF  Error_status NE 0. THEN BEGIN
    
     print, ' '
     print, 'Oops the following went wrong in FAT_HANNING:'
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
  SBR=SBRin
  SBRout=dblarr(n_elements(SBR))
  case 1 of
     n_elements(SBR) LT 3:return,SBR
     n_elements(SBR) LT 15:points=5
     else:points=7
  endcase
  window=dblarr(points)
  xaxis=findgen(points)
  window=0.5*(1-Cos(!pi*2.*xaxis[*]/(points-1)))
  window=window[*]/TOTAL(window)
  for i=fix((points-1)/2.),n_elements(SBR)-(fix((points-1)/2.)+1) do begin
     SBRout[i]=TOTAL(SBR[i-fix((points-1)/2.):i+fix((points-1)/2.)]*window)
  endfor
                                ;We'll just use a truncated edge in the center
  for i=0,fix((points-1)/2.)-1 do begin
     SBRout[i]=TOTAL(SBR[0:i+fix((points-1)/2.)]*window[fix((points-1)/2.)-i:points-1])
  endfor
                                ;at the outer edge we can use
                                ;zer0-padding
  for i=n_elements(SBR)-fix((points-1)/2.),n_elements(SBR)-1 do begin
    SBRout[i]=TOTAL([SBR[i-fix((points-1)/2.):n_elements(SBR)-1],replicate(0.,i-(n_elements(SBR)-(fix((points-1)/2.)+1)))]*window)
  endfor
  return,SBRout
end
