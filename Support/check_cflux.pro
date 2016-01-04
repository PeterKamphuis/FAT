Pro check_cflux,nopoints,tirificfirst,tirificfirstvars,cfluxadjusted,log=log

;+
; NAME:
;       CHECK_CFLUX
;
; PURPOSE:
;       Routine to adjust the cflux such that the number of points in
;       the models is between 0.5E6 and 2.2E6
;
; CATEGORY:
;       Support
; 
; CALLING SEQUENCE:
;       CHECK_CFLUX,nopoints,tirificfirst,tirificfirstvars,cfluxadjusted,log=log
;
;
; INPUTS:
;       nopoints = The number of points in the model
;       tirificfirst = String array with the template of a tirific def file 
;       tirificfirstvars = String array with the variables in the previous template
;
; OPTIONAL INPUTS:
;       -
;
; KEYWORD PARAMETERS:
;       log= general fitting log to document the action the routine has taken
;
; OUTPUTS:
;       cfluxadjusted = a trigger that indicates whether the cflux is
;       adjusted 0 no, 1 yes
;
; OPTIONAL OUTPUTS:
;       -
; 
; PROCEDURES CALLED:
;       LINENUMBER(),STRTRIM(),STRCOMPRESS(),STRING(),STR_SEP()
;
; EXAMPLE:
;
;
; MODIFICATION HISTORY:
;     23-07-2015 added a upper limit of 0.005 and a check that CFLUX
;     is never set to 0.
;     Written by P.Kamphuis 01-01-2015 
;
; NOTE:
;-
  COMPILE_OPT IDL2 
  cfluxadjusted=0.
  for j=0,n_elements(nopoints)-1 do begin
     case (1) of
        nopoints[j] LT 0.5E6: begin
           IF j EQ 0 then begin
              tmppos=where('CFLUX' EQ tirificfirstvars)
              currentcflux=str_sep(strtrim(strcompress(tirificfirst[tmppos]),2),'=')
              IF double(string(double(currentcflux[1])*nopoints[j]/1e6)) NE 0. then $
                 tirificfirst[tmppos]='CFLUX= '+string(double(currentcflux[1])*nopoints[j]/1e6) else $
                    tirificfirst[tmppos]='CFLUX= 1e-5'
           ENDIF ELSE BEGIN
              tmppos=where('CFLUX_'+strtrim(strcompress(fix(j+1)),2) EQ tirificfirstvars)
              currentcflux=str_sep(strtrim(strcompress(tirificfirst[tmppos]),2),'=')
              IF double(string(double(currentcflux[1])*nopoints[j]/1e6)) NE 0. then $
                 tirificfirst[tmppos]=currentcflux[0]+'= '+string(double(currentcflux[1])*nopoints[j]/1e6) else $
                    tirificfirst[tmppos]=currentcflux[0]+'= 1e-5'
           ENDELSE
           if j eq n_elements(nopoints)-1 then begin
              IF size(log,/TYPE) EQ 7 then begin
                 openu,66,log,/APPEND
                 printf,66,linenumber()+"CFLUX is scaled down to"+STRJOIN([string(double(currentcflux[1])*nopoints[0]/1e6),string(double(currentcflux[1])*nopoints[1]/1e6)],' and ')
                 close,66
              endif
           endif
           cfluxadjusted=1.
        end
        nopoints[j] GT 2.2E6: begin
           IF j EQ 0 then begin
              tmppos=where('CFLUX' EQ tirificfirstvars)
              currentcflux=str_sep(strtrim(strcompress(tirificfirst[tmppos]),2),'=')
              IF double(string(double(currentcflux[1])*nopoints[j]/1e6)) LT 0.005 then $
                 tirificfirst[tmppos]='CFLUX= '+string(double(currentcflux[1])*nopoints[j]/1E6) else $
                    tirificfirst[tmppos]='CFLUX= 0.005'
           ENDIF ELSE BEGIN
              tmppos=where('CFLUX_'+strtrim(strcompress(fix(j+1)),2) EQ tirificfirstvars)
              currentcflux=str_sep(strtrim(strcompress(tirificfirst[tmppos]),2),'=')
              IF double(string(double(currentcflux[1])*nopoints[j]/1e6)) LT 0.005 then $
                 tirificfirst[tmppos]=currentcflux[0]+'= '+string(double(currentcflux[1])*nopoints[j]/1E6) else $
                    tirificfirst[tmppos]=currentcflux[0]+'= 0.005'
           ENDELSE
           if j eq n_elements(nopoints)-1 then begin
              IF size(log,/TYPE) EQ 7 then begin
                 openu,66,log,/APPEND
                 printf,66,linenumber()+"CFLUX is scaled up to"+STRJOIN([string(double(currentcflux[1])*nopoints[0]/1e6),string(double(currentcflux[1])*nopoints[1]/1e6)],' and ')
                 close,66
              endif
           endif
        end  
        else:begin
           if j eq n_elements(nopoints)-1 then begin
              IF size(log,/TYPE) EQ 7 then begin
                 openu,66,log,/APPEND
                 printf,66,linenumber()+"CFLUX is ok."
                 close,66
              endif
           endif
        end
     ENDCASE
  ENDFOR
END  
