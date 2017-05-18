Pro check_cflux,nopoints,norings,tirificfirst,tirificfirstvars,cfluxadjusted,log=log

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
;     06-03-2017 P. Kamphuis: Added a condition that allows for higher
;     number of point sources for larger models. Additionally also
;     check that nopoints is finite.  
;     23-07-2015 added a upper limit of 0.005 and a check that CFLUX
;     is never set to 0.
;     Written by P.Kamphuis 01-01-2015 
;
; NOTE:
;-
  COMPILE_OPT IDL2 
  cfluxadjusted=0.
 
  for j=0,n_elements(nopoints)-1 do begin
 
     if FINITE(nopoints[j]) EQ 0 then begin
        IF size(log,/TYPE) EQ 7 then begin
           openu,66,log,/APPEND
           printf,66,linenumber()+"CHECK_CFLUX: We had an infinite number of points somehow. We stop the code"
           close,66
        endif
        stop
     ENDIF ELSE BEGIN
        IF norings[0] LT 15 then begin
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
                    currentcflux1=str_sep(strtrim(strcompress(tirificfirst[tmppos]),2),'=')
                    IF double(string(double(currentcflux1[1])*nopoints[j]/1e6)) NE 0. then $
                       tirificfirst[tmppos]=currentcflux1[0]+'= '+string(double(currentcflux1[1])*nopoints[j]/1e6) else $
                          tirificfirst[tmppos]=currentcflux1[0]+'= 1e-5'
                 ENDELSE
                 if j eq n_elements(nopoints)-1 then begin
                    IF n_elements(currentcflux) EQ 0 then begin
                       tmppos=where('CFLUX' EQ tirificfirstvars)
                       currentcflux=str_sep(strtrim(strcompress(tirificfirst[tmppos]),2),'=')
                    ENDIF
                       
                    IF size(log,/TYPE) EQ 7 then begin
                       openu,66,log,/APPEND
                       printf,66,linenumber()+"CHECK_CFLUX: CFLUX is scaled down to"+STRJOIN([string(double(currentcflux[1])*nopoints[0]/1e6),string(double(currentcflux1[1])*nopoints[1]/1e6)],' and ')
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
                    currentcflux1=str_sep(strtrim(strcompress(tirificfirst[tmppos]),2),'=')
                    IF double(string(double(currentcflux1[1])*nopoints[j]/1e6)) LT 0.005 then $
                       tirificfirst[tmppos]=currentcflux1[0]+'= '+string(double(currentcflux1[1])*nopoints[j]/1E6) else $
                          tirificfirst[tmppos]=currentcflux1[0]+'= 0.005'
                 ENDELSE
                 if j eq n_elements(nopoints)-1 then begin
                    IF n_elements(currentcflux) EQ 0 then begin
                       tmppos=where('CFLUX' EQ tirificfirstvars)
                       currentcflux=str_sep(strtrim(strcompress(tirificfirst[tmppos]),2),'=')
                    ENDIF
                    IF size(log,/TYPE) EQ 7 then begin
                       openu,66,log,/APPEND
                       printf,66,linenumber()+"CHECK_CFLUX: CFLUX is scaled up to"+STRJOIN([string(double(currentcflux[1])*nopoints[0]/1e6),string(double(currentcflux1[1])*nopoints[1]/1e6)],' and ')
                       close,66
                    endif
                 endif
              end  
              else:begin
                 if j eq n_elements(nopoints)-1 then begin
                    IF size(log,/TYPE) EQ 7 then begin
                       openu,66,log,/APPEND
                       printf,66,linenumber()+"CHECK_CFLUX: CFLUX is ok."
                       close,66
                    endif
                 endif
              end
           ENDCASE
        ENDIF ELSE BEGIN
           factor=(norings[0]/15.)^1.5
             case (1) of
              nopoints[j] LT 0.5E6*factor: begin
                 IF j EQ 0 then begin
                    tmppos=where('CFLUX' EQ tirificfirstvars)
                    currentcflux=str_sep(strtrim(strcompress(tirificfirst[tmppos]),2),'=')
                    IF double(string(double(currentcflux[1])*nopoints[j]/(1e6*factor))) NE 0. then $
                       tirificfirst[tmppos]='CFLUX= '+string(double(currentcflux[1])*nopoints[j]/(1e6*factor)) else $
                          tirificfirst[tmppos]='CFLUX= 1e-5'
                 ENDIF ELSE BEGIN
                    tmppos=where('CFLUX_'+strtrim(strcompress(fix(j+1)),2) EQ tirificfirstvars)
                    currentcflux1=str_sep(strtrim(strcompress(tirificfirst[tmppos]),2),'=')
                    IF double(string(double(currentcflux1[1])*nopoints[j]/(1e6*factor))) NE 0. then $
                       tirificfirst[tmppos]=currentcflux1[0]+'= '+string(double(currentcflux1[1])*nopoints[j]/(1e6*factor)) else $
                          tirificfirst[tmppos]=currentcflux1[0]+'= 1e-5'
                 ENDELSE
                 if j eq n_elements(nopoints)-1 then begin
                    IF n_elements(currentcflux) EQ 0 then begin
                       tmppos=where('CFLUX' EQ tirificfirstvars)
                       currentcflux=str_sep(strtrim(strcompress(tirificfirst[tmppos]),2),'=')
                    ENDIF
                    IF size(log,/TYPE) EQ 7 then begin
                       openu,66,log,/APPEND
                       printf,66,linenumber()+"CHECK_CFLUX: CFLUX is scaled down to"+STRJOIN([string(double(currentcflux[1])*nopoints[0]/(1e6*factor)),string(double(currentcflux1[1])*nopoints[1]/(1e6*factor))],' and ')
                       close,66
                    endif
                 endif
                 cfluxadjusted=1.
              end
              nopoints[j] GT 2.2E6*factor: begin
                 IF j EQ 0 then begin
                    tmppos=where('CFLUX' EQ tirificfirstvars)
                    currentcflux=str_sep(strtrim(strcompress(tirificfirst[tmppos]),2),'=')
                    IF double(string(double(currentcflux[1])*nopoints[j]/(1e6*factor))) LT 0.005 then $
                       tirificfirst[tmppos]='CFLUX= '+string(double(currentcflux[1])*nopoints[j]/(1E6*factor)) else $
                          tirificfirst[tmppos]='CFLUX= 0.005'
                 ENDIF ELSE BEGIN
                    tmppos=where('CFLUX_'+strtrim(strcompress(fix(j+1)),2) EQ tirificfirstvars)
                    currentcflux1=str_sep(strtrim(strcompress(tirificfirst[tmppos]),2),'=')
                    IF double(string(double(currentcflux1[1])*nopoints[j]/(1e6*factor))) LT 0.005 then $
                       tirificfirst[tmppos]=currentcflux1[0]+'= '+string(double(currentcflux1[1])*nopoints[j]/(1E6*factor)) else $
                          tirificfirst[tmppos]=currentcflux1[0]+'= 0.005'
                 ENDELSE
                 if j eq n_elements(nopoints)-1 then begin
                    IF n_elements(currentcflux) EQ 0 then begin
                       tmppos=where('CFLUX' EQ tirificfirstvars)
                       currentcflux=str_sep(strtrim(strcompress(tirificfirst[tmppos]),2),'=')
                    ENDIF
                    IF size(log,/TYPE) EQ 7 then begin
                       openu,66,log,/APPEND
                       printf,66,linenumber()+"CHECK_CFLUX: CFLUX is scaled up to"+STRJOIN([string(double(currentcflux[1])*nopoints[0]/(1e6*factor)),string(double(currentcflux1[1])*nopoints[1]/(1e6*factor))],' and ')
                       close,66
                    endif
                 endif
              end  
              else:begin
                 if j eq n_elements(nopoints)-1 then begin
                    IF size(log,/TYPE) EQ 7 then begin
                       openu,66,log,/APPEND
                       printf,66,linenumber()+"CHECK_CFLUX: CFLUX is ok."
                       close,66
                    endif
                 endif
              end
           ENDCASE
        ENDELSE
     ENDELSE
  ENDFOR
END  
