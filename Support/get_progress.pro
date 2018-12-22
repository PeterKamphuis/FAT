Pro get_progress,name,AC,nopoints,loops,models

;+
; NAME:
;       GET_PROGRESS
;
; PURPOSE:
;       This routine analyses the progress log from tirific and
;       extracts the acceppted parameter, number of points, loops and
;       total amount of models 
;
; CATEGORY:
;       Support
; 
; CALLING SEQUENCE:
;       GET_PROGRESS,name,AC,nopoints,loops,models
;
;
; INPUTS:
;       name = name of the progress file 
;       
; OPTIONAL INPUTS:
;       - 
;
; KEYWORD PARAMETERS:
;       -
;
; OUTPUTS:
;       AC = value of the acceptance parameter
;       nopoints = number of point sources in the model. 2D array for
;       both sides
;       loops = amount of completed big loops
;       models = number of models produced
;
; OPTIONAL OUTPUTS:
;       -
; 
; PROCEDURES CALLED:
;       STRUPCASE(), STRTRIM() STRCOMPRESS(), STR_SEP()
;
; EXAMPLE:
;      
;
; MODIFICATION HISTORY:
;       Written by P.Kamphuis 01-01-2015 
;
; NOTE:
;     
;-
  COMPILE_OPT IDL2
  AC=dblarr(1)
  nopoints=dblarr(1)
  openr,1,name
  h=' '
                                ;Just need the first line
  readf,1,h
  tmp=str_sep(strtrim(strcompress(h),2),' ')
  close,1
   ;AC is the ninth item (In my modified Tirific)
  for i=0,n_elements(tmp)-1 do begin
     tmp1=str_sep(strtrim(strcompress(tmp[i]),2),':')
     IF STRUPCASE(tmp1[0]) EQ 'N' then begin
        tmp2=str_sep(strtrim(strcompress(tmp1[1]),2),'/')
        nopoints=dblarr(n_elements(tmp2))
        nopoints=double(tmp2)
     ENDIF
     IF STRUPCASE(tmp1[0]) EQ 'AC' then AC=double(tmp1[1])
     IF STRUPCASE(tmp1[0]) EQ 'BL' then loops=double(tmp1[1])
     IF STRUPCASE(tmp1[0]) EQ 'TM' then models=double(tmp1[1])
  endfor
                                ;fitmode two does not give all this output so for now
  AC= 1
  loops=10.
  models=1000.
  
end
