Pro read_template, name, array, variables,SOFIA=sofia

;+
; NAME:
;       READ_TEMPLATE
;
; PURPOSE:
;       Routine to read in a specific template file with = defined commands.
;
; CATEGORY:
;       Support
;
; CALLING SEQUENCE:
;       READ_TEMPLATE, name, array, variables, /SOFIA
;
;
; INPUTS:
;       Name = name of the file to be read.
;
; OPTIONAL INPUTS:
;       -
;
; KEYWORD PARAMETERS:
;       /SOFIA - trigger to indicate it is a sofia input template not tirific
;
; OUTPUTS:
;       array = An array with each line of the input file
;       variable = an array with all the variables (i.e. the words
;       appaearing before =) of the file
; OPTIONAL OUTPUTS:
;       -
;
; PROCEDURES CALLED:
;       FILE_LINES(), STR_SEP(), STRTRIM(), STRCOMPRESS()
;
; EXAMPLE:
;
;
; MODIFICATION HISTORY:
;       Written 01-01-2015 P.Kamphuis v1.0
;
; NOTE:
;
;-
  COMPILE_OPT IDL2
  IF not keyword_set(sofia) then begin
     openr,1,name
     filelength=FILE_LINES(name)
     h=' '
     array=strarr(filelength)
     variables=strarr(filelength)
     for j=0,filelength-1 do begin
        readf,1,h
        tmp=str_sep(strtrim(strcompress(h),2),'=')
        IF n_elements(tmp) GT 1 then begin
           variables[j]=tmp[0]
        endif
        array[j]=h
     endfor
     close,1
  ENDIF ELSE BEGIN
     openr,1,name
     filelength=FILE_LINES(name)
     h=' '
     array=strarr(filelength)
     variables=intarr(7)
     for j=0,filelength-1 do begin
           readf,1,h
           tmp=strtrim(strcompress(str_sep(h,'=')),2)
           IF n_elements(tmp) GT 1 then begin
              case tmp[0] of
                 'import.inFile':variables[0]=j
                 'steps.doReliability':variables[1]=j
                 'parameters.dilatePixMax':variables[2]=j
                 'SCfind.threshold':variables[3]=j
                 'SCfind.kernels':variables[4]=j
                 'steps.doFlag':variables[5]=j
                 'flag.regions':variables[6]=j
                 else:
              endcase
           ENDIF
           array[j]=h
        endfor
     close,1
  ENDELSE
end
