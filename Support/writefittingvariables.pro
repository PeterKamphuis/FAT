Pro writefittingvariables,inputarray,v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,v13,v14,v15, $
                          v16,v17,v18,v19,v20,v21,v22,v23,v24,v25,v26,v27,v28,v29,v30,$
                          v31,v32,v33,v34,v35,v36,v37,v38,v39,v40

;+
; NAME:
;       WRITEFITTINGVARIABLES
;
; PURPOSE:
;       Routine to write the fitting variables to the tirific template
;       array
;
; CATEGORY:
;       Support
; 
; CALLING SEQUENCE:
;       WRITEFITTINGVARIABLES,inputarray,v1[,v2,v3,...v40]
;
;
; INPUTS:
;       inputarray = tirific template
;       v1-v40 = Arrays with the various fitting parameters that
;       should be written to the template
;
; OPTIONAL INPUTS:
;       - 
;
; KEYWORD PARAMETERS:
;       -
;
; OUTPUTS:
;       inputarray = updated tirific template
;
; OPTIONAL OUTPUTS:
;       -
; 
; PROCEDURES CALLED:
;       STR_SEP(), STRTRIM(), STRCOMPRESS(), STRING()
;
; EXAMPLE:
;      
;
; MODIFICATION HISTORY:
;       21-10-2018 P.Kamphuis; Modified to work with Fitmode 2 in tirific  
;       Written 01-01-2015 P.Kamphuis v1.0
;
; NOTE:
;      Fitting VARIABLE input
;  order= [Vary, Parmax,parmin,delstart,delend,satdelt,mindelt,moderate,ITESTART,ITEEND,VARINDX]
;      or
;  order= [Vary, Parmax,parmin,delstart,delend,mindelt,moderate,VARINDX]
;-
  COMPILE_OPT IDL2 
  On_error,2                    ;Return to caller
;We need at least three input parameters
  IF N_params() LT 2 then begin
     print,'WRITEFITTINGVARIABLES: You need to change at least one variable'
     return
  endif
  
  nvar = N_params() - 1         ;Number of variables to deal with
;Let's check wether the amount of limits and variables match

;let's build  some arrays from the input
  strings=strarr(8)
  strings[0]='VARY= '    
  strings[1]='PARMAX= '
  strings[2]='PARMIN= '
  strings[3]='MODERATE= '
  strings[4]='DELSTART= '
  strings[5]='DELEND= '
  ;strings[6]='ITESTART= '
  ;strings[7]='ITEEND= '
  ;strings[8]='SATDELT= '
  strings[6]='MINDELTA= '             
  strings[7]='VARINDX= '
  vv = 'v' + strtrim( indgen(nvar)+1, 2)


  for i=0,nvar-1 do begin
     limits=SCOPE_VARFETCH(vv[i],LEVEL=0)
                                ;Allow for old writing still into
                                ;proper new mode such that not all
                                ;fitting arrays need to be modfied immediately
     IF n_elements(limits) EQ 10 then begin
        strings[0]=strings[0]+strtrim(strcompress(string(limits[0])))+','
        strings[1]=strings[1]+strtrim(strcompress(string(limits[1])))+' '
        strings[2]=strings[2]+strtrim(strcompress(string(limits[2])))+' '
        strings[3]=strings[3]+strtrim(strcompress(string(limits[7])))+' '
        strings[4]=strings[4]+strtrim(strcompress(string(limits[3])))+' '
        strings[5]=strings[5]+strtrim(strcompress(string(limits[4])))+' '
        ;strings[6]=strings[6]+strtrim(strcompress(string(limits[8])))+' '
        ;strings[7]=strings[7]+strtrim(strcompress(string(limits[9])))+' '
        ;strings[8]=strings[8]+strtrim(strcompress(string(limits[5])))+' '
        strings[6]=strings[6]+strtrim(strcompress(string(limits[6])))+' '
     ENDIF
     IF n_elements(limits) EQ 11 then begin
        strings[0]=strings[0]+strtrim(strcompress(string(limits[0])))+','
        strings[1]=strings[1]+strtrim(strcompress(string(limits[1])))+' '
        strings[2]=strings[2]+strtrim(strcompress(string(limits[2])))+' '
        strings[3]=strings[3]+strtrim(strcompress(string(limits[7])))+' '
        strings[4]=strings[4]+strtrim(strcompress(string(limits[3])))+' '
        strings[5]=strings[5]+strtrim(strcompress(string(limits[4])))+' '
        ;strings[6]=strings[6]+strtrim(strcompress(string(limits[8])))+' '
        ;strings[7]=strings[7]+strtrim(strcompress(string(limits[9])))+' '
        ;strings[8]=strings[8]+strtrim(strcompress(string(limits[5])))+' '
        strings[6]=strings[6]+strtrim(strcompress(string(limits[6])))+' '
        strings[7]=strings[7]+strtrim(strcompress(string(limits[10])))+' '
     ENDIF
     ;New writing start here  strings=strarr(8)
     IF n_elements(limits) EQ 7 then begin order= [Vary, Parmax,parmin,delstart,delend,mindelt,moderate,VARINDX]
        strings[0]=strings[0]+strtrim(strcompress(string(limits[0])))+','
        strings[1]=strings[1]+strtrim(strcompress(string(limits[1])))+' '
        strings[2]=strings[2]+strtrim(strcompress(string(limits[2])))+' '
        strings[3]=strings[3]+strtrim(strcompress(string(limits[6])))+' '
        strings[4]=strings[4]+strtrim(strcompress(string(limits[3])))+' '
        strings[5]=strings[5]+strtrim(strcompress(string(limits[4])))+' '
        strings[6]=strings[6]+strtrim(strcompress(string(limits[5])))+' '
     ENDIF
     IF n_elements(limits) EQ 8 then begin
        strings[0]=strings[0]+strtrim(strcompress(string(limits[0])))+','
        strings[1]=strings[1]+strtrim(strcompress(string(limits[1])))+' '
        strings[2]=strings[2]+strtrim(strcompress(string(limits[2])))+' '
        strings[3]=strings[3]+strtrim(strcompress(string(limits[6])))+' '
        strings[4]=strings[4]+strtrim(strcompress(string(limits[3])))+' '
        strings[5]=strings[5]+strtrim(strcompress(string(limits[4])))+' '
        strings[6]=strings[6]+strtrim(strcompress(string(limits[5])))+' '
        strings[7]=strings[7]+strtrim(strcompress(string(limits[7])))+' '  
     ENDIF
  endfor
  strings[0]=STRMID(strings[0], 0 , STRLEN(strings[0])-1)



  for i=0,n_elements(inputarray)-1 do begin
     tmp=str_sep(strtrim(strcompress(inputarray[i]),2),'=')
     case tmp[0] of
        'VARY':begin  
           inputarray[i]=strings[0]
        end
        'PARMAX':begin  
           inputarray[i]=strings[1]
        end
        'PARMIN':begin  
           inputarray[i]=strings[2]
        end
        'MODERATE':begin  
           inputarray[i]=strings[3]
        end
        'DELSTART':begin  
           inputarray[i]=strings[4]
        end
        'DELEND':begin  
           inputarray[i]=strings[5]
        end
;        'ITESTART':begin  
;           inputarray[i]=strings[6]
;        end
;        'ITEEND':begin  
;           inputarray[i]=strings[7]
;        end
;        'SATDELT':begin  
;           inputarray[i]=strings[8]
;        end
        'MINDELTA':begin  
           inputarray[i]=strings[6]
        end
        'VARINDX':begin  
           inputarray[i]=strings[7]
        end

        else:begin
        end

     endcase
  endfor

end

