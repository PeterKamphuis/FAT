function isnumeric,input

;+
; NAME:
;       ISNUMERIC
;
; PURPOSE:
;       Routine to determine wether a scalar is numeric
;
; CATEGORY:
;       Support
; 
; CALLING SEQUENCE:
;       Result = isnumeric(input)
;
; INPUTS:
;       input = the value to be checked
;
; OPTIONAL INPUTS:
;       - 
;
; KEYWORD PARAMETERS:
;       -
;
; OUTPUTS:
;       result= binary to check wether the value is numeric 0 is no 1
;       is yes.
;
; OPTIONAL OUTPUTS:
;       -
; 
; PROCEDURES CALLED:
;       -
;
; EXAMPLE:
;      
;
; MODIFICATION HISTORY:
;       15-11-2016 P.Kamphuis; Added full string length capability to
;                              check against numeric excepting '+-.ED'
;                              as otherwise a single digit is enough
;                              to make the string numeric.  
;       Written by ??????
;
; NOTE:
;     
;-
  COMPILE_OPT IDL2
                                ;if double returns an error then the string is not numeric
  on_ioerror, false
                                ;Get the length of the string
  tmp=strlen(input)
                                ;Keep a set of counters for numeric signs
  pluscount=0
  minuscount=0
  Ecount=0
  Dcount=0
  dotcount=0
                                ;Loop through each element of the string to check it is numeric
  for i=0,tmp-1 do begin
     tmpcheck=STRMID(input,i,1)
     case 1 of
        STRUPCASE(tmpcheck) EQ '+' AND pluscount LT 1:pluscount++
        STRUPCASE(tmpcheck) EQ '-' AND minuscount LT 2:minuscount++
        STRUPCASE(tmpcheck) EQ '.' AND dotcount LT 2:dotcount++
        STRUPCASE(tmpcheck) EQ 'E' AND Ecount LT 1:Ecount++
        STRUPCASE(tmpcheck) EQ 'D' AND Dcount LT 1:Dcount++
        else: test = double(tmpcheck)
     endcase
  endfor
                                ;Make sure the total is numeric
  test=double(input)
                                ;Return 1 on succes 0 on error
  return, 1
  false: return, 0
end
