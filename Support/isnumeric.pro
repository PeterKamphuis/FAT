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
;       Written by ??????
;
; NOTE:
;     
;-
  COMPILE_OPT IDL2 
  on_ioerror, false
  test = double(input)
  return, 1
  false: return, 0
end
