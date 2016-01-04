Pro rename,stringin,stringreplace

  
;+
; NAME:
;       RENAME
;
; PURPOSE:
;      Rename Tirific Output that is produced while the pipeline runs
;
; CATEGORY:
;       Support
; 
; CALLING SEQUENCE:
;       RENAME,stringin,stringreplace,files
;
;
; INPUTS:
;      stringin = the string to be replaced
;      stringreplace = The string that will be put in
;
; OPTIONAL INPUTS:
;
; KEYWORD PARAMETERS:
;       -
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;       -
; 
; PROCEDURES CALLED:
;       SPAWN
;
; EXAMPLE:
;      
;
; MODIFICATION HISTORY:
;       Written 04-01-2016 P.Kamphuis v1.0
;
; NOTE: This only works for the set of files produced by tirific this
; is not a generic rename
;     
;-

extensions=['def','log','ps','fits']

for i=0,n_elements(extensions)-1 do begin
   IF FILE_TEST(stringreplace+extensions[i]) then spawn,'rm -f '+stringreplace+extensions[i]
   IF FILE_TEST(stringin+extensions[i]) then spawn,'mv '+stringin+extensions[i]+' '+stringreplace+extensions[i]
endfor

end

   
