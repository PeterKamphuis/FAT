Pro changeradii,tirifictemplate,numberofrings

;+
; NAME:
;       CHANGERADII
;
; PURPOSE:
;       Program to reduce the number of rings in a tirific template by 1 ring
;
; CATEGORY:
;       Support
; 
; CALLING SEQUENCE:
;       CHANGERADII,tirifictemplate,numberofrings
;
;
; INPUTS:
;       tirifictemplate =  a string array with all the input of a
;       tirific def file. 
;       numberofrings = the new numberof ring requested
;
; OPTIONAL INPUTS:
;       -
;
; KEYWORD PARAMETERS:
;       -
;
; OUTPUTS:
;       tirifictemplate =  the modified string array with all the
;       input of a tirific def file
;
; OPTIONAL OUTPUTS:
;       -
;
; PROCEDURES CALLED:
;       STRTRIM(),STRCOMPRESS(),STRING(),STR_SEP(),WHERE(),DOUBLE(),FLOOR()
;
; EXAMPLE:
;
;
; MODIFICATION HISTORY:
;     Written by P.Kamphuis 01-01-2015
; 
; NOTE:
;
;-
  COMPILE_OPT IDL2 

possiblevars=['VROT','Z0','SBR','INCL','PA','XPOS','YPOS','VSYS','AZ1W','AZ1P','SDIS']     

for i=0,n_elements(tirifictemplate)-1 do begin
   tmp=str_sep(strtrim(strcompress(tirifictemplate[i]),2),'=')
   IF tmp[0] EQ 'NUR' then begin
      oldrings=tmp[1]
      tirifictemplate[i]=tmp[0]+'='+strtrim(strcompress(string(numberofrings)),1)
      break
   ENDIF
endfor

if oldrings GT numberofrings then oldrings=numberofrings

for i=0,n_elements(tirifictemplate)-1 do begin
                                ;first we split at the is sign
   tmp=str_sep(strtrim(strcompress(tirifictemplate[i]),2),'=')    
   IF n_elements(tmp) GT 1 then begin
                                ;then we check wether we have a disk variable
      checkdisk=str_sep(strtrim(strcompress(tmp[0]),2),'_')
      IF n_elements(checkdisk) GT 1 then begin
         currentvar=checkdisk[0]
      endif else begin
         currentvar=tmp[0]
      ENDELSE
      
      found=WHERE(currentvar EQ possiblevars)
      IF found[0] NE -1 then begin
         tmp2=str_sep(strtrim(strcompress(tmp[1]),2),' ')
         outputstring=strarr(1)
         outputstring=tmp[0]+'='
         IF n_elements(tmp2) GT 1 then begin
            IF  n_elements(tmp2) LT numberofrings then amount=n_elements(tmp2)-1 else amount=numberofrings-1
            for j=0,amount do begin
               outputstring=outputstring+tmp2[j]+' '
            endfor
            tirifictemplate[i]=outputstring
         ENDIF ELSE begin
            tirifictemplate[i]=tmp[0]+'='+tmp[1]
         ENDELSE
         
      ENDIF     
      IF tmp[0] EQ 'RADI' then begin
         tmp2=str_sep(strtrim(strcompress(tmp[1]),2),' ')
         outputstring=strarr(1)
         outputstring=tmp[0]+'='
         IF n_elements(tmp2) LT numberofrings then begin
            for j=0,n_elements(tmp2)-1 do begin
               outputstring=outputstring+tmp2[j]+' '
            endfor
            for j=n_elements(tmp2),numberofrings-1 do begin
               outputstring=outputstring+strtrim(strcompress(string(double(tmp2[n_elements(tmp2)-1])+(j-double(n_elements(tmp2))+1)*(double(tmp2[n_elements(tmp2)-1])-double(tmp2[n_elements(tmp2)-2])))),1)+' '
            endfor
         ENDIF ELSE BEGIN
            for j=0,numberofrings-1 do begin
               outputstring=outputstring+tmp2[j]+' '
            endfor
         ENDELSE
         tirifictemplate[i]=outputstring
      ENDIF
      IF tmp[0] EQ 'VARY' then begin
         tirifictemplate[i]=tmp[0]+'= '
         tmpsep1=str_sep(strtrim(strcompress(tmp[1]),2),',')
         for j=0,n_elements(tmpsep1)-1 do begin
            tmpsep2=str_sep(strtrim(strcompress(tmpsep1[j]),2),' ')
            for x=0,n_elements(tmpsep2)-1 do begin 
               tmpsep3=str_sep(strtrim(strcompress(tmpsep2[x]),2),':')
               IF tmpsep3[0] EQ -1 then begin
                  IF isnumeric(tmpsep2[x]) then begin
                     IF Double(tmpsep2[x]) GE double(oldrings) then tmpsep2[x]=strtrim(strcompress(string(numberofrings,format='(F7.4)')),1)
                  ENDIF
               endif ELSE BEGIN
                  for y=0,n_elements(tmpsep3)-1 do begin
                     IF isnumeric(tmpsep3[y]) then begin
                        IF Double(tmpsep3[y]) GE double(oldrings) then tmpsep3[y]=strtrim(strcompress(string(numberofrings,format='(F7.4)')),1)
                     ENDIF
                     IF y EQ 0 then tmpsep2[x]=tmpsep3[y] else tmpsep2[x]=tmpsep2[x]+':'+tmpsep3[y]
                  endfor
               ENDELSE
              
               IF x EQ 0 then tmpsep1[j]=tmpsep2[x] else tmpsep1[j]=tmpsep1[j]+' '+tmpsep2[x]
            endfor
            IF j EQ 0 then tirifictemplate[i]=tirifictemplate[i]+tmpsep1[j] else tirifictemplate[i]=tirifictemplate[i]+','+tmpsep1[j]
            
         endfor
         
      ENDIF 
;this part should be adjusted to accomodate any number of
;regularization parameters
      IF tmp[0] EQ 'REGNUME' then begin
         tirifictemplate[i]='REGNUME= '+string(floor((numberofrings)/2.))+','+string(floor((numberofrings)/2.))
      endif
      IF tmp[0] EQ 'REGDENO' then begin
                                ;first we need to get the low number
         
         tmpnow=str_sep(strtrim(strcompress(tmp[1]),2),',')   
         tmplow=str_sep(strtrim(strcompress(tmpnow[0]),2),':')
         if double(tmplow[0]) EQ floor(numberofrings/2.) then tmplow[0]=string(tmplow-1)
         tirifictemplate[i]='REGDENO= '+tmplow[0]+':'+strtrim(strcompress(string(floor((numberofrings)/2.))),2)+','+tmplow[0]+':'+strtrim(strcompress(string(floor((numberofrings)/2.))),2)
      endif
   ENDIF
ENDFOR


end
