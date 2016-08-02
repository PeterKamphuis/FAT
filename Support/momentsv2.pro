Pro momentsv2,Cube,Momentmap,header,map

;+
; NAME:
;       MOMENTSV2
;
; PURPOSE:
;       Program to calculate the moment zero map of a 3D fits cube  and update the header.
;
; CATEGORY:
;       Support
; 
; CALLING SEQUENCE:
;       momentsv2,Cube,Momentmap,header,map
;
;
; INPUTS:
;       Cube = A 3D-array containing the values of the cube 
;       header = The fits header of the data cube
;       map = 1 or 0 indicating the requested moment
;
; OPTIONAL INPUTS:
;       - 
;
; KEYWORD PARAMETERS:
;       -
;
; OUTPUTS:
;       momentmap = the calaculated moment map
;       header = the header will be updated to reflect the reduction
;       in dimension and proper units.
;
; OPTIONAL OUTPUTS:
;       -
; 
; PROCEDURES CALLED:
;       STR_SEP(), STRTRIM(), STRCOMPRESS(), TOTAL(), SXADDPAR(),
;       SXPAR(), SXDELPAR()
;
; EXAMPLE:
;      
;
; MODIFICATION HISTORY:
;       01-06-2016 P.Kamphuis; Added a condition to check that datamax
;                              and datamin are finite.   
;       07-01-2016 P.Kamphuis; Replaced SUM commands with the proper
;       TOTAL commands reducing the need for outside routines.  
;       Modified to deal with existing CUNIT3 without error message 12-08-2015 P. Kamphuis v2
;       Modified to deal with missing CUNIT3  29-07-2015 P. Kamphuis v2
;       Modified to use SUM which is much faster 25-05-2015 P. Kamphuis v2
;       Written 15-07-2010 P.Kamphuis v1.0
;
; NOTE:
;     
;-




IF map EQ 0 then begin
   blank=WHERE(FINITE(Cube) NE 1.)
   IF blank[0] NE -1 then Cube[blank]=0
   Momentmap=fltarr(n_elements(Cube[*,0,0]),n_elements(Cube[0,*,0]))
   Momentmap[*,*]=TOTAL(Cube,3)*ABS(sxpar(header,'CDELT3'))
   IF isnumeric(sxpar(header,'CUNIT3')) then begin
      IF sxpar(header,'CDELT3') GT 500. then sxaddpar,header,'CUNIT3','M/S' else sxaddpar,header,'CUNIT3','KM/S'
   ENDIF
   IF STRUPCASE(strtrim(sxpar(header,'CUNIT3'),2)) EQ 'M/S' then begin
      sxaddpar,header,'CUNIT3','KM/S'
      momentmap=momentmap/1000.
      print,linenumber()+'MOMENTSV2: We have converted the units to Jy/Beam x Km/s'
   ENDIF
   sxaddpar,header,'BUNIT',strtrim(strcompress(sxpar(header,'BUNIT')),2)+'.'+strtrim(strcompress(sxpar(header,'CUNIT3')),2)
   maxmap=MAX(momentmap,MIN=minmap)
   if FINITE(maxmap) then sxaddpar,header,'DATAMAX',maxmap else sxaddpar,header,'DATAMAX',MAX(Cube)
   if FINITE(minmap) then sxaddpar,header,'DATAMIN',minmap else sxaddpar,header,'DATAMIN',0.

endif
IF map EQ 1 then begin
   buildaxii,header,xaxis,yaxis,zaxis=zaxis
   blank=WHERE(FINITE(Cube) NE 1.)
   IF blank[0] NE -1 then Cube[blank]=0  
   IF isnumeric(sxpar(header,'CUNIT3')) then begin
      IF sxpar(header,'CDELT3') GT 500. then sxaddpar,header,'CUNIT3','M/S' else sxaddpar,header,'CUNIT3','KM/S'
   ENDIF
   IF STRUPCASE(strtrim(sxpar(header,'CUNIT3'),2)) EQ 'M/S' then begin
      print,linenumber()+'MOMENTSV2: We are converting to KM/S'
      zaxis=zaxis/1000.
   ENDIF
   Momentmap=fltarr(n_elements(Cube[*,0,0]),n_elements(Cube[0,*,0]))
   c=rebin(reform(zaxis,1,1,n_elements(zaxis)),n_elements(Cube[*,0,0]),n_elements(Cube[0,*,0]),n_elements(Cube[0,0,*]))
   Momentmap=TOTAL(c*Cube,3)/TOTAL(Cube,3)
   maxmap=MAX(momentmap,MIN=minmap)
   if FINITE(maxmap) then sxaddpar,header,'DATAMAX',maxmap else sxaddpar,header,'DATAMAX',zaxis[n_elements(zaxis)-1]
   if FINITE(minmap) then sxaddpar,header,'DATAMIN',minmap else sxaddpar,header,'DATAMIN',zaxis[0]
   sxaddpar,header,'BUNIT','KM/S'
   c=0
ENDIF

sxdelpar,header,'CUNIT3'
sxdelpar,header,'CTYPE3'
sxdelpar,header,'CRVAL3'
sxdelpar,header,'CDELT3'
sxdelpar,header,'NAXIS3'
sxdelpar,header,'CRPIX3'
sxaddpar,header,'NAXIS',2

end


