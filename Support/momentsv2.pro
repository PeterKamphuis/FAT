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
;       STR_SEP(), STRTRIM(), STRCOMPRESS(), SUM(), SXADDPAR(),
;       SXPAR(), SXDELPAR()
;
; EXAMPLE:
;      
;
; MODIFICATION HISTORY:
;       Modified to use SUM which is much faster 25-05-2015 P. Kamphuis v2
;       Written 15-07-2010 P.Kamphuis v1.0
;
; NOTE:
;     
;-




IF map EQ 0 then begin
   blank=WHERE(FINITE(Cube) NE 1.)
   IF blank[0] NE -1 then Cube[blank]=0
   Momentmap=dblarr(n_elements(Cube[*,0,0]),n_elements(Cube[0,*,0]))
   Momentmap[*,*]=SUM(Cube,2)*ABS(sxpar(header,'CDELT3'))
   IF STRUPCASE(strtrim(sxpar(header,'CUNIT3'),2)) EQ 'M/S' then begin
      sxaddpar,header,'CUNIT3','KM/S'
      momentmap=momentmap/1000.
      print,'We have converted the units to Jy/Beam x Km/s'
   ENDIF
   sxaddpar,header,'BUNIT',strtrim(strcompress(sxpar(header,'BUNIT')),2)+'.'+strtrim(strcompress(sxpar(header,'CUNIT3')),2)
   sxaddpar,header,'DATAMAX',MAX(momentmap,MIN=minmap)
   sxaddpar,header,'DATAMIN',minmap
endif
IF map EQ 1 then begin
   buildaxii,header,xaxis,yaxis,zaxis=zaxis
   blank=WHERE(FINITE(Cube) NE 1.)
   IF blank[0] NE -1 then Cube[blank]=0
   Momentmap=dblarr(n_elements(Cube(*,0,0)),n_elements(Cube(0,*,0)))

   c=rebin(reform(zaxis,1,1,n_elements(zaxis)),n_elements(Cube[*,0,0]),n_elements(Cube[0,*,0]),n_elements(Cube[0,0,*]))
   Momentmap=SUM(c*Cube,2)/sum(Cube,2)
   IF strtrim(sxpar(header,'CUNIT3'),2) EQ 'M/S' or strtrim(sxpar(header,'CUNIT3'),2) eq 'm/s' then begin
      print,'We are converting to KM/S'
      Momentmap=Momentmap/1000.
   ENDIF
   sxaddpar,header,'DATAMAX',MAX(momentmap,MIN=minmap)
   sxaddpar,header,'DATAMIN',minmap
   sxaddpar,header,'BUNIT','KM/S'
ENDIF

sxdelpar,header,'CUNIT3'
sxdelpar,header,'CTYPE3'
sxdelpar,header,'CRVAL3'
sxdelpar,header,'CDELT3'
sxdelpar,header,'NAXIS3'
sxdelpar,header,'CRPIX3'
sxaddpar,header,'NAXIS',2

end


