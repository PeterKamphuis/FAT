Pro momentsv2,Cube,Momentmap,header,map,BLANK_VALUE=blanked,gdlidl=gdlidl


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
;       gdlidl = indicator of running gdl or idl. Is 1 when running GDL
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
;       20-08-2017 P.Kamphuis; Rebin is broken on the mac 0.9.7 hence
;                              we have added a condition to rebin in
;                              loop when running GDL. This is only
;                              used for velocity fields so no need to
;                              add to intergrated moment map
;       12-05-2017 P.Kamphuis; In the  maps the blanks should
;                              propagate.
;       16-11-2016 P.Kamphuis; Now dealing with the (non-)presence of
;                              CUNIT3 properly
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
     ;In the case of moment 0 we want blanks to propagate through.
   ;blank=WHERE(FINITE(Cube) NE 1.)
   ;IF blank[0] NE -1 then Cube[blank]=0
   Momentmap=fltarr(n_elements(Cube[*,0,0]),n_elements(Cube[0,*,0]))
   Momentmap[*,*]=TOTAL(Cube,3)*ABS(sxpar(header,'CDELT3'))
   IF ~(sxpar(header,'CUNIT3')) then begin
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
IF map GE 1 then begin
   buildaxii,header,xaxis,yaxis,zaxis=zaxis
   ;blank=WHERE(FINITE(Cube) NE 1.)
   ;IF blank[0] NE -1 then Cube[blank]=0
   IF ~(sxpar(header,'CUNIT3')) then begin
      IF sxpar(header,'CDELT3') GT 500. then sxaddpar,header,'CUNIT3','M/S' else sxaddpar,header,'CUNIT3','KM/S'
   ENDIF
   IF STRUPCASE(strtrim(sxpar(header,'CUNIT3'),2)) EQ 'M/S' then begin
      print,linenumber()+'MOMENTSV2: We are converting to KM/S'
      zaxis=zaxis/1000.
   ENDIF
   Momentmap=fltarr(n_elements(Cube[*,0,0]),n_elements(Cube[0,*,0]))

   if gdlidl then begin
                                ;this is necessary for gdl as rebin is
                                ;broken in gdl on the mac. This seems
                                ;an issue of 0.9.7. We should check on ubuntu and
                                ;issu a ticket or some such.
      c=dblarr(n_elements(Cube[*,0,0]),n_elements(Cube[0,*,0]),n_elements(Cube[0,0,*]))
      for i=0,n_elements(Cube[0,0,*])-1 do begin
         c[*,*,i]=zaxis[i]
      endfor
   endif else begin
      c=rebin(reform(zaxis,1,1,n_elements(zaxis)),n_elements(Cube[*,0,0]),n_elements(Cube[0,*,0]),n_elements(Cube[0,0,*]))
   endelse
   Momentmap=TOTAL(c*Cube,3)/TOTAL(Cube,3)
   IF map GE 2 then begin
      if gdlidl then begin
         b=dblarr(n_elements(Cube[*,0,0]),n_elements(Cube[0,*,0]),n_elements(Cube[0,0,*]))
         for i=0,n_elements(Cube[0,0,*])-1 do begin
            b[*,*,i]=Momentmap[*,*]
         endfor
      endif else begin
         b=rebin(reform(Momentmap,n_elements(Momentmap[0,*]),n_elements(Momentmap[*,0]),1),n_elements(Cube[*,0,0]),n_elements(Cube[0,*,0]),n_elements(Cube[0,0,*]))
      endelse
      d=c-b
      b=1.
      c=1.
      Momentmap=SQRT(TOTAL(Cube*(d)^2,3)/TOTAL(Cube,3))
      d=0.
   ENDIF
   maxmap=MAX(momentmap[where(FINITE(momentmap) eq 1)],MIN=minmap)
   if FINITE(maxmap) then sxaddpar,header,'DATAMAX',maxmap else sxaddpar,header,'DATAMAX',zaxis[n_elements(zaxis)-1]
   if FINITE(minmap) then sxaddpar,header,'DATAMIN',minmap else sxaddpar,header,'DATAMIN',zaxis[0]
   sxaddpar,header,'BUNIT','KM/S'
   c=0
ENDIF
if n_elements(blanked) NE 0  then begin
   tmp=WHERE(momentmap EQ 0.)
   IF tmp[0] NE -1 then momentmap[tmp]=blanked
ENDIF
sxdelpar,header,'CUNIT3'
sxdelpar,header,'CTYPE3'
sxdelpar,header,'CRVAL3'
sxdelpar,header,'CDELT3'
sxdelpar,header,'NAXIS3'
sxdelpar,header,'CRPIX3'
sxaddpar,header,'NAXIS',2

end
