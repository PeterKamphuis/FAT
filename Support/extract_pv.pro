Pro extract_pv,Cube,header,pa,xv,CENTER=center,XVHEADER=new_header

;+
; NAME:
;       EXTRACT_PV
;
; PURPOSE:
;       Program to extract a PV diagram along a Position Angle from a
;       Line emission cube. The strip width along which to extract the
;       XV-diagram is the beam's major axis FWHM  if stated in
;       the header else a single pixel is used. 
;
; CATEGORY:
;       Support
; 
; CALLING SEQUENCE:
;       EXTRACT_PV,Cube,header,pa,xv,CENTER=center,XVHEADER=new_header
;
;
; INPUTS:
;       Cube = The array containing the Data Cube
;       header = The header of the data cube
;       pa = The position angle along which to extract this runs from
;       north east wards
;
; OPTIONAL INPUTS:
;       - 
;
; KEYWORD PARAMETERS:
;       CENTER = if set the rotation will happen around this center,
;       if blank it will rotate around the center set in the header of
;       the cube
;       XVHEADER = A new header to write the xv array to a fits files
;
; OUTPUTS:
;       xv = a 2-dimensional array with the XV - Diagram
;
; OPTIONAL OUTPUTS:
;       -
; 
; PROCEDURES CALLED:
;       REVERSE(),ROT(),SXADDPAR,SXPAR(),SXDELPAR,TOTAL()
;
; EXAMPLE:
;      
;
; MODIFICATION HISTORY:
;       11-07-2018 P.Kamphuis; Replaced bitwise operator not with
;                              logical operator ~
;       07-01-2016 P.Kamphuis; Replaced SUM commands with the proper
;       TOTAL commands reducing the need for outside routines.  
;       05-01-2016 P.Kamphuis; Updated NAXIS1 in new header  
;       Written by P.Kamphuis 01-01-2015 
;
; NOTE:
;     
;-
  COMPILE_OPT IDL2 
  inheader=header
  if n_elements(center) EQ 0 then center=[sxpar(inheader,'CRVAL1'), sxpar(inheader,'CRVAL2')]
  if n_elements(width) EQ 0 then begin
     IF ~(sxpar(inheader,'BMAJ')) then begin
        width=0.
     ENDIF else width=sxpar(inheader,'BMAJ')
  endif    
  xv=dblarr(n_elements(Cube[*,0,0]),n_elements(Cube[0,0,*]))
  ypix=sxpar(inheader,'CRPIX2')+(center[1]-sxpar(inheader,'CRVAL2'))/sxpar(inheader,'CDELT2')
  xpix=sxpar(inheader,'CRPIX1')+(center[0]-sxpar(inheader,'CRVAL1'))/sxpar(inheader,'CDELT1')*COS(center[1]*!DtoR)
  xrange=[-xpix,sxpar(inheader,'NAXIS1')-xpix]
  IF xpix GT sxpar(inheader,'NAXIS1')-xpix-1 then xsize=floor(sxpar(inheader,'NAXIS1')-xpix-1) else xsize=floor(xpix)
  xv=dblarr(fix(2*xsize),n_elements(Cube[0,0,*]))
  newCube=dblarr(n_elements(Cube[*,0,0]),n_elements(Cube[0,*,0]),n_elements(Cube[0,0,*]))
  for j=0,n_elements(Cube[0,0,*])-1 do begin
     Newcube[*,*,j]=ROT(Cube[*,*,j],pa-90.,1.0,xpix,ypix,missing=!values.f_nan,cubic=-1,/PIVOT) ;For HI
  endfor
  If width/(sxpar(inheader,'CDELT2')) LT 1 then begin
     xv[*,*]=Newcube[fix(xpix-xsize):fix(xpix+xsize-1),round(ypix),*]
  ENDIF ELSE BEGIN
     centrpix=xpix
     istart=fix(xpix-xsize)
     
     IF istart lt 0 then begin
        istart=0
     ENDIF
     IF istart Ne 0 then centrpix=xpix-istart
     iend=fix(xpix+xsize-1)
     IF iend-istart GT n_elements(xv[*,0])-1 then iend=n_elements(xv[*,0])-1
     xv[0:iend-istart,*]=TOTAL(Newcube[istart:iend,fix(ypix-(width)/(2.*(sxpar(inheader,'CDELT2')))):fix(ypix+(width)/(2.*(sxpar(inheader,'CDELT2')))),*],2)/n_elements(Newcube[0,fix(ypix-(width)/(2.*(sxpar(inheader,'CDELT2')))):fix(ypix+(width)/(2.*(sxpar(inheader,'CDELT2')))),0])
  ENDELSE
  IF sxpar(inheader,'CDELT1') LT 0 then BEGIN
     sxaddpar,inheader,'CDELT1',ABS(sxpar(inheader,'CDELT1'))
     xv=REVERSE(xv)
  ENDIF

  new_header=inheader
  sxaddpar,new_header,'NAXIS',2
  sxaddpar,new_header,'CDELT1',sxpar(inheader,'CDELT1')
  sxaddpar,new_header,'CRVAL1',0.
  sxaddpar,new_header,'CRPIX1',xsize+1
  sxaddpar,new_header,'CTYPE1','ANGLE'
  sxaddpar,new_header,'NAXIS1',fix(2*xsize)
  sxaddpar,new_header,'CUNIT1','DEGREE',AFTER='CRPIX1'
  sxaddpar,new_header,'PA',pa,AFTER='CTYPE2'
  sxaddpar,new_header,'RA POS',center[0],AFTER='PA'
  sxaddpar,new_header,'DEC POS',center[1],AFTER='RA POS'
  IF width/(sxpar(inheader,'CDELT2')) GT 1 then sxaddpar,new_header,'Strip Width',width,AFTER='PA'
  sxaddpar,new_header,'CDELT2',sxpar(inheader,'CDELT3')
  sxaddpar,new_header,'CRVAL2',sxpar(inheader,'CRVAL3')
  sxaddpar,new_header,'CRPIX2',sxpar(inheader,'CRPIX3')
  sxaddpar,new_header,'CTYPE2',sxpar(inheader,'CTYPE3')
  sxaddpar,new_header,'NAXIS2',sxpar(inheader,'NAXIS3') 
  sxdelpar,new_header,['NAXIS3','CDELT3','CRPIX3','CRVAL3','CTYPE3','LTYPE']
  IF sxpar(inheader,'CUNIT3') then begin
     sxaddpar,new_header,'CUNIT2',sxpar(inheader,'CUNIT3'),after='CTYPE2'
     sxdelpar,new_header,'CUNIT3'
  ENDIF
  IF ~(sxpar(new_header,'CUNIT2')) then begin
     IF sxpar(new_header,'CDELT2') GT 500. then sxaddpar,new_header,'CUNIT2','M/S',after='CTYPE2' else sxaddpar,new_header,'CUNIT2','KM/S',after='CTYPE2'
  ENDIF
  IF STRUPCASE(strtrim(sxpar(new_header,'CUNIT2'),2)) EQ 'M/S' then begin
     print,linenumber()+'EXTRACT_PV: We are converting to KM/S'
     sxaddpar,new_header,'CDELT2',sxpar(new_header,'CDELT2')/1000.
     sxaddpar,new_header,'CRVAL2',sxpar(new_header,'CRVAL2')/1000.
     sxaddpar,new_header,'CUNIT2','KM/S'
  ENDIF
   
end






