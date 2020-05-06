Pro buildaxii,header,xaxis,yaxis,ZAXIS=zaxis

;+
; NAME:
;       BUILDAXII
;
; PURPOSE:
;       Program to build the axii of a fits file plot just give the header of the
;       fits file and it will build the axii (2 or 3 dimensions). No rotation included.
;
; CATEGORY:
;       Support
;
; CALLING SEQUENCE:
;       BUILDAXII,header,xaxis,yaxis,ZAXIS=zaxis
;
;
; INPUTS:
;       header = header of the fits file for which the axii need to be
;       build
;
; OPTIONAL INPUTS:
;       -
;
; KEYWORD PARAMETERS:
;       -
;
; OUTPUTS:
;       xaxis = the variable that will hold the output x-axis
;       yaxis = the variable that will hold the output y-axis
;
;
; OPTIONAL OUTPUTS:
;       ZAXIS= the variable that will hold the zaxis if a zaxis is required
; PROCEDURES CALLED:
;       SXPAR()
;
;
; EXAMPLE:
;
;
; MODIFICATION HISTORY:
;     8-2-2011 modified so that it now builds a proper RA axis if present
;     18-11-2010 added axis from cd matrix
;     Written by P.Kamphuis 01-01-2010
; NOTE:
;     BUILDAXII will always build 3 axii in case of a 3 dimensional
;     cube and 2 in case of a 2 dimensional cube. It is just optional
;     whether you want the third axis returned.
;-
  COMPILE_OPT IDL2
  CATCH,Error_status
  IF  Error_status NE 0. THEN BEGIN
     print, 'BUILDAXII: Oops the following went wrong:'
     print, !ERROR_STATE.MSG
     print, 'BUILDAXII: Use buildaxii in this way'
     print, 'BUILDAXII: Build the axii based on a header red from a fits file   '
     print, 'BUILDAXII: CALLING SEQUENCE: buildaxii,header,xaxis,yaxis,ZAXIS=zaxis'
     print, 'BUILDAXII:    Keywords:'
     print, 'BUILDAXII:     header       array with the header keywords from the fits file'
     print, 'BUILDAXII:      xaxis       name of the first axis '
     print, 'BUILDAXII:      yaxis       name of the second axis '
     print, 'BUILDAXII:      ZAXIS=      name of a third axis '
     goto,ending
  ENDIF
  numaxii=sxpar(header,'NAXIS')
  xsize=sxpar(header,'NAXIS1')
  cdelt1=sxpar(header,'CDELT1')
  ;crval1=sxpar(header,'CRVAl1')
  crpix1=sxpar(header,'CRPIX1')
  ysize=sxpar(header,'NAXIS2')
  ;crval2=sxpar(header,'CRVAl2')
  crpix2=sxpar(header,'CRPIX2')
  cdelt2=sxpar(header,'CDELT2')
  ctype1=sxpar(header,'CTYPE1')
  ctype2=sxpar(header,'CTYPE2')
  ;As we do not do projections accurately let's make sure our crpix is set to the accurate center of the image
  crpix1=xsize/2.
  crpix2=ysize/2.
  xyad,header,crpix1,crpix2,crval1,crval2
  tmp=str_sep(strtrim(strcompress(ctype1),2),'---')
  ctype1=tmp[0]
  tmp=str_sep(strtrim(strcompress(ctype2),2),'---')
  ctype2=tmp[0]
  IF cdelt1 EQ 0 OR cdelt2 EQ 0 then begin
     cd=dblarr(2,2)
     cd[0,0]=sxpar(header,'CD1_1')
     cd[0,1]=sxpar(header,'CD1_2')
     cd[1,0]=sxpar(header,'CD2_1')
     cd[1,1]=sxpar(header,'CD2_2')
     cdelt1=cd[0,0]/(COS(ATAN(cd[1,0]/cd[0,0])))
     cdelt2=cd[1,1]/(COS(ATAN(-cd[0,1]/cd[1,1])))
     CROTA2=ATAN(cd[1,0]/cd[0,0])*(180./!pi)
     IF ABS(crota2) GT 2 then begin
        print,'BUILDAXII: You have significant rotation (crota2='+string(crota2)+'from the cd matrix. Buildaxii is not made for that aborting'
        goto,ending
     endif
  endif

  IF numaxii GE 3. then begin
     zsize=sxpar(header,'NAXIS3')
     cdelt3=sxpar(header,'CDELT3')
     crval3=sxpar(header,'CRVAl3')
     crpix3=sxpar(header,'CRPIX3')
     ctype3=sxpar(header,'CTYPE3')
     cunit3=sxpar(header,'CUNIT3')
  endif

  xaxis=dblarr(xsize)
  bla=dblarr(xsize)
  bla=findgen(xsize)
  IF strupcase(ctype1) EQ 'RA' then begin
     crval2rad=(crval2/360.)*2.*!pi
     cdelt1=(cdelt1/(cos(crval2rad)))
  ENDIF
  xaxis[*]=crval1+(bla[*]-crpix1+1)*(cdelt1)
  yaxis=dblarr(ysize)
  bla=dblarr(ysize)
  bla=findgen(ysize)
  IF strupcase(ctype2) EQ 'RA' then begin
     crval1rad=(crval1/360.)*2.*!pi
     cdelt2=(cdelt2/(cos(crval1rad)))
  ENDIF
  yaxis[*]=crval2+(bla[*]-crpix2+1)*(cdelt2)
  IF numaxii GE 3. then begin
     zaxis=dblarr(zsize)
     bla=dblarr(zsize)
     bla=findgen(zsize)
     zaxis[*]=crval3+(bla[*]-crpix3+1)*(cdelt3)
  endif

ending:
end
