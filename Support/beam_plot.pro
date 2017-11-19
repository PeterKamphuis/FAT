Pro beam_plot, BMAJ,BMIN,BPA=bpa,center=center,Fill=Fill,_EXTRA=ex,Transparent=transparent

;+
; NAME:
;       BEAM_PLOT
;
; PURPOSE:
;       Program to plot a beam based on the FWHM's
;
; CATEGORY:
;       Support
; 
; CALLING SEQUENCE:	
;       BEAM_PLOT, BMAJ,BMIN,BPA=bpa,center=center,Fill=Fill,_EXTRA=ex,Transparent=transparent 
;
; INPUTS:
;         BMAJ = Beam major axis FWHM
;         BMIN = Beam minor axis FWHM
; OPTIONAL INPUTS:
;         center = Central position of the beam. Default 0,0
;         BPA= beam PA. Default =0
;  
; KEYWORD PARAMETERS:
;         FILL - If set the beam will be filled by lines under an angle       -
;      Transparent - IF set the beam will be overplotted
;                    directly. Normally a white background is drawn
;                    below the beam
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;       -
; 
; PROCEDURES CALLED:
;       FINDGEN(), COS(), SIN(), MAX(), POLYFILL, ROTTAB TVLCT
;
; MODIFICATION HISTORY:
;       Written 01-01-2009 P.Kamphuis v1.0
;
; NOTE:
;     
;-

  COMPILE_OPT IDL2
  if n_elements(center) EQ 0 then center=[0,0]
  if n_elements(bpa) EQ 0 then bpa=0.
  A = FINDGEN(100) * (!PI*2/99.)
  xcoord=(COS(A)*(BMIN/2.))/3600.
  ycoord=(SIN(A)*(BMAJ/2.))/3600.
  rottab,xcoord,ycoord,bpa

  ycoordmax=MAX(ycoord+center[1],min=ycoordmin)
  xcoordmax=MAX(xcoord/cos((center[1]/360.)*2.*!pi)+center[0],min=xcoordmin)
  bufferx=(xcoordmax-xcoordmin)/5.
  buffery=(ycoordmax-ycoordmin)/5.

  IF ~keyword_set(transparent) then begin
     TVLCT, Rin, Gin, Bin, /GET
     loadct,0,/SILENT
     POLYFILL,[xcoordmin-bufferx,xcoordmin-bufferx,xcoordmax+bufferx,xcoordmax+bufferx,xcoordmin-bufferx], [ycoordmin-buffery,ycoordmax+buffery,ycoordmax+buffery,ycoordmin-buffery,ycoordmin-buffery],color=254
     TVLCT, Rin, Gin, Bin
  endif
  oplot,xcoord/cos((center[1]/360.)*2.*!pi)+center[0],ycoord+center[1], _EXTRA=ex
  IF keyword_set(FILL) then begin
     POLYFILL,xcoord/cos((center[1]/360.)*2.*!pi)+center[0], ycoord +center[1],/LINE_FILL,ORIENTATION=45,spacing=0.1 , _EXTRA=ex
  ENDIF
end

