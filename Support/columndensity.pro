Pro columndensity,levels,vsys,beam,VWIDTH=vwidth,NCOLUMN=ncolumn,DENSITIES=densities,ARCSQUARE=arcsquare,MSOLAR=Msolar

;+
; NAME:
;       COLUMNDENSITY
;
; PURPOSE:
;       Program to calculate the column density based on the given flux level 
;
; CATEGORY:
;       Support
; 
; CALLING SEQUENCE:
;       COLUMNDENSITY,levels,vsys,beam,VWIDTH=vwidth,/NCOLUMN,DENSITIES=densities,/ARCSQUARE
;
;
; INPUTS:
;       levels = Is an array or scalar with the flux levels (mJy/beam) that are
;       to be transformed to column densities
;       vsys = is scalar with the systemic velocity of the system in km/s
;       beam = the beam of the observations in arcsec if it is 1 parameter then a
;       circular beam is assumed
;
; OPTIONAL INPUTS:
;       - 
;
; KEYWORD PARAMETERS:
;       vwidth = the channel width of the observations must be given when
;       flux is in mJy/beam if not given it is assumed that the flux level
;       is mJy/beam*km/s
;       /NCOLUMN - if set it is assumed that level is column densities
;       /ARCSQUARE - if set it is assumed that levels is in
;                    mJy/arcsec^2 this is specifically useful
;                    for converting tirific's SBR to
;                    columndensities program to calculate the column
;                    density
;        /MSOLAR - Give output in M_solar pc^2,  if ncolumn is set it
;                  is assumed the input is in M_solar pc^2  
;
; OUTPUTS:
;       levels = the transformed units if densities is unset
;
; OPTIONAL OUTPUTS:
;       DENSITIES = if given the result will be published there
;       otherwise levels will be converted. IF /NCOLUM is given and
;       densities is given it will be assumed that densities is the
;       input and the levels output that have to be transforemed to flux levels in mJy/beam
; 
; PROCEDURES CALLED:
;       -
;
; EXAMPLE:
;      
;
; MODIFICATION HISTORY:
;       16-12-2016 P.Kamphuis; Added the option to have the levels
;       converted to Msolar pc^2  
;       Written by P.Kamphuis 01-01-2015 
;
; NOTE:
;      A more rough shod way would be n(HI)=1e+21 F(mJy/beam km/s)/bmaj(arcsec)/bmin(arcsec)
;-
  COMPILE_OPT IDL2 
  f0=1.420405751786E9           ;Hz rest freq
  c=299792.458                  ;light speed in km/s 
  pc=3.086e+18                  ;parsec in cm
  solarmass=1.98855e30          ;Solar mass in kg
  mHI=1.6737236e-27             ;neutral hydrogen mass in kg
  IF vsys GT 10000 then vsys=vsys/1000.
  f=f0*(1-(vsys/c))             ;Systemic frequency
  IF keyword_set(arcsquare) then begin
                                ;we do not want to use the beam but
                                ;need to correct for that with a
                                ;factor the difference is a factor
                                ;1.1330900354567984..., which comes
                                ;from the normalisation to the unit
                                ;"beam". It is the Integral of a 2D
                                ;Gaussian with unity HPBW. From
                                ;Josh's e-mail.
     HIconv=605.7383*1.823E18*(2.*!pi/(ALOG(256.))) 
     IF keyword_set(ncolumn) then begin
                                ;I can't find the redshift dependence
                                ;so let's leave it out in this
                                ;one note that this one gives you
                                ;column densities in mJy/arcsec^2
                             
        IF n_elements(densities) GE 1 then begin
           levels=dblarr(n_elements(densities))
           levels[*]=densities[*]
        ENDIF
                                ;if msolar is set we want to convert
                                ;back to column densities
        if keyword_set(msolar) then begin
           levels[*]=levels[*]*solarmass/(mHI*pc^2)
        ENDIF
        IF n_elements(vwidth) NE 0. then  levels[*]=levels[*]/(HIconv*vwidth) else  levels[*]=levels[*]/(HIconv)
     ENDIF ELSE BEGIN
        NH=dblarr(n_elements(levels))
        IF n_elements(vwidth) NE 0. then NH=HIconv*levels*vwidth else NH=HIconv*levels

        IF n_elements(densities) GE 1 then begin
           densities=dblarr(n_elements(levels))
           densities[*]=NH[*]
        ENDIF else begin
           levels[*]=NH[*]
        endelse
     ENDELSE
     
  endif else begin
     IF n_elements(beam) EQ 1 then begin
        tmp=beam
        beam=dblarr(2)
        beam[*]=tmp
     ENDIF
     b=beam[0]*beam[1]
     IF keyword_set(ncolumn) then begin
        IF keyword_set(densities) then begin
           levels=dblarr(n_elements(densities))
           levels[*]=densities[*]
        ENDIF 
                                ;if msolar is set we want to convert
                                ;back to column densities
        if keyword_set(msolar) then begin
           levels[*]=levels[*]*solarmass/(mHI*pc^2)
        ENDIF
      
        TK=dblarr(n_elements(levels))
        IF n_elements(vwidth) NE  0 then TK=levels[*]/(1.823E18*vwidth) else  TK=levels[*]/(1.823E18)
        levels[*]=TK[*]/(((605.7383)/(b))*(f0/f)^2)  
     endif else begin   
        TK=dblarr(n_elements(levels))
        TK=((605.7383)/(b))*(f0/f)^2*levels[*]
        NH=dblarr(n_elements(levels))
        IF n_elements(vwidth) NE 0. then NH=1.823E18*Tk*vwidth else NH=1.823E18*Tk
        IF n_elements(densities) GE 1 then begin
           densities=dblarr(n_elements(levels))
           densities[*]=NH[*]
        ENDIF else begin
           levels[*]=NH[*]
        endelse
     endelse
  endelse
  
  IF ~keyword_set(ncolumn)  and keyword_set(msolar) then begin
     levels[*]=levels[*]*mHI*pc^2/solarmass
  ENDIF
end
