Pro columndensity,levels,vsys,beam,vwidth=vwidth,NCOLUMN=ncolumn,DENSITIES=densities,ARCSQUARE=arcsquare,MSOLAR=Msolar

;Program to calculate the column density based on the given flux level 
;    levels=Is the flux levels that are to be transformed to column
;    densities in mJy/beam
;    vsys= is the systemic velocity of the system in km/s
;      beam=the beam of the observations in arcsec if it is 1 parameter then a
;      circular beam is assumed
;      vwidth=the channel width of the observations must be given when
;      flux in mJy/beam if not given it is assumed that the flux level
;      is mJy/beam*km/s
;        /NCOLUMN if set it is assumed that level is column densities
;      if /MSOLAR is set then the output will be in M_solar pc^2     
;IF DENSITIES  is given the result will be published there otherwise
;levels will be converted. IF /NCOLUM is given and densities is given
;it will be assumed that densities is the input and the levels output
;        that have to be transforemed to flux levels in mJy/beam
;    IF /arcsquare is set it is assumed that levels is in mJy/arcsec^2
;    this is specifically useful for converting tirific's SBR to columndensities
;program to calculate the column density
;paolo actually uses n(HI)=1e+21 F(mJy/beam km/s)/bmaj(arcsec)/bmin(arcsec)
  f0=1.420405751786E9           ;Hz rest freq
  
  c=299792.458
  pc=3.086e+18                  ;cm
  solarmass=1.98855e30          ;kg
  mHI=1.6737236e-27             ;neutral hydrogen mass in kg
IF vsys GT 10000 then vsys=vsys/1000.
f=f0*(1-(vsys/c)) ;Systemic frequency
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
                                ;one note that this one gives youck
                                ;column densities in mJy/arcsec^2
      ;if msolar is set we first want to return to column densities
      if keyword_set(msolar) then begin
         levels[*]=levels[*]*solarmass/(mHI*pc^2)
      ENDIF
      IF n_elements(densities) GE 1 then begin
         levels=dblarr(n_elements(densities))
         levels[*]=densities[*]
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
