pro dec_names, values, tick_value = tv, tick_name = tn, n_ticks = n_ticks, $
               increment = increment,force = force
;+
; NAME:
;   DEC_NAMES
; PURPOSE:
;   To generate the locations and names of appropriately places DEC
;   ticks on an axis of declination in a graph.
;
; CALLING SEQUENCE:
;   DEC_NAMES, values, TICK_VAL = TICK_VAL, TICK_NAME = TICK_NAME
;
; INPUTS:
;   VALUES -- The values of DEC along the axis.
;   FORCE -- force the given tick marks instead of converting to nice values
;
; KEYWORD PARAMETERS:
;   N_TICKS -- Forces labelling to this number of ticks.
;   INCREMENT -- Spacing between tick marks
;   
; OUTPUTS:
;   TICK_VAL -- The value (in decimal degrees) of the location of the ticks.
;   TICK_NAME -- An array of strings naming the ticks in the
;                appropriate format.
; MODIFICATION HISTORY:
;       18-07-2011
;       Modified to correctly produce the degree symbol when using
;       True Type and postscript fonts 
;       P.Kamphuis
;
;       15-09-2009
;       Modified to correctly handle southern declinations.
;       P.Kamphuis
;
;       Modified to correctly handle southern declinations.????
;       Mon Mar 3, 2003, Nate McCrady <nate@astro>
;
;       Written in a deeply fatigued state.
;       Wed Jul 24 02:13:37 2002, Erik Rosolowsky <eros@cosmic>
;
; NOTE: P. Kamphuis, GDL does not recognize the degree symbol not does it accept
;       superscripts in the tickmarks of plot. This leads to bad output in
;       the GDL plot.
;  -

; Set up some information about the coordinate axis.
  values = double(values)
  nelts = n_elements(values) 
  start = values[0]
  finish = values[nelts-1]
  range = finish-start
; Establish the array of legitimate values for an increment.
  legit = [[60, 30, 15, 10, 5, 2, 1], 1/60d0*[30, 15, 10, 5, 2, 1], $
           1/3.6d3*[30, 15, 10, 5, 2, 1], 1/3.6d4*[5, 2, 1], $
           1/3.6d5*[5, 2, 1], 1/3.6d6*[5, 2, 1]]

; Determine the spacing.  If number of ticks not specified, assume
; between 3 and 6 ticks will do.
  if n_elements(n_ticks) gt 0 then begin
    incr = range/n_ticks 
    nticks = n_ticks+2
    if not keyword_set(force) then begin
; Choose ideal increment as largest increment smaller than required
; increment.  Measure spacing logarithmically. 
       diff = alog10(incr[0])-alog10(legit)
       min_diff = min(diff[where(diff gt 0)])
       ind = where(diff eq min_diff)
       increment = legit[ind]
    endif else increment = incr
  endif else begin    
; As above but choose number of ticks as well as the ideal increment.
    incr = range/(dindgen(4)+3)
    legit_array = (fltarr(4)+1) # legit
    incr_array = incr # (fltarr(n_elements(legit))+1) 
    diff = alog10(abs(incr_array))-alog10(legit_array)

    ind = where(diff gt 0, ct)

    if ct gt 0 then begin 
      min_diff = min(diff[ind])
      target = where(diff eq min_diff)
      nticks = max((target mod 4)+3)+2
      increment = min(legit[target/4])*incr/abs(incr)
    endif else begin
      message, 'This condition not yet implemented. I Suck.', /con
      return
    endelse  
  endelse
  increment = increment[0]

  if increment ge 0 then begin
    tv = start+increment*(start gt 0 and increment gt 0)-(start mod increment)
    while tv[n_elements(tv)-1] lt finish do $
      tv = [tv, tv[n_elements(tv)-1]+increment]
  endif else begin
    tv = start+increment*(start gt 0 and increment gt 0)-(start mod increment)
    while tv[n_elements(tv)-1] gt finish do $
      tv = [tv, tv[n_elements(tv)-1]+increment]
  endelse
;  tv = dindgen(nticks)*increment+start+increment*(start mod increment gt 0)-(start mod increment)

                     ; Separate into nice Deg, Min, Sec.  This includes some
                     ; black magic to avoid irritating floor results with the DOUBLE values.
  
  deg = floor(tv)*(tv ge 0) + ceil(tv)*(tv lt 0)
  minute = floor(float((tv - deg)*60))*(tv ge 0)+$
    (tv lt 0)*floor(float(abs((tv - deg))*60))
  sec = (double((double((tv-deg)*60)-minute)*60))*(tv ge 0)+$
    abs((double((double((tv-deg)*60)+minute)*60)))*(tv lt 0)

  sec = double(sec)
; Carrying
  
  c_sec = sec ge 60
  sec = sec-60.*c_sec
  minute = minute+c_sec
  c_min = minute ge 60
  minute = minute-60*c_min
  deg = deg+c_min
  if abs(increment*3600) lt 1 then $
    decs = ceil(-alog10(abs(increment)*3600)) else decs = 0  

  a = string("047B)
  b = string("042B)
                                ; Now, the tedium of checking all the
                                ; formatting and only marking
                                ; differences.  

  deg_name = strarr(n_elements(tv))
  min_name = deg_name
  sec_name = deg_name  
  IF !P.FONT EQ -1 then begin
     degreesymbol='!9%!X'
  endif else begin
     degreesymbol='!9'+string("260B)+'!X'
  endelse
  if tv[0] LT 0. then begin
     if fix(tv[0]) EQ 0 then begin
        deg_name[0] = '-'+strcompress(string(deg[0]), /rem)+degreesymbol
     endif else begin
        deg_name[0] = strcompress(string(deg[0]), /rem)+degreesymbol
     endelse
  endif else begin
     deg_name[0] =  strcompress(string(deg[0]), /rem)+degreesymbol
  endelse
;  deg_name[0] = strcompress(string(deg[0]), /rem)+"!9%!X"
  min_name[0] = string(minute[0],  format = '(i2.2)')+a
;  sec_name[0] = decimals(sec[0], decs)+b 
  sec_name[0] = string(sec[0], format = '(i2.2)')+"."+string(decs, format = '(i1.1)' )+b
  tracksec=0.
  for i = 1, n_elements(tv)-1 do begin
     if deg[i]-deg[i-1] ne 0 then begin
        if tv[i] LT 0. then begin
           if fix(tv[i]) EQ 0 then begin
              deg_name[i] = '-'+strcompress(string(deg[i]), /rem)+degreesymbol
           endif else begin
              deg_name[i] = strcompress(string(deg[i]), /rem)+degreesymbol
           endelse
       endif else begin
           deg_name[i] =  strcompress(string(deg[i]), /rem)+degreesymbol
        endelse
     endif

    if minute[i]-minute[i-1] ne 0 then $
      min_name[i] = string(minute[i],  format = '(i2.2)')+a
    if sec[i]-sec[i-1] ne 0 then $
;      sec_name[i] = decimals(sec[i], decs)+b
      sec_name[i] = string(sec[i], format = '(i2.2)')+"."+string(decs, format = '(i1.1)' )+b
    tracksec=tracksec+double(string(sec[i], format = '(i2.2)')+"."+string(decs, format = '(i1.1)' ))
 endfor

  if total(sec) eq 0 then sec_name = strarr(n_elements(tv))
  if total(minute) eq 0 then min_name = strarr(n_elements(tv))
  if tracksec eq 0. then tn = deg_name+min_name else tn = deg_name+min_name+sec_name
  return
end
