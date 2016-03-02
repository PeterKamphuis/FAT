pro ra_names, values, tick_value = tv, tick_name = tn, n_ticks = n_ticks, $
              increment = increment, force = force
;+
; NAME:
;   RA_NAMES
; PURPOSE:
;   To generate the locations and names of appropriately places RA
;   ticks on an axis of declination in a graph.
;
; CALLING SEQUENCE:
;   RA_NAMES, values, TICK_VAL = TICK_VAL, TICK_NAME = TICK_NAME
;
; INPUTS:
;   VALUES -- The values of DEC along the axis.
;
; KEYWORD PARAMETERS:
;   N_TICKS -- Forces labelling to this number of ticks.
;   INCREMENT -- Spacing between tickmarks
;   FORCE -- force the given  number of ticks instead of
;            converting to nice values 
;
; OUTPUTS:
;   TICK_VAL -- The value (in decimal degrees) of the location of the ticks.
;   TICK_NAME -- An array of strings naming the ticks in the
;                appropriate format.
; MODIFICATION HISTORY:
;       28-02-2016 P.Kamphuis; Replaced !U with !E for GDL
;       compatibility. Which didn't work.
;       august 2009 Decimals command unknown changed to string with decimals
;       Written in a deeply fatigued state.
;       Wed Jul 24 02:13:37 2002, Erik Rosolowsky <eros@cosmic>
;
; NOTE: P. Kamphuis, GDL does not accept superscripts in the tickmarks
;       of plot. This leads to bad output in the GDL plot.  
;-

; Set up some information about the coordinate axis.
  values = double(values)
  nelts = n_elements(values) 
  start = values[0]
  finish = values[nelts-1]
  range = finish-start

; Establish the array of legitimate values for an increment.
  legit = [6, 4, 2, 1,  1/60d0*[30, 15, 10, 5, 2, 1], $
           1/3.6d3*[30, 15, 10, 5, 2, 1], 1/3.6d4*[5, 2, 1], $
           1/3.6d5*[5, 2, 1], 1/3.6d6*[5, 2, 1]]*15

; Determine the spacing.  If number of ticks not specified, assume
; between 3 and 6 ticks will do.
  if n_elements(n_ticks) gt 0 then begin
    incr = range/n_ticks 
    nticks = n_ticks+2
    if not keyword_set(force) then begin
; Choose ideal increment as largest increment smaller than required
; increment.  Measure spacing logarithmically. 
       diff = alog10(abs(incr))-alog10(legit)
       min_diff = min(diff[where(diff gt 0)])
       ind = where(diff eq min_diff)
       increment = legit[ind]*incr/abs(incr)
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
      increment = min(legit[target/4]*incr/abs(incr))
    endif else begin
      message, 'This condition not yet implemented. I Suck.', /con
      return
    endelse  
  endelse
  increment = increment[0]
  if increment ge 0 then begin
    tv = start+increment*(start gt 0)-(start mod increment)
    while tv[n_elements(tv)-1] lt finish do $
      tv = [tv, tv[n_elements(tv)-1]+increment]
  endif else begin
    tv = start-(start mod increment)
    while tv[n_elements(tv)-1] gt finish do $
      tv = [tv, tv[n_elements(tv)-1]+increment]
  endelse  

; Separate into nice Hour, Min, Sec.  This includes some
; black magic to  avoid irritating floor results with the DOUBLE values.

  tv_hold = tv/15
  ind = where(tv_hold lt 0, ct)
  if ct gt 0 then tv_hold[ind] = tv_hold[ind]+24
  hour = floor(tv_hold)
  minute = floor(double((tv_hold - hour)*60d))
  sec =  float(((double((tv_hold-hour)*60d)-minute)*60))

;Carrying
  c_sec = sec ge 60
  sec = sec-60*c_sec
  minute = minute+c_sec
  c_min = minute ge 60
  minute = minute-60*c_min
  hour = hour+c_min

  if abs(increment/15*3600) lt 1 then $
    decs = ceil(-alog10(abs(increment/15)*3600)) else decs = 0
  a = string("140B)
  b = string("042B)
                                ; Now, the tedium of checking all the
                                ; formatting and only marking
                                ; differences.  

  hour_name = strarr(n_elements(tv))
  min_name = hour_name
  sec_name = hour_name
  hour_name[0] = string(hour[0],  format = '(i2.2)')+"!Eh!N"
  min_name[0] = string(minute[0],  format = '(i2.2)')+"!Em!N"
  sec_name[0] = string(sec[0], format = '(i2.2)')+"."+string(decs, format = '(i1.1)' )+"!Es!N"
  for i = 1, n_elements(tv)-1 do begin
    if hour[i]-hour[i-1] ne 0 then $
      hour_name[i] =  string(hour[i],  format = '(i2.2)')+"!Eh!N"
    if minute[i]-minute[i-1] ne 0 then $
      min_name[i] = string(minute[i],  format = '(i2.2)')+"!Em!N"
    if sec[i]-sec[i-1] ne 0 then $
;      sec_name[i] = decimals(sec[i], decs)+"!Us!N"
       sec_name[i] = string(sec[i], format = '(i2.2)')+"."+string(decs, format = '(i1.1)' )+"!Es!N"

  endfor
  if total(sec) eq 0 then sec_name = strarr(n_elements(tv))
  if total(minute) eq 0 then min_name = strarr(n_elements(tv))
  if total(sec) eq 0 then  tn = hour_name+min_name else tn = hour_name+min_name+sec_name
  return
end

