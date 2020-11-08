function fat_spline,x,y,x_new

;+
; NAME:
;       fat_spline
;
; PURPOSE:
;       Perform an akima spline interpolation.
;
; CATEGORY:
;       Support
; 
; CALLING SEQUENCE:
;       FAT_SPLINE(x,y,new_x)	
;
;
; INPUTS:
;         x = xaxis
;         y = yaxis
;     new_x = the new xaxis
; OPTIONAL INPUTS:
;  
; KEYWORD PARAMETERS:
;  
; OUTPUTS:
;    new_y = the new interpolated values for new_x
; OPTIONAL OUTPUTS:
;       -
; 
; PROCEDURES CALLED:
;
; MODIFICATION HISTORY:
;       Written 27-08-2019 P.Kamphuis v1.0
;
; NOTE:
;     Not sure how wel this works 
;-


  y_new=dblarr(n_elements(x_new))
  if n_elements(x) LT 3 then begin
     print,linenumber()+"FAT_SPLINE: The input array is too small."
     stop
  endif
  if n_elements(x) NE n_elements(y) then begin
     print,linenumber()+"FAT_SPLINE: The input arrays must have the same size."
     stop
  endif
  dx=dblarr(n_elements(x)-1)
  val1=dblarr(n_elements(x)-1)
  for i=0,n_elements(x)-2 do begin
     if x[i] GT x[i+1] then begin
        print,linenumber()+"FAT_SPLINE: the axis must be increasing."
        stop
     endif
     dx[i]=x[i+1]-x[i]
     val1[i] = (y[i+1]-y[i])/dx[i]
  endfor
 
  tmp =WHERE(x_new LT x[0] OR x_new GT x[n_elements(x)-1])
  if tmp[0] NE -1 then begin
     print,linenumber()+"FAT_SPLINE: the new axis cannot be larger than the input axis."
     stop
  endif

  val2=2.*val1[0]-val1[1]
  val3 = 2.*val2 - val1[0]
  endval1= 2.*val1[n_elements(val1)-2]-val1[n_elements(val1)-3]
  endval2= 2.*endval1-val1[n_elements(val1)-2]

  newval1=[val3,val2,val1,endval1,endval2]
  absnewval1=dblarr(n_elements(newval1)-1)
  for i=0,n_elements(newval1)-2 do begin
     absnewval1[i]=abs(newval1[i+1]-newval1[i])
  endfor

  g1 = absnewval1[2:n_elements(absnewval)-1]
  g2 = absnewval1[0:n_elements(x)]
  g12 = g1+g2

  inds = WHERE(g12 > 1e-9*max(g12))
  d1 = newval1[1:n_elements(newval1)-2]

  d1[inds]=(g1[inds]*newval1[inds+1]+g2[inds]*newval1[inds+2])/g12[inds]
  r = (3.*val1-2.0*d1[0:n_elements(x)-1]-d1[1:n_elements(x)])/dx
  klm = (d1[0:n_elements(x)-1]+d1[1:n_elements(x)]-2.*val1)/dx^2
  bins = dblarr(n_elements(x_new))
  counter =0
  for i=0,n_elements(x_new)-1 do begin
     IF x_new[i] LT x[counter+1] then bins[i] = counter else begin
        counter++
        bins[i]=counter
     endelse
  endfor
  luft= x_new-x[bins]
  y_new= (luft*klm[bins]+r[bins])*luft+d1[bins]*luft+y[bins]
  print,y_new
  print,'Got to end'
  return,y_new
end


  
  
  
  
