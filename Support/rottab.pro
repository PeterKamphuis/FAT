Pro rottab,x,y,rotangle,xrotated=xrot,yrotated=yrot,center=center

;+
; NAME:
;       ROTTAB
;
; PURPOSE:
;       Program to rotate tables consisting of x positions and ypositions;      
;
; CATEGORY:
;       Support
; 
; CALLING SEQUENCE:	
;       ROTTAB,x,y,rotangle,xrotated=xrot,yrotated=yrot,center=center
;
; INPUTS:
;         x = Array of x positions
;         y = Array of y positions  
;  rotangle = Required angle of rotation
;  
; OPTIONAL INPUTS:
;        xrotated = array that contains the rotated values. If not
;        present then x will contain the output.
;        yrotated = array that contains the rotated values. If not
;        present then y will contain the output
;        center = Point of rotation, Default 0,0
;  
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;       -
; 
; PROCEDURES CALLED:
;      COS(), SIN()
;
; MODIFICATION HISTORY:
;       Written 01-01-2009 P.Kamphuis v1.0
;
; NOTE:
;     
;-

  COMPILE_OPT IDL2
  IF n_elements(center) GE 1 then begin
     tempx=dblarr(n_elements(x))
     tempx=dblarr(n_elements(y))
     tempx=x-center[0]
     tempy=y-center[1]
  endif else begin
     tempx=dblarr(n_elements(x))
     tempx=dblarr(n_elements(y))
     tempx=x
     tempy=y
  ENDELSE
  tempx2=dblarr(n_elements(x))
  tempy2=dblarr(n_elements(y))
  for i=0,n_elements(tempx)-1 do begin
     tempx2[i]=(COS((rotangle*!pi/180.))*tempx[i]-SIN((rotangle*!pi/180.))*tempy[i])
     tempy2[i]=(SIN((rotangle*!pi/180.))*tempx[i]+COS((rotangle*!pi/180.))*tempy[i])
  endfor
  IF n_elements(center) GE 1 then begin
     tempx2[*]=tempx2[*]+center[0]
     tempy2[*]=tempy2[*]+center[1]
  ENDIF
  IF n_elements(xrot) LT 1 then begin
     x=tempx2
  endif else begin
     xrot=dblarr(n_elements(x))
     xrot=tempx2
  endelse
  IF n_elements(yrot) LT 1 then begin
     y=tempy2
  endif else begin
     yrot=dblarr(n_elements(y))
     yrot=tempy2
  endelse
end


