Pro Interpolate,Values,radii,output=output,newradii=newradii


;+
; NAME:
;       INTERPOLATE
;
; PURPOSE:
;       Routine to interpolate values to a new set of nodes
;
; CATEGORY:
;       Support
; 
; CALLING SEQUENCE:
;       INTERPOLATE,Values,radii,OUTPUT=output,NEWRADII=newradii
;
;
; INPUTS:
;       Values = the original values 
;       radii =  the original nodes 
;
; OPTIONAL INPUTS:
;       - 
;
; KEYWORD PARAMETERS:
;       OUTPUT = the values at the new node locations
;       NEWRADII = the new node locations. If unset it will be half
;       the original nodes
;
; OUTPUTS:
;       -
;
; OPTIONAL OUTPUTS:
;       -
; 
; PROCEDURES CALLED:
;       -
;
; EXAMPLE:
;      
;
; MODIFICATION HISTORY:
;       Written by P.Kamphuis 01-01-2015 
;
; NOTE:
;     
;-
  COMPILE_OPT IDL2 

  IF  n_elements(radii) LT 2 OR n_elements(radii) NE n_elements(Values) then begin
     print,'INTERPOLATE: Please provide at least two values and make sure that x and y values have the same amount of elements'
     goto,breakthis
  ENDIF
                                ;let's make the new radii half the old if not given
  halved=0.
  IF n_elements(newradii) EQ 0. then begin
     newradii=findgen( n_elements(radii)*2)*(radii[1]-radii[0])/2.+radii[0]
     IF n_elements(radii) GT 2 then halved=1. else halved=2
  ENDIF
  addzerorad=0.
  IF n_elements(newradii) EQ 1 then halved=2.
  IF radii[0] NE newradii[0] then begin
     tmp=newradii
     newradii=dblarr(n_elements(newradii)+1.)
     newradii=[radii[0],tmp]
     addzerorad=1.
  ENDIF
  
  newparameter=dblarr(n_elements(newradii))
  newparameter[0]=Values[0]
  for i=1,n_elements(newradii)-1 do begin
     x3=newradii[i]
     case 1 of
        x3 GT radii[n_elements(radii)-1]: begin
           y1=Values[n_elements(radii)-2]
           y2=Values[n_elements(radii)-1]    
           x1=radii[n_elements(radii)-2]
           x2=radii[n_elements(radii)-1]
           newparameter[i]=y1-(((x3-x1)*(y1-y2))/(x2-x1))
           
        end
        
        x3 LT radii[0] :begin
           y1=Values[0]
           y2=Values[1]    
           x1=radii[0]
           x2=radii[1]
           newparameter[i]=y1-(((x3-x1)*(y1-y2))/(x2-x1))
        end
        else:begin
           for j=0,n_elements(radii)-1 do begin
              if radii[n_elements(radii)-1-j] LE x3 then begin
                 
                 pos1=n_elements(radii)-1-j
                 goto,found1
              endif
           endfor
           found1:
           
           for j=0,n_elements(radii)-1 do begin
              if radii[j] GE x3 then begin
                 pos2=j
                 goto,found2
              endif
           endfor
           found2:
           y1=Values[pos1]
           y2=Values[pos2]    
           x1=radii[pos1]
           x2=radii[pos2]
           
           newparameter[i]=y1-(((x3-x1)*(y1-y2))/(x2-x1))
           
           IF x1-x2 EQ 0. then newparameter[i]=Values[pos1]
        end
     endcase
  endfor
  

                                ;Cleaning up 
  IF addzerorad EQ 1. then begin
     tmp=newradii
     tmp2=newparameter
     newradii=dblarr(n_elements(tmp)-1)
     newparameter=dblarr(n_elements(tmp2)-1)
     newradii=tmp[1:n_elements(tmp)-1]
     newparameter=tmp2[1:n_elements(tmp2)-1]
  ENDIF
  maxx=MAX([radii,newradii],Min=minx)
  maxy=MAX([newparameter,Values],Min=miny)

  case halved of
     0:begin
        IF n_elements(output) EQ 0 then begin
           Values=dblarr(n_elements(newparameter))
           Values=newparameter
        ENDIF ELSE begin
           output=dblarr(n_elements(newparameter))
           output=newparameter
        ENDELSE
     end
     1:begin
        radii=dblarr(n_elements(newradii))
        radii=newradii
        IF n_elements(output) EQ 0 then begin
           Values=dblarr(n_elements(newparameter))
           Values=newparameter
        ENDIF ELSE begin
           output=dblarr(n_elements(newparameter))
           output=newparameter
        ENDELSE
     end
     2:begin
        radii=dblarr(1)
        radii=newradii[0]
        IF n_elements(output) EQ 0 then begin
           Values=dblarr(1)
           Values=newparameter[0]
        ENDIF ELSE begin
           output=dblarr(1)
           output=newparameter[0]
        ENDELSE
     end
     else:begin
     end
  endcase

  
  
  breakthis:
end
