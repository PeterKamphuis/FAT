PRO fat_ploterror,plotradii,plotVariable,xerr,yerr,ERRCOLOR = errcolor, ERRTHICK= errthick, _EXTRA = extra,over_plot=oplot, HATLENGTH= hatlength

;+
; NAME:
;       FAT_PLOTERROR
;
; PURPOSE:
;       plot a set of variables with errors
;
; CATEGORY:
;       Support
; 
; CALLING SEQUENCE:
;       FAT_PLOTERROR,plotradii,plotVariable,xerr,yerr,ERRCOLOR = errcolor, ERRTHICK= errthick, _EXTRA = extra,over_plot=oplot
;
;
; INPUTS:
;         PLOTRADII = Distance to the galaxy
;         PLOTVARIABLE = an indicator whether we are running idl or
;         XERR= Error in the x direction
;         YERR=Error in the Y direction
;  
; OPTIONAL INPUTS:
;    ERRCOLOR =  color of errors
;    ERRTHICK=  Thickness of errors
; 
; KEYWORD PARAMETERS:
;       -
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;       -
; 
; PROCEDURES CALLED:
;        OPLOT, PLOT
;
; MODIFICATION HISTORY:
;  
;       Written 01-06-2016 P.Kamphuis v1.0
;
; NOTE: Written because in the NASA library ploterror they suddenly
; introduced the setdefault values from the fanning library which woul
; require more libraries to be installed. Also helps with the move to
; GDL.
;     
;-

  compile_opt idl2

  IF n_elements(errcolor) EQ 0 then errcolor=0
  IF n_elements(errthick) EQ 0 then errthick=!p.thick
  IF keyword_set(oplot) then begin
     oplot,plotradii,plotVariable,_EXTRA=extra
  endif else begin
     plot,plotradii,plotVariable,_EXTRA=extra
  endelse
  IF n_elements(hatlength) EQ 0 then hatlength=[(!Y.CRANGE[1]-!Y.CRANGE[0])/100.,(!X.CRANGE[1]-!X.CRANGE[0])/100.]
  for i=0,n_elements(plotVariable)-1 do begin
     IF xerr[i] NE 0 then begin
        oplot,[plotradii[i]-xerr[i],plotradii[i]+xerr[i]],[PlotVariable[i],PlotVariable[i]],thick=errthick,color=errcolor
                                ;plot a hat
        oplot,[plotradii[i]-xerr[i],plotradii[i]-xerr[i]],[PlotVariable[i]-hatlength[0],PlotVariable[i]+hatlength[0]],thick=errthick,color=errcolor
        oplot,[plotradii[i]+xerr[i],plotradii[i]+xerr[i]],[PlotVariable[i]-hatlength[0],PlotVariable[i]+hatlength[0]],thick=errthick,color=errcolor
     ENDIF
     IF yerr[i] NE 0 then begin       
        oplot,[plotradii[i],plotradii[i]],[PlotVariable[i]-yerr[i],PlotVariable[i]+yerr[i]],thick=errthick,color=errcolor
                                ;plot a hat
        oplot,[plotradii[i]-hatlength[1],plotradii[i]+hatlength[1]],[PlotVariable[i]-yerr[i],PlotVariable[i]-yerr[i]],thick=errthick,color=errcolor
        oplot,[plotradii[i]-hatlength[1],plotradii[i]+hatlength[1]],[PlotVariable[i]+yerr[i],PlotVariable[i]+yerr[i]],thick=errthick,color=errcolor
     ENDIF
  endfor
end

  
