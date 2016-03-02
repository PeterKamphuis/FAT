PRO colour_bar, x, y, minval, maxval, TITLE=title, VERTICAL=dummy1,$
                    BLACK=blackw,COLOR=legcolor,_EXTRA=ex,MULTIPLE_TICKS=multiple_ticks,OPPOSITE_LABEL=opposite_label, $
                xtickname=xtn,ytickname=ytn,ytickformat=ytf,xtickformat=xtf,supress=supress
  compile_opt idl2
  ;+
; NAME:
;       COLOUR_BAR
;
; PURPOSE:
;       Program to make colour bars, Any size or shape or position
;
; CATEGORY:
;       Support
; 
; CALLING SEQUENCE:
;       colour_bar, x, y, minval, maxval, TITLE=title, VERTICAL=dummy1,
;                   BLACK=blackw,COLOR=legcolor,_EXTRA=ex,MULTIPLE_TICKS=multiple_ticks,
;                   OPPOSITE_LABEL=opposite_label,
;                   xtickname=xtn,ytickname=ytn,
;                   ytickformat=ytf,xtickformat=xtf,supress=supress
;
;
; INPUTS:
;               x = horizontal positions of the plot 
;               y = vertical positions of the plot
;          minval = minimum value of color range
;          maxval = maximum value of color range
;  
; OPTIONAL INPUTS:
;           COLOR = Color to be used for labelling
;          _EXTRA = Any plot keywords
;           TITLE = Units of the color table
;       xtickname = name of x tickmarks
;     xtickformat = format of x tickmarks
;       ytickname = name of y tickmarks
;     ytickformat = format of y tickmarks 
; 
; KEYWORD PARAMETERS:
;          /BLACK - Set this keyword to ensure black outlining
;                   regardless of the colour table.
; /MULTIPLE_TICKS - Have multiple ticks on your color bar instead of
;                   just minimum and maximum  
; /OPPOSITE_LABEL - Puts the label on the opposite side of the bar
;                   (i.e. on top or right side)
;        /SUPRESS - Supresses labeling?
;       /VERTICAL - Create a vertical colour bar instead of a
;                   horizontal one.;  
; OUTPUTS:
;       Result = the smoothed map 
;
; OPTIONAL OUTPUTS:
;       -
;
; PROCEDURES CALLED:
;       AXIS,LOADCT,PLOT,POLYFILL,TVLCT,XYOUTS
;
; EXAMPLE:
;      
;
; MODIFICATION HISTORY: 
;       8-10-2009 P.Kamphuis;  Minimum and maximum value have to be given as a string or plotted
;                 with 8 digits. Give the string length yourself
;       Written 22-09-2009 by P.Kamphuis
;
; NOTE:
;  
;-

  CATCH,Error_status
  IF  Error_status NE 0. THEN BEGIN
     print, 'Oops the following went wrong:'
     print, !ERROR_STATE.MSG
     print, 'Use Colour_bar in this way:'
     print,"CALLING SEQUENCE:colour_bar,x,y,'minval','maxval', EXTRA=extra"
     print, '              x = the left and right positions of the colourbar in page units'
     print, '              y = the lower and upper positions of the colourbar in page units'
     print, '         minval = the minimum value of the plotted color'
     print, '         maxval = the maximum value of the plotted color'
     print,' Give these as strings for perfect control on what is plotted'
     print, '          TITLE = title of the colorbar'
     print, '       /VERTICAL= if set the colour bar will run vertical instead of horizontal'
     print, '         _EXTRA = extra keywords. Colourbarfits accepts the normal '
     print, '                  plot keywords     '
     print, '         /BLACK = Set this keyword for black outline if not '
     print, '                  available in the colortable'
     print, '          COLOR = color of the outline and legends'
     goto,ending
  endif


if n_elements(legcolor) LE 0 then legcolor=0.
POS=[X[0],Y[0],X[1],Y[1]]

;drawing the colors

if not keyword_set(dummy1) then begin
    FOR z=0,255 DO BEGIN
        y1 = pos[0] + z/256.0*(pos[2]-pos[0])
        y2 = y1 + 1.0/256.0*(pos[2]-pos[0])
        POLYFILL,  [y1, y2, y2, y1,y1],[pos[1], pos[1], pos[3], pos[3], pos[1]], COLOR=z+0.5, /NORMAL
    ENDFOR
endif else begin
 FOR z=0,255 DO BEGIN
        y1 = pos[1] + z/256.0*(pos[3]-pos[1])
        y2 = y1 + 1.0/256.0*(pos[3]-pos[1])
        POLYFILL, [pos[0], pos[0], pos[2], pos[2], pos[0]],  [y1, y2, y2, y1,y1], COLOR=z+0.5, /NORMAL
    ENDFOR
endelse
if keyword_set(blackw) then begin
    TVLCT, Rin, Gin, Bin, /GET
    loadct,40,/SILENT
    legcolor=0.
endif

;plotting a title if requested
if n_elements(TITLE) GT 0. then begin
   if not keyword_set(multiple_ticks) then begin
      IF not keyword_set(opposite_label) then begin
         if not keyword_set(dummy1) then begin
            XYOUTS,[(X[0]+x[1])/2.],[Y[0]-0.03],title,alignment=0.5, $
                   /NORMAL, COLOR=legcolor,_EXTRA=ex
         endif else begin
            XYOUTS,[X[1]-0.02],[(Y[0]+Y[1])/2.],title,alignment=0.5, $
                   /NORMAL, COLOR=legcolor,_EXTRA=ex,ORIENTATION=90
         endelse
      ENDIF ELSE BEGIN
          if not keyword_set(dummy1) then begin
            XYOUTS,[(X[0]+x[1])/2.],[Y[1]+0.006],title,alignment=0.5, $
                   /NORMAL, COLOR=legcolor,_EXTRA=ex
         endif else begin
            XYOUTS,[X[1]+0.02],[(Y[0]+Y[1])/2.],title,alignment=0.5, $
                   /NORMAL, COLOR=legcolor,_EXTRA=ex,ORIENTATION=90
         endelse
      ENDELSE
   endif else begin
      IF not keyword_set(opposite_label) then begin
         if not keyword_set(dummy1) then begin
            XYOUTS,[(X[0]+x[1])/2.],[Y[0]-0.03],title,alignment=0.5, $
                   /NORMAL, COLOR=legcolor,_EXTRA=ex
         endif else begin
            XYOUTS,[X[0]-0.06],[(Y[0]+Y[1])/2.],title,alignment=0.5, $
                   /NORMAL, COLOR=legcolor,_EXTRA=ex,ORIENTATION=90
         endelse
      endif else begin
         if not keyword_set(dummy1) then begin
            XYOUTS,[(X[0]+x[1])/2.],[Y[1]+0.025],title,alignment=0.5, $
                   /NORMAL, COLOR=legcolor,_EXTRA=ex
         endif else begin
            XYOUTS,[X[1]+0.06],[(Y[0]+Y[1])/2.],title,alignment=0.5, $
                   /NORMAL, COLOR=legcolor,_EXTRA=ex,ORIENTATION=90
         endelse
      endelse
   endelse
endif
;Writing the minimum and maximum value
 a=[double(minval),double(maxval)]
 b=[0,1]
if not keyword_set(multiple_ticks) then begin
   plot,a,b,position=pos,color=legcolor, _EXTRA=ex,/NODATA,/NOERASE,ytickname=[' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' '],yticks=1.,xtickname=[' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' '],xticks=1.
   IF not keyword_set(supress) then begin
      IF not keyword_set(opposite_label) then begin
         if not keyword_set(dummy1) then begin
            XYOUTS,[X[0]],[Y[0]-0.03],/NORMAL, strtrim(string(minval),2), COLOR=legcolor,alignment=0.5, _EXTRA=ex
            XYOUTS,[X[1]],[Y[0]-0.03],/NORMAL, strtrim(string(maxval),2), COLOR=legcolor,alignment=0.5, _EXTRA=ex     
         endif else begin
            XYOUTS,[X[0]-0.02],[Y[0]], minval,alignment=0.5, $
                   /NORMAL, ORIENTATION=90, COLOR=legcolor,_EXTRA=ex
            XYOUTS,[X[0]-0.02],[Y[1]],/NORMAL, maxval,alignment=0.5, $
                   ORIENTATION=90, COLOR=legcolor,_EXTRA=ex
         endelse
      ENDIF ELSE BEGIN
         if not keyword_set(dummy1) then begin
            XYOUTS,[X[0]],[Y[1]+0.006],/NORMAL, strtrim(string(minval),2), COLOR=legcolor,alignment=0.5, _EXTRA=ex
            XYOUTS,[X[1]],[Y[1]+0.006],/NORMAL, strtrim(string(maxval),2), COLOR=legcolor,alignment=0.5, _EXTRA=ex     
         endif else begin
            XYOUTS,[X[1]+0.02],[Y[0]], minval,alignment=0.5, $
                   /NORMAL, ORIENTATION=90, COLOR=legcolor,_EXTRA=ex
            XYOUTS,[X[1]+0.02],[Y[1]],/NORMAL, maxval,alignment=0.5, $
                   ORIENTATION=90, COLOR=legcolor,_EXTRA=ex
         endelse
      ENDELSE
   ENDIF

endif else begin
  
   if  not keyword_set(dummy1) then begin
      if  not keyword_set(opposite_label) AND not keyword_set(supress) then begin
         plot,a,b,position=pos,color=legcolor, _EXTRA=ex,/NODATA,/NOERASE,ytickname=[' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' '],yticks=1.
      endif else begin
         plot,a,b,position=pos,color=legcolor, _EXTRA=ex,/NODATA,/NOERASE,ytickname=[' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' '],yticks=1.,xtickname=[' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' '],xticks=1
         if not keyword_set(supress) then begin
            if not n_elements(xtf) then begin
               IF not n_elements(xtn) then begin
                  Axis,XAXIS=1,_EXTRA=ex
               endif else begin
                  AXIS,XAXIS=1,_EXTRA=ex,xtickv=double(xtn),xticks=n_elements(xtn)-1
               endelse
            endif else begin
               IF not n_elements(xtn) then begin
                  Axis,XAXIS=1,_EXTRA=ex,xtickformat=xtf
               endif else begin
                  AXIS,XAXIS=1,_EXTRA=ex,xtickv=double(xtn),xticks=n_elements(xtn)-1,xtickformat=xtf
               endelse
            endelse
         Endif
      endelse

   endif else begin
      if  not keyword_set(opposite_label) and not keyword_set(supress) then begin
         plot,b,a,position=pos,color=legcolor, _EXTRA=ex,/NODATA,/NOERASE,xtickname=[' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' '],xticks=1.
         
      endif else begin
         plot,b,a,position=pos,color=legcolor, _EXTRA=ex,/NODATA,/NOERASE,yrange=[double(minval),double(maxval)],ystyle=1,xticks=1,yticks=1,ytickname=[' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' '],xtickname=[' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ']
         if not keyword_set(supress) then begin
            if not n_elements(ytf) then begin
               IF not n_elements(ytn) then begin
                  AXIS,YAXIS=1,_EXTRA=ex,COLOR=legcolor,ystyle=1
               endif else begin
                  AXIS,YAXIS=1,_EXTRA=ex,ytickv=double(ytn),yticks=n_elements(ytn)-1,COLOR=legcolor
               endelse
            endif else begin
               IF not n_elements(ytn) then begin
                  AXIS,YAXIS=1,_EXTRA=ex,ytickformat=ytf,COLOR=legcolor,ystyle=1
               endif else begin
                  AXIS,YAXIS=1,_EXTRA=ex,ytickname=ytn,ytickv=double(ytn),yticks=n_elements(ytn)-1,COLOR=legcolor,ytickformat=ytf
               endelse
            endelse
         endif
      endelse
      
   endelse
endelse


if keyword_set(blackw) then begin
   TVLCT, Rin, Gin, Bin
endif
ending: 
END
