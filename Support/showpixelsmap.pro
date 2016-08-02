pro showpixelsmap, xaxis, yaxis, int, BLANK_VALUE=blankval,RANGE=range, PLOT_ALL=plotall,PIXELSIZE=pixelsize,WCS=wcs, BLACK=blackw,_EXTRA=ex, $
                   XRANGE=xrng,YRANGE=yrng,XTICKS=n_ticksx,YTICKS=n_ticksy, FORCE= force,XTICKNAME=xtn,YTICKNAME=ytn,HEX_COLOR=hex_color
  compile_opt idl2
;+
; NAME:
;       SHOWPIXELSMAP
;
; PURPOSE:
;       This routine will plot a fits file into a colormap on the
;       current DEVICE. It is up to the user to ensure the right ratio
;       of the plot vs fits file if square pixels are to plotted.
;
; CATEGORY:
;       Support
; 
; CALLING SEQUENCE:
;       SHOWPIXELSMAP, xaxis, yaxis, int, BLANK_VALUE=blankval,RANGE=range, PLOT_ALL=plotall,PIXELSIZE=pixelsize,WCS=wcs, BLACK=blackw,_EXTRA=ex,
;                      XRANGE=xrng,YRANGE=yrng,XTICKS=n_ticksx,YTICKS=n_ticksy, FORCE= force,XTICKNAME=xtn,YTICKNAME=ytn
;
;
; INPUTS:
;              int = values of the pixels
;            xaxis = the xaxis of the map to be displayed (location of
;                    pixels centre)
;            yaxis = the yaxis of the map to be displayed (location of
;                    pixels centre) 
; 
; OPTIONAL INPUTS:
;       BLANK_VALUE = value to treated as a blank 
;            _EXTRA = Additional plotting parameters
;         PIXELSIZE = The size of a pixel in units of xaxis/yaxis
;             RANGE = Intensity range over which the colormap should
;                     vary. Values outside this range will not be
;                     plotted unless /PLOT_ALL is set.
;            XRANGE = Range of xaxis (The size of the map is
;                     irrelevant)
;            XTICKS = number of xtick marks
;         XTICKNAME = name of the xtick marks  
;            YRANGE = Range of yaxis (The size of the map is
;                     irrelevant)
;            YTICKS = number of ytick marks
;         YTICKNAME = name of the ytick marks     
;           
; KEYWORD PARAMETERS:
;            /BLACK - Outlines should be in black regardless of the colortable
;            /FORCE - forces the exact tickmarks in WCS mode
;         /PLOT_ALL - Set this keyword to ensure that all values are
;                     plotted.
;              /WCS - Axis labels are WCS coordinates
;  
; OUTPUTS:
;       The map (Likely fits image) defined by xaxis,yaxis,int is
;       plotted onto the current DEVICE
;
; OPTIONAL OUTPUTS:
;       -
;
; PROCEDURES CALLED:
;       DEC_NAMES, KEYWORD_SET(), N_ELEMENTS(), MAX(), POLYFILL, PLOT, RA_NAMES
;
; EXAMPLE:
;      
;
; MODIFICATION HISTORY:
;       08-02-2011 P.Kamphuis; Can now deal with input map larger than range and larger
;                  range than input map in wcs coord. Non finite and blankvalues are not plotted.
;       22-09-2009 P.Kamphuis; Added an error message, enabled
;                  axis' that are larger than the cube.  
;       03-09-2009 Written by P.Kamphuis; 
;
; NOTE:
;  

  CATCH,Error_status
  IF  Error_status NE 0. THEN BEGIN
     print, 'Oops the following went wrong:'
     print, !ERROR_STATE.MSG
     print, 'Use Showpixelsmap in this way:'
     print,'CALLING SEQUENCE:showpixelsmap,xaxis,yaxis, int, EXTRA=extra'
     
     print, '      xaxis = the x central positions of the cells'
     print, '      yaxis = the y central positions of the cells'
     print, '        int = the values of the image. Has to be a 2D array'
     print, 'BLANK_VALUE = value of blanks  [-987654]'
     print, '      RANGE = color range values outside are blanked automatically'
     print, '  /PLOT_ALL = if set then all values are plotted'
     print, '  PIXELSIZE = Size of the pixels plotted default is calculated from the axii '
     print, '       /WCS = if keyword is set it is assumed that x axis is right'
     print, '              ascension in degrees and y axis is DEC in degrees'
     print, '     /BLACK = if keyword is set axis and labels will be in black' 
     print, '      /FORCE = forces the exact tickmarks in WCS mode'
     goto,ending
  endif
  if n_elements(int) EQ 0. then begin
     print, 'SHOWPIXELSMAP: Your input cube has no values'  
     print, 'Use Showpixelsmap in this way:'
     print,'CALLING SEQUENCE:showpixelsmap,xaxis,yaxis,int, EXTRA=extra'
     goto, ending
  endif
  if n_elements(blankval) LT 1 then blankval=-987654
  if n_elements(xrng) LT 1 then xrng=[xaxis[0],xaxis[n_elements(xaxis)-1]]
  if n_elements(yrng) LT 1 then yrng=[yaxis[0],yaxis[n_elements(yaxis)-1]]

  nxax = n_elements(xaxis)
  nyax = n_elements(yaxis)
  nxCube=n_elements(int[*,0])
  nyCube=n_elements(int[0,*])
  IF nxCube GT nxax then nx=nxax else nx=nxCube
  IF nyCube GT nyax then ny=nyax else ny=nyCube


  if n_elements(pixelsize) NE 0. then begin
     pixelsizex =abs(pixelsize)
     pixelsizey =abs(pixelsize)
  endif  else begin
     pixelsizex = ABS(xaxis[fix(nx/2)]-xaxis[fix(nx/2)+1]) 
     pixelsizey = ABS(yaxis[fix(ny/2)]-yaxis[fix(ny/2)+1])
  endelse
  maxx = max(xrng, MIN=minx)
  maxy = max(yrng, MIN=miny)
  if keyword_set(blackw) then begin
     TVLCT, Rin, Gin, Bin, /GET
     loadct,0,/SILENT 
  endif
  
  IF xrng[0] LT xrng[1] then begin 
     xrng=[minx-(0.5*pixelsizex), maxx+(0.5*pixelsizex)]
  ENDIF ELSE  xrng=[ maxx+(0.5*pixelsizex),minx-(0.5*pixelsizex)]
  IF yrng[0] LT yrng[1] then begin 
     yrng=[miny-(0.5*pixelsizey), maxy+(0.5*pixelsizey)]
  ENDIF ELSE  yrng=[ maxy+(0.5*pixelsizey),miny-(0.5*pixelsizey)]

  If not keyword_set(wcs) then begin 
     PLOT, xrng, yrng, /NODATA,/XSTYLE,/YSTYLE, _EXTRA=ex,COLOR=0,xrange=xrng,yrange=yrng,$
           xticks =n_ticksx, yticks =n_ticksy,xtickname=xtn,ytickname=ytn
  endif else begin
     if n_elements(xtn) GT 0 then begin
        xtntmp=xtn
     endif

     if n_elements(ytn) GT 0 then ytntmp=ytn
     if keyword_set(force) then begin
        ra_names, xrng, tick_name = xtn, tick_value = xtv, $
                  n_ticks = n_ticksx, incr = incrx,/FORCE
        dec_names,yrng, tick_name = ytn, tick_value = ytv, $
                  n_ticks = n_ticksy, incr = incry,/FORCE

     ENDIF else begin
        
        ra_names, xrng, tick_name = xtn, tick_value = xtv, $
                  n_ticks = n_ticksx, incr = incrx
        dec_names,yrng, tick_name = ytn, tick_value = ytv, $
                  n_ticks = n_ticksy, incr = incry
     ENDELSE
     IF n_elements(xtntmp) GT 0 then xtn=xtntmp
     IF n_elements(ytntmp) GT 0 then ytn=ytntmp
                                ; print,xtn,ytn
     PLOT, xrng,yrng,/NODATA,/XSTYLE,/YSTYLE, _EXTRA=ex, xtickname = xtn,        $
           ytickname = ytn, xticks = n_elements(xtv)-2,           $
           yticks = n_elements(ytv)-2,xminor = 4, yminor = 4,     $
           ytickv = ytv, xtickv = xtv,COLOR=0,xrange=xrng,yrange=yrng
  endelse

  if keyword_set(blackw) then begin
     TVLCT, Rin, Gin, Bin
  endif

  if n_elements(range) GT 0. then begin
     intmin = range[0]
  endif else begin
     tmp1=where(int NE blankval AND FINITE(int) NE 0.)
     intmax = MAX(int[tmp1],MIN=intmin)
     range=dblarr(2)
     range[0] = intmin
     range[1] = intmax
  endelse
  colorRange = range[1] - range[0]
  xblock=pixelsizex*[-0.5, -0.5, +0.5, +0.5, -0.5]
  yblock=pixelsizey*[+0.5, -0.5, -0.5, +0.5, +0.5]
  
  for i=0, nx-1 do begin
     FOR j=0, ny-1 DO BEGIN        
        IF xaxis[i] GE minx AND xaxis[i] LE maxx $
           AND yaxis[j] GE miny AND yaxis[j] LE maxy then begin
           
           xpix = xaxis[i]+xblock
           ypix = yaxis[j]+yblock
           if keyword_set(plotall) then begin
              IF int[i,j] NE blankval AND FINITE(int[i,j]) then begin
                 POLYFILL, xpix, ypix, COLOR=(int[i,j]-intmin)/colorrange*255<255>0
              endif
           endif else begin
              IF int[i,j] NE blankval AND int[i,j] GE Range[0] $
                 AND int[i,j] LE Range[1] AND FINITE(int[i,j]) then begin
                 POLYFILL, xpix, ypix, COLOR=(int[i,j]-intmin)/colorrange*255<255>0
              endif
           endelse
        endif
     ENDFOR
  endfor


  if keyword_set(blackw) then begin
     TVLCT, Rin, Gin, Bin, /GET
     loadct,0,/SILENT
  endif
  If not keyword_set(wcs) then begin 
     PLOT, xrng, yrng, /NODATA,/XSTYLE,/YSTYLE,/NOERASE,COLOR=1,_EXTRA=ex,$
           xrange=xrng,yrange=yrng, xticks =n_ticksx, yticks =n_ticksy,xtickname=xtn,ytickname=ytn
  endif else begin
     PLOT,xrng, yrng,/NODATA,/NOERASE,/XSTYLE,/YSTYLE, _EXTRA=ex, xtickname = xtn,$
           ytickname = ytn, xticks = n_elements(xtv)-2,           $
           yticks = n_elements(ytv)-2,xminor = 4, yminor = 4,     $
           ytickv = ytv, xtickv = xtv,COLOR=0,xrange=xrng,yrange=yrng
  endelse
  if keyword_set(blackw) then begin
     TVLCT, Rin, Gin, Bin
  endif
ending:
END
