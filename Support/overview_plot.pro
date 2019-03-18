Pro overview_plot,distance,gdlidl,noise=noise,finishafter = finishafter,filenames=filenames,splined=splined,version=version,fixedpars=fixedpars

;+
; NAME:
;       overview_plot
;
; PURPOSE:
;       Create an overview plot to easily see how well the model fits the data.
;
; CATEGORY:
;       Support
; 
; CALLING SEQUENCE:
;       OVERVIEW_PLOT	
;
;
; INPUTS:
;         distance = Distance to the galaxy
;         gdlidl = an indicator whether we are running idl or gdl
; OPTIONAL INPUTS:
;         Noise = the noise in the cube.
;   finishafter = key for what kind of fitting was done.
;   filenames   = nams of the created files   
; KEYWORD PARAMETERS:
;       /SPLINED - Plot with spline interpolation between the points.
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:00B4FF
;       -
; 
; PROCEDURES CALLED:
;       AXIS,COLORMAPS,COLOUR_BAR,COLUMNDENSITY,CONTOUR,DEVICE,FILE_TEST(),LOADCT,MAX(),N_ELEMENTS(),OPLOT,PLOT,FAT_PLOTERR,READFITS(),SET_PLOT,SHOWPIXELSMAP,WRITENEWTOTEMPLATE,XYOUTS
;
; MODIFICATION HISTORY:
;       10-10-2018 P.Kamphuis; An existing Overview.png will be moved
;                              to Overview_Prev.png.  
;       10-10-2018 P.Kamphuis; Addition of SDIS plotting as radial parameter.  
;       14-11-2017 P.Kamphuis; Increased and improved checking for convert      
;       22-03-2017 P.Kamphuis; Increased checking for Imagick convert
;                              in order to trim the png. Improved
;                              plotting limits on velocity contours
;                              and PV-Diagram color scale.  
;       15-11-2016 P.Kamphuis; Fixed a bounding Box issue in GDL  
;       15-11-2016 P.Kamphuis; Introduced hexadecimal plotting for
;                              color in idl and gdl. In idl and gdl
;                              this is coded as RGB colors converting
;                              to hex(B)+hex(G)+hex(R).   
;       08-11-2016 P.Kamphuis; Added the splined keyword to have the
;                              variables plotted in a spline
;                              interpolation
;       08-11-2016 P.Kamphuis; Removed the plotting through the X
;                              widget for GDL as it cannot be run in
;                              screen in that case. However this
;                              results in black and white plots. Need
;                              to move to hexidecimal plotting   
;       01-08-2016 P.Kamphuis; Updated the GDL plotting to use png as
;                              well. This now happens through the X
;                              buffer which unfortunately means that
;                              it is limited to the amount of pixels
;                              in the users screen. This is because
;                              the Z-buffer in GDL does not accept the
;                              decomposed keyword.
;       09-06-2016 P.Kamphuis; Case matched the call to Finalmodel as
;                              unix is case sensitive.  
;       23-03-2016 P.Kamphuis; Added filename and finishafter
;                              extension to make sure that the right files are picked up and
;                              to provide info about the fit  
;       Written 20-02-2016 P.Kamphuis v1.0
;
; NOTE:
;     
;-

  
  arrays=1.
  IF gdlidl then SET_PLOT,'PS' else SET_PLOT, 'Z'
  plotpara=['RADI','SBR','SBR_2','VROT','VROT_ERR','PA','PA_ERR','PA_2','PA_2_ERR','INCL','INCL_ERR','INCL_2','INCL_2_ERR','BMAJ','SDIS','XPOS','YPOS','VSYS','SDIS_ERR']
  plotstart=[[1,3,14,5,9],[2,3,14,7,11],[0,1,4,1,1]]
  Template=1.
  WriteNewToTemplate,Template,'Finalmodel/Finalmodel.def',ARRAYS=Arrays,VARIABLECHANGE=plotpara,/EXTRACT
  IF FILE_TEST('ModelInput.def') then begin
     Writenewtotemplate,Template,'ModelInput.def',ARRAYS=ModArrays,VARIABLECHANGE=plotpara,/EXTRACT
  ENDIF
  varunits=strarr(n_elements(plotpara))
  for i=0,n_elements(plotpara)-1 do begin
     tmp=str_sep(strtrim(strcompress(plotpara[i]),2),'_')
     case tmp[0] of
        'VROT':Varunits[i]='(km s!E-1!N)'
        'SDIS':Varunits[i]='(km s!E-1!N)'
        'PA' :Varunits[i]='(Degrees)'
        'INCL':Varunits[i]='(Degrees)'
        'SBR':Varunits[i]='(Jy km s!E-1!N arcsec!E-2!N)'
        'Z0':Varunits[i]='(Arcsec)'
        else:Varunits[i]=''
     endcase
  endfor
  maxvar=dblarr(5)
  minvar=dblarr(5)
  buffer=dblarr(5)
  minvar[*]=100000
  maxvar[*]=-10000
  for i=0,4 do begin
     tmpvals=[Arrays[*,plotstart[i,0]],Arrays[*,plotstart[i,1]]]
     tmplocs=WHERE(tmpvals NE 0.)
     tmpmax=MAX(tmpvals[tmplocs],MIN=tmpmin)
     IF tmpmax GT Maxvar[i] then Maxvar[i]=tmpmax
     IF tmpmin LT Minvar[i] then Minvar[i]=tmpmin
     buffer[i]=(ABS(Maxvar[i])+ABS(minvar[i]))/20.
  endfor
  tmp=WHERE(plotpara EQ 'XPOS')
  RA=double(Arrays[0,tmp])
  tmp=WHERE(plotpara EQ 'YPOS')
  DEC=double(Arrays[0,tmp])
  convertradec,RA,DEC
  tmp=WHERE(plotpara EQ 'VSYS')
  vsys=strtrim(string(double(Arrays[0,tmp]),format='(F10.1)'),2)
  tmp=WHERE(plotpara EQ 'SDIS')
  disper=strtrim(string(double(Arrays[0,tmp]),format='(F10.1)'),2)
  tmp=WHERE(plotpara EQ 'BMAJ')
  majbeam=strtrim(string(double(Arrays[0,tmp]),format='(F10.1)'),2)
  tmp=WHERE(plotpara EQ 'RADI')
  out_ringsize=strtrim(string(double(Arrays[n_elements(Arrays[*,0])-1,tmp]-Arrays[n_elements(Arrays[*,0])-2,tmp]),format='(F10.1)'),2)
  in_ringsize=strtrim(string(double(Arrays[2,tmp]-Arrays[1,tmp]),format='(F10.1)'),2)
  tmp=WHERE(plotpara EQ 'INCL')
  ceninc=Arrays[0,tmp]
  tmp=WHERE(plotpara EQ 'VROT')
  maxvrot=MAX(Arrays[*,tmp])
  ysize=0.1
  !x.style=1.5
  !y.style=1.5
  !p.charsize=3.7
  !p.charthick=2.15569
  !p.thick=6
  !P.FONT=1
  !p.background=0
  thick=2
  scrdim = [8.27*300.,11.69*300]
  A = FIndGen(16) * (!PI*2/16.) 
  UserSym, cos(A), sin(A), /fill
  spawn,'mv Overview.png Overview_Prev.png'
  IF gdlidl then begin

;                                ;Currently GDL does not recognize true
                                ;type fonts yet. This leads to errors
                                ;in using the degree symbol. It also
                                ;does not yet recognize superscript
                                ;commands in tickmarks. Additonally as
                                ;we cannot use the Z-buffer as it does
                                ;not accept the decomposed keyword and
                                ;because we do not want to be screen
                                ;size depndent we need a widget.
     !p.charsize=0.4
     !p.symsize=0.3
     charthick=0.7
     ssize=0.3
     DEVICE,xsize=scrdim[0]/200.,ysize=scrdim[1]/200,FILENAME='Overview.ps',/color,/PORTRAIT,/DECOMPOSED,/ENCAPSULATED
  endif else begin
     !p.charsize=3.7
     charthick=2.15569
     !p.symsize=2.5
     ssize=2.5
     DEVICE,SET_FONT='Times',/TT_FONT,SET_RESOLUTION=[scrdim[0],scrdim[1]],/DECOMPOSED,SET_PIXEL_DEPTH=24
  endelse
  tmp=WHERE(plotpara EQ 'RADI')
  plotradii=Arrays[*,tmp]
  tmppos=WHERE(plotpara EQ 'SBR')
  tmp=WHERE(Arrays[*,tmppos[0]] GT 1.1E-16)
  tmppos=WHERE(plotpara EQ 'SBR_2')
  tmp2=WHERE(Arrays[*,tmppos[0]] GT 1.1E-16)
  
  maxradii=MAX([plotradii[tmp],plotradii[tmp2]])+(plotradii[n_elements(plotradii)-1]-plotradii[n_elements(plotradii)-2])/2.
  
  for i=0,4 do begin    
     IF i EQ 0 then begin
        tmppos=WHERE(plotpara EQ 'SBR')
        plotvariable=Arrays[tmp,tmppos[0]]
        loadct,0,/silent
        xerr=dblarr(n_elements(plotVariable))
       
        ;plot,plotradii,plotVariable,position=[0.15,0.95-5.*ysize,0.55,0.95-4.*ysize],xtitle='Radius (arcmin)',$
         ;    xrange=[0.,maxradii],yrange=[minvar[i]-buffer[i],maxvar[i]+buffer[i]],ytickname=[' ',' ',' ',' ',' ',' ',' ',' ',' ',' '],xtickname=[' ',' ',' ',' ',' ',' ',' ',' ',' ',' '],xticklayout=1,background='ffffff'x,color='ffffff'x,/nodata
        fat_ploterror,plotradii,plotVariable,xerr,xerr,position=[0.15,0.95-(5.-i)*ysize,0.55,0.95-(4.-i)*ysize],$
                         xrange=[0.,maxradii],yrange=[minvar[i]-buffer[i],maxvar[i]+buffer[i]],xthick=xthick,ythick=ythick,xtickname=[' ',' ',' ',' ',' ',' ',' ',' ',' ',' '],ytickname=[' ',' ',' ',' ',' ',' ',' ',' ',' ',' '],xticklayout=1,yticklayout=1,charthick=charthick,thick=thick,charsize=charsize,linestyle=0,$
                         color='000000'x,ERRCOLOR = '000000'x, ERRTHICK=!p.thick*0.4,psym=8,symsize=ssize
        if keyword_set(splined) then begin
           newrad=dblarr((n_elements(plotradii)-1)*10.+1)
           for h=1,n_elements(plotradii)-1 do begin
              newrad[(h-1)*10:h*10-1]=findgen(10)*(plotradii[h]-plotradii[h-1])/10.+plotradii[h-1]
           endfor
           newrad[(h-1)*10]=plotradii[h-1]
           newvar=spline(plotradii,plotVariable,newrad)
           oplot,newrad,newvar,color='000000'x,linestyle=0
        ENDIF ELSE oplot,plotradii,plotVariable,color='000000'x,linestyle=0     
        levelsrange=[minvar[i]-buffer[i],maxvar[i]+buffer[i]]*1000.
        ;fat_ploterror,plotradii,plotVariable,xerr,xerr,thick=lthick,color='000000'x,linestyle=2,ERRCOLOR = '000000'x, ERRTHICK=!p.thick*0.4,/over_plot,psym=8,symsize=ssize
        
        columndensity,levelsrange,double(vsys),[1.,1.],vwidth=1.,/arcsquare
        levelsranges=levelsrange/1e20
        reset=1e20
        levelsranges[0]=ceil(levelsranges[0])
        levelsranges[1]=floor(levelsranges[1])
        adst='x10!E20!N'
        if levelsranges[0] EQ levelsranges[1] then begin
           levelsranges=levelsrange/1e19
           levelsranges[0]=ceil(levelsranges[0])
           levelsranges[1]=floor(levelsranges[1])
           adst='x10!E19!N'
           reset=1e19
        endif
        midrange=levelsranges[0]+(levelsranges[1]-levelsranges[0])/2.
        IF fix(midrange) NE midrange then begin
           levelsranges[1]=levelsranges[1]-1 
           midrange=levelsranges[0]+(levelsranges[1]-levelsranges[0])/2.
        endif
        newlevels=[levelsranges[0],midrange,levelsranges[1]]
        jynewlevels=newlevels*reset
        columndensity,jynewlevels,double(vsys),[1.,1.],vwidth=1.,/NCOLUMN,/arcsquare
        AXIS,YAXIS=0,charthick=charthick,xthick=xthick,ythick=ythick,charsize=charsize,color='000000'x
        AXIS,YAXIS=1,charthick=charthick,xthick=xthick,ythick=ythick,charsize=charsize,ytickv=jynewlevels/1000.,ytickname=[strtrim(strcompress(string(newlevels[0],format='(I3)')),2),strtrim(strcompress(string(fix(newlevels[1]),format='(I2)')),2),strtrim(strcompress(string(fix(newlevels[2]),format='(I2)')),2)],yticks=2,yminor=3,color='000000'x
        XYOUTs,0.60,0.95-4.5*ysize,'N!IH',/NORMAL,alignment=0.5,ORIENTATION=90, CHARTHICK=charthick,charsize=!p.charsize*1.25,color='000000'x
        XYOUTs,0.63,0.95-4.5*ysize,'('+adst+' cm!E-2!N)' ,/NORMAL,alignment=0.5,ORIENTATION=90, CHARTHICK=charthick,charsize=charsize,color='000000'x
        AXIS,XAXIS=0,charthick=charthick,xthick=xthick,ythick=ythick,charsize=charsize,color='000000'x ,XTITLE='Radius (arcsec)'
        AXIS,XAXIS=1,charthick=charthick,xthick=xthick,ythick=ythick,charsize=charsize,XRANGE = convertskyanglefunction(!X.CRANGE,distance),xtickname=[' ',' ',' ',' ',' ',' ',' ',' ',' ',' '],color='000000'x 
        loadct,40,/silent
        tmppos=WHERE(plotpara EQ 'SBR_2')
        xerr=dblarr(n_elements(Arrays[tmp2,tmppos[0]]))
        fat_ploterror,plotradii,Arrays[tmp2,tmppos[0]],xerr,xerr,thick=lthick,color='0000FF'x,linestyle=2,ERRCOLOR = '0000FF'x, ERRTHICK=!p.thick*0.4,/over_plot,psym=8,symsize=ssize
        if keyword_set(splined) then begin
           newvar=spline(plotradii,Arrays[tmp2,tmppos[0]],newrad)
           oplot,newrad,newvar,color='0000FF'x,linestyle=2
        ENDIF ELSE oplot,plotradii,Arrays[tmp2,tmppos[0]],thick=lthick,color='0000FF'x,linestyle=2
      
        XYOUTs,0.05,0.95-4.5*ysize,plotpara[plotstart[i,0]],/NORMAL,alignment=0.5,ORIENTATION=90,charsize=!p.charsize*1.25,color='000000'x,charthick=charthick
        XYOUTs,0.08,0.95-4.5*ysize,varunits[plotstart[i,0]],/NORMAL,alignment=0.5,ORIENTATION=90,color='000000'x,charthick=charthick
        
        IF FILE_TEST('ModelInput.def') then begin
           oplot,ModArrays[*,0],ModArrays[*,1],thick=lthick,color='FF0010'x
           oplot,ModArrays[*,0],ModArrays[*,1],psym=8,color='FF0010'x,symsize=ssize
           oplot,ModArrays[*,0],ModArrays[*,2],thick=lthick,color='00B4FF'x,linestyle=2
           oplot,ModArrays[*,0],ModArrays[*,2],psym=8,color='00B4FF'x,linestyle=2,symsize=ssize
        ENDIF
     ENDIF ELSE begin
        IF plotstart[i,0] NE plotstart[i,1] then begin
           plotvariable=Arrays[tmp,plotstart[i,0]]
           plotVariableErr=Arrays[tmp,plotstart[i,0]+plotstart[i,2]]
        ENDIF ELSE BEGIN
           IF tmp[n_elements(tmp)-1] GT tmp2[n_elements(tmp2)-1] then begin
              plotvariable=Arrays[tmp,plotstart[i,0]]
              plotVariableErr=Arrays[tmp,plotstart[i,0]+plotstart[i,2]]
           ENDIF ELSE BEGIN
              plotvariable=Arrays[tmp2,plotstart[i,0]]
              plotVariableErr=Arrays[tmp2,plotstart[i,0]+plotstart[i,2]]
           ENDELSE
        ENDELSE
        loadct,0,/silent
        xerr=dblarr(n_elements(plotVariableErr))
        IF TOTAL(plotVariableErr) NE 0. then begin
           fat_ploterror,plotradii,plotVariable,xerr,plotVariableErr,position=[0.15,0.95-(5.-i)*ysize,0.55,0.95-(4.-i)*ysize],$
                     xrange=[0.,maxradii],yrange=[minvar[i]-buffer[i],maxvar[i]+buffer[i]],xthick=xthick,ythick=ythick,xtickname=[' ',' ',' ',' ',' ',' ',' ',' ',' ',' '],xticklayout=1,charthick=charthick,thick=thick,charsize=charsize,linestyle=0,$
                     /noerase,color='000000'x,ERRCOLOR = '000000'x, ERRTHICK=!p.thick*0.4,psym=8,symsize=ssize
        ENDIF ELSE BEGIN
           fat_ploterror,plotradii,plotVariable,xerr,xerr,position=[0.15,0.95-(5.-i)*ysize,0.55,0.95-(4.-i)*ysize],$
                         xrange=[0.,maxradii],yrange=[minvar[i]-buffer[i],maxvar[i]+buffer[i]],xthick=xthick,ythick=ythick,xtickname=[' ',' ',' ',' ',' ',' ',' ',' ',' ',' '],xticklayout=1,charthick=charthick,thick=thick,charsize=charsize,linestyle=0,$
                         /noerase,color='000000'x,ERRCOLOR = '000000'x, ERRTHICK=!p.thick*0.4,psym=8,symsize=ssize

           
        ENDELSE
        if keyword_set(splined) then begin
           newvar=spline(plotradii,plotVariable,newrad)
           oplot,newrad,newvar,color='000000'x,linestyle=0,thick=lthick
        ENDIF ELSE oplot,plotradii,plotVariable,thick=lthick,color='000000'x,linestyle=0
      
        IF i EQ 4 then begin
           AXIS,XAXIS=1,charthick=charthick,xthick=xthick,ythick=ythick,charsize=charsize,XRANGE = convertskyanglefunction(!X.CRANGE,distance),XTITLE='Radius (kpc)',color='000000'x 
        endif else begin
           AXIS,XAXIS=1,charthick=charthick,xthick=xthick,ythick=ythick,charsize=charsize,XRANGE = convertskyanglefunction(!X.CRANGE,distance),xtickname=[' ',' ',' ',' ',' ',' ',' ',' ',' ',' '],color='000000'x 
        endelse
        XYOUTs,0.05,0.95-(4.5-i)*ysize,plotpara[plotstart[i,0]],/NORMAL,alignment=0.5,ORIENTATION=90,charsize=!p.charsize*1.25,color='000000'x,charthick=charthick
        XYOUTs,0.08,0.95-(4.5-i)*ysize,varunits[plotstart[i,0]],/NORMAL,alignment=0.5,ORIENTATION=90,color='000000'x,charthick=charthick
        IF plotpara[plotstart[i,0]] EQ 'PA' and fixedpars[0] EQ 0. then XYOUTs,0.53,0.93-(4.-i)*ysize,'Fixed',/NORMAL,alignment=1.0,charsize=!p.charsize*1.25,color='000000'x,charthick=charthick
        IF plotpara[plotstart[i,0]] EQ 'INCL' and fixedpars[1] EQ 0. then XYOUTs,0.53,0.93-(4.-i)*ysize,'Fixed',/NORMAL,alignment=1.0,charsize=!p.charsize*1.25,color='000000'x,charthick=charthick
        IF plotpara[plotstart[i,0]] EQ 'SDIS' and fixedpars[2] EQ 0. then XYOUTs,0.53,0.93-(4.-i)*ysize,'Fixed',/NORMAL,alignment=1.0,charsize=!p.charsize*1.25,color='000000'x,charthick=charthick
       loadct,40,/silent
        IF plotstart[i,0] NE plotstart[i,1] then begin
           plotvariable=Arrays[tmp2,plotstart[i,1]]
           plotVariableErr=Arrays[tmp2,plotstart[i,1]+plotstart[i,2]]          
           IF TOTAL(plotVariableErr) NE 0. then begin
             xerr=dblarr(n_elements(plotVariableErr))
              fat_ploterror,plotradii,plotVariable,xerr,plotVariableErr,thick=lthick,color='0000FF'x,linestyle=2,ERRCOLOR = '0000FF'x, ERRTHICK=!p.thick*0.4,/over_plot,psym=8,symsize=ssize
           ENDIF ELSE BEGIN
              xerr=dblarr(n_elements(plotVariableErr))
              fat_ploterror,plotradii,plotVariable,xerr,xerr,thick=lthick,color='0000FF'x,linestyle=2,ERRCOLOR = '0000FF'x, ERRTHICK=!p.thick*0.4,/over_plot,psym=8,symsize=ssize
           ENDELSE
           if keyword_set(splined) then begin
              newvar=spline(plotradii,plotVariable,newrad)
              oplot,newrad,newvar,color='0000FF'x,linestyle=2
           ENDIF ELSE  oplot,plotradii,plotVariable,thick=lthick,color='0000FF'x,linestyle=2
      
        ENDIF
        IF FILE_TEST('ModelInput.def') then begin
           oplot,ModArrays[*,0],ModArrays[*,plotstart[i,0]],thick=lthick,color='FF0010'x
           oplot,ModArrays[*,0],ModArrays[*,plotstart[i,0]],psym=8,color='FF0010'x,symsize=ssize
           IF plotstart[i,0] NE plotstart[i,1] then begin
              oplot,ModArrays[*,0],ModArrays[*,plotstart[i,1]],thick=lthick,color='00B4FF'x,linestyle=2
              oplot,ModArrays[*,0],ModArrays[*,plotstart[i,1]],psym=8,color='00B4FF'x,linestyle=2,symsize=ssize
           ENDIF
        ENDIF
     ENDELSE
  endfor
  IF FILE_TEST('ModelInput.def') then begin
     tmp=WHERE(plotpara EQ 'XPOS')
     RAmod=double(ModArrays[0,tmp])
     tmp=WHERE(plotpara EQ 'YPOS')
     DECmod=double(ModArrays[0,tmp])
     convertradec,RAmod,DECmod
     tmp=WHERE(plotpara EQ 'VSYS')
     vsysmod=strtrim(string(double(ModArrays[0,tmp]),format='(F10.1)'),2)
     dispermod=strtrim(string(double(ModArrays[0,n_elements(plotpara)-4]),format='(F10.1)'),2)
     XYOUTS,0.60,0.89,'Systemic Velocity= '+vsys+' ('+vsysmod+') km s!E-1',/normal,alignment=0.,charthick=charthick,color='000000'x
     XYOUTS,0.60,0.87,'R.A.= '+RA+' ('+RAmod+')',/normal,alignment=0.,charthick=charthick,color='000000'x
     XYOUTS,0.60,0.85,'DEC.= '+DEC+' ('+DECmod+')',/normal,alignment=0.,charthick=charthick,color='000000'x
     XYOUTS,0.60,0.83,'Black lines: approaching side parameters.',/normal,alignment=0.,charthick=charthick,color='000000'x
     XYOUTS,0.60,0.81,'Red lines: receding side parameters.',/normal,alignment=0.,charthick=charthick,color='000000'x
     XYOUTS,0.60,0.79,'Blue lines: approaching side input model parameters.',/normal,alignment=0.,charthick=charthick,color='000000'x
     XYOUTS,0.60,0.77,'Yellow lines: receding side input model parameters.',/normal,alignment=0.,charthick=charthick,color='000000'x
   
     XYOUTS,0.60,0.75,'The major FWHM beam is '+majbeam+'" ',/normal,color='000000'x
     IF in_ringsize EQ out_ringsize then XYOUTS,0.60,0.73,'We used rings of size '+in_ringsize+'"',/normal,color='000000'x else begin
        XYOUTS,0.60,0.73,'We used rings of size '+in_ringsize+'" in the inner part',/normal,color='000000'x
        XYOUTS,0.60,0.71,'We used rings of size '+out_ringsize+'" in the outer part',/normal,color='000000'x
     ENDELSE
  ENDIF ELSE BEGIN
     XYOUTS,0.60,0.89,'Systemic Velocity= '+vsys+' km s!E-1',/normal,alignment=0.,charthick=charthick,color='000000'x
     XYOUTS,0.60,0.87,'R.A.= '+RA,/normal,alignment=0.,charthick=charthick,color='000000'x
     XYOUTS,0.60,0.85,'DEC.= '+DEC,/normal,alignment=0.,charthick=charthick,color='000000'x
     XYOUTS,0.60,0.83,'Black lines: approaching side parameters.',/normal,alignment=0.,charthick=charthick,color='000000'x
     XYOUTS,0.60,0.81,'Red lines: receding side parameters.',/normal,alignment=0.,charthick=charthick,color='000000'x
     IF in_ringsize EQ out_ringsize then XYOUTS,0.60,0.79,'We used rings of size '+in_ringsize+'"',/normal,color='000000'x else begin
        XYOUTS,0.60,0.79,'We used rings of size '+in_ringsize+'" in the inner part',/normal,color='000000'x
        XYOUTS,0.60,0.77,'We used rings of size '+out_ringsize+'" in the outer part',/normal,color='000000'x
     ENDELSE
  ENDELSE
                                ;Currently GDL does not recognize true
                                ;type fonts yet. This leads to errors
                                ;in using the degree symbol. It also
                                ;does not yet recognize superscript
                                ;commands in tickmarks.
  tmpstr=strtrim(strsplit(filenames[1],'/',/extract),2)
  IF strupcase(tmpstr[0]) EQ 'MOMENTS' then begin
     mom0=readfits(filenames[1]+'.fits',mom0hed,/SILENT)
  ENDIF else begin
     mom0=readfits('Moments/'+filenames[1]+'.fits',mom0hed,/SILENT)
  ENDELSE
  mom0mod=readfits('Moments/Finalmodel_mom0.fits',mom0hedmod,/SILENT)
  mapmax=MAX(mom0,min=mapmin)
  buildaxii,mom0hed,xaxis,yaxis
  colormaps,'heat',/invert
  showpixelsmap,xaxis,yaxis,mom0,position=[0.15,0.1,0.35,0.1+0.2*scrdim[0]/scrdim[1]],/WCS, xtitle='RA (J2000)',ytitle='DEC (J2000)',BLANK_VALUE=0.,range=[0.,mapmax],/NOERASE,charthick=charthick,thick=thick,/black,/hex_color
  levels=[1E20, 2E20, 4E20, 8E20,16E20,32E20]
  beam=[sxpar(mom0hed,'BMAJ')*3600.,sxpar(mom0hed,'BMIN')*3600.]
  IF sxpar(mom0hed,'BPA') then bpa=sxpar(mom0hed,'BPA') else bpa=0
  columndensity,levels,double(vsys),beam,vwidth=1.,/NCOLUMN
  levels=levels/1000.
  loadct,0,/SILENT
  Contour,mom0,xaxis,yaxis,levels=levels,/overplot,c_colors=['000000'x]
  loadct,40,/SILENT
  Contour,mom0mod,xaxis,yaxis,levels=levels,/overplot,c_colors=['0000FF'x]
  beam_plot,beam[0],beam[1],bpa=bpa,center=[xaxis[0]-beam[0]/3600.,yaxis[0]+beam[0]/3600.],/fill,color='000000'x,/transparent
  colormaps,'heat',/invert
  colour_bar,[0.37,0.39],[0.12,0.1+0.2*scrdim[0]/scrdim[1]-0.02],strtrim(string(0,format='(F10.1)'),2),strtrim(string(mapmax,format='(F10.1)'),2),/OPPOSITE_LABEL,/BLACK,TITLE='(Jy bm!E-1!N x km s!E-1!N)',/VERTICAL,charthick=charthick,/hex_color
  loadct,0,/SILENT
  XYOUTS,0.45,0.01+0.2*scrdim[0]/scrdim[1],'Velocity Field, Moment0 and PV-Diagram along the major axis.',color='000000'x,/normal,charthick=charthick
  XYOUTS,0.45,0.01+0.2*scrdim[0]/scrdim[1]-0.02,'Black Contours: Data, White/Red Contours: Final Model',color='000000'x,/normal,charthick=charthick
  XYOUTS,0.45,0.01+0.2*scrdim[0]/scrdim[1]-0.04,'Moment 0 Contours are at 1, 2, 4, 8, 16, 32 x 10!E20!N cm!E-2',color='000000'x,/normal,charthick=charthick

;Velocity Field
  tmpstr=strtrim(strsplit(filenames[2],'/',/extract),2)
  IF strupcase(tmpstr[0]) EQ 'MOMENTS' then begin
     mom0=readfits(filenames[2]+'.fits',mom0hed,/SILENT)
  ENDIF else begin
     mom0=readfits('Moments/'+filenames[2]+'.fits',mom0hed,/SILENT)
  ENDELSE

  mom0mod=readfits('Moments/Finalmodel_mom1.fits',mom0hedmod,/SILENT)
  tmp=WHERE(FINITE(mom0mod))
  ;Too low inclination can result in a too small range
  ceninc=ceninc+12.5
  IF ceninc GT 90 then ceninc=90.
  velext=1.25*maxvrot*SIN((ceninc)*!DtoR)
  mapmax=double(vsys[0]+velext[0])
  mapmin=double(vsys[0]-velext[0])
  buildaxii,mom0hed,xaxis,yaxis
  colormaps,'sauron_colormap'
  showpixelsmap,xaxis,yaxis,mom0,position=[0.15,0.1+0.2*scrdim[0]/scrdim[1],0.35,0.1+0.4*scrdim[0]/scrdim[1]],/WCS,ytitle='DEC (J2000)',BLANK_VALUE=0.,range=[mapmin,mapmax],/NOERASE,charthick=charthick,thick=thick,xtickname=[' ',' ',' ',' ',' ',' ',' ',' ',' ',' '],/hex_color

  velostep=fix((fix(mapmax-mapmin)-fix(mapmax-mapmin)/10.)/10.)
  IF velostep LT 1 then velostep=1
  levels=(findgen(9)+0.5)*velostep+mapmin
  loadct,0,/SILENT
  Contour,mom0,xaxis,yaxis,levels=levels,/overplot,c_colors=['000000'x]
  Contour,mom0mod,xaxis,yaxis,levels=levels,/overplot,c_colors=['ffffff'x]
  beam_plot,beam[0],beam[1],bpa=bpa,center=[xaxis[0]-beam[0]/3600.,yaxis[0]+beam[0]/3600.],/fill,/transparent,color='000000'x
  colormaps,'sauron_colormap'
  colour_bar,[0.37,0.39],[0.1+0.2*scrdim[0]/scrdim[1]+0.02,0.1+0.4*scrdim[0]/scrdim[1]-0.02],strtrim(string(mapmin,format='(F10.1)'),2),strtrim(string(mapmax,format='(F10.1)'),2),/OPPOSITE_LABEL,/BLACK,TITLE='(km s!E-1!N)',/VERTICAL,charthick=charthick,/hex_color
  loadct,0,/SILENT
  XYOUTS,0.45,0.01+0.2*scrdim[0]/scrdim[1]-0.06,'Velocity Field Contours start at '+strtrim(string(levels[0],format='(F10.1)'),2)+' km s!E-1!N and increase with '+strtrim(string(velostep,format='(I10)'),2)+' km s!E-1!N.',color='000000'x,/normal,charthick=charthick
  
  
;PV Diagram along major axis
  IF fix(version) EQ version then begin
     spawn,'ls -1 PV-Diagrams/'+filenames[0]+'_2_xv.fits',mom0name
  ENDIF ELSE spawn,'ls -1 PV-Diagrams/'+filenames[0]+'_1_xv.fits',mom0name
  mom0=readfits(mom0name[n_elements(mom0name)-1],mom0hed,/SILENT)
  mom0mod=readfits('PV-Diagrams/Finalmodel_xv.fits',mom0hedmod,/SILENT)
  velbuf=(2.*velext)*0.2+disper/2.+2.*sxpar(mom0hed,'CDELT2')

  mapinmax=MAX(mom0[WHERE(FINITE(mom0))],min=mapinmin)
  if ABS(mapinmin/mapinmax) GT 0.2 then mapinmin=-1*mapinmax/5
  buildaxii,mom0hed,xaxis,yaxis
  IF mapmax+velbuf GT yaxis[n_elements(yaxis)-1] then tmp=yaxis[n_elements(yaxis)-1] else tmp= mapmax+velbuf
  IF mapmin-velbuf LT yaxis[0] then tmpmin=yaxis[0] else tmpmin=mapmin-velbuf
  yrange=[tmpmin,tmp]
  colormaps,'heat',/invert
  showpixelsmap,xaxis*3600.,yaxis,mom0,position=[0.65,0.1+0.2*scrdim[0]/scrdim[1],0.85,0.1+0.4*scrdim[0]/scrdim[1]], xtitle='Offset (arcsec)',ytitle='Velocity (km s!E-1!N)',BLANK_VALUE=0.,range=[mapinmin,mapinmax],/NOERASE,charthick=charthick,thick=thick,/black,yrange=yrange,/hex_color
  if n_elements(noise) LT 1 then noise=STDDEV(mom0[0:10,0:10])
  levels=[1,2,4,8,16,32,64,128]*1.5*noise
  levelsneg=([-2,-1])*1.5*noise
  loadct,0,/SILENT
  Contour,mom0,xaxis*3600.,yaxis,levels=levels,/overplot,c_colors=['000000'x]
  Contour,mom0,xaxis*3600,yaxis,levels=levelsneg,/overplot,c_colors=['999999'x],c_linestyle=2
  loadct,40,/silent
  Contour,mom0mod,xaxis*3600,yaxis,levels=levels,/overplot,c_colors=['0000FF'x]
  colormaps,'heat',/invert
  colour_bar,[0.87,0.89],[0.1+0.2*scrdim[0]/scrdim[1]+0.02,0.1+0.4*scrdim[0]/scrdim[1]-0.02],strtrim(string(mapinmin,format='(F10.4)'),2),strtrim(string(mapinmax,format='(F10.4)'),2),/OPPOSITE_LABEL,/BLACK,TITLE='(Jy bm!E-1!N)',/VERTICAL,charthick=charthick,/hex_color
  loadct,0,/SILENT
  XYOUTS,0.45,0.01+0.2*scrdim[0]/scrdim[1]-0.08,'PV-Diagram Contours start are at -3, -1.5, 1.5, 3, 6, 12, 24 x rms.',color='000000'x,/normal,charthick=charthick
  XYOUTS,0.45,0.01+0.2*scrdim[0]/scrdim[1]-0.1,'rms = '+strtrim(string(noise,format='(F10.5)'),2)+' Jy bm!E-1!N.',color='000000'x,/normal,charthick=charthick
  XYOUTS,0.45,0.01+0.2*scrdim[0]/scrdim[1]-0.12,'The distance used for conversions = '+strtrim(string(Distance,format='(F10.1)'),2)+' Mpc',color='000000'x,/normal,charthick=charthick
  IF ~(gdlidl) then image = tvrd(/true)
  DEVICE,/CLOSE  
  IF ~(gdlidl) then  write_png,'Overview.png',image
  IF gdlidl then begin
     spawn,"sed -i -- 's/1185 669/506 679/g' Overview.ps"
     spawn,"sed -i -- 's/1186 669/506 679/g' Overview.ps"
     IF FILE_TEST('Overview.ps--') then spawn,'rm -f Overview.ps--'
     spawn,'gs -help',result
     if n_elements(result) GT 1 then begin
        spawn,'gs -r300  -dEPSCrop -dTextAlphaBits=4 -sDEVICE=png16m -sOutputFile="Overview.png" -dBATCH -dNOPAUSE "Overview.ps"'
        spawn,'rm -f Overview.ps'
     ENDIF
  ENDIF
  converted=0
  spawn,'convert --help',result,notfound
  IF n_elements(result) GT 1 then begin
     gf=strsplit(result[0],' ',/extract)
     IF n_elements(gf) GT 1 then begin
        IF strtrim(gf[1],2) EQ 'ImageMagick' then begin
           converted=1
           spawn,'convert Overview.png -trim Overview.png'
        ENDIF
     ENDIF
  ENDIF
  IF converted LT 1 then begin
     spawn,'/usr/bin/convert --help',result,notfound
     IF n_elements(result) GT 1 then begin
        gf=strsplit(result[0],' ',/extract)
        IF n_elements(gf) GT 1 then begin
           IF strtrim(gf[1],2) EQ 'ImageMagick' then begin
              converted=1
              spawn,'/usr/bin/convert Overview.png -trim Overview.png'
           ENDIF
        ENDIF
     ENDIF
  ENDIF
  IF converted LT 1 then begin
     spawn,'imconvert --help',result,notfound
     IF n_elements(result) GT 1 then begin
        gf=strsplit(result[0],' ',/extract)
        IF n_elements(gf) GT 1 then begin
           IF strtrim(gf[1],2) EQ 'ImageMagick' then begin
              converted=1
              spawn,'imconvert Overview.png -trim Overview.png'
           ENDIF
        ENDIF
     ENDIF
  ENDIF
  IF converted LT 1 then begin
     spawn,'convert-im6 --help',result,notfound
     IF n_elements(result) GT 1 then begin
        gf=strsplit(result[0],' ',/extract)
        IF n_elements(gf) GT 1 then begin
           IF strtrim(gf[1],2) EQ 'ImageMagick' then begin
              converted=1
              spawn,'convert-im6 Overview.png -trim Overview.png'
           ENDIF
        ENDIF
     ENDIF
  ENDIF

  
end
