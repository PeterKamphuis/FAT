Pro int_profilev2,inputor,profile,PA=paor,MINOR=minor,HEADER=header,XCENTER=xcenter,YCENTER=ycenter,RANGE=range,AXIS=axis,ROTIMAGE=rotimage,NOAVERAGE=noaverage

;+
; NAME:
;       INT_PROFILEV2
;
; PURPOSE:
;       Program to make a radial or vertical profile from an image
;
; CATEGORY:
;       Support
; 
; CALLING SEQUENCE:
;       INT_PROFILEV2,inputor,profile, PA = paor, /MINOR,
;       HEADER = header, XCENTER = xcenter, YCENTER = ycenter, 
;       RANGE = range, AXIS = axis, ROTIMAGE = rotimage, /NOAVERAGE
;
;
; INPUTS:
;       inputor = 2D array with the values of the image pixels
;
; OPTIONAL INPUTS:
;       - 
;
; KEYWORD PARAMETERS:
;       /MINOR - extract the minor axis profile instead of the major axis
;       PA = position angle of the galaxy be aware that that PA is
;       defined from North so a galaxy aligned on the y-axis has a PA of 90
;       HEADER = Image's header can be given to rotate the image
;       around the images center or to give range in map units.
;       XCENTER= rotation center on xaxis in pixels. Supersedes the
;       header. Default is half the xaxis
;       YCENTER= rotation center on yaxis in pixels. Supersedes the
;       header default is half the yaxis
;       RANGE = the range over which to collapse the map/cube without
;       a header in pixels when header is given it assumed to be in
;       map units if not given the whole axis will be collapsed, when
;       2x2 array is assumed to be [xrange,yrange] 
;       AXIS = can be set to return a xaxis for your profile, if
;       header is given it will be based on header otherwise it is pixels offset from the center.
;       /NOAVERAGE - Set to not get an average profile but an
;                    integrated one 
;       ROTIMAGE = If set the rotated image will be returned in this array
;
; OUTPUTS:
;       profile = the extracted profile
;
; OPTIONAL OUTPUTS:
;       -
; 
; PROCEDURES CALLED:
;       ROT(), SIZE(), STRUPCASE(), STRTRIM() STRCOMPRESS(),
;       STR_SEP(), SXPAR(), TOTAL()
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
  input=inputor 
  IF n_elements(paor) GT 0 then pa=paor else pa=0
  arsize=SIZE(input)
  IF arsize[0] GT 2 then begin
     moments,input,inputmap
     input=dblarr(arsize[1],arsize[2])
     input=inputmap
     arsize=SIZE(input)
  ENDIF
  IF keyword_set(minor) then pa=pa+90

  IF n_elements(xcenter) EQ 0 then begin
     if n_elements(header) EQ 0 then begin
        xcenter=fix(arsize[1]/2.)
     endif else begin
        xcenter=sxpar(header,'CRPIX1')
     ENDELSE
  ENDIF
  IF n_elements(ycenter) EQ 0 then begin
     if n_elements(header) EQ 0 then begin
        ycenter=fix(arsize[2]/2.)
     endif else begin
        ycenter=sxpar(header,'CRPIX2')
     ENDELSE
  ENDIF

  if n_elements(header) NE 0 then begin
     xcenterdeg=(xcenter-sxpar(header,'CRPIX1'))*sxpar(header,'CDELT1')+sxpar(header,'CRVAL1')
     ycenterdeg=(ycenter-sxpar(header,'CRPIX2'))*sxpar(header,'CDELT2')+sxpar(header,'CRVAL2')
  ENDIF

  IF pa NE 90 then begin
     result=ROT(input, pa-90,1.0,xcenter,ycenter,CUBIC=-0.5,/PIVOT)
     input=result
  ENDIF


  IF n_elements(header) GT 0 then begin  
     buildaxii,header,xaxis,yaxis
  ENDIF
  pixpos=dblarr(2,2)
  IF n_elements(range) GT 0 then begin 
     rangesize=SIZE(range)
     if n_elements(header) GT 0 then begin  
        IF rangesize[0] EQ 1 then begin
           
           IF keyword_set(minor) then begin
              pixpos[0,0]=0
              pixpos[1,0]=n_elements(input[0,*])-1
              IF xaxis[0]-xaxis[1] LT 0 then begin
                 tmp=WHERE(xaxis GT range[0]+xcenterdeg)
                 pixpos[0,1]=tmp[0]
                 tmp=WHERE(xaxis LT range[1]+xcenterdeg)
                 pixpos[1,1]=tmp[n_elements(tmp)-1]
              endif else begin
                 tmp=WHERE(xaxis LT range[0]+xcenterdeg)
                 pixpos[1,1]=tmp[0]
                 tmp=WHERE(xaxis GT range[1]+xcenterdeg)
                 pixpos[0,1]=tmp[n_elements(tmp)-1]
              endelse


           ENDIF ELSE BEGIN
              pixpos[0,0]=0
              pixpos[1,0]=n_elements(input[*,1])-1
              tmp=WHERE(yaxis GT range[0]+ycenterdeg)
              pixpos[0,1]=tmp[0]
              tmp=WHERE(yaxis LT range[1]+ycenterdeg)
              
              pixpos[1,1]=tmp[n_elements(tmp)-1]
              
           ENDELSE
           
        ENDIF
        IF rangesize[0] EQ 2 then begin
           
           IF keyword_set(minor) then begin
              IF xaxis[0]-xaxis[1] LT 0 then begin
                 tmp=WHERE(xaxis GT range[0,0]+xcenterdeg)
                 pixpos[0,1]=tmp[0]
                 tmp=WHERE(xaxis LT range[1,0]+xcenterdeg)
                 pixpos[1,1]=tmp[n_elements(tmp)-1]
              endif else begin
                 tmp=WHERE(xaxis LT range[0,0]+xcenterdeg)
                 pixpos[1,1]=tmp[0]
                 tmp=WHERE(xaxis GT range[1,0]+xcenterdeg)
                 pixpos[0,1]=tmp[n_elements(tmp)-1]
              endelse
              tmp=WHERE(yaxis GT range[0,1]+ycenterdeg)
              pixpos[0,0]=tmp[0]
              tmp=WHERE(yaxis LT range[1,1]+ycenterdeg)
              pixpos[1,0]=tmp[n_elements(tmp)-1]
           ENDIF ELSE BEGIN
              IF xaxis[0]-xaxis[1] LT 0 then begin
                 tmp=WHERE(xaxis GT range[0,0]+xcenterdeg)
                 pixpos[0,0]=tmp[0]
                 tmp=WHERE(xaxis LT range[1,0]+xcenterdeg)
                 pixpos[1,0]=tmp[n_elements(tmp)-1]
              endif else begin
                 tmp=WHERE(xaxis LT range[0,0]+xcenterdeg)
                 pixpos[1,0]=tmp[0]
                 tmp=WHERE(xaxis GT range[1,0]+xcenterdeg)
                 pixpos[0,0]=tmp[n_elements(tmp)-1]
              endelse
              tmp=WHERE(yaxis GT range[0,1]+ycenterdeg)
              pixpos[0,1]=tmp[0]
              tmp=WHERE(yaxis LT range[1,1]+ycenterdeg)
              pixpos[1,1]=tmp[n_elements(tmp)-1]
           ENDELSE
           
        ENDIF
     ENDIF ELSE BEGIN
        IF rangesize[0] EQ 1 then begin
           pixpos[0,0]=0
           pixpos[1,0]=n_elements(input[*,1])-1
           pixpos[0,1]=range[0]
           pixpos[1,1]=range[1]
        ENDIF
        IF rangesize[0] EQ 2 then begin
           IF keyword_set(minor) then begin
              pixpos[0,1]=range[0,0]
              pixpos[1,1]=range[1,0]
              pixpos[0,0]=range[0,1]
              pixpos[1,0]=range[1,1]
           ENDIF ELSE BEGIN
              pixpos[0,0]=range[0,0]
              pixpos[1,0]=range[1,0]
              pixpos[0,1]=range[0,1]
              pixpos[1,1]=range[1,1]
           ENDELSE
        ENDIF
     ENDELSE
  ENDIF ELSE BEGIN
     
     pixpos[0,*]=0
     pixpos[1,0]=n_elements(input[*,1])-1
     pixpos[1,1]=n_elements(input[1,*])-1
  ENDELSE
  IF pixpos[0,0] GT pixpos[1,0] then pixpos[*,0]=reverse(pixpos[*,0])
  IF pixpos[0,1] GT pixpos[1,1] then pixpos[*,1]=reverse(pixpos[*,1])


  profile=dblarr(pixpos[1,0]+1-pixpos[0,0])
  clear=dblarr(n_elements(input[*,0]),n_elements(input[0,*]))
  tmp2=WHERE(FINITE(input))
  clear[tmp2]=input[tmp2]
  tmp=dblarr(fix(pixpos[1,1]-pixpos[0,1]))

  for i=pixpos[0,0],pixpos[1,0] do begin
     
     tmp=clear[i,pixpos[0,1]:pixpos[1,1]]
     tmp2=WHERE(tmp EQ 0)
     IF tmp2[0] EQ -1 then begin
        IF keyword_set(noaverage) then begin
           profile[i-pixpos[0,0]]=TOTAL(clear[i,pixpos[0,1]:pixpos[1,1]])/(n_elements(clear[i,pixpos[0,1]:pixpos[1,1]]))
        ENDIF ELSE BEGIN         
           profile[i-pixpos[0,0]]=TOTAL(clear[i,pixpos[0,1]:pixpos[1,1]])/(n_elements(clear[i,pixpos[0,1]:pixpos[1,1]]))
        ENDELSE
     ENDIF ELSE BEGIN
        IF keyword_set(noaverage) then begin
           profile[i-pixpos[0,0]]=TOTAL(clear[i,pixpos[0,1]:pixpos[1,1]])
        ENDIF ELSE BEGIN         
           profile[i-pixpos[0,0]]=TOTAL(clear[i,pixpos[0,1]:pixpos[1,1]])/(n_elements(clear[i,pixpos[0,1]:pixpos[1,1]])-n_elements(tmp2))
        ENDELSE
     ENDELSE
  endfor
  tmpindex=where(FINITE(profile) EQ 0. )
  IF tmpindex[0] NE -1 then profile[where(FINITE(profile) EQ 0. )]=0.
  axis=dblarr(n_elements(profile))
  IF n_elements(header) GT 0 then begin
     IF keyword_set(minor) then axis[*]=yaxis[pixpos[0,0]:pixpos[1,0]]-ycenterdeg else axis[*]=xaxis[pixpos[0,0]:pixpos[1,0]]-xcenterdeg
  ENDIF ELSE begin
     if keyword_set(minor) then axis=findgen(n_elements(axis))+pixpos[0,0]-ycenter else  axis=findgen(n_elements(axis))+pixpos[0,0]-xcenter
  ENDELSE
  if n_elements(rotimage) GT 0 then begin
     rotimage=dblarr(n_elements(input[*,0]),n_elements(input[0,*]))
     rotimage=input
  ENDIF

end


