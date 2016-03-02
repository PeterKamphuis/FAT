Function FAT_SMOOTH,map,sigma,GAUSSIAN=Gaussian,BOX=box

;+
; NAME:
;       FAT_SMOOTH
;
; PURPOSE:
;       Routine to apply a gaussian or box averaging smooth kernel to an image
;
; CATEGORY:
;       Support
; 
; CALLING SEQUENCE:
;       FAT_SMOOTH,map,sigma
;
;
; INPUTS:
;       map = moment 1 map of the galaxy
;     sigma = the standard deviation of the gaussian to be used, the
;             gaussian kernel will have the size of 10 times the FWHM defined
;             by this sigma.
; OPTIONAL INPUTS:
;       - 
;
; KEYWORD PARAMETERS:
;       /GAUSSIAN - Gaussian Smoothing
;       /BOX      - Box averaging  
;
; OUTPUTS:
;       Result = the smoothed map 
;
; OPTIONAL OUTPUTS:
;       -
;
; PROCEDURES CALLED:
;       
;
; EXAMPLE:
;      
;
; MODIFICATION HISTORY: 
;       Written 29-12-2016 by P.Kamphuis
;
; NOTE:
;     
;-
  
mapor=map

if not keyword_set(gaussian) and not keyword_set(box) then gaussian=1

if keyword_set(gaussian) then begin
   xgauss=fltarr(fix(3*2.3548*sigma))
   ygauss=fltarr(fix(3*2.3548*sigma))
   xgauss=findgen(fix(3*2.3548*sigma))-fix(3*2.3548*sigma/2.)
   ygauss=findgen(fix(3*2.3548*sigma))-fix(3*2.3548*sigma/2.)
   Gauss=fltarr(fix(3*2.3548*sigma),fix(3*2.3548*sigma))
   for i=0,fix(3*2.3548*sigma)-1 do begin
                                ;From http://mathworld.wolfram.com/GaussianFunction.html
      Gauss[i,*]=(1./(2.*!pi*sigma^2))*EXP(-((0.5*xgauss[i]^2)/(sigma^2))-((0.5*ygauss[*]^2)/(sigma^2)))
   endfor
   tmp_smoothed_map=dblarr(fix(3*2.3548*sigma),fix(3*2.3548*sigma))
   startpixx=fix(3*2.3548*sigma/2.)
   endpixx=n_elements(map[*,0])-startpixx-1
   startpixy=fix(3*2.3548*sigma/2.)
   endpixy=n_elements(map[0,*])-startpixy-1
   for xax=startpixx,endpixx do begin
      startx=xax-startpixx
      endx=xax+startpixx
      for yax=startpixy,endpixy do begin
         starty=yax-startpixy
         endy=yax+startpixy
         map[xax,yax]=TOTAL(mapor[startx:endx,starty:endy]*Gauss[*,*])
      endfor
   endfor
   return,map
endif
if keyword_set(box) then begin
   inf=WHERE(FINITE(mapor) EQ 0)
   mapor[inf]=0.
   startpixx=fix(sigma/2.)
   endpixx=n_elements(map[*,0])-startpixx-1
   startpixy=fix(sigma/2.)
   endpixy=n_elements(map[0,*])-startpixy-1
   for xax=startpixx,endpixx do begin
      startx=xax-startpixx
      endx=xax+startpixx
      for yax=startpixy,endpixy do begin
         starty=yax-startpixy
         endy=yax+startpixy
         zero=WHERE(mapor[startx:endx,starty:endy] EQ 0)
         if zero[0] EQ -1 then map[xax,yax]=TOTAL(mapor[startx:endx,starty:endy])/n_elements(mapor[startx:endx,starty:endy]) else begin
            IF TOTAL(mapor[startx:endx,starty:endy]) EQ 0. then map[xax,yax]=!values.f_nan else begin
               map[xax,yax]=TOTAL(mapor[startx:endx,starty:endy])/(n_elements(mapor[startx:endx,starty:endy])-n_elements(zero))
            ENDELSE
         ENDELSE
      endfor
   endfor
   return,map
endif
   
end
