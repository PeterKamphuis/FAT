Pro extract_pv,Cube,header,pa,xv,center=center,xvheader=new_header,width=width

;width should be in the same units as the pixel sizes



inheader=header
if n_elements(center) EQ 0 then center=[sxpar(inheader,'CRVAL1'), sxpar(inheader,'CRVAL2')]
xv=dblarr(n_elements(Cube[*,0,0]),n_elements(Cube[0,0,*]))
ypix=sxpar(inheader,'CRPIX2')+(center[1]-sxpar(inheader,'CRVAL2'))/sxpar(inheader,'CDELT2')
xpix=sxpar(inheader,'CRPIX1')+(center[0]-sxpar(inheader,'CRVAL1'))/sxpar(inheader,'CDELT1')*COS(center[1]*!DtoR)
xrange=[-xpix,sxpar(inheader,'NAXIS1')-xpix]
print,'it is here',xpix,sxpar(inheader,'NAXIS1')-xpix-1
IF xpix GT sxpar(inheader,'NAXIS1')-xpix-1 then xsize=floor(sxpar(inheader,'NAXIS1')-xpix-1) else xsize=floor(xpix)
print,xsize
xv=dblarr(fix(2*xsize),n_elements(Cube[0,0,*]))
newCube=dblarr(n_elements(Cube[*,0,0]),n_elements(Cube[0,*,0]),n_elements(Cube[0,0,*]))
for j=0,n_elements(Cube[0,0,*])-1 do begin
   Newcube[*,*,j]=ROT(Cube[*,*,j],pa-90.,1.0,xpix,ypix,missing=!values.f_nan,cubic=-1,/PIVOT) ;For HI
endfor
print,n_elements(xv[*,0]),n_elements(xv[*,0]),xpix,xsize,n_elements(Cube[*,0,0])
print,fix(xpix-xsize),fix(xpix+xsize-1)
If width/(sxpar(inheader,'CDELT2')) LT 1 then begin
   xv[*,*]=Newcube[fix(xpix-xsize):fix(xpix+xsize-1),round(ypix),*]
   print,'We should not be doing this'
ENDIF ELSE BEGIN
   centrpix=xpix
   istart=fix(xpix-xsize)
   
   IF istart lt 0 then begin
    ;  centrpix=xpix-istart
      istart=0
      print,'here we look',xpix,centrpix
      stop
   ENDIF
   IF istart Ne 0 then centrpix=xpix-istart
   iend=fix(xpix+xsize-1)
   IF iend GT n_elements(xv[*,0])-1 then iend=n_elements(xv[*,0])-1
   print,istart,iend,xpix,centrpix,'is here',xsize,n_elements(xv[*,0])
;  for i=istart,iend do begin
  ;    for j=0,n_elements(Cube[0,0,*])-1 do begin
    ;     print,n_elements(xv[*,0]),n_elements(xv[0,*]),i,j
   
         xv[0:iend-istart,*]=SUM(Newcube[istart:iend,fix(ypix-(width)/(2.*(sxpar(inheader,'CDELT2')))):fix(ypix+(width)/(2.*(sxpar(inheader,'CDELT2')))),*],1)/n_elements(Newcube[0,fix(ypix-(width)/(2.*(sxpar(inheader,'CDELT2')))):fix(ypix+(width)/(2.*(sxpar(inheader,'CDELT2')))),0])

      ;   xv[i,]=TOTAL(Newcube[i,fix(ypix-(width)/(2.*(sxpar(inheader,'CDELT2')))):fix(ypix+(width)/(2.*(sxpar(inheader,'CDELT2')))),j])/n_elements(Newcube[i,fix(ypix-(width)/(2.*(sxpar(inheader,'CDELT2')))):fix(ypix+(width)/(2.*(sxpar(inheader,'CDELT2')))),j])
   ;   endfor
  ; endfor
ENDELSE
IF sxpar(inheader,'CDELT1') LT 0 then BEGIN
   sxaddpar,inheader,'CDELT1',ABS(sxpar(inheader,'CDELT1'))
   xv=REVERSE(xv)
ENDIF

new_header=inheader
sxaddpar,new_header,'NAXIS',2
sxaddpar,new_header,'NAXIS1',n_elements(xv[*,0])
sxaddpar,new_header,'NAXIS2',n_elements(xv[0,*])
sxaddpar,new_header,'CDELT1',sxpar(inheader,'CDELT1')
sxaddpar,new_header,'CRVAL1',0.
sxaddpar,new_header,'CRPIX1',centrpix
sxaddpar,new_header,'CTYPE1','ANGLE'
sxaddpar,new_header,'CUNIT1','DEGREE',AFTER='CRPIX1'
sxaddpar,new_header,'CDELT2',sxpar(inheader,'CDELT3')
sxaddpar,new_header,'CRVAL2',sxpar(inheader,'CRVAL3')
sxaddpar,new_header,'CRPIX2',sxpar(inheader,'CRPIX3')
sxaddpar,new_header,'CTYPE2',sxpar(inheader,'CTYPE3')
sxaddpar,new_header,'NAXIS2',sxpar(inheader,'NAXIS3')
sxaddpar,new_header,'PA',pa,AFTER='CTYPE2'
IF width/(sxpar(inheader,'CDELT2')) GT 1 then sxaddpar,new_header,'Strip Width',width,AFTER='PA'
sxaddpar,new_header,'RA POS',center[0],AFTER='PA'
sxaddpar,new_header,'DEC POS',center[1],AFTER='RA POS'
;sxaddpar,new_header,'BUNIT','Jy/Beam',AFTER='RA POS'


sxdelpar,new_header,['NAXIS3','CDELT3','CRPIX3','CRVAL3','CTYPE3','LTYPE']


IF sxpar(inheader,'CUNIT3') then begin
   sxaddpar,new_header,'CUNIT2',sxpar(inheader,'CUNIT3'),after='CTYPE2'
   sxdelpar,new_header,'CUNIT3'
ENDIF

;SET_PLOT, 'X', /COPY
;WINDOW, 0, XSIZE = 400, YSIZE = 400
;TVIMAGE,bytscl(xv),pos=[0.2,0.2,0.9,0.9]




end






