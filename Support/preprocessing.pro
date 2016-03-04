Pro preprocessing,cube,header,writecube,log=log,catalogue=outputcatalogue,noise=noise,name=name,directory=dir


;+
; NAME:
;       PREPROCESSING
;
; PURPOSE:
;       Clean up the cube  and make sure it has all the right
;       characteristics and blanks
;
; CATEGORY:
;       Support
; 
; CALLING SEQUENCE:
;        PREPROCESSING,cube,headerwritecube,log=log
;
;
; INPUTS:
;         cube = the array containing the input cube.  
;       header = the header of the cube.
;    writecube = trigger to determine whether the cube is modified and
;                we should write it back to disk. (0=no, 1=yes)
;          log = the logging file.
;
; OPTIONAL INPUTS:
;       
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
;       FILE_TEST(),LINENUMBER(),WHERE(),FINITE(),SXPAR(),SXADDPAR(),SIGMA(),N_ELEMENTS()
;
; EXAMPLE:
;      
;
; MODIFICATION HISTORY:
;       18-02-2016 P.Kamphuis; Replaced sigma with STDDEV 
;       Written 04-01-2016 P.Kamphuis v1.0
;
; NOTE:
;     
;-  

                                ; check for 0 values, they
                                ; shouldn't be present and mess
                                ; things up
                                ;IF present assume they are blanks and
                                ;add warning

  IF sxpar(header,'CDELT3') LT 0 then begin
     IF size(log,/TYPE) EQ 7 then begin
        openu,66,log,/APPEND
        printf,66,linenumber()+'PREPROCESSING: Your velocity axis is declining with increasing channels'
        printf,66,linenumber()+'PREPROCESSING: We reversed the velocity axis'   
        close,66
     ENDIF ELSE BEGIN
        print,linenumber()+'PREPROCESSING: Your velocity axis is declining with increasing channels'
        print,linenumber()+'PREPROCESSING: We reversed the velocity axis'    
     ENDELSE
     sxaddpar,header,'CDELT3',ABS(sxpar(header,'CDELT3'))
     cube=reverse(cube,3)
     writecube=1
  ENDIF
  
  
  tmp=WHERE(cube EQ 0.d)
  IF tmp[0] NE -1 then begin
     IF size(log,/TYPE) EQ 7 then begin
        openu,66,log,/APPEND
        printf,66,linenumber()+'PREPROCESSING: Your cube contains values exactly 0. If this is padding these should be blanks.'
        printf,66,linenumber()+'PREPROCESSING: We have changed them to blanks.'
        close,66
     ENDIF ELSE BEGIN
        print,linenumber()+'PREPROCESSING: Your cube contains values exactly 0. If this is padding these should be blanks.'
        print,linenumber()+'PREPROCESSING: We have changed them to blanks.'
     ENDELSE
     cube[tmp]=!values.f_nan
     IF writecube EQ 0 then writecube=1
  ENDIF
                                ;first we check wether the cube has some proper noise statistics by
                                ;comparing the last and first channel, if they differ to much we cut
                                ;off the high noise channel
  difference = 10.
  changedcube=0.
  firstcompprev=0.
  secondcompprev=0.
  firstcut=0
  secondcut=0
  WHILE difference GT 1 do begin
     tmpnoblank=cube[*,*,0]
     wherefinite=WHERE(FINITE(tmpnoblank))
     WHILE wherefinite[0] EQ -1 DO begin
        tmp=fltarr(n_elements(cube[*,0,0]),n_elements(cube[0,*,0]),n_elements(cube[0,0,*])-1)
        tmp[*,*,0:n_elements(tmp[0,0,*])-1]=cube[*,*,1:n_elements(cube[0,0,*])-1]
        sxaddpar,header,'CRPIX3',sxpar(header,'CRPIX3')-1.
        cube=fltarr(n_elements(tmp[*,0,0]),n_elements(tmp[0,*,0]),n_elements(tmp[0,0,*]))
        IF n_elements(cube[0,0,*]) LT 5 then begin
           IF size(log,/TYPE) EQ 7 then begin
              openu,66,log,/APPEND
              printf,66,linenumber()+'PREPROCESSING: '+dir+'/'+name+' has too many blanked channels.'
              close,66
           ENDIF ELSE BEGIN
              printf,linenumber()+'PREPROCESSING: '+dir+'/'+name+'has too many blanked channels.'
           ENDELSE         
           openu,1,outputcatalogue,/APPEND
           printf,1,format='(A60,A90)', dir,'The Cube has too many blanked channels'
           close,1
           writecube=2
           goto,finishup
        ENDIF
        IF size(log,/TYPE) EQ 7 then begin
           openu,66,log,/APPEND
           printf,66,linenumber()+'PREPROCESSING: We are cutting the cube as the first channel is completely blank.'
           close,66
        ENDIF
        sxaddpar,header,'NAXIS3',fix(sxpar(header,'NAXIS3')-1)
        IF writecube EQ 0 then writecube=1
        cube=tmp
        tmpnoblank=cube[*,*,0]
        wherefinite=WHERE(FINITE(tmpnoblank))
     ENDWHILE
     rmsfirstchannel=STDDEV(tmpnoblank[WHERE(FINITE(tmpnoblank))])          
     tmpnoblank=cube[*,*,n_elements(cube[0,0,*])-1]
     wherefinite=WHERE(FINITE(tmpnoblank))
     WHILE wherefinite[0] EQ -1 DO begin
        tmp=fltarr(n_elements(cube[*,0,0]),n_elements(cube[0,*,0]),n_elements(cube[0,0,*])-1)
        tmp[*,*,0:n_elements(tmp[0,0,*])-1]=cube[*,*,0:n_elements(cube[0,0,*])-2]
        cube=fltarr(n_elements(tmp[*,0,0]),n_elements(tmp[0,*,0]),n_elements(tmp[0,0,*]))
        IF n_elements(cube[0,0,*]) LT 5 then begin
           IF size(log,/TYPE) EQ 7 then begin
              openu,66,log,/APPEND
              printf,66,linenumber()+'PREPROCESSING: '+dir+'/'+name+' has too many blanked channels.'
              close,66
           ENDIF ELSE BEGIN
              printf,linenumber()+'PREPROCESSING: '+dir+'/'+name+'has too many blanked channels.'
           ENDELSE         
           openu,1,outputcatalogue,/APPEND
           printf,1,format='(A60,A90)', dir,'The cube has too many blanked channels'
           close,1    
           writecube=2
           goto,finishup
        ENDIF
        IF size(log,/TYPE) EQ 7 then begin
           openu,66,log,/APPEND
           printf,66,linenumber()+'PREPROCESSING: We are cutting the cube as the last channel is completely blank.'
           close,66
        ENDIF
        sxaddpar,header,'NAXIS3',fix(sxpar(header,'NAXIS3')-1)
        changedcube=1
        cube=tmp
        tmpnoblank=cube[*,*,n_elements(cube[0,0,*])-1]
        wherefinite=WHERE(FINITE(tmpnoblank))
     ENDWHILE
     rmslastchannel=STDDEV(tmpnoblank[WHERE(FINITE(tmpnoblank))])
     IF rmsfirstchannel EQ 0. then rmsfirstchannel=2.*rmslastchannel
     IF rmslastchannel EQ 0. then rmslastchannel=2.*rmsfirstchannel
     tmpnoblank=cube[0:5,0:5,*]
     rmsbottoml=STDDEV(tmpnoblank[WHERE(FINITE(tmpnoblank))])
     tmpnoblank=cube[n_elements(cube[*,0,0])-6:n_elements(cube[*,0,0])-1,0:5,*]
     rmsbottomr=STDDEV(tmpnoblank[WHERE(FINITE(tmpnoblank))])
     tmpnoblank=cube[n_elements(cube[*,0,0])-6:n_elements(cube[*,0,0])-1,n_elements(cube[0,*,0])-6:n_elements(cube[0,*,0])-1,*]
     rmstopr=STDDEV(tmpnoblank[WHERE(FINITE(tmpnoblank))])
     tmpnoblank=cube[0:5,n_elements(cube[0,*,0])-6:n_elements(cube[0,*,0])-1,*]
     rmstopl=STDDEV(tmpnoblank[WHERE(FINITE(tmpnoblank))])
     rmschan=(rmsfirstchannel+rmslastchannel)/2.
     rmscorn=(rmsbottoml+rmsbottomr+rmstopl+rmstopr)/4.
     diff=ABS((rmsfirstchannel-rmslastchannel)/rmsfirstchannel)
     diff2=ABS((rmschan-rmscorn)/rmschan)
     IF diff LT 0.2 AND FINITE(diff) AND diff2 LT 0.2 then difference=0. else begin
        IF size(log,/TYPE) EQ 7 then begin
           openu,66,log,/APPEND
           printf,66,linenumber()+'PREPROCESSING: We are cutting the cube as clearly the noise statistics are off.'
           printf,66,linenumber()+'PREPROCESSING: Noise in the first channel is '+string(rmsfirstchannel)
           printf,66,linenumber()+'PREPROCESSING: Noise in the last channel is '+string(rmslastchannel)
           printf,66,linenumber()+'PREPROCESSING: Noise in the corners is '+string(rmscorn)
           close,66
        ENDIF ELSE BEGIN
           print,linenumber()+'PREPROCESSING: We are cutting the cube as clearly the noise statistics are off.'
           print,linenumber()+'PREPROCESSING: Noise in the first channel is '+string(rmsfirstchannel)
           print,linenumber()+'PREPROCESSING: Noise in the last channel is '+string(rmslastchannel)       
        ENDELSE
        IF writecube EQ 0 then writecube=1.
        tmp=fltarr(n_elements(cube[*,0,0]),n_elements(cube[0,*,0]),n_elements(cube[0,0,*])-1)
        firstcomp=ABS((rmsfirstchannel-rmscorn)/rmscorn)
        secondcomp=ABS((rmslastchannel-rmscorn)/rmscorn)
        IF (firstcomp GT secondcomp AND ABS((firstcomp-firstcompprev)/firstcomp) GT ABS((secondcomp-secondcompprev)/secondcomp) AND firstcut LT 8.) OR NOT FINITE(rmsfirstchannel) then begin
           firstcompprev=firstcomp
           firstcut++
           tmp[*,*,0:n_elements(tmp[0,0,*])-1]=cube[*,*,1:n_elements(cube[0,0,*])-1]
           sxaddpar,header,'CRPIX3',sxpar(header,'CRPIX3')-1.
        ENDIF ELSE BEGIN
           secondcompprev=secondcomp
           secondcut++
           tmp[*,*,0:n_elements(tmp[0,0,*])-1]=cube[*,*,0:n_elements(cube[0,0,*])-2]
        endelse
        IF firstcut GE 8 AND secondcut GE 8 then begin
           difference=0.
           IF size(log,/TYPE) EQ 7 then begin
              openu,66,log,/APPEND
              printf,66,linenumber()+'PREPROCESSING: '+dir+'/'+name+' has non-uniform noise statistics.'
              close,66
           ENDIF ELSE BEGIN
              print,linenumber()+'PREPROCESSING: '+dir+'/'+name+'  has non-uniform noise statistics.'
           ENDELSE     
        ENDIF
        sxaddpar,header,'NAXIS3',fix(sxpar(header,'NAXIS3')-1)
        cube=fltarr(n_elements(tmp[*,0,0]),n_elements(tmp[0,*,0]),n_elements(tmp[0,0,*]))
        cube=tmp
        IF n_elements(cube[0,0,*]) LT 5 then begin
           IF size(log,/TYPE) EQ 7 then begin
              openu,66,log,/APPEND
              printf,66,linenumber()+'PREPROCESSING: '+dir+'/'+name+' has noise statistics that cannot be dealt with.'
              close,66
           ENDIF ELSE BEGIN
              print,linenumber()+'PREPROCESSING: '+dir+'/'+name+' has noise statistics that cannot be dealt with.'
           ENDELSE         
           openu,1,outputcatalogue,/APPEND
           printf,1,format='(A60,A90)', dir,'The Cube has noise statistics that cannot be dealt with'
           close,1    
           writecube=2
           goto,finishup
        ENDIF
     ENDELSE
  ENDWHILE
  noise=rmscorn
  IF writecube EQ 1. then begin        
     IF size(log,/TYPE) EQ 7 then begin
        openu,66,log,/APPEND
        printf,66,linenumber()+'PREPROCESSING: The cube was cut due to incorrect noise statistics.'
        close,66
     ENDIF ELSE BEGIN
        print,linenumber()+'PREPROCESSING: The cube was cut due to incorrect noise statistics.'    
     ENDELSE
     writecube=3
  ENDIF
                                ; And last we want to cut out central absorption parts. We will cut
                                ; below negative 10 sigma. If this messes up the process there is
                                ; something seriously wrong with the cube                            
  tmp=WHERE(cube LT -10.*noise)
  IF tmp[0] NE -1 then begin
     cube[tmp]=!values.f_nan
     writecube=3
     IF size(log,/TYPE) EQ 7 then begin
        openu,66,log,/APPEND
        printf,66,linenumber()+'PREPROCESSING: Your cube had values below -10*sigma. If you do not have a central absorption source there is something seriously wrong with the cube'
        printf,66,linenumber()+'PREPROCESSING: We blanked these values.'
        close,66
     ENDIF ELSE BEGIN
        print,linenumber()+'PREPROCESSING: Your cube had values below -10.*sigma. If you do not have a central absorption source there is something seriously wrong with the cube.'        
     ENDELSE
  ENDIF
                                ;if the moment maps or the proper binmask do not yet exist we create
                                ;them here
  IF size(log,/TYPE) EQ 7 then begin
     openu,66,log,/APPEND
     printf,66,linenumber()+'PREPROCESSING: The noise in the cube is '+string(noise)
     close,66
  ENDIF ELSE BEGIN
     print,linenumber()+'PREPROCESSING: The noise in the cube is '+string(noise)
  ENDELSE   
finishup:
end
