Pro clean_header,header,writecube,beam,log=log,catalogue=outputcatalogue,directory=dir

;+
; NAME:
;       CLEAN_HEADER
;
; PURPOSE:
;       Clean up the cube header and make sure it has all the right
;       variables that we require in the process of fitting
;
; CATEGORY:
;       Support
; 
; CALLING SEQUENCE:
;       CLEAN_HEADER,header,log=log
;
;
; INPUTS:
;       header = the header of the cube
;    writecube = trigger to determine whether the cube is modified and
;                we should write it back to disk. (0=no, 1=yes)
;          log = the logging file
;
; OPTIONAL INPUTS:
;       
;
; KEYWORD PARAMETERS:
;       -
;
; OUTPUTS:
;          beam = the beam recorded in the header. If no beam is
;          present it will return [NaN,Nan] and FAT will abort
; OPTIONAL OUTPUTS:
;       -
; 
; PROCEDURES CALLED:
;       SXPAR(),SXDELPAR(),SXADDPAR(),STRUPCASE()
;
; EXAMPLE:
;      
;
; MODIFICATION HISTORY:
;       11-07-2018 P. Kamphuis; Replace the bitwise not with the
;                               logical negation ~ to avoid problems
;                               with even integers. Additionally added
;                               a check that there are more than 2
;                               pixels per beam and that the beam is
;                               smaller than the cube.       
;       12-06-2016 P. Kamphuis; Fixed the condition for detecting a
;       frequency axis to exit cleanly and check both unit and type.  
;       Written 04-01-2016 P.Kamphuis v1.0
;
; NOTE:
;     
;-

  howmanyaxis=sxpar(header,'NAXIS')
  writecube=0
  IF howmanyaxis GT 3 then begin
     sxaddpar,header,'NAXIS',3
                                ;whenever we change something we want to rewrite the cube
     writecube=1
  ENDIF
  IF sxpar(header,'NAXIS4') then begin
     sxdelpar,header,'NAXIS4'
     writecube=1
  ENDIF
  IF ~sxpar(header,'EPOCH') then begin
     IF sxpar(header,'EQUINOX') then begin
        IF size(log,/TYPE) EQ 7 then begin
           openu,66,log,/APPEND
           printf,66,linenumber()+'CLEAN_HEADER: Your cube has no EPOCH keyword but we found EQUINOX.'
           printf,66,linenumber()+'CLEAN_HEADER: We have set EPOCH to '+string(sxpar(header,'EQUINOX'))
           close,66
        ENDIF ELSE BEGIN
           print,linenumber()+'CLEAN_HEADER: Your cube has no EPOCH keyword but we found EQUINOX.'
           print,linenumber()+'CLEAN_HEADER: We have set EPOCH to '+string(sxpar(header,'EQUINOX'))   
        ENDELSE
        sxaddpar,header,'EPOCH',sxpar(header,'EQUINOX')
     ENDIF ELSE BEGIN
         IF size(log,/TYPE) EQ 7 then begin
           openu,66,log,/APPEND
           printf,66,linenumber()+'CLEAN_HEADER: Your cube has no EPOCH keyword'
           printf,66,linenumber()+'CLEAN_HEADER: We assumed J2000'
           close,66
        ENDIF ELSE BEGIN
           print,linenumber()+'CLEAN_HEADER: Your cube has no EPOCH keyword'
           print,linenumber()+'CLEAN_HEADER: We assumed J2000'
        ENDELSE
        sxaddpar,header,'EPOCH',2000.
     ENDELSE  
     writecube=1
  ENDIF
  
 
  channelwidth=ABS(sxpar(header,'CDELT3'))   
  veltype=strtrim(strcompress(sxpar(header,'CUNIT3')))
  velproj=sxpar(header,'CTYPE3')	
  IF STRUPCASE(veltype) EQ 'HZ' OR  STRUPCASE(strtrim(velproj,2)) EQ 'FREQ' then begin
     IF size(log,/TYPE) EQ 7 then begin
        openu,66,log,/APPEND
        printf,66,linenumber()+'CLEAN_HEADER: FREQUENCY IS NOT A SUPPORTED VELOCITY AXIS.'          
        close,66
     ENDIF ELSE BEGIN
        print,linenumber()+'CLEAN_HEADER: FREQUENCY IS NOT A SUPPORTED VELOCITY AXIS.'    
     ENDELSE
     openu,1,outputcatalogue,/APPEND
     printf,1,format='(A60,2A12,A120)',Dir,0.,0.,'The Cube has frequency as a velocity axis this is not supported'
     close,1
     writecube=2 
     goto,finishup	 
  ENDIF
  IF ~(sxpar(header,'CUNIT3')) then begin
     IF channelwidth GT 100. then begin
        veltype='M/S'
        sxaddpar,header,'CUNIT3','M/S',after='CDELT3'
        writecube=1
     endif else begin 
        veltype='KM/S'
        sxaddpar,header,'CUNIT3','KM/S',after='CDELT3'
        writecube=1
     endelse
     IF size(log,/TYPE) EQ 7 then begin
        openu,66,log,/APPEND
        printf,66,linenumber()+'CLEAN_HEADER: Your header did not have a unit for the third axis, that is bad policy.'          
        printf,66,linenumber()+'CLEAN_HEADER: We have set it to '+veltype+'. Please ensure that is correct.'          
        close,66
     ENDIF ELSE BEGIN
        print,linenumber()+'CLEAN_HEADER: Your header did not have a unit for the third axis, that is bad policy.'          
        print,linenumber()+'CLEAN_HEADER: We have set it to '+veltype+'. Please ensure that is correct.'     
     ENDELSE
  ENDIF
 
  IF STRUPCASE(strtrim(velproj,2)) NE 'VELO-HEL' AND $
     STRUPCASE(strtrim(velproj,2)) NE 'VELO-LSR' AND $
     STRUPCASE(strtrim(velproj,2)) NE 'FELO-HEL' AND $
     STRUPCASE(strtrim(velproj,2)) NE 'FELO-LSR' AND $
     STRUPCASE(strtrim(velproj,2)) NE 'VELO' AND $
     STRUPCASE(strtrim(velproj,2)) NE 'FREQ' then begin
     checkax=strtrim(strsplit(velproj,'-',/extract),2)
     IF STRUPCASE(checkax[0]) EQ 'DEC' or STRUPCASE(checkax[0]) EQ 'RA' then begin
        IF size(log,/TYPE) EQ 7 then begin
           openu,66,log,/APPEND
           printf,66,linenumber()+'CLEAN_HEADER: Your zaxis is a spatial axis not a velocity axis.'
           printf,66,linenumber()+'CLEAN_HEADER: Please arrange your cube logically'
           close,66
        ENDIF ELSE BEGIN
           print,linenumber()+'CLEAN_HEADER: Your zaxis is a spatial axis not a velocity axis.'
           print,linenumber()+'CLEAN_HEADER: Please arrange your cube logically'
        ENDELSE
        openu,1,outputcatalogue,/APPEND
        printf,1,format='(A60,2A12,A120)',Dir,0.,0.,'The Cube is not arranged properly'
        close,1
        writecube=2
        goto,finishup
     ENDIF   
     sxaddpar,header,'CTYPE3','VELO',after='CUNIT3'
     writecube=1
     IF size(log,/TYPE) EQ 7 then begin
        openu,66,log,/APPEND
        printf,66,linenumber()+'CLEAN_HEADER: Your velocity projection is not standard. The keyword is changed to VELO (relativistic definition). This might be dangerous.'          
        close,66
     ENDIF ELSE BEGIN
        print,linenumber()+'CLEAN_HEADER: Your velocity projection is not standard. The keyword is changed to VELO (relativistic definition). This might be dangerous.'   
     ENDELSE
  ENDIF
  IF STRUPCASE(veltype) EQ 'KM/S' then begin
     IF size(log,/TYPE) EQ 7 then begin
        openu,66,log,/APPEND
        printf,66,linenumber()+'CLEAN_HEADER: The channels in your input cube are in km/s. This sometimes leads to problems with wcs lib, hence we change it to m/s.'          
        close,66
     ENDIF ELSE BEGIN
        printf,66,linenumber()+'CLEAN_HEADER: The channels in your input cube are in km/s. This sometimes leads to problems with wcs lib, hence we change it to m/s.'   
     ENDELSE
     sxaddpar,header,'CDELT3',sxpar(header,'CDELT3')*1000.
     sxaddpar,header,'CRVAL3',sxpar(header,'CRVAL3')*1000.
     sxaddpar,header,'CUNIT3','M/S'
  ENDIF
                                ;Let's check for presence of the beam in the header. IF present
                                ;supersede the input. !!!!Be careful with smoothed data. If not let's add
                                ;it from the file; This is important if you do it incorrectly the
                                ;pipeline will use a incorrect pixel size.
                                ;   writecube=0

  IF ~(sxpar(header,'BMAJ')) then begin
     IF sxpar(header,'BMMAJ') then begin
        sxaddpar,header,'BMAJ',sxpar(header,'BMMAJ')/3600.
         writecube=1
      endif else begin
         IF n_elements(sxpar(header,'HISTORY')) GT 0 then begin
            history=sxpar(header,'HISTORY')
            found=0.
            for j=n_elements(history)-1,0, -1 do begin
               tmp=strtrim(strsplit(history[j],' ',/extract),2)
               tmp2=WHERE(strupcase(tmp) EQ 'BMAJ=')
               IF tmp2[0] NE -1 then begin
                  sxaddpar,header,'BMAJ',double(tmp[tmp2[0]+1])
                  found =1.
                  writecube=1
               ENDIF
               tmp2=WHERE(strupcase(tmp) EQ 'BMIN=')
               IF tmp2[0] NE -1 then begin
                  sxaddpar,header,'BMIN',double(tmp[tmp2[0]+1])
               ENDIF
               tmp2=WHERE(strupcase(tmp) EQ 'BPA=')
               IF tmp2[0] NE -1 then begin
                  sxaddpar,header,'BPA',double(tmp[tmp2[0]+1])
               ENDIF
               IF found then break
            endfor
         ENDIF
         IF ~(found) then begin
            IF beam[0] NE 1 then begin
               sxaddpar,header,'BMAJ',double(beam[0]/3600.)
               If beam[1] NE 1 then  sxaddpar,header,'BMIN',double(beam[1]/3600.) else  sxaddpar,header,'BMIN',double(beam[0]/3600.)
               writecube=1
            ENDIF ELSE BEGIN
               IF size(log,/TYPE) EQ 7 then begin
                  openu,66,log,/APPEND
                  printf,66,linenumber()+'CLEAN_HEADER: WE CANNOT FIND THE MAJOR AXIS FWHM IN THE HEADER'          
                  close,66
               ENDIF ELSE BEGIN
                  print,linenumber()+'CLEAN_HEADER: WE CANNOT FIND THE MAJOR AXIS FWHM IN THE HEADER'   
               ENDELSE
               openu,1,outputcatalogue,/APPEND
               printf,1,format='(A60,2A12,A120)',Dir,0.,0.,'The Cube has no major axis FWHM in the header.'
               close,1
               writecube=2
               goto,finishup
            ENDELSE
         ENDIF
     ENDELSE
    
  ENDIF 
  IF ~(sxpar(header,'BMIN')) then begin
     IF sxpar(header,'BMMIN') then begin
        sxaddpar,header,'BMIN',sxpar(header,'BMMIN')/3600.
        writecube=1
     endif else begin
        IF sxpar(header,'BMAJ') then begin
           sxaddpar,header,'BMIN', sxpar(header,'BMAJ')
           writecube=1
           IF size(log,/TYPE) EQ 7 then begin
              openu,66,log,/APPEND
              printf,66,linenumber()+'CLEAN_HEADER: We cannot find the minor axis FWHM. Assuming a circular beam.'          
              close,66
           ENDIF ELSE BEGIN
              print,linenumber()+'CLEAN_HEADER: We cannot find the minor axis FWHM. Assuming a circular beam.'
           ENDELSE
        ENDIF ELSE BEGIN
           IF size(log,/TYPE) EQ 7 then begin
              openu,66,log,/APPEND
              printf,66,linenumber()+'CLEAN_HEADER: WE CANNOT FIND THE MINOR AXIS FWHM IN THE HEADER'          
              close,66
           ENDIF ELSE BEGIN
              print,linenumber()+'CLEAN_HEADER: WE CANNOT FIND THE MINOR AXIS FWHM IN THE HEADER'   
           ENDELSE
           
        ENDELSE
     ENDELSE
  ENDIF
 
  IF n_elements(sxpar(header,'HISTORY')) GT 10. then begin
     sxdelpar,header,'HISTORY'
     writecube=1
     IF size(log,/TYPE) EQ 7 then begin
        openu,66,log,/APPEND
        printf,66,linenumber()+'CLEAN_HEADER: Your cube has a significant history attached we are removing it for easier interpretation.'          
        close,66
     ENDIF ELSE BEGIN
        print,linenumber()+'CLEAN_HEADER:  Your cube has a significant history attached we are removing it for easier interpretation.'    
     ENDELSE
  ENDIF
  IF ABS(sxpar(header,'BMAJ')/sxpar(header,'CDELT1')) LT 2. then begin
     IF size(log,/TYPE) EQ 7 then begin
        openu,66,log,/APPEND
        printf,66,linenumber()+'CLEAN_HEADER: !!!!!!!!!!Your cube has less than two pixels per beam major axis.!!!!!!!!!!!!!!!!!'
        printf,66,linenumber()+'CLEAN_HEADER: !!!!!!!!!!           This will lead to bad results.              !!!!!!!!!!!!!!!!'     
        close,66
     ENDIF
     print,linenumber()+'CLEAN_HEADER: !!!!!!!!!!Your cube has less than two pixels per beam major axis.!!!!!!!!!!!!!!!!!'
     print,linenumber()+'CLEAN_HEADER: !!!!!!!!!!           This will lead to bad results.              !!!!!!!!!!!!!!!!'     
  ENDIF
  IF ABS(sxpar(header,'BMAJ')/sxpar(header,'CDELT1')) GT (sxpar(header,'NAXIS1')+sxpar(header,'NAXIS1'))/2. then begin
     IF size(log,/TYPE) EQ 7 then begin
        openu,66,log,/APPEND
        printf,66,linenumber()+'CLEAN_HEADER: !!!!!!!!!!Your cube is smaller than the beam major axis. !!!!!!!!!!!!!!!!!'
        printf,66,linenumber()+'CLEAN_HEADER: !!!!!!!!!!         This will not work.          !!!!!!!!!!!!!!!!'     
        close,66
     ENDIF
     print,linenumber()+'CLEAN_HEADER: !!!!!!!!!!Your cube is smaller than the beam major axis. !!!!!!!!!!!!!!!!!'
     print,linenumber()+'CLEAN_HEADER: !!!!!!!!!!         This will not work.          !!!!!!!!!!!!!!!!'
     openu,1,outputcatalogue,/APPEND
     printf,1,format='(A60,2A12,A120)',Dir,0.,0.,'The Cube is not arranged properly'
     close,1
     writecube=2
     goto,finishup
  ENDIF   

  
  beam=[sxpar(header,'BMAJ')*3600,sxpar(header,'BMIN')*3600.]
  finishup:
end
