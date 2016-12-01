Pro set_sbr,SBRinput1,SBRinput2,SBRinput3,SBRinput4,SBRinput5,SBRinput6,SBRarr,cutoff,norings,finishafter,log=log,initial=initial,doubled=doubled

;+
; NAME:
;       SET_SBR
;
; PURPOSE:
;       Routine to update the sbr fitting parameters based on the cutoff and the ring numbers
;
; CATEGORY:
;       Support
; 
; CALLING SEQUENCE:
;       SET_SBR,SBRinput1,SBRinput2,SBRinput3,SBRinput4,SBRinput5,SBRinput6,SBRarr,cutoff,norings,finishafter,LOG=log,/INITIAL
;
;
; INPUTS:
;       SBRinput1 = template fitting setting
;       cutoff = cutoff values for the nodes
;       norings = number of rings in model
;       finishafter = indicator of type of fit
;       
;
; OPTIONAL INPUTS:
;       LOG = name of the tracing log 
;
; KEYWORD PARAMETERS:
;       /INITIAL - indicates it is the first time parameters are being set
;
; OUTPUTS:
;       SBRinput1-6= Fitting parameters for different ranges of the model
;
; OPTIONAL OUTPUTS:
;       -
; 
; PROCEDURES CALLED:
;       STR_SEP(), STRTRIM(), STRCOMPRESS(), STRING()
;
; EXAMPLE:
;      
;
; MODIFICATION HISTORY:
;       Written 01-01-2015 P.Kamphuis v1.0
;
; NOTE:
;     
;-
  COMPILE_OPT IDL2 

                           
                                ;first define some stepsizes based on
                                ;the cutoff values
  minstep=strtrim(strcompress(string(cutoff[norings[0]]/2.,format='(E8.1)')),1)
  startstep=strtrim(strcompress(string(5.*cutoff[norings[0]-1],format='(E8.1)')),1)
  if keyword_set(initial) then satlevel=strtrim(strcompress(string(5.*double(startstep),format='(E8.1)')),1) else satlevel=strtrim(strcompress(string(2.*double(startstep),format='(E8.1)')),1)
  SBRinput1=['!SBR '+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+':3','1',strtrim(strcompress(string(cutoff[n_elements(cutoff)-1],format='(E12.5)')),1),startstep,minstep,satlevel,minstep,'3','70','70']
  IF doubled then sbrinput1[2]=strtrim(strcompress(string(cutoff[norings[0]]/1.2,format='(E8.1)')),1)
                                ;Then copy the array into multiple arrays
  SBRinput2=SBRinput1
  SBRinput3=SBRinput1
  SBRinput4=SBRinput1
                                ;and adjust them properly
  SBRinput2[0]='!SBR_2 '+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+':3'
  SBRinput3[0]='SBR  1 SBR_2 1 '
  SBRinput3[1]=strtrim(strcompress(string(SBRarr[1],format='(E12.5)')),1)
  SBRinput4[0]='SBR 2 SBR_2 2 '
                                ;If the disk is too small
  IF norings[0] LE 4 OR finishafter EQ 1.1  then begin
     SBRinput1[2]=strtrim(strcompress(string(cutoff[n_elements(cutoff)-1]/4.,format='(E12.5)')),1)
     SBRinput2[2]=strtrim(strcompress(string(cutoff[n_elements(cutoff)-1]/4.,format='(E12.5)')),1)
     SBRinput3[2]=strtrim(strcompress(string(cutoff[n_elements(cutoff)-1],format='(E12.5)')),1)
     SBRinput4[2]=strtrim(strcompress(string(cutoff[n_elements(cutoff)-1]*2.,format='(E12.5)')),1)
  ENDIF
                                ;IF the disk is large than treat inner and outer parts differently
  IF norings[0] GT 6 and not keyword_set(initial) then begin
     SBRinput5=SBRinput1
     SBRinput6=SBRinput2
     halfrings=floor((norings[0]-2)/2.)
     IF halfrings LE 4 then begin
        SBRinput1[0]='!SBR '+strtrim(strcompress(string(norings[0]-2,format='(F7.4)')),1)+':3'
        SBRinput1[2]=strtrim(strcompress(string(cutoff[norings[0]-1],format='(E12.5)')),1)
        SBRinput5[0]='!SBR '+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+':'+strtrim(strcompress(string(norings[0]-1,format='(F7.4)')),1)
        SBRinput5[2]=strtrim(strcompress(string(cutoff[norings[0]-1]/4.,format='(E12.5)')),1)
        SBRinput2[0]='!SBR_2 '+strtrim(strcompress(string(norings[0]-2,format='(F7.4)')),1)+':3'
        SBRinput2[2]=strtrim(strcompress(string(cutoff[norings[0]-1],format='(E12.5)')),1)
        SBRinput6[0]='!SBR_2 '+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+':'+strtrim(strcompress(string(norings[0]-1,format='(F7.4)')),1)
        SBRinput6[2]=strtrim(strcompress(string(cutoff[norings[0]-1]/4.,format='(E12.5)')),1)
     ENDIF ELSE BEGIN
        SBRinput1[0]='!SBR '+strtrim(strcompress(string(halfrings+2,format='(F7.4)')),1)+':3'
        SBRinput1[2]=strtrim(strcompress(string(cutoff[norings[0]-1],format='(E12.5)')),1)
        SBRinput5[0]='!SBR '+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+':'+strtrim(strcompress(string(halfrings+3,format='(F7.4)')),1)
        SBRinput5[2]=strtrim(strcompress(string(cutoff[norings[0]-1]/4.,format='(E12.5)')),1)
        SBRinput2[0]='!SBR_2 '+strtrim(strcompress(string(halfrings+2,format='(F7.4)')),1)+':3'
        SBRinput2[2]=strtrim(strcompress(string(cutoff[norings[0]-1],format='(E12.5)')),1)
        SBRinput6[0]='!SBR_2 '+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+':'+strtrim(strcompress(string(halfrings+3,format='(F7.4)')),1)
        SBRinput6[2]=strtrim(strcompress(string(cutoff[norings[0]-1]/4.,format='(E12.5)')),1)
     ENDELSE
  ENDIF
end
