Pro set_sdis,sdisinput1,SDISarr,velconstused,sdismax,sdismin,norings,channelwidth,AVINNER=avinner,START=start,CENTRALEXCLUDE=centralexclude,FINISH_AFTER=finish_after,slope=slope,debug=debug

;+
; NAME:
;       SET_SDIS
;
; PURPOSE:
;       Routine to update the SDIS fitting parameters, this is a
;       modified copy of set_sdisv6
;
; CATEGORY:
;       Support
; 
; CALLING SEQUENCE:
;       SET_SDISV6, sdisinput1, SDISarr, velconstused, sdismax,
;       sdismin, norings, channelwidth, AVINNER=avinner, START=start,
;       CENTRALEXCLUDE=centralexclude, FINISH_AFTER=finish_after
;
;
; INPUTS:
;       sdisinput1 = template fitting setting
;       SDISarr = currently fitted sdis
;       velconstused = number of rings that are fitted as slope at the
;       outer radii
;       sdismax = maximum allowed value
;       sdismin = minimum allowed value
;       norings = number of rings in the model
;       channelwidth = the width of a channel in the cube
;
; OPTIONAL INPUTS:
;       FINISH_AFTER = indicator of type of fit
;       START = indicator from which node fitting should start
;       CENTRALEXCLUDE = trigger to exclude the inner ring from
;       fitting
;       SLOPE = if 0 then fit slope if 1 then fit flat.  
;
; KEYWORD PARAMETERS:
;       /INITIAL - indicates it is the first time parameters are being set
;
; OUTPUTS:
;       sdisinput1= Updated fitting parameters for SDIS
;
; OPTIONAL OUTPUTS:
;       AVINNER =  the average value of the outer part of the rotation curve
; 
; PROCEDURES CALLED:
;       STR_SEP(), STRTRIM(), STRCOMPRESS(), STRING()
;
; EXAMPLE:
;      
;
; MODIFICATION HISTORY:
;       14-09-2018 P.Kamphuis; Start of set_sdis  
;       10-03-2016 P.Kamphuis; Added an option to not fit a slope but
;                              fit a flat extension of the last reliable ring.  
;       Written 01-01-2015 P.Kamphuis v1.0
;
; NOTE:
;     
;-
  COMPILE_OPT IDL2 
  
                                ;If the model is big we always want to
                                ;fit a slope to the outer parts
  IF n_elements(slope) EQ 0 then slope=0
 
  IF double(norings[0]) GT 15. then begin
     IF n_elements(finish_after) EQ 0 then finish_after=2.
     IF finish_after EQ 2.1 then begin
        IF velconstused GT norings[0]-ceil(norings[0]/10.) then velconstused=norings[0]-ceil(norings[0]/10.)
     endif else begin
        IF velconstused GT norings[0]-ceil(norings[0]/5.) then velconstused=norings[0]-ceil(norings[0]/5.)
     Endelse
  ENDIF
  if n_elements(start) EQ 0 then begin
     start=floor(norings[0]/4.)
  endif
  IF start LT 4 then start=4
  IF n_elements(centralexclude) EQ 0 then centralexclude=0 
                                ;If we want to fit everything as slope
                                ;then just leave the most inner ring free
  IF norings[0]-velconstused LT 2 then velconstused=norings[0]-1
  IF velconstused LT 5 then velconstused=5
  IF velconstused LT start+1 then velconstused=start+1
  avinner=0.
 
                                ;set the fitting parameters for the rotation curve
  case (1) of
                                ;If a small galaxy we fit the
                                ;dispersion as 1
     norings[0] LE 4 OR finish_after EQ 1.1:begin
  ;      SDISinput1=['SDIS 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
  ;               ' SDIS_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1),$
  ;               '25','5','1','0.1','0.5','0.05','3','70','70']    
        string1='SDIS '+strtrim(strcompress(string(norings[0],format='(I3)')),1)+':'+strtrim(strcompress(string(1,format='(I3)')),1)+' SDIS_2 '+strtrim(strcompress(string(norings[0],format='(I3)')),1)+':'+strtrim(strcompress(string(1,format='(I3)')),1)
        string2=string(SDISmax)
        string3=string(SDISmin)
        string4=string(0.5*channelwidth)
        string5=string(0.01*channelwidth)
        string6=string(0.5*channelwidth)
        string7=string(0.01*channelwidth)
        string8='3'
        string9='70'
        string10=' '
     end
                                ; IF only the second to last ring is sloped
     velconstused GE norings[0]-1:begin
        avinner=SDISarr[n_elements(SDISarr)-2]
        case (1) of
           avinner LT SDISmin:begin
              tmp=where(SDISarr LT sdismin)
              IF tmp[0] NE -1 then SDISarr[tmp]=sdismin*1.2
              avinner=sdismin*1.2
           end
           avinner LT channelwidth: avinner=channelwidth 
           else:begin
           end
        endcase
        string1='SDIS '+strtrim(strcompress(string(norings[0],format='(I3)')),1)+$
                ':'+strtrim(strcompress(string(norings[0]-2,format='(I3)')),1)+$
                ' SDIS_2 '+strtrim(strcompress(string(norings[0],format='(I3)')),1)+$
                ':'+strtrim(strcompress(string(norings[0]-2,format='(I3)')),1)+$
                ',!SDIS '+strtrim(strcompress(string(norings[0]-3,format='(I3)')),1)+':'$
                +strtrim(strcompress(string(start+1,format='(I3)')),1)+' SDIS_2 '+$
                strtrim(strcompress(string(norings[0]-3,format='(I3)')),1)+$
                ':'+strtrim(strcompress(string(start+1,format='(I3)')),1)+', SDIS 1:'+strtrim(strcompress(string(start,format='(I3)')),1)+' SDIS_2 1:'+strtrim(strcompress(string(start,format='(I3)')),1)
        string2=string(SDISmax)+' '+string(SDISmax)+' '+string(SDISmax)
        string3=string(SDISmin)+' '+string(SDISmin)+' '+string(avinner)
        string4=string(1.5*channelwidth)+' '+string(1.5*channelwidth)+' '+string(1.5*channelwidth)
        string5=string(0.1*channelwidth)+' '+string(0.1*channelwidth)+' '+string(0.1*channelwidth)
        string6=string(0.75*channelwidth)+' '+string(0.75*channelwidth)+' '+string(0.75*channelwidth)
        string7=string(0.1*channelwidth)+' '+string(0.01*channelwidth)+' '+string(0.1*channelwidth)
        string8='3 3 3'
        string9='70 70 70'
;        string10=' SDIS '+strtrim(strcompress(string(norings[0]-1,format='(I3)')),1)+' '+strtrim(strcompress(string(norings[0],format='(I3)')),1)+' SDIS_2 '+strtrim(strcompress(string(norings[0]-1,format='(I3)')),1)+' '+strtrim(strcompress(string(norings[0],format='(I3)')),1)
        string10=' SDIS '+strtrim(strcompress(string(norings[0]-1,format='(I3)')),1)+' SDIS_2 '+strtrim(strcompress(string(norings[0]-1,format='(I3)')),1)
     End
                                ;All other cases
     else:begin
        avinner=SDISarr[velconstused]
        case (1) of
           avinner LT SDISmin:begin
              tmp=where(SDISarr LT sdismin)
              IF tmp[0] NE -1 then SDISarr[tmp]=sdismin*1.2
              avinner=sdismin*1.2
           end
           avinner LT channelwidth: avinner=channelwidth
           avinner LT SDISarr[n_elements(SDISarr)-1]: avinner= SDISarr[n_elements(SDISarr)-1]*1.1
           else:begin
           end
        endcase
     
        string1='SDIS '+strtrim(strcompress(string(norings[0],format='(I3)')),1)+':'+$
                strtrim(strcompress(string(velconstused,format='(I3)')),1)+' SDIS_2 '$
                +strtrim(strcompress(string(norings[0],format='(I3)')),1)+':'$
                +strtrim(strcompress(string(velconstused,format='(I3)')),1)+$
                ', !SDIS '+strtrim(strcompress(string(velconstused-1,format='(I3)')),1)+$
                ':'+strtrim(strcompress(string(start+1,format='(I3)')),1)+' SDIS_2 '$
                +strtrim(strcompress(string(velconstused-1,format='(I3)')),1)+$
                ':'+strtrim(strcompress(string(start+1,format='(I3)')),1)+', SDIS 1:'+strtrim(strcompress(string(start,format='(I3)')),1)+' SDIS_2 1:'+strtrim(strcompress(string(start,format='(I3)')),1)
        string2=string(SDISmax)+' '+string(SDISmax)+' '+string(SDISmax)
        string3=string(avinner)+' '+string(SDISmin)+' '+string(SDISmin)
        string4=string(1.5*channelwidth)+' '+string(1.5*channelwidth)+' '+string(1.5*channelwidth)
        string5=string(0.1*channelwidth)+' '+string(0.1*channelwidth)+' '+string(0.1*channelwidth)
        string6=string(0.75*channelwidth)+' '+string(0.75*channelwidth)+' '+string(0.75*channelwidth)
        string7=string(0.1*channelwidth)+' '+string(0.01*channelwidth)+' '+string(0.1*channelwidth)
        string8='3 3 3'
        string9='70 70 70'
        string10=' SDIS '+strtrim(strcompress(string(norings[0]-1,format='(I3)')),1)+$
                 ':'+strtrim(strcompress(string(velconstused+1,format='(I3)')),1)+$
                 ' SDIS_2 '+strtrim(strcompress(string(norings[0]-1,format='(I3)')),1)+$
                 ':'+strtrim(strcompress(string(velconstused+1,format='(I3)')),1)
     end
  endcase
                                ;make an array
  SDISinput1=[string1,string2,string3,string4,string5,string6,string7,string8,string9,string9]   
end
