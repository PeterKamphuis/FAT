Pro set_vrotv6,vrotinput1,VROTarr,velconstused,vrotmax,vrotmin,norings,channelwidth,AVINNER=avinner,START=start,CENTRALEXCLUDE=centralexclude,FINISH_AFTER=finish_after,slope=slope

;+
; NAME:
;       SET_VROTV6
;
; PURPOSE:
;       Routine to update the VROT fitting parameters
;
; CATEGORY:
;       Support
; 
; CALLING SEQUENCE:
;       SET_VROTV6, vrotinput1, VROTarr, velconstused, vrotmax,
;       vrotmin, norings, channelwidth, AVINNER=avinner, START=start,
;       CENTRALEXCLUDE=centralexclude, FINISH_AFTER=finish_after
;
;
; INPUTS:
;       vrotinput1 = template fitting setting
;       VROTarr = currently fitted vrot
;       velconstused = number of rings that are fitted as slope at the
;       outer radii
;       vrotmax = maximum allowed value
;       vrotmin = minimum allowed value
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
;       vrotinput1= Updated fitting parameters for VROT
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
;       14-09-2018 P.Kamphuis; When sloping the outer ring were done
;                              later not as the first ones. Swapped these  
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
  if n_elements(start) EQ 0 then start=2
  IF start LT 1 then start=2
  IF n_elements(centralexclude) EQ 0 then centralexclude=0 
                                ;If we want to fit everything as slope
                                ;then just leave the most inner ring free
  IF norings[0]-velconstused LT 2 then velconstused=norings[0]-1
  IF velconstused LT 5 then velconstused=5
  IF velconstused LT start+1 then velconstused=start+1
  avinner=0.
                                ;set the fitting parameters for the rotation curve
  case (1) of
                                ;If a small galaxy then no slope fitting
     norings[0] LE 4:begin
        string1='!VROT '+strtrim(strcompress(string(norings[0],format='(I3)')),1)+':'+strtrim(strcompress(string(start,format='(I3)')),1)+' VROT_2 '+strtrim(strcompress(string(norings[0],format='(I3)')),1)+':'+strtrim(strcompress(string(start,format='(I3)')),1)
        string2=string(VROTmax)
        string3=string(VROTmin)
        string4=string(channelwidth)
        string5=string(0.01*channelwidth)
        string6=string(0.5*channelwidth)
        string7=string(0.01*channelwidth)
        string8='3'
        string9='70'
        string10=' '
     end
                                ; IF only the second to last ring is sloped
     velconstused GE norings[0]-1:begin
        avinner=VROTarr[n_elements(VROTarr)-2]*0.9
        case (1) of
           avinner GT VROTmax:begin
              tmp=where(VROTarr GT vrotmax)
              IF tmp[0] NE -1 then VROTarr[tmp]=vrotmax*0.8
              avinner=vrotmax*0.8
           end
           avinner GT 150.: avinner=150. 
           else:begin
           end
        endcase
        string1='!VROT '+strtrim(strcompress(string(norings[0],format='(I3)')),1)+':'+$
                strtrim(strcompress(string(norings[0]-2,format='(I3)')),1)+$
                ' VROT_2 '+strtrim(strcompress(string(norings[0],format='(I3)')),1)+':'+$
                strtrim(strcompress(string(norings[0]-2,format='(I3)')),1)+$
                ',!VROT '+strtrim(strcompress(string(norings[0]-3,format='(I3)')),1)+':'$
                +strtrim(strcompress(string(start,format='(I3)')),1)+' VROT_2 '+$
                strtrim(strcompress(string(norings[0]-3,format='(I3)')),1)+$
                ':'+strtrim(strcompress(string(start,format='(I3)')),1)
        string2=string(VROTmax)+' '+string(VROTmax)
        string3=string(avinner)+' '+string(VROTmin)
        string4=string(0.5*channelwidth)+' '+string(channelwidth)
        string5=string(0.01*channelwidth)+' '+string(0.01*channelwidth)
        string6=string(0.5*channelwidth)+' '+string(0.5*channelwidth)
        string7=string(0.01*channelwidth)+' '+string(0.01*channelwidth)
        string8='3 3'
        string9='70 70'
        case (1) of
           centralexclude EQ 1 AND slope EQ 0:begin
              string10=' VROT 2 '+strtrim(strcompress(string(norings[0]-1,format='(I3)')),1)+' VROT_2 2 '+strtrim(strcompress(string(norings[0]-1,format='(I3)')),1)
           end
           centralexclude EQ 1 AND slope EQ 1:begin
              string10=' VROT 2 '+strtrim(strcompress(string(norings[0]-1,format='(I3)')),1)+' '+$
                       strtrim(strcompress(string(norings[0],format='(I3)')),1)+$
                       ' VROT_2 2 '+strtrim(strcompress(string(norings[0]-1,format='(I3)')),1)+' '+$
                       strtrim(strcompress(string(norings[0],format='(I3)')),1)
           end
           
           centralexclude EQ 0 AND slope EQ 1:begin
              string10=' VROT '+strtrim(strcompress(string(norings[0]-1,format='(I3)')),1)+' '+$
                       strtrim(strcompress(string(norings[0],format='(I3)')),1)+$
                       ' VROT_2 '+strtrim(strcompress(string(norings[0]-1,format='(I3)')),1)+' '+$
                       strtrim(strcompress(string(norings[0],format='(I3)')),1)

           end
           else:begin
              string10=' VROT '+strtrim(strcompress(string(norings[0]-1,format='(I3)')),1)+' VROT_2 '+strtrim(strcompress(string(norings[0]-1,format='(I3)')),1)
           end
        ENDCASE
     End
                                ;All other cases
     else:begin
        avinner=VROTarr[velconstused]*0.8
        case (1) of
           avinner GT VROTmax:begin
              tmp=where(VROTarr GT vrotmax)
              IF tmp[0] NE -1 then VROTarr[tmp]=vrotmax*0.8
              avinner=vrotmax*0.8
           end
           avinner GT 150.: avinner=150. 
           avinner GT VROTarr[n_elements(VROTarr)-1]: avinner= VROTarr[n_elements(VROTarr)-1]*0.9
           else:begin
           end
        endcase
        string1='!VROT '+strtrim(strcompress(string(norings[0],format='(I3)')),1)+':'+$
                strtrim(strcompress(string(velconstused,format='(I3)')),1)+' VROT_2 '+$
                strtrim(strcompress(string(norings[0],format='(I3)')),1)+':'$
                +strtrim(strcompress(string(velconstused,format='(I3)')),1)+', !VROT '$
                +strtrim(strcompress(string(velconstused-1,format='(I3)')),1)+':'$
                +strtrim(strcompress(string(start,format='(I3)')),1)+' VROT_2 '+$
                strtrim(strcompress(string(velconstused-1,format='(I3)')),1)+$
                ':'+strtrim(strcompress(string(start,format='(I3)')),1)
        string2=string(VROTmax)+' '+string(VROTmax)
        string3=string(VROTmin)+' '+string(avinner)
        string4=string(channelwidth)+' '+string(0.1*channelwidth)
        string5=string(0.01*channelwidth)+' '+string(0.01*channelwidth)
        string6=string(0.5*channelwidth)+' '+string(channelwidth)
        string7=string(0.01*channelwidth)+' '+string(0.01*channelwidth)
        string8='3 3'
        string9='70 70'


        case (1) of
           centralexclude EQ 1 AND slope EQ 0:begin
              string10=' VROT 2 '+strtrim(strcompress(string(norings[0]-1,format='(I3)')),1)+$
                       ':'+strtrim(strcompress(string(velconstused+1,format='(I3)')),1)+$
                       ' VROT_2 2 '+strtrim(strcompress(string(norings[0]-1,format='(I3)')),1)+$
                       ':'+strtrim(strcompress(string(velconstused+1,format='(I3)')),1)
           end
           centralexclude EQ 1 AND slope EQ 1:begin
              string10=' VROT 2 '+strtrim(strcompress(string(norings[0],format='(I3)')),1)+$
                       ':'+strtrim(strcompress(string(velconstused+1,format='(I3)')),1)+$
                       ' VROT_2 2 '+strtrim(strcompress(string(norings[0],format='(I3)')),1)+$
                       ':'+strtrim(strcompress(string(velconstused+1,format='(I3)')),1)
           end
           
           centralexclude EQ 0 AND slope EQ 1:begin
              string10=' VROT '+strtrim(strcompress(string(norings[0],format='(I3)')),1)+$
                       ':'+strtrim(strcompress(string(velconstused+1,format='(I3)')),1)+$
                       ' VROT_2 '+strtrim(strcompress(string(norings[0],format='(I3)')),1)+$
                       ':'+strtrim(strcompress(string(velconstused+1,format='(I3)')),1)

           end
           else:begin
              string10=' VROT '+strtrim(strcompress(string(norings[0]-1,format='(I3)')),1)+$
                       ':'+strtrim(strcompress(string(velconstused+1,format='(I3)')),1)+$
                       ' VROT_2 '+strtrim(strcompress(string(norings[0]-1,format='(I3)')),1)+$
                       ':'+strtrim(strcompress(string(velconstused+1,format='(I3)')),1)
           end
        ENDCASE

     end
  endcase
                                ;make an array
  VROTinput1=[string1,string2,string3,string4,string5,string6,string7,string8,string9,string9,string10]   
end
