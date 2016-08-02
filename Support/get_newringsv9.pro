Pro get_newringsv9,SBR1in,SBR2in,cutoffin,newrings,INDIVIDUAL=individual,debug=debug

;+
; NAME:
;       GET_NEWRINGSV8
;
; PURPOSE:
;       This routine compares the SBR profiles and decides what the new number of ring should be
;
; CATEGORY:
;       Support
; 
; CALLING SEQUENCE:
;       GET_NEWRINGSV9,SBR1in,SBR2in,cutoffin,newrings,individual=individual
;
;
; INPUTS:
;       SBR1in = The sbr profile of the approaching side
;       SBR2in = The sbr profile of the receding side
;       cutoffin = Array with the cutoff values
;       
; OPTIONAL INPUTS:
;       - 
;
; KEYWORD PARAMETERS:
;       /INDIVIDUAL - Set this keyword to get an independent ring for
;                     each sides.
;       /DEBUG      - Set this keyword to get printed output during
;                     the running  
;
; OUTPUTS:
;       newrings = the new amount of rings. a 2D array when
;       /INDIVIDUAL is set
;
; OPTIONAL OUTPUTS:
;       -
; 
; PROCEDURES CALLED:
;       MAX(), FLOOR()
;
; EXAMPLE:
;      
;
; MODIFICATION HISTORY:
;       30-05-2016 P.Kamphuis; Added debug keyword for testing  
;       17-05-2016 P.Kamphuis; Added a condition that when the last
;       ring is too low but the second too last is 5* cutoff or
;       bigger then no cutting.  
;       12-03-2016 P.Kamphuis; Removed the +1 in the multiple
;                              detection of toolow.
;       Written by P.Kamphuis 01-01-2015 
;
; NOTE:
;     
;-
  COMPILE_OPT IDL2
  SBR1=SBR1in
  SBR2=SBR2in
  cutoff=cutoffin


  newrings=n_elements(SBR1)
  toolow=WHERE(SBR1 LT cutoff)
  toolow2=WHERE(SBR2 LT cutoff)
  If keyword_set(individual) then begin
     if keyword_set(debug) then print,'Running in individual mode'
     newrings=dblarr(2)
     newrings=[n_elements(SBR1),n_elements(SBR2)]
     IF n_elements(toolow) GT 1 then begin
        if keyword_set(debug) then print,'Multiple rings below cutoff in first array'
        IF toolow[n_elements(toolow)-1] EQ n_elements(SBR1)-1 OR  toolow[n_elements(toolow)-1] EQ n_elements(SBR1)-2 then begin
           for j=n_elements(toolow)-1,1,-1 do begin
              if toolow[j] - toolow[j-1] NE 1 then break
           endfor
           last1=toolow[j]
        ENDIF ELSE last1=n_elements(SBR1)
     ENDIF ELSE BEGIN
        if keyword_set(debug) then print,'Single ring below cutoff in first array. Ring='+string(toolow)
        case (1) of
                                ; first case where the elements are
                                ; the last elements of the array 
           toolow[0] EQ n_elements(SBR1)-1: begin
                                ;If the secondring is more then 5
                                ;times the cutoff this is acceptable
                                ;and we do not cut
              IF SBR1[toolow[0]-1] GE 5.*cutoff[toolow[0]-1] then last1=toolow[0]+1 else last1=toolow[0]
           end
                                ;The second case deals with the last
                                ;ring being risen above the threshold
                                ;but the second to last is below. This
                                ;does require a check to make sure
                                ;that the last ring does not actually
                                ;warrant growing the model.
           toolow[0] EQ n_elements(SBR1)-2 : begin
              IF (SBR1[toolow[0]+1] LE 5.*cutoff[toolow[0]+1]) then last1= n_elements(SBR1)-2 else last1=n_elements(SBR1)
           end
                                ;if the low rings are in some other
                                ;area of the model we do not want to cut
           else: begin
              last1=n_elements(SBR1)
           end
        endcase
     ENDELSE
     IF n_elements(toolow2) GT 1 then begin
         if keyword_set(debug) then print,'Multiple rings below cutoff in second array.'
        IF toolow2[n_elements(toolow2)-1] EQ n_elements(SBR2)-1 OR  toolow2[n_elements(toolow2)-1] EQ n_elements(SBR2)-2 then begin
           for j=n_elements(toolow2)-1,1,-1 do begin
              if toolow2[j] - toolow2[j-1] NE 1 then break
           endfor
           last2=toolow2[j]
        ENDIF ELSE last2=n_elements(SBR2)
     ENDIF ELSE BEGIN
         if keyword_set(debug) then print,'Single ring below cutoff in second array. Ring='+string(toolow)
        case (1) of
                                ; first case where the elements are
                                ; the last elements of the array 
           toolow2[0] EQ n_elements(SBR2)-1: begin
                                ;If the secondring is more then 5
                                ;times the cutoff this is acceptable
                                ;and we do not cut
              IF SBR2[toolow2[0]-1] GE 5.*cutoff[toolow2[0]-1] then last2=toolow2[0]+1 else last2=toolow2[0]
           end 
                                ;The second case deals with the last
                                ;ring being risen above the threshold
                                ;but the second to last is below. This
                                ;does require a check to make sure
                                ;that the last ring does not actually
                                ;warrant growing the model.
           toolow2[0] EQ n_elements(SBR2)-2 : begin
              IF (SBR2[toolow2[0]+1] LE 5.*cutoff[toolow2[0]+1]) then last2= n_elements(SBR1)-2 else last2=n_elements(SBR2)
           end
                                ;if the low rings are in some other
                                ;area of the model we do not want to cut
           else: begin
              last2=n_elements(SBR2)
           end
        endcase
     ENDELSE
     newrings[*]=[last1,last2]

     IF newrings[0] LT 3 then newrings[0]=3
     IF newrings[1] LT 3 then newrings[1]=3 
                                ; we want make sure that the last
                                ; rings are not very bright as
                                ; we want make sure that the last
                                ; rings are not very bright as
     check=0                    ; we'll just end up adding rings 
     WHILE check EQ 0 do begin
        IF (SBR1[newrings[0]-1] GT 25*cutoff[newrings[0]-1] AND newrings[0] LT n_elements(SBR1))  then newrings[0]=newrings[0]+1. else begin 
           check=1
        ENDELSE
     ENDWHILE 
     check=0                    ; we'll just end up adding rings 
     WHILE check EQ 0 do begin
        IF (SBR2[newrings[1]-1] GT 25*cutoff[newrings[1]-1] AND newrings[1] LT n_elements(SBR2))  then newrings[1]=newrings[1]+1. else begin 
           check=1
        ENDELSE
     ENDWHILE 
     IF newrings[0] GT  n_elements(SBR1) then newrings[0]=n_elements(SBR1)
     IF newrings[1] GT  n_elements(SBR2) then newrings[1]=n_elements(SBR2)
     
     
  ENDIF else begin
     if keyword_set(debug) then print,'Running in non-individual mode'
     IF toolow[0] NE -1 AND toolow2[0] NE -1 then begin
        IF n_elements(toolow) GT 1 AND n_elements(toolow2) GT 1 then begin
           if keyword_set(debug) then begin
              print,'Multiple rings below cutoff in Both arrays'
              print,'Array 1:'+STRJOIN(toolow,' ')
              print,'Array 2:'+STRJOIN(toolow2,' ')
           endif
           IF toolow[n_elements(toolow)-1] EQ n_elements(SBR1)-1 OR  toolow[n_elements(toolow)-1] EQ n_elements(SBR1)-2 then begin
              for j=n_elements(toolow)-1,1,-1 do begin
                 if toolow[j] - toolow[j-1] NE 1 then break
              endfor
              last1=toolow[j]
           ENDIF ELSE last1=n_elements(SBR1)
           IF toolow2[n_elements(toolow2)-1] EQ n_elements(SBR2)-1 OR  toolow2[n_elements(toolow2)-1] EQ n_elements(SBR2)-2 then begin
              for j=n_elements(toolow2)-1,1,-1 do begin
                 if toolow2[j] - toolow2[j-1] NE 1 then break
              endfor
              last2=toolow2[j]
           ENDIF ELSE last2=n_elements(SBR2)
           newrings=MAX([last1,last2])

           IF newrings EQ n_elements(SBR1) then goto,checkadd
           IF newrings LT 3 then newrings=3
           
                                ; we want make sure that the last
                                ; rings are not very bright as
                                ; we'll just end up adding rings 
           
        ENDIF ELSE begin
           if keyword_set(debug) then print,'One of the arrays has only one array below the cutoff'
           IF n_elements(toolow) GT 1 then toolow=toolow[n_elements(toolow)-1]
           IF n_elements(toolow2) GT 1 then toolow2=toolow2[n_elements(toolow2)-1]
            if keyword_set(debug) then begin
              print,'Ring below cutoff in Both arrays'
              print,'Array 1:'+STRJOIN(toolow,' ')
              print,'Array 2:'+STRJOIN(toolow2,' ')
              print,'Array extent = '+string(n_elements(SBR1)-1)
           endif
           
           case (1) of
                                ; first case where the elements are
                                ; the last elements of the array 
              toolow[0] EQ n_elements(SBR1)-1 $
                 AND toolow2[0] EQ n_elements(SBR2)-1 $
                 : begin
                                ;If the secondring is more then 5
                                ;times the cutoff this is acceptable
                                ;and we do not cut
                  IF SBR2[toolow2[0]-1] GE 5.*cutoff[toolow2[0]-1] OR SBR1[toolow[0]-1] GE 5.*cutoff[toolow[0]-1] then newrings=toolow[0]+1 else newrings=toolow[0]                            
              end
                                ;The second case deals with the last
                                ;ring being risen above the threshold
                                ;but the second to last is below. This
                                ;does require a check to make sure
                                ;that the last ring does not actually
                                ;warrant growing the model.
              toolow[0] EQ n_elements(SBR1)-2 AND toolow2[0] EQ n_elements(SBR2)-1 : begin
                 IF (SBR1[toolow[0]+1] LE 5.*cutoff[toolow[0]+1]) then newrings= n_elements(SBR2)-1
              end
                                ;case 3 is the reverse of case 2
              toolow2[0] EQ n_elements(SBR2)-2 AND toolow[0] EQ n_elements(SBR1)-1 : begin
                 IF (SBR2[toolow2[0]+1] LE 5.*cutoff[toolow2[0]+1]) then newrings= n_elements(SBR1)-1
              end
                                ;the final case is where in both case
                                ;the second to last ring is less
              toolow[0] EQ n_elements(SBR1)-2 AND toolow2[0] EQ n_elements(SBR2)-2 : begin
                 IF (SBR2[toolow2[0]+1] LE 5.*cutoff[toolow2[0]+1] AND SBR1[toolow[0]+1] LE 5.*cutoff[toolow[0]+1]) then begin
                    newrings=n_elements(SBR1)-2
                 ENDIF ELSE goto,checkadd

                 
                 
                 
              end
                                ;if the low rings are in some other
                                ;area of the model we do not want to cut
              else: begin
                 goto,checkadd
              end
           endcase
           
        ENDELSE 
                                ; we want make sure that the last
                                ; rings are not very bright as
        check=0                 ; we'll just end up adding rings 
        WHILE check EQ 0 do begin
           IF (SBR1[newrings-1] GT 25*cutoff[newrings-1] OR  SBR2[newrings-1] GT 25*cutoff[newrings-1]) AND newrings LT n_elements(SBR1)  then newrings=newrings+1. else begin 
              check=1
              IF newrings EQ n_elements(SBR1)-1 then goto,checkadd
           ENDELSE
        ENDWHILE 
     ENDIF ELSE begin
        checkadd:
        IF SBR1[n_elements(SBR1)-1] GT cutoff[n_elements(SBR1)-1] then begin
           IF (SBR1[n_elements(SBR1)-2] GT 5.*cutoff[n_elements(SBR1)-2] AND SBR1[n_elements(SBR1)-1] GT 3.*cutoff[n_elements(SBR1)-1]) OR SBR1[n_elements(SBR1)-1] GT 5.*cutoff[n_elements(SBR1)-1]  then newrings=n_elements(SBR1)+1
        ENDIF
        IF SBR2[n_elements(SBR2)-1] GT cutoff[n_elements(SBR2)-1] then begin
           IF (SBR2[n_elements(SBR2)-2] GT 5.*cutoff[n_elements(SBR2)-2] AND SBR2[n_elements(SBR2)-1] GT 3.*cutoff[n_elements(SBR2)-1]) OR SBR2[n_elements(SBR2)-1] GT 5.*cutoff[n_elements(SBR2)-1] then newrings=n_elements(SBR2)+1
        ENDIF
     ENDELSE
  ENDELSE
end





