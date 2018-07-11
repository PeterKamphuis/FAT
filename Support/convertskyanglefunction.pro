function convertskyanglefunction,angle,distance,UNIT=unit,DISTANCE_UNIT=unitdistance,PHYSICAL=inv

;+
; NAME:
;       CONVERTSKYANGLEFUNCTION
;
; PURPOSE:
;       Program to convert sky angles in arcsecond arcminute or degree
;       to physical size and vice versa
;
; CATEGORY:
;       Support
; 
; CALLING SEQUENCE:
;       Result = CONVERTSKYANGLEFUNCTION(angle,distance,UNIT=unit,DISTANCE_UNIT=unitdistance,/PHYSICAL)
;
;
; INPUTS:
;       Angle = The sky angle to be converted in arcsec
;       Distance = The distance to the object in Mpc
;
; OPTIONAL INPUTS:
;       - 
;
; KEYWORD PARAMETERS:
;       /PHYSICAL - if set it is assumed that angle is in kpc and has
;                   to be converted to arcsec
;       UNIT = unit of the angle, arcsec, arcmin and degree are recognized
;       DISTANCE_UNIT = unit of the distance, Mpc, Kpc and pc are
;       recognized
;
; OUTPUTS:
;       Result = the transformed angle
;
; OPTIONAL OUTPUTS:
; 
; PROCEDURES CALLED:
;       STRTRIM(),STRCOMPRESS(),STRING(),STR_SEP(),FLOOR()
;
; EXAMPLE:
;      
;
; MODIFICATION HISTORY:
;       11-12-2009 added a error message displaying usage (P. Kamphuis)
;       Written by P.Kamphuis 03-09-2009 
;
; NOTE:
;     
;-
  COMPILE_OPT IDL2 
  CATCH,Error_status
  IF  Error_status NE 0. THEN BEGIN
     print, 'CONVERTSKYANGLEFUNCTION: Oops the following went wrong:'
     print, !ERROR_STATE.MSG
     print, 'CONVERTSKYANGLEFUNCTION: Use convertskyangle in this way:'
     print,'CONVERTSKYANGLEFUNCTION: CALLING SEQUENCE: convertskyangle,angle,distance'
     print,"CONVERTSKYANGLEFUNCTION: UNIT = unit of the angle ('arcsec')"
     print,"CONVERTSKYANGLEFUNCTION: DISTANCE_UNIT = unit of the distance ('Mpc')"
     print,"CONVERTSKYANGLEFUNCTION: !!!! Don't forget the apostrophes when using"
     print,"CONVERTSKYANGLEFUNCTION: /PHYSICAL : Keyword to go from kpc to arcsec default unit becomes kpc"  
     goto,ending
  endif

  if n_elements(unit) EQ 0 then begin
     if not keyword_set(inv) then unit='arcsec' else unit='kpc'
  endif

  if n_elements(unitdistance) EQ 0 then unitdistance='Mpc'
  case unitdistance of
     'Mpc': begin
        Dkpc=distance*10E2
     end
     'kpc': begin
        Dkpc=distance
     end
     'pc': begin
        Dkpc=distance/1E3
     end
     'mpc': begin
        Dkpc=distance*10E2
     end
     'Kpc': begin
        Dkpc=distance
     end
     'Pc': begin
        Dkpc=distance/1E3
     end
     else: begin
        print,'CONVERTSKYANGLEFUNCTION: '+unitdistance+' is an unknown unit to convertskyangle'
        print,'CONVERTSKYANGLEFUNCTION: please use Mpc, kpc or pc'
        goto,ending
     end
  endcase
  if not keyword_set(inv) then begin
     case unit of
        'arcsec': begin
           radians=dblarr(n_elements(angle))
           radians=(angle/3600.)*((2.*!pi)/360.)      
        end
        'arcmin': begin
           radians=dblarr(n_elements(angle))
           radians=(angle/60.)*((2.*!pi)/360.) 
        end
        'degree': begin
           radians=dblarr(n_elements(angle))
           radians=(angle)*((2.*!pi)/360.) 
        end
        'Arcsec': begin
           radians=dblarr(n_elements(angle))
           radians=(angle/3600.)*((2.*!pi)/360.) 
        end
        'Arcmin': begin
           radians=dblarr(n_elements(angle))
           radians=(angle/60.)*((2.*!pi)/360.) 
        end
        'Degree': begin
           radians=dblarr(n_elements(angle))
           radians=(angle)*((2.*!pi)/360.) 
        end
        else: begin
           print,'CONVERTSKYANGLEFUNCTION: '+unit+' is an unknown unit to convertskyangle'
           print,'CONVERTSKYANGLEFUNCTION: please use arcsec, arcmin or degree'
           goto,ending
        end
        
     endcase
     kpc=2.*(Dkpc*TAN(radians/2.))
     out=kpc
  endif

  if keyword_set(inv) then begin
     case unit of
        'kpc': begin
           kpc=dblarr(n_elements(angle))
           kpc=angle
        end
        'mpc': begin
           kpc=dblarr(n_elements(angle))
           kpc=angle*1E3
        end
        'pc': begin
           kpc=dblarr(n_elements(angle))
           kpc=angle/1E3
        end
        'Kpc': begin
           kpc=dblarr(n_elements(angle))
           kpc=angle
        end
        'Mpc': begin
           kpc=dblarr(n_elements(angle))
           kpc=angle*1E3
        end
        'Pc': begin
           kpc=dblarr(n_elements(angle))
           kpc=angle/1E3
        end
        else: begin
           print,'CONVERTSKYANGLEFUNCTION: '+unit+' is an unknown unit to convertskyangle'
           print,'CONVERTSKYANGLEFUNCTION: please use kpc, Mpc or pc'
           goto,ending
        end
        
     endcase
     angsec=dblarr(n_elements(kpc))
     radians=dblarr(n_elements(kpc))
     radians[*]=2.*ATAN((kpc[*])/(2.*Dkpc))
     angsec[*]=(radians[*]*(360./(2.*!pi)))*3600.

     out=angsec
  endif




  
  return,out
ending:
  return,-999999
end
