Pro pre_ran_sofia,names,new_dir,supportdirchecked,pixfwhm,header,errormessage,VSYSpix,RApix,DECpix,Totflux,log=log
;+
; NAME:
;       PRE_RAN_SOFIA
;
; PURPOSE:
;       Set parameters from a pre ran sofia. This assumes that all
;       sources belong to galaxy
;
; CATEGORY:
;       Support
; 
; CALLING SEQUENCE:
;       PRE_RAN_SOFIA,names,new_dir,currentfitcube,catcatalogname,supportdirchecked,pixfwhm,header,errormessage,VSYSpix,RApix,DECpix,Totflux,log=log 
;
;
; INPUTS:
;            names  = nameofmask,nameofcatalog,nameofmom0,nameofmom1
;           new_dir = Current working directory.
;    currentfitcube = Currently used cube.
;     catcatlogname = Name of the catalog with SoFiA output.          
; supportdirchecked = Name of support directory without spaces.
;           pixfwhm = FWHM of major axis beam in pixels.
;            header = header of the input cube.
;      errormessage = 2 elements string array to contain the
;                     errormessage and bookkeeping value in case of failure.
;           VSYSpix = 3-elements array containg the systemic velocity pixel
;                     and its boundaries.
;             RApix = 3-elements array containg the central RA pixel
;                     and its boundaries.
;            DECpix = 3-elements array containg the central DEC pixel
;                     and its boundaries.
;           Totflux = variable containing the total fluyx found in the source.  
;     
; OPTIONAL INPUTS:
;       LOG = name of the tracing log 
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
;       FILE_TEST(),READ_FITS(),READ_TEMPLATE,SPAWN
;
; EXAMPLE:
;      
;
; MODIFICATION HISTORY:
;      Written 27-04-2017 P.Kamphuis v1.0
;
; NOTE:
;     
;-
  errormessage=['0','All Good']
  print,'We are using the following files'
  print,'Mask = '+names[0]+'.fits'
  print,'Catalog = '+names[1]
  allnew=2
  cd,new_dir,CURRENT=old_dir
                                ;We need to read the sofia output file     
  openr,1,names[1]
                                ;get the number of lines
  paralines=FILE_LINES(names[1])
                                ;define some variables
  counter=1
  line=''
  ID=0.
  trig=0.
  catvals=0.
                                ;Define the parameters that need to be present and that we want.
  sofia_parameters=['id','x_geo','x_min','x_max','y_geo','y_min','y_max','z_geo','z_min','z_max','x','y','z','f_int']
                                ;Check whether they are present and
                                ;which index number they have
  readf,1,line
  readf,1,line
  sofia_column_ids=strtrim(strcompress(str_sep(strtrim(strcompress(line),2),' ')),2)
  sofia_locations=intarr(n_elements(sofia_parameters))
  for j=0,n_elements(sofia_parameters)-1 do begin
     sofia_locations[j]=WHERE(sofia_parameters[j] EQ sofia_column_ids)-1
     IF sofia_locations[j] EQ -2 then begin
        print,'We cannot find the required column for '+sofia_parameters[i]+' in the sofia catalogue'
        print,'This can happen because a) you have tampered with the sofiainput.txt file in the directory, '+supportdirchecked
        print,'b) you are using an updated version of SoFiA where the names have changed and FAT is not yet updated.'
        print,'   In this case please file a bug report at https://github.com/PeterKamphuis/FAT/issues/'
        print,'c) You are using pre processed SoFiA output of your own and do not have all the output'
        print,'   Required output is '+STRJOIN(sofia_parameters,',')
        stop
     ENDIF
  endfor
                                ;reset file pointer to the beginning of the file
  point_lun,1,0
  
                                ;find the number 1 source
  WHILE ID NE 1 and counter LT paralines do begin
                                ;create a skip array
     skiparray=strarr(paralines-counter)
                                ;skip the header
     readf,1,skiparray
                                ;read the line
     readf,1,line
                                ;break the line
     vals=strtrim(strcompress(str_sep(strtrim(strcompress(line),2),' ')),2)
                                ;set ID
     IF vals[0] EQ '#' then goto,nogame
     ID=fix(vals[sofia_locations[0]])
     IF ID eq 1 then begin
        IF trig EQ 0. then catvals=double(vals) ELSE begin
           tmp=catvals
           catvals=dblarr(n_elements(vals),n_elements(tmp[0,*])+1)
           catvals[*,0:n_elements(tmp[0,*])-1]=tmp[*,*]
           catvals[*,n_elements(tmp[0,*])]=double(vals[*])
        endelse
     ENDIF
     IF ID GT 1 then begin
        IF n_elements(catvals) LT 3 then begin
           catvals=vals
           trig=1
        ENDIF else begin
           tmp=catvals
           catvals=dblarr(n_elements(vals),n_elements(tmp[0,*])+1)
           catvals[*,0:n_elements(tmp[0,*])-1]=tmp[*,*]
           catvals[*,n_elements(tmp[0,*])]=vals[*]
           trig=1
        ENDELSE
     ENDIF
     nogame:
                                ;reset file pointer to the beginning of the file
     point_lun,1,0
                                ;increase the counter
     counter++
  ENDWHILE
  close,1
  IF n_elements(catvals[sofia_locations[0],*]) GT 1 then begin
     vals=strarr(n_elements(catvals[*,0]))
     tmp=WHERE(FINITE(catvals[sofia_locations[13],*]) EQ 0)
     for j=0,n_elements(vals)-1 do begin
        vals[j]=string(MEAN(catvals[j,tmp]))
     endfor
  ENDIF ELSE vals=catvals
                                ;If sofia didn't work and there is no preprocessed stuff then abort
  IF size(log,/TYPE) EQ 7 then begin
     openu,66,log,/APPEND
     printf,66,linenumber()+"We picked the "+string(vals[sofia_locations[0]])+" object of the parameter list."
     close,66
  ENDIF
  
                                ;Read the centre values and make sure
                                ;they fall inside the cube
     
  RApix=[double(vals[sofia_locations[1]]),double(vals[sofia_locations[2]]),double(vals[sofia_locations[3]])]
                                ;IF we read the small cube this needs to be corrected
  
                                ;RApixboun=[double(vals[sofia_locations[2]]),double(vals[sofia_locations[3]])]
  IF allnew EQ 3 then begin
     newsize=sxpar(header,'NAXIS1')-1
     Old=RApix[0]
     RApix[0]=floor(newsize/2.)+(Rapix[0]-floor(RApix[0]))
     pixshift=Old-RApix[0]
     RApix[1:2]=RApixboun[1:2]-pixshift
  ENDIF    
  IF RApix[0] GT sxpar(header,'NAXIS1')-10 then BEGIN
     RApix[0]=double(sxpar(header,'NAXIS1')/2.)
     RApix[1:2]=[double(sxpar(header,'NAXIS1')/4.),double(sxpar(header,'NAXIS1')/4.)*3.]
  ENDIF
  DECpix=[double(vals[sofia_locations[4]]),double(vals[sofia_locations[5]]),double(vals[sofia_locations[6]])]
   ;  DECpixboun=[double(vals[sofia_locations[5]]),double(vals[sofia_locations[6]])]
  IF allnew EQ 3 then begin
     newsize=sxpar(header,'NAXIS1')-1
     Old=DECpix[0]
     DECpix[0]=floor(newsize/2.)+(DECpix[0]-floor(DECpix[0]))
     pixshift=Old-DECpix[0]
     DECpix[1:2]=DECpix[1:2]-pixshift
  ENDIF    
  IF DECpix[0] GT sxpar(header,'NAXIS2')-10 then begin
     DECpix[0]=fix( sxpar(header,'NAXIS2')/2.)
     DECpix[1:2]=[double(sxpar(header,'NAXIS2')/4.),double(sxpar(header,'NAXIS2')/4.)*3.]
  ENDIF
  VSYSpix=[double(vals[sofia_locations[7]]),double(vals[sofia_locations[8]]),double(vals[sofia_locations[9]])]
                                ;   ROTpixboun=[double(vals[sofia_locations[8]]),double(vals[sofia_locations[9]])]
  IF VSYSpix[0] GT sxpar(header,'NAXIS3')-5 then BEGIN
     VSYSpix[0]=fix( sxpar(header,'NAXIS3')/2.)
     VSYSpix[1:2]=[double(sxpar(header,'NAXIS3')/4.),double(sxpar(header,'NAXIS3')/4.)*3.]
  ENDIF
  IF RApix[0] LT RApix[1] OR RApix[0] GT RApix[2] $
     OR DECpix[0] LT DECpix[1] OR DECpix[0] GT DECpix[2] $
     OR VSYSpix[0] LT VSYSpix[1] OR VSYSpix[0] GT VSYSpix[2] then begin
     RApix[0]=double(vals[sofia_locations[10]])
     DECpix[0]=double(vals[sofia_locations[11]])
     VSYSpix[0]=double(vals[sofia_locations[12]])
  ENDIF
  Totflux=[double(vals[sofia_locations[13]])] ;Jy/beam
  
  
  cd,old_dir
end
