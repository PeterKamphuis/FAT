Pro run_sofia,allnew,new_dir,currentfitcube,catcatalogname,supportdirchecked,pixfwhm,header,errormessage,VSYSpix,RApix,DECpix,Totflux,log=log,flagchannels=flagchannels,beams_in_cube =beams_in_cube,channels_in_cube=channels_in_cube
;+
; NAME:
;       RUN_SOFIA
;
; PURPOSE:
;       Run SoFiA to create a mask and find initial central coordinates
;
; CATEGORY:
;       Support
;
; CALLING SEQUENCE:
;       RUN_SOFIA,allnew,new_dir,currentfitcube,catcatalogname,supportdirchecked,pixfwhm,header,errormessage,VSYSpix,RApix,DECpix,Totflux,log=log
;
;
; INPUTS:
;            allnew = trigger determining whether to run SoFiA from scratch
;                     or use provided output.
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
;      flagchannels = channels to be flagged before running SoFiA
;     beams_in_cube = The amount of beams that fit in a cube. Used as a proxy for the galaxy size.
;  channels_in_cube = The amount of the channels in a cube. Used as a proxy for the velocity width of the cube.
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
;       07-03-2017 P.Kamphuis; Setup such that SoFiA output now goes
;                              to the Log file instead of the screen.
;       06-03-2017 P.Kamphuis; Made sure that when there is a central
;                              source and nan is found in the output
;                              file that source is selected.
;       Written 22-02-2017 P.Kamphuis v1.0
;
; NOTE:
;
;-
  errormessage=['0','All Good']
  if allnew GE 2 then goto,skipallsofia
  read_template,supportdirchecked+'/sofiainput.txt',sofia,sofiatriggers,/SOFIA
  CD,new_dir,CURRENT=old_dir
  ;velkern=['0,'+string(39B)+'b'+string(39B),'2,'+string(39B)+'b'+string(39B),'4,'+string(39B)+'b'+string(39B),$
  ;         '8,'+string(39B)+'b'+string(39B),'16,'+string(39B)+'b'+string(39B)]
  threshold= 5.
  counter=3.
  ;We want to take out any fully blanked channels
  ;if flagchannels[0] NE -1 then begin
  ;  sofia[sofiatriggers[5]]='steps.doFlag = True'
;    flag_regions = '['
;    for j=0,n_elements(flagchannels)-1 do begin
;      flag_regions = flag_regions+"[0,'',0,'',"+STRTRIM(string(flagchannels[j]))+','+STRTRIM(string(flagchannels[j]))+'],'
;    ENDFOR
;    flag_regions = STRMID(flag_regions,0,strlen(flag_regions)-1)+']'

;    sofia[sofiatriggers[6]]='flag.regions = '+flag_regions
;  ENDIF

  lowthresagain:
  sofia[sofiatriggers[0]]='import.inFile = '+currentfitcube+'.fits'
  sofia[sofiatriggers[1]]='steps.doReliability             =       false'
  sofia[sofiatriggers[2]]='parameters.dilatePixMax	= '+string(fix(pixfwhm/4*beams_in_cube/10.*counter))
  sofia[sofiatriggers[3]]='SCfind.threshold	=' +string(threshold)
  sofia_kernels = "[ 0, 0, 0,'g'],[ 0, 0, 2,'g'],[ 0, 0, 4,'g'],[ 0, 0, 8,'g'],"+$
                  '['+strtrim(string(fix(pixfwhm)),2)+','+strtrim(string(fix(pixfwhm)),2)+",0,'g'],"+$
                  '['+strtrim(string(fix(pixfwhm)),2)+','+strtrim(string(fix(pixfwhm)),2)+",3,'g'],"+$
                  '['+strtrim(string(fix(pixfwhm)),2)+','+strtrim(string(fix(pixfwhm)),2)+",6,'g'],"+$
                  '['+strtrim(string(fix(pixfwhm)),2)+','+strtrim(string(fix(pixfwhm)),2)+",12,'g'],"+$
                  '['+strtrim(string(fix(pixfwhm)*2),2)+','+strtrim(string(fix(pixfwhm)*2),2)+",0,'g'],"+$
                  '['+strtrim(string(fix(pixfwhm)*2),2)+','+strtrim(string(fix(pixfwhm)*2),2)+",3,'g'],"+$
                  '['+strtrim(string(fix(pixfwhm)*2),2)+','+strtrim(string(fix(pixfwhm)*2),2)+",6,'g'],"+$
                  '['+strtrim(string(fix(pixfwhm)*2),2)+','+strtrim(string(fix(pixfwhm)*2),2)+",12,'g']"
  ;if the cube is very large we add an additional 3 beams scale
  If beams_in_cube GT 30. then BEGIN
    IF size(log,/TYPE) EQ 7 then begin
      openu,66,log,/APPEND
      printf,66,linenumber()+"RUN_SOFIA: Adding an extra kernel scale as the cube is more than 30 beams "
      close,66
    ENDIF
    sofia_kernels = sofia_kernels+',['+strtrim(string(fix(pixfwhm)*3),2)+','+strtrim(string(fix(pixfwhm)*3),2)+",0,'g'],"+$
                  '['+strtrim(string(fix(pixfwhm)*3),2)+','+strtrim(string(fix(pixfwhm)*3),2)+",3,'g'],"+$
                  '['+strtrim(string(fix(pixfwhm)*3),2)+','+strtrim(string(fix(pixfwhm)*3),2)+",6,'g'],"+$
                  '['+strtrim(string(fix(pixfwhm)*3),2)+','+strtrim(string(fix(pixfwhm)*3),2)+",12,'g']"
  endif
  If channels_in_cube GT 52 then begin
    IF size(log,/TYPE) EQ 7 then begin
      openu,66,log,/APPEND
      printf,66,linenumber()+"RUN_SOFIA: Adding an extra kernel scale as the cube has more than 52 channels "
      close,66
    ENDIF
    sofia_kernels = sofia_kernels+',['+strtrim(string(fix(pixfwhm)*0),2)+','+strtrim(string(fix(pixfwhm)*0),2)+",16,'g'],"+$
                  '['+strtrim(string(fix(pixfwhm)*1),2)+','+strtrim(string(fix(pixfwhm)*1),2)+",16,'g'],"+$
                  '['+strtrim(string(fix(pixfwhm)*2),2)+','+strtrim(string(fix(pixfwhm)*2),2)+",16,'g']"
    If beams_in_cube GT 30. then BEGIN
      sofia_kernels = sofia_kernels+',['+strtrim(string(fix(pixfwhm)*3),2)+','+strtrim(string(fix(pixfwhm)*3),2)+",16,'g']"
    ENDIF
  ENDIF

  sofia[sofiatriggers[4]]="SCfind.kernels= ["+sofia_kernels+"]"


                                ;We run sofia to get a bunch of initial estimates
  openw,1,'sofia_input.txt'
  for index=0,n_elements(sofia)-1 do begin
     printf,1,sofia[index]
  endfor
  close,1
  IF FILE_TEST(currentfitcube+'_cat.ascii') then spawn,'rm -f '+currentfitcube+'_cat.ascii'
  IF FILE_TEST(currentfitcube+'_cont.png') then spawn,'rm -f '+currentfitcube+'_cont.png'
  IF FILE_TEST(currentfitcube+'_mask.fits') then spawn,'rm -f '+currentfitcube+'_mask.fits'
  IF FILE_TEST(currentfitcube+'_scat.png') then spawn,'rm -f '+currentfitcube+'_scat.png'
  print,linenumber()+'RUN_SOFIA: Running SoFiA.'
  spawn,'python '+supportdirchecked+'/sofia_pipeline.py sofia_input.txt ',sofiaoutput
  IF size(log,/TYPE) EQ 7 then begin
     openu,66,log,/APPEND
     printf,66,linenumber()+'RUN_SOFIA: The Sofia Output'
     for j=0,n_elements(sofiaoutput)-1 do begin
        printf,66,sofiaoutput[j]
     endfor
     close,66
  ENDIF

  IF not FILE_TEST(currentfitcube+'_cat.ascii') then begin
     IF threshold GT 4 then begin
        IF size(log,/TYPE) EQ 7 then begin
           openu,66,log,/APPEND
           printf,66,linenumber()+"RUN_SOFIA: We did not find a source at a threshold of "+string(threshold)
           printf,66,linenumber()+"RUN_SOFIA: Lowering the threshold and trying again."
           close,66
        ENDIF
        IF counter EQ 1 then counter=3 else counter++
        threshold--
        goto,lowthresagain
     ENDIF ELSE BEGIN
        errormessage=['5','Cannot produce a SoFiA mask']
        goto,endSofiA
     ENDELSE
  endif
  catCatalogname=currentfitcube+'_cat.ascii'
  catmaskname=currentfitcube+'_mask' ;first clean up the cat
  skipallsofia:
                                ;We need to read the sofia output file
  openr,1,catCatalogname
                                ;get the number of lines
  paralines=FILE_LINES(catCatalogname)
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
        print,'RUN_SOFIA: We cannot find the required column for '+sofia_parameters[i]+' in the sofia catalogue'
        print,'RUN_SOFIA: This can happen because a) you have tampered with the sofiainput.txt file in the directory, '+supportdirchecked
        print,'RUN_SOFIA: b) you are using an updated version of SoFiA where the names have changed and FAT is not yet updated.'
        print,'RUN_SOFIA:    In this case please file a bug report at https://github.com/PeterKamphuis/FAT/issues/'
        print,'RUN_SOFIA: c) You are using pre processed SoFiA output of your own and do not have all the output'
        print,'RUN_SOFIA:    Required output is '+STRJOIN(sofia_parameters,',')
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
     tmp=WHERE(FINITE(catvals[sofia_locations[13],*]) EQ 0)
     IF tmp[0] NE -1 then begin
        dummy=readfits(currentfitcube+'.fits',hedtmp)
        ;,/NOSCALE,/SILENT)
        For j=0,n_elements(tmp)-1 do begin
           IF FINITE(dummy[double(catvals[sofia_locations[1],tmp[j]]),double(catvals[sofia_locations[4],tmp[j]]),double(catvals[sofia_locations[7],tmp[j]])]) EQ 0 then $
              catvals[sofia_locations[13],tmp[j]]=1e9-j else begin
              IF ABS(catvals[sofia_locations[1],tmp[j]]-sxpar(hedtmp,'NAXIS1')/2.) LT pixfwhm AND ABS(catvals[sofia_locations[4],tmp[j]]-sxpar(hedtmp,'NAXIS2')/2.) LT pixfwhm AND    $
                 ABS(catvals[sofia_locations[7],tmp[j]]-sxpar(hedtmp,'NAXIS3')/2.) LT 2 then $
                    catvals[sofia_locations[13],tmp[j]]=1e9-j else catvals[sofia_locations[13],tmp[j]]=0.
           ENDELSE
        ENDFOR

     ENDIF
     voxratio=dblarr(n_elements(catvals[0,*]))
     maxvoxel=MAX(catvals[sofia_locations[13],*])
     for j=0,n_elements(catvals[0,*])-1 do begin
        voxratio[j]=maxvoxel/catvals[sofia_locations[13],j]
     endfor

     rmp=WHERE(voxratio LT 3 AND voxratio GT 0)


     IF rmp[0] EQ -1 then begin
        vals=catvals[*,n_elements(catvals[0,*])-1]
     ENDIF ELSE BEGIN
        diff=dblarr(n_elements(rmp))
        for j=0,n_elements(rmp)-1 do begin
           diff[j]=SQRT((sxpar(header,'CRPIX1')-catvals[sofia_locations[1],rmp[j]])^2+(sxpar(header,'CRPIX2')-catvals[sofia_locations[4],rmp[j]])^2)
        endfor
        maxdiff=MAX(diff,MIN=whatwelookfor)
        vals=catvals[*,rmp[WHERE(diff EQ whatwelookfor)]]
     ENDELSE
  ENDIF ELSE vals=catvals
                                ;If sofia didn't work and there is no preprocessed stuff then abort
  IF double(vals[sofia_locations[0]]) EQ 0. AND NOT allnew GE 2 then begin
     errormessage=['5','Sofia failed, the original cube is too small.']
     goto,endSofiA
  ENDIF
  IF size(log,/TYPE) EQ 7 then begin
     openu,66,log,/APPEND
     printf,66,linenumber()+"RUN_SOFIA: We picked the "+string(vals[sofia_locations[0]])+" object of the parameter list."
     close,66
  ENDIF

                                ;If not using preprocessed stuff then we make logical names
  IF allnew LT 2 then begin
     nummask=readfits(catmaskname+'.fits',hedmasknotproper,/SILENT,/NOSCALE)
     tmp=where(nummask NE double(vals[sofia_locations[0]]))
     IF tmp[0] NE -1 then begin
        nummask[tmp]=0.
     ENDIF
     IF vals[0] NE 1. then begin
        tmp=where(nummask EQ double(vals[sofia_locations[0]]))
        nummask[tmp]=1.
     ENDIF
                                ;We want to expand this mask as it is
                                ;usually a bit tight
     ;pixelsizeRA=ABS(sxpar(header,'CDELT1'))

                                ;newmask=fat_smooth(nummask,catmajbeam[i]/(pixelsizeRA*3600.),/MASK)


     sxaddpar,header,'BITPIX',-32
     writefits,currentfitcube+'_binmask.fits',float(nummask),header
     nummask=0.
     newmask=0.
     catmaskname=currentfitcube+'_binmask'
     IF FILE_TEST(currentfitcube+'_cont.png') then spawn,'rm -f '+currentfitcube+'_cont.png'
     IF FILE_TEST(currentfitcube+'_mask.fits') then spawn,'rm -f '+currentfitcube+'_mask.fits'
     IF FILE_TEST(currentfitcube+'_scat.png') then spawn,'rm -f '+currentfitcube+'_scat.png'
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

  endSoFiA:
  cd,old_dir
end
