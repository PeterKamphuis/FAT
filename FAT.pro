Pro FAT,SUPPORT=supportdir,CONFIGURATION_FILE=configfile

;+
; NAME:
;      FAT 
; PURPOSE:
;      Fit Tilted Ring Models with Tirific in a fully automated manner
; CATEGORY:
;      Main for fitting galaxies. Tirific still requires interactive fitting this code attempts
;      to remedy that 
;
; CALLING SEQUENCE:
;      FAT,support='supportdir',configuration_file='configfile' 
;

; INPUTS:
;      - 
; OPTIONAL INPUT KEYWORDS:
;      SUPPORT  = path to the directory where FAT's support
;      routines are located. The default location is ./Support/
;      CONFIGURATION_FILE = A configuration file for FAT. This file
;      should contain the locations of the galaxies to be fitted. See
;      readme for more detailed info.
; OPTIONAL KEYWORD OUTPUT:
;      -
;
; OUTPUTS:
;     See Readme or just run the code 
;
; EXAMPLE:
;     IDL> FAT,CONFIGURATION_FILE='~/Myfittingdirectory/FAT.config'
;
; PROCEDURES CALLED
;  BOOK_KEEPING, BUILDAXII, CALC_EDGE, CHANGERADII, CHECK_CFLUX, CLEANUP, COLUMNDENSITY,
;  CONVERTRADEC, CONVERTSKYANGLEFUNCTION(), EXTRACT_PV, FIT_ELLIPSE(),FAT_SMOOTH()
;  GETDHI, GET_FIXEDRINGSV8, GET_NEWRINGSV8, GET_PROGRESS, 
;  INT_PROFILEV2, INTERPOLATE, ISNUMERIC(), LINENUMBER(), MOMENTSV2,
;  OBTAIN_INCLINATIONV8, OBTAIN_PAV2, OBTAIN_VELPA, OBTAIN_W50, 
;  PARAMETERREGUV87, READ_TEMPLATE, SBR_CHECK, SET_SBR, SET_VROTV6,
;  SET_WARP_SLOPEV2, WRITEFITTINGVARIABLES, WRITENEWTOTEMPLATE,
;  RESOLVE_ROUTINE, STRLOWCASE, and likely more.
;
; MODIFICATION HISTORY:
;      04-01-2016 P.Kamphuis; Added fail conditions when final model
;      is outside the boundaries set in Kamphuis et al. 2015. 
;      04-01-2016 P.Kamphuis; Removed the usage of the rename command
;      and replaced it with a tirific specific idl routine.
;      04-01-2016 P.Kamphuis; Extensive editing of output messages to
;      keep log clear and readable.
;      04-01-2016 P.Kamphuis; Minor bug fixes.
;      29-12-2016 P.Kamphuis: Replaced GAUSS_SMOOTH and SMOOTH with FAT_SMOOTH
;      for increased compatibility with IDL 7.0 and GDL 
;      29-12-2015 P.Kamphuis; Updated interaction to the latest
;      version of SoFiA. Also introduced column recognition in order
;      to facilitate future changes more easily. Additionally this
;      allows the user to add extra output to the SofiA file. S+C
;      threshold is changed from 4.0 to 5.5 to get similar masks dues
;      to SoFiA moving back to their initial way of defining the
;      threshold. SNR values no longer present in SoFiA.
;      29-12-2015 P.Kamphuis; There is an issue with the usage of
;      NAXIS3 as the commit history show going back and forth between
;      making this an integer or Double precision. I do not remember
;      the reason for going to double precision but it has been
;      reversed now because SoFiA requires and integer (As it should be)  
;      12-08-2015 P. Kamphuis; Improved Bookkeeping for failed fits,
;      improved noise determination  and fixed an issue with spaces in
;      the config file
;      10-08-2015 P. Kamphuis; On line 4012 it was possible to set the
;      rotation curve to zero by having the inner INCL and PA hit the
;      border. Made sure that if taking ring 0 for PA and INCL that
;      the mean is taken for VROT.
;      30-07-2015 P. Kamphuis; Added a set of routines for
;      bookkeeping. Still need to implement proper removal.  
;      24-07-2015 P. Kamphuis; Added routine for the creation of residual files
;      24-07-2015 P. Kamphuis; Replaced the usage of SNR from Sofia
;      with a proper flux from the integrated intensity map. Improved
;      clean up and cflux handling. Also improved treatment of small
;      galaxies and a better determination of CONDISP.  
;      Written by P.Kamphuis 01-01-2015 
; NOTES:
                                ; please note that this pipeline uses the following command line
                                ; commands. Ensure that these work
                                ; under IDL
                                ;the original unix rename working as
                                ;---> rename string replacestring  inputfiles
  
                                ; Standard unix commands pwd,mkdir,rm,
                                ; cp, ls, python
  
                                ; Tirific (any standalone version) should be properly installed
                                ; and callable from the command line
                                ; as tirific.

                                ; SoFIA (The pipeline is updated to be
                                ; compatible with the latest version,
                                ; if there is any problem please let
                                ; me know.) should be properly installed and the
                                ; python file sofia_pipeline.py should
                                ; be copied or linked into the Support
                                ; directory

                                ;A standard up to date IDL (v7.0 or
                                ;higher) installation with ASTRO LIB
                                ;installed should be enough to run the
                                ;program

                                ;All Fully Automated Tirific (FAT) specific routines are availaible in the
                                ;Support directory

                                ;The directory where the support
                                ;routines for the program are located
                                ;Python needs this to be a full path
                                ;on my mac so we'll set that
                                ;with a pwd just in case


                          
  COMPILE_OPT IDL2 

  version='v3.1'
  if n_elements(supportdir) EQ 0 then supportdir='Support'
  CD,supportdir,CURRENT=old_dir
  spawn,'pwd',supportdir
  cd,old_dir
  supportdir=supportdir[0]+'/'
                                ;Extension name for a basic information (e.g HI diameter/mass central
                                ;position) file
  basicinfo='BasicInfo'
  spacecheck=strtrim(str_sep(supportdir,' '),2)
  supportdirchecked=STRJOIN(spacecheck,'\ ')
                                ;Location of the input configuration file

  IF n_elements(configfile) EQ 0 then configfile='FAT_INPUT.config'

                                ;dumb resolve routine doesn't accept paths so cd to support dir

  spawn,'pwd',originaldir
  CD,supportdir,current=old_dir
  RESOLVE_ROUTINE, 'book_keeping'
  RESOLVE_ROUTINE, 'buildaxii'
  RESOLVE_ROUTINE, 'calc_edge'
  RESOLVE_ROUTINE, 'changeradii'
  RESOLVE_ROUTINE, 'check_cflux'
  RESOLVE_ROUTINE, 'cleanup'
  RESOLVE_ROUTINE, 'columndensity'
  RESOLVE_ROUTINE, 'convertradec'
  RESOLVE_ROUTINE, 'convertskyanglefunction',/IS_FUNCTION
  RESOLVE_ROUTINE, 'create_residuals'
  RESOLVE_ROUTINE, 'extract_pv'
  RESOLVE_ROUTINE, 'fat_smooth',/IS_FUNCTION
  RESOLVE_ROUTINE, 'fit_ellipse',/IS_FUNCTION
  RESOLVE_ROUTINE, 'getdhi'
  RESOLVE_ROUTINE, 'get_fixedringsv8'
  RESOLVE_ROUTINE, 'get_newringsv8'
  RESOLVE_ROUTINE, 'get_progress'
  RESOLVE_ROUTINE, 'int_profilev2'
  RESOLVE_ROUTINE, 'interpolate'
  RESOLVE_ROUTINE, 'isnumeric',/IS_FUNCTION
  RESOLVE_ROUTINE, 'linenumber',/IS_FUNCTION
  RESOLVE_ROUTINE, 'momentsv2'
  RESOLVE_ROUTINE, 'obtain_inclinationv8'
  RESOLVE_ROUTINE, 'obtain_pav2'
  RESOLVE_ROUTINE, 'obtain_velpa'
  RESOLVE_ROUTINE, 'obtain_w50'
  RESOLVE_ROUTINE, 'organize_output'
  RESOLVE_ROUTINE, 'parameterreguv87'
  RESOLVE_ROUTINE, 'read_template'
  RESOLVE_ROUTINE, 'rename'
  RESOLVE_ROUTINE, 'sbr_check'
  RESOLVE_ROUTINE, 'set_sbr'
  RESOLVE_ROUTINE, 'set_vrotv6'
  RESOLVE_ROUTINE, 'set_warp_slopev2'
  RESOLVE_ROUTINE, 'sum',/IS_FUNCTION
  RESOLVE_ROUTINE, 'writefittingvariables'
  RESOLVE_ROUTINE, 'writenewtotemplate'
  CD,old_dir

tryconfigagain:
                                ;open the configuration file and get the proper parameters
  fileexist=FILE_TEST(configfile)
                                ;however sometimes the config is not found
  IF strlen(configfile) NE 0 AND fileexist EQ 0 then begin
     print,"You have provided a config file but it can't be found"
     print,"If you want to provide a config file please give the name"
     print,"Else press enter to continue"
     disp=''
     WHILE (disp EQ '') DO BEGIN
        read,disp
        IF disp EQ '' then begin
           configfile=''
           disp='hh'
        ENDIF else begin
           configfile=disp
        ENDELSE
     ENDWHILE
     goto,tryconfigagain
  ENDIF
  IF fileexist EQ 0 then goto,noconfig
  close,1
  h=' '
  openr,1,configfile
  filelength=FILE_LINES(configfile)-1
  WHILE ~ EOF(1) DO BEGIN
     readf,1,h
     tmp=str_sep(strtrim(strcompress(h),2),'=')
     IF n_elements(tmp) GT 1 then begin
        case STRLOWCASE(strtrim(tmp[0],2)) of
                                ;thefirst galaxy to be fitted
           'startgalaxy':startgalaxy=double(tmp[1])
                                ;thelast galaxy to be fitted
                                ;if set to -1 the whole catalogue will be fitted
           'endgalaxy':endgalaxy=double(tmp[1])
                                ;parameters to skip certain fits for testing
                                ;Files need to be present as it will only skip the tirific fitting part
                                ;if you want to run only part of the code use testing remainder and will go on
                                ;to the next galaxy
           'testing':testing=double(tmp[1])
                                ;Remove all previous pre-processing as cutting the cube 0 no 1 yes 2
                                ;use provided preprocessing
           'allnew':allnew=double(tmp[1])
                                ;Parameters for finishing the fitting process early
                                ;IF set to one the program finishes after this loop (The first 1 it encounters)
           'finishafter':finishafter=double(tmp[1])
                                ;Input catalogue for the pipeline
           'catalogue':catalogue=tmp[1]
                                ;Directory with all the directories of the galaxies to be fitted
           'maindir':maindir=tmp[1]
                                ;Output file for the fitting results
           'outputcatalogue':outputcatalogue=tmp[1]
                                ;Output file for the log of the fit
           'outputlog':olog=tmp[1]
                                ;trigger for creating new output
                                ;results file 'y' or append the old one 'n'
           'new_output':newresult=tmp[1]
                                ;trigger for creating a new output
                                ;log file 'y' or append the old one 'n'
           'new_log':newlog=tmp[1]
                                ;Optimal number of pixels per maj axis beam
           'opt_pixelbeam':optpixelbeam=fix(tmp[1])
                                ;Hanning smoothed or not (y=1, n=0)
           'velocity_resolution':vresolution=double(tmp[1])
                                ;Clean up after the fitting is complete.
           'maps_output':bookkeeping=double(tmp[1])
           else:begin
           end
        endcase
     endif
  ENDWHILE
  close,1
 
noconfig:
                                ;Make idiot failsafe standards for the config files 
  IF size(maindir,/TYPE) NE 7 then begin
     print,'You have not provided the maindirectory for the data cubes'
     disp=''
     WHILE (disp EQ '') DO BEGIN
        print,'Please provide the maindir for the data (Default = home)'
        read,disp
        IF disp EQ '' then disp='~'
        filetest=FILE_TEST(disp)
        IF filetest EQ 0 then begin
           print,'That is not a valid directory please try again'
           disp=''
        ENDIF
     ENDWHILE
     maindir=disp
  ENDIF
  IF size(catalogue,/TYPE) NE 7 then begin
     print,'You have not provided a catalogue with directories and initial estimates'
     disp=''
     WHILE (disp EQ '') DO BEGIN
        print,'Please provide the catalogue with initial estimates (No Default)'
        read,disp
        filetest=FILE_TEST(disp)
        IF filetest EQ 0 then begin
           print,'That is not a valid file, please try again'
           disp=''
        ENDIF
     ENDWHILE
     catalogue=disp
  ENDIF
  IF size(outputcatalogue,/TYPE) NE 7 then begin
     disp=''
     WHILE (disp EQ '') DO BEGIN
        print,'Please provide the file for results (Default = '+maindir+'/output_dir/fit_result.txt)'
        read,disp
        if disp EQ '' then begin
           wehaveadir=FILE_TEST(maindir+'/output_dir')
           IF wehaveadir EQ 0 then begin
              spawn,'mkdir '+maindir+'/output_dir/',isthere
           ENDIF
           disp=maindir+'output_dir/fit_results.txt'
        endif
     ENDWHILE
     outputcatalogue=disp
  ENDIF
                                ;defaults for configuration file
  IF size(finishafter,/TYPE) EQ 0 then finishafter=2
  IF size(testing,/TYPE) EQ 0 then testing=0
  IF size(startgalaxy,/TYPE) EQ 0 then startgalaxy=0
  IF size(endgalaxy,/TYPE) EQ 0 then endgalaxy=-1
  IF size(newlog,/TYPE) EQ 0 then newlog='y'
  IF size(newresult,/TYPE) EQ 0 then newresult='y'
  IF n_elements(vresolution) EQ 0 then vresolution=1.
  IF n_elements(optpixelbeam) EQ 0 then optpixelbeam=4.
  IF n_elements(allnew) EQ 0 then allnew=1
  IF n_elements(bookkeeping) EQ 0 then bookkeeping=3
 
  IF bookkeeping EQ 5. then bookkeeping=4

  bookkeepingin=bookkeeping
  finishafterold=finishafter
                                ;start a log with current program name
  help,calls=cc
  print,linenumber()+"This is version "+version[0]+" of the program"
                                ;Create a file to write the results to
  If not FILE_TEST(outputcatalogue) OR strupcase(newresult) EQ 'Y' then begin
     openw,1,outputcatalogue
     printf,1,format='(A60,2A12)','Name','AC1','AC2'
     close,1
  endif
                                ;read the input catalogue
  openr,1,catalogue
  filelength=FILE_LINES(catalogue)-1
  catnumber=strarr(filelength)
  catRA=strarr(filelength)
  catDEC=strarr(filelength)
  epoch=strarr(filelength)
  catPA=dblarr(filelength)
  catPAdev=dblarr(filelength)
  catinc=dblarr(filelength)
  catincdev=dblarr(filelength)
  catvsys=dblarr(filelength)
  catmaxrot=dblarr(filelength)
  catmaxrotdev=dblarr(filelength)
  catDistance=dblarr(filelength)
  catnoise=dblarr(filelength)
  catmajbeam=dblarr(filelength)
  catminbeam=dblarr(filelength)
  catbpa=dblarr(filelength)
  catDirname=strarr(filelength)
  catCubename=strarr(filelength)
  catMom0name=strarr(filelength)
  catMom1name=strarr(filelength)
  catMaskname=strarr(filelength)
  catCatalogname=strarr(filelength)
  h=' '
  readf,1,h
  tmp=str_sep(strtrim(strcompress(h),2),'|')
                                ;IF we have a new catalog we need to
                                ;determine what each column is
                                ;referring to
  IF n_elements(tmp) LT 16 then begin
                                ;When using the Sofia input it is not
                                ;necessary to have a specific ordering
                                ;of the input catalogue just the
                                ;correct column names
     tmpnumber=WHERE(STRLOWCASE(tmp) EQ 'number')
     tmpdistance=WHERE(STRLOWCASE(tmp) EQ 'distance')
     tmpdir=WHERE(STRLOWCASE(tmp) EQ 'directoryname')
     tmpcube=WHERE(STRLOWCASE(tmp) EQ 'cubename')
     IF n_elements(tmp) GT 4 then begin
        test=WHERE(STRLOWCASE(tmp) EQ 'basename')
        IF test[0] EQ -1 then begin
           tmpmom1=WHERE(STRLOWCASE(tmp) EQ 'moment1')
           tmpmom0=WHERE(STRLOWCASE(tmp) EQ 'moment0')
           tmpmask=WHERE(STRLOWCASE(tmp) EQ 'mask')
           tmpcatalog=WHERE(STRLOWCASE(tmp) EQ 'catalog')
                                ;If we want to use preprocessed input
                                ;than we need moment0 moment1 a mask
                                ;and a catalog
           IF allnew GT 1 AND (tmpmom1[0] EQ -1 OR tmpmom0[0] EQ -1 OR tmpmask[0] EQ -1 OR tmpcatalog[0] EQ -1) then begin
              print,'Your catalog file is not configured properly. Aborting'
              stop
           ENDIF
        ENDIF
     ENDIF ELSE BEGIN
                                ;If we want to use preprocessed input
                                ;than we need moment0 moment1 a mask
                                ;and a catalog
        IF allnew GT 1 then begin
           print,'Your catalog file is not configured properly. Aborting'
           stop
        ENDIF
     ENDELSE
  ENDIF
                                ;Then read the actual input from the catalo
  h=' '
  for i=0,filelength-1 do begin
     readf,1,h
     tmp=str_sep(strtrim(strcompress(h),2),'|')
     IF n_elements(tmp) GE 16 then begin 
                                ;Old fashioned input of initial
                                ;estimates from the literature
                                ;Better to us the Sofia estimates as
                                ;the program is fine tuned to them
        catnumber[i]=tmp[0]
        catRA[i]=tmp[1]
        catDEC[i]=tmp[2]
        epoch[i]=tmp[3]
        catPA[i]=double(tmp[4])
        catinc[i]=double(tmp[5])
        catincdev[i]=double(tmp[6])
        catvsys[i]=double(tmp[7])
        catmaxrot[i]=double(tmp[8])
        catDistance[i]=double(tmp[9])
        IF catDistance[i] LT 1. THEN catDistance[i]=1.
        catnoise[i]=double(tmp[10])
        catmajbeam[i]=double(tmp[11])
        catminbeam[i]=double(tmp[12])
        catbpa[i]=double(tmp[13])
        catDirname[i]=tmp[14]
        tmpagain=str_sep(tmp[15],'.')        
        IF n_elements(tmpagain) GT 1 and STRLOWCASE(tmpagain[n_elements(tmpagain)-1]) EQ 'fits' then begin           
           IF n_elements(tmpagain) GT 2 then catCubename[i]=STRJOIN(tmpagain[0:n_elements(tmpagain)-2],'.') else catCubename[i]=tmpagain[0]
        endif else catCubename[i]=tmp[15]
        IF n_elements(tmp) GT 16 then begin
           catMom0name[i]=tmp[16]
           catMom1name[i]=tmp[17]
        ENDIF ELSE BEGIN
           catMom0name[i]='youshouldnotnameyourfilesthisway'
           catMom1name[i]='youshouldnotnameyourfilesthisway'
        ENDELSE
     ENDIF ELSE BEGIN
                                ;IF we have a new catalog we use the determined order
        catnumber[i]=tmp[tmpnumber]
        catDirname[i]=tmp[tmpdir]
        tmpagain=str_sep(tmp[tmpcube],'.')
        IF n_elements(tmpagain) GT 1 and STRLOWCASE(tmpagain[n_elements(tmpagain)-1]) EQ 'fits' then begin
           IF n_elements(tmpagain) GT 2 then catCubename[i]=STRJOIN(tmpagain[0:n_elements(tmpagain)-2],'.') else catCubename[i]=tmpagain[0]
        endif else catCubename[i]=tmp[tmpcube]
        catDistance[i]=double(tmp[tmpdistance])
                                ;IF we have preprocessed data we read the file names
        IF allnew GT 1 then begin
           IF test[0] EQ -1 then begin
              catMom0name[i]=tmp[tmpmom0]
              catMom1name[i]=tmp[tmpmom1]
              catMaskname[i]=tmp[tmpmask]
              catCatalogname[i]=tmp[tmpcatalog]
              tmpagain=str_sep(tmp[tmpcatalog],'.')
              IF n_elements(tmpagain) GT 1 and STRLOWCASE(tmpagain[n_elements(tmpagain)-1]) NE 'ascii' then catCatalogname[i]=tmp[tmpcatalog]+'.ascii'
           ENDIF ELSE BEGIN
              catMom0name[i]=tmp[test]+'_mom0'
              catMom1name[i]=tmp[test]+'_mom1'
              catMaskname[i]=tmp[test]+'_mask'
              catCatalogname[i]=tmp[test]+'_cat.ascii'
           ENDELSE
                                ; else we make sure they have a
                                ; silly name so they don't
                                ; exist in the directory
        ENDIF ELSE BEGIN
           catMom0name[i]='youshouldnotnameyourfilesthisway'
           catMom1name[i]='youshouldnotnameyourfilesthisway'
        ENDELSE
     ENDELSE

  endfor
  close,1
                                ;if our end galaxy is -1 then we want to fit the full catalog
  IF endgalaxy EQ -1 then begin
     endgalaxy=filelength-1
  ENDIF
                                ;start the main galaxy fitting loop
  for i=startgalaxy,endgalaxy do begin 
                                ;To ensure using the right template we
                                ;read them everytime we start a new galaxy
     finishafter=finishafterold
     bookkeeping=bookkeepingin
     noisemapname=catCubename[i]+'_6.0_noisemap'
     close,1
                                ;Start a galaxy specific fitting log
     IF size(olog,/TYPE) EQ 7 then begin
        log=maindir+'/'+catdirname[i]+'/'+olog
        IF not FILE_TEST(log) OR strupcase(newlog) EQ 'Y' then begin
           openw,66,log
           printf,66,linenumber()+"This file is a log of the fitting process run at "+systime()
           printf,66,linenumber()+"This is version "+version[0]+" of the program."
           close,66
        ENDIF ELSE BEGIN
           openu,66,log,/APPEND
           printf,66,linenumber()+"This file is a log of the fitting process continued at "+systime()
           printf,66,linenumber()+"This is version "+version[0]+" of the program."
           close,66
        ENDELSE
     ENDIF
                                ;Read the template for the first
                                ;tirific fit 
     read_template,supportdir+'/1stfit.def', tirificfirst, tirificfirstvars
                                ;Read the template for the second
                                ;tirific fit 
     read_template,supportdir+'/2ndfit.def', tirificsecond, tirificsecondvars
                                ;Read the template sofia input file
     IF allnew LT 2 then begin
        read_template,supportdir+'/sofiainput.txt',sofia,sofiatriggers,/SOFIA
     ENDIF
                                ;print which galaxy we are at
     print,linenumber()+"We're at galaxy number "+strtrim(string(i,format='(I10)'),2)+". Which is catalogue id number "+catnumber[i]
     IF size(log,/TYPE) EQ 7 then begin
        openu,66,log,/APPEND
        printf,66,linenumber()+"We're at galaxy number "+strtrim(string(i,format='(I10)'),2)+". Which is catalogue id number "+catnumber[i]
        close,66
     ENDIF
                                ;set the name of the cube to be fitted
     currentfitcube=catcubename[i]
                                ;first some extension detection In
                                ;case somebody is still using gipsy file
     cubeext=' '
     new_dir=maindir+'/'+catdirname[i]+'/'
   
     CD,new_dir,CURRENT=old_dir
                                ;To avoid confusion we remove all
                                ;initial fits and subsequent fits that
                                ;are present if we want to start all
                                ;over again
     IF allnew EQ 1 then begin
        IF testing LT 1 then cleanup,catcubename[i]               
     ENDIF
     IF testing LT 1 then spawn,'rm -f 1stfit*',isthere
     IF testing LT 2 then spawn,'rm -f 2ndfit*',isthere
     IF testing LT 3 then spawn,'rm -f 3rdfit*',isthere
                                ;Let's see what kind of input cube we have
     spawn,'ls '+catcubename[i]+'.*',filenameext
     CD,OLD_DIR
     IF (filenameext[0]) then begin
        for j=0,n_elements(filenameext)-1 do begin
           tmp=str_sep(strtrim(strcompress(filenameext[j]),2),'.')
           IF n_elements(tmp) LT 3 then begin
              cubeext='.fits'
              goto,extset
           ENDIF
           for g=n_elements(tmp)-1,n_elements(tmp)-3,-1 do begin
              IF STRUPCASE(tmp[g]) EQ 'FITS' then begin
                 case g of
                    n_elements(tmp)-1:begin
                       nametmp=tmp[0]
                       for z=1,n_elements(tmp)-2 do begin
                          nametmp=nametmp+'.'+tmp[z]
                       endfor
                       IF nametmp EQ catcubename[i] then begin
                          cubeext='.'+tmp[g]
                          goto,extset
                       ENDIF
                    end
                    n_elements(tmp)-2:begin
                       IF STRUPCASE(tmp[g+1]) EQ 'GZ' then begin
                          nametmp=tmp[0]
                          for z=1,n_elements(tmp)-3 do begin
                             nametmp=nametmp+'.'+tmp[z]
                          endfor
                          IF nametmp EQ catcubename[i] then begin
                             cubeext='.'+tmp[g]+'.'+tmp[g+1]
                             goto,extset
                          ENDIF
                          
                       ENDIF
                    END
                 endcase
              endif
           endfor
        endfor
        extset:
     endif

                                ;See if we already cut the initial cube to something more accesible
     smallexists=FILE_TEST(maindir+'/'+catdirname[i]+'/'+catCubename[i]+'_small.fits')
     IF smallexists then begin
        cubeext='.fits'
        currentfitcube=catcubename[i]+'_small'
        catcubename[i]=catcubename[i]+'_small'
     ENDIF
                                ;and we make sure that all the naming
                                ;stuff adds up
     fitsexists=FILE_TEST(maindir+'/'+catdirname[i]+'/'+catcubename[i]+cubeext)


                                ;Let's go to the right directory
     IF size(log,/TYPE) EQ 7 then begin
        openu,66,log,/APPEND
        printf,66,linenumber()+"Working in Directory "+new_dir
        close,66
     ENDIF ELSE BEGIN
        print,linenumber()+"Working in Directory "+new_dir
     ENDELSE
                                ;Break if the cube cannot be found
     if fitsexists EQ 0 then begin
        openu,1,outputcatalogue,/APPEND
        printf,1,format='(A60,A100)',catDirname[i],' This galaxy has no fits cube to work with, it is skipped.'
        close,1   
        IF size(log,/TYPE) EQ 7 then begin
           openu,66,log,/APPEND
           printf,66,linenumber()+catDirname[i]+" This galaxy has no fits cube to work with, it is skipped."
           close,66
        ENDIF ELSE BEGIN
           print,linenumber()+catDirname[i]+" This galaxy has no fits to work with, it is skipped."
        ENDELSE
        bookkeeping=5
        optimized=0
        goto,finishthisgalaxy
     endif
                                ;set some triggers at 0
     maxrings=0.
     propermom0=0
     propermom1=0
     optimized=0
     sofiacount=0
     sofiafail=0
     wroteblank=0
                                ;read in the cube and check it's header
     dummy=readfits(maindir+'/'+catdirname[i]+'/'+catcubename[i]+cubeext,hed,/NOSCALE)
     howmanyaxis=sxpar(hed,'NAXIS')
     writecube=0
     IF howmanyaxis GT 3 then begin
        sxaddpar,hed,'NAXIS',3
                                ;whenever we change something we want to rewrite the cube
        writecube=1
     ENDIF
     IF sxpar(hed,'NAXIS4') then begin
        sxdelpar,hed,'NAXIS4'
        writecube=1
     ENDIF
     cubecrvalRA=sxpar(hed,'CRVAL1')
     cubecdelt=(ABS(sxpar(hed,'CDELT1'))+ABS(sxpar(hed,'CDELT2')))/2.
     cubecrvalDEC=sxpar(hed,'CRVAL2')
     channelwidth=ABS(sxpar(hed,'CDELT3'))   
     veltype=strtrim(strcompress(sxpar(hed,'CUNIT3')))
     IF STRUPCASE(veltype) EQ 'HZ' then begin
        IF size(log,/TYPE) EQ 7 then begin
           openu,66,log,/APPEND
           printf,66,linenumber()+'FREQUENCY IS NOT A SUPPORTED VELOCITY AXIS.'          
           close,66
        ENDIF ELSE BEGIN
           print,linenumber()+'FREQUENCY IS NOT A SUPPORTED VELOCITY AXIS.'    
        ENDELSE
        bookkeeping=5
        GOTO,FINISHTHISGALAXY
     ENDIF
     IF isnumeric(veltype) then begin
        IF channelwidth GT 100. then begin
           veltype='M/S'
           sxaddpar,hed,'CUNIT3','M/S',after='CDELT3'
           writecube=1
        endif else begin 
           veltype='KM/S'
           sxaddpar,hed,'CUNIT3','KM/S',after='CDELT3'
           writecube=1
        endelse
        IF size(log,/TYPE) EQ 7 then begin
           openu,66,log,/APPEND
           printf,66,linenumber()+'Your header did not have a unit for the third axis, that is bad policy.'          
           printf,66,linenumber()+'We have set it to '+veltype+'. Please ensure that is correct.'          
           close,66
        ENDIF ELSE BEGIN
           print,linenumber()+'Your header did not have a unit for the third axis, that is bad policy.'          
           print,linenumber()+'We have set it to '+veltype+'. Please ensure that is correct.'     
        ENDELSE
     ENDIF
     velproj=sxpar(hed,'CTYPE3')
     IF STRUPCASE(strtrim(velproj,2)) NE 'VELO-HEL' AND $
        STRUPCASE(strtrim(velproj,2)) NE 'VELO-LSR' AND $
        STRUPCASE(strtrim(velproj,2)) NE 'FELO-HEL' AND $
        STRUPCASE(strtrim(velproj,2)) NE 'FELO-LSR' AND $
        STRUPCASE(strtrim(velproj,2)) NE 'VELO' AND $
        STRUPCASE(strtrim(velproj,2)) NE 'FREQ' then begin
        sxaddpar,hed,'CTYPE3','VELO',after='CUNIT3'
        writecube=1
        IF size(log,/TYPE) EQ 7 then begin
           openu,66,log,/APPEND
           printf,66,linenumber()+'Your velocity projection is not standard. The keyword is changed to VELO (relativistic definition). This might be dangerous.'          
           close,66
        ENDIF ELSE BEGIN
           print,linenumber()+'Your velocity projection is not standard. The keyword is changed to VELO (relativistic definition). This might be dangerous.'   
        ENDELSE
     ENDIF
     IF STRUPCASE(veltype) EQ 'KM/S' then begin
         IF size(log,/TYPE) EQ 7 then begin
           openu,66,log,/APPEND
           printf,66,linenumber()+'The channels in your input cube are in km/s. This sometimes leads to problems with wcs lib, hence we change it to m/s.'          
           close,66
        ENDIF ELSE BEGIN
           printf,66,linenumber()+'The channels in your input cube are in km/s. This sometimes leads to problems with wcs lib, hence we change it to m/s.'   
        ENDELSE
        sxaddpar,hed,'CDELT3',sxpar(hed,'CDELT3')*1000.
        sxaddpar,hed,'CRVAL3',sxpar(hed,'CRVAL3')*1000.
        sxaddpar,hed,'CUNIT3','M/S'
        channelwidth=ABS(sxpar(hed,'CDELT3'))   
        veltype=strtrim(strcompress(sxpar(hed,'CUNIT3')))
     ENDIF
        
     sizex=sxpar(hed,'NAXIS1')
     sizey=sxpar(hed,'NAXIS2')
                                ;Let's check for presence of the beam in the header. IF present
                                ;supersede the input. !!!!Be careful with smoothed data. If not let's add
                                ;it from the file; This is important if you do it incorrectly the
                                ;pipeline will use a incorrect pixels size.
                                ;   writecube=0
     IF NOT sxpar(hed,'BMAJ') then begin
        IF sxpar(hed,'BMMAJ') then begin
           sxaddpar,hed,'BMAJ',sxpar(hed,'BMMAJ')/3600. 
           catmajbeam[i]=sxpar(hed,'BMMAJ')
        endif  else sxaddpar,hed,'BMAJ',double(catmajbeam[i]/3600.)
        writecube=1
     ENDIF else catmajbeam[i]=sxpar(hed,'BMAJ')*3600.
     IF NOT sxpar(hed,'BMIN') then begin
        IF sxpar(hed,'BMMIN') then begin
           sxaddpar,hed,'BMIN',sxpar(hed,'BMMIN')/3600. 
           catminbeam[i]=sxpar(hed,'BMMIN')
        endif else sxaddpar,hed,'BMIN',double(catminbeam[i]/3600.)
        writecube=1
     ENDIF else catminbeam[i]=sxpar(hed,'BMIN')*3600.
                                ; check for 0 values, they
                                ; shouldn't be present and mess
                                ; things up
                                ;IF present assume they are blanks and
                                ;add warning
  
     tmp=WHERE(dummy EQ 0.d)
     IF tmp[0] NE -1 then begin
        IF size(log,/TYPE) EQ 7 then begin
           openu,66,log,/APPEND
           printf,66,linenumber()+'Your cube contains values exactly 0. If this is padding these should be blanks.'
           printf,66,linenumber()+'We have changed them to blanks.'
           close,66
        ENDIF ELSE BEGIN
           print,linenumber()+'Your cube contains values exactly 0. If this is padding these should be blanks.'
           print,linenumber()+'We have changed them to blanks.'
        ENDELSE
        dummy[tmp]=!values.f_nan
        writecube=1
     ENDIF
                                ;If we changed something we rewrite the cube
     IF writecube EQ 1 then writefits,maindir+'/'+catdirname[i]+'/'+catcubename[i]+cubeext,dummy,hed
     writecube=0.
     IF sxpar(hed,'BPA') then catbpa[i]=sxpar(hed,'BPA')
                                ;We will now optimize the cube for speed such that there are n pixels
                                ;per FWHM on the minor axis
     
                                ;first we check wether the cube has some proper noise statistics by
                                ;comparing the last and first channel, if they differ to much we cut
                                ;off the high noise channel
     difference = 10.
     changedcube=0.
     WHILE difference GT 1 do begin
       
        tmpnoblank=dummy[*,*,0]
        wherefinite=WHERE(FINITE(tmpnoblank))
        WHILE wherefinite[0] EQ -1 DO begin
           tmp=fltarr(n_elements(dummy[*,0,0]),n_elements(dummy[0,*,0]),n_elements(dummy[0,0,*])-1)
           tmp[*,*,0:n_elements(tmp[0,0,*])-1]=dummy[*,*,1:n_elements(dummy[0,0,*])-1]
           sxaddpar,hed,'CRPIX3',sxpar(hed,'CRPIX3')-1.
           dummy=fltarr(n_elements(tmp[*,0,0]),n_elements(tmp[0,*,0]),n_elements(tmp[0,0,*]))
           IF n_elements(dummy[0,0,*]) LT 5 then begin
              IF size(log,/TYPE) EQ 7 then begin
                 openu,66,log,/APPEND
                 printf,66,linenumber()+catdirname[i]+'/'+catcubename[i]+' has too many blanked channels.'
                 close,66
              ENDIF ELSE BEGIN
                 printf,linenumber()+catdirname[i]+'/'+catcubename[i]+'has too many blanked channels.'
              ENDELSE         
              openu,1,outputcatalogue,/APPEND
              printf,1,format='(A60,A90)', catDirname[i],'The Cube has too many blanked channels'
              close,1    
              bookkeeping=5
              goto,finishthisgalaxy
           ENDIF
           IF size(log,/TYPE) EQ 7 then begin
              openu,66,log,/APPEND
              printf,66,linenumber()+' We are cutting the cube as the first channel is completely blank.'
              close,66
           ENDIF
           sxaddpar,hed,'NAXIS3',fix(sxpar(hed,'NAXIS3')-1)
           changedcube=1
           dummy=tmp
           tmpnoblank=dummy[*,*,0]
           wherefinite=WHERE(FINITE(tmpnoblank))
        ENDWHILE
        rmsfirstchannel=SIGMA(tmpnoblank[WHERE(FINITE(tmpnoblank))])          
        tmpnoblank=dummy[*,*,n_elements(dummy[0,0,*])-1]
        wherefinite=WHERE(FINITE(tmpnoblank))
        WHILE wherefinite[0] EQ -1 DO begin
           tmp=fltarr(n_elements(dummy[*,0,0]),n_elements(dummy[0,*,0]),n_elements(dummy[0,0,*])-1)
           tmp[*,*,0:n_elements(tmp[0,0,*])-1]=dummy[*,*,0:n_elements(dummy[0,0,*])-2]
           dummy=fltarr(n_elements(tmp[*,0,0]),n_elements(tmp[0,*,0]),n_elements(tmp[0,0,*]))
           IF n_elements(dummy[0,0,*]) LT 5 then begin
              IF size(log,/TYPE) EQ 7 then begin
                 openu,66,log,/APPEND
                 printf,66,linenumber()+catdirname[i]+'/'+catcubename[i]+' has too many blanked channels.'
                 close,66
              ENDIF ELSE BEGIN
                 printf,linenumber()+catdirname[i]+'/'+catcubename[i]+'has too many blanked channels.'
              ENDELSE         
              openu,1,outputcatalogue,/APPEND
              printf,1,format='(A60,A90)', catDirname[i],'The cube has too many blanked channels'
              close,1    
              bookkeeping=5
              goto,finishthisgalaxy
           ENDIF
           IF size(log,/TYPE) EQ 7 then begin
              openu,66,log,/APPEND
              printf,66,linenumber()+'We are cutting the cube as the last channel is completely blank.'
              close,66
           ENDIF
           sxaddpar,hed,'NAXIS3',fix(sxpar(hed,'NAXIS3')-1)
           changedcube=1
           dummy=tmp
           tmpnoblank=dummy[*,*,n_elements(dummy[0,0,*])-1]
           wherefinite=WHERE(FINITE(tmpnoblank))
        ENDWHILE
        rmslastchannel=SIGMA(tmpnoblank[WHERE(FINITE(tmpnoblank))])
        IF rmsfirstchannel EQ 0. then rmsfirstchannel=2.*rmslastchannel
        IF rmslastchannel EQ 0. then rmslastchannel=2.*rmsfirstchannel
        tmpnoblank=dummy[0:5,0:5,*]
        rmsbottoml=SIGMA(tmpnoblank[WHERE(FINITE(tmpnoblank))])
        tmpnoblank=dummy[n_elements(dummy[*,0,0])-6:n_elements(dummy[*,0,0])-1,0:5,*]
        rmsbottomr=SIGMA(tmpnoblank[WHERE(FINITE(tmpnoblank))])
        tmpnoblank=dummy[n_elements(dummy[*,0,0])-6:n_elements(dummy[*,0,0])-1,n_elements(dummy[0,*,0])-6:n_elements(dummy[0,*,0])-1,*]
        rmstopr=SIGMA(tmpnoblank[WHERE(FINITE(tmpnoblank))])
        tmpnoblank=dummy[0:5,n_elements(dummy[0,*,0])-6:n_elements(dummy[0,*,0])-1,*]
        rmstopl=SIGMA(tmpnoblank[WHERE(FINITE(tmpnoblank))])
        rmschan=(rmsfirstchannel+rmslastchannel)/2.
        rmscorn=(rmsbottoml+rmsbottomr+rmstopl+rmstopr)/4.
        diff=ABS((rmsfirstchannel-rmslastchannel)/rmsfirstchannel)
        diff2=ABS((rmschan-rmscorn)/rmschan)
        IF diff LT 0.2 AND FINITE(diff) AND diff2 LT 0.2 then difference=0. else begin
           IF size(log,/TYPE) EQ 7 then begin
              openu,66,log,/APPEND
              printf,66,linenumber()+'We are cutting the cube as clearly the noise statistics are off.'
              printf,66,linenumber()+'Noise in the first channel is '+string(rmsfirstchannel)
              printf,66,linenumber()+'Noise in the last channel is '+string(rmslastchannel)
              printf,66,linenumber()+'Noise in the corners is '+string(rmscorn)
              close,66
           ENDIF ELSE BEGIN
              print,linenumber()+'We are cutting the cube as clearly the noise statistics are off.'
              print,linenumber()+'Noise in the first channel is '+string(rmsfirstchannel)
              print,linenumber()+'Noise in the last channel is '+string(rmslastchannel)       
           ENDELSE
           changedcube=1.
           tmp=fltarr(n_elements(dummy[*,0,0]),n_elements(dummy[0,*,0]),n_elements(dummy[0,0,*])-1)
           firstcomp=ABS((rmsfirstchannel-rmscorn)/rmscorn)
           secondcomp=ABS((rmslastchannel-rmscorn)/rmscorn)
           IF firstcomp GT secondcomp OR NOT FINITE(rmsfirstchannel) then begin
              tmp[*,*,0:n_elements(tmp[0,0,*])-1]=dummy[*,*,1:n_elements(dummy[0,0,*])-1]
              sxaddpar,hed,'CRPIX3',sxpar(hed,'CRPIX3')-1.
           ENDIF ELSE BEGIN
              tmp[*,*,0:n_elements(tmp[0,0,*])-1]=dummy[*,*,0:n_elements(dummy[0,0,*])-2]
           endelse
           sxaddpar,hed,'NAXIS3',fix(sxpar(hed,'NAXIS3')-1)
           dummy=fltarr(n_elements(tmp[*,0,0]),n_elements(tmp[0,*,0]),n_elements(tmp[0,0,*]))
           dummy=tmp
           IF n_elements(dummy[0,0,*]) LT 5 then begin
           IF size(log,/TYPE) EQ 7 then begin
              openu,66,log,/APPEND
              printf,66,linenumber()+catdirname[i]+'/'+catcubename[i]+' has noise statistics that cannot be dealt with.'
              close,66
           ENDIF ELSE BEGIN
              printf,linenumber()+catdirname[i]+'/'+catcubename[i]+' has noise statistics that cannot be dealt with.'
           ENDELSE         
           openu,1,outputcatalogue,/APPEND
           printf,1,format='(A60,A90)', catDirname[i],'The Cube has noise statistics that cannot be dealt with'
           close,1    
           bookkeeping=5
           goto,finishthisgalaxy
        ENDIF
        ENDELSE
     ENDWHILE
     catnoise[i]=rmscorn
     IF changedcube EQ 1. then begin        
        currentfitcube=catcubename[i]+'_cut'
        catcubename[i]=catcubename[i]+'_cut'      
        writefits,maindir+'/'+catdirname[i]+'/'+catcubename[i]+cubeext,dummy,hed
        IF size(log,/TYPE) EQ 7 then begin
           openu,66,log,/APPEND
           printf,66,linenumber()+'The cube was cut due to incorrect noise statistics.'
           printf,66,linenumber()+'The new cube is '+catcubename[i]+cubeext
           close,66
        ENDIF ELSE BEGIN
           print,linenumber()+'The cube was cut due to incorrect noise statistics.'
           print,linenumber()+'The new cube is '+catcubename[i]+cubeext      
        ENDELSE
        propermom0=-1
        propermom1=-1
     ENDIF
                                ; And last we want to cut out central absorption parts. We will cut
                                ; below negative 5 sigma. If this messes up the process there is
                                ; something seriously wrong with the cube                            
     tmp=WHERE(dummy LT -10.*catnoise[i])
     pythoncube=catcubename[i]
     possiblecentral=0
     IF tmp[0] NE -1 then begin
        extbl=fltarr(n_elements(dummy[*,0,0]),n_elements(dummy[0,*,0]),n_elements(dummy[0,0,*]))
        extbl=dummy
        extbl[tmp]=!values.f_nan
        pythoncube=catcubename[i]+'_extbl'
        writefits,maindir+'/'+catdirname[i]+'/'+pythoncube+cubeext,extbl,hed
        tmp=WHERE(dummy LT -50.*catnoise[i])
        IF tmp[0] NE -1 and wroteblank ne 1 then begin
           dummy[tmp]=!values.f_nan
           currentfitcube=catcubename[i]+'_bl'
           catcubename[i]=catcubename[i]+'_bl'      
           writefits,maindir+'/'+catdirname[i]+'/'+catcubename[i]+cubeext,dummy,hed
           wroteblank=1
        ENDIF   
        possiblecentral=1
        IF size(log,/TYPE) EQ 7 then begin
           openu,66,log,/APPEND
           printf,66,linenumber()+'Your cube had values below -10*sigma. If you do not have a central absorption source there is something seriously wrong with the cube.'
           
           close,66
        ENDIF ELSE BEGIN
           print,linenumber()+'Your cube had values below -10.*sigma. If you do not have a central absorption source there is something seriously wrong with the cube.'
           
           
        ENDELSE
     ENDIF
                                ;if the moment maps or the proper binmask do not yet exist we create
                                ;them here
     IF size(log,/TYPE) EQ 7 then begin
        openu,66,log,/APPEND
        printf,66,linenumber()+'The noise in the cube is '+string(catnoise[i])
        close,66
     ENDIF ELSE BEGIN
        print,linenumber()+'The noise in the cube is '+string(catnoise[i])
     ENDELSE   
     createmoment=0
                                ;if we have preprepared sofia files we will skip all the sofia steps
     IF allnew EQ 2 then goto,skipallsofia
     CD,new_dir,CURRENT=old_dir
     IF propermom0 NE -1 then propermom0=FILE_TEST(catMom0name[i]+'.fits')
     IF propermom0 EQ 0 AND FILE_TEST(catcubename[i]+'_6.0_mom0.fits') then begin
        propermom0=1
        catMom0name[i]=catcubename[i]+'_6.0_mom0'
        mom0ext='.fits'
        createmoment=1
     ENDIF
     IF propermom1 NE -1 then propermom1=FILE_TEST(catMom1name[i]+'.fits')
     IF propermom1 EQ 0 AND FILE_TEST(catcubename[i]+'_6.0_mom1.fits') then begin
        propermom1=1
        catMom1name[i]=catcubename[i]+'_6.0_mom1'
        mom1ext='.fits'
     ENDIF
     CD,old_dir
     mismatchedmoments:
     CD,new_dir,CURRENT=old_dir
     propermask=FILE_TEST(catcubename[i]+'_6.0_binmask.fits')
     sofia[sofiatriggers[0]]='import.inFile = '+pythoncube+'.fits'
     sofia[sofiatriggers[1]]='steps.doReliability             =       false'
                                ;We run sofia to get a bunch of initial estimates
     openw,1,'sofia_input.txt'
     for index=0,n_elements(sofia)-1 do begin
        printf,1,sofia[index]
     endfor
     close,1
     spawn,'rm -f '+pythoncube+'_cat.ascii '+pythoncube+'_cont.png '+pythoncube+'_mask.fits '+pythoncube+'_scat.png',isthere
     spawn,'python '+supportdirchecked+'/sofia_pipeline.py sofia_input.txt '
     IF not FILE_TEST(pythoncube+'_cat.ascii') then begin
        IF size(log,/TYPE) EQ 7 then begin
           openu,66,log,/APPEND
           printf,66,linenumber()+"Can't make a mask, aborting."
           close,66
        ENDIF
        openu,1,outputcatalogue,/APPEND
        printf,1,format='(A60,A90)', catDirname[i],'Cannot produce a mask'
        close,1 
        CD,old_dir
        bookkeeping=5
        goto,finishthisgalaxy
     endif
     catCatalogname[i]=pythoncube+'_cat.ascii'
     catmaskname[i]=pythoncube+'_mask' ;first clean up the cat
     skipallsofia:
     IF allnew EQ 2 then CD,new_dir,CURRENT=old_dir
                                ;We need to read the sofia output file     
     openr,1,catCatalogname[i]
                                ;get the number of lines
     paralines=FILE_LINES(catCatalogname[i])
                                ;define some variables
     counter=1
     line=''
     ID=0.
     trig=0.
     catvals=0.
                                ;Define the parameters that need to be present and that we want.
     sofia_parameters=['id','x','x_min','x_max','y','y_min','y_max','z','z_min','z_max','x_geo','y_geo','z_geo','f_int']
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
           print,'This can happen because a) you have tampered with the sofiainput.txt file in the directory '+supportdirchecked
           print,'or b) you are using an updated version of SoFiA where the names have changed and FAT is not yet updated.'
           print,'In this case please file a bug report at https://github.com/PeterKamphuis/FAT/issues/'
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
        voxratio=dblarr(n_elements(catvals[0,*]))
        maxvoxel=MAX(catvals[sofia_locations[13],*])
        for j=0,n_elements(catvals[0,*])-1 do begin
           voxratio[j]=maxvoxel/catvals[sofia_locations[13],j]
        endfor
        rmp=WHERE(voxratio LT 3)
        IF rmp[0] EQ -1 then begin
           vals=catvals[*,n_elements(catvals[0,*])-1]
        ENDIF ELSE BEGIN        
           diff=dblarr(n_elements(rmp))
           for j=0,n_elements(rmp)-1 do begin
              diff[j]=SQRT((sxpar(hed,'CRPIX1')-catvals[sofia_locations[1],rmp[j]])^2+(sxpar(hed,'CRPIX2')-catvals[sofia_locations[4],rmp[j]])^2)
           endfor
           maxdiff=MAX(diff,MIN=whatwelookfor)
           vals=catvals[*,rmp[WHERE(diff EQ whatwelookfor)]]
        ENDELSE
     ENDIF ELSE vals=catvals
                                ;If sofia didn't work and there is no preprocessed stuff then abort
     IF double(vals[sofia_locations[0]]) EQ 0. AND allnew NE 2 then begin
        IF size(log,/TYPE) EQ 7 then begin
           openu,66,log,/APPEND
           printf,66,linenumber()+"Sofia failed, the original cube is too small."
           close,66
        ENDIF
        openu,1,outputcatalogue,/APPEND
        printf,1,format='(A60,A90)', catDirname[i],'Sofia failed on this cube.'
        close,1 
        CD,old_dir
        bookkeeping=5
        goto,finishthisgalaxy
     ENDIF
     IF size(log,/TYPE) EQ 7 then begin
        openu,66,log,/APPEND
        printf,66,linenumber()+"We picked the "+string(vals[sofia_locations[0]])+" object of the parameter list."
        close,66
     ENDIF
                                ;If not using preprocessed stuff then we make logical names 
     IF allnew NE 2 then begin
        nummask=readfits(catmaskname[i]+'.fits',hedmasknotproper)
        tmp=where(nummask NE double(vals[sofia_locations[0]])) 
        IF tmp[0] NE -1 then begin
           nummask[tmp]=0.
        ENDIF
        IF vals[0] NE 1. then begin
           print,vals[sofia_locations[0]]
           tmp=where(nummask EQ double(vals[sofia_locations[0]]))
           nummask[tmp]=1.
        ENDIF
        writefits,pythoncube+'_6.0_binmask.fits',nummask,hed
        catmaskname[i]=pythoncube+'_6.0_binmask'
        spawn,'rm -f '+pythoncube+'_cont.png '+pythoncube+'_mask.fits '+pythoncube+'_scat.png',isthere
     ENDIF
                                ;Read the centre values and make sure
                                ;they fall inside the cube
     RApix=double(vals[sofia_locations[1]])
     RApixboun=[double(vals[sofia_locations[2]]),double(vals[sofia_locations[3]])]   
     IF RApix GT sxpar(hed,'NAXIS1')-10 then BEGIN
        RApix=double(sxpar(hed,'NAXIS1')/2.)
        RApixboun=[double(sxpar(hed,'NAXIS1')/4.),double(sxpar(hed,'NAXIS1')/4.)*3.]
     ENDIF
     DECpix=double(vals[sofia_locations[4]])
     DECpixboun=[double(vals[sofia_locations[5]]),double(vals[sofia_locations[6]])]
     IF DECpix GT sxpar(hed,'NAXIS2')-10 then begin
        DECpix=fix( sxpar(hed,'NAXIS2')/2.)
        DECpixboun=[double(sxpar(hed,'NAXIS2')/4.),double(sxpar(hed,'NAXIS2')/4.)*3.]
     ENDIF
     VSYSpix=double(vals[sofia_locations[7]])
     ROTpixboun=[double(vals[sofia_locations[8]]),double(vals[sofia_locations[9]])]
     IF VSYSpix GT sxpar(hed,'NAXIS3')-5 then BEGIN
        VSYSpix=fix( sxpar(hed,'NAXIS3')/2.)
        ROTpixboun=[double(sxpar(hed,'NAXIS3')/4.),double(sxpar(hed,'NAXIS3')/4.)*3.]
     ENDIF
     IF RApix LT RApixboun[0] OR RApix GT RApixboun[1] $
        OR DECpix LT DECpixboun[0] OR DECpix GT DECpixboun[1] $
        OR VSYSpix LT ROTpixboun[0] OR VSYSpix GT ROTpixboun[1] then begin
        RApix=double(vals[sofia_locations[10]])
        DECpix=double(vals[sofia_locations[11]])
        VSYSpix=double(vals[sofia_locations[12]])
     ENDIF
     Totflux=[double(vals[sofia_locations[13]])] ;Jy/beam   
                                ;let's check whether we have flux in this central pixel
   

                                ;Convert pixel coordinates to Degrees
     DECdeg=sxpar(hed,'CRVAL2')+(DECpix-sxpar(hed,'CRPIX2'))*sxpar(hed,'CDELT2')
     RAdeg=sxpar(hed,'CRVAL1')+(RApix-sxpar(hed,'CRPIX1'))*sxpar(hed,'CDELT1')/COS(DECdeg*(!pi/180.))
     DECboundeg=sxpar(hed,'CRVAL2')+(DECpixboun-sxpar(hed,'CRPIX2'))*sxpar(hed,'CDELT2')
     RAboundeg=sxpar(hed,'CRVAL1')+((RApixboun-sxpar(hed,'CRPIX1'))*sxpar(hed,'CDELT1'))/COS(DECboundeg*(!pi/180.))
     RAboundeg=RAboundeg[SORT(RAboundeg)]
     DECboundeg=DECboundeg[SORT(DECboundeg)]
     IF strupcase(veltype) EQ 'M/S' then begin
        catVSYS[i]=sxpar(hed,'CRVAL3')/1000.+(VSYSpix-sxpar(hed,'CRPIX3'))*sxpar(hed,'CDELT3')/1000.
        ROTboun=sxpar(hed,'CRVAL3')/1000.+(ROTpixboun-sxpar(hed,'CRPIX3'))*sxpar(hed,'CDELT3')/1000.
        ROTboun=ROTboun[SORT(ROTboun)]
     ENDIF else begin
        catVSYS[i]=sxpar(hed,'CRVAL3')+(VSYSpix-sxpar(hed,'CRPIX3'))*sxpar(hed,'CDELT3')
        ROTboun=sxpar(hed,'CRVAL3')+(ROTpixboun-sxpar(hed,'CRPIX3'))*sxpar(hed,'CDELT3')
        ROTboun=ROTboun[SORT(ROTboun)]
     ENDELSE
                                ;if we do not have a moment 0 map and
                                ;a moment one map then we are going to
                                ;make them now
     IF propermom0 LE 0 OR propermom1 LE 0 and allnew NE 2 THEN BEGIN
        print,linenumber()+'We are removing the moment maps.'
                                ;We are not remake merely one as the
                                ;velocity field and the moment 0 map
                                ;should be consistent
        spawn,'rm -f '+catcubename[i]+'_6.0_mom0.fits',isthere
        spawn,'rm -f '+catcubename[i]+'_6.0_mom1.fits',isthere
                                ;read the mask 
        mask=readfits(maindir+'/'+catdirname[i]+'/'+catmaskname[i]+'.fits',hedmask,/NOSCALE)        
                                ;mask the data cube
        tmpmask=fltarr(n_elements(dummy[*,0,0]),n_elements(dummy[0,*,0]),n_elements(dummy[0,0,*]))
        tmpmask[WHERE(mask GT 0.)]=dummy[WHERE(mask GT 0.)]
        hedmap=hed
                                ;Same as in gipsy
        momentsv2,tmpmask,moment0map,hedmap,0.
        writefits,catcubename[i]+'_6.0_mom0.fits',moment0map,hedmap
        hedvel=hed
        momentsv2,tmpmask,moment1map,hedvel,1.
        writefits,catcubename[i]+'_6.0_mom1.fits',moment1map,hedvel
        createmoment=1
                                ;reset the mom0 map to the one just created
        catMom0name[i]=catcubename[i]+'_6.0_mom0'
                                ;reset the mom1 map to the one just created
        catMom1name[i]=catcubename[i]+'_6.0_mom1'
                                ;the extension is fits
        mom0ext='.fits'
        mom1ext='.fits'
                                ;Make a noisemap as well
        IF veltype EQ 'M/S' then channelwidthkm=channelwidth/1000. else channelwidthkm=channelwidth
        noisemap=fltarr(n_elements(mask[*,0,0]),n_elements(mask[*,0,0]))
        for j=0,n_elements(mask[0,0,*])-1 do begin
           noisemap=noisemap+(mask[*,*,j]*catnoise[i]*channelwidthkm)^2
        endfor
        noisemap=SQRT(noisemap)
        writefits,catcubename[i]+'_6.0_noisemap.fits',noisemap,hedmap
        noisemapname=catCubename[i]+'_6.0_noisemap'
        tmp=WHERE(noisemap NE 0.)
        momnoise=TOTAL(noisemap[tmp])/n_elements(tmp)
        
                                ;Print the log what we have done
        IF size(log,/TYPE) EQ 7 then begin
           openu,66,log,/APPEND
           printf,66,linenumber()+"We have created the moment maps."
           close,66
        ENDIF ELSE BEGIN
           print,linenumber()+"We have created the moment maps."
        Endelse    
     ENDIF ELSE BEGIN
                                ; if we are using preprocessed stuff then read and adapt them
        IF allnew EQ 2 then begin
           mask=readfits(maindir+'/'+catdirname[i]+'/'+catmaskname[i]+'.fits',hedmask,/NOSCALE)
           tmp=WHERE(mask GT 1)
           IF tmp[0] NE -1 then begin
              tmp=WHERE(mask NE 0.)
              mask[tmp]=1.
           ENDIF
           noisemap=fltarr(n_elements(mask[*,0,0]),n_elements(mask[*,0,0]))
           IF veltype EQ 'M/S' then channelwidthkm=channelwidth/1000. else channelwidthkm=channelwidth
           for j=0,n_elements(mask[0,0,*])-1 do begin
              noisemap=noisemap+(mask[*,*,j]*catnoise[i]*channelwidthkm)^2
           endfor
           writefits,maindir+'/'+catdirname[i]+'/'+catcubename[i]+'_noisemap.fits',noisemap,hedmap
           noisemapname=catCubename[i]+'_noisemap'
        endif else begin
           noisemap=readfits(catcubename[i]+'_6.0_noisemap.fits',hednoisemap)
           noisemapname=catCubename[i]+'_6.0_noisemap'
           tmp=WHERE(noisemap NE 0.)
           momnoise=TOTAL(noisemap[tmp])/n_elements(tmp)
        ENDELSE
     ENDELSE
     CD,old_dir
                                ;let's check the size of this cube to see whether it is reasonable 
     IF maxrings EQ 0. then begin
        mask=readfits(maindir+'/'+catdirname[i]+'/'+catmaskname[i]+'.fits',hedmask,/NOSCALE)
        foundmin=0.
        minx=intarr(sxpar(hedmask,'NAXIS1'))
        maxx=intarr(sxpar(hedmask,'NAXIS1'))
        countmin=-1
        for pix=0,n_elements(mask[*,0,0])-1 do begin
           tmp=WHERE(mask[pix,*,*] GT 0.)
           IF tmp[0] NE -1 AND foundmin EQ 0. then begin
              countmin++
              minx[countmin]=pix
              foundmin=1
           ENDIF
           IF tmp[0] EQ -1 AND foundmin EQ 1. then begin
              maxx[countmin]=pix
              foundmin=0.
           ENDIF
        endfor
        IF maxx[countmin] EQ 0 then maxx[countmin]=n_elements(mask[*,0,0])-1
        foundmin=0.
        countmin=-1
        miny=intarr(sxpar(hedmask,'NAXIS2'))
        maxy=intarr(sxpar(hedmask,'NAXIS2'))
        for pix=0,n_elements(mask[0,*,0])-1 do begin
           tmp=WHERE(mask[*,pix,*] GT 0.)
           IF tmp[0] NE -1 AND foundmin EQ 0. then begin
              countmin++
              miny[countmin]=pix
              foundmin=1
           ENDIF
           IF tmp[0] EQ -1 AND foundmin EQ 1. then begin
              maxy[countmin]=pix
              foundmin=0.
           ENDIF
        endfor
        IF maxy[countmin] EQ 0 then maxy[countmin]=n_elements(mask[0,*,0])-1
        xhandle=MAX(maxx-minx)
        yhandle=MAX(maxy-miny)
        maxrings=round(SQRT((xhandle/2.)^2+(yhandle/2.)^2)/(catmajbeam[i]/(ABS(sxpar(hed,'cdelt1'))*3600.))+3.)
        IF xhandle[0] GT yhandle[0] then begin
           IF xhandle[0]/(catmajbeam[i]/(ABS(sxpar(hed,'cdelt1'))*3600.)) LE 6 then newsize=fix(round(xhandle[0]+(10.+sofiacount)*catmajbeam[i]/(ABS(sxpar(hed,'cdelt1'))*3600.))) else $
              newsize=fix(round(xhandle[0]+(8.+sofiacount)*catmajbeam[i]/(ABS(sxpar(hed,'cdelt1'))*3600.)))
        ENDIF ELSE BEGIN
           IF yhandle[0]/(catmajbeam[i]/(ABS(sxpar(hed,'cdelt2'))*3600.)) LE 6 then newsize=fix(round(yhandle[0]+(10.+sofiacount)*catmajbeam[i]/(ABS(sxpar(hed,'cdelt2'))*3600.))) else $
              newsize=fix(round(yhandle[0]+(8.+sofiacount)*catmajbeam[i]/(ABS(sxpar(hed,'cdelt2'))*3600.)))
        ENDELSE
        IF floor(RApix-newsize/2.) GE 0. AND floor(RApix+newsize/2.) LT  sxpar(hed,'NAXIS1') AND $
           floor(DECpix-newsize/2.) GE 0. AND floor(DECpix+newsize/2.) LT sxpar(hed,'NAXIS2') AND $
           not smallexists then begin
           newdummy=fltarr(newsize+1,newsize+1,n_elements(dummy[0,0,*]))
           newdummy[*,*,*]=dummy[floor(RApix-newsize/2.):floor(RApix+newsize/2.),floor(DECpix-newsize/2.):floor(DECpix+newsize/2.),*]
           IF size(log,/TYPE) EQ 7 then begin
              openu,66,log,/APPEND
              printf,66,linenumber()+"In order to speed up the fitting process significantly we cut the cube to "+string(newsize)+ " pixels." 
              close,66
           ENDIF          
           sxaddpar,hed,'NAXIS1',newsize
           sxaddpar,hed,'NAXIS2',newsize
           sxaddpar,hed,'CRPIX1',sxpar(hed,'CRPIX1')-floor(RApix-newsize/2.)
           sxaddpar,hed,'CRPIX2',sxpar(hed,'CRPIX2')-floor(DECpix-newsize/2.)
           writefits,maindir+'/'+catdirname[i]+'/'+catCubename[i]+'_small.fits',newdummy,hed
           dummy=newdummy
           newdummy[*,*,*]=0.
           newdummy[*,*,*]=mask[floor(RApix-newsize/2.):floor(RApix+newsize/2.),floor(DECpix-newsize/2.):floor(DECpix+newsize/2.),*]
           writefits,maindir+'/'+catdirname[i]+'/'+catmaskname[i]+'_small.fits',newdummy,hed
           mask=newdummy
           moment0map=readfits(maindir+'/'+catdirname[i]+'/'+catMom0Name[i]+'.fits',hedmap,/SILENT)
           newdummy=fltarr(newsize+1,newsize+1)
           sxaddpar,hedmap,'NAXIS1',newsize
           sxaddpar,hedmap,'NAXIS2',newsize
           sxaddpar,hedmap,'CRPIX1',sxpar(hedmap,'CRPIX1')-floor(RApix-newsize/2.)
           sxaddpar,hedmap,'CRPIX2',sxpar(hedmap,'CRPIX2')-floor(DECpix-newsize/2.)
           newdummy[*,*]=moment0map[floor(RApix-newsize/2.):floor(RApix+newsize/2.),floor(DECpix-newsize/2.):floor(DECpix+newsize/2.)]
           writefits,maindir+'/'+catdirname[i]+'/'+catMom0name[i]+'_small.fits',newdummy,hedmap
           moment1map=readfits(maindir+'/'+catdirname[i]+'/'+catMom1Name[i]+'.fits',hedvel,/SILENT)
           newdummy=fltarr(newsize+1,newsize+1)
           sxaddpar,hedvel,'NAXIS1',newsize
           sxaddpar,hedvel,'NAXIS2',newsize
           sxaddpar,hedvel,'CRPIX1',sxpar(hedvel,'CRPIX1')-floor(RApix-newsize/2.)
           sxaddpar,hedvel,'CRPIX2',sxpar(hedvel,'CRPIX2')-floor(DECpix-newsize/2.)
           newdummy[*,*]=moment1map[floor(RApix-newsize/2.):floor(RApix+newsize/2.),floor(DECpix-newsize/2.):floor(DECpix+newsize/2.)]
           writefits,maindir+'/'+catdirname[i]+'/'+catMom1name[i]+'_small.fits',newdummy,hedvel
           noisemap=readfits(maindir+'/'+catdirname[i]+'/'+catCubename[i]+'_6.0_noisemap.fits',hednoise,/SILENT)
           newdummy=fltarr(newsize+1,newsize+1)
           sxaddpar,hednoise,'NAXIS1',newsize
           sxaddpar,hednoise,'NAXIS2',newsize
           sxaddpar,hednoise,'CRPIX1',sxpar(hednoise,'CRPIX1')-floor(RApix-newsize/2.)
           sxaddpar,hednoise,'CRPIX2',sxpar(hednoise,'CRPIX2')-floor(DECpix-newsize/2.)
           newdummy[*,*]=noisemap[floor(RApix-newsize/2.):floor(RApix+newsize/2.),floor(DECpix-newsize/2.):floor(DECpix+newsize/2.)]
           writefits,maindir+'/'+catdirname[i]+'/'+catCubename[i]+'_6.0_noisemap_small.fits',newdummy,hednoise
           noisemapname=catCubename[i]+'_6.0_noisemap_small'
           cubeext='.fits'
           CD, new_dir, CURRENT=old_dir
           spawn,'rm -f '+pythoncube+'_cont.png '+pythoncube+'_scat.png',isthere
           IF allnew NE 2 then begin
              spawn,'rm -f '+catmaskname[i]+'.fits '+catMom1name[i]+'.fits '+catMom0name[i]+'.fits '+catCubename[i]+'_6.0_noisemap.fits',isthere
           ENDIF
           CD,old_dir
           currentfitcube=catcubename[i]+'_small'
           catCubename[i]=catCubename[i]+'_small'
           smallexists=1
           catMom0name[i]=catMom0name[i]+'_small'
           catMom1name[i]=catMom1name[i]+'_small'
           catmaskname[i]=catmaskname[i]+'_small'
           RApix=floor(newsize/2.)+(Rapix-floor(RApix))
           DECpix=floor(newsize/2.)+(DECpix-floor(DECpix))
        ENDIF
     ENDIF
     centralflux=dummy[fix(RApix),fix(DECpix),fix(VSYSpix)]
     IF size(log,/TYPE) EQ 7 then begin
        openu,66,log,/APPEND
        printf,66,linenumber()+"We have this value in the central pixel "+string(centralflux)+" at the pixel (x,y,z) "+string(fix(RApix))+","+string(fix(DECpix))+","+string(fix(VSYSpix))
        IF not FINITE(centralflux) then           printf,66,linenumber()+"Central exclude = 1" else   printf,66,linenumber()+"Central exclude = 0"
        close,66
     ENDIF
     IF not FINITE(centralflux) then centralexclude=1 else centralexclude=0.
                                ; Print the source finder results
     IF size(log,/TYPE) EQ 7 then begin
        openu,66,log,/APPEND
        printf,66,linenumber()+"The source finder found the following center in pixels."
        printf,66,linenumber()+"DEC center="+string(DECpix)+"DEC Boundaries"+string(DECpixboun[0])+","+string(DECpixboun[1])
        printf,66,linenumber()+"RA center="+string(RApix)+"RA Boundaries"+string(RApixboun[0])+","+string(RApixboun[1])
        printf,66,linenumber()+"V_sys center="+string(VSYSpix)+"V_sys Boundaries"+string(ROTpixboun[0])+","+string(ROTpixboun[1])
        printf,66,linenumber()+"Which converts to real coordinates as."      
        printf,66,linenumber()+"DEC center="+string(DECdeg)+"DEC Boundaries"+string(DECboundeg[0])+","+string(DECboundeg[1])
        printf,66,linenumber()+"RA center="+string(RAdeg)+"RA Boundaries"+string(RAboundeg[0])+","+string(RAboundeg[1])
        printf,66,linenumber()+"V_sys center="+string(catVSYS[i])+"V_sys Boundaries"+string(ROTboun[0])+","+string(ROTboun[1])
        close,66
     ENDIF     
                                ;if we have not created moment0 and
                                ;moment1 then we want to make sure we
                                ;have a proper extension.
     IF NOT createmoment then begin
        mom0ext=' '
        CD,new_dir,CURRENT=old_dir
        spawn,'ls '+catMom0name[i]+'.*',filenameext
        CD,OLD_DIR
        IF (filenameext[0]) then begin
           for j=0,n_elements(filenameext)-1 do begin
              tmp=str_sep(strtrim(strcompress(filenameext[j]),2),'.')
              for g=n_elements(tmp)-1,n_elements(tmp)-3,-1 do begin
                 
                 IF STRUPCASE(tmp[g]) EQ 'FITS' then begin
                    case g of
                       n_elements(tmp)-1:begin
                          nametmp=tmp[0]
                          for z=1,n_elements(tmp)-2 do begin
                             nametmp=nametmp+'.'+tmp[z]
                          endfor
                          
                          IF nametmp EQ catMom0name[i] then begin
                             mom0ext='.'+tmp[g]
                             goto,extmom0set
                          ENDIF
                       end
                       n_elements(tmp)-2:begin
                          IF STRUPCASE(tmp[g+1]) EQ 'GZ' then begin
                             nametmp=tmp[0]
                             for z=1,n_elements(tmp)-3 do begin
                                nametmp=nametmp+'.'+tmp[z]
                             endfor
                             
                             IF nametmp EQ catMom0name[i] then begin
                                mom0ext='.'+tmp[g]+'.'+tmp[g+1]
                                goto,extmom0set
                             ENDIF
                             
                          ENDIF
                       END
                    endcase
                 endif
              endfor
           endfor
           extmom0set:
        endif
        mom1ext=' '
        CD,new_dir,CURRENT=old_dir
        spawn,'ls '+catMom1name[i]+'.*',filenameext
        CD,OLD_DIR
        IF (filenameext[0]) then begin
           for j=0,n_elements(filenameext)-1 do begin
              tmp=str_sep(strtrim(strcompress(filenameext[j]),2),'.')
              for g=n_elements(tmp)-1,n_elements(tmp)-3,-1 do begin
                 IF STRUPCASE(tmp[g]) EQ 'FITS' then begin
                    case g of
                       n_elements(tmp)-1:begin
                          nametmp=tmp[0]
                          for z=1,n_elements(tmp)-2 do begin
                             nametmp=nametmp+'.'+tmp[z]
                          endfor
                          IF nametmp EQ catMom1name[i] then begin
                             mom1ext='.'+tmp[g]
                             goto,extmom1set
                          ENDIF
                       end
                       n_elements(tmp)-2:begin
                          IF STRUPCASE(tmp[g+1]) EQ 'GZ' then begin
                             nametmp=tmp[0]
                             for z=1,n_elements(tmp)-3 do begin
                                nametmp=nametmp+'.'+tmp[z]
                             endfor
                             IF nametmp EQ catMom1name[i] then begin
                                mom1ext='.'+tmp[g]+'.'+tmp[g+1]
                                goto,extmom1set
                             ENDIF
                          ENDIF
                       END
                    endcase
                 endif
              endfor
           endfor
        endif
        extmom1set:     
     ENDIF
                                ;We check whether the galaxy is bright
                                ;enough
     tmp=WHERE(mask EQ 1.)
     maxbright=MAX(dummy[tmp])
     maxSN=maxbright/catnoise[i]
     IF maxSN LT 5. then begin
        openu,1,outputcatalogue,/APPEND
        printf,1,format='(A60,A90)', catDirname[i],' The maximum Signal to Noise in this cube is '+string(MaxSN)+' that is not enough for a fit.'
        close,1    
        IF size(log,/TYPE) EQ 7 then begin
           openu,66,log,/APPEND
           printf,66,linenumber()+catDirname[i]+' The maximum Signal to Noise in this cube is '+string(MaxSN)+' that is not enough for a fit.'
           close,66
        ENDIF ELSE BEGIN
           print,linenumber()+catDirname[i]+' The maximum Signal to Noise in this cube is '+string(MaxSN)+' that is not enough for a fit.'
        ENDELSE
        bookkeeping=5 
        goto,finishthisgalaxy
     ENDIF 
     pixelarea=ABS(3600.^2*sxpar(hed,'CDELT2')*sxpar(hed,'CDELT1'))
     pixelsizeRA=ABS(sxpar(hed,'CDELT1'))
     pixelsizeDEC=ABS(sxpar(hed,'CDELT2'))
                                ;Let's get the size of the cube
     imagesize=sxpar(hed,'CDELT2')*(sxpar(hed,'NAXIS2')/2.)*3600.
                                ;IF it is not in km/s convert
     IF strupcase(veltype) EQ 'M/S' then channelwidth=channelwidth/1000.
   
     beamarea=(!pi*ABS(double(catminbeam[i])*double(catmajbeam[i])))/(4.*ALOG(2.))
     pixperbeam=beamarea/(ABS(pixelsizeRA*3600.)*ABS(pixelsizeDEC*3600.))
                                ;Obtain the cutoff values without
                                ;inclination corrections.
     rad=[0.,((findgen(maxrings+2.))*catmajbeam[i]+catmajbeam[i]/5.)]
     calc_edge,catnoise[i],rad,[catmajbeam[i],catminbeam[i]],cutoffor
                                ;And print the values
   
     IF size(log,/TYPE) EQ 7 then begin
        openu,66,log,/APPEND
        printf,66,linenumber()+"We use an intrinsic cutoff value of: "
        printf,66,linenumber()+"Radius           Minimum SBR"
        for j=0,n_elements(rad)-1 do begin
           printf,66,rad[j],cutoffor[j]
        endfor
        close,66
     ENDIF ELSE begin
        print,linenumber()+"We calculated the cutoff  values. "
     endelse
                                ;create the rings to be fitted with major beam
                                ;Let's calculate how many rings
                                ;we need to investigate the full HI
                                ;image. But let's reduce it by 1 rings
                                ;as to not run into edge problems 
     moment0map=readfits(maindir+'/'+catdirname[i]+'/'+catMom0name[i]+mom0ext,hedmap,/NOSCALE,/SILENT) 
     tmpfinite=WHERE(FINITE(moment0map) EQ 0)
     if tmpfinite[0] NE -1 then moment0map[tmpfinite]=0.
     totflux=TOTAL(moment0map)
    If NOT createmoment then begin
        IF sxpar(hedmap,'CDELT1') NE sxpar(hed,'CDELT1') OR sxpar(hedmap,'CDELT2') NE sxpar(hed,'CDELT2') then begin
           IF size(log,/TYPE) EQ 7 then begin
              openu,66,log,/APPEND
              printf,66,linenumber()+"the pixel scale in your cube and moment 0 map do not correspond."
              printf,66,linenumber()+"recreating the moment maps."
              close,66
           ENDIF 
           propermom0=0.
           propermom1=0.
           goto,mismatchedmoments
        ENDIF
        IF NOT sxpar(hedmap,'BMAJ') then sxaddpar,hedmap,'BMAJ',catmajbeam[i]
        IF NOT sxpar(hedmap,'BMIN') then sxaddpar,hedmap,'BMIN',catminbeam[i]
     ENDIF
                                ;Initial size guess is
     norings=fix(round((imagesize)/(catmajbeam[i]))-2)
                                ; If there is lt 1.5 beams in the cube we can stop here                     
     IF norings LE 4 then begin    
        IF norings LT 1.5 then begin
           openu,1,outputcatalogue,/APPEND
           printf,1,format='(A60,A80)', catDirname[i],'This Cube is too small'
           close,1
           bookkeeping=5
           goto,finishthisgalaxy
        ENDIF                   
     ENDIF
                       ;Now let's obtain a initial pa and inclination for the galaxy    
     obtain_pav2,moment0map,newPA,inclination=newinclination,center=[RApix,DECpix],noise=momnoise,iterations=10.
     moment1map=readfits(maindir+'/'+catdirname[i]+'/'+catMom1name[i]+mom1ext,hedvel,/NOSCALE,/SILENT) 
                                ;Let's check wether the PA
                                ;lines up with the kinematical
                                ;or needs a 180  degree shift
     dummyvel=moment1map
     tmpfinite=WHERE(FINITE(dummyvel) EQ 0)
     if tmpfinite[0] GT -1 then dummyvel[tmpfinite]=0.
     twobeampixel=3.*catmajbeam[i]/(ABS(sxpar(hedmap,'cdelt1'))*3600.)
     case (1) of
        (sxpar(hedvel,'cdelt1') LT 0):begin
           elowRA=fix(RApix-twobeampixel)
           IF elowRA LT 0 then elowRA=0
           ehighRA=fix(RApix)
           IF ehighRA GT n_elements(dummyvel[*,0])-1 then ehighRA=n_elements(dummyvel[*,0])-1
           elowDEC=fix(DECpix-twobeampixel)
           IF elowDEC LT 0 then elowDEC=0
           ehighDEC=fix(DECpix+twobeampixel)
           IF ehighDEC GT n_elements(dummyvel[0,*])-1 then ehighDEC=n_elements(dummyvel[0,*])-1
           IF elowRA GE ehighRA then elowRA=0
           IF elowDEC GE ehighDEC then elowDEC=0
           wlowRA=fix(RApix)
           IF wlowRA LT 0 then wlowRA=0
           whighRA=fix(RApix+twobeampixel)
           IF whighRA GT n_elements(dummyvel[*,0])-1 then whighRA=n_elements(dummyvel[*,0])-1
           wlowDEC=fix(DECpix-twobeampixel)
           IF wlowDEC LT 0 then wlowDEC=0
           whighDEC=fix(DECpix+twobeampixel)
           IF whighDEC GT n_elements(dummyvel[0,*])-1 then whighDEC=n_elements(dummyvel[0,*])-1
           IF wlowRA GE whighRA then whighRA=n_elements(dummyvel[*,0])-1
           IF wlowDEC GE whighDEC then whighDEC=n_elements(dummyvel[0,*])-1
        end
        (sxpar(hedvel,'cdelt1') GT 0):begin
           elowRA=fix(RApix)
           IF elowRA LT 0 then elowRA=0
           ehighRA=fix(RApix+twobeampixel)
           IF ehighRA GT n_elements(dummyvel[*,0])-1 then ehighRA=n_elements(dummyvel[*,0])-1
           elowDEC=fix(DECpix-twobeampixel)
           IF elowDEC LT 0 then elowDEC=0
           ehighDEC=fix(DECpix+twobeampixel)
           IF ehighDEC GT n_elements(dummyvel[0,*])-1 then ehighDEC=n_elements(dummyvel[0,*])-1
           IF elowRA GE ehighRA then ehighRA=n_elements(dummyvel[*,0])-1
           IF elowDEC GE ehighDEC then ehighDEC=n_elements(dummyvel[0,*])-1
           wlowRA=fix(RApix-twobeampixel)
           IF wlowRA LT 0 then wlowRA=0
           whighRA=fix(RApix)
           IF whighRA GT n_elements(dummyvel[*,0])-1 then whighRA=n_elements(dummyvel[*,0])-1
           wlowDEC=fix(DECpix-twobeampixel)
           IF wlowDEC LT 0 then wlowDEC=0
           whighDEC=fix(DECpix+twobeampixel)
           IF whighDEC GT n_elements(dummyvel[0,*])-1 then whighDEC=n_elements(dummyvel[0,*])-1
           IF wlowRA GE whighRA then wlowRA=0
           IF wlowDEC GE whighDEC then wlowDEC=0
        end
        else:begin
           elowRA=0
           elowDEC=0
           ehighDEC=n_elements(dummyvel[*,0])-1
           ehighRA=n_elements(dummyvel[0,*])-1
           wlowRA=0
           wlowDEC=0
           whighDEC=n_elements(dummyvel[*,0])-1
           whighRA=n_elements(dummyvel[0,*])-1
        end
     endcase
     IF elowRA gt ehighRA then elowRA=0.
     IF elowDEC gt ehighDEC then elowDEC=0.
     tmpvelmap=dummyvel[elowRA:ehighRA,elowDEC:ehighDEC]
     nozervel=WHERE(tmpvelmap NE 0)
     East=TOTAL(tmpvelmap[nozervel])/n_elements(nozervel)
     tmpvelmap=dummyvel[wlowRA:whighRA,wlowDEC:whighDEC]
     nozervel=WHERE(tmpvelmap NE 0)
     West=TOTAL(tmpvelmap[nozervel])/n_elements(nozervel)
     IF ABS(EAST-WEST) GT channelwidth*(1.+vresolution) then begin
        IF East GT West then begin
           IF newPA[0] GT 180. then begin
              newPA[0]=newPA[0]-180.
              IF size(log,/TYPE) EQ 7 then begin
                 openu,66,log,/APPEND
                 printf,66,linenumber()+"We flipped the PA because the East side had higher velocties East="+string(East)+' West='+string(West)
                 close,66
              ENDIF 
           ENDIF
           IF newPA[0] LT 170 and newPA[0] GT 10 and ABS(east-west) GT 20. then PAfixboun='e' else PAfixboun='n'
        endif else begin
           IF newPA[0] lT 180. then begin
              newPA[0]=newPA[0]+180.
              IF size(log,/TYPE) EQ 7 then begin
                 openu,66,log,/APPEND
                 printf,66,linenumber()+"We flipped the PA because the West side had higher velocties East="+string(East)+' West='+string(West) 
                 close,66
              ENDIF 
              
           ENDIF
           IF newPA[0] LT 350 and newPA[0] GT 190 AND ABS(east-west) GT 20. then PAfixboun='w'  else PAfixboun='n'
        endelse
     ENDIF ELSE BEGIN
        newPA[1]=360.
        PAfixboun='n'
     ENDELSE
     avPA=dblarr(2)
                                ;get the kinematical PA
     obtain_velpa,moment1map,velPA,center=[RApix,DECpix]
     IF ABS(velPA[0]-newPA[0]) LT 25. then begin
        avPA[0]=(velPA[0]+newPA[0])/2. 
        avPA[1]=SQRT(velPA[1]^2+newPA[1]^2)/2. 
     endif else begin
        IF FINITE(velPA[0]) then begin
           avPA[0]=velPA[0]
           avPA[1]=velPA[1]
        ENDIF else avPA=newPA
     ENDELSE
                                ;recalculate the inclination
    
     obtain_inclinationv8,moment0map,avPA,newinclination,[RApix,DECpix],extend=noringspix,noise=momnoise,beam=catmajbeam[i]/(pixelsizeRA*3600.)
     obtain_w50,dummy,mask,hed,W50   
                                ;We want to adjust the cutoff values
                                ;with the inclination as edge-on
                                ;galaxies are naturally integrated
                                ;The test galaxy in Jozsa 2007 has an
                                ;inclination of 75 deg, hence above
                                ;that the cutoff valu should be
                                ;declined and below that increased
     cutoffcorrection=SIN(75.*!DtoR)/SIN(newinclination[0]*!DtoR)
     IF newinclination[0] GT 60 then cutoffcorrection=1.
     IF newinclination[0] LT 50 AND newinclination[0] GT 40  then cutoffcorrection=1.+(50-newinclination[0])*0.05
     IF cutoffcorrection GT 2.5 then cutoffcorrection=2.5
     ringcorrection=cutoffcorrection
 ;    ringcorrection=0.
     WHILE maxrings-4-floor(ringcorrection/2.) LT 4 AND ringcorrection GT 1.03 do ringcorrection=ringcorrection-0.25
    
     IF size(log,/TYPE) EQ 7 then begin
        openu,66,log,/APPEND
        printf,66,linenumber()+"Cutoff values will be adjusted for inclination by multiplying them with"+string(cutoffcorrection)
        printf,66,linenumber()+"And the rings are reduced by "+string(floor(ringcorrection/2.))
        close,66
     ENDIF
     cutoff=cutoffor*cutoffcorrection
   
                                ;Also our initial ring estimates
                                ;should take this inclination into
                                ;account
     IF sofiafail then begin
        case 1 of
           (cutoffcorrection GT 1.03):norings=fix(round((noringspix*ABS(sxpar(hedmap,'cdelt1'))*3600.)/catmajbeam[i]))-floor(ringcorrection/2.)
           (cutoffcorrection LT 0.981):norings=fix(round((noringspix*ABS(sxpar(hedmap,'cdelt1'))*3600.)/catmajbeam[i])+1)
           else:norings=fix(round((noringspix*ABS(sxpar(hedmap,'cdelt1'))*3600.)/catmajbeam[i]))
        endcase
     endif else begin
        case 1 of
           cutoffcorrection GT 1.03 :norings=maxrings-4-floor(ringcorrection/2.)
           cutoffcorrection LT 0.981:norings=maxrings-4+1          
           else:norings=maxrings-4
        endcase
        IF newinclination[0] LT 60 then norings=norings[0]-1.
        noringspix=norings[0]*catmajbeam[i]/(ABS(sxpar(hedmap,'cdelt1'))*3600.)
     ENDELSE
     IF size(log,/TYPE) EQ 7 then begin
        openu,66,log,/APPEND
        printf,66,linenumber()+'We find a W50 of'+string(W50)
        printf,66,linenumber()+'We find an extend of '+strtrim(string(fix(norings)),2)+' at this point from the inclination estimates.'
        close,66
     ENDIF
     IF size(log,/TYPE) EQ 7 then begin
        openu,66,log,/APPEND
        printf,66,linenumber()+"The original estimates are:"
        printf,66,linenumber()+"maxrings:"+string(fix(maxrings))
        printf,66,linenumber()+"norings:"+strtrim(string(fix(norings)),2)
        printf,66,linenumber()+"Kinematical PA:"+string(velPA[0])+'+/-'+string(velPA[1])
        printf,66,linenumber()+"Morphological PA:"+string(newPA[0])+'+/-'+string(newPA[1])
        printf,66,linenumber()+"noringspix:"+string(noringspix)
        printf,66,linenumber()+"newinclination:"+string(newinclination[0])+'+/-'+string(newinclination[1])
        printf,66,linenumber()+"PA:"+string(avPA[0])+'+/-'+string(avPA[1])
        close,66
     ENDIF 
     newPA=avPA
     setfinishafter=0.
                                ; if the galaxy is too small we will
                                ; only fit a flat disk
     if norings[0] LE 3. then begin
        norings[0]=(maxrings-3)*2.
        noringspix=norings[0]*(catmajbeam[i]/2.)/(ABS(sxpar(hedmap,'cdelt1'))*3600.)
        rad=[0.,((findgen(maxrings*2.))*(catmajbeam[i]/2)+catmajbeam[i]/10.)]
        calc_edge,catnoise[i],rad,[catmajbeam[i],catminbeam[i]],cutoffor
        cutoffor=cutoffor/norings[0]
        cutoff=cutoffor*cutoffcorrection
        IF norings LT 3. then begin
           openu,1,outputcatalogue,/APPEND
           printf,1,format='(A60,A80)', catDirname[i],'This should not happen'
           close,1
           IF size(log,/TYPE) EQ 7 then begin
              openu,66,log,/APPEND
              printf,66,linenumber()+"Finished "+catDirname[i]+"in galaxy # "+strtrim(string(fix(i)),2)+" at "+systime()
              close,66
           ENDIF 
           finishafter=0.
           setfinishafter=1
        ENDIF ELSE begin
                                ;we don't want to do individual ring fits on this cube
           IF finishafter GT 1. then finishafter=1.1
           maxrings=maxrings*2.
                                ;we also want to make sure that we are in the cube    
           IF size(log,/TYPE) EQ 7 then begin
              openu,66,log,/APPEND
              printf,66,linenumber()+"Being a small galaxy we have halved the rings sizes "+catDirname[i]+"in galaxy # "+strtrim(string(fix(i)),2)+" at "+systime()
              printf,66,linenumber()+"The new cutoff values are: "
              printf,66,linenumber()+"Radius           Minimum SBR"
              for j=0,n_elements(rad)-1 do begin
                 printf,66,rad[j],cutoff[j]
              endfor
              close,66
           ENDIF   
        ENDELSE
     ENDIF ELSE BEGIN
        if norings[0] GT maxrings then begin
           norings[0]=maxrings
           noringspix=norings[0]*catmajbeam[i]/(ABS(sxpar(hedmap,'cdelt1'))*3600.)
        ENDIF
       
           
     ENDELSE
     
     IF finishafter EQ 1.1 then begin
        rings=(findgen(norings[0]))*(catmajbeam[i]/2.)+(catmajbeam[i]/10.)
        norings[0]=norings[0]+1.
        IF size(log,/TYPE) EQ 7 then begin
           openu,66,log,/APPEND
           printf,66,linenumber()+"The number of rings for the fitting = "+strtrim(string(fix(norings[0])),2)
           printf,66,linenumber()+"Fitting with half the beam major axis FWHM as ringsize."
           close,66
        ENDIF ELSE BEGIN
           print,linenumber()+"The number of rings for the fitting = "+strtrim(string(fix(norings[0])),2)
           print,linenumber()+"Fitting with half the beam major axis FWHM as ringsize."
        Endelse  
     ENDIF ELSE BEGIN
        IF maxrings GT 25 then begin
           tmpring=norings[0]-10.
           norings=10.+fix(tmpring/2.)
           rings=dblarr(norings[0])
           rings[0:9]=(findgen(10))*catmajbeam[i]+catmajbeam[i]/5.
           rings[10:norings[0]-1]=(findgen(fix(tmpring/2.)))*catmajbeam[i]*2+catmajbeam[i]/5.+11.*catmajbeam[i]
           norings[0]=norings[0]+1.
           rad=[0.,rings[0:9],(findgen(fix((maxrings-10.)/2.)))*catmajbeam[i]*2+catmajbeam[i]/5.+11.*catmajbeam[i]]
           finishafter=2.1
           maxrings=10.+fix((maxrings-10.)/2.)
           calc_edge,catnoise[i],rad,[catmajbeam[i],catminbeam[i]],cutoffor
           cutoff=cutoffor*cutoffcorrection
           IF size(log,/TYPE) EQ 7 then begin
              openu,66,log,/APPEND
              printf,66,linenumber()+"Being a large galaxy we have doubled the outer the ring sizes "+catDirname[i]+"in galaxy # "+strtrim(string(fix(i)),2)+" at "+systime()
              printf,66,linenumber()+"The new cutoff values are: "
              printf,66,linenumber()+"Radius           Minimum SBR"
              for j=0,n_elements(rad)-1 do begin
                 printf,66,rad[j],cutoff[j]
              endfor
              close,66
           ENDIF   
           IF size(log,/TYPE) EQ 7 then begin
              openu,66,log,/APPEND
              printf,66,linenumber()+"The number of rings for the fitting = "+strtrim(string(fix(norings[0])),2)
              printf,66,linenumber()+"Fitting with the twice the beam major axis FWHM beyond ring 10."
              close,66
           ENDIF ELSE BEGIN
              print,linenumber()+"The number of rings for the fitting = "+strtrim(string(fix(norings[0])),2)
              print,linenumber()+"Fitting with the twice the beam major axis FWHM beyond ring 10."
           Endelse  


        ENDIF ELSE BEGIN
           rings=(findgen(norings[0]))*catmajbeam[i]+catmajbeam[i]/5.
           norings[0]=norings[0]+1.
           IF size(log,/TYPE) EQ 7 then begin
              openu,66,log,/APPEND
              printf,66,linenumber()+"The number of rings for the fitting = "+strtrim(string(fix(norings[0])),2)
              printf,66,linenumber()+"Fitting with the beam major axis FWHM as ringsize."
              close,66
           ENDIF ELSE BEGIN
              print,linenumber()+"The number of rings for the fitting = "+strtrim(string(fix(norings[0])),2)
              print,linenumber()+"Fitting with the beam major axis FWHM as ringsize."
           Endelse  
        ENDELSE
     ENDELSE



                                ;since we have velocities we need to
                                ;check which half our PA should be
                                ;IF negative rotation values are in
                                ;the east the PA should be between 180
                                ;to 360 else between 0-180 
                                ;This is so that the rotation curve is always positive
     tmpfinite=WHERE(FINITE(dummyvel) EQ 0)
     if tmpfinite[0] GT -1 then dummyvel[tmpfinite]=0.
                                ;If the moment 0 is provided check that
                                ;it has the same properties as the cube
     If NOT createmoment then begin
        IF sxpar(hedvel,'CDELT1') NE sxpar(hed,'CDELT1') OR sxpar(hedvel,'CDELT2') NE sxpar(hed,'CDELT2') then begin
           IF size(log,/TYPE) EQ 7 then begin
              openu,66,log,/APPEND
              printf,66,linenumber()+"the pixel scale in your cube and moment 0 map do not correspond."
              printf,66,linenumber()+"recreating the moment maps."
              close,66
           ENDIF 
           propermom0=0.
           propermom1=0.
           goto,mismatchedmoments
        ENDIF
        IF NOT sxpar(hedvel,'BMAJ') then sxaddpar,hedvel,'BMAJ',catmajbeam[i]
        IF NOT sxpar(hedvel,'BMIN') then sxaddpar,hedvel,'BMIN',catminbeam[i]
     ENDIF
                                ;print the inclination and PA we found to the output log
     IF size(log,/TYPE) EQ 7 then begin
        openu,66,log,/APPEND
        printf,66,linenumber()+"We are using an inclination of "+string(newinclination[0])+" +/- "+string(newinclination[1])
        printf,66,linenumber()+"We are using an PA of "+string(newPA[0])+" +/- "+string(newPA[1])      
        close,66
     ENDIF 
                                ; store the original estimates to compare against later on
     IF catPA[i] NE 0 AND catinc[i] NE 0 then begin
        UndatedPA=catPA[i]
        UndatedINCL=catinc[i]
        UndatedINCLdev=catincdev[i]
        Undatedmaxrot=catmaxrot[i]
     ENDIF ELSE BEGIN
        UndatedPA=newPA[0]
        UndatedINCL=newinclination[0]
        UndatedINCLdev=newinclination[1]
        IF catinc[i] LT 40 then begin
           Undatedmaxrot=W50/2./SIN(ABS(newinclination[0]-2.5)*!pi/180.)
        ENDIF ELSE Undatedmaxrot=W50/2./SIN(newinclination[0]*!pi/180.)
        
     ENDELSE
                                ;Write the new estimates as catalogue value
     catPA[i]=newPA[0]
     catPAdev[i]=newPA[1]
     catinc[i]=newinclination[0]
    
     catincdev[i]=newinclination[1]
                                ;calculate the vmax based on inclination and W50
     IF catinc[i] LT 40 then begin
        catmaxrot[i]=W50/2./SIN(ABS(newinclination[0]+5.)*!pi/180.)
     ENDIF ELSE catmaxrot[i]=W50/2./SIN(newinclination[0]*!pi/180.)
                                ;Let's optimize the cube for fitting 
     IF catminbeam[i]/(cubecdelt*3600.) GT optpixelbeam  then begin        
        newcdelt=catminbeam[i]/double(optpixelbeam)
        ratio=double(newcdelt/double((cubecdelt*3600.)))
        new_hed=hed
        new_Cube=CONGRID(dummy,fix(sxpar(hed,'NAXIS1')/ratio),fix(sxpar(hed,'NAXIS2')/ratio),sxpar(hed,'NAXIS3'),/CENTER,/MINUS_ONE)
        sxaddpar,new_hed,'CDELT1',sxpar(hed,'CDELT1')*(double(sxpar(hed,'NAXIS1'))/(fix(sxpar(hed,'NAXIS1')/ratio)-1)),format='(E20.12)'
        sxaddpar,new_hed,'CDELT2',sxpar(hed,'CDELT2')*(double(sxpar(hed,'NAXIS2'))/(fix(sxpar(hed,'NAXIS2')/ratio)-1)),format='(E20.12)'
        sxaddpar,new_hed,'CRPIX1',sxpar(hed,'CRPIX1')/(double(sxpar(hed,'NAXIS1'))/(fix(sxpar(hed,'NAXIS1')/ratio))),format='(E20.12)'
        sxaddpar,new_hed,'CRPIX2',sxpar(hed,'CRPIX2')/(double(sxpar(hed,'NAXIS2'))/(fix(sxpar(hed,'NAXIS2')/ratio))),format='(E20.12)'
        sxaddpar,new_hed,'NAXIS1',fix(sxpar(hed,'NAXIS1')/ratio)
        sxaddpar,new_hed,'NAXIS2',fix(sxpar(hed,'NAXIS2')/ratio)
        CD,new_dir,CURRENT=old_dir
        currentfitcube=catcubename[i]+'_opt'
        catcubename[i]=catcubename[i]+'_opt'      
        writefits,maindir+'/'+catdirname[i]+'/'+catcubename[i]+cubeext,new_Cube,new_hed
        cd,old_dir
        dummy=new_Cube
        hed=new_hed
        optimized=1
     ENDIF
                                ;Let's write a file with the basic initial parameters
     oldRA=RAdeg
     oldDEC=DECdeg
     oldVSYS=catVSYS[i]
     RAhr=RAdeg
     RAdiff=TOTAL(ABS(RADeg-RAboundeg))/2.*3600./15.
     DEChr=DECdeg
     DECdiff=TOTAL(ABS(DECDeg-DECboundeg))/2.*3600.
     VSYSdiff=TOTAL(ABS(catVSYS[i]-ROTboun))
     totflux[0]=totflux[0]/pixperbeam
     HIMASS=2.36E5*catDistance[i]^2*totflux*ABS(channelwidth)
     convertradec,RAhr,DEChr
     getDHI,moment0map,hedmap,catPA[i],[RADeg,DECdeg,catVSYS[i]],DHI
     catmaxrotdev[i]=(VSYSdiff/2.)/SIN((newinclination[0])*!pi/180.)
     IF catmaxrot[i] LT 20*channelwidth AND catmaxrotdev[i] LT 25.*channelwidth then  catmaxrotdev[i]=25.*channelwidth
     IF catmaxrotdev[i] LT 4*channelwidth then catmaxrotdev[i]=4*channelwidth
     IF size(log,/TYPE) EQ 7 then begin
        openu,66,log,/APPEND
        printf,66,linenumber()+"We are using maxrot of "+string(catmaxrot[i])+" +/- "+string(catmaxrotdev[i])   
        close,66
     ENDIF 
     basicinfofile=maindir+'/'+catdirname[i]+'/'+basicinfo+'_'+catcubename[i]+'.txt'
     openw,1,basicinfofile
     printf,1,"#This file contains the basic parameters of the Galaxy"
     printf,1,"#The total flux is determined by the source finder in the initial estimates and the total emission in the masked cube of the final model in the fits" 
     printf,1,"#D_HI is determined as the diameter where the major axis with given PA cuts the 10^20 column of the moment0 map" 
     printf,1,format='(A25,A25,A20,A20,A20,A20,A20,A20,A20,A20,A20,A20)', 'RA','DEC','VSYS','PA','Inclination','Max VRot','V_mask', 'Tot FLux','D_HI', 'Distance','HI Mass','D_HI'
     printf,1,format='(A25,A25,A20,A20,A20,A20,A20,A20,A20,A20,A20,A20)', 'hh:mm:ss','dd:mm:ss','km/s','deg','deg','km/s','km/s', 'Jy','arcsec', 'Mpc','M_solar','kpc'
     printf,1,"#The initial input"
     printf,1,format='(A25,A25,A20,A20,A20,A20,A20,A20,A20,A20,A20,A20)',string(RAhr+'+/-'+strtrim(strcompress(string(RAdiff,format='(F6.1)')),2)),$
            string(DEChr+'+/-'+strtrim(strcompress(string(DECdiff,format='(F6.1)')),2)),$
            string(strtrim(strcompress(string(catVSYS[i])),2)+'+/-'+strtrim(strcompress(string(VSYSdiff,format='(F6.1)')),2)),$
            string(strtrim(strcompress(string(catPA[i])),2)+'+/-'+strtrim(strcompress(string(catPAdev[i],format='(F6.1)')),2)),$
            string(strtrim(strcompress(string(catinc[i])),2)+'+/-'+strtrim(strcompress(string(catincdev[i],format='(F6.1)')),2)),$
            string(strtrim(strcompress(string(catmaxrot[i])),2)+'+/-'+strtrim(strcompress(string(catmaxrotdev[i],format='(F6.1)')),2)),$
            string(strtrim(strcompress(string(VSYSdiff)),2)+'+/-'+strtrim(strcompress(string(channelwidth*(1.+vresolution),format='(F6.1)')),2)),$
            string(strtrim(strcompress(string(Totflux[0],format='(E10.3)')),2)),$
            string(strtrim(strcompress(string(DHI,format='(F8.1)')),2)),$
            string(strtrim(strcompress(string(catDistance[i])),2)),$
            string(strtrim(strcompress(string(HIMass[0],format='(E10.3)')),2)),$
            string(strtrim(strcompress(string(convertskyanglefunction(DHI,catDistance[i]),format='(F8.1)')),2))
     close,1
     
                                ;and we make a pv-diagram based on these parameters
     IF optimized then begin
        noptname=str_sep(strtrim(strcompress(currentfitcube),2),'_opt')
        nooptcube=readfits(maindir+'/'+catdirname[i]+'/'+noptname[0]+'.fits',nopthed,/NOSCALE)
        extract_pv,nooptcube,nopthed,catpa[i],xv,center=[RAdeg,DECdeg],xvheader=new_header
        writefits,maindir+'/'+catdirname[i]+'/'+noptname[0]+'_0_xv.fits',xv,new_header
     ENDIF ELSE BEGIN
        extract_pv,dummy,hed,catpa[i],xv,center=[RAdeg,DECdeg],xvheader=new_header
        writefits,maindir+'/'+catdirname[i]+'/'+catcubename[i]+'_0_xv.fits',xv,new_header
     ENDELSE
                                ;break if we do not want to use tirific
     IF finishafter EQ 0. then begin
        if setfinishafter NE 1 then begin 
           openu,1,outputcatalogue,/APPEND
           printf,1,format='(A60,A12,A80)', catDirname[i],'You have chosen to skip the fitting process after all preparations for the fit'
           close,1
        ENDIF
        bookkeeping=5
        goto,finishthisgalaxy
     ENDIF
                                ;build up moment 0 axis
     buildaxii,hedmap,xaxmom0,yaxmom0
                                ;some counters for keeping track
     prevmodification=0.
     overwrite=0.
     lastsbr=100
     counter=0.
     countsbr=0.
     ringmodifier=0
     plus2beamshift=0
     shiftcentercounter=0.
                                ;CD to the correct directory
     CD, new_dir, CURRENT=old_dir
     newrings=fix(norings[0])
                                ;Calculate the beam in pixels  and get an estimate for the peak sbr
     beaminpixels=fix(catmajbeam[i]/((ABS(pixelsizeRA)+ABS(pixelsizeDEC))/2.*3600.))
     centralarea=moment0map[fix(RApix-beaminpixels):fix(RApix+beaminpixels),fix(DECpix-beaminpixels):fix(DECpix+beaminpixels)]
     cenav=TOTAL(centralarea[WHERE(FINITE(centralarea))])/n_elements(centralarea[WHERE(FINITE(centralarea))])
     peaksbr=cenav*channelwidth/(catmajbeam[i]*catminbeam[i])     
                                ; at high inclinations rings overlap
                                ; so the peaksbr needs to be reduced
                                ; according to the inclination and the
                                ; physical size of the galaxy
     IF catinc[i] GT 70. then peaksbr=peaksbr/(convertskyanglefunction(DHI,catDistance[i])*SIN(catinc[i]*!DtoR)*0.5)   
     IF size(log,/TYPE) EQ 7 then begin
        openu,66,log,/APPEND
        printf,66,linenumber()+"We estimate the peaksbr to be "+string(peaksbr)   
        close,66
     ENDIF 
                                ;Calculate an initial surface
                                ;brightness profile based on the
                                ;central brightness and an exponential
     velaverage=catmaxrot[i]
     signvelaverage=catmaxrot[i]
     stringsbr='SBR= '+string(peaksbr)+' '+string(peaksbr)
     scalelength=convertskyanglefunction(7.5,catDistance[i],/PHYSICAL)
     for j=1,norings[0]-2 do begin
        exposbr=peaksbr*exp(-(rings[j]-rings[0])/(scalelength))
        stringsbr=stringsbr+string(exposbr)+' ' 
     endfor
     tmppos=where('SBR' EQ tirificfirstvars)
     tirificfirst[tmppos]=stringsbr
                                ;setting the values for the 2nd disk
     tmppos=where('SBR_2' EQ tirificfirstvars)
     tmp=str_sep(strtrim(strcompress(stringsbr),2),'=')
     tirificfirst[tmppos]='SBR_2='+tmp[1]
     tmp2=str_sep(strtrim(strcompress(tmp[1]),2),' ')
     

                                ;but that does mean we'd like
                                ;to change SBR a lot so we set quite
                                ;large starting steps
     string1='SBR '+strtrim(strcompress(string(norings[0],format='(I3)')),1)+' SBR_2 '+strtrim(strcompress(string(norings[0],format='(I3)')),1)
     string2='1'
     string3=strtrim(strcompress(string(cutoff[fix(norings[0])/2.],format='(E12.5)')),1)
     string4='7.5E-5'
     string5='2E-6'
     string6='5E-5'
     string7='1E-6'
     string8='3'
     string9='70'
                                ;The minimum is based on the cutoff
                                ;values and therefore we need to set
                                ;every ring in the tirific file
     for j=norings[0]-1,2,-1 do begin
        string1=string1+','+'SBR '+strtrim(strcompress(string(j,format='(I3)')),1)+' SBR_2 '+strtrim(strcompress(string(j,format='(I3)')),1)
        string2=string2+' 1' 
        string3=string3+' '+strtrim(strcompress(string(cutoff[fix(j-1)]/2.,format='(E12.5)')),1)
        string4=string4+' 7.5E-5'
        string5=string5+' 2E-6'
        string6=string6+' 5E-5'
        string7=string7+' 1E-6'
        string8=string8+' 3'
        string9=string9+' 70'
     endfor
     SBRinput1=[string1,string2,string3,string4,string5,string6,string7,string8,string9,string9]
     SBRinput2=[' SBR 1 SBR_2 1',$
                strtrim(strcompress(string(peaksbr,format='(E12.5)'))),'0','1E-5','1E-6','5E-5','1E-6','3','70','70']
                                ;Now setting some general values
     tmppos=where('INSET' EQ tirificfirstvars)
     tirificfirst[tmppos]='INSET='+strtrim(strcompress(string(currentfitcube+'.fits')))
     tmppos=where('BMIN' EQ tirificfirstvars)
     tirificfirst[tmppos]='BMIN='+strtrim(strcompress(string(catminbeam[i])))
     tmppos=where('BMAJ' EQ tirificfirstvars)
     tirificfirst[tmppos]='BMAJ='+strtrim(strcompress(string(catmajbeam[i])))
     tmppos=where('BPA' EQ tirificfirstvars)
     tirificfirst[tmppos]='BPA='+strtrim(strcompress(string(catbpa[i])))
     tmppos=where('SDIS' EQ tirificfirstvars)
     tirificfirst[tmppos]='SDIS=  8.'
     tmppos=where('SDIS_2' EQ tirificfirstvars)
     tirificfirst[tmppos]='SDIS_2=  8.'
     tmppos=where('CONDISP' EQ tirificfirstvars)
     IF vresolution EQ 0. then begin
        tirificfirst[tmppos]='CONDISP=  '+string(channelwidth*1.2/(2*SQRT(2*ALOG(2.))))
     ENDIF ELSE BEGIN
        tirificfirst[tmppos]='CONDISP=  '+string((channelwidth*(1.+vresolution))/(2*SQRT(2*ALOG(2.))))
     ENDELSE
     tmppos=where('RMS' EQ tirificfirstvars) 
     tirificfirst[tmppos]='RMS='+strtrim(strcompress(string(catnoise[i])))
     tmppos=where('DISTANCE' EQ tirificfirstvars)
     tirificfirst[tmppos]='DISTANCE='+strtrim(strcompress(string(catDistance[i])))
     prevrings=norings[0]
                                ; defining some flags that we will use
                                ; for checks on what is going on in
                                ; the fitting
                                ;the First estimate of Z0 is 0.2 kpc
                                ; if 0.2 kpc is less
                                ;than 1/4th of the beam we want to
                                ;make it a quarter of the beam at high inclinations
     If catdistance[i] EQ 1. then begin
        IF catinc[i] GT 80. then inpzval=MAX([(catmajbeam[i]*norings[0])/150.,catmajbeam[i]/4.]) else  inpzval=(catmajbeam[i]*norings[0])/150.
        tmppos=where('Z0' EQ tirificfirstvars)
        tirificfirst[tmppos]='Z0='+strtrim(strcompress(string((inpzval))))
        tmppos=where('Z0_2' EQ tirificfirstvars)
        tirificfirst[tmppos]='Z0_2='+strtrim(strcompress(string((inpzval))))
     ENDIF ELSE BEGIN
        IF catinc[i] GT 80. then inpzval=MAX([convertskyanglefunction(0.2,double(catDistance[i]),/PHYSICAL),catmajbeam[i]/4.]) else inpzval=convertskyanglefunction(0.2,double(catDistance[i]),/PHYSICAL)
        tmppos=where('Z0' EQ tirificfirstvars)
        tirificfirst[tmppos]='Z0='+strtrim(strcompress(string(inpzval)))
        tmppos=where('Z0_2' EQ tirificfirstvars)
        tirificfirst[tmppos]='Z0_2='+strtrim(strcompress(string(inpzval)))     
     ENDELSE
     constring=0.
     forcedring=0.
     lastcutrings=0.
     lastaddrings=0.
     secondtime=0.
     sofiarings=norings[0]
     stringvelocities='VROT= 0. '+string(ABS(double(catmaxrot[i])))
                                ;If we change the amount of rings we need to come back here
     sbrshift:
     case finishafter of
        1.1:rings=(findgen(norings[0]))*(catmajbeam[i]/2.)+(catmajbeam[i]/10.)
        2.1:begin
           tmpring=norings[0]-10.
           rings=dblarr(norings[0])
           rings[0:9]=(findgen(10))*catmajbeam[i]+catmajbeam[i]/5.
           rings[10:norings[0]-1]=(findgen(fix(tmpring)))*catmajbeam[i]*2+catmajbeam[i]/5.+11.*catmajbeam[i]
        end
       else:rings=(findgen(norings[0]))*catmajbeam[i]+catmajbeam[i]/5.   
     endcase
     
     
     IF size(log,/TYPE) EQ 7 then begin
        openu,66,log,/APPEND
        printf,66,linenumber()+"The number of rings used = "+strtrim(string(fix(norings[0])),2)
        close,66
     ENDIF ELSE begin
        print,linenumber()+"The number of rings used = "+strtrim(string(fix(norings[0])),2)
     endelse
     fluxadjust=0. 
                                ;let's write the number of rings to the tirific file
     tmppos=where('NUR' EQ tirificfirstvars)
     tirificfirst[tmppos]='NUR='+strtrim(strcompress(string(norings[0])),1) 
                                ;and cflux
     tmppos=where('CFLUX' EQ tirificfirstvars)
     IF float(string((totflux[0])/7.5E5)) NE 0. then tirificfirst[tmppos]='CFLUX= '+string((totflux[0])/7.5E5) else tirificfirst[tmppos]='CFLUX= 1e-5'
     tmppos=where('CFLUX_2' EQ tirificfirstvars)
      IF float(string((totflux[0])/7.5E5)) NE 0. then tirificfirst[tmppos]='CFLUX_2= '+string((totflux[0])/7.5E5) else tirificfirst[tmppos]='CFLUX_2= 1e-5'

                                ;Now we write the radii of our rings
     stringring='RADI=0.0 '
     for j=0,norings[0]-2 do begin   
        stringring=stringring+string(rings[j])+' ' 
     endfor
     tmppos=where('RADI' EQ tirificfirstvars)
     tirificfirst[tmppos]=stringring
                                ;Using the parameters from sofia or a
                                ;previous fit
     tmppos=where('VROT' EQ tirificfirstvars)
     tirificfirst[tmppos]=stringvelocities
     tmp=str_sep(strtrim(strcompress(stringvelocities),2),'=')
     tmppos=where('VROT_2' EQ tirificfirstvars)
     tirificfirst[tmppos]='VROT_2='+tmp[1]
     tmppos=where('INCL' EQ tirificfirstvars)
     tirificfirst[tmppos]='INCL='+strtrim(strcompress(string(catinc[i]))) 
     tmppos=where('PA' EQ tirificfirstvars)
     tirificfirst[tmppos]='PA='+strtrim(strcompress(string(catPA[i])))
     tmppos=where('XPOS' EQ tirificfirstvars)
     tirificfirst[tmppos]='XPOS='+strtrim(strcompress(string(RAdeg)))
     tmppos=where('YPOS' EQ tirificfirstvars)
     tirificfirst[tmppos]='YPOS='+strtrim(strcompress(string(DECdeg)))
     tmppos=where('VSYS' EQ tirificfirstvars)
     tirificfirst[tmppos]='VSYS='+strtrim(strcompress(string(catvsys[i])))
     tmppos=where('INCL_2' EQ tirificfirstvars)
     tirificfirst[tmppos]='INCL_2='+strtrim(strcompress(string(catinc[i]))) 
     tmppos=where('PA_2' EQ tirificfirstvars)
     tirificfirst[tmppos]='PA_2='+strtrim(strcompress(string(catPA[i]))) 
     tmppos=where('XPOS_2' EQ tirificfirstvars)
     tirificfirst[tmppos]='XPOS_2='+strtrim(strcompress(string(RAdeg)))
     tmppos=where('YPOS_2' EQ tirificfirstvars)
     tirificfirst[tmppos]='YPOS_2='+strtrim(strcompress(string(DECdeg)))
     tmppos=where('VSYS_2' EQ tirificfirstvars)
     tirificfirst[tmppos]='VSYS_2='+strtrim(strcompress(string(catvsys[i])))
                             
                                ; Set the tirific fitting parameters for the inclination
                                ;If the inclination is above 75 the
                                ;initial estimates must be good so we
                                ;want less change
     case 1 of      
        catinc[i] LT 40: begin
           INCLinput1=['INCL 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                       ' INCL_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1),$
                       string(catinc[i]+catincdev[i]+10),string(catinc[i]-catincdev[i]-10),string(0.2),string(0.1),string(5.0),string(0.1),'5','70','70']
           ;We put this inside might not work then take it out again
           IF INCLinput1[2] LT 0.1 then INCLinput1[2]=0.1
           IF INCLinput1[2] GT 60 then INCLinput1[2]=60
        end
        catinc[i] LT 75: begin
           INCLinput1=['INCL 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                       ' INCL_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1),$
                       string(catinc[i]+catincdev[i]+10),string(catinc[i]-catincdev[i]-10),string(1.),string(0.1),string(0.5),string(0.1),'3','70','70']
          
           IF INCLinput1[2] LT 0.1 then INCLinput1[2]=0.1
           IF INCLinput1[2] GT 60 then INCLinput1[2]=60
        end
        else:INCLinput1=['INCL 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                         ' INCL_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1),$
                         string(catinc[i]+catincdev[i]),string(catinc[i]-catincdev[i]),string(0.2),string(0.01),string(0.5),string(0.001),'3','70','70']
     endcase
                                ; have to ensure that the parameters are within the limits
     IF INCLinput1[1] GT 90. then INCLinput1[1]=90.
     IF INCLinput1[1] LT 30.then INCLinput1[1]=30.
   
                                ;Set the input for the PA
     
     PAinput1=['PA 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
               ' PA_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1),$
               string((catPA[i])+catPAdev[i]+10.),string((catPA[i])-catPAdev[i]-10),string(1.),string(0.1),string(1.),string(0.1),'3','70','70']
     
                                ;If the estimate PA isn't close
                                ;to 180 or 360 we want to block
                                ;changes that can flip the rotation curve 
     IF PAfixboun NE 'n' then begin
        IF PAfixboun EQ 'w' then begin
           IF catPA[i]+catPAdev[i]+10. GT 360 then PAinput1[1]='360'
           IF catPA[i]-catPAdev[i]-10. LT 180 then PAinput1[2]='180'
        ENDIF
        IF PAfixboun EQ 'e' then begin
           IF catPA[i]+catPAdev[i]+10. GT 180. then PAinput1[1]='180'
           IF catPA[i]-catPAdev[i]-10. LT 0. then PAinput1[2]='0'
        ENDIF
     ENDIF
                                ;Also the minimum and maximum of Z0 which have to be based on physical
                                ;values 
     Z0input1=['Z0 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
               ' Z0_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1),$
               string(convertskyanglefunction(0.5,double(catDistance[i]),/PHYSICAL)),$
               string(convertskyanglefunction(0.075,double(catDistance[i]),/PHYSICAL)),$
               string(convertskyanglefunction(0.1,double(catDistance[i]),/PHYSICAL)),$
               string(convertskyanglefunction(0.01,double(catDistance[i]),/PHYSICAL)),$
               string(convertskyanglefunction(1.0,double(catDistance[i]),/PHYSICAL)),$
               string(convertskyanglefunction(0.01,double(catDistance[i]),/PHYSICAL)),'3','70','70']
                                ;if the inclination is high we want to
                                ;make sure that the scale height can
                                ;be at least half a beam or 1 kpc
                                ;whichever is larger because otherwise
                                ;it might be pushing the inclination
                                ;down. 

     IF catinc[i] GT 80 then begin
        Z0input1[1]=string(MAX([convertskyanglefunction(1.0,double(catDistance[i])),catmajbeam[i]/2.]))
     ENDIF
     IF catDistance[i] EQ 1. then begin
                                ;ok if we do not know the distance
                                ;then let us assume each disk is about
                                ;30 kpc which means that
                                ;0.5=norings[0]/60. and so on
        IF catinc[i] GT 80 then begin
           Z0input1[1]=string(MAX([catmajbeam[i]*norings[0]/66.,catmajbeam[i]/2.]))
        ENDIF ELSE  Z0input1[1]=string(catmajbeam[i]*norings[0]/66.)
        Z0input1[2]='0.'
        Z0input1[3]=string(-1*catmajbeam[i]*norings[0]/1.5E5)
        Z0input1[4]=string(catmajbeam[i]*norings[0]/1.5E6)
        Z0input1[5]=string(catmajbeam[i]*norings[0]/1500)
        Z0input1[6]=string(catmajbeam[i]*norings[0]/1.5E6)
     ENDIF
                             
  



                                ;Now The rotation
    
     VROTmax=catmaxrotdev[i]
     VROTmin=channelwidth
                                ;ensure reasonable vrotmax dependent
                                ;on inclination as more unsure at low inclination
     IF VROTmax LT 80. then VROTmax=80.                           
     IF VROTmax GT 600. then VROTmax=600.
    
     string1='!VROT '+strtrim(strcompress(string(norings[0],format='(I3)')),1)+':2 VROT_2 '+strtrim(strcompress(string(norings[0],format='(I3)')),1)+':2'
     string2=string(VROTmax)
     string3=string(VROTmin)
     string4=string(channelwidth)
                                ;At low inclinations we want to be
                                ;very accepting of the rotational
                                ;values as they are very unsure.
     case 1 of
        catinc[i] LT 30:begin
           string4=string(2*channelwidth)
           string5=string(0.02*channelwidth)
           string6=string(0.2*channelwidth)
           string7=string(0.02*channelwidth)
        end
        else:begin
           string4=string(channelwidth)
           string5=string(0.01*channelwidth)
           string6=string(0.1*channelwidth)
           string7=string(0.01*channelwidth)
        end
     ENDCASE
     string8='3'
     string9='70'
                                ;If the model is large enough we fit
                                ;only a slope to the outer quarter of
                                ;the model at low inclinations we want
                                ;this slope to be at least half of the
                                ;rings and be flat
     IF centralexclude then begin
        IF catinc[i] LT 40 then begin
           string10='VROT 2 '+strtrim(strcompress(string(norings[0]-1,format='(I3)')),1)+':'+strtrim(strcompress(string(fix(norings[0]*0.5+2),format='(I3)')),1)+' VROT_2 2 '+strtrim(strcompress(string(norings[0]-1,format='(I3)')),1)+':'+strtrim(strcompress(string(fix(norings[0]*0.5+2),format='(I3)')),1)
        ENDIF else begin
           string10='VROT 2 '+strtrim(strcompress(string(norings[0]-1,format='(I3)')),1)+':'+strtrim(strcompress(string(fix(norings[0]*3./4.),format='(I3)')),1)+' VROT_2 2 '+strtrim(strcompress(string(norings[0]-1,format='(I3)')),1)+':'+strtrim(strcompress(string(fix(norings[0]*3./4.),format='(I3)')),1)
        ENDELSE
     ENDIF ELSE BEGIN
        IF catinc[i] LT 40 then begin
           string10='VROT '+strtrim(strcompress(string(norings[0]-1,format='(I3)')),1)+':'+strtrim(strcompress(string(fix(norings[0]*0.5+2),format='(I3)')),1)+' VROT_2 '+strtrim(strcompress(string(norings[0]-1,format='(I3)')),1)+':'+strtrim(strcompress(string(fix(norings[0]*0.5+2),format='(I3)')),1)
        ENDIF else begin
           string10='VROT '+strtrim(strcompress(string(norings[0]-1,format='(I3)')),1)+':'+strtrim(strcompress(string(fix(norings[0]*3./4.),format='(I3)')),1)+' VROT_2 '+strtrim(strcompress(string(norings[0]-1,format='(I3)')),1)+':'+strtrim(strcompress(string(fix(norings[0]*3./4.),format='(I3)')),1)
        ENDELSE
     ENDELSE
     IF norings[0] GT 4 then begin
        VROTinput1=[string1,string2,string3,string4,string5,string6,string7,string8,string9,string9,string10]
     ENDIF ELSE BEGIN
        VROTinput1=[string1,string2,string3,string4,string5,string6,string7,string8,string9,string9,' ']
     ENDELSE  
     vsysinput1=[ ' VSYS 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                  ' VSYS_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1),+$
                  strtrim(strcompress(string(catvsys[i]+100.)),1),strtrim(strcompress(string(catvsys[i]-100.)),1),'10','0.01','0.5','0.01','3','70','70']
     
                                ;IF we have a small cube we will only accept very small changes to the
                                ;central position                              
     IF norings[0] LT 4 then begin
        maxxpos=strcompress(string(RAdeg+catmajbeam[i]/3600.))
        minxpos=strcompress(string(RAdeg-catmajbeam[i]/3600.))
        maxypos=strcompress(string(DECdeg+catmajbeam[i]/3600.))
        minypos=strcompress(string(DECdeg-catmajbeam[i]/3600.))
        xposinput1=[ ' XPOS 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                     ' XPOS_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1), maxxpos,minxpos,'5E-6','1E-6','1E-4','1E-6','3','70','70']
        yposinput1=[ ' YPOS 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                     ' YPOS_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1), maxypos,minypos,'5E-6','1E-6','1E-4','1E-6','3','70','70']
     ENDIF ELSE BEGIN
        maxxpos=strcompress(string(RAdeg+3.*catmajbeam[i]/3600.))
        minxpos=strcompress(string(RAdeg-3.*catmajbeam[i]/3600.))
        maxypos=strcompress(string(DECdeg+3.*catmajbeam[i]/3600.))
        minypos=strcompress(string(DECdeg-3.*catmajbeam[i]/3600.))
        xposinput1=[ ' XPOS 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                     ' XPOS_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1), maxxpos,minxpos,'1E-5','1E-6','1E-4','1E-6','3','70','70']
        yposinput1=[ ' YPOS 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                     ' YPOS_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1), maxypos,minypos,'1E-5','1E-6','1E-4','1E-6','3','70','70']
     ENDELSE
     tryone=0.
                                ;   IF we have a low inclination we first want to fit the PA by itself
     PAest=0.
     INCLest=0.
     IF catinc[i]  LT 50 AND counter EQ 0. then begin
                                ;If we have a small number of beams
                                ;acros the minor axis the inclination
                                ;is very unsure and we first want to
                                ;fit the inclination
        IF ceil(norings[0]*COS(catinc[i]*!DtoR)) LE 5 then begin
           goto,notfornow
           INCLest=1
           INCLinput1=['INCL 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                       ' INCL_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1),$
                       '90.',string(catinc[i]-catincdev[i]),string(0.5),string(0.1),string(1.0),string(0.1),'3','70','70']  
           PAinput1=['PA 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                     ' PA_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1),$
                     string((catPA[i])+catPAdev[i]+40.),string((catPA[i])-catPAdev[i]-40),string(6),string(0.1),string(0.05),string(0.01),'3','70','70']

           IF PAfixboun NE 'n' then begin
              IF PAfixboun EQ 'w' then begin
                 IF catPA[i]+catPAdev[i]+40. GT 360 then PAinput1[1]='360'
                 IF catPA[i]-catPAdev[i]-40. LT 180 then PAinput1[2]='180'
              ENDIF
              IF PAfixboun EQ 'e' then begin
                 IF catPA[i]+catPAdev[i]+40. GT 180. then PAinput1[1]='180'
                 IF catPA[i]-catPAdev[i]-40. LT 0. then PAinput1[2]='0'
              ENDIF
           ENDIF
           VROTinputINCL=['VROT 2:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                          ' VROT_2 2:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)$
                          ,string(VROTmax),string(0.),$
                          string(channelwidth),string(0.1*channelwidth),string(channelwidth),string(0.01*channelwidth),'3','70','70']
        
           Writefittingvariables,tirificfirst,inclinput1,VROTinputINCL,painput1,sbrinput1,sbrinput2
           againINCLestimate:
           IF size(log,/TYPE) EQ 7 then begin
              openu,66,log,/APPEND
              printf,66,linenumber()+"Because there are only a few beams across the minor axis we first adjust the inclination "+catDirname[i]+" which is galaxy # "+strtrim(string(fix(i)),2)+" at "+systime()
              close,66
           ENDIF 
           openw,1,maindir+'/'+catdirname[i]+'/tirific.def'
           for index=0,n_elements(tirificfirst)-1 do begin
              printf,1,tirificfirst[index]
           endfor
           close,1
           gipsyfirst=strarr(1)
           IF size(log,/TYPE) EQ 7 then begin
              openu,66,log,/APPEND
              printf,66,linenumber()+"Starting tirific the INCL estimate in  "+catDirname[i]+" which is galaxy # "+strtrim(string(fix(i)),2)+" at "+systime()
              close,66
           ENDIF
           
           IF testing GE 1 then goto,testing1INCL
           print,linenumber()+"Starting tirific the INCL estimate in  "+catDirname[i]+" which is galaxy # "+strtrim(string(fix(i)),2)+" at "+systime()
           gipsyfirst='tirific DEFFILE=tirific.def ACTION=1'
           spawn,gipsyfirst,isthere2
           ; Let's check how we did in the fit
           get_progress,maindir+'/'+catdirname[i]+'/progress1.txt',AC1,nopoints,loops,toymodels
           ;If not accepted we try again
           IF AC1 EQ 0. and INCLest LT 5 then begin
              IF size(log,/TYPE) EQ 7 then begin
                 openu,66,log,/APPEND
                 printf,66,linenumber()+"The INCL estimate was not accepted try again. Tries: "+string(INCLest) 
                 close,66
              ENDIF 
              VariablesWanted=['PA','PA_2','INCL','INCL_2','VROT','VROT_2']
              firstfitvalues=0.
              writenewtotemplate,tirificfirst,maindir+'/'+catdirname[i]+'/1stfit.def',Arrays=firstfitvalues,VariableChange=VariablesWanted
              catPA[i]=firstfitvalues[0,0]
              catinc[i]=firstfitvalues[0,2]
              catmaxrot[i]=firstfitvalues[1,4]
              PAinput1=['PA 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                        ' PA_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1),$
                        string((catPA[i])+catPAdev[i]+30.),string((catPA[i])-catPAdev[i]-30),string(5),string(0.1),string(0.05),string(0.01),'3','70','70']
              IF PAfixboun NE 'n' then begin
                 IF PAfixboun EQ 'w' then begin
                    IF catPA[i]+catPAdev[i]+40. GT 360 then PAinput1[1]='360'
                    IF catPA[i]-catPAdev[i]-40. LT 180 then PAinput1[2]='180'
                 ENDIF
                 IF PAfixboun EQ 'e' then begin
                    IF catPA[i]+catPAdev[i]+40. GT 180. then PAinput1[1]='180'
                    IF catPA[i]-catPAdev[i]-40. LT 0. then PAinput1[2]='0'
                 ENDIF
              ENDIF
              INCLinput1=['INCL 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                          ' INCL_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1),$
                          '90.',string(catinc[i]-catincdev[i]-20),string(3),string(0.1),string(0.5),string(0.1),'3','70','70']  
              Writefittingvariables,tirificfirst,painput1,inclinput1,vrotinputINCL
              INCLest=INCLest+1.
              goto,againINCLestimate
           ENDIF ELSE BEGIN
              
              ;When accepted update our values
              VariablesWanted=['PA','PA_2','INCL','INCL_2','VROT','VROT_2']
              firstfitvalues=0.
              writenewtotemplate,tirificfirst,maindir+'/'+catdirname[i]+'/1stfit.def',Arrays=firstfitvalues,VariableChange=VariablesWanted
              catPA[i]=firstfitvalues[0,0]
              catinc[i]=firstfitvalues[0,2]
              catmaxrot[i]=firstfitvalues[1,4]
           ENDELSE
           notfornow:
        ENDIF
        testing1INCL:
        INCLinput1=['INCL 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                    ' INCL_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1),$
                    string(catinc[i]+catincdev[i]+20),string(catinc[i]-catincdev[i]-20),string(0.5),string(0.1),string(5.),string(0.1),'3','70','70']
        IF INCLinput1[1] GT 90. then INCLinput1[1]=90
        IF INCLinput1[1] LT 30 then INCLinput1[1]=30
        IF INCLinput1[2] LT 0 then INCLinput1[2]=0
        IF INCLinput1[2] GT 60 then INCLinput1[2]=60
        PAinput1=['PA 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                  ' PA_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1),$
                  string((catPA[i])+catPAdev[i]+30.),string((catPA[i])-catPAdev[i]-30),string(2.5),string(0.1),string(0.05),string(0.01),'3','70','70']
        IF PAfixboun NE 'n' then begin
           IF PAfixboun EQ 'w' then begin
              IF catPA[i]+catPAdev[i]+30. GT 360 then PAinput1[1]='360'
              IF catPA[i]-catPAdev[i]-30. LT 180 then PAinput1[2]='180'
           ENDIF
           IF PAfixboun EQ 'e' then begin
              IF catPA[i]+catPAdev[i]+30. GT 180. then PAinput1[1]='180'
              IF catPA[i]-catPAdev[i]-30. LT 0. then PAinput1[2]='0'
           ENDIF
        ENDIF
        VROTinputPA=['VROT 2:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                     ' VROT_2 2:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)$
                     ,string(VROTmax),string(0.),$
                     string(channelwidth),string(0.01*channelwidth),string(channelwidth),string(0.01*channelwidth),'3','70','70']
        Writefittingvariables,tirificfirst,painput1,sbrinput1,sbrinput2,VROTinputPA
        
        againPAestimate:

        openw,1,maindir+'/'+catdirname[i]+'/tirific.def'
        for index=0,n_elements(tirificfirst)-1 do begin
           printf,1,tirificfirst[index]
        endfor
        close,1
        gipsyfirst=strarr(1)
        IF size(log,/TYPE) EQ 7 then begin
           openu,66,log,/APPEND
           printf,66,linenumber()+"Starting tirific the PA estimate in  "+catDirname[i]+" which is galaxy # "+strtrim(string(fix(i)),2)+" at "+systime()
           close,66
        ENDIF        
        IF testing GE 1 then goto,testing1PA
        print,linenumber()+"Starting tirific the PA estimate in "+catDirname[i]+" which is galaxy # "+strtrim(string(fix(i)),2)+" at "+systime()
        gipsyfirst='tirific DEFFILE=tirific.def ACTION=1'
        spawn,gipsyfirst,isthere2
                                ;get the results
        get_progress,maindir+'/'+catdirname[i]+'/progress1.txt',AC1,nopoint,loops,toymodels     
                                ;If failed try again
        IF AC1 EQ 0. and PAest LT 5 then begin
           IF size(log,/TYPE) EQ 7 then begin
              openu,66,log,/APPEND
              printf,66,linenumber()+"The Pa estimate was not accepted try again. Tries: "+string(PAest) 
              close,66
           ENDIF 
           VariablesWanted=['PA','PA_2']
           firstfitvalues=0.
           writenewtotemplate,tirificfirst,maindir+'/'+catdirname[i]+'/1stfit.def',Arrays=firstfitvalues,VariableChange=VariablesWanted
           catPA[i]=firstfitvalues[0,0]
           PAinput1=['PA 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                     ' PA_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1),$
                     string((catPA[i])+catPAdev[i]+30.),string((catPA[i])-catPAdev[i]-30),string(5),string(0.1),string(0.05),string(0.01),'3','70','70']
           IF PAfixboun NE 'n' then begin
              IF PAfixboun EQ 'w' then begin
                 IF catPA[i]+catPAdev[i]+30. GT 360 then PAinput1[1]='360'
                 IF catPA[i]-catPAdev[i]-30. LT 180 then PAinput1[2]='180'
              ENDIF
              IF PAfixboun EQ 'e' then begin
                 IF catPA[i]+catPAdev[i]+30. GT 180. then PAinput1[1]='180'
                 IF catPA[i]-catPAdev[i]-30. LT 0. then PAinput1[2]='0'
              ENDIF
           ENDIF
           Writefittingvariables,tirificfirst,painput1,vrotinputPA
           PAest=PAest+1.
           goto,againPAestimate
        ENDIF
        testing1PA:
                                ;get the values from these adapted
                                ;models and write them into the
                                ;tirific array
        VariablesWanted=['PA','PA_2','INCL','INCL_2','VROT','VROT_2']
        firstfitvalues=0.
        writenewtotemplate,tirificfirst,maindir+'/'+catdirname[i]+'/1stfit.def',Arrays=firstfitvalues,VariableChange=VariablesWanted
                                ;if the change is large we want to
                                ;adjust our inc and vrot estimates as
                                ;well
        newPA=[firstfitvalues[0,0],ABS(firstfitvalues[0,0]-newPA[0])+newPA[1]]
        IF INCLest GT 0. then begin
           newinclination=[firstfitvalues[0,2],20.]
           newrot=firstfitvalues[1,4]
        ENDIF else begin
           obtain_inclinationv8,moment0map,newPA,newinclination,[RApix,DECpix],extend=noringspix,noise=momnoise,beam=catmajbeam[i]/(pixelsizeRA*3600.) 
           
           IF newinclination[0] LT 40 then newrot=W50/2./SIN(ABS(newinclination[0]+5.0)*!pi/180.) else newrot=W50/2./SIN(newinclination[0]*!pi/180.) 
        ENDELSE
                                ;when we get a new inclination we need
                                ;to reset our cutoff values.
        cutoffcorrection=SIN(75.*!DtoR)/SIN(newinclination[0]*!DtoR)
        IF newinclination[0] GT 50 then cutoffcorrection=1.
        IF newinclination[0] LT 50 AND newinclination[0] GT 40  then cutoffcorrection=1.+(50-newinclination[0])*0.05
        IF cutoffcorrection GT 2.5 then cutoffcorrection=2.5
       
        IF size(log,/TYPE) EQ 7 then begin
           openu,66,log,/APPEND
           printf,66,linenumber()+"Cutoff values will be adjusted for inclination by multiplying them with"+string(cutoffcorrection)
           close,66
        ENDIF
        cutoff=cutoffor*cutoffcorrection   
    
        IF INCLest GT 0. then begin
           IF size(log,/TYPE) EQ 7 then begin
              openu,66,log,/APPEND
              printf,66,linenumber()+"Because of a low number of beams across the minor axis,"
              printf,66,linenumber()+"We have adjusted the PA from to "+string(catPA[i])+'+/-'+string(newPA[1])
              printf,66,linenumber()+"The inclination to "+string(catInc[i])+'+/-'+string(newinclination[1])
              printf,66,linenumber()+"Maxrot to "+string(newrot[0])+'+/-'+string(catmaxrotdev[i])
              printf,66,linenumber()+"We have adjusted the fit setting of the inclination and VROT." 
              close,66
           ENDIF
        ENDIF ELSE BEGIN
           IF size(log,/TYPE) EQ 7 then begin
              openu,66,log,/APPEND
              printf,66,linenumber()+"Because of low inclination "+string(catinc[i])
              printf,66,linenumber()+"We have adjusted the PA from "+string(catPA[i])+" to "+string(firstfitvalues[0,0])+'+/-'+string(newPA[1])
              printf,66,linenumber()+"The inclination from "+string(catInc[i])+" to "+string(newinclination[0])+'+/-'+string(newinclination[1])
              printf,66,linenumber()+"We have adjusted the fit setting of the inclination and VROT." 
              close,66
           ENDIF
        ENDELSE
        catPAdev[i]=15.
        catinc[i]=newinclination[0]
        catincdev[i]=newinclination[1]
        catPA[i]=firstfitvalues[0,0]
                                ;adapt the inc values
        tmppos=where('INCL' EQ tirificfirstvars)
        tirificfirst[tmppos]='INCL='+strtrim(strcompress(string(catinc[i]))) 
        tmppos=where('INCL_2' EQ tirificfirstvars)
        tirificfirst[tmppos]='INCL_2='+strtrim(strcompress(string(catinc[i]))) 
        INCLinput1=['INCL 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                    ' INCL_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1),$
                    string(catinc[i]+catincdev[i]+10),string(catinc[i]-catincdev[i]-10),string(0.2),string(0.1),string(1.),string(0.1),'3','70','70']

        IF INCLinput1[1] GT 90. then INCLinput1[1]=90
        IF INCLinput1[1] LT 30 then INCLinput1[1]=30
        
        IF INCLinput1[2] LT 0.1 then INCLinput1[2]=0.1
        IF INCLinput1[2] GT 60 then INCLinput1[2]=60
                                ;and the vrot values
        
        VROTmax=catmaxrotdev[i]
        VROTmin=channelwidth
                                ;WE want to ensure a reasonable vrot max
        IF VROTmax LT 80 then VROTmax=80.
        IF VROTmax GT 600. then VROTmax=600.       
        string1='!VROT '+strtrim(strcompress(string(norings[0],format='(I3)')),1)+':2 VROT_2 '+strtrim(strcompress(string(norings[0],format='(I3)')),1)+':2'
        string2=string(VROTmax)
        string3=string(VROTmin)
        string4=string(2*channelwidth)
        string5=string(0.01*channelwidth)
        string6=string(0.2*channelwidth)
        string7=string(0.01*channelwidth)
        string8='3'
        string9='70'
        IF centralexclude then begin
           IF catinc[i] LT 40 then begin
              string10='VROT 2 '+strtrim(strcompress(string(norings[0]-1,format='(I3)')),1)+':'+strtrim(strcompress(string(fix(norings[0]*0.5+2),format='(I3)')),1)+' VROT_2 2 '+strtrim(strcompress(string(norings[0]-1,format='(I3)')),1)+':'+strtrim(strcompress(string(fix(norings[0]*0.5+2),format='(I3)')),1)
           ENDIF else begin
              string10='VROT 2 '+strtrim(strcompress(string(norings[0]-1,format='(I3)')),1)+':'+strtrim(strcompress(string(fix(norings[0]*3./4.),format='(I3)')),1)+' VROT_2 2 '+strtrim(strcompress(string(norings[0]-1,format='(I3)')),1)+':'+strtrim(strcompress(string(fix(norings[0]*3./4.),format='(I3)')),1)
           ENDELSE
        ENDIF ELSE BEGIN
           IF catinc[i] LT 40 then begin
              string10='VROT '+strtrim(strcompress(string(norings[0]-1,format='(I3)')),1)+':'+strtrim(strcompress(string(fix(norings[0]*0.5+2),format='(I3)')),1)+' VROT_2 '+strtrim(strcompress(string(norings[0]-1,format='(I3)')),1)+':'+strtrim(strcompress(string(fix(norings[0]*0.5+2),format='(I3)')),1)
           ENDIF else begin
              string10='VROT '+strtrim(strcompress(string(norings[0]-1,format='(I3)')),1)+':'+strtrim(strcompress(string(fix(norings[0]*3./4.),format='(I3)')),1)+' VROT_2 '+strtrim(strcompress(string(norings[0]-1,format='(I3)')),1)+':'+strtrim(strcompress(string(fix(norings[0]*3./4.),format='(I3)')),1)
           ENDELSE
        ENDELSE
   
        IF norings[0] GT 4 then begin
           VROTinput1=[string1,string2,string3,string4,string5,string6,string7,string8,string9,string9,string10]
        ENDIF ELSE BEGIN
           VROTinput1=[string1,string2,string3,string4,string5,string6,string7,string8,string9,string9,' ']
        ENDELSE
        tmppos=where('VROT' EQ tirificfirstvars)
        tirificfirst[tmppos]='VROT= 0. '+string(catmaxrot[i])
        tmppos=where('VROT_2' EQ tirificfirstvars)
        tirificfirst[tmppos]='VROT_2= 0. '+string(catmaxrot[i])
     ENDIF
                                ;Since our first SBR estimate is very
                                ;crude we want to adapt that by itself first
     IF counter EQ 0 then begin  
        Writefittingvariables,tirificfirst,sbrinput1,sbrinput2,painput1,vrotinput1
        openw,1,maindir+'/'+catdirname[i]+'/tirific.def'
        for index=0,n_elements(tirificfirst)-1 do begin
           printf,1,tirificfirst[index]
        endfor
        close,1
        IF size(log,/TYPE) EQ 7 then begin
           openu,66,log,/APPEND
           printf,66,linenumber()+"Starting tirific SBR estimate in  "+catDirname[i]+" which is galaxy # "+strtrim(string(fix(i)),2)+" at "+systime()
           close,66
        ENDIF 
        print,linenumber()+"Starting tirific SBR estimate in  "+catDirname[i]+" which is galaxy # "+strtrim(string(fix(i)),2)+" at "+systime()
        gipsyfirst=strarr(1)
        IF testing GE 1 then goto,testing1SBR
        gipsyfirst='tirific DEFFILE=tirific.def ACTION=1'
        testing1SBR:
        spawn,gipsyfirst,isthere2
        VariablesWanted=['RADI','SBR','SBR_2','PA','PA_2','VROT','VROT_2']
        firstfitvalues=0.
        get_progress,maindir+'/'+catdirname[i]+'/progress1.txt',AC1,nopoints,loops,toymodels
        check_cflux,nopoints,tirificfirst,tirificfirstvars,cfluxadjusted,log=log
        IF cfluxadjusted then  VariablesWanted=['RADI','SBR','SBR_2','VROT','VROT_2']
        writenewtotemplate,tirificfirst,maindir+'/'+catdirname[i]+'/1stfit.def',Arrays=firstfitvalues,VariableChange=VariablesWanted
        IF not cfluxadjusted then begin
           newPA=[firstfitvalues[0,4],ABS(firstfitvalues[0,4]-newPA[0])+newPA[1]]
        ENDIF ELSE BEGIN
           tmp=firstfitvalues
           firstfitvalues=dblarr(n_elements(tmp[*,0]),n_elements(tmp[0,*])+2)
           firstfitvalues[*,0:2]=tmp[*,0:2]
           firstfitvalues[*,3:4]=newPA[0]
           firstfitvalues[*,5:6]=tmp[*,3:4]
           VariablesWanted=['RADI','SBR','SBR_2','PA','PA_2','VROT','VROT_2']
        ENDELSE

        
        IF INCLest LT 1 then begin
           obtain_inclinationv8,moment0map,newPA,newinclination,[RApix,DECpix],extend=noringspix,noise=momnoise,beam=catmajbeam[i]/(pixelsizeRA*3600.)
                                ;when we get a new inclination we need
                                ;to reset our cutoff values.
           cutoffcorrection=SIN(75.*!DtoR)/SIN(newinclination[0]*!DtoR)
           IF newinclination[0] GT 50 then cutoffcorrection=1.
           IF newinclination[0] LT 50 AND newinclination[0] GT 40  then cutoffcorrection=1.+(50-newinclination[0])*0.05
           IF cutoffcorrection GT 2.5 then cutoffcorrection=2.5
           IF size(log,/TYPE) EQ 7 then begin
              openu,66,log,/APPEND
              printf,66,linenumber()+"Cutoff values will be adjusted for inclination by multiplying them with"+string(cutoffcorrection)
              close,66
           ENDIF
           cutoff=cutoffor*cutoffcorrection
        ENDIF
        tmppos=where('SBR' EQ VariablesWanted)
        SBRarr=firstfitvalues[*,tmppos]
        tmppos=where('SBR_2' EQ VariablesWanted)
        SBRarr2=firstfitvalues[*,tmppos]
                                ;checking the surface brightness
        sbr_check,tirificfirst, tirificfirstvars,sbrarr,sbrarr2,cutoff
        IF newinclination[0] LT 40 then newrot=W50/2./SIN(ABS(newinclination[0]+5)*!pi/180.) else newrot=W50/2./SIN(newinclination[0]*!pi/180.)
        catmaxrotdev[i]=(VSYSdiff/2.)/SIN(newinclination[0]*!pi/180.)
        IF catmaxrotdev[i] LT 3*channelwidth then catmaxrotdev[i]=3*channelwidth
        IF size(log,/TYPE) EQ 7 then begin
           openu,66,log,/APPEND
           printf,66,linenumber()+"Because of changing the SBR,"
           printf,66,linenumber()+"We have adjusted the PA from "+string(catPA[i])+" to "+string(firstfitvalues[0,4])
           printf,66,linenumber()+"The inclination from "+string(catInc[i])+" to "+string(newinclination[0])+'+/-'+string(newinclination[1])
           printf,66,linenumber()+"Maxrot from "+string(catmaxrot[i])+" to "+string(newrot[0])+'+/-'+string(catmaxrotdev[i])
           close,66
        ENDIF
                                ;write new values to overall arrays
        catPAdev[i]=15.
        catinc[i]=newinclination[0]
        catincdev[i]=newinclination[1]
        catmaxrot[i]=newrot         
        catPA[i]=firstfitvalues[0,4]
                                ;adapt the inc values
        tmppos=where('INCL' EQ tirificfirstvars)
        tirificfirst[tmppos]='INCL='+strtrim(strcompress(string(catinc[i]))) 
        tmppos=where('INCL_2' EQ tirificfirstvars)
        tirificfirst[tmppos]='INCL_2='+strtrim(strcompress(string(catinc[i]))) 
                                ;Inclination is tricky and the settings
                                ;depend greatly on whether high or low inclination
        case 1 of
           catinc[i] GE 75:BEGIN
              INCLinput1=['INCL 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                          ' INCL_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1),$
                          string(catinc[i]+catincdev[i]),string(catinc[i]-catincdev[i]),string(0.1),string(0.01),string(0.5),string(0.01),'3','70','70']
           END
           catinc[i] GE 40:BEGIN
              INCLinput1=['INCL 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                          ' INCL_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1),$
                          string(catinc[i]+catincdev[i]+20),string(catinc[i]-catincdev[i]-20),string(0.2),string(0.1),string(1.),string(0.1),'3','70','70']
           END
           catinc[i] LT 30:begin
              INCLinput1=['INCL 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                          ' INCL_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1),$
                          string(catinc[i]+catincdev[i]+20),string(catinc[i]-catincdev[i]-20),string(0.2),string(0.05),string(3),string(0.05),'5','70','70']
           end
           catinc[i] LT 40:begin
              INCLinput1=['INCL 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                          ' INCL_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1),$
                          string(catinc[i]+catincdev[i]+20),string(catinc[i]-catincdev[i]-20),string(0.2),string(0.05),string(3),string(0.05),'5','70','70']
           end
           ELSE:BEGIN
              INCLinput1=['INCL 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                          ' INCL_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1),$
                          string(catinc[i]+catincdev[i]+20),string(catinc[i]-catincdev[i]-20),string(0.2),string(0.1),string(1.),string(0.1),'3','70','70']
           END
        endcase
        IF INCLinput1[1] GT 90. then INCLinput1[1]=90
        IF INCLinput1[1] LT 30 then INCLinput1[1]=30
        IF catinc[i] LT 75 then begin
           IF INCLinput1[2] LT 0.1 then INCLinput1[2]=0.1
           IF INCLinput1[2] GT 60 then INCLinput1[2]=60
        endif
                                ;and the vrot values
        VROTmax=catmaxrotdev[i]
        VROTmin=channelwidth
                                ;We want to ensure that we have a
                                ;reasonable Vrot max 
        IF VROTmax LT 80. then VROTmax=80.
        IF VROTmax GT 600. then VROTmax=600.
        string1='!VROT '+strtrim(strcompress(string(norings[0],format='(I3)')),1)+':2 VROT_2 '+strtrim(strcompress(string(norings[0],format='(I3)')),1)+':2'
        string2=string(VROTmax)
        string3=string(VROTmin)
        
        IF catinc[i] LT 40 then begin
           string4=string(channelwidth)
           string5=string(0.1*channelwidth)
           string6=string(channelwidth)
           string7=string(0.1*channelwidth)
        ENDIF ELSE BEGIN
           string4=string(2*channelwidth)
           string5=string(0.01*channelwidth)
           string6=string(0.2*channelwidth)
           string7=string(0.01*channelwidth)
        ENDELSE
        string8='3'
        string9='70'
        IF centralexclude then begin
           IF catinc[i] LT 40 then begin
              string10='VROT 2 '+strtrim(strcompress(string(norings[0]-1,format='(I3)')),1)+':'+strtrim(strcompress(string(fix(norings[0]*0.5+2),format='(I3)')),1)+' VROT_2 2 '+strtrim(strcompress(string(norings[0]-1,format='(I3)')),1)+':'+strtrim(strcompress(string(fix(norings[0]*0.5+2),format='(I3)')),1)
           ENDIF else begin
              string10='VROT 2 '+strtrim(strcompress(string(norings[0]-1,format='(I3)')),1)+':'+strtrim(strcompress(string(fix(norings[0]*3./4.),format='(I3)')),1)+' VROT_2 2 '+strtrim(strcompress(string(norings[0]-1,format='(I3)')),1)+':'+strtrim(strcompress(string(fix(norings[0]*3./4.),format='(I3)')),1)
           ENDELSE
        ENDIF ELSE BEGIN
           IF catinc[i] LT 40 then begin
              string10='VROT '+strtrim(strcompress(string(norings[0]-1,format='(I3)')),1)+':'+strtrim(strcompress(string(fix(norings[0]*0.5+2),format='(I3)')),1)+' VROT_2 '+strtrim(strcompress(string(norings[0]-1,format='(I3)')),1)+':'+strtrim(strcompress(string(fix(norings[0]*0.5+2),format='(I3)')),1)
           ENDIF else begin
              string10='VROT '+strtrim(strcompress(string(norings[0]-1,format='(I3)')),1)+':'+strtrim(strcompress(string(fix(norings[0]*3./4.),format='(I3)')),1)+' VROT_2 '+strtrim(strcompress(string(norings[0]-1,format='(I3)')),1)+':'+strtrim(strcompress(string(fix(norings[0]*3./4.),format='(I3)')),1)
           ENDELSE
        ENDELSE
      
        IF norings[0] GT 4 then begin 
           VROTinput1=[string1,string2,string3,string4,string5,string6,string7,string8,string9,string9,string10]
        ENDIF ELSE BEGIN
           VROTinput1=[string1,string2,string3,string4,string5,string6,string7,string8,string9,string9,' ']
        ENDELSE
        tmppos=where('VROT' EQ VariablesWanted)
        VROTarr=firstfitvalues[*,tmppos]
        tmppos=where('VROT' EQ tirificfirstvars)
        tirificfirst[tmppos]='VROT= 0. '+STRJOIN(VROTarr[1:n_elements(VROTarr)-1],' ')
        tmppos=where('VROT_2' EQ tirificfirstvars)
        tirificfirst[tmppos]='VROT_2= 0. '+STRJOIN(VROTarr[1:n_elements(VROTarr)-1],' ')
     ENDIF
                                ;Setting some SBR values for the fitting
     minstep=strtrim(strcompress(string(cutoff[norings[0]-1]/3.,format='(E8.1)')),1)
     IF countsbr EQ 0. then begin
        startstep='1E-5' 
        minstep='1E-6'
     ENDIF else begin
        startstep=strtrim(strcompress(string(5.*cutoff[norings[0]-1],format='(E8.1)')),1)
        minstep=strtrim(strcompress(string(cutoff[norings[0]-1]/2.,format='(E8.1)')),1)
     ENDELSE

     satlevel=strtrim(strcompress(string(5.*double(startstep),format='(E8.1)')),1)
     string1='SBR '+strtrim(strcompress(string(fix(norings[0]),format='(I3)')),1)+' SBR_2 '+strtrim(strcompress(string(fix(norings[0]),format='(I3)')),1)
     string2='1'
     string3=strtrim(strcompress(string(cutoff[fix(norings[0])]/2.,format='(E12.5)')),1)
     string4=startstep
     string5=minstep
     string6=strtrim(strcompress(string(5.*double(startstep),format='(E8.1)')),1)
     string7=minstep
     string8='3'
     string9='70'

     for j=norings[0]-1,2,-1 do begin
        string1=string1+','+'SBR '+strtrim(strcompress(string(j,format='(I3)')),1)+' SBR_2 '+strtrim(strcompress(string(j,format='(I3)')),1)
        string2=string2+' 1' 
        string3=string3+' '+strtrim(strcompress(string(cutoff[fix(j-1)]/2.,format='(E12.5)')),1)
        string4=string4+' '+startstep
        string5=string5+' '+minstep
        string6=string6+' '+strtrim(strcompress(string(5.*double(startstep),format='(E8.1)')),1)
        string7=string7+' '+minstep
        string8=string8+' 3'
        string9=string9+' 70'
     endfor
     SBRinput1=[string1,string2,string3,string4,string5,string6,string7,string8,string9,string9]
                                ;The inner rings should be fitted as one
     SBRinput2=[' SBR 1 SBR_2 1',$
                strtrim(strcompress(string(ABS(SBRarr[1]),format='(E12.5)'))) ,'0','1E-5','1E-6','5E-5','1E-6','3','70','70']
     smoothrotation=0.
     keepcenter=0.
     fixedcenter=0.
     shiftcenter:
     notacceptedone:    
     IF smoothrotation then vrotinput1=vrotinputoriginal
     smoothrotation=0.
     smoothingone:
                                ;write the paramter to the array based
                                ;on inclination
     if not keepcenter and not fixedcenter and not cfluxadjusted then begin
        case 1 of
           catinc[i] LT 30: Writefittingvariables,tirificfirst,xposinput1,yposinput1,vsysinput1,painput1,vrotinput1,sbrinput1,sbrinput2,inclinput1
           catinc[i] LT 50:Writefittingvariables,tirificfirst,xposinput1,yposinput1,vsysinput1,painput1,vrotinput1,inclinput1,sbrinput1,sbrinput2
           catinc[i] GT 75:Writefittingvariables,tirificfirst,xposinput1,yposinput1,vsysinput1,painput1,vrotinput1,sbrinput1,sbrinput2,inclinput1,z0input1
           else:Writefittingvariables,tirificfirst,xposinput1,yposinput1,vsysinput1,painput1,inclinput1,vrotinput1,sbrinput1,sbrinput2
        endcase
     ENDIF ELSE BEGIN
        case 1 of
           catinc[i] LT 30:Writefittingvariables,tirificfirst,painput1,vrotinput1,sbrinput1,sbrinput2,inclinput1
           catinc[i] LT 50:Writefittingvariables,tirificfirst,painput1,vrotinput1,inclinput1,sbrinput1,sbrinput2
           catinc[i] GT 75:Writefittingvariables,tirificfirst,painput1,vrotinput1,sbrinput1,sbrinput2,inclinput1,z0input1
           else:Writefittingvariables,tirificfirst,painput1,inclinput1,vrotinput1,sbrinput1,sbrinput2
        endcase
     ENDELSE
                                ;If we tried fifty times and have not
                                ;found a reasonable solution then we
                                ;give up
     IF counter GT 50 then begin
        IF size(log,/TYPE) EQ 7 then begin
           openu,66,log,/APPEND
           printf,66,linenumber()+"We have rerun "+string(counter)+" times and still haven't found an acceptable solution."
           printf,66,linenumber()+"We have modified the rings "+string(countsbr)+" times."
           printf,66,linenumber()+"We are aborting this galaxy."
           close,66
        ENDIF
        openu,1,outputcatalogue,/APPEND
        printf,1,format='(A60,A80)', catDirname[i],'We could not find a proper center.'
        close,1 
        bookkeeping=0
        goto,finishthisgalaxy
     ENDIF
     ;Write the current array to the def file
     openw,1,maindir+'/'+catdirname[i]+'/tirific.def'
     for index=0,n_elements(tirificfirst)-1 do begin
        printf,1,tirificfirst[index]
     endfor
     close,1
                                ;And starting the first fit while
                                ;moving any previous fits to old.
                                ;if this is the second fit we do want
                                ;to keep the first one not being old
     IF testing GE 1 then goto,testing1
     IF tryone EQ 1 then begin
        rename,'1stfit.','1stfitall.'
     endif
     counter++
     IF size(log,/TYPE) EQ 7 then begin
        openu,66,log,/APPEND
        printf,66,linenumber()+"Starting tirific first estimate in  "+catDirname[i]+"  which is galaxy #  "+strtrim(string(fix(i)),2)+" at "+systime()
        printf,66,linenumber()+"We have rerun "+string(counter)+" times."
        printf,66,linenumber()+"We have modified the rings "+string(countsbr)+" times."
        close,66
     ENDIF 
     print,linenumber()+"Starting tirific first estimate in  "+catDirname[i]+"  which is galaxy #  "+strtrim(string(fix(i)),2)+" at "+systime()
     ;Add this line in to follow the evolution of the first fit.
;     spawn,'mv 1stfit.def 1stfit_'+strtrim(strcompress(string(fix(counter))),2)+'.def'

     rename,'1stfit.','1stfitold.'
     gipsyfirst=strarr(1)
     gipsyfirst='tirific DEFFILE=tirific.def ACTION=1'
     spawn,gipsyfirst,isthere2
     
     testing1:
                                ;Let's check if it is accepted 
     get_progress,maindir+'/'+catdirname[i]+'/progress1.txt',AC1,nopoints,loops,toymodels
     IF size(log,/TYPE) EQ 7 then begin
        openu,66,log,/APPEND
        IF AC1 EQ 1 then begin
           printf,66,linenumber()+"The first estimate was accepted in this run." 
           printf,66,linenumber()+"Tirific ran through "+strtrim(strcompress(string(loops)))+" loops and produced "+strtrim(strcompress(string(toymodels)))+" models."
           printf,66,linenumber()+"Disk 1 had "+strtrim(strcompress(string(nopoints[0])))+" point sources and Disk 2 "+strtrim(strcompress(string(nopoints[1])))
        ENDIF ELSE BEGIN
           printf,66,linenumber()+"The First estimate was not accepted."
           printf,66,linenumber()+"Tirific ran through "+strtrim(strcompress(string(loops)))+" loops and produced "+strtrim(strcompress(string(toymodels)))+" models."
           printf,66,linenumber()+"Disk 1 had "+strtrim(strcompress(string(nopoints[0])))+" point sources and Disk 2 "+strtrim(strcompress(string(nopoints[1])))
        ENDELSE
        close,66
     ENDIF   
     IF AC1 EQ 1 then begin
        print,linenumber()+"The first estimate is accepted."      
     ENDIF ELSE begin
        print,linenumber()+"The first estimate is not accepted." 
     ENDELSE     
                                ;Check the number of points in the models
     check_cflux,nopoints,tirificfirst,tirificfirstvars,cfluxadjusted,log=log
     IF  cfluxadjusted then goto,notacceptedone
     
     IF testing GE 1 then begin
        maxchangeRA=ABS(0.15*catmajbeam[i])/3600.
        maxchangeDEC=ABS(0.15*catmajbeam[i])/3600.
        IF maxchangeRA LT ABS(0.5*pixelsizeRA) then maxchangeRA=ABS(pixelsizeRA)
        IF maxchangeDEC LT ABS(0.5*pixelsizeDEC) then maxchangeDEC=ABS(pixelsizeDEC)
        maxchangevel=ABS(0.5*channelwidth)
        IF maxchangeRA LT  1./3600. then maxchangeRA=1./3600.
        IF maxchangeDEC LT  1./3600. then maxchangeDEC=1./3600.
        IF maxchangevel LT 2.5 then  maxchangevel=2.5
        goto,testing1skip
     ENDIF
                               ;Let's extract some fitting parameters from the file we just fitted
     VariablesWanted=['INCL','PA','SBR','XPOS','YPOS','VSYS','CFLUX','VROT','SBR_2','RADI','Z0']
     firstfitvalues=0.
     writenewtotemplate,tirificfirst,maindir+'/'+catdirname[i]+'/1stfit.def',Arrays=firstfitvalues,VariableChange=VariablesWanted,/EXTRACT
     tmppos=where('SBR' EQ VariablesWanted)
     SBRarr=firstfitvalues[*,tmppos]
     tmppos=where('SBR_2' EQ VariablesWanted)
     SBRarr2=firstfitvalues[*,tmppos]
     tmppos=where('VROT' EQ VariablesWanted)
     VROTarr=firstfitvalues[*,tmppos]
     stringvelocities='VROT= 0. '+STRJOIN(VROTarr[1:n_elements(VROTarr)-1],' ')
     IF size(log,/TYPE) EQ 7 then begin
        openu,66,log,/APPEND
        printf,66,linenumber()+"SBRs before they are checked."
        printf,66,SBRarr,SBRarr2
        close,66
     endif
     sbr_check,tirificfirst, tirificfirstvars,sbrarr,sbrarr2,cutoff
     IF size(log,/TYPE) EQ 7 then begin
        openu,66,log,/APPEND
        printf,66,linenumber()+"SBRs after they are checked."
        printf,66,SBRarr,SBRarr2
        close,66
     endif
                                ; Let's see if the values are too close to the boundaries or not
                                ;First the inclination
     diffmax=ABS(firstfitvalues[0,0]-double(inclinput1[1]))
     diffmin=ABS(firstfitvalues[0,0]-double(inclinput1[2]))                            
     oldinc=catinc[i]
     oldpa=catpa[i]
     paraised=0.
     inclraised=0
     IF (diffmax GT 0.1*catincdev[i] and diffmin GT 0.1*catincdev[i]) then begin 
                                ; if the  value is  more than 10% of
                                ; the error removed from the boundary
                                ; values it is ok and we update the
                                ; inclination or we have high inclination
        catinc[i]=ABS(firstfitvalues[0,0])
        inclraised=1
        IF catinc[i] GT 90. then catinc[i]=90.-(catinc[i]-90.)
        tmppos=where('INCL' EQ tirificfirstvars)
        tirificfirst[tmppos]='INCL= '+strtrim(string(catinc[i]))
        tmppos=where('INCL_2' EQ tirificfirstvars)
        tirificfirst[tmppos]='INCL_2= '+strtrim(string(catinc[i]))                      
        catmaxrot[i]=W50/2./SIN(catinc[i]*!pi/180.) 
        catmaxrotdev[i]=(VSYSdiff/2.)/SIN(catinc[i]*!pi/180.)
        tmppos=where('VROT' EQ tirificfirstvars)
        tirificfirst[tmppos]=stringvelocities
        tmp=str_sep(strtrim(strcompress(stringvelocities),2),'=')
        tmppos=where('VROT_2' EQ tirificfirstvars)
        tirificfirst[tmppos]='VROT_2='+tmp[1]
        IF size(log,/TYPE) EQ 7 then begin
           openu,66,log,/APPEND
           printf,66,linenumber()+"The incl was adjusted to"+string(catinc[i])
           close,66
        ENDIF ELSE BEGIN
           print,linenumber()+"The incl was adjusted to"+string(catinc[i])
        ENDELSE
     ENDIF ELSE Begin
        ; Otherwise we will adjust the boundaries if they are insecure
        IF ceil(norings[0]*COS(firstfitvalues[0,0]*!DtoR)) GT 4 OR (norings[0] GT 6. AND catinc[i] GT 80.) then begin
                                ;If we are on the boundary we reject
                                ;the new value unless we are in very
                                ;uncertain conditions
           IF size(log,/TYPE) EQ 7 then begin
              openu,66,log,/APPEND
              printf,66,linenumber()+"The incl was on the boundary"+string(firstfitvalues[0,0])+'+'+string(inclinput1[1])+'-'+string(inclinput1[2])
              close,66
           ENDIF ELSE BEGIN
              print,linenumber()+"The incl was on the boundary"+string(firstfitvalues[0,0])+'+'+string(inclinput1[1])+'-'+string(inclinput1[2])
           ENDELSE
                                ;However if this is due to high
                                ;inclination we still want to reset vrot
           IF  (norings[0] GT 6. AND catinc[i] GT 80.) then begin
              catmaxrot[i]=W50/2./SIN(catinc[i]*!pi/180.) 
              catmaxrotdev[i]=(VSYSdiff/2.)/SIN(catinc[i]*!pi/180.)
              tmppos=where('VROT' EQ tirificfirstvars)
              tirificfirst[tmppos]=stringvelocities
              tmp=str_sep(strtrim(strcompress(stringvelocities),2),'=')
              tmppos=where('VROT_2' EQ tirificfirstvars)
              tirificfirst[tmppos]='VROT_2='+tmp[1]     
          ENDIF ELSE AC1=0
        ENDIF ELSE BEGIN
                                ; Otherwise we will adjust the boundaries if they are insecure
           inclraised=1
           catinc[i]=ABS(firstfitvalues[0,0])
           AC1=0
           IF catinc[i] GT 90. then catinc[i]=90.-(catinc[i]-90.)
           tmppos=where('INCL' EQ tirificfirstvars)
           tirificfirst[tmppos]='INCL= '+strtrim(string(catinc[i]))
           tmppos=where('INCL_2' EQ tirificfirstvars)
           tirificfirst[tmppos]='INCL_2= '+strtrim(string(catinc[i]))
           catmaxrot[i]=W50/2./SIN(catinc[i]*!pi/180.)
           catmaxrotdev[i]=(VSYSdiff/2.)/SIN(catinc[i]*!pi/180.)
           tmppos=where('VROT' EQ tirificfirstvars)
           tirificfirst[tmppos]=stringvelocities
           tmp=str_sep(strtrim(strcompress(stringvelocities),2),'=')
           tmppos=where('VROT_2' EQ tirificfirstvars)
           tirificfirst[tmppos]='VROT_2='+tmp[1]
           IF size(log,/TYPE) EQ 7 then begin
              openu,66,log,/APPEND
              printf,66,linenumber()+"Although the incl hit its boundary we adjusted to"+string(catinc[i])+"due to the small ring size. "
              close,66
           ENDIF ELSE BEGIN
              print,linenumber()+"Although the incl hit its boundary we adjusted to"+string(catinc[i])+"due to the small ring size. "
           ENDELSE
           VROTmax=catmaxrotdev[i]
           INCLmin=catinc[i]-catincdev[i]-20
                                ;We base the minimum vrot max on
                                ;inclination however the maximum of
                                ;600 is a physical limit.
           IF VROTmax LT 80 then VROTmax=80.
           IF VROTmax GT 600. then VROTmax=600.
           VROTinput1[1]=strtrim(strcompress(string(VROTmax)),2)
           INCLinput1=['INCL 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                       ' INCL_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1),$
                       string(90.),string(INCLmin),string(1.0),string(0.1),string(1.),string(0.1),'5','70','70']
        ENDELSE
     ENDELSE
                                ;when we get a new inclination we need
                                ;to reset our cutoff values.
     cutoffcorrection=SIN(75.*!DtoR)/SIN(catinc[i]*!DtoR)
     IF newinclination[0] GT 50 then cutoffcorrection=1.
     IF newinclination[0] LT 50 AND newinclination[0] GT 40  then cutoffcorrection=1.+(50-newinclination[0])*0.05
     IF cutoffcorrection GT 2.5 then cutoffcorrection=2.5
   
     IF size(log,/TYPE) EQ 7 then begin
        openu,66,log,/APPEND
        printf,66,linenumber()+"Cutoff values will be adjusted for inclination by multiplying them with "+string(cutoffcorrection)
        close,66
     ENDIF
     cutoff=cutoffor*cutoffcorrection
                                ;Then we look at the PA
     diffmax=ABS(firstfitvalues[0,1]-double(painput1[1]))
     diffmin=ABS(firstfitvalues[0,1]-double(painput1[2])) 
     IF diffmax GT 0.1*catpadev[i] and diffmin GT 0.1*catpadev[i] then begin 
        catpa[i]=firstfitvalues[0,1]
        paraised=1
        ;if not between -360 and 360 adapt
        WHILE ABS(catpa[i])/360. GT 1. DO BEGIN
           IF catpa[i] GT 360. then catpa[i]=catpa[i]-360.
           IF catpa[i] LT -360. then catpa[i]=catpa[i]+360.
        ENDWHILE
        IF catpa[i] LT 0. then catpa[i]=catpa[i]+360. 
        tmppos=where('PA' EQ tirificfirstvars)
        tirificfirst[tmppos]='PA= '+strtrim(string(catpa[i]))
        tmppos=where('PA_2' EQ tirificfirstvars)
        tirificfirst[tmppos]='PA_2= '+strtrim(string(catpa[i]))
        IF size(log,/TYPE) EQ 7 then begin
           openu,66,log,/APPEND
           printf,66,linenumber()+"The pa was adjusted to "+string(catpa[i])
           close,66
        ENDIF ELSE BEGIN
           print,linenumber()+"The pa was adjusted to "+string(catpa[i])
        ENDELSE
     ENDIF ELSE begin 
        IF shiftcentercounter EQ 1. then begin
           IF UndatedPA LT 180 then begin
              firstfitvalues[*,1]=firstfitvalues[*,1]+180.
              catpa[i]=catpa[i]+180.
              UndatedPA=UndatedPA+180. 
           endif else begin
              firstfitvalues[*,1]=firstfitvalues[*,1]-180.
              catpa[i]=catpa[i]-180.
              UndatedPA=UndatedPA-180.
           endelse
           
        ENDIF ELSE BEGIN
           IF oldinc LT 50 AND UndatedINCL LT 65 AND catinc[i] LT 55 then begin  
              
              paraised=1
              catpa[i]=firstfitvalues[0,1]
              AC1=0
              WHILE ABS(catpa[i])/360. GT 1. DO BEGIN
                 IF catpa[i] GT 360. then catpa[i]=catpa[i]-360.
                 IF catpa[i] LT -360. then catpa[i]=catpa[i]+360.
              ENDWHILE
              IF catpa[i] LT 0. then catpa[i]=catpa[i]+360. 
              tmppos=where('PA' EQ tirificfirstvars)
              tirificfirst[tmppos]='PA= '+strtrim(string(catpa[i]))
              tmppos=where('PA_2' EQ tirificfirstvars)
              tirificfirst[tmppos]='PA_2= '+strtrim(string(catpa[i]))
              IF size(log,/TYPE) EQ 7 then begin
                 openu,66,log,/APPEND
                 printf,66,linenumber()+"Although the pa hit its boundary we adjusted to"+string(catpa[i])+"because of low inclination."
                 close,66
              ENDIF ELSE BEGIN
                 print,linenumber()+"Although the pa hit its boundary we adjusted to"+string(catpa[i])+"because of low inclination."
              ENDELSE
              IF PAest LT 5 then begin
                 PAinput1=['PA 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                           ' PA_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1),$
                           string((catPA[i])+catPAdev[i]+30.),string((catPA[i])-catPAdev[i]-30),string(5),string(0.1),string(0.05),string(0.01),'3','70','70']
                 IF PAfixboun NE 'n' then begin
                    IF PAfixboun EQ 'w' then begin
                       IF catPA[i]+catPAdev[i]+30. GT 360 then PAinput1[1]='360'
                       IF catPA[i]-catPAdev[i]-30. LT 180 then PAinput1[2]='180'
                    ENDIF
                    IF PAfixboun EQ 'e' then begin
                       IF catPA[i]+catPAdev[i]+30. GT 180. then PAinput1[1]='180'
                       IF catPA[i]-catPAdev[i]-30. LT 0. then PAinput1[2]='0'
                    ENDIF
                 ENDIF
                 Writefittingvariables,tirificfirst,painput1
                 PAest++
                 goto,againPAestimate 
              ENDIF else begin
                 PAinput1=['PA 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                           ' PA_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1),$
                           string((catPA[i])+catPAdev[i]+30.),string((catPA[i])-catPAdev[i]-30),string(2.5),string(0.1),string(0.05),string(0.01),'3','70','70'] 
                 IF PAfixboun NE 'n' then begin
                    IF PAfixboun EQ 'w' then begin
                       IF catPA[i]+catPAdev[i]+30. GT 360 then PAinput1[1]='360'
                       IF catPA[i]-catPAdev[i]-30. LT 180 then PAinput1[2]='180'
                    ENDIF
                    IF PAfixboun EQ 'e' then begin
                       IF catPA[i]+catPAdev[i]+30. GT 180. then PAinput1[1]='180'
                       IF catPA[i]-catPAdev[i]-30. LT 0. then PAinput1[2]='0'
                    ENDIF
                 ENDIF
              ENDELSE
           ENDIF
        ENDELSE
        IF size(log,/TYPE) EQ 7 then begin
           openu,66,log,/APPEND
           printf,66,linenumber()+"The pa was on the boundary"+string(firstfitvalues[0,1])+'+'+string(painput1[1])+'-'+string(painput1[2])
           close,66
        ENDIF ELSE BEGIN
           print,linenumber()+"The pa was on the boundary"+string(firstfitvalues[0,1])+'+'+string(painput1[1])+'-'+string(painput1[2])
        ENDELSE
        AC1=0
     ENDELSE
     lastsbr=firstfitvalues[n_elements(firstfitvalues[*,2])-1,2]
     avlastsbr=(TOTAL(firstfitvalues[n_elements(firstfitvalues[*,2])-2:n_elements(firstfitvalues[*,2])-1,2]))/2.
     lastsbr2=firstfitvalues[n_elements(firstfitvalues[*,8])-1,2]
     avlastsbr2=(TOTAL(firstfitvalues[n_elements(firstfitvalues[*,2])-2:n_elements(firstfitvalues[*,8])-1,2]))/2.
     tmppos=where('Z0' EQ tirificfirstvars)
     tirificfirst[tmppos]='Z0= '+strtrim(string(firstfitvalues[0,10]))
     tmppos=where('Z0_2' EQ tirificfirstvars)
     tirificfirst[tmppos]='Z0_2= '+strtrim(string(firstfitvalues[0,10]))
                                ; let's check wether the center
                                ; shifted a lot since this affects everything
                                ;if there is a major change in central
                                ;position as it determines many parameters
     
     newxpos=firstfitvalues[0,3] 
     newypos=firstfitvalues[0,4]
     newvsys=firstfitvalues[0,5]
     maxchangeRA=ABS(0.15*catmajbeam[i])/3600.
     maxchangeDEC=ABS(0.15*catmajbeam[i])/3600.
     IF maxchangeRA LT ABS(0.5*pixelsizeRA) then maxchangeRA=ABS(pixelsizeRA)
     IF maxchangeDEC LT ABS(0.5*pixelsizeDEC) then maxchangeDEC=ABS(pixelsizeDEC)
     maxchangevel=ABS(0.5*channelwidth)
     IF maxchangeRA LT  1./3600. then maxchangeRA=1./3600.
     IF maxchangeDEC LT  1./3600. then maxchangeDEC=1./3600.
     IF maxchangevel LT 2.5 then  maxchangevel=2.5
                                ;if the change between this fit and
                                ;the previous one is bigger that the
                                ;limits, refit
     IF ABS(RADeg-newxpos) GT maxchangeRA OR ABS(DECDeg-newypos) GT maxchangeDEC OR ABS(catvsys[i]-newvsys) GT maxchangevel OR fixedcenter EQ 1. then begin
        ; If the center was fixed refit
        IF fixedcenter EQ 1 then begin
           IF size(log,/TYPE) EQ 7 then begin
              openu,66,log,/APPEND
              printf,66,linenumber()+"The center was not fitted."     
              close,66
           ENDIF ELSE BEGIN
              print,linenumber()+"The center was not fitted."
           ENDELSE
           fixedcenter=0
           goto,shiftcenter
        ENDIF ELSE BEGIN
                                ;if the shift is more than two beams
                                ;from the initial guess something went wrong
           IF    ABS(RADeg-newxpos) GT catmajbeam[i]/1800. OR ABS(DECDeg-newypos) GT catmajbeam[i]/1800. then begin
              IF size(log,/TYPE) EQ 7 then begin
                 openu,66,log,/APPEND
                 printf,66,linenumber()+"The center shifted more than 2 major beams. Not applying this shift."
                 close,66
              ENDIF ELSE BEGIN
                 print,linenumber()+"The center shifted than 2 major beams. Not applying this shift."
              ENDELSE
              IF paraised then begin
                 tmppos=where('PA' EQ tirificfirstvars)
                 tirificfirst[tmppos]='PA= '+strtrim(string(oldpa))
                 tmppos=where('PA_2' EQ tirificfirstvars)
                 tirificfirst[tmppos]='PA_2= '+strtrim(string(oldpa))
              ENDIF
              IF inclraised then begin
                 tmppos=where('INCL' EQ tirificfirstvars)
                 tirificfirst[tmppos]='INCL= '+strtrim(string(oldinc))
                 tmppos=where('INCL_2' EQ tirificfirstvars)
                 tirificfirst[tmppos]='INCL_2= '+strtrim(string(oldinc))
              ENDIF
              fixedcenter=1
              tryone=0.
              plus2beamshift++
                                ; if the center has tried to shift by 2 beams
                                ; more than 10 times we are just going
                                ; to keep this center
              IF plus2beamshift GT 10 then begin
                 keepcenter=1.
                 fixedcenter=0.
                 IF size(log,/TYPE) EQ 7 then begin
                    openu,66,log,/APPEND
                    printf,66,linenumber()+"The center keeps trying to change unreasonably. Accepting the current center and fixing it."
                    close,66
                 ENDIF ELSE BEGIN
                    print,linenumber()+"The center keeps trying to change unreasonably. Accepting the current center and fixing it."
                 ENDELSE
              ENDIF
              goto,shiftcenter
           ENDIF
                                ;If the shift in center is jumping
                                ;back and forth between two poistions
                                ;accept the first/last position
           IF ABS(oldRA-newxpos) LT maxchangeRA AND ABS(oldDEC-newypos) LT maxchangeDEC AND ABS(oldVSYS-newvsys) LT maxchangevel then begin
              IF size(log,/TYPE) EQ 7 then begin
                 openu,66,log,/APPEND
                 printf,66,linenumber()+"The center shifted back to the previous fit. Accepting this center and fixing it."
                 close,66
              ENDIF ELSE BEGIN
                 print,linenumber()+"The center shifted back to the previous fit. Accepting this center and fixing it."
              ENDELSE
              keepcenter=1.
              goto,shiftback
           ENDIF
                                ;If the shift was to large try again
           IF size(log,/TYPE) EQ 7 then begin
              openu,66,log,/APPEND
              printf,66,linenumber()+"The center shifted too much trying again with new center."
              printf,66,linenumber()+"The RA has shifted from "+string(RAdeg)+" to "+string(newxpos)+" which is a difference of"$
                     +string(ABS(RADeg-newxpos))+" needed ="+string(maxchangeRA)
              printf,66,linenumber()+"The DEC has shifted from "+string(DECdeg)+" to "+string(newypos)+" which is a difference of"$
                     +string(ABS(DECDeg-newypos))+" needed ="+string(maxchangeDEC)
              printf,66,linenumber()+"The systemic has shifted from "+string(catvsys[i])+" to "+string(newvsys)+" which is a difference of"$
                     +string(ABS(catvsys[i]-newvsys))+" needed ="+string(maxchangevel)
              close,66
           ENDIF ELSE BEGIN
              print,linenumber()+"The center shifted too much trying again with new center."
           ENDELSE
        ENDELSE
        oldRA=RADeg
        oldDEC=DECdeg
        oldVSYS=catvsys[i]
        RADeg=newxpos
        DECdeg=newypos
        catvsys[i]=newvsys
        VariablesWanted=['XPOS','XPOS_2','YPOS','YPOS_2','VSYS','VSYS_2']
        firstfitvalues=0.
        writenewtotemplate,tirificfirst,maindir+'/'+catdirname[i]+'/1stfit.def',Arrays=firstfitvalues,VariableChange=VariablesWanted
        tryone=0.
        goto,shiftcenter
     ENDIF  ELSE BEGIN
        IF newxpos GT RAboundeg[0] AND newxpos LT RAboundeg[1] AND newypos GT DECboundeg[0] AND newypos LT DECboundeg[1] AND $
           newvsys GT ROTboun[0] AND newvsys LT ROTboun[1] AND testing LT 1 then begin 
           shiftback:
           RADeg=newxpos
           DECdeg=newypos
           catvsys[i]=newvsys
           if fixedcenter EQ 1. then begin
              IF size(log,/TYPE) EQ 7 then begin
                 openu,66,log,/APPEND
                 printf,66,linenumber()+"The center was accepted because it was fixed."
                 close,66
              ENDIF
              goto,shiftcenter
           endif        
           ;let's see if the model has the right size
           IF secondtime then norings=newrings else get_newringsv8,SBRarr,SBRarr2,cutoff,newrings
           if newrings LT 3 then begin
              IF size(log,/TYPE) EQ 7 then begin
                 openu,66,log,/APPEND
                 printf,66,linenumber()+"Cause the newrings was "+strtrim(string(fix(newrings)),2)+" we set it to 3."
                 close,66
              ENDIF         
              newrings=3
           ENDIF
                                ;see if the newsize is not too big
           IF newrings GT maxrings then begin
              IF size(log,/TYPE) EQ 7 then begin
                 openu,66,log,/APPEND
                 printf,66,linenumber()+"Cause the newrings was "+strtrim(string(fix(newrings)),2)+" we set it to maxrings."
                 close,66
              ENDIF       
              newrings=maxrings 
           ENDIF
           
           IF size(log,/TYPE) EQ 7 then begin
              openu,66,log,/APPEND
              printf,66,linenumber()+"We found this as rings="+strtrim(string(fix(norings[0])),2)+"  new="+strtrim(string(fix(newrings)),2)
              close,66
           ENDIF
           IF NOT sofiafail then begin
              IF newrings GT sofiarings+2 OR newrings LT sofiarings-2 then begin        
                 IF size(log,/TYPE) EQ 7 then begin
                    openu,66,log,/APPEND
                    printf,66,linenumber()+"The new amount of rings ("+strtrim(string(fix(newrings)),2)+") deviates too much from the sofia estimate ("+string(sofiarings)+")."
                    close,66
                 ENDIF
                 newrings=norings[0] 
              ENDIF
           ENDIF
           IF newrings LT norings[0] AND constring NE 1 then begin 
              IF newrings LE forcedring then begin
                 IF norings[0] GT forcedring then begin
                    newrings=forcedring
                    constring=1
                    prevmodification=0.
                 ENDIF ELSE BEGIN
                    IF size(log,/TYPE) EQ 7 then begin
                       openu,66,log,/APPEND
                       printf,66,linenumber()+"As the cut would modify the rings to less than the maximum forced addition of a ring we do not apply the cut."
                       printf,66,linenumber()+"new rings "+strtrim(string(fix(newrings)),2)+" old rings "+strtrim(string(fix(norings[0])),2)+"forcedring="+string(forcedring)
                       Close,66
                    ENDIF
                    goto,nocutfirst
                 ENDELSE
              ENDIF    
                                ;check that we do not go back and forth
              IF prevmodification EQ 1 then begin
                 IF newrings GT norings[0]-1 then begin
                    prevmodification=-1
                 ENDIF else BEGIN
                    IF newrings EQ norings[0]-1 then begin
                       IF size(log,/TYPE) EQ 7 then begin
                          openu,66,log,/APPEND
                          printf,66,linenumber()+"As the last modification was the addition of a ring we will fix this ring number."
                          printf,66,linenumber()+"new rings "+strtrim(string(fix(newrings)),2)+" old rings "+strtrim(string(fix(norings[0])),2)
                          Close,66
                       ENDIF
                       constring=1
                       secondtime=1
                       newrings=newrings+1
                    ENDIF 
                    prevmodification=-2
                ENDELSE
              ENDIF 
              oldrings=norings[0]
                                ;If no change then go on with the fitting process
              IF newrings EQ norings[0] then goto,nocutfirst else norings[0]=newrings
                                ;check if we fitted this size before
              tmp=WHERE(prevrings EQ newrings)
              IF tmp[0] NE -1 AND newrings NE oldrings then begin
                 IF size(log,/TYPE) EQ 7 then begin
                    openu,66,log,/APPEND
                    printf,66,linenumber()+"Cause the newrings was "+strtrim(string(fix(newrings)),2)+" which we have fitted before we gonna fix the rings at this value."
                    close,66
                 ENDIF     
                 secondtime=1
              ENDIF
              IF size(log,/TYPE) EQ 7 then begin
                 openu,66,log,/APPEND
                 printf,66,linenumber()+"We cut a ring! lastcut"+string(lastcutrings)+" lastadd "+string(lastaddrings)+" new "+strtrim(string(fix(newrings)),2)+" old "+string(oldrings)
                 Close,66
              ENDIF ELSE begin            
                 print,linenumber()+"We cut a ring! lastcut"+string(lastcutrings)+" lastadd "+string(lastaddrings)+" new "+strtrim(string(fix(newrings)),2)+" old "+string(oldrings)
              ENDELSE
              lastcutrings=norings[0]
              fluxadjust=0.             
              tmppos=where('SBR' EQ tirificfirstvars)
              tirificfirst[tmppos]='SBR='+STRJOIN(SBRarr[0:norings-1],' ')
              tmppos=where('SBR_2' EQ tirificfirstvars)
              tirificfirst[tmppos]='SBR_2='+STRJOIN(SBRarr2[0:norings-1],' ')
              countsbr++
              overwrite=0.
              prevmodification=-1
              ringmodifier=ringmodifier+1
              prevrings=[prevrings,norings[0]]
              ;refit
              goto,sbrshift
           ENDIF
           nocutfirst:
           overwrite=0
           ;if the last rings are really bright force add an addition
           IF (SBRarr[newrings-2] GT 7.5*cutoff[newrings-2] OR SBRarr2[newrings-2] GT 7.5*cutoff[newrings-2]) AND newrings GT norings[0] AND prevmodification EQ -1 then begin
              prevmodification=0.
              IF size(log,/TYPE) EQ 7 then begin
                 openu,66,log,/APPEND
                 printf,66,linenumber()+"We force added a ring to the profile."
                 printf,66,linenumber()+"Newrings"+strtrim(string(fix(newrings)),2)+", Rings"+strtrim(string(fix(norings[0])),2)
                 Close,66
              ENDIF
              IF newrings GT forcedring then forcedring=newrings
              overwrite=1
           ENDIF
                                ;If we previously subtracted we will not add
           IF prevmodification EQ -1 AND overwrite NE 1 AND newrings GT norings[0] then begin
              IF size(log,/TYPE) EQ 7 then begin
                 openu,66,log,/APPEND
                 printf,66,linenumber()+"We wanted to add a ring to the model but the previous step was a subtraction."
                 Close,66
              ENDIF
           ENDIF
                                ;If we want to add a ring we need to
                                ;set some parameters and make estimate
           IF newrings GT norings[0] AND (prevmodification NE -1 OR overwrite EQ 1) AND secondtime NE 1 then begin        
              tmp=WHERE(prevrings EQ newrings)
              IF tmp[0] NE -1 AND newrings NE norings[0] then begin
                 IF size(log,/TYPE) EQ 7 then begin
                    openu,66,log,/APPEND
                    printf,66,linenumber()+"Cause the newrings was "+strtrim(string(fix(newrings)),2)+" which we have fitted before we gonna fix the rings at this value."
                    close,66
                 ENDIF     
                 secondtime=1
              ENDIF
              norings[0]=newrings
              lastaddrings=norings[0]
              IF size(log,/TYPE) EQ 7 then begin
                 openu,66,log,/APPEND
                 printf,66,linenumber()+"We added a ring to the model. "+string(secondtime)
                 Close,66
              ENDIF ELSE BEGIN
                 print,linenumber()+"We added a ring to the model."
              ENDELSE
              prevmodification=1.
              ringmodifier=ringmodifier+1
              countsbr++
              prevrings=[prevrings,norings[0]]
              goto,sbrshift
           ENDIF
        ENDIF ELSE BEGIN
           ; If we really go crazy then just end the fitting process
           IF size(log,/TYPE) EQ 7 then begin
              openu,66,log,/APPEND
              printf,66,linenumber()+"The fit diverged out of the boundaries set by the parameterization. That can't be good."
              printf,66,linenumber()+"Boundaries are; DEC current="+string(newypos)+"DEC Boundaries"+string(DECboundeg[0])+","+string(DECboundeg[1])
              printf,66,linenumber()+"RA current="+string(newxpos)+"RA Boundaries"+string(RAboundeg[0])+","+string(RAboundeg[1])
              printf,66,linenumber()+"V_sys current="+string(newvsys)+"V_sys Boundaries"+string(ROTboun[0])+","+string(ROTboun[1])
              close,66
           ENDIF ELSE BEGIN
              print,linenumber()+"The fit diverged out of the boundaries set by the parameterization. That can't be good."
           ENDELSE  
           openu,1,outputcatalogue,/APPEND
           printf,1,format='(A60,A80)',catDirname[i],'This galaxy diverged out of the set boundaries.'
           close,1   
           bookkeeping=0
           goto,finishthisgalaxy
        ENDELSE
     ENDELSE
     IF AC1 NE 1 And tryone LT 10. then begin
                                ;If not accepted let's try again
                                
        
        IF size(log,/TYPE) EQ 7 then begin
           openu,66,log,/APPEND
           printf,66,linenumber()+"The fit is not accepted trying once more."
           close,66
        ENDIF ELSE BEGIN
           print,linenumber()+"The fit is not accepted trying once more."
        ENDELSE
                                ;but from the newvalues
        writenewtotemplate,tirificfirst,maindir+'/'+catdirname[i]+'/1stfit.def'
        tryone=tryone+1.
        goto,notacceptedone
     endif 
                                ;let's do a final check with a
                                ;smoothed rotation curve when all is well
     IF NOT smoothrotation then begin 
        smoothrotation=1
        VROTarr=firstfitvalues[*,7]
        SBRarr=(firstfitvalues[*,2]+firstfitvalues[*,8])/2.
        VROTarr[0]=0.
        vmaxdev=MAX([30,7.5*channelwidth*(1.+vresolution)])
        verror=MAX([5.,channelwidth*(1.+vresolution)/SQRT(sin(catinc[i]*!DtoR))])
        IF VROTarr[1] LT 120. AND VROTarr[2] LT 120. then begin
           parameterreguv87,VROTarr,SBRarr, firstfitvalues[*,9],/REVERSE,fixedrings=norings[0]-fix(norings[0]*3./4.),difference=verror,cutoff=cutoff,arctan=0.,order=polorder,max_par=VROTmax,min_par=channelwidth
        ENDIF ELSE Begin
           parameterreguv87,VROTarr,SBRarr, firstfitvalues[*,9],/REVERSE,fixedrings=norings[0]-fix(norings[0]*3./4.),difference=verror,cutoff=cutoff,arctan=0.,order=polorder,max_par=VROTmax,min_par=channelwidth,/NOCENTRAL
        ENDELSE
        tmp0check=WHERE(VROTarr LT 0)
        IF tmp0check[0] NE -1 then VROTarr[tmp0check]=3.*channelwidth
        tmppos=where('VROT' EQ tirificfirstvars)
        tirificfirst[tmppos]='VROT= 0. '+STRJOIN(strtrim(string(VROTarr[1:n_elements(VROTarr[*])-1]),2),' ')
        tmppos=where('VROT_2' EQ tirificfirstvars)
        tirificfirst[tmppos]='VROT_2= 0. '+STRJOIN(strtrim(string(VROTarr[1:n_elements(VROTarr[*])-1]),2),' ')
        vrotinputoriginal=vrotinput1
        VROTinput1=['VROT 2:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                    ' VROT_2 2:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1),$
                    string(VROTmax),string(VROTmin),string(channelwidth),string(0.1*channelwidth),string(channelwidth),string(0.1*channelwidth),'3','70','70', ' ']
        IF size(log,/TYPE) EQ 7 then begin
           openu,66,log,/APPEND
           printf,66,linenumber()+"We have smoothed the rotation curve in the first fit."
           close,66
        ENDIF ELSE BEGIN
           print,linenumber()+"We have smoothed the rotation curve in the first fit."
        ENDELSE

        goto,smoothingone
     ENDIF
     IF size(log,/TYPE) EQ 7 then begin
        openu,66,log,/APPEND
        printf,66,linenumber()+"We finished the first fit of "+catDirname[i]+" at "+systime()
        printf,66,linenumber()+"We have rerun "+string(counter)+" times."
        printf,66,linenumber()+"We have modified the rings "+string(countsbr)+" times."
        IF AC1 EQ 0 then  printf,66,linenumber()+"The fit was not accepted."
        IF AC1 EQ 1 then  printf,66,linenumber()+"The fit was accepted."
        close,66
     ENDIF 
     testing1skip:
                                ;Reading out the fit values 
     Basicinfovars=['XPOS','YPOS','VSYS','PA','INCL','VROT']
     writenewtotemplate,tirificfirst,maindir+'/'+catdirname[i]+'/1stfit.def',Arrays=Basicinfovalues,VariableChange=Basicinfovars,Variables=tirificfirstvars,/EXTRACT
     RAhr=Basicinfovalues[0,0]
     RAdiff=maxchangeRA*3600./15.
     DEChr=Basicinfovalues[0,1]
     DECdiff=maxchangeDEC*3600.

                               ;If the cube is optimized then make a model at full resolution 
     IF optimized then begin
        tmp=str_sep(strtrim(strcompress(currentfitcube),2),'_opt')
        IF n_elements(tmp) EQ 2 then begin
           currentfitcube=tmp[0]
           tmppos=where('INSET' EQ tirificfirstvars)
           tirificfirst[tmppos]='INSET=  '+strtrim(strcompress(string(currentfitcube+'.fits')))
           writenewtotemplate,tirificfirst,maindir+'/'+catdirname[i]+'/1stfit.def'
           tmppos=where('LOOPS' EQ tirificfirstvars)
           tirificfirst[tmppos]='LOOPS=  0'
           openw,1,maindir+'/'+catdirname[i]+'/tirific.def'
           for index=0,n_elements(tirificfirst)-1 do begin
              printf,1,tirificfirst[index]
           endfor
           close,1                       
           IF testing GE 1 then goto,skiporver        
           rename,'1stfit.','1stfit_opt.'
           gipsyfirst=strarr(1)
           gipsyfirst='tirific DEFFILE=tirific.def ACTION=1'
           spawn,gipsyfirst,isthere2
           skiporver:
           currentfitcube=currentfitcube+'_opt'
        ENDIF
     ENDIF





     tmpcube=readfits(maindir+'/'+catdirname[i]+'/1stfit.fits',hedtmp1stcube)    
                                ;and we make a pv-diagram based on these parameters
     IF optimized then begin
        extract_pv,nooptcube,nopthed,Basicinfovalues[0,3],xv,center=[RAhr,DEChr],xvheader=new_header
        writefits,maindir+'/'+catdirname[i]+'/'+noptname[0]+'_1_xv.fits',xv,new_header
     ENDIF ELSE BEGIN
        extract_pv,dummy,hed,Basicinfovalues[0,3],xv,center=[RAhr,DEChr],xvheader=new_header
        writefits,maindir+'/'+catdirname[i]+'/'+catcubename[i]+'_1_xv.fits',xv,new_header
     ENDELSE
     extract_pv,tmpcube,hedtmp1stcube,Basicinfovalues[0,3],xv,center=[RAhr,DEChr],xvheader=new_header
     writefits,maindir+'/'+catdirname[i]+'/1stfit_xv.fits',xv,new_header
                                ;We check that the model actually has enough flux to be reasonable
     hedtmp1st=hedtmp1stcube
     tmpix=WHERE(tmpcube GT catnoise[i])
     IF tmpix[0] EQ -1 then begin
        openu,1,outputcatalogue,/APPEND
        printf,1,format='(A60,A12,A80)', catDirname[i],AC1,'The first fit does not have flux above the noise level, this means a mis fit.'
        close,1
        IF size(log,/TYPE) EQ 7 then begin
           openu,66,log,/APPEND
           printf,66,linenumber()+"No flux in final fit, aborting."
           close,66
        ENDIF 
        bookkeeping=0
        goto,finishthisgalaxy
     endif
     IF size(log,/TYPE) EQ 7 then begin
        openu,66,log,/APPEND
        printf,66,linenumber()+'we use this mask for the moment maps'+catmaskname[i]
        close,66
     ENDIF 
     mask=readfits(maindir+'/'+catdirname[i]+'/'+catmaskname[i]+'.fits',hedmask,/NOSCALE)        
                                ;mask the data cube
     tmpmask=fltarr(n_elements(tmpcube[*,0,0]),n_elements(tmpcube[0,*,0]),n_elements(tmpcube[0,0,*]))
     tmpmask[WHERE(mask GT 0.)]=tmpcube[WHERE(mask GT 0.)]
     momentsv2,tmpmask,tmpmap,hedtmp1st,0.
     writefits,maindir+'/'+catdirname[i]+'/1stfit_mom0.fits',tmpmap,hedtmp1st
     hedtmp1stv=hedtmp1stcube
     momentsv2,tmpmask,tmpmapv,hedtmp1stv,1.
     writefits,maindir+'/'+catdirname[i]+'/1stfit_mom1.fits',tmpmapv,hedtmp1stv
     getDHI,tmpmap,hedtmp1st,Basicinfovalues[0,3],[RAhr,DEChr,Basicinfovalues[0,4]],DHI
     totflux=[TOTAL(tmpcube[tmpix])/pixperbeam,(TOTAL(2.*cutoff[0:norings[0]-1]))/(n_elements(tmpix)/pixperbeam)]
     VSYSdiff=maxchangevel
     HIMASS=2.36E5*catDistance[i]^2*totflux*ABS(channelwidth)
     convertradec,RAhr,DEChr
 
                                ;write our general overview parameters
     openu,1,basicinfofile,/APPEND
     printf,1,"#After the first fit"
     printf,1,format='(A25,A25,A20,A20,A20,A20,A20,A20,A20,A20,A20,A20,A20)',string(RAhr+'+/-'+strtrim(strcompress(string(RAdiff,format='(F6.1)')),2)),$
            string(DEChr+'+/-'+strtrim(strcompress(string(DECdiff,format='(F6.1)')),2)),$
            string(strtrim(strcompress(string(Basicinfovalues[0,2])),2)+'+/-'+strtrim(strcompress(string(VSYSdiff,format='(F6.1)')),2)),$
            string(strtrim(strcompress(string(Basicinfovalues[0,3])),2)+'+/-'+strtrim(strcompress(string(catPAdev[i],format='(F6.1)')),2)),$
            string(strtrim(strcompress(string(Basicinfovalues[0,4])),2)+'+/-'+strtrim(strcompress(string(catincdev[i],format='(F6.1)')),2)),$
            string(strtrim(strcompress(string(Basicinfovalues[n_elements(Basicinfovalues[*,5])-1,5])),2)+'+/-'+strtrim(strcompress(string(catmaxrotdev[i],format='(F6.1)')),2)),$
            string(strtrim(strcompress(string(W50)),2)+'+/-'+strtrim(strcompress(string(channelwidth,format='(F6.1)')),2)),$
            string(strtrim(strcompress(string(Totflux[0],format='(E10.3)')),2)),$
            string(strtrim(strcompress(string(DHI,format='(F8.1)')),2)),$
            string(strtrim(strcompress(string(catDistance[i])),2)),$
            string(strtrim(strcompress(string(HIMass[0],format='(E10.3)')),2)),$
            string(strtrim(strcompress(string(convertskyanglefunction(DHI,catDistance[i]),format='(F8.1)')),2))
     close,1
                       
     IF finishafter EQ 1 then begin
        openu,1,outputcatalogue,/APPEND
        printf,1,format='(A60,A12,A80)', catDirname[i],AC1,'You have chosen to skip the fitting process after the first fit'
        close,1
        bookkeeping=bookkeeping+0.5
        goto,finishthisgalaxy
     ENDIF
     prevmodification=0
;******************************************This is the end of the first fit******************************

     sigmapa1=0.
     sigmapa2=0.
     sigmaincl1=0.
     sigmaincl2=0.
     sigmarot=0.
     velfixrings=1
     tmppos=where('INSET' EQ tirificsecondvars)
     tirificsecond[tmppos]='INSET=  '+strtrim(strcompress(string(currentfitcube+'.fits')))
                                ;then open the previous fit 
     firstfitvaluesnames=0.
     writenewtotemplate,tirificsecond,maindir+'/'+catdirname[i]+'/1stfit.def',Arrays=firstfitvalues,VariableChange=firstfitvaluesnames,Variables=tirificsecondvars
     tmppos=where('RADI' EQ firstfitvaluesnames)
     RADarr=firstfitvalues[*,tmppos]
     tmppos=where('VROT' EQ firstfitvaluesnames)
     VROTarr=dblarr(n_elements(firstfitvalues[*,tmppos]))
     Vfrom1=dblarr(n_elements(firstfitvalues[*,tmppos]))
     Vfrom1=firstfitvalues[*,tmppos]
     half=fix(n_elements(firstfitvalues[*,tmppos])/2.)
     VROTarr[*]=firstfitvalues[*,tmppos]
     tmppos=where('SBR' EQ firstfitvaluesnames)
     SBRarr=firstfitvalues[*,tmppos]
     tmppos=where('INCL' EQ firstfitvaluesnames)
     INCLang=firstfitvalues[*,tmppos]
     tmppos=where('PA' EQ firstfitvaluesnames)
     PAang=firstfitvalues[*,tmppos]
     VROTarr2=VROTarr
     tmppos=where('SBR_2' EQ firstfitvaluesnames)
     SBRarr2=firstfitvalues[*,tmppos]
     tmppos=where('INCL_2' EQ firstfitvaluesnames)
     INCLang2=firstfitvalues[*,tmppos]
     IF catinc[i] GT 80. then begin
        INCLang=catinc[i]
        INCLang2=catinc[i]
     ENDIF
     tmppos=where('PA_2' EQ firstfitvaluesnames)
     PAang2=firstfitvalues[*,tmppos]    
     tmppos=where('NUR' EQ firstfitvaluesnames)
     norings=firstfitvalues[0,tmppos]
     IF norings LE 4 then begin     
        IF norings GE 3 and finishafter NE 1.1 then begin
           finishafter=1.1 
           norings[0]=(norings[0]-1)*2.
           maxrings=maxrings*2.
           noringspix=norings[0]*(catmajbeam[i]/2.)/(ABS(sxpar(hedmap,'cdelt1'))*3600.)
           rad=[0.,((findgen(maxrings))*(catmajbeam[i]/2.)+catmajbeam[i]/10.)]
           calc_edge,catnoise[i],rad,[catmajbeam[i],catminbeam[i]],cutoffor
           cutoff=cutoffor*cutoffcorrection
           tmp=rad[0:norings[0]-1]
           rad=tmp
           tmp=1
           interpolate,VROTarr,RADarr,newradii=rad,output=tmp
           VROTarr=tmp
           VROTarr2=tmp
           tmp=1
           interpolate,SBRarr,RADarr,newradii=rad,output=tmp
           SBRarr=tmp
           tmp=1
           interpolate,SBRarr2,RADarr,newradii=rad,output=tmp
           SBRarr2=tmp
           RADarr=rad
           tmppos=where('RADI' EQ tirificsecondvars)
           tirificsecond[tmppos]='RADI= '+STRJOIN(strtrim(strcompress(string(RADarr))),' ')
           tmppos=where('SBR' EQ tirificsecondvars)
           tirificsecond[tmppos]='SBR= '+STRJOIN(strtrim(strcompress(string(SBRarr))),' ')
           tmppos=where('SBR_2' EQ tirificsecondvars)
           tirificsecond[tmppos]='SBR_2= '+STRJOIN(strtrim(strcompress(string(SBRarr2))),' ')
           tmppos=where('NUR' EQ tirificsecondvars)
           tirificsecond[tmppos]='NUR= '+STRJOIN(strtrim(strcompress(string(norings[0]))),' ')
          
        endif else begin
           openu,1,outputcatalogue,/APPEND
           printf,1,format='(A60,A12,A80)',catDirname[i],AC1,'The first fit model is too small to fit variations.'
           close,1
           IF size(log,/TYPE) EQ 7 then begin
              openu,66,log,/APPEND
              printf,66,linenumber()+"Finished "+catDirname[i]+" which is galaxy #  "+strtrim(string(fix(i)),2)+" at "+systime()
              close,66
           ENDIF 
           bookkeeping=0
           goto,finishthisgalaxy
        ENDELSE
     ENDIF
;Let's see how many of the inner rings we want to fix
     tmppos=where('VSYS' EQ firstfitvaluesnames)
     vsys=firstfitvalues[0,tmppos]
     levels=(sbrarr+sbrarr2)/2.*1000.
     Columndensity,levels,vsys,[catmajbeam[i],catminbeam[i]],/ARCSQUARE
     tmp=WHERE(levels GT 2E20)
     IF tmp[0] NE -1 then innerfix=floor(tmp[n_elements(tmp)-1]/1.5)-1. else innerfix=4
     IF innerfix LT 4 OR innerfix GE norings[0] OR finishafter EQ 1.1 then innerfix=4


     IF sofiafail then begin
        lowring=3
        highring=maxrings
        IF finishafter EQ 1.1 then highring=(maxrings-2)*2
     ENDIF else begin
        lowring=norings[0]-2
        highring=norings[0]+2
     ENDELSE
     tmppos=where('VROT' EQ tirificsecondvars)
     tirificsecond[tmppos]='VROT= 0.'+STRJOIN(strtrim(strcompress(string(VROTarr[1:n_elements(VROTarr)-1]))),' ')
     tmppos=where('VROT_2' EQ tirificsecondvars)
     tirificsecond[tmppos]='VROT_2= 0. '+STRJOIN(strtrim(strcompress(string(VROTarr[1:n_elements(VROTarr)-1]))),' ')
     tmppos=where('DISTANCE' EQ tirificsecondvars)
     tirificsecond[tmppos]='DISTANCE='+strtrim(strcompress(string(catDistance[i])))
                                ;check the sbr
     sbr_check,tirificsecond, tirificsecondvars,sbrarr,sbrarr2,cutoff
  
     IF norings[0] LT 4 then begin
        norings[0]=4.
        changeradii,tirificsecond,norings[0]
        IF size(log,/TYPE) EQ 7 then begin
           openu,66,log,/APPEND
           printf,66,linenumber()+"We added a ring to the model cause there were too few the new number of rings = "+strtrim(string(fix(norings[0])),2)
           Close,66
        ENDIF ELSE BEGIN
           print,linenumber()+"We added a ring to the model cause there were too few the new number of rings = "+strtrim(string(fix(norings[0])),2)
        ENDELSE
     ENDIF
    
     tmppos=where('VSYS' EQ firstfitvaluesnames)
     catvsys[i]=firstfitvalues[0,tmppos]
     newrings=norings[0]
                                ;Let's set some limits for the
                                ;first make the best case scenario arrays
                                ;INCL
     INCLtir=(TOTAL(INCLang)+TOTAL(INCLang2))/(n_elements(INCLang)+n_elements(INCLang2))
     INCLmax=INCLtir+30
     IF INCLmax GT 90. then INCLmax=90
     INCLmin=INCLtir-30
     IF INCLmin LT 0.1 then INCLmin=0.1
     INCLinput1=['INCL 1:'+strtrim(strcompress(string(innerfix,format='(F7.4)')),1)+' INCL_2 1:'+strtrim(strcompress(string(innerfix,format='(F7.4)')),1),string(INCLmax),string(INCLmin),string(1.),string(0.1),string(0.5),string(0.1),'3','70','70']
                                ;PA
     PAtir=(TOTAL(PAang)+TOTAL(PAang2))/(n_elements(PAang)+n_elements(PAang2))
     PAinput1=['PA 1:'+strtrim(strcompress(string(innerfix,format='(F7.4)')),1)+' PA_2 1:'+strtrim(strcompress(string(innerfix,format='(F7.4)')),1),string(PAtir+40),string(PAtir-40),string(0.5),string(0.1),string(1.),string(0.1),'3','70','70']
                                ;VROT
     VROTmax=MAX(Vfrom1[1:n_elements(Vfrom1)-1],Min=VROTmin)
     VROTmax=VROTmax*1.5
     VROTmin=VROTmin/4.
                                ;the minimum can never be less than
                                ;the channelwidth
     IF VROTmin LT channelwidth then VROTmin=channelwidth
                                ; The maximum should not be less than
                                ; 80 /sin(inclination) as it becomes
                                ; more uncartain at lower inclination
     IF VROTmax LT 60. then VROTmax=60.
     IF vrotmax GT 600. then VROTmax=600.
                                ; See how much of the rotation curve we want to fit as a slope
     get_newringsv8,SBRarr,SBRarr2,3.5*cutoff,velconstused
     IF norings[0] GT 8 AND not finishafter EQ 2.1 then velconstused=velconstused-1
     set_vrotv6,vrotinput1,VROTarr,velconstused,vrotmax,vrotmin,norings,channelwidth,avinner=avinner,centralexclude=centralexclude,finish_after=finishafter
                                ;set the surface brightness values
     set_sbr,SBRinput1,SBRinput2,SBRinput3,SBRinput4,SBRinput5,SBRinput6,SBRarr,cutoff,norings,finishafter,/initial
                                ;SDIS
     SDISinput1=['SDIS 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                 ' SDIS_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1),$
                 '25','5','1','0.1','0.5','0.05','3','70','70']    
                                ;Z0
     Z0input1=['Z0 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
               ' Z0_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1),$
               string(convertskyanglefunction(0.5,double(catDistance[i]),/PHYSICAL)),$
               string(convertskyanglefunction(0.075,double(catDistance[i]),/PHYSICAL)),$
               string(convertskyanglefunction(0.005,double(catDistance[i]),/PHYSICAL)),$
               string(convertskyanglefunction(0.0005,double(catDistance[i]),/PHYSICAL)),$
               string(convertskyanglefunction(0.5,double(catDistance[i]),/PHYSICAL)),$
               string(convertskyanglefunction(0.0005,double(catDistance[i]),/PHYSICAL)),'3','70','70']
     IF catinc[i] GT 80 then begin
        Z0input1[1]=string(MAX([convertskyanglefunction(1.0,double(catDistance[i])),catmajbeam[i]/2.]))
     ENDIF
     IF catDistance[i] EQ 1. then begin
                                ;ok if we do not know the distance
                                ;then let us assume each disk is about
                                ;30 kpc which means that
                                ;0.5=norings[0]/60. and so on
        IF catinc[i] GT 80 then begin
           Z0input1[1]=string(MAX([catmajbeam[i]*norings[0]/60.,catmajbeam[i]/2.]))
        ENDIF ELSE  Z0input1[1]=string(catmajbeam[i]*norings[0]/60.)
        Z0input1[2]='0.'
        Z0input1[3]=string(-1*catmajbeam[i]*norings[0]/6E4)
        Z0input1[4]=string(catmajbeam[i]*norings[0]/6E5)
        Z0input1[5]=string(catmajbeam[i]*norings[0]/60.)
        Z0input1[6]=string(catmajbeam[i]*norings[0]/6E5)
     ENDIF
                                ;And then make the other string input variables with the same
                                ;parameters and an additional string
                                ;where we can put slope fitting for
                                ;low SBR rings
     INCLinput2=[INCLinput1,' '] 
     INCLinput3=[INCLinput1,' ']
                                ;If we have a high inclination we put stringent limits on the inner
                                ;inclination
     IF catinc[i] GT 80 then begin
        INCLinput1[1]=catinc[i]+catincdev[i]
        IF INCLinput1[1] GT 90. then INCLinput1[1]=90
        INCLinput1[2]=catinc[i]-catincdev[i] 
        IF INCLinput1[2] LT 0.1 then INCLinput1[2]=0.1
     ENDIF
                                ;PA
     PAinput2=[PAinput1,' ']
     PAinput3=[PAinput1,' ']
                                ;IF we have a decent amount of rings
                                ;and we are not fitting a flat disk
                                ;with half beams than we let the rings
                                ;beyond 4 free                          
     IF norings[0] GT 4 AND finishafter NE 1.1 then begin
                                ;And some different fitting parameters for the free rings
        INCLinput2[5]='0.5'
        INCLinput2[3]='1.0'
        INCLinput3[5]='0.5'
        INCLinput3[3]='1.0'
        PAinput2[5]='2.0'
        PAinput2[3]='0.5'
        PAinput3[5]='2.0'
        PAinput3[3]='0.5'
                                                          ;update INCL and PA
        set_warp_slopev2,SBRarr,SBRarr2,cutoff,INCLinput2,PAinput2,INCLinput3,PAinput3,norings,log=log,innerfix=innerfix
     endif else begin
                                ;Otherwise we just fit  a single flat disk
        INCLinput1[0]='INCL 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                      ' INCL_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)
        PAinput1[0]='PA 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                    ' PA_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)
        SBRinput1[2]=strtrim(strcompress(string(cutoff[n_elements(cutoff)-1]/4.,format='(E12.5)')),1)
     endelse  
                       
                                ;If the first fit is accepted we only change the central position minimally
     IF AC1 EQ 1 then begin     
        xposinput1=[ ' XPOS 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                     ' XPOS_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1), '360','0','1E-5','1E-6','1E-4','1E-6','3','70','70']
        yposinput1=[ ' YPOS 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                     ' YPOS_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1), '90','-90','1E-5','1E-6','1E-4','1E-6','3','70','70']
        vsysinput1=[ ' VSYS 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                     ' VSYS_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1),strtrim(strcompress(string(catvsys[i]+100.)),1),$
                     strtrim(strcompress(string(catvsys[i]-100.)),1),'0.5','0.01','2','0.5','3','70','70']
        IF norings[0] LE 4 OR finishafter EQ 1.1 then begin
           Writefittingvariables,tirificsecond,inclinput1,painput1,vrotinput1,sbrinput1,sbrinput2,sbrinput4,sbrinput3,z0input1,$
                                 xposinput1,yposinput1,vsysinput1,sdisinput1
        ENDIF ELSE BEGIN
           Writefittingvariables,tirificsecond,inclinput1,inclinput2,inclinput3,painput1,painput2,painput3,vrotinput1,sbrinput1,sbrinput2,sbrinput4,sbrinput3,z0input1,$
                                 xposinput1,yposinput1,vsysinput1,sdisinput1
        ENDELSE
     endif else begin
                                ;Else we allow more variation in the center
        xposinput1=[ ' XPOS 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                     ' XPOS_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1), '360','0','1E-4','1E-5','1E-4','1E-5','3','70','70']
        yposinput1=[ ' YPOS 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                     ' YPOS_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1), '90','-90','1E-4','1E-5','1E-4','1E-5','3','70','70']
        vsysinput1=[ ' VSYS 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                     ' VSYS_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1),strtrim(strcompress(string(catvsys[i]+100.)),1),strtrim(strcompress(string(catvsys[i]-100.)),1),'2','0.5','2','0.5','3','70','70']
        IF norings[0] LE 4 or finishafter EQ 1.1 then begin
           Writefittingvariables,tirificsecond,xposinput1,yposinput1,vsysinput1,inclinput1,painput1,vrotinput1,$
                                 sbrinput1,sbrinput2,sbrinput4,sbrinput3,z0input1,sdisinput1
        ENDIF ELSE BEGIN
           Writefittingvariables,tirificsecond,xposinput1,yposinput1,vsysinput1,inclinput1,inclinput2,inclinput3,painput1,painput2,painput3,vrotinput1,$
                                 sbrinput1,sbrinput2,sbrinput4,sbrinput3,z0input1,sdisinput1
        ENDELSE
     endelse
                                ;Setting a bunch of fit tracking variables
     forcedring=0.
     letvrotvary=0.
     counter=0.
     trytwo=0.
     constring=0.
     sbrmodify=0.
     lastcutrings=0.
     noacceptbeforesmooth=0.
     lastaddrings=0.
     prevrings=norings[0]
     secondtime=0.
                                ;If the second fit is not accepted we
                                ;come back to this point to do the refitting 
     notacceptedtwo:
                                ;If we have optimized the cube we want
                                ;to set the _opt extension to the
                                ;input cubes name. Otherwise the input
                                ;cube should have properly copied from
                                ;the 1st fit  
     IF optimized then begin
        tmp=str_sep(strtrim(strcompress(currentfitcube),2),'_opt')
        if n_elements(tmp) LT 2 then begin
           currentfitcube=tmp[0]+'_opt'
           tmppos=where('INSET' EQ tirificsecondvars)
           tirificsecond[tmppos]='INSET=  '+strtrim(strcompress(string(currentfitcube+'.fits')))
        ENDIF
     ENDIF
                                ;If we are testing than skip the fitting
     IF testing GE 2 then goto,testing2
                                ;Write the tirific file and update log
     openw,1,maindir+'/'+catdirname[i]+'/tirific.def'
     for index=0,n_elements(tirificsecond)-1 do begin
        printf,1,tirificsecond[index]
     endfor
     close,1
     IF size(log,/TYPE) EQ 7 then begin
        openu,66,log,/APPEND
        printf,66,linenumber()+"Starting tirific Second estimate in "+catDirname[i]+" at "+systime()
        printf,66,linenumber()+"We have changed the ring numbers "+string(sbrmodify)+" times."
        printf,66,linenumber()+"We have changed the fitting paramaters "+string(trytwo)+" times."
        close,66
     ENDIF 
     print,linenumber()+"Starting tirific Second estimate in "+catDirname[i]+" at "+systime()
                                ;rename the files
     counter++
     
;     spawn,'mv 2ndfit.def 2ndfit_'+strtrim(strcompress(string(fix(counter))),2)+'.def'
     rename,'2ndfit.','2ndfitold.'
     gipsyfirst=strarr(1)
     gipsyfirst='tirific DEFFILE=tirific.def ACTION=1'
     spawn,gipsyfirst,isthere2
     testing2:
                                ;Then we check wether the fit is accepted if not we see if the PA and
                                ;INCL can be accepted IF the positions were not accepted before we try
                                ;then as well
     nopoints=0
     get_progress,maindir+'/'+catdirname[i]+'/progress2.txt',AC2,nopoints,loops,toymodels
     IF size(log,/TYPE) EQ 7 then begin
        openu,66,log,/APPEND
        IF AC2 EQ 1 then begin
           printf,66,linenumber()+"The Second estimate was accepted in this run." 
           printf,66,linenumber()+"Tirific ran through "+strtrim(strcompress(string(loops)))+" loops and produced "+strtrim(strcompress(string(toymodels)))+" models"
           printf,66,linenumber()+"Disk 1 had "+strtrim(strcompress(string(nopoints[0])))+" point sources and Disk 2 "+strtrim(strcompress(string(nopoints[1])))
        ENDIF ELSE BEGIN
           printf,66,linenumber()+"The Second estimate was not accepted."
           printf,66,linenumber()+"Tirific ran through "+strtrim(strcompress(string(loops)))+" loops and produced "+strtrim(strcompress(string(toymodels)))+" models"
           printf,66,linenumber()+"Disk 1 had "+strtrim(strcompress(string(nopoints[0])))+" point sources and Disk 2 "+strtrim(strcompress(string(nopoints[1])))
        ENDELSE
        close,66
     ENDIF   
                                ;let's moderate the cflux to ensure between 1.5 and 3 million points
     check_cflux,nopoints,tirificsecond,tirificsecondvars,cfluxadjusted,log=log
     IF AC2 EQ 1 then print,linenumber()+"The second estimate was accepted in this run." $
     else print,linenumber()+"The second estimate was not accepted."
     secondfitvaluesnames=0.
     secondfitvalues=0.
     finalsmooth=0.
                                ; If the fit is fully accepted and we
                                ; only want to do a smoothing we come
                                ; back to this point
     lastadjust:
 
                                ;We get the previous fit and read it into the fitting template
     writenewtotemplate,tirificsecond,maindir+'/'+catdirname[i]+'/2ndfit.def',Arrays=secondfitvalues,VariableChange=secondfitvaluesnames,Variables=tirificsecondvars
     tmppos=where('RADI' EQ secondfitvaluesnames)
     RADarr=secondfitvalues[*,tmppos]
     tmppos=where('VROT' EQ secondfitvaluesnames)
     VROTarr=secondfitvalues[*,tmppos]
     tmppos=where('SBR' EQ secondfitvaluesnames)
     SBRarr=secondfitvalues[*,tmppos]
     SBRarror=SBRarr
     tmppos=where('INCL' EQ secondfitvaluesnames)
     INCLang=secondfitvalues[*,tmppos]
     tmppos=where('PA' EQ secondfitvaluesnames)
     PAang=secondfitvalues[*,tmppos]
     tmppos=where('VROT_2' EQ secondfitvaluesnames)
     VROTarr2=secondfitvalues[*,tmppos]
     tmppos=where('SBR_2' EQ secondfitvaluesnames)
     SBRarr2=secondfitvalues[*,tmppos]
     SBRarr2or=SBRarr2
     tmppos=where('INCL_2' EQ secondfitvaluesnames)
     INCLang2=secondfitvalues[*,tmppos]
     tmppos=where('PA_2' EQ secondfitvaluesnames)
     PAang2=secondfitvalues[*,tmppos]    
     IF PAang2[0] NE PAang[0] then begin
        IF size(log,/TYPE) EQ 7 then begin
           openu,66,log,/APPEND
           printf,66,linenumber()+"This is where the inner PA starts to deviate."
           Close,66
        ENDIF        
     ENDIF

;Let's see how many of the inner rings we want to fix
     tmppos=where('VSYS' EQ firstfitvaluesnames)
     vsys=firstfitvalues[0,tmppos]
     levels=(sbrarr+sbrarr2)/2.*1000.
     Columndensity,levels,vsys,[catmajbeam[i],catminbeam[i]],/ARCSQUARE
     tmp=WHERE(levels GT 2E20)
     IF tmp[0] NE -1 then innerfix=floor(tmp[n_elements(tmp)-1]/1.5)-1.
     IF innerfix LT 4 OR innerfix GE norings[0] OR finishafter EQ 1.1 then innerfix=4
     tmppos=where('NUR' EQ secondfitvaluesnames)
     norings[0]=secondfitvalues[0,tmppos]   
     sbr_check,tirificsecond, tirificsecondvars,sbrarr,sbrarr2,cutoff     
                                ;We get the rings for which we only
                                ;want to fit a slope in the rotation curve
     get_newringsv8,SBRarr,SBRarr2,3.5*cutoff,velconstused
     IF norings[0] GT 8 AND not finishafter EQ 2.1 then velconstused=velconstused-1
     set_vrotv6,vrotinput1,VROTarr,velconstused,vrotmax,vrotmin,norings,channelwidth,avinner=avinner,finish_after=finishafter
                                ;Then set the surface brighness profile parameters
     set_sbr,SBRinput1,SBRinput2,SBRinput3,SBRinput4,SBRinput5,SBRinput6,SBRarr,cutoff,norings,finishafter
                                ;Update the rings to be fitted in the other parameters

     IF norings[0] LE 4 or finishafter EQ 1.1 then begin                 
        INCLinput1[0]='INCL 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                      ' INCL_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)
        PAinput1[0]='PA 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                    ' PA_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)         
     ENDIF ELSE BEGIN
                                ;Need to update the arrays otherwise might use
                                ;them later with wrong ring numbers
        INCLinput1[0]='INCL 1:'+strtrim(strcompress(string(innerfix,format='(F7.4)')),1)+' INCL_2 1:'+strtrim(strcompress(string(innerfix,format='(F7.4)')),1)
        PAinput1[0]='PA 1:'+strtrim(strcompress(string(innerfix,format='(F7.4)')),1)+' PA_2 1:'+strtrim(strcompress(string(innerfix,format='(F7.4)')),1)
                                ;update INCL and PA
        set_warp_slopev2,SBRarr,SBRarr2,cutoff,INCLinput2,PAinput2,INCLinput3,PAinput3,norings,log=log,innerfix=innerfix
     ENDELSE

                                ;additionallly we should updat the min
                                ;and max in the fitting if more than
                                ;two rings reach it
                                ;This doesn't happen when the
                                ;fit is accepted 
                                ;when this happens all rings after are
                                ;unreliable and should be set to the
                                ;last ring that isn't affected
     boundaryadjustment=0.
     lastreliablerings=n_elements(INCLang)
     IF double(INCLinput2[1]) LT 90. OR  double(INCLinput3[1]) LT 90.  OR double(INCLinput1[1]) LT 90. AND norings[0] GT 4 then begin
        tmp=WHERE(INCLang GE (double(INCLinput2[1])-double(INCLinput2[3])))
        tmp2=WHERE(INCLang2 GE (double(INCLinput3[1])-double(INCLinput3[3])))      
        IF n_elements(tmp) GE 2 OR n_elements(tmp2) GT 2 then begin
           boundaryadjustment=1
           IF double(INCLinput1[1]+5.) LT 90 then INCLinput1[1]=strtrim(strcompress(string(double(INCLinput1[1]+5.))),2) else INCLinput1[1]='90.'
           IF double(INCLinput2[1]+5.) LT 90. then INCLinput2[1]=strtrim(strcompress(string(double(INCLinput2[1]+5.))),2) else INCLinput2[1]='90.'
           IF double(INCLinput3[1]+5.) LT 90. then INCLinput3[1]=strtrim(strcompress(string(double(INCLinput3[1]+5.))),2) else INCLinput3[1]='90.'
           IF tmp[0] EQ -1 then tmp[0]=n_elements(INCLang)
           IF tmp2[0] EQ -1 then tmp2[0]=n_elements(INCLang)
           tmp=MAX([tmp,tmp2],min=lastreliablerings)
        ENDIF
     ENDIF
     IF double(INCLinput1[2]) GT 0. OR double(INCLinput2[2]) GT 0 OR double(INCLinput3[2]) GT 0 then begin
        IF norings[0] LE 4 then begin
           tmp=WHERE(INCLang LE (double(INCLinput1[2])+double(INCLinput1[3])))
           tmp2=WHERE(INCLang2 LE (double(INCLinput1[2])+double(INCLinput1[3])))
        ENDIF ELSE BEGIN
           tmp=WHERE(INCLang LE (double(INCLinput2[2])+double(INCLinput2[3])))
           tmp2=WHERE(INCLang2 LE (double(INCLinput2[2])+double(INCLinput2[3])))
        ENDELSE
        IF n_elements(tmp) GE 2 OR n_elements(tmp2) GT 2 then begin
           boundaryadjustment=1
           IF double(INCLinput1[2]-5.) GT 0. then INCLinput1[2]=strtrim(strcompress(string(double(INCLinput1[2]-5.))),2) else INCLinput1[2]='0.'
           IF double(INCLinput2[2]-5.) GT 0. then INCLinput2[2]=strtrim(strcompress(string(double(INCLinput2[2]-5.))),2) else INCLinput2[2]='0.'
           IF double(INCLinput3[2]-5.) GT 0. then INCLinput3[2]=strtrim(strcompress(string(double(INCLinput3[2]-5.))),2) else INCLinput3[2]='0.'
           IF tmp[0] EQ -1 then tmp[0]=n_elements(INCLang)
           IF tmp2[0] EQ -1 then tmp2[0]=n_elements(INCLang)
           tmp=MAX([tmp,tmp2,lastreliablerings],min=lastreliablerings)
        ENDIF
     ENDIF
     IF ABS(double(PAinput1[1])-double(PAinput1[2])) LT 400. OR  ABS(double(PAinput2[1])-double(PAinput2[2])) LT 400 OR  ABS(double(PAinput3[1])-double(PAinput3[2])) LT 400 AND norings[0] GT 4 then begin
        tmp=WHERE(PAang GE (double(PAinput2[1])-double(PAinput2[3])))
        tmp2=WHERE(PAang2 GE (double(PAinput2[1])-double(PAinput2[3])))
        IF n_elements(tmp) GE 2 OR n_elements(tmp2) GT 2 then begin
           boundaryadjustment=1
           PAinput1[1]=strtrim(strcompress(string(double(PAinput1[1]+10.))),2)
           PAinput2[1]=strtrim(strcompress(string(double(PAinput2[1]+10.))),2)
           PAinput3[1]=strtrim(strcompress(string(double(PAinput3[1]+10.))),2)
        ENDIF
        IF norings[0] LE 4 or finishafter EQ 1.1 then begin
           tmp=WHERE(PAang LE (double(PAinput1[2])+double(PAinput1[3])))
           tmp2=WHERE(PAang2 LE (double(PAinput1[2])+double(PAinput1[3])))
        ENDIF ELSE BEGIN
           tmp=WHERE(PAang LE (double(PAinput2[2])+double(PAinput2[3])))
           tmp2=WHERE(PAang2 LE (double(PAinput2[2])+double(PAinput2[3])))
        ENDELSE
        IF n_elements(tmp) GE 2 OR n_elements(tmp2) GT 2 then begin
           boundaryadjustment=1
           PAinput1[2]=strtrim(strcompress(string(double(PAinput1[2]-10.))),2)
           PAinput2[2]=strtrim(strcompress(string(double(PAinput2[2]-10.))),2)
           PAinput3[2]=strtrim(strcompress(string(double(PAinput3[2]-10.))),2)
           IF tmp[0] EQ -1 then tmp[0]=n_elements(INCLang)
           IF tmp2[0] EQ -1 then tmp2[0]=n_elements(INCLang)
           tmp=MAX([tmp,tmp2,lastreliablerings],min=lastreliablerings)
        ENDIF
     ENDIF
                                ;If we determined that the boundaries should be updated than we need
                                ;to check SBR and reset the values, as
                                ;well as smooth the INCL, PAand VROT
     IF boundaryadjustment then begin
        for j=n_elements(SBRarr)-1,0,-1 do begin
           IF SBRarr[j] LE cutoff[j] then begin
              SBRarr[j]=cutoff[j]*2.
           ENDIF else break
        endfor
        IF lastreliablerings EQ 0 then lastreliablerings=1
        IF lastreliablerings LT n_elements(SBRarr)-1 then begin
           INCLang[lastreliablerings:n_elements(SBRarr)-1]=INCLang[lastreliablerings-1]
           PAang[lastreliablerings:n_elements(SBRarr)-1]=PAang[lastreliablerings-1]
           IF lastreliablerings EQ 1 then VROTarr[lastreliablerings:n_elements(SBRarr)-1]=MEAN(VROTarr[1:n_elements(VROTarr)-1]) else VROTarr[lastreliablerings:n_elements(SBRarr)-1]=VROTarr[lastreliablerings-1]
           parameterreguv87,PAang,SBRarror,RADarr,fixedrings=3,difference=1.0,cutoff=cutoff,arctan=1,order=polorder,max_par=PAinput2[1],min_par=PAinput2[2],accuracy=1/4.,error=sigmapa1
           parameterreguv87,INCLang,SBRarror,RADarr,fixedrings=3,difference=8.*exp(-catinc[i]^2.5/10^3.5)+2.,cutoff=cutoff,arctan=1,order=polorder,max_par=PAinput2[1],min_par=PAinput2[2],accuracy=1,error=sigmaincl1
           stringINCL='INCL= '+STRJOIN(string(INCLang),' ')
           tmppos=where('INCL' EQ tirificsecondvars)
           tirificsecond[tmppos]=stringINCL      
           stringINCL='PA= '+STRJOIN(string(PAang),' ')
           tmppos=where('PA' EQ tirificsecondvars)
           tirificsecond[tmppos]=stringINCL
        ENDIF
        for j=n_elements(SBRarr2)-1,0,-1 do begin
           IF SBRarr2[j] LE cutoff[j] then begin
              SBRarr2[j]=cutoff[j]*2.
           ENDIF else break
        endfor
       
        IF lastreliablerings LT n_elements(SBRarr)-1 then begin
           INCLang2[lastreliablerings:n_elements(SBRarr2)-1]=INCLang2[lastreliablerings-1]
           PAang2[lastreliablerings:n_elements(SBRarr2)-1]=PAang2[lastreliablerings-1]
           IF lastreliablerings EQ 1 then VROTarr2[lastreliablerings:n_elements(SBRarr)-1]=MEAN(VROTarr2[1:n_elements(VROTarr2)-1]) else VROTarr2[lastreliablerings:n_elements(SBRarr2)-1]=VROTarr2[lastreliablerings-1]
           parameterreguv87,PAang2,SBRarr2or,RADarr,fixedrings=3,difference=1.0,cutoff=cutoff,arctan=1,order=polorder,max_par=PAinput2[1],min_par=PAinput2[2],accuracy=1/4.,error=sigmapa2
           parameterreguv87,INCLang2,SBRarr2or,RADarr,fixedrings=3,difference=8.*exp(-catinc[i]^2.5/10^3.5)+2.,cutoff=cutoff,arctan=1,order=polorder,max_par=PAinput2[1],min_par=PAinput2[2],accuracy=1,error=sigmaincl2
           stringINCL='INCL_2= '+STRJOIN(string(INCLang2),' ')
           tmppos=where('INCL_2' EQ tirificsecondvars)
           tirificsecond[tmppos]=stringINCL
           stringINCL='PA_2= '+STRJOIN(string(PAang2),' ')
           tmppos=where('PA_2' EQ tirificsecondvars)
           tirificsecond[tmppos]=stringINCL
        ENDIF
        VROTarr=(VROTarr+VROTarr2)/2.
        tmpSBR=(SBRarror+SBRarr2or)/2.
        verror=MAX([5.,channelwidth*(1.+vresolution)/SQRT(sin(catinc[i]*!DtoR))])
        IF double(VROTarr[1]) LT 120. AND  double(VROTarr[2]) LT 120. then begin
           parameterreguv87,VROTarr,tmpSBR, RADarr,/REVERSE,fixedrings=velfixrings,difference=verror,cutoff=cutoff,arctan=2,order=polorder,error=sigmarot
        ENDIF ELSE BEGIN
           parameterreguv87,VROTarr,tmpSBR, RADarr,/REVERSE,fixedrings=velfixrings,difference=verror,cutoff=cutoff,arctan=2,/NOCENTRAL,error=sigmarot
        ENDELSE

        stringVROT='VROT= 0. '+STRJOIN(string(VROTarr[1:n_elements(VROTarr)-1]),' ') 
        tmppos=where('VROT' EQ tirificsecondvars)
        tirificsecond[tmppos]=stringVROT
        stringVROT='VROT_2= 0. '+STRJOIN(string(VROTarr[1:n_elements(VROTarr)-1]),' ') 
        tmppos=where('VROT_2' EQ tirificsecondvars)
        tirificsecond[tmppos]=stringVROT
                                ;If we do this we also need to check the vrot min
        if VROTarr[lastreliablerings-1] LT avinner then vrotinput1[2]=string(VROTmin)+' '+string(VROTarr[lastreliablerings-1]*0.6)
        stringSBR='SBR= '+STRJOIN(string(SBRarr),' ') 
        tmppos=where('SBR' EQ tirificsecondvars)
        tirificsecond[tmppos]=stringSBR
        stringSBR='SBR_2= '+STRJOIN(string(SBRarr2),' ') 
        tmppos=where('SBR_2' EQ tirificsecondvars)
        tirificsecond[tmppos]=stringSBR
                                ;Write the fitting variables to files
        IF norings[0] LE 4 or finishafter EQ 1.1 then begin
           Writefittingvariables,tirificsecond,inclinput1,painput1,vrotinput1,$
                                 sbrinput1,sbrinput2,sbrinput4,sbrinput3,z0input1,sdisinput1
        ENDIF ELSE BEGIN
           Writefittingvariables,tirificsecond,inclinput1,inclinput2,inclinput3,painput1,painput2,painput3,vrotinput1,$
                                 sbrinput1,sbrinput2,sbrinput4,sbrinput3,z0input1,sdisinput1
        ENDELSE
        IF size(log,/TYPE) EQ 7 then begin
           openu,66,log,/APPEND
           printf,66,linenumber()+"We changed the boundaries. These are the new boundaries: "
           printf,66,linenumber()+"PA upper limits="+PAinput1[1]+" "+PAinput2[1]+" "+PAinput2[1]+" "
           printf,66,linenumber()+"PA lower limits="+PAinput1[2]+" "+PAinput2[2]+" "+PAinput2[2]+" "
           printf,66,linenumber()+"INCL upper limits="+INCLinput1[1]+" "+INCLinput2[1]+" "+INCLinput2[1]+" "
           printf,66,linenumber()+"INCL lower limits="+INCLinput1[2]+" "+INCLinput2[2]+" "+INCLinput2[2]+" "
           
           Close,66
        ENDIF
                                ;goto the point where the the second
                                ;fit is not accepted and try again
        goto,notacceptedtwo
     ENDIF                 
                                ;If we are smoothing with polynomials
                                ;we want to make sure that the minimum
                                ;sbr is actually of some value for the
                                ;weighing in the smoothing.
     IF finalsmooth EQ 1 then begin
        for j=n_elements(SBRarr)-1,0,-1 do begin
           IF SBRarr[j] LE cutoff[j] then begin
              SBRarr[j]=cutoff[j]*0.5
           ENDIF else break
        endfor      
        tmp=WHERE(SBRarr LT 0.)
        IF tmp[0] NE -1 then SBRarr[tmp]=cutoff[tmp]*0.95
        for j=n_elements(SBRarr2)-1,0,-1 do begin
           IF SBRarr2[j] LE cutoff[j] then begin
              SBRarr2[j]=cutoff[j]*0.5
           ENDIF else break
        endfor
        tmp=WHERE(SBRarr2 LT 0.)
        IF tmp[0] NE -1 then SBRarr2[tmp]=cutoff[tmp]*0.95
     ENDIF    
                                ;We will correct the output by fitting a polynomial
                                ;first we set the amount of rings that
                                ;are clearly one value
     get_fixedringsv8,[[PAang],[PAang2],[INCLang],[INCLang2]],fixedrings
     IF fixedrings GE norings[0] then begin
        IF norings[0] GT 5 then fixedrings=norings[0]-2 else fixedrings=norings[0]-1
     ENDIF
                                ;if we are not smoothing with a
                                ;polynomial we still want to simply
                                ;smooth the profile by averaging such
                                ;that the jigsaw pattern
                                ;doesn't amplify every iteration
     IF finalsmooth EQ 1 then prefunc=0 else prefunc=1
     centralincl=(INCLang[0]+INCLang2[0])/2.
     IF norings[0] GT 15 then accuracy=0.5+COS(centralincl*!DtoR)*0.5*15./norings[0] else accuracy=1
     IF finalsmooth LE 1 then begin
        sigmapa1=0.
        sigmapa2=0.
        sigmaincl1=0.
        sigmaincl2=0.
        sigmarot=0.
     ENDIF
                                ;Smoothing PA_1
     if finalsmooth LE 1 AND norings[0] GT 4 and finishafter NE 1.1 then begin     
        parameterreguv87,PAang,SBRarror,RADarr,fixedrings=fixedrings,difference=3.5-(ATAN((n_elements(norings[0])-5.)*2.))*5.25/!pi,cutoff=cutoff,arctan=prefunc,order=polorder,max_par=PAinput2[1],min_par=PAinput2[2],accuracy=accuracy/4.,error=sigmapa1
        IF size(log,/TYPE) EQ 7 then begin
           openu,66,log,/APPEND
           printf,66,linenumber()+"We used a polynomial of the   "+string(polorder)+" order to smooth the PA profile."
           printf,66,linenumber()+"The profile is forced to remain between "+PAinput2[1]+" and "+PAinput2[2]
           printf,66,linenumber()+"The inner" +string(fixedrings)+" rings are fixed."
           Close,66
        ENDIF
     ENDIF 
     padiff=PAang[*]-PAang[0]
                                ;Smoothing PA_2
     if finalsmooth LE 1 AND norings[0] GT 4 and finishafter NE 1.1 then begin
        parameterreguv87,PAang2,SBRarr2or,RADarr,fixedrings=fixedrings,difference=3.5-(ATAN((n_elements(norings[0])-5.)*2.))*5.25/!pi,cutoff=cutoff,arctan=prefunc,order=polorder,max_par=PAinput3[1],min_par=PAinput3[2],accuracy=accuracy/4.,error=sigmapa2
        IF size(log,/TYPE) EQ 7 then begin
           openu,66,log,/APPEND
           printf,66,linenumber()+"We used a polynomial of the   "+string(polorder)+" order to smooth the PA_2 profile."
           printf,66,linenumber()+"The profile is forced to remain between "+PAinput3[1]+" and "+PAinput3[2]
           printf,66,linenumber()+"The inner" +string(fixedrings)+" rings  are fixed."
           Close,66
        ENDIF
     ENDIF
                                ;making sure that the rings that are
                                ;fixed are the average of both sides
                                ;and the same
     PAang[0:fixedrings-1]=(double(PAang[0])+double(PAang2[0]))/2.
     PAang2[0:fixedrings-1]=PAang[0]
                                ;And writing them to the file template
     stringPA='PA= '+STRJOIN(string(PAang),' ')
     tmppos=where('PA' EQ tirificsecondvars)
     tirificsecond[tmppos]=stringPA
     stringPA='PA_2= '+STRJOIN(string(PAang2),' ')
     tmppos=where('PA_2' EQ tirificsecondvars)
     tirificsecond[tmppos]=stringPA
                                ;inclination
                                ;Smoothing INCL_1
     polorder1=0.
     if finalsmooth LE 1 AND norings[0] GT 4 and finishafter NE 1.1 then begin
        parameterreguv87,INCLang,SBRarror,RADarr,fixedrings=fixedrings,difference=8.*exp(-catinc[i]^2.5/10^3.5)+2.,cutoff=cutoff,arctan=prefunc,order=polorder1,max_par=INCLinput2[1],min_par=INCLinput2[2],accuracy=accuracy,error=sigmaincl1
        IF size(log,/TYPE) EQ 7 then begin
           openu,66,log,/APPEND
           printf,66,linenumber()+"We used a polynomial of the   "+string(polorder)+" order to smooth the INCL profile."
           printf,66,linenumber()+"The profile is forced to remain between "+INCLinput2[1]+" and "+INCLinput2[2]
           printf,66,linenumber()+"The inner" +string(fixedrings)+" rings are fixed."
           Close,66
        ENDIF
     ENDIF 
                                ;Smoothing INCL_2
     polorder2=0.
     if finalsmooth LE 1 AND norings[0] GT 4 and finishafter NE 1.1 then begin
        parameterreguv87,INCLang2,SBRarr2or,RADarr,fixedrings=fixedrings,difference=8.*exp(-catinc[i]^2.5/10^3.5)+2.,cutoff=cutoff,arctan=prefunc,order=polorder2,max_par=INCLinput3[1],min_par=INCLinput3[2],accuracy=accuracy,error=sigmaincl2
        IF size(log,/TYPE) EQ 7 then begin
           openu,66,log,/APPEND
           printf,66,linenumber()+"We used a polynomial of the   "+string(polorder)+" order to smooth the INCL_2 profile."
           printf,66,linenumber()+"The profile is forced to remain between "+INCLinput3[1]+" and "+INCLinput3[2]
           printf,66,linenumber()+"The inner" +string(fixedrings)+" rings are fixed."
           Close,66
        ENDIF
     ENDIF
                                ;making sure that the rings that are
                                ;fixed are the average of both sides
                                ;and the same  
     INCLang[0:fixedrings-1]=(double(INCLang[0])+double(INCLang2[0]))/2.
     INCLang2[0:fixedrings-1]=(double(INCLang[0])+double(INCLang2[0]))/2.
     stringINCL='INCL= '+STRJOIN(string(INCLang),' ')
     tmppos=where('INCL' EQ tirificsecondvars)
     tirificsecond[tmppos]=stringINCL
     stringINCL='INCL_2= '+STRJOIN(string(INCLang2),' ')
     tmppos=where('INCL_2' EQ tirificsecondvars)
     tirificsecond[tmppos]=stringINCL
                                ;rotation curve 
                                ;THis is a bit more tricky since we
                                ;want the flat parts on the outside we
                                ;need to inverse them
                                ;First let's see if the size of the
                                ;galaxy is reasonable 
                                ;As we fix the rotation curve to be
                                ;the same for both sides we want to
                                ;use an average surface brightness
                                ;profile
                                ;Not dividing the profile by two seems
                                ;like a mistake but as this
                                ;combination gives good results
                                ;let's leave it like this.
     SBRav=(SBRarr+SBRarr2)
     get_newringsv8,SBRav,SBRav,3.5*cutoff,velconstused    
     IF double(norings[0]) GT 15. then begin
        IF velconstused GT norings[0]-ceil(norings[0]/10.) then velconstused=norings[0]-ceil(norings[0]/10.)
     ENDIF
     prefunc=0. 
     IF norings[0]-velconstused LT 2 then velconstused=norings[0]-1
     IF norings[0] GT 8 AND not finishafter EQ 2.1 then velconstused=velconstused-1
     if finalsmooth LE 1 then begin
                                ;it is important to use the original
                                ;non centrallly modified profiles here
                                ;for the error weighing. For PA and
                                ;INCL it doesn't matter as
                                ;those errors are set to 1
        IF finalsmooth EQ 1 AND finishafter NE 1.1 AND norings[0] GT 4 then prefunc=0 else prefunc=2
                                ; IF finalsmooth EQ 1 then we first
                                ; want to refit the rotation curve
                                ; with the smoothed inclination as it
                                ; often drives eachother 
        IF finalsmooth EQ 1 then begin
           tmppos=where('VARINDX' EQ tirificsecondvars)
           tmp=strsplit(tirificsecond[tmppos],': ',/extract)
           IF n_elements(tmp) GT 3 then begin
              IF isnumeric(tmp[3]) then velconstused=tmp[3]-1 else velconstused=norings[0]
           ENDIF else velconstused=norings[0]
        ENDIF
        IF finalsmooth EQ 1 AND velconstused LT norings[0] AND (norings[0] LT 15 OR (polorder1 LT 4. AND polorder2 LT 4)) then begin
           vrotslopinput=strarr(10)
           start=3.+fix((norings[0]-3.)/5.)
           set_vrotv6,vrotslopinput,VROTarr,velconstused,vrotmax,vrotmin,norings,channelwidth,avinner=avinner,start=start,centralexclude=centralexclude,finish_after=finishafter 
           INCLinputall=['INCL 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+' '+$
                    'INCL_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1),$
                    '90','0',string(0.5),string(0.1),string(0.5),string(0.1),'3','70','70'] 
           Writefittingvariables,tirificsecond,INCLinputall,vrotslopinput
           openw,1,maindir+'/'+catdirname[i]+'/tirific.def'
           for index=0,n_elements(tirificsecond)-1 do begin
              printf,1,tirificsecond[index]
           endfor
           close,1
                                ;Perform tirific check with only changes as a whole to the parameters
           IF size(log,/TYPE) EQ 7 then begin
              openu,66,log,/APPEND
              printf,66,linenumber()+"Starting slope adjustment in "+catDirname[i]+" at "+systime()
              close,66
           ENDIF 
           print,linenumber()+"Starting slope adjustment in "+catDirname[i]+" at "+systime()
           rename,'2ndfit.','2ndfittmp.'
           gipsyfirst=strarr(1)
           gipsyfirst='tirific DEFFILE=tirific.def ACTION=1'
           spawn,gipsyfirst,isthere2
           writenewtotemplate,tirificsecond,maindir+'/'+catdirname[i]+'/2ndfit.def',Arrays=secondfitvalues,VariableChange=secondfitvaluesnames,Variables=tirificsecondvars
           rename,'2ndfit.','2ndfitslop.'
           rename,'2ndfittmp.','2ndfit.'
           tmppos=where('VROT' EQ secondfitvaluesnames)
           VROTarr=secondfitvalues[*,tmppos]
        ENDIF

;this is where we end this experiment
        tmpSBR=(SBRarror+SBRarr2or)/2.
        tmphigh=str_sep(strtrim(strcompress(VROTinput1[1]),2),' ')
        tmplow=str_sep(strtrim(strcompress(VROTinput1[2]),2),' ')
        velfixrings=norings[0]-velconstused
        IF velfixrings LT 1 then velfixrings=1
        IF velfixrings EQ 1 AND norings[0] GT 5 then velfixrings=2
        vmaxdev=MAX([30.,7.5*channelwidth*(1.+vresolution)])
        verror=MAX([5.,channelwidth*(1.+vresolution)/SQRT(sin(INCLang[0]*!DtoR))])
        IF size(log,/TYPE) EQ 7 then begin
           openu,66,log,/APPEND
           printf,66,linenumber()+"This is VROT."
           printf,66,VROTarr
           printf,66,linenumber()+"This is SBR."
           printf,66,tmpSBR
           printf,66,"The number of rings= "+strtrim(string(fix(norings[0])),2)
           printf,66,"The maximum/minimum= "+strtrim(string(tmphigh[0]),2)+"/"+strtrim(string(tmplow[0]),2)
           printf,66,"The number of fixed rings= "+strtrim(string(fix(velfixrings)),2)
           printf,66,"The velocity error = "+strtrim(string(verror),2)
           Close,66
        ENDIF
                                ;If the central parameters are LT 120
                                ;we want to take the central ring into
                                ;account otherwise not
        IF double(VROTarr[1]) LT 120. AND  double(VROTarr[2]) LT 120. then begin
           parameterreguv87,VROTarr,tmpSBR, RADarr,/REVERSE,fixedrings=velfixrings,difference=verror,cutoff=cutoff,arctan=prefunc,order=polorder,max_par=tmphigh[0],min_par=tmplow[0],accuracy=accuracy,error=sigmarot
        ENDIF ELSE BEGIN
           parameterreguv87,VROTarr,tmpSBR, RADarr,/REVERSE,fixedrings=velfixrings,difference=verror,cutoff=cutoff,arctan=prefunc,/NOCENTRAL,order=polorder,max_par=tmphigh[0],min_par=tmplow[0],accuracy=accuracy,error=sigmarot
        ENDELSE
        tmp0check=WHERE(VROTarr LT 0)
        IF tmp0check[0] NE -1 then VROTarr[tmp0check]=tmplow[0]
        IF size(log,/TYPE) EQ 7 then begin
           openu,66,log,/APPEND
           printf,66,linenumber()+"This is VROT after smoothing."
           printf,66,VROTarr
           Close,66
        ENDIF
        VROTarr[0]=0.
        IF size(log,/TYPE) EQ 7 then begin
           openu,66,log,/APPEND
           printf,66,linenumber()+"We used a polynomial of the "+string(polorder)+" order to smooth the VROT profile."
           printf,66,linenumber()+"And fixed "+string(velfixrings)+" rings."
           Close,66
        ENDIF
     ENDIF
                                ;As the SBR used is the average and
                                ;the rotation curve is the same for
                                ;both side we only need to fit one
                                ;side and can copy the
                                ;other. Velconstused also determines
                                ;the amount of rings only fitted in a slope
     IF size(log,/TYPE) EQ 7 then begin
        openu,66,log,/APPEND
        printf,66,linenumber()+"We fit a slope to the rotation curve from ring "+string(velconstused)+" to ring"+strtrim(string(fix(norings[0])),2)
        for x=0,n_elements(SBRarr)-1 do begin
           printf,66,SBRarr[x],SBRarr2[x],3.5*cutoff[x]
        endfor
        Close,66
     ENDIF ELSE BEGIN
        print,linenumber()+"We fit a slope to the rotation curve from ring "+string(velconstused)+" to ring"+strtrim(string(fix(norings[0])),2)
     ENDELSE
     VROTarr[0]=0.
     VROTarr2=VROTarr  
     IF velconstused LT norings[0]-1 AND finalsmooth LT 1 then begin
        VROTarr[velconstused-1:norings[0]-1]=VROTarr[velconstused-2]
        VROTarr2=VROTarr
     ENDIF
                                ;Write all to the template
     stringVROT='VROT= 0. '+STRJOIN(string(VROTarr[1:n_elements(VROTarr)-1]),' ') 
     tmppos=where('VROT' EQ tirificsecondvars)
     tirificsecond[tmppos]=stringVROT
     stringVROT='VROT_2= 0. '+STRJOIN(string(VROTarr[1:n_elements(VROTarr)-1]),' ') 
     tmppos=where('VROT_2' EQ tirificsecondvars)
     tirificsecond[tmppos]=stringVROT
     newindrings=[n_elements(SBRarr),n_elements(SBRarr2)]
                                ;IF we are on the last smoothing we
                                ;want to set the rings in the SBR
                                ;profile that are below the cutoff to
                                ;a very small number so they do not
                                ;affect the fit
     IF finalsmooth EQ 1 then begin
        get_newringsv8,SBRarr,SBRarr2,cutoff,newindrings,/individual
        IF newindrings[0] EQ n_elements(SBRarr)-2 and SBRarr[n_elements(SBRarr)-1] GT cutoff[n_elements(SBRarr)-1] then newindrings[0]=n_elements(SBRarr)
        IF newindrings[0] LT n_elements(SBRarr) then begin
           SBRarr[newindrings[0]-1:n_elements(SBRarr)-1]=1E-16
        ENDIF
        stringSBR='SBR= '+STRJOIN(string(SBRarr),' ') 
        tmppos=where('SBR' EQ tirificsecondvars)
        tirificsecond[tmppos]=stringSBR
        IF newindrings[1] EQ n_elements(SBRarr2)-2 and SBRarr2[n_elements(SBRarr2)-1] GT cutoff[n_elements(SBRarr2)-1] then newindrings[1]=n_elements(SBRarr2)
        IF newindrings[1] LT n_elements(SBRarr2) then begin
           SBRarr2[newindrings[1]-1:n_elements(SBRarr2)-1]=1E-16
        ENDIF
        stringSBR='SBR_2= '+STRJOIN(string(SBRarr2),' ') 
        tmppos=where('SBR_2' EQ tirificsecondvars)
        tirificsecond[tmppos]=stringSBR      
     ENDIF
     IF size(log,/TYPE) EQ 7 then begin
        openu,66,log,/APPEND
        printf,66,linenumber()+"The last SBR ring 1 is "+string(newindrings[0])
        printf,66,linenumber()+"The last SBR ring 2 is "+string(newindrings[1])
        printf,66,linenumber()+"The amount of rings is "+strtrim(string(fix(norings[0])),2)     
        close,66
     ENDIF
                                ;IF we have not done too many
                                ;iterations and we are not on the
                                ;final smoothing than we want to check
                                ;the extend of the models
     IF sbrmodify LT 10. AND finalsmooth NE 1 AND secondtime NE 1  then begin
                                ;check the SBR profiles to see if we
                                ;need to cut/add a ring.
        get_newringsv8,SBRarr,SBRarr2,cutoff,newrings 
        IF size(log,/TYPE) EQ 7 then begin
           openu,66,log,/APPEND
           printf,66,linenumber()+"Newrings should always be bigger or equal to velconstused "+strtrim(string(fix(newrings)),2)+" =>"+string(velconstused)
           close,66
        ENDIF 
                                ;Checking the new amount of rings against various criteria
                                ;Not too small
        if newrings LT 4 then begin
           IF size(log,/TYPE) EQ 7 then begin
              openu,66,log,/APPEND
              printf,66,linenumber()+"Cause the newrings was "+strtrim(string(fix(newrings)),2)+" we set it to 4."
              close,66
           ENDIF         
           newrings=4
        ENDIF
                                ;not too big
        IF newrings GT maxrings then begin
           IF size(log,/TYPE) EQ 7 then begin
              openu,66,log,/APPEND
              printf,66,linenumber()+"Cause the newrings was "+strtrim(string(fix(newrings)),2)+" we set it to maxrings."
              close,66
           ENDIF       
           newrings=maxrings 
        ENDIF
        IF size(log,/TYPE) EQ 7 then begin
           openu,66,log,/APPEND
           printf,66,linenumber()+"We found this as rings = "+strtrim(string(fix(norings[0])),2)+"  new="+strtrim(string(fix(newrings)),2)
           close,66
        ENDIF
                                ;Not previously fitted. If it has been
                                ;than fix to this number of rings.
        tmp=WHERE(prevrings EQ newrings)
        IF tmp[0] NE -1 AND newrings NE norings[0] then begin
           IF finalsmooth EQ 2 then begin 
              IF size(log,/TYPE) EQ 7 then begin
                 openu,66,log,/APPEND
                 printf,66,linenumber()+"Cause the newrings was "+strtrim(string(fix(newrings)),2)+" which we have fitted before and we are smoothing we are going to fix the rings at "+strtrim(string(fix(norings[0])),2)
                 close,66
              ENDIF  
              newrings=norings[0]
           ENDIF ELSE BEGIN
              IF size(log,/TYPE) EQ 7 then begin
                 openu,66,log,/APPEND
                 printf,66,linenumber()+"Cause the newrings was "+strtrim(string(fix(newrings)),2)+" which we have fitted before we gonna fix the rings at this value."
                 close,66
              ENDIF     
           ENDELSE
           secondtime=1
        ENDIF
                                ;IF we want to cut a ring we check
                                ;what we did previously and whether it
                                ;is acceptable.
        IF newrings LT norings[0] AND constring NE 1 AND newrings GT lowring then begin 
           IF newrings LE forcedring then begin
              IF norings[0] GT forcedring then begin
                 newrings=forcedring
                 constring=1
                 prevmodification=0.
              ENDIF ELSE BEGIN
                 IF size(log,/TYPE) EQ 7 then begin
                    openu,66,log,/APPEND
                    printf,66,linenumber()+"As the cut would modify the rings to less than the maximum forced addition of a ring we do not apply the cut."
                    printf,66,linenumber()+"new rings "+strtrim(string(fix(newrings)),2)+" old rings "+strtrim(string(fix(norings[0])),2)+" forcedring="+string(forcedring)
                    Close,66
                 ENDIF
                 goto,nocut
              ENDELSE
           ENDIF    
           IF prevmodification EQ 1 then begin
              IF newrings GE norings[0]-1 then begin
                 prevmodification=-1
              ENDIF else BEGIN
                 prevmodification=-2
              ENDELSE
           ENDIF 
           IF finalsmooth EQ 2 then begin
              newrings=norings[0]-1
              IF size(log,/TYPE) EQ 7 then begin
                 openu,66,log,/APPEND
                 printf,66,linenumber()+"As we are smoothing we will only subtract one ring."
                 printf,66,linenumber()+"new rings "+strtrim(string(fix(newrings)),2)+" old rings "+strtrim(string(fix(norings[0])),2)
                 Close,66
              ENDIF
           ENDIF
           oldrings=norings[0]
           IF newrings EQ norings[0] then goto,nocut else norings[0]=newrings
           
           IF size(log,/TYPE) EQ 7 then begin
              openu,66,log,/APPEND
              printf,66,linenumber()+"We cut a ring! lastcut"+string(lastcutrings)+" lastadd "+string(lastaddrings)+" new "+strtrim(string(fix(newrings)),2)+" old "+string(oldrings)
              Close,66
           ENDIF ELSE begin            
              print,linenumber()+"We cut a ring! lastcut"+string(lastcutrings)+" lastadd "+string(lastaddrings)+" new "+strtrim(string(fix(newrings)),2)+" old "+string(oldrings)
           ENDELSE
                                ;change the template to have a ring less
           changeradii,tirificsecond,norings[0]  
           lastcutrings=norings[0] 
                                ;Set the fitting parameters       
           set_vrotv6,vrotinput1,VROTarr,velconstused,vrotmax,vrotmin,norings,channelwidth,avinner=avinner,centralexclude=centralexclude,finish_after=finishafter    
                                ;Adapt the other fitting parameters
                                ;first the SBRarr
           set_sbr,SBRinput1,SBRinput2,SBRinput3,SBRinput4,SBRinput5,SBRinput6,SBRarr,cutoff,norings,finishafter,log=log
           IF norings[0] LE 4 OR finishafter EQ 1.1 then begin
              INCLinput1[0]='INCL 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                            ' INCL_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)              
              PAinput1[0]='PA 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                          ' PA_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)
              z0input1[0]='Z0 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                          ' Z0_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)
              xposinput1[0]=' XPOS 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                            ' XPOS_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)
              yposinput1[0]=' YPOS 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                            ' YPOS_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)
              vsysinput1[0]=' VSYS 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                            ' VSYS_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)
              sdisinput1[0]=' SDIS 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                            ' SDIS_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)
              Writefittingvariables,tirificsecond,inclinput1,painput1,$
                                    vrotinput1,sbrinput1,sbrinput2,sbrinput4,sbrinput3,z0input1,$
                                    xposinput1,yposinput1,vsysinput1,sdisinput1              
           ENDIF ELSE BEGIN
                                ;Need to update the arrays otherwise might use
                                ;them later with wrong ring numbers
              INCLinput1[0]='INCL 1:'+strtrim(strcompress(string(innerfix,format='(F7.4)')),1)+' INCL_2 1:'+strtrim(strcompress(string(innerfix,format='(F7.4)')),1)
              PAinput1[0]='PA 1:'+strtrim(strcompress(string(innerfix,format='(F7.4)')),1)+' PA_2 1:'+strtrim(strcompress(string(innerfix,format='(F7.4)')),1)
                                ;update INCL and PA
              set_warp_slopev2,SBRarr,SBRarr2,cutoff,INCLinput2,PAinput2,INCLinput3,PAinput3,norings,log=log,innerfix=innerfix
                                ;And the central parameters
              z0input1[0]='Z0 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                          ' Z0_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)
              xposinput1[0]=' XPOS 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                            ' XPOS_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)
              yposinput1[0]=' YPOS 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                            ' YPOS_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)
              vsysinput1[0]=' VSYS 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                            ' VSYS_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)
              sdisinput1[0]=' SDIS 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                            ' SDIS_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)
              IF norings[0] LE 6 then begin
                 Writefittingvariables,tirificsecond,inclinput1,inclinput2,inclinput3,painput1,painput2,painput3,$
                                       vrotinput1,sbrinput1,sbrinput2,sbrinput4,sbrinput3,z0input1,$
                                       xposinput1,yposinput1,vsysinput1,sdisinput1
              ENDIF ELSE BEGIN
                 Writefittingvariables,tirificsecond,inclinput1,inclinput2,inclinput3,painput1,painput2,painput3,$
                                       vrotinput1,sbrinput5,sbrinput1,sbrinput6,sbrinput2,sbrinput4,sbrinput3,z0input1,$
                                       xposinput1,yposinput1,vsysinput1,sdisinput1
              ENDELSE              
           ENDELSE
                                ;Update counters and go to a new
                                ;iteration of the model.
           sbrmodify++
           finalsmooth=0.
           overwrite=0.
           prevrings=[prevrings,norings[0]]
           goto,notacceptedtwo
        ENDIF
        nocut:
        overwrite=0
                                ;Under certain conditions we want to
                                ;add a previously cut ring back on
        IF (SBRarr[newrings-2] GT 7.5*cutoff[newrings-2] OR SBRarr2[newrings-2] GT 7.5*cutoff[newrings-2]) AND newrings GT norings[0] AND prevmodification EQ -1 then begin
           prevmodification=0.
           IF size(log,/TYPE) EQ 7 then begin
              openu,66,log,/APPEND
              printf,66,linenumber()+"We force added a ring to the profile."
              printf,66,linenumber()+"Newrings "+strtrim(string(fix(newrings)),2)+", Rings "+strtrim(string(fix(norings[0])),2)
              Close,66
           ENDIF
           IF newrings GT forcedring then forcedring=newrings
           overwrite=1 
        ENDIF
                                ;If we do not have those conditions
                                ;and the previous iteration was a cut
                                ;we do not add a ring
        IF prevmodification EQ -1 AND overwrite NE 1 then begin
           IF size(log,/TYPE) EQ 7 then begin
              openu,66,log,/APPEND
              printf,66,linenumber()+"We wanted to add a ring to the model but the previous step was a subtraction."
              Close,66
           ENDIF
        ENDIF
                                ;If all conditions for adding a ring
                                ;are satified then update the tirific file
        IF newrings GT norings[0] AND (prevmodification NE -1 OR overwrite EQ 1) AND newrings LT highring then begin        
           norings[0]=newrings
           addring:
           changeradii,tirificsecond,norings[0]            
                                ;we need to estimate some values for the added rings
           prefunc=0.
           tmppos=where('RADI' EQ tirificsecondvars)
           tmp=str_sep(strtrim(strcompress(tirificsecond[tmppos]),2),'=')
           RADarr=double(str_sep(strtrim(strcompress(tmp[1]),2),' '))
           tmp=[PAang,PAang[n_elements(PAang)-1]]
           tmpsbr=[SBRarror,cutoff[n_elements(SBRarr)]]
                                ;PA_1
           parameterreguv87,tmp,tmpsbr,RADarr,fixedrings=fixedrings,difference=2.,cutoff=cutoff,arctan=prefunc,order=polorder,max_par=PAinput2[1],min_par=PAinput2[2],/EXTENDING
           IF size(log,/TYPE) EQ 7 then begin
              openu,66,log,/APPEND
              printf,66,linenumber()+"We used a polynomial of the   "+string(polorder)+" order to extend the PA profile."
              Close,66
           ENDIF
           tmp2=PAang
           PAang=dblarr(n_elements(RADarr))
           PAang[0:n_elements(PAang)-2]=tmp2
           PAang[n_elements(PAang)-1]=tmp[n_elements(tmp)-1]
                                ;PA_2
           tmp=[PAang2,PAang2[n_elements(PAang2)-1]]
           tmpsbr2=[SBRarr2or,cutoff[n_elements(SBRarr2)]]
           parameterreguv87,tmp,tmpsbr2,RADarr,fixedrings=fixedrings,difference=2.,cutoff=cutoff,arctan=prefunc,order=polorder,max_par=PAinput3[1],min_par=PAinput3[2],/EXTENDING
           IF size(log,/TYPE) EQ 7 then begin
              openu,66,log,/APPEND
              printf,66,linenumber()+"We used a polynomial of the   "+string(polorder)+" order to extend the PA_2 profile."
              Close,66
           ENDIF 
           tmp2=PAang2
           PAang2=dblarr(n_elements(RADarr))
           PAang2[0:n_elements(PAang2)-2]=tmp2
           PAang2[n_elements(PAang2)-1]=tmp[n_elements(tmp)-1]
           PAang[0:fixedrings]=(PAang[0]+PAang2[0])/2.
           PAang2[0:fixedrings]=(PAang[0]+PAang2[0])/2.
           stringPA='PA= '+STRJOIN(string(PAang),' ')
           tmppos=where('PA' EQ tirificsecondvars)
           tirificsecond[tmppos]=stringPA
           stringPA='PA_2= '+STRJOIN(string(PAang2),' ')
           tmppos=where('PA_2' EQ tirificsecondvars)
           tirificsecond[tmppos]=stringPA
                                ;inclination
           tmp=[INCLang,INCLang[n_elements(INCLang)-1]]
           parameterreguv87,tmp,tmpsbr,RADarr,fixedrings=fixedrings,difference=6.*exp(-catinc[i]^2.5/10^3.5)+4.,cutoff=cutoff,arctan=prefunc,order=polorder,max_par=INCLinput2[1],min_par=INCLinput2[2],/EXTENDING
           IF size(log,/TYPE) EQ 7 then begin
              openu,66,log,/APPEND
              printf,66,linenumber()+"We used a polynomial of the   "+string(polorder)+" order to extend the INCL profile."
              Close,66
           ENDIF 
           tmp2=INCLang
           INCLang=dblarr(n_elements(RADarr))
           INCLang[0:n_elements(INCLang)-2]=tmp2
           INCLang[n_elements(INCLang)-1]=tmp[n_elements(tmp)-1]
                                ;INCL_2
           tmp=[INCLang2,INCLang2[n_elements(INCLang2)-1]]
           parameterreguv87,tmp,tmpsbr2,RADarr,fixedrings=fixedrings,difference=6.*exp(-catinc[i]^2.5/10^3.5)+4.,cutoff=cutoff,arctan=prefunc,order=polorder,max_par=INCLinput3[1],min_par=INCLinput3[2],/EXTENDING
           IF size(log,/TYPE) EQ 7 then begin
              openu,66,log,/APPEND
              printf,66,linenumber()+"We used a polynomial of the   "+string(polorder)+" order to extend the INCL_2 profile."
              Close,66
           ENDIF
           tmp2=INCLang2
           INCLang2=dblarr(n_elements(RADarr))
           INCLang2[0:n_elements(INCLang2)-2]=tmp2
           INCLang2[n_elements(INCLang2)-1]=tmp[n_elements(tmp)-1]
           INCLang[0:fixedrings]=(INCLang[0]+INCLang2[0])/2.
           INCLang2[0:fixedrings]=(INCLang[0]+INCLang2[0])/2.
           tmp=WHERE(Inclang GT 90.)
           IF tmp[0] NE -1 then INCLang[tmp]=90.
           tmp=WHERE(Inclang LT 0.)
           IF tmp[0] NE -1 then INCLang[tmp]=0.
           tmp=WHERE(Inclang2 GT 90.)
           IF tmp[0] NE -1 then INCLang2[tmp]=90.
           tmp=WHERE(Inclang2 LT 0.)
           IF tmp[0] NE -1 then INCLang2[tmp]=0.
           stringINCL='INCL= '+STRJOIN(string(INCLang),' ')
           tmppos=where('INCL' EQ tirificsecondvars)
           tirificsecond[tmppos]=stringINCL
           stringINCL='INCL_2= '+STRJOIN(string(INCLang2),' ')
           tmppos=where('INCL_2' EQ tirificsecondvars)
           tirificsecond[tmppos]=stringINCL
                                ;VROT
           IF n_elements(VROTarr) GE 4 then begin
              vrotext=string(strtrim(strcompress(TOTAL([VROTarr[n_elements(VROTarr)-4:n_elements(VROTarr)-1],VROTarr2[n_elements(VROTarr2)-4:n_elements(VROTarr2)-1]])/8.),2))
           ENDIF ELSE BEGIN
              vrotext=string(strtrim(strcompress(TOTAL([VROTarr,VROTarr2])/(n_elements(VROTarr)+n_elements(VROTarr))),2))
              vrotext2=string(strtrim(strcompress(TOTAL(VROTarr2)/n_elements(VROTarr2)),2))
           ENDELSE
           stringVROT='VROT= 0. '+STRJOIN(string(VROTarr[1:n_elements(VROTarr)-1]),' ')+' '+vrotext 
           tmppos=where('VROT' EQ tirificsecondvars)
           tirificsecond[tmppos]=stringVROT
           stringVROT='VROT_2= 0. '+STRJOIN(string(VROTarr[1:n_elements(VROTarr)-1]),' ')+' '+vrotext
           tmppos=where('VROT_2' EQ tirificsecondvars)
           tirificsecond[tmppos]=stringVROT
                                ;Update fitting setting and Log file
           lastaddrings=norings[0]
           IF size(log,/TYPE) EQ 7 then begin
              openu,66,log,/APPEND
              printf,66,linenumber()+"We added a ring to the model."
              Close,66
           ENDIF ELSE BEGIN
              print,linenumber()+"We added a ring to the model."
           ENDELSE
           prevmodification=1.
           INCLinput1[0]='INCL 1:'+strtrim(strcompress(string(innerfix,format='(F7.4)')),1)+' INCL_2 1:'+strtrim(strcompress(string(innerfix,format='(F7.4)')),1)           
           PAinput1[0]='PA 1:'+strtrim(strcompress(string(innerfix,format='(F7.4)')),1)+' PA_2 1:'+strtrim(strcompress(string(innerfix,format='(F7.4)')),1)       
                                ;Need to update the arrays otherwise might use
                                ;them later with wrong ring numbers
           IF norings[0] LE 4 OR finishafter EQ 1.1 then begin
              INCLinput1[0]='INCL 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+' INCL_2 1:' +strtrim(strcompress(string(norings[0],format='(F7.4)')),1)        
              PAinput1[0]='PA 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+' PA_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)              
           ENDIF ELSE BEGIN
                                ;update INCL and PA
              set_warp_slopev2,SBRarr,SBRarr2,cutoff,INCLinput2,PAinput2,INCLinput3,PAinput3,norings,log=log,innerfix=innerfix           
           ENDELSE
                                ;Updating surface brightness settings
           set_sbr,SBRinput1,SBRinput2,SBRinput3,SBRinput4,SBRinput5,SBRinput6,SBRarr,cutoff,norings,finishafter,log=log
                                ;Central parameter setting and Z0    
           z0input1[0]='Z0 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                       ' Z0_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)
           xposinput1[0]=' XPOS 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                         ' XPOS_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)
           yposinput1[0]=' YPOS 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                         ' YPOS_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)
           vsysinput1[0]=' VSYS 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                         ' VSYS_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)
           sdisinput1[0]=' SDIS 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                         ' SDIS_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)
                                ;Rotation settings
           set_vrotv6,vrotinput1,VROTarr,velconstused,vrotmax,vrotmin,norings,channelwidth,avinner=avinner,centralexclude=centralexclude,finish_after=finishafter
                                ;Write the parameters to the tirific array
           IF norings[0] LE 4 or finishafter EQ 1.1 then begin
              Writefittingvariables,tirificsecond,inclinput1,painput1,vrotinput1,$
                                    sbrinput1,sbrinput2,sbrinput4,sbrinput3,z0input1,xposinput1,yposinput1,vsysinput1 ,sdisinput1 
           ENDIF ELSE BEGIN
              IF norings[0] LE 6 then begin
                 Writefittingvariables,tirificsecond,inclinput1,inclinput2,inclinput3,painput1,painput2,painput3,vrotinput1,$
                                       sbrinput1,sbrinput2,sbrinput4,sbrinput3,z0input1,xposinput1,yposinput1,vsysinput1 ,sdisinput1 
              ENDIF else begin
                 Writefittingvariables,tirificsecond,inclinput1,inclinput2,inclinput3,painput1,painput2,painput3,vrotinput1,$
                                       sbrinput5,sbrinput1,sbrinput6,sbrinput2,sbrinput4,sbrinput3,z0input1,xposinput1,yposinput1,vsysinput1 ,sdisinput1 
              ENDELSE
           ENDELSE
                                ;Update counters and go back to the fitting proces
           sbrmodify++
           finalsmooth=0.
           prevrings=[prevrings,norings[0]]
           goto,notacceptedtwo
        ENDIF
     ENDIF
     noadd:
                                ;If we have tried the second fit three
                                ;times without changing the rings or
                                ;boundaries than accept that the
                                ;second fit failed.
     IF AC2 EQ 1 and trytwo EQ 3. then begin
        AC2=0
        AC1=1
     endif
                                ;IF the second fit is not accepted and
                                ;we have not tried three times without
                                ;changing anything already then update
                                ;parameters and try again
     IF AC2 NE 1  And trytwo LT 3.  then begin                
        INCLinput1[0]='INCL 1:'+strtrim(strcompress(string(innerfix,format='(F7.4)')),1)+' INCL_2 1:'+strtrim(strcompress(string(innerfix,format='(F7.4)')),1)         
        PAinput1[0]='PA 1:'+strtrim(strcompress(string(innerfix,format='(F7.4)')),1)+' PA_2 1:'+strtrim(strcompress(string(innerfix,format='(F7.4)')),1)       
                                ;Need to update the arrays otherwise might use
                                ;them later with wrong ring numbers
        IF norings[0] LE 4 OR finishafter EQ 1.1  then begin
           INCLinput1[0]='INCL 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+' INCL_2 1:' +strtrim(strcompress(string(norings[0],format='(F7.4)')),1)
           PAinput1[0]='PA 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+' PA_2 1:'+strtrim(strcompress(string(norings[0]-1,format='(F7.4)')),1)  
        ENDIF ELSE BEGIN
                                ;update INCL and PA
           set_warp_slopev2,SBRarr,SBRarr2,cutoff,INCLinput2,PAinput2,INCLinput3,PAinput3,norings,log=log,innerfix=innerfix
        ENDELSE
                                ;set the SBR fitting parameters
        set_sbr,SBRinput1,SBRinput2,SBRinput3,SBRinput4,SBRinput5,SBRinput6,SBRarr,cutoff,norings,finishafter,log=log
                                ;and some other parameters
        z0input1[0]='Z0 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                    ' Z0_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)
        xposinput1[0]=' XPOS 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                      ' XPOS_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)
        yposinput1[0]=' YPOS 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                      ' YPOS_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)
        vsysinput1[0]=' VSYS 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                      ' VSYS_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)
        sdisinput1[0]=' SDIS 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                      ' SDIS_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)
        set_vrotv6,vrotinput1,VROTarr,velconstused,vrotmax,vrotmin,norings,channelwidth,avinner=avinner,centralexclude=centralexclude,finish_after=finishafter 
                                ;But we need to reset all the fitting parameters
        IF trytwo LT 2. then begin
           IF AC1 NE 1 then begin
              ;If the first estimate is not accepted either 
              IF norings[0] LE 4 OR finishafter EQ 1.1  then begin
                 Writefittingvariables,tirificsecond,inclinput1,painput1,z0input1,xposinput1,yposinput1,$
                                       vsysinput1,vrotinput1,sbrinput1,sbrinput2,sbrinput4,sbrinput3,sdisinput1  
              ENDIF ELSE BEGIN
                 IF norings[0] LE 6 then begin
                    Writefittingvariables,tirificsecond,inclinput1,inclinput2,inclinput3,painput1,painput2,painput3,z0input1,xposinput1,yposinput1,$
                                          vsysinput1,vrotinput1,sbrinput1,sbrinput2,sbrinput4,sbrinput3,sdisinput1  
                 ENDIF else begin
                    Writefittingvariables,tirificsecond,inclinput1,inclinput2,inclinput3,painput1,painput2,painput3,z0input1,xposinput1,yposinput1,$
                                          vsysinput1,vrotinput1,sbrinput5,sbrinput1,sbrinput6,sbrinput2,sbrinput4,sbrinput3,sdisinput1 
                 ENDELSE
              ENDELSE
              IF size(log,/TYPE) EQ 7 then begin
                 openu,66,log,/APPEND
                 printf,66,linenumber()+"Both the first and the second estimate weren't accepted, retrying."
                 close,66
              ENDIF             
           endif else begin
                                ;IF only the second model is not
                                ;accepted try without fitting the center
              IF norings[0] LE 4  OR finishafter EQ 1.1  then begin
                 Writefittingvariables,tirificsecond,inclinput1,painput1,vrotinput1,sbrinput1,sbrinput2,sbrinput4,sbrinput3 ,sdisinput1 
              ENDIF ELSE BEGIN
                 IF norings[0] LE 6 then begin
                    Writefittingvariables,tirificsecond,inclinput1,inclinput2,inclinput3,painput1,painput2,painput3,vrotinput1,sbrinput1,sbrinput2,sbrinput4,sbrinput3,sdisinput1
                 ENDIF ELSE begin
                    Writefittingvariables,tirificsecond,inclinput1,inclinput2,inclinput3,painput1,painput2,painput3,vrotinput1,sbrinput5,sbrinput1,sbrinput6,sbrinput2,sbrinput4,sbrinput3,sdisinput1
                 ENDELSE
              ENDELSE
              IF size(log,/TYPE) EQ 7 then begin
                 openu,66,log,/APPEND                
                 printf,66,linenumber()+"The second estimate wasn't accepted. Running with fixed center."
                 close,66
              ENDIF
           endelse
        endif else begin
                                ;if we are on the third try
           IF AC1 NE 1 then begin
                                ;if the center is not yet accepted try
                                ;just changing the center
              xposinput1[3]='1E3'
              yposinput1[3]='1E3'
              vsysinput1[3]='4' 
              Writefittingvariables,tirificsecond,xposinput1,yposinput1,vsysinput1  
              IF size(log,/TYPE) EQ 7 then begin
                 openu,66,log,/APPEND
                 printf,66,linenumber()+"Last try on the second fit."               
                 close,66
              ENDIF              
           endif else begin
                                ;Else wrap up a failed fit
              goto,noretrytwo
           endelse
        endelse
                                ; update counter and try again
        trytwo++
        goto,notacceptedtwo
     endif
     noretrytwo:
                                ;if we get here without polynomial
                                ;smoothing then set smoothing and run parameters
     IF finalsmooth LT 1 then begin
        finalsmooth=1
        IF trytwo EQ 3 and AC2 NE 1 then noacceptbeforesmooth=1
        IF size(log,/TYPE) EQ 7 then begin
           openu,66,log,/APPEND
           printf,66,linenumber()+"We have a succesful model and we will smooth with polynomials now." 
           close,66
        END
        goto,lastadjust
     ENDIF 
                                ;IF the fitting rings are very small set them to the minimum of 4
     IF newindrings[0] LT 4 then newindrings[0]=4
     IF newindrings[1] LT 4 then newindrings[1]=4
     IF size(log,/TYPE) EQ 7 then begin
        openu,66,log,/APPEND
        printf,66,linenumber()+"We'll fit from "+strtrim(strcompress(string(newindrings,format='(F7.4)')),1)
        close,66
     ENDIF
                                ;set the parameters for single fitting
     SBRinput1=['!SBR '+strtrim(strcompress(string(newindrings[0],format='(F7.4)')),1)+':3','1',strtrim(strcompress(string(cutoff[n_elements(cutoff)-1]/2.,format='(E12.5)')),1),'1E-6','1E-7','1E-5','1E-7','3','70','70']
     SBRinput2=SBRinput1
     SBRinput2[0]='!SBR_2 '+strtrim(strcompress(string(newindrings[1],format='(F7.4)')),1)+':3'
     SBRinput3=SBRinput1
     SBRinput3[0]='SBR 1 SBR_2 1'
     SBRinput3[1]=strtrim(strcompress(string(SBRarr[1],format='(E12.5)')),1)
     SBRinput3[2]=strtrim(strcompress(string(cutoff[n_elements(cutoff)-1],format='(E12.5)')),1)
     SBRinput4=SBRinput1
     SBRinput4[0]='SBR 2 SBR_2 2'
     SBRinput4[2]=strtrim(strcompress(string(cutoff[n_elements(cutoff)-1],format='(E12.5)')),1)
     IF norings GT 4 and finishafter NE 1.1 then begin 
                                ;PA
        PAinput1=['PA 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+' '+$
                  'PA_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1),$
                  '360','0',string(0.5),string(0.1),string(0.5),string(0.1),'3','70','70']
                                ;INCL
        INCLinput1=['INCL 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+' '+$
                    'INCL_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1),$
                    '90','0',string(0.5),string(0.1),string(0.5),string(0.1),'3','70','70']
 
        VROTinput1=['VROT 2:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+' VROT_2 2:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)$
                    ,'500',string(channelwidth),string(channelwidth*0.5),string(0.1*channelwidth),string(0.5*channelwidth),string(0.1*channelwidth),'3','70','70',' ']
                                ;SDIS
        SDISinput1=['SDIS 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+' SDIS_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1),$
                    '40','3','1','0.1','0.5','0.05','3','70','70']    
                                ;Z0
        Z0input1=['Z0 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                  ' Z0_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1),$
                  string(convertskyanglefunction(0.4,double(catDistance[i]),/PHYSICAL)),$
                  string(convertskyanglefunction(0.075,double(catDistance[i]),/PHYSICAL)),$
                  string(convertskyanglefunction(0.005,double(catDistance[i]),/PHYSICAL)),$
                  string(convertskyanglefunction(0.0005,double(catDistance[i]),/PHYSICAL)),$
                  string(convertskyanglefunction(0.5,double(catDistance[i]),/PHYSICAL)),$
                  string(convertskyanglefunction(0.0005,double(catDistance[i]),/PHYSICAL)),'3','70','70']

        IF catinc[i] GT 80 then begin
           Z0input1[1]=string(MAX([convertskyanglefunction(1.0,double(catDistance[i])),catmajbeam[i]/2.]))
        ENDIF
        IF catDistance[i] EQ 1. then begin
                                ;ok if we do not know the distance
                                ;then let us assume each disk is about
                                ;30 kpc which means that
                                ;0.5=norings[0]/60. and so on
           IF catinc[i] GT 80 then begin
              Z0input1[1]=string(MAX([catmajbeam[i]*norings[0]/60.,catmajbeam[i]/2.]))
           ENDIF ELSE  Z0input1[1]=string(catmajbeam[i]*norings[0]/60.)
           Z0input1[2]='0.'
           Z0input1[3]=string(-1*catmajbeam[i]*norings[0]/6E4)
           Z0input1[4]=string(catmajbeam[i]*norings[0]/6E5)
           Z0input1[5]=string(catmajbeam[i]*norings[0]/60.)
           Z0input1[6]=string(catmajbeam[i]*norings[0]/6E5)    
        ENDIF
     ENDIF ELSE BEGIN
        VROTinput1=['VROT 2:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+' VROT_2 2:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)$
                    ,'500',string(channelwidth),string(channelwidth*0.5),string(0.1*channelwidth),string(0.5*channelwidth),string(0.1*channelwidth),'3','70','70',' ']
     ENDELSE
                                ;write the parameters to the tirific array
     Writefittingvariables,tirificsecond,PAinput1,INCLinput1,VROTinput1,SBRinput1,SBRinput2,SBRinput4,SDISinput1,Z0input1 
     IF testing GE 2 OR finalsmooth EQ 2. then goto,testing2check
     tmppos=where('PA' EQ tirificsecondvars)
     stringPA=tirificsecond[tmppos]
     tmppos=where('PA_2' EQ tirificsecondvars)
     stringPA=tirificsecond[tmppos]
                                ;Write to file
     openw,1,maindir+'/'+catdirname[i]+'/tirific.def'
     for index=0,n_elements(tirificsecond)-1 do begin
        printf,1,tirificsecond[index]
     endfor
     close,1
                                ;Perform tirific check with only changes as a whole to the parameters
     IF size(log,/TYPE) EQ 7 then begin
        openu,66,log,/APPEND
        printf,66,linenumber()+"Starting tirific check of the smoothed second estimate in "+catDirname[i]+" at "+systime()
        printf,66,linenumber()+"We have changed the ring numbers "+string(sbrmodify)+" times."
        printf,66,linenumber()+"We have changed the fitting paramaters "+string(trytwo)+" times."
        close,66
     ENDIF 
     print,linenumber()+"Starting tirific check of second estimate in "+catDirname[i]+" at "+systime()
     rename,'2ndfit.','2ndfitunsmooth.'
     gipsyfirst=strarr(1)
     gipsyfirst='tirific DEFFILE=tirific.def ACTION=1'
     spawn,gipsyfirst,isthere2
     finalsmooth=2.
     overwrite=0.
     get_progress,maindir+'/'+catdirname[i]+'/progress2.txt',AC2,nopoints,loops,toymodels
     IF size(log,/TYPE) EQ 7 then begin
        openu,66,log,/APPEND
        IF AC2 EQ 1 then begin
           printf,66,linenumber()+"The smoothed second estimate was accepted in this run." 
           printf,66,linenumber()+"Tirific ran through "+strtrim(strcompress(string(loops)))+" loops and produced "+strtrim(strcompress(string(toymodels)))+" models."
           printf,66,linenumber()+"Disk 1 had "+strtrim(strcompress(string(nopoints[0])))+" point sources and Disk 2 "+strtrim(strcompress(string(nopoints[1])))
        ENDIF ELSE BEGIN
           printf,66,linenumber()+"The smoothed second estimate was not accepted."
           printf,66,linenumber()+"Tirific ran through "+strtrim(strcompress(string(loops)))+" loops and produced "+strtrim(strcompress(string(toymodels)))+" models."
           printf,66,linenumber()+"Disk 1 had "+strtrim(strcompress(string(nopoints[0])))+" point sources and Disk 2 "+strtrim(strcompress(string(nopoints[1])))
        ENDELSE
        close,66
     ENDIF  
                                ;Go and see wether all parameters are satisfied
     goto,lastadjust
    
     testing2check:


                
  

                           ;Clean up the model        
     IF finalsmooth LT 3 then begin
        firstfitvaluesnames=0.
        writenewtotemplate,tirificsecond,maindir+'/'+catdirname[i]+'/2ndfit.def',Arrays=firstfitvalues,VariableChange=firstfitvaluesnames,Variables=tirificsecondvars
        tmp=WHERE(firstfitvaluesnames EQ 'SBR')
        SBRarr=firstfitvalues[*,tmp]
        tmp=WHERE(firstfitvaluesnames EQ 'SBR_2')
        SBRarr2=firstfitvalues[*,tmp]
        tmpSBR=(SBRarr+SBRarr2)/2.
        IF finishafter EQ 1.1 then begin
           get_newringsv8,tmpSBR,tmpSBR,cutoff,newend
           IF newend LT 4 then newend=4
           rename,'2ndfit.','2ndfituncor.'
           tmp=WHERE(firstfitvaluesnames EQ 'NUR')
           rings=firstfitvalues[0,tmp]
           IF newend LT rings then begin
              changeradii,tirificsecond,newend
              INCLinput1=['INCL 1','90','89',string(0.5),string(0.1),string(0.5),string(0.1),'3','70','70',' ']
              writefittingvariables,tirificsecond,INCLinput1
              newrings=newend
           ENDIF ELSE newrings=rings
        ENDIF
        IF optimized then begin
           tmp=str_sep(strtrim(strcompress(currentfitcube),2),'_opt')
           IF n_elements(tmp) EQ 2 then begin
              currentfitcube=tmp[0]
              tmppos=where('INSET' EQ tirificsecondvars)
              tirificsecond[tmppos]='INSET=  '+strtrim(strcompress(string(currentfitcube+'.fits')))
           ENDIF
        ENDIF
        tmp1=WHERE(SBRarr LT 1.1E-16)
        tmp2=WHERE(SBRarr2 LT 1.1E-16)
        IF n_elements(tmp1) GT 1 AND n_elements(tmp2) GT 1 then begin       
           IF tmp1[0] GT tmp2[0] then newrings=tmp1[0]+1 else newrings=tmp2[0]+1
           IF newrings LT 4 then newrings=4
           changeradii,tirificsecond,newrings
        ENDIF ELSE BEGIN
           IF finishafter NE 1.1 then newrings=n_elements(SBRarr)
        ENDELSE
        tmppos=where('LOOPS' EQ tirificsecondvars)
        tirificsecond[tmppos]='LOOPS=  0'
        errorsadd=['VROT','VROT_2','PA','PA_2','INCL','INCL_2']
        tmpfile=strarr(n_elements(tirificsecond)+6)
        added=0
                                ;Add the errors to the final tirific
                                ;file
        for j=0,n_elements(tirificsecond)-1 do begin
           tmp=str_sep(strtrim(strcompress(tirificsecond[j]),2),'=')
           tmpex=WHERE(tmp[0] EQ errorsadd)
           if tmpex[0] NE -1 then begin
              tmpfile[j+added]=tirificsecond[j]
              added++
              case tmpex[0] of
                 0:begin
                    IF n_elements(sigmarot) LT newrings then sigmarot=replicate(channelwidth,newrings) 
                    tmpinc=WHERE(firstfitvaluesnames EQ 'INCL')
                    incarr1=double(firstfitvalues[*,tmpinc])
                    errs=dblarr(n_elements(sigmarot))
                    errs=sigmarot/SIN(incarr1*!DtoR)
                    errs[0]=double(channelwidth)
                    tmpposv=WHERE(firstfitvaluesnames EQ 'VROT')
                    avrot=TOTAL(firstfitvalues[1:n_elements(firstfitvalues[*,tmpposv])-1,tmpposv])/n_elements(firstfitvalues[1:n_elements(firstfitvalues[*,tmpposv])-1,tmpposv])  
                    tmpavrot=WHERE(errs GT avrot/2.)
                    
                    IF tmpavrot[0] NE -1 AND avrot NE 0 then errs[WHERE(errs GT avrot/2.)]=avrot/2.
                    tmpfile[j+added]='# '+tmp[0]+'_ERR='+STRJOIN(strtrim(strcompress(string(errs[0:newrings-1]))),' ')
                 end
                 1: begin
                    IF n_elements(sigmarot) LT newrings then sigmarot=replicate(channelwidth,newrings)  
                    tmpinc=WHERE(firstfitvaluesnames EQ 'INCL_2')
                    incarr2=firstfitvalues[*,tmpinc]
                    errs=sigmarot/SIN(incarr2*!DtoR)
                    errs[0]=double(channelwidth)
                    tmpposv=WHERE(firstfitvaluesnames EQ 'VROT_2')
                    avrot=TOTAL(firstfitvalues[1:n_elements(firstfitvalues[*,tmpposv])-1,tmpposv])/n_elements(firstfitvalues[1:n_elements(firstfitvalues[*,tmpposv])-1,tmpposv])
                    tmpavrot=WHERE(errs GT avrot/2.)
                    
                    IF tmpavrot[0] NE -1 AND avrot NE 0 then errs[WHERE(errs GT avrot/2.)]=avrot/2.
                    tmpfile[j+added]='# '+tmp[0]+'_ERR='+STRJOIN(strtrim(strcompress(string(errs[0:newrings-1]))),' ')
                 end
                 2: begin
                    IF n_elements(sigmapa1) LT newrings then sigmapa1=replicate(2,newrings) 
                    tmpfile[j+added]='# '+tmp[0]+'_ERR='+STRJOIN(strtrim(strcompress(string(sigmapa1[0:newrings-1]))),' ')
                 end
                 3: begin
                     IF n_elements(sigmapa2) LT newrings then  sigmapa2=replicate(2,newrings)
                    tmpfile[j+added]='# '+tmp[0]+'_ERR='+STRJOIN(strtrim(strcompress(string(sigmapa2[0:newrings-1]))),' ')
                 end
                 4: begin
                     IF n_elements(sigmaincl1) LT newrings then  sigmaincl1=replicate(4,newrings)
                     tmpfile[j+added]='# '+tmp[0]+'_ERR='+STRJOIN(strtrim(strcompress(string(sigmaincl1[0:newrings-1]))),' ')
                  end
                 5: begin
                     IF n_elements(sigmaincl2) LT newrings then   sigmaincl2=replicate(4,newrings)
                    tmpfile[j+added]='# '+tmp[0]+'_ERR='+STRJOIN(strtrim(strcompress(string(sigmaincl2[0:newrings-1]))),' ')
                 end
                 else:print,'odd'
              endcase
           endif else tmpfile[j+added]=tirificsecond[j]
        endfor
        tirificsecond=tmpfile
                                ;write the final file and produce the final model
        openw,1,maindir+'/'+catdirname[i]+'/tirific.def'
        for index=0,n_elements(tirificsecond)-1 do begin
           printf,1,tirificsecond[index]
        endfor
        close,1
        if optimized then begin
           rename,'2ndfit.','2ndfit_opt.'
        endif
        gipsyfirst=strarr(1)
        gipsyfirst='tirific DEFFILE=tirific.def ACTION=1'
        spawn,gipsyfirst,isthere2
        spawn,'cp tirific.def 2ndfit.def',isthere
        finalsmooth=3
        secondtime=1
        openu,66,log,/APPEND
        printf,66,linenumber()+"This should be the very last thing and not happening again." 
        close,66
     ENDIF
                                ;Write the final info and pv-diagrams
     Basicinfovars=['XPOS','YPOS','VSYS','PA','INCL','VROT']
     writenewtotemplate,tirificsecond,maindir+'/'+catdirname[i]+'/2ndfit.def',Arrays=Basicinfovalues,VariableChange=Basicinfovars,Variables=tirificsecondvars,/EXTRACT
     RAhr=Basicinfovalues[0,0]
     RAdiff=maxchangeRA*3600./15.
     DEChr=Basicinfovalues[0,1]
     DECdiff=maxchangeDEC*3600.
     tmpcube=readfits(maindir+'/'+catdirname[i]+'/2ndfit.fits',hedtmp1stcube)
     IF optimized then begin
        extract_pv,nooptcube,nopthed,Basicinfovalues[0,3],xv,center=[RAhr,DEChr],xvheader=new_header
        writefits,maindir+'/'+catdirname[i]+'/'+noptname[0]+'_2_xv.fits',xv,new_header
     ENDIF ELSE BEGIN
        extract_pv,dummy,hed,Basicinfovalues[0,3],xv,center=[RAhr,DEChr],xvheader=new_header
        writefits,maindir+'/'+catdirname[i]+'/'+catcubename[i]+'_2_xv.fits',xv,new_header
     ENDELSE
     extract_pv,tmpcube,hedtmp1stcube,Basicinfovalues[0,3],xv,center=[RAhr,DEChr],xvheader=new_header
     writefits,maindir+'/'+catdirname[i]+'/2ndfit_xv.fits',xv,new_header
     hedtmp1st=hedtmp1stcube
     tmpix=WHERE(tmpcube GT catnoise[i])
     totflux=[TOTAL(tmpcube[tmpix])/pixperbeam,(2.*cutoff*norings[0])/(n_elements(tmpix)/pixperbeam)]
     mask=readfits(maindir+'/'+catdirname[i]+'/'+catmaskname[i]+'.fits',hedmask,/NOSCALE)        
                                ;mask the data cube
     tmpmask=fltarr(n_elements(tmpcube[*,0,0]),n_elements(tmpcube[0,*,0]),n_elements(tmpcube[0,0,*]))
     tmpmask[WHERE(mask GT 0.)]=tmpcube[WHERE(mask GT 0.)]
     momentsv2,tmpmask,tmpmap,hedtmp1st,0.
     writefits,maindir+'/'+catdirname[i]+'/2ndfit_mom0.fits',tmpmap,hedtmp1st
     hedtmp1stv=hedtmp1stcube
     momentsv2,tmpmask,tmpmapv,hedtmp1stv,1.
     writefits,maindir+'/'+catdirname[i]+'/2ndfit_mom1.fits',tmpmapv,hedtmp1stv
     getDHI,tmpmap,hedtmp1st,Basicinfovalues[0,3],[RAhr,DEChr,Basicinfovalues[0,4]],DHI
     VSYSdiff=maxchangevel
     HIMASS=2.36E5*catDistance[i]^2*totflux*ABS(channelwidth)
     convertradec,RAhr,DEChr
     openu,1,basicinfofile,/APPEND
     printf,1,"#After the second fit"
     printf,1,format='(A25,A25,A20,A20,A20,A20,A20,A20,A20,A20,A20,A20,A20)',string(RAhr+'+/-'+strtrim(strcompress(string(RAdiff,format='(F6.1)')),2)),$
            string(DEChr+'+/-'+strtrim(strcompress(string(DECdiff,format='(F6.1)')),2)),$
            string(strtrim(strcompress(string(Basicinfovalues[0,2])),2)+'+/-'+strtrim(strcompress(string(VSYSdiff,format='(F6.1)')),2)),$
            string(strtrim(strcompress(string(Basicinfovalues[0,3])),2)+'+/-'+strtrim(strcompress(string(catPAdev[i],format='(F6.1)')),2)),$
            string(strtrim(strcompress(string(Basicinfovalues[0,4])),2)+'+/-'+strtrim(strcompress(string(catincdev[i],format='(F6.1)')),2)),$
            string(strtrim(strcompress(string(Basicinfovalues[n_elements(Basicinfovalues[*,5])-1,5])),2)+'+/-'+strtrim(strcompress(string(catmaxrotdev[i],format='(F6.1)')),2)),$
            string(strtrim(strcompress(string(W50)),2)+'+/-'+strtrim(strcompress(string(channelwidth,format='(F6.1)')),2)),$
            string(strtrim(strcompress(string(totflux[0],format='(E10.3)')),2)),$
            string(strtrim(strcompress(string(DHI,format='(F8.1)')),2)),$
            string(strtrim(strcompress(string(catDistance[i])),2)),$
            string(strtrim(strcompress(string(HIMass[0],format='(E10.3)')),2)),$
            string(strtrim(strcompress(string(convertskyanglefunction(DHI,catDistance[i]),format='(F8.1)')),2))
     close,1
     
     IF size(log,/TYPE) EQ 7 then begin
        openu,66,log,/APPEND
        printf,66,linenumber()+"We finished the second fit of "+catDirname[i]+"  at "+systime()
        case prevmodification of
           0:lastringstr=" no change"
           -1:lastringstr=" a subtraction"
           1:lastringstr=" adding a ring"
           else: lastringstr=" something weird happened"
        endcase        
        printf,66,linenumber()+"We have rerun "+string(counter)+" times and the last modification was "+lastringstr
        close,66
     ENDIF 
     IF noacceptbeforesmooth EQ 1 then begin
        IF AC2 EQ 1 then AC2=2
     ENDIF
     IF finishafter EQ 2 OR finishafter EQ 2.1 then begin
                                ;Below 8 beams in diameter FAT is not
                                ;reliable so if in that range or when
                                ;the inclination is below 20 the fit
                                ;is failed
        case 1 of
           norings[0] LT 5:begin
              openu,1,outputcatalogue,/APPEND
              printf,1,format='(A60,2A12,A120)',catDirname[i],AC1,0.,' The finalmodel is less than 8 beams in diameter. FAT is not reliable in this range.'
              close,1
           end
           Basicinfovalues[0,4] LT 20:begin
              openu,1,outputcatalogue,/APPEND
              printf,1,format='(A60,2A12,A120)',catDirname[i],AC1,0.,' The final inclination is below 20 degrees. FAT is not reliable in this range.'
              close,1
           end
           else:begin
              openu,1,outputcatalogue,/APPEND
              printf,1,format='(A60,2A12,A120)',catDirname[i],AC1,AC2,' You have chosen to skip the fitting process after the second fit'
              close,1  
           end
        endcase
        goto,finishthisgalaxy
     ENDIF ELSE begin
        IF finishafter EQ 1.1 then begin
                                ;Below 8 beams in diameter FAT is not
                                ;reliable so if in that range or when
                                ;the inclination is below 20 the fit
                                ;is failed
           case 1 of
              norings[0] LT 9:begin
                 openu,1,outputcatalogue,/APPEND
                 printf,1,format='(A60,2A12,A120)',catDirname[i],AC1,0.,'The finalmodel is less than 8 beams in diameter. FAT is not reliable in this range.'
                 close,1
              end
              Basicinfovalues[0,4] LT 20:begin
                 openu,1,outputcatalogue,/APPEND
                 printf,1,format='(A60,2A12,A120)',catDirname[i],AC1,0.,'The final inclination is below 20 degrees. FAT is not reliable in this range.'
                 close,1
              end
              else:begin
                 openu,1,outputcatalogue,/APPEND
                 printf,1,format='(A60,2A12,A120)',catDirname[i],AC1,AC2,'As the galaxy radius is already halved we stop here'
                 close,1
              end
           endcase
           goto,finishthisgalaxy
        ENDIF ELSE BEGIN
           IF size(log,/TYPE) EQ 7 then begin
              openu,66,log,/APPEND
              IF AC2 EQ 1 then printf,66,linenumber()+"The Second estimate was accepted." $
              else printf,66,linenumber()+"The Second estimate was not accepted."
              close,66
           ENDIF ELSE BEGIN
              IF AC2 EQ 1 then print,linenumber()+"The Second estimate was accepted." $
              else print,linenumber()+"The Second estimate was not accepted."
           ENDELSE
        ENDELSE
     ENDELSE
     IF AC1 EQ 0 AND AC2 EQ 0 then begin
        openu,1,outputcatalogue,/APPEND
        printf,1,format='(A60,2A12,A80)',catDirname[i],AC1,AC2,'Both AC1 and AC2 could not get accepted. Aborting the fit.'
        close,1
        goto,finishthisgalaxy
     ENDIF

;!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!The current Code ends here!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                ;This is where the third estimate
                                ;starts in previous versions. Here we change to 0.5 beam
                                ;size radii and we try to determine
                                ;the start of the warp fitting those
                                ;parameters as 1 for greater accuracy
                                ;first the standard changes

     
     finishthisgalaxy:
   
     IF optimized then begin
        catcubename[i]= noptname[0]
     ENDIF
     cd, new_dir
     names=[catcubename[i],catMom0name[i],catMom1name[i],catmaskname[i],noisemapname,catCatalogname[i],basicinfo]
     book_keeping,names,bookkeeping
     
     IF size(log,/TYPE) EQ 7 then begin
        openu,66,log,/APPEND
        printf,66,linenumber()+"Finished "+catDirname[i]+" which is galaxy # "+strtrim(string(fix(i)),2)+" at "+systime()
        close,66
     ENDIF
     CD,old_dir 
  endfor
  CD,originaldir
  close,3
end

