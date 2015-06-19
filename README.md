# Fully Automated TiRiFiC (FAT)
=====

Introduction
------------

The Fully Automated TiRiFiC is an IDL wrapper around the tilted ring fitting code  (TiRiFiC) that aims to fully automize the process of fitting tilted ring models to line emission cubes. This code is still very much in the development phase and hence errors and bugs are to be expected. Nevertheless the code has extensively been tested and the results and a more extensive description of the code are documented in Kamphuis et al. 2015 

Requirements
------------
The code requires full installation of:

    IDL 7.0 or higher with astrolib
    TiRiFiC v2.2.3 or higher (http://gigjozsa.github.io/tirific/download_and_installation.html)
    SoFiA v 0.4.0 or higher  (https://github.com/SoFiA-Admin/SoFiA)
    Standard unix commands pwd, mkdir, rm, cp, ls, python, rename
    
IDL needs to be able to exececute tirific, sofia, rename and the standard unix commands from a spawn command. All other dependencies should be in IDL and available with the normal IDL distributions. 

Installation
------------

Unpack the zip file in a desired directory and you are ready to run FAT from this directory under IDL. 
The rename command might have to be aliased to rename -s, as this depends on which exact rename command and their are many versions available the code excutes the command "rename originalstring replacestring filesonwhichtoexecute" make sure that this is what the rename command does in the shell that is run by IDL.

You will also have to make a softlink in the Support directory to file sofia_pipeline.py in the sofia distribution i.e.:

        cd Support/
        ln -s pathtosofiainstallation/sofia_pipeline.py sofia_pipeline.py

Running FAT
-----------
FAT is currently run under IDL. It is called as a regular IDL program, i.e. in IDL:

    IDL >.r FAT_vx.x

where x.x is the version number

    IDL >FAT,configfile='pathtodir/configfile.config',support='pathtosupportfilesdir'
    
All information that the code needs about output directories fitting steps and input parameters are taken from the configfile.
If a config file is not given it will look for the file 'FAT_INPUT.config' in the directory from which FAT is run.
The default support directory is /Support however you can specify an alternative directory with the keyword support.

Configuration File
------

A configuration file will require the following parameters:

        catalogue=Path_to_catalog_dir/Catalog.txt

The code requires a catalog with input sources to know which cubes to fit and where they are (See Below). There is no default for this.

        maindir=Path_to_dir_with_input/

This is the path where the directories for all galaxies are stored. As FAT is still very much in a development state it still produces quite some output (e.g. Models for each step, xvdiagrams, Fitting logs). In order to keep this managable each galaxy requires its own directory. There is no default for this.

        outputcatalogue=Path/nameofresult.txt

The code will write a summary of the succes of the fit for each galaxy in this file. No default

        new_output='y'

Do you want to write a new output catalogue. If set to 'n'  the existing log file will be appendended. Default ='y'

        startgalaxy=0

The catalog number at which the code should start. The default is 0 which is the first line

        endgalaxy=-1

The catalog number at which the code should stop. The default is -1 which means that it should run until the end of the catalog.

        outputlog=fitlog.txt

The name of a log file that will trace the iterations and steps that the code is executing for each galaxy. This file is written into the galaxy directory. If left out no log file will be written and additional output will be printed to the terminal.

        new_log='y'

Do you want to write a new log file. If set to 'n'  the existing log file will be appendended. Default='y'

        velocity_resolution=1

The velocity resolution of the data cubes. If set to zero the code assume that the instrumental dispersion is equal to a (1.2 x channel)/(2 x SQRT(2ln2)) otherwise (1+vresolution) x channel/(2 x SQRT(2ln2)). That is, if set to 1 it assumes Hanning smoothing. Default=1.

        maps_output = 2

A parameter to control the amount of outpur created by FAT.  0.= all possible output (This is a lot), 1= all steps model + maps + def files, 2 = Final model + maps + def files for steps + logs, 3 = Only final model def + logs. Default = 2


The following parameters are predominantly set for testing the code and normally best left to their defaults, i.e. unset.

        testing=0

parameter to skip certain steps of the fitting for testing. If set to 1 it will skip the flat disk fitting and assume all necessary files for step 2 are in place. Default= 0

        allnew=1.

Parameter for setting the type of input for the initial guesses. possible setting are 0, 1, 2 
0) use initial guesse produced by FAT in a previous run. This is not recommended but can slightly speed up the start of the fitting. However it can lead to mismatches if slight changes have occured between runs.
1) FAT produces all initial input. Default
2) Pre-produced SoFiA output should be used for the initial guesses. These need to be specified in the input catalog (See input catalog)

        finishafter=2.

Parameter for finishing the fitting process early. if set to one the program finishes after fitting the flat disk. Default = 2

        opt_pixelbeam=4.
        
The amount of pixels in the FWHM of the minor axis. Default = 4.

    

A default config file (FAT_INPUT.config) is included in the distribution.

Input Catalog
-----------

The input catalog should have at least 4 columns named as 

        number|distance|directoryname|cubename

and seperated by |
The number is an easy identifier to keep track of which galaxy is being fitted.
the distance is the distance to the galaxy in Mpc. This is used to make some initial guesses for the structure of the galaxy. If it is unknown it should be set to 1.
The directory name is the name of the directory of the galaxy to be fitted. This directory should be located in the specified maindir in the config file.
cubename is the name of the cube to be fitted. This should be without the fits extension.

An example catalog is included in the distribution. This also gives examples for how to set up a cataog when using pre-made sofia input, i.e. allnew=2



