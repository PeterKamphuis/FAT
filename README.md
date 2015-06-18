# Fully Automated TiRiFiC (FAT)
=====

Introduction
------------

The Fully Automated TiRiFiC is an IDL wrapper around the tilted ring fitting code  (TiRiFiC) that aims to fully automize the process of fitting tilted ring models to line emission cubes. 

Requirements
------------
The code requires full installation of:

    IDL 7.0 or higher with astrolib
    TiRiFiC v2.2.3 or higher (http://gigjozsa.github.io/tirific/download_and_installation.html)
    SoFiA v 0.4.0 or higher  (https://github.com/SoFiA-Admin/SoFiA)

IDL needs to be able to exececute tirific and sofia from a spawn command. All other dependencies should be in IDL and available with the normal IDL distributions.

Installation
------------

Unpack the zip file in a desired directory and you are ready to run FAT from this directory under IDL. 

Running FAT
-----------
FAT is currently run under IDL. It is called as a regular IDL program, i.e. in IDL:

    IDL >.r FATvx.x
    IDL >FAT,configfile='hereismydir/configfile.config'
    
All information that the code needs about output directories fitting steps and input parameters are taken from the configfile.


  
