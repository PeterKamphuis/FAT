import numpy as np
import os.path
import glob
from cycler import cycler
from writenewtotemplate import writenewtotemplate
from convertradec import convertradec
from plot_fits import plot_fits
from collections import Counter
import matplotlib.pyplot as plt
import matplotlib as mpl
def overview_plot(distance=1.0, deffiles=None,galaxy='Plot'):
   """
    NAME:
          overview_plot
   
    PURPOSE:
          Create an overview plot to easily see how well the model fits the data.
   
    CATEGORY:
          Support
   
    CALLING SEQUENCE:
          OVERVIEW_PLOT
   
   
    INPUTS:
            distance = Distance to the galaxy
            gdlidl = an indicator whether we are running idl or gdl
    OPTIONAL INPUTS:
            Noise = the noise in the cube.
      finishafter = key for what kind of fitting was done.
      filenames   = nams of the created files
    KEYWORD PARAMETERS:
          /SPLINED - Plot with spline interpolation between the points.
   
    OUTPUTS:
   
    OPTIONAL OUTPUTS:
          -
   
    PROCEDURES CALLED:
          AXIS,COLORMAPS,COLOUR_BAR,COLUMNDENSITY,CONTOUR,DEVICE,FILE_TEST(),LOADCT,MAX(),N_ELEMENTS(),OPLOT,PLOT,FAT_PLOTERR,READFITS(),SET_PLOT,SHOWPIXELSMAP,WRITENEWTOTEMPLATE,XYOUTS
   
    MODIFICATION HISTORY:
          14-11-2017 P.Kamphuis; Increased and improved checking for convert
          22-03-2017 P.Kamphuis; Increased checking for Imagick convert
                                 in order to trim the png. Improved
                                 plotting limits on velocity contours
                                 and PV-Diagram color scale.
          15-11-2016 P.Kamphuis; Fixed a bounding Box issue in GDL
          15-11-2016 P.Kamphuis; Introduced hexadecimal plotting for
                                 color in idl and gdl. In idl and gdl
                                 this is coded as RGB colors converting
                                 to hex(B)+hex(G)+hex(R).
          08-11-2016 P.Kamphuis; Added the splined keyword to have the
                                 variables plotted in a spline
                                 interpolation
          08-11-2016 P.Kamphuis; Removed the plotting through the X
                                 widget for GDL as it cannot be run in
                                 screen in that case. However this
                                 results in black and white plots. Need
                                 to move to hexidecimal plotting
          01-08-2016 P.Kamphuis; Updated the GDL plotting to use png as
                                 well. This now happens through the X
                                 buffer which unfortunately means that
                                 it is limited to the amount of pixels
                                 in the users screen. This is because
                                 the Z-buffer in GDL does not accept the
                                 decomposed keyword.
          09-06-2016 P.Kamphuis; Case matched the call to Finalmodel as
                                 unix is case sensitive.
          23-03-2016 P.Kamphuis; Added filename and finishafter
                                 extension to make sure that the right files are picked up and
                                 to provide info about the fit
          Written 20-02-2016 P.Kamphuis v1.0
   
    NOTE:
   
   """
   plotpara = ['RADI', 'SBR', 'SBR_2', 'VROT', 'VROT_ERR', 'PA', 'PA_ERR', 'PA_2', 'PA_2_ERR', 'INCL', 'INCL_ERR', 'INCL_2', 'INCL_2_ERR', 'BMAJ', 'SDIS', 'XPOS', 'YPOS', 'VSYS']
   plotstart = np.array([[1, 3, 5, 9], [2, 3, 7, 11],[0, 1, 1, 1]])
   if len(deffiles) == 1:
      arrays = writenewtotemplate(Newfilename=deffiles[0], Variablechange=plotpara, Extract=True)
   else:
      for i in range(len(deffiles)):
         print(deffiles[i])
         tmp = writenewtotemplate(Newfilename=deffiles[i], Variablechange=plotpara, Extract=True)
         if i == 0:
            arrays = tmp
         elif i == 1:
            tmp2 = arrays
            arrays= np.zeros([max([len(tmp[:,0]),len(tmp2[:,0])]),len(plotpara),i+1])
            arrays[0:len(tmp2[:,0]),:,0] = tmp2[:,:]
            arrays[0:len(tmp[:,0]),:,1]=tmp[:,:]
         else:
            tmp2 = arrays
            arrays= np.zeros([max([len(tmp[:,0]),len(tmp2[:,0,0])]),len(plotpara),i+1])
            arrays[0:len(tmp2[:,0]),:,0:i] = tmp2[:,:,:]
            arrays[0:len(tmp[:,0]),:,i]=tmp[:,:]
         
   if os.path.isfile('ModelInput.def'):               
      modarrays= writenewtotemplate(Newfilename='ModelInput.def', Variablechange=plotpara, Extract=True)
   varunits=[]
   for varia in plotpara:
      tmp = varia.split('_')
      if tmp[0] == 'VROT' or tmp[0] == 'SDIS':
         varunits.append('(km s^{-1})')
      elif tmp[0] == 'PA' or tmp[0] == 'INCL':   
         varunits.append('(Degrees)')
      elif tmp[0] == 'SBR':   
         varunits.append('(Jy km s!E-1!N arcsec!E-2!N)')
      elif tmp[0] == 'Z0':   
         varunits.append('(Arcsec)')
      else:   
         varunits.append('')
  # minvar = np.array([100000.,100000.,100000.,100000.])
  # maxvar = np.array([-100000.,-100000.,-100000.,-100000.])
  # buffer = np.zeros(4)
  # for i in np.arange(4):
  #    try:
  #       tmpvals = np.array([arrays[:,plotstart[0,i],:], arrays[:,plotstart[1,i],:]])
  #    except IndexError:
  #       tmpvals = np.array([arrays[:,plotstart[0,i]], arrays[:,plotstart[1,i]]])
  #    tmplocs = np.where(tmpvals != 0.)
  #    if tmpvals[tmplocs].max() > maxvar[i]:   
  #       maxvar[i] = tmpvals[tmplocs].max()
  #    if tmpvals[tmplocs].min() < minvar[i]:   
  #       minvar[i] = tmpvals[tmplocs].min()
  #    buffer[i] = (abs(maxvar[i]) + abs(minvar[i])) / 20.
   try:
      ra = arrays[0,-3,:]
      dec = arrays[0,-2,:]
      vsys=["{:10.1f}".format(itm) for itm in arrays[0,-1,:]]
      disper=["{:10.1f}".format(itm) for itm in arrays[0,-4,:]]
      majbeam=["{:10.1f}".format(itm) for itm in arrays[0,-5,:]]
      indexgh = np.where(arrays[:,0,:] != 0.)
      last = 1
      while len(np.where(indexgh[:][0] == last)[0]) == len(vsys):
         last += 1  
      last -= 1
      ringsize=["{:10.1f}".format(itm) for itm in list(arrays[last,0,:]-arrays[last-1,0,:])]
      ceninc=arrays[0,plotpara.index("INCL"),:]
      maxvrot=np.max(arrays[:,plotpara.index("VROT"),:],axis=0)     
   except IndexError:
      ra = arrays[0,-3]
      dec = arrays[0,-2]
      vsys="{:10.1f}".format(arrays[0,-1])
      disper="{:10.1f}".format(arrays[0,-4])
      majbeam="{:10.1f}".format(arrays[0,-5])
      ringsize="{:10.1f}".format(arrays[-1,0]-arrays[-2,0])
      ceninc=arrays[0,plotpara.index("INCL")]
      maxvrot=np.max(arrays[:,plotpara.index("VROT")]) 
   ra,dec = convertradec(ra,dec)
   # create a figure
 #  mpl.font_manager._rebuild()
   mpl.rcParams['axes.prop_cycle'] = cycler('color',['k','k','r', 'r','b', 'b', 'g','g'])
   mpl.rc('lines', linewidth=1)
   mpl.rc('xtick', color='k', labelsize='medium', direction='in')
   mpl.rc('ytick', color='k', labelsize='medium', direction='in')
   labelfont= {'family':'Times New Roman',
               'weight':'normal',
               'size':12}
   mpl.rc('font',**labelfont) 
   overview = plt.figure(2, figsize=(8.27, 11.69), dpi=300, facecolor='w', edgecolor='k')
   sbr = overview.add_axes([0.1, 0.5, 0.6, 0.125]) #left, bottom, width, height
   #parameters = overview.add_subplot(4, 1, 4)
    
   # first we plot the paramters for the galaxy
   try:
      tmp = arrays[0,-3,:]
      for i in range(len(vsys)):
         nozero=np.where(arrays[:,plotpara.index("SBR"),i] != 0)
         sbr.plot(arrays[nozero[0],plotpara.index("RADI"),i],arrays[nozero[0],plotpara.index("SBR"),i],'-o')
         sbr.plot(arrays[nozero[0],plotpara.index("RADI"),i],arrays[nozero[0],plotpara.index("SBR_2"),i],'--o')
   except IndexError:
      sbr.plot(arrays[:,plotpara.index("RADI")],arrays[:,plotpara.index("SBR")],'-ko')
      sbr.plot(arrays[:,plotpara.index("RADI")],arrays[:,plotpara.index("SBR_2")],'--ko')
   ymin, ymax = plt.ylim()
   bufferpl=(abs(ymax)-abs(ymin))/20.
   sbr.set_ylim(ymin-bufferpl,ymax+bufferpl)
   sbr.set_xlabel('Radius (arcsec)',**labelfont)
   sbr.set_ylabel('SBR (Jy km s$^{-1}$ arcsec$^{-2}$)',**labelfont)
   vrot = overview.add_axes([0.1, 0.625, 0.6, 0.125])
   miniy = 1000.
   try:
      tmp = arrays[0,-3,:]
      for i in range(len(vsys)):
         nozero=np.where(arrays[:,plotpara.index("VROT"),i] != 0)[0]
         vrot.plot(arrays[0:nozero[-1],plotpara.index("RADI"),i],arrays[0:nozero[-1],plotpara.index("VROT"),i],'-o')
         vrot.plot(arrays[0:nozero[-1],plotpara.index("RADI"),i],arrays[0:nozero[-1],plotpara.index("VROT"),i],'--o')
         tmp = min(arrays[nozero,plotpara.index("VROT"),i])
         if tmp < miniy:
            miniy = tmp
   except IndexError:
      vrot.plot(arrays[:,plotpara.index("RADI")],arrays[:,plotpara.index("VROT")],'-o')
      nozero=np.where(arrays[:,plotpara.index("VROT")] != 0)[0]
      tmp = min(arrays[nozero,plotpara.index("VROT")])
      if tmp < miniy:
         miniy = tmp
  
   ymin, ymax = plt.ylim()
   ymin= miniy
   bufferpl=(abs(ymax)-abs(ymin))/20.
   vrot.set_ylim(ymin-bufferpl,ymax+bufferpl)
   vrot.set_ylabel('V$_{rot}$ (km s$^{-1}$)',**labelfont)
   vrot.set_xticklabels([])

   pa = overview.add_axes([0.1, 0.75, 0.6, 0.125])
   try:
      tmp = arrays[0,-3,:]
      for i in range(len(vsys)):
         nozero=np.where(arrays[:,plotpara.index("PA"),i] != 0)
         pa.plot(arrays[nozero[0],plotpara.index("RADI"),i],arrays[nozero[0],plotpara.index("PA"),i],'-o')
         pa.plot(arrays[nozero[0],plotpara.index("RADI"),i],arrays[nozero[0],plotpara.index("PA_2"),i],'--o')
   except IndexError:
      pa.plot(arrays[:,plotpara.index("RADI")],arrays[:,plotpara.index("PA")],'-ko')
      pa.plot(arrays[:,plotpara.index("RADI")],arrays[:,plotpara.index("PA_2")],'--ko')
   ymin, ymax = plt.ylim()
   bufferpl=(abs(ymax)-abs(ymin))/20.
   pa.set_ylim(ymin-bufferpl,ymax+bufferpl)
   pa.set_ylabel('PA (Degrees)',**labelfont)
   pa.set_xticklabels([])

   incl = overview.add_axes([0.1, 0.875, 0.6, 0.125])
   try:
      tmp = arrays[0,-3,:]
      for i in range(len(vsys)):
         nozero=np.where(arrays[:,plotpara.index("INCL"),i] != 0)
         incl.plot(arrays[nozero[0],plotpara.index("RADI"),i],arrays[nozero[0],plotpara.index("INCL"),i],'-o')
         incl.plot(arrays[nozero[0],plotpara.index("RADI"),i],arrays[nozero[0],plotpara.index("INCL_2"),i],'--o')
   except IndexError:
      incl.plot(arrays[:,plotpara.index("RADI")],arrays[:,plotpara.index("INCL")],'-ko')
      incl.plot(arrays[:,plotpara.index("RADI")],arrays[:,plotpara.index("INCL_2")],'--ko')
   ymin, ymax = plt.ylim()
   bufferpl=(abs(ymax)-abs(ymin))/20.
   incl.set_ylim(ymin-bufferpl,ymax+bufferpl)
   incl.set_ylabel('INCL (Degrees)',**labelfont)
   incl.set_xticklabels([])
   
   filecolors=['black','red','blue','green']
   # then we write the different filese central coordinates and some input parameters
   try:
      tmp = arrays[0,-3,:]
      for i in range(len(deffiles)):
         plt.gcf().text(0.72, 0.98-i*0.1, "This plot shows {} in {}".format(deffiles[i],filecolors[i]) , fontsize=12)
         plt.gcf().text(0.72, 0.95-i*0.1, "With RA= {}, DEC = {}, V$_{{sys}}$ = {} ".format(ra[i],dec[i],vsys[i]) , fontsize=12)
         plt.gcf().text(0.72, 0.92-i*0.1,"The Dispersion = {} km s$^{{-1}}$, FWHM = {} arcsec and the ringsize = {} arcsec".format(disper[i],majbeam[i], ringsize[i]), fontsize=12)
         wrkdir=deffiles[i].split('/')[0]
         mom0files=glob.glob(wrkdir+'/Moments/*preprocessed_mom0*.fits')
         mom0model=glob.glob(wrkdir+'/Moments/Finalmodel_mom0.fits')
         plot_fits(mom0files[0],overview,figure_coordinates=[0.05+(i*(0.3*11.69/8.27)+i*0.075),0.2,0.3*11.69/8.27,0.3*8.27/11.69],cmap='hot_r',contours=mom0model[0],log=True)
         plt.gcf().text(0.1+(i*(0.3*11.69/8.27)+i*0.075)+0.05*11.69/8.27, 0.21+0.3*8.27/11.69, "Mom0 {}".format(wrkdir) , fontsize=12)
         mom1files=glob.glob(wrkdir+'/Moments/*preprocessed_mom1*.fits')
         mom1model=glob.glob(wrkdir+'/Moments/Finalmodel_mom1.fits')
         plot_fits(mom1files[0],overview,figure_coordinates=[0.05+(i*(0.3*11.69/8.27)+i*0.075),-0.1,0.3*11.69/8.27,0.3*8.27/11.69],cmap='jet',contours=mom1model[0])
         plt.gcf().text(0.1+(i*(0.3*11.69/8.27)+i*0.075)+0.05*11.69/8.27, -0.09+0.3*8.27/11.69, "Mom1 {}".format(wrkdir) , fontsize=12)
         xvfiles=glob.glob(wrkdir+'/PV-Diagrams/*_2_*.fits')
         print(xvfiles)
         xvmodel=glob.glob(wrkdir+'/PV-Diagrams/Finalmodel_xv.fits')
         plot_fits(xvfiles[0],overview,figure_coordinates=[0.1+(i*(0.3*11.69/8.27)+i*0.1),-0.4,0.275*11.69/8.27,0.3*8.27/11.69],cmap='hot_r',contours=xvmodel[0],aspect='auto')
         plt.gcf().text(0.1+(i*(0.3*11.69/8.27)+i*0.075)+0.05*11.69/8.27, -0.39+0.3*8.27/11.69, "PV-Diagram {}".format(wrkdir) , fontsize=12)
   except IndexError:
      i = 1
      plt.gcf().text(0.72, 0.98-i*0.1, "This plot shows {} in {}".format(deffiles,filecolors[0]) , fontsize=12)
      plt.gcf().text(0.72, 0.95-i*0.1, "With RA= {}, DEC = {}, V$_{{sys}}$ = {} ".format(ra,dec,vsys) , fontsize=12)
      plt.gcf().text(0.72, 0.92-i*0.1,"The Dispersion = {} km s$^{{-1}}$, FWHM = {} arcsec and the ringsize = {} arcsec".format(disper,majbeam, ringsize), fontsize=12)
      wrkdir=deffiles[0].split('/')[0]
      mom0files=glob.glob(wrkdir+'/Moments/*preprocessed_mom0*.fits')
      mom0model=glob.glob(wrkdir+'/Moments/Finalmodel_mom0.fits')
      plot_fits(mom0files[0],overview,figure_coordinates=[0.05+(i*(0.3*11.69/8.27)+i*0.075),0.2,0.3*11.69/8.27,0.3*8.27/11.69],cmap='hot_r',contours=mom0model[0],log=True)
      plt.gcf().text(0.1+(i*(0.3*11.69/8.27)+i*0.075)+0.05*11.69/8.27, 0.21+0.3*8.27/11.69, "Mom0 {}".format(wrkdir) , fontsize=12)
      mom1files=glob.glob(wrkdir+'/Moments/*preprocessed_mom1*.fits')
      mom1model=glob.glob(wrkdir+'/Moments/Finalmodel_mom1.fits')
      plot_fits(mom1files[0],overview,figure_coordinates=[0.05+(i*(0.3*11.69/8.27)+i*0.075),-0.1,0.3*11.69/8.27,0.3*8.27/11.69],cmap='jet',contours=mom1model[0])
      plt.gcf().text(0.1+(i*(0.3*11.69/8.27)+i*0.075)+0.05*11.69/8.27, -0.09+0.3*8.27/11.69, "Mom1 {}".format(wrkdir) , fontsize=12)
      xvfiles=glob.glob(wrkdir+'/PV-Diagrams/*_2_*.fits')
      print(xvfiles)
      xvmodel=glob.glob(wrkdir+'/PV-Diagrams/Finalmodel_xv.fits')
      plot_fits(xvfiles[0],overview,figure_coordinates=[0.1+(i*(0.3*11.69/8.27)+i*0.1),-0.4,0.275*11.69/8.27,0.3*8.27/11.69],cmap='hot_r',contours=xvmodel[0],aspect='auto')
      plt.gcf().text(0.1+(i*(0.3*11.69/8.27)+i*0.075)+0.05*11.69/8.27, -0.39+0.3*8.27/11.69, "PV-Diagram {}".format(wrkdir) , fontsize=12)
    
 
   #finally we must plot the different velocity fields and PV-Diagrams

    
#   plt.show()
   plt.savefig(galaxy+'_Overview.png', bbox_inches='tight')
   plt.close()
   """0.3*11.69/8.27
   
   
   device(set_font='Times', tt_font=True, set_resolution=concatenate([scrdim[0], scrdim[1]]), decomposed=True, set_pixel_depth=24)
   plotradii = arrays[0,:]
   tmp = where(ravel(arrays[1,:] > 1.1e-16))[0]
   tmp2 = where(ravel(arrays[2,:] > 1.1e-16))[0]
   
   maxradii = max(concatenate([plotradii[tmp], plotradii[tmp2]])) + (plotradii[array(plotradii, copy=0).nelements() - 1] - plotradii[array(plotradii, copy=0).nelements() - 2]) / 2.
   
   for i in arange(0, 4):
      if i == 0:   
         
         plotvariable = arrays[1,tmp]
         loadct(0, silent=True)
         plot(plotradii, plotvariable, position=concatenate([0.15, 0.9 - 4 * ysize, 0.55, 0.9 - 3 * ysize]), xtitle='Radius (arcmin)', xrange=concatenate([0., maxradii]), yrange=concatenate([minvar[i] - buffer[i], maxvar[i] + buffer[i]]), ytickname=concatenate([' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ']), xtickname=concatenate([' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ']), xticklayout=1, background=0xffffff, color=0xffffff, nodata=True)
         if (splined is not None):   
            newrad = dblarr((array(plotradii, copy=0).nelements() - 1) * 10. + 1)
            for h in arange(1, (array(plotradii, copy=0).nelements() - 1)+(1)):
               newrad[(h - 1) * 10:(h * 10 - 1)+1] = findgen(10) * (plotradii[h] - plotradii[h - 1]) / 10. + plotradii[h - 1]
            newrad[(h - 1) * 10] = plotradii[h - 1]
            newvar = spline(plotradii, plotvariable, newrad)
            oplot(newrad, newvar, color=0x000000, linestyle=0, symsize=ssize)
         else:   
            oplot(plotradii, plotvariable, color=0x000000, linestyle=0, symsize=ssize)
         levelsrange = concatenate([minvar[i] - buffer[i], maxvar[i] + buffer[i]]) * 1000.
         oplot(plotradii, plotvariable, psym=8, color=0x000000, linestyle=2, symsize=ssize)
         columndensity(levelsrange, arrays[14,0], concatenate([1., 1.]), vwidth=1., arcsquare=True)
         levelsranges = levelsrange / 1e20
         reset = 1e20
         levelsranges[0] = ceil(levelsranges[0])
         levelsranges[1] = floor(levelsranges[1])
         adst = 'x10!E20!N'
         if levelsranges[0] == levelsranges[1]:   
            levelsranges = levelsrange / 1e19
            levelsranges[0] = ceil(levelsranges[0])
            levelsranges[1] = floor(levelsranges[1])
            adst = 'x10!E19!N'
            reset = 1e19
         midrange = levelsranges[0] + (levelsranges[1] - levelsranges[0]) / 2.
         if array(midrange, copy=0).astype(Int32) != midrange:   
            levelsranges[1] = levelsranges[1] - 1
            midrange = levelsranges[0] + (levelsranges[1] - levelsranges[0]) / 2.
         newlevels = concatenate([levelsranges[0], midrange, levelsranges[1]])
         jynewlevels = newlevels * reset
         columndensity(jynewlevels, array(vsys, copy=0).astype(Float64), concatenate([1., 1.]), vwidth=1., ncolumn=True, arcsquare=True)
         axis(yaxis=0, charthick=charthick, xthick=xthick, ythick=ythick, charsize=charsize, color=0x000000)
         axis(yaxis=1, charthick=charthick, xthick=xthick, ythick=ythick, charsize=charsize, ytickv=jynewlevels / 1000., ytickname=concatenate([strtrim(strcompress(string(newlevels[0], format='(I3)')), 2), strtrim(strcompress(string(array(newlevels[1], copy=0).astype(Int32), format='(I2)')), 2), strtrim(strcompress(string(array(newlevels[2], copy=0).astype(Int32), format='(I2)')), 2)]), yticks=2, yminor=3, color=0x000000)
         xyouts(0.60, 0.9 - 3.5 * ysize, 'N!IH', normal=True, alignment=0.5, orientation=90, charthick=charthick, charsize=_sys_p.charsize * 1.25, color=0x000000)
         xyouts(0.63, 0.9 - 3.5 * ysize, '(' + adst + ' cm!E-2!N)', normal=True, alignment=0.5, orientation=90, charthick=charthick, charsize=charsize, color=0x000000)
         axis(xaxis=0, charthick=charthick, xthick=xthick, ythick=ythick, charsize=charsize, color=0x000000, xtitle='Radius (arcsec)')
         axis(xaxis=1, charthick=charthick, xthick=xthick, ythick=ythick, charsize=charsize, xrange=convertskyanglefunction(_sys_x.crange, distance), xtickname=concatenate([' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ']), color=0x000000)
         loadct(40, silent=True)
         if (splined is not None):   
            newvar = spline(plotradii, arrays[2,tmp2], newrad)
            oplot(newrad, newvar, color=0x0000ff, linestyle=2, symsize=ssize)
         else:   
            oplot(plotradii, arrays[2,tmp2], thick=lthick, color=0x0000ff, linestyle=2)
         oplot(plotradii, arrays[2,tmp2], psym=8, color=0x0000ff, linestyle=2, symsize=ssize)
         xyouts(0.05, 0.9 - 3.5 * ysize, plotpara[plotstart[0,i]], normal=True, alignment=0.5, orientation=90, charsize=_sys_p.charsize * 1.25, color=0x000000, charthick=charthick)
         xyouts(0.08, 0.9 - 3.5 * ysize, varunits[plotstart[0,i]], normal=True, alignment=0.5, orientation=90, color=0x000000, charthick=charthick)
         
         if file_test('ModelInput.def'):   
            oplot(modarrays[0,:], modarrays[1,:], thick=lthick, color=0xff0010)
            oplot(modarrays[0,:], modarrays[1,:], psym=8, color=0xff0010, symsize=ssize)
            oplot(modarrays[0,:], modarrays[2,:], thick=lthick, color=0x00b4ff, linestyle=2)
            oplot(modarrays[0,:], modarrays[2,:], psym=8, color=0x00b4ff, linestyle=2, symsize=ssize)
      else:   
         if plotstart[0,i] != plotstart[1,i]:   
            plotvariable = arrays[plotstart[0,i],tmp]
            plotvariableerr = arrays[plotstart[0,i] + plotstart[2,i],tmp]
         else:   
            if tmp[array(tmp, copy=0).nelements() - 1] > tmp2[array(tmp2, copy=0).nelements() - 1]:   
               plotvariable = arrays[plotstart[0,i],tmp]
               plotvariableerr = arrays[plotstart[0,i] + plotstart[2,i],tmp]
            else:   
               plotvariable = arrays[plotstart[0,i],tmp2]
               plotvariableerr = arrays[plotstart[0,i] + plotstart[2,i],tmp2]
         loadct(0, silent=True)
         if total(plotvariableerr) != 0.:   
            xerr = dblarr(array(plotvariableerr, copy=0).nelements())
            fat_ploterror(plotradii, plotvariable, xerr, plotvariableerr, position=concatenate([0.15, 0.9 - (4 - i) * ysize, 0.55, 0.9 - (3 - i) * ysize]), xrange=concatenate([0., maxradii]), yrange=concatenate([minvar[i] - buffer[i], maxvar[i] + buffer[i]]), xthick=xthick, ythick=ythick, xtickname=concatenate([' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ']), xticklayout=1, charthick=charthick, thick=thick, charsize=charsize, linestyle=0, noerase=True, color=0x000000, errcolor=0x000000, errthick=_sys_p.thick * 0.4, psym=8, symsize=ssize)
         else:   
            plot(plotradii, plotvariable, position=concatenate([0.15, 0.9 - (4 - i) * ysize, 0.55, 0.9 - (3 - i) * ysize]), xrange=concatenate([0., maxradii]), yrange=concatenate([minvar[i] - buffer[i], maxvar[i] + buffer[i]]), xthick=xthick, ythick=ythick, xtickname=concatenate([' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ']), xticklayout=1, charthick=charthick, thick=thick, charsize=charsize, linestyle=0, noerase=True, color=0x000000, psym=8, symsize=ssize)
         if (splined is not None):   
            newvar = spline(plotradii, plotvariable, newrad)
            oplot(newrad, newvar, color=0x000000, linestyle=0, symsize=ssize)
         else:   
            oplot(plotradii, plotvariable, thick=lthick, color=0x000000, linestyle=0)
         
         if i == 3:   
            axis(xaxis=1, charthick=charthick, xthick=xthick, ythick=ythick, charsize=charsize, xrange=convertskyanglefunction(_sys_x.crange, distance), xtitle='Radius (kpc)', color=0x000000)
         else:   
            axis(xaxis=1, charthick=charthick, xthick=xthick, ythick=ythick, charsize=charsize, xrange=convertskyanglefunction(_sys_x.crange, distance), xtickname=concatenate([' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ']), color=0x000000)
         xyouts(0.05, 0.9 - (3.5 - i) * ysize, plotpara[plotstart[0,i]], normal=True, alignment=0.5, orientation=90, charsize=_sys_p.charsize * 1.25, color=0x000000, charthick=charthick)
         xyouts(0.08, 0.9 - (3.5 - i) * ysize, varunits[plotstart[0,i]], normal=True, alignment=0.5, orientation=90, color=0x000000, charthick=charthick)
         loadct(40, silent=True)
         if plotstart[0,i] != plotstart[1,i]:   
            plotvariable = arrays[plotstart[1,i],tmp2]
            plotvariableerr = arrays[plotstart[1,i] + plotstart[2,i],tmp2]
            if total(plotvariableerr) != 0.:   
               xerr = dblarr(array(plotvariableerr, copy=0).nelements())
               fat_ploterror(plotradii, plotvariable, xerr, plotvariableerr, thick=lthick, color=0x0000ff, linestyle=2, errcolor=0x0000ff, errthick=_sys_p.thick * 0.4, over_plot=True, psym=8, symsize=ssize)
            else:   
               oplot(plotradii, plotvariable, thick=lthick, color=0x0000ff, linestyle=2, psym=8, symsize=ssize)
            if (splined is not None):   
               newvar = spline(plotradii, plotvariable, newrad)
               oplot(newrad, newvar, color=0x0000ff, linestyle=2, symsize=ssize)
            else:   
               oplot(plotradii, plotvariable, thick=lthick, color=0x0000ff, linestyle=2)
            
         if file_test('ModelInput.def'):   
            oplot(modarrays[0,:], modarrays[plotstart[0,i],:], thick=lthick, color=0xff0010)
            oplot(modarrays[0,:], modarrays[plotstart[0,i],:], psym=8, color=0xff0010, symsize=ssize)
            if plotstart[0,i] != plotstart[1,i]:   
               oplot(modarrays[0,:], modarrays[plotstart[1,i],:], thick=lthick, color=0x00b4ff, linestyle=2)
               oplot(modarrays[0,:], modarrays[plotstart[1,i],:], psym=8, color=0x00b4ff, linestyle=2, symsize=ssize)
   if file_test('ModelInput.def'):   
      ramod = array(modarrays[array(plotpara, copy=0).nelements() - 3,0], copy=0).astype(Float64)
      decmod = array(modarrays[array(plotpara, copy=0).nelements() - 2,0], copy=0).astype(Float64)
      convertradec(ramod, decmod)
      vsysmod = strtrim(string(array(modarrays[array(plotpara, copy=0).nelements() - 1,0], copy=0).astype(Float64), format='(F10.1)'), 2)
      dispermod = strtrim(string(array(modarrays[array(plotpara, copy=0).nelements() - 4,0], copy=0).astype(Float64), format='(F10.1)'), 2)
      xyouts(0.60, 0.89, 'Systemic Velocity= ' + vsys + ' (' + vsysmod + ') km s!E-1', normal=True, alignment=0., charthick=charthick, color=0x000000)
      xyouts(0.60, 0.87, 'R.A.= ' + ra + ' (' + ramod + ')', normal=True, alignment=0., charthick=charthick, color=0x000000)
      xyouts(0.60, 0.85, 'DEC.= ' + dec + ' (' + decmod + ')', normal=True, alignment=0., charthick=charthick, color=0x000000)
      xyouts(0.60, 0.83, 'Dispersion= ' + disper + ' (' + dispermod + ') km s!E-1', normal=True, alignment=0., charthick=charthick, color=0x000000)
      xyouts(0.60, 0.81, 'Black lines: approaching side parameters.', normal=True, alignment=0., charthick=charthick, color=0x000000)
      xyouts(0.60, 0.79, 'Red lines: receding side parameters.', normal=True, alignment=0., charthick=charthick, color=0x000000)
      xyouts(0.60, 0.77, 'Blue lines: approaching side input model parameters.', normal=True, alignment=0., charthick=charthick, color=0x000000)
      xyouts(0.60, 0.75, 'Yellow lines: receding side input model parameters.', normal=True, alignment=0., charthick=charthick, color=0x000000)
      if bitwise_or(array(finishafter, copy=0).astype(Int32) / finishafter != 1, finishafter == 1):   
         xyouts(0.60, 0.73, 'The inclination and PA were not allowed to vary', normal=True, color=0x000000)
         xyouts(0.60, 0.71, 'The major FWHM beam is ' + majbeam + ' arcsec', normal=True, color=0x000000)
         xyouts(0.60, 0.69, 'We used rings of size ' + ringsize + ' arcsec', normal=True, color=0x000000)
      else:   
         xyouts(0.60, 0.73, 'The major FWHM beam is ' + majbeam + ' arcsec', normal=True, color=0x000000)
         xyouts(0.60, 0.71, 'We used rings of size ' + ringsize + ' arcsec', normal=True, color=0x000000)
   else:   
      xyouts(0.60, 0.89, 'Systemic Velocity= ' + vsys + ' km s!E-1', normal=True, alignment=0., charthick=charthick, color=0x000000)
      xyouts(0.60, 0.87, 'R.A.= ' + ra, normal=True, alignment=0., charthick=charthick, color=0x000000)
      xyouts(0.60, 0.85, 'DEC.= ' + dec, normal=True, alignment=0., charthick=charthick, color=0x000000)
      xyouts(0.60, 0.83, 'Dispersion= ' + disper + ' km s!E-1', normal=True, alignment=0., charthick=charthick, color=0x000000)
      xyouts(0.60, 0.81, 'Black lines: approaching side parameters.', normal=True, alignment=0., charthick=charthick, color=0x000000)
      xyouts(0.60, 0.79, 'Red lines: receding side parameters.', normal=True, alignment=0., charthick=charthick, color=0x000000)
      if bitwise_or(array(finishafter, copy=0).astype(Int32) / finishafter != 1, finishafter == 1):   
         xyouts(0.60, 0.77, 'The inclination and PA were not allowed to vary', normal=True, charthick=charthick, color=0x000000)
         xyouts(0.60, 0.75, 'The major FWHM beam is ' + majbeam + ' arcsec', normal=True, charthick=charthick, color=0x000000)
         xyouts(0.60, 0.73, 'We used rings of size ' + ringsize + ' arcsec', normal=True, charthick=charthick, color=0x000000)
      else:   
         xyouts(0.60, 0.77, 'The major FWHM beam is ' + majbeam + ' arcsec', normal=True, charthick=charthick, color=0x000000)
         xyouts(0.60, 0.75, 'We used rings of size ' + ringsize + ' arcsec', normal=True, charthick=charthick, color=0x000000)
   #Currently GDL does not recognize true
   #type fonts yet. This leads to errors
   #in using the degree symbol. It also
   #does not yet recognize superscript
   #commands in tickmarks.
   tmpstr = strtrim(strsplit(filenames[1], '/', extract=True), 2)
   if strupcase(tmpstr[0]) == 'MOMENTS':   
      mom0 = readfits(filenames[1] + '.fits', mom0hed, silent=True)
   else:   
      mom0 = readfits('Moments/' + filenames[1] + '.fits', mom0hed, silent=True)
   mom0mod = readfits('Moments/Finalmodel_mom0.fits', mom0hedmod, silent=True)
   mapmax = max(mom0, min=mapmin)
   buildaxii(mom0hed, xaxis, yaxis)
   colormaps('heat', invert=True)
   showpixelsmap(xaxis, yaxis, mom0, position=concatenate([0.15, 0.1, 0.35, 0.1 + 0.2 * scrdim[0] / scrdim[1]]), wcs=True, xtitle='RA (J2000)', ytitle='DEC (J2000)', blank_value=0., range=concatenate([0., mapmax]), noerase=True, charthick=charthick, thick=thick, black=True, hex_color=True)
   levels = concatenate([1e20, 2e20, 4e20, 8e20, 16e20, 32e20])
   beam = concatenate([sxpar(mom0hed, 'BMAJ') * 3600., sxpar(mom0hed, 'BMIN') * 3600.])
   if sxpar(mom0hed, 'BPA'):   
      bpa = sxpar(mom0hed, 'BPA')
   else:   
      bpa = 0
   columndensity(levels, array(vsys, copy=0).astype(Float64), beam, vwidth=1., ncolumn=True)
   levels = levels / 1000.
   loadct(0, silent=True)
   contour(mom0, xaxis, yaxis, levels=levels, overplot=True, c_colors=concatenate([0x000000]))
   loadct(40, silent=True)
   contour(mom0mod, xaxis, yaxis, levels=levels, overplot=True, c_colors=concatenate([0x0000ff]))
   beam_plot(beam[0], beam[1], bpa=bpa, center=concatenate([xaxis[0] - beam[0] / 3600., yaxis[0] + beam[0] / 3600.]), fill=True, color=0x000000, transparent=True)
   colormaps('heat', invert=True)
   colour_bar(concatenate([0.37, 0.39]), concatenate([0.12, 0.1 + 0.2 * scrdim[0] / scrdim[1] - 0.02]), strtrim(string(0, format='(F10.1)'), 2), strtrim(string(mapmax, format='(F10.1)'), 2), opposite_label=True, black=True, title='(Jy bm!E-1!N x km s!E-1!N)', vertical=True, charthick=charthick, hex_color=True)
   loadct(0, silent=True)
   xyouts(0.45, 0.01 + 0.2 * scrdim[0] / scrdim[1], 'Velocity Field, Moment0 and PV-Diagram along the major axis.', color=0x000000, normal=True, charthick=charthick)
   xyouts(0.45, 0.01 + 0.2 * scrdim[0] / scrdim[1] - 0.02, 'Black Contours: Data, White/Red Contours: Final Model', color=0x000000, normal=True, charthick=charthick)
   xyouts(0.45, 0.01 + 0.2 * scrdim[0] / scrdim[1] - 0.04, 'Moment 0 Contours are at 1, 2, 4, 8, 16, 32 x 10!E20!N cm!E-2', color=0x000000, normal=True, charthick=charthick)
   
   #Velocity Field
   tmpstr = strtrim(strsplit(filenames[2], '/', extract=True), 2)
   if strupcase(tmpstr[0]) == 'MOMENTS':   
      mom0 = readfits(filenames[2] + '.fits', mom0hed, silent=True)
   else:   
      mom0 = readfits('Moments/' + filenames[2] + '.fits', mom0hed, silent=True)
   
   mom0mod = readfits('Moments/Finalmodel_mom1.fits', mom0hedmod, silent=True)
   tmp = where(ravel(finite(mom0mod)))[0]
   #Too low inclination can result in a too small range
   ceninc = ceninc + 12.5
   if ceninc > 90:   
      ceninc = 90.
   velext = 1.25 * maxvrot * sin((ceninc) * _sys_dtor)
   mapmax = array(vsys[0] + velext[0], copy=0).astype(Float64)
   mapmin = array(vsys[0] - velext[0], copy=0).astype(Float64)
   buildaxii(mom0hed, xaxis, yaxis)
   colormaps('sauron_colormap')
   showpixelsmap(xaxis, yaxis, mom0, position=concatenate([0.15, 0.1 + 0.2 * scrdim[0] / scrdim[1], 0.35, 0.1 + 0.4 * scrdim[0] / scrdim[1]]), wcs=True, ytitle='DEC (J2000)', blank_value=0., range=concatenate([mapmin, mapmax]), noerase=True, charthick=charthick, thick=thick, xtickname=concatenate([' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ']), hex_color=True)
   
   velostep = array((array(mapmax - mapmin, copy=0).astype(Int32) - array(mapmax - mapmin, copy=0).astype(Int32) / 10.) / 10., copy=0).astype(Int32)
   if velostep < 1:   
      velostep = 1
   levels = (findgen(9) + 0.5) * velostep + mapmin
   loadct(0, silent=True)
   contour(mom0, xaxis, yaxis, levels=levels, overplot=True, c_colors=concatenate([0x000000]))
   contour(mom0mod, xaxis, yaxis, levels=levels, overplot=True, c_colors=concatenate([0xffffff]))
   beam_plot(beam[0], beam[1], bpa=bpa, center=concatenate([xaxis[0] - beam[0] / 3600., yaxis[0] + beam[0] / 3600.]), fill=True, transparent=True, color=0x000000)
   colormaps('sauron_colormap')
   colour_bar(concatenate([0.37, 0.39]), concatenate([0.1 + 0.2 * scrdim[0] / scrdim[1] + 0.02, 0.1 + 0.4 * scrdim[0] / scrdim[1] - 0.02]), strtrim(string(mapmin, format='(F10.1)'), 2), strtrim(string(mapmax, format='(F10.1)'), 2), opposite_label=True, black=True, title='(km s!E-1!N)', vertical=True, charthick=charthick, hex_color=True)
   loadct(0, silent=True)
   xyouts(0.45, 0.01 + 0.2 * scrdim[0] / scrdim[1] - 0.06, 'Velocity Field Contours start at ' + strtrim(string(levels[0], format='(F10.1)'), 2) + ' km s!E-1!N and increase with ' + strtrim(string(velostep, format='(I10)'), 2) + ' km s!E-1!N.', color=0x000000, normal=True, charthick=charthick)
   
   
   #PV Diagram along major axis
   if array(version, copy=0).astype(Int32) == version:   
      spawn('ls -1 PV-Diagrams/' + filenames[0] + '_2_xv.fits', mom0name)
   else:   
      spawn('ls -1 PV-Diagrams/' + filenames[0] + '_1_xv.fits', mom0name)
   mom0 = readfits(mom0name[array(mom0name, copy=0).nelements() - 1], mom0hed, silent=True)
   mom0mod = readfits('PV-Diagrams/Finalmodel_xv.fits', mom0hedmod, silent=True)
   velbuf = (2. * velext) / (ceninc * 0.2) + disper + 2. * sxpar(mom0hed, 'CDELT2')
   
   mapinmax = max(mom0[where(ravel(finite(mom0)))[0]], min=mapinmin)
   if abs(mapinmin / mapinmax) > 0.2:   
      mapinmin = -1 * mapinmax / 5
   buildaxii(mom0hed, xaxis, yaxis)
   if mapmax + velbuf > yaxis[array(yaxis, copy=0).nelements() - 1]:   
      tmp = yaxis[array(yaxis, copy=0).nelements() - 1]
   else:   
      tmp = mapmax + velbuf
   if mapmin - velbuf < yaxis[0]:   
      tmpmin = yaxis[0]
   else:   
      tmpmin = mapmin - velbuf
   yrange = concatenate([tmpmin, tmp])
   colormaps('heat', invert=True)
   showpixelsmap(xaxis * 3600., yaxis, mom0, position=concatenate([0.65, 0.1 + 0.2 * scrdim[0] / scrdim[1], 0.85, 0.1 + 0.4 * scrdim[0] / scrdim[1]]), xtitle='Offset (arcsec)', ytitle='Velocity (km s!E-1!N)', blank_value=0., range=concatenate([mapinmin, mapinmax]), noerase=True, charthick=charthick, thick=thick, black=True, yrange=yrange, hex_color=True)
   if array(noise, copy=0).nelements() < 1:   
      noise = stddev(mom0[0:11,0:11])
   levels = concatenate([1, 2, 4, 8, 16, 32, 64, 128]) * 1.5 * noise
   levelsneg = (concatenate([-2, -1])) * 1.5 * noise
   loadct(0, silent=True)
   contour(mom0, xaxis * 3600., yaxis, levels=levels, overplot=True, c_colors=concatenate([0x000000]))
   contour(mom0, xaxis * 3600, yaxis, levels=levelsneg, overplot=True, c_colors=concatenate([0x999999]), c_linestyle=2)
   loadct(40, silent=True)
   contour(mom0mod, xaxis * 3600, yaxis, levels=levels, overplot=True, c_colors=concatenate([0x0000ff]))
   colormaps('heat', invert=True)
   colour_bar(concatenate([0.87, 0.89]), concatenate([0.1 + 0.2 * scrdim[0] / scrdim[1] + 0.02, 0.1 + 0.4 * scrdim[0] / scrdim[1] - 0.02]), strtrim(string(mapinmin, format='(F10.4)'), 2), strtrim(string(mapinmax, format='(F10.4)'), 2), opposite_label=True, black=True, title='(Jy bm!E-1!N)', vertical=True, charthick=charthick, hex_color=True)
   loadct(0, silent=True)
   xyouts(0.45, 0.01 + 0.2 * scrdim[0] / scrdim[1] - 0.08, 'PV-Diagram Contours start are at -3, -1.5, 1.5, 3, 6, 12, 24 x rms.', color=0x000000, normal=True, charthick=charthick)
   xyouts(0.45, 0.01 + 0.2 * scrdim[0] / scrdim[1] - 0.1, 'rms = ' + strtrim(string(noise, format='(F10.5)'), 2) + ' Jy bm!E-1!N.', color=0x000000, normal=True, charthick=charthick)
   xyouts(0.45, 0.01 + 0.2 * scrdim[0] / scrdim[1] - 0.12, 'The distance used for conversions = ' + strtrim(string(distance, format='(F10.1)'), 2) + ' Mpc', color=0x000000, normal=True, charthick=charthick)
   if logical_not((gdlidl)):   
      image = tvrd(true=True)
   device(close=True)
   if logical_not((gdlidl)):   
      write_png('Overview.png', image)
   if gdlidl:   
      spawn("sed -i -- 's/1185 669/506 679/g' Overview.ps")
      spawn("sed -i -- 's/1186 669/506 679/g' Overview.ps")
      if file_test('Overview.ps--'):   
         spawn('rm -f Overview.ps--')
      spawn('gs -help', result)
      if array(result, copy=0).nelements() > 1:   
         spawn('gs -r300  -dEPSCrop -dTextAlphaBits=4 -sDEVICE=png16m -sOutputFile="Overview.png" -dBATCH -dNOPAUSE "Overview.ps"')
         spawn('rm -f Overview.ps')
   converted = 0
   spawn('convert --help', result, notfound)
   if array(result, copy=0).nelements() > 1:   
      gf = strsplit(result[0], ' ', extract=True)
      if array(gf, copy=0).nelements() > 1:   
         if strtrim(gf[1], 2) == 'ImageMagick':   
            converted = 1
            spawn('convert Overview.png -trim Overview.png')
   if converted < 1:   
      spawn('/usr/bin/convert --help', result, notfound)
      if array(result, copy=0).nelements() > 1:   
         gf = strsplit(result[0], ' ', extract=True)
         if array(gf, copy=0).nelements() > 1:   
            if strtrim(gf[1], 2) == 'ImageMagick':   
               converted = 1
               spawn('/usr/bin/convert Overview.png -trim Overview.png')
   if converted < 1:   
      spawn('imconvert --help', result, notfound)
      if array(result, copy=0).nelements() > 1:   
         gf = strsplit(result[0], ' ', extract=True)
         if array(gf, copy=0).nelements() > 1:   
            if strtrim(gf[1], 2) == 'ImageMagick':   
               converted = 1
               spawn('imconvert Overview.png -trim Overview.png')
   if converted < 1:   
      spawn('convert-im6 --help', result, notfound)
      if array(result, copy=0).nelements() > 1:   
         gf = strsplit(result[0], ' ', extract=True)
         if array(gf, copy=0).nelements() > 1:   
            if strtrim(gf[1], 2) == 'ImageMagick':   
               converted = 1
               spawn('convert-im6 Overview.png -trim Overview.png')
   
   
   
   return _ret()

"""
