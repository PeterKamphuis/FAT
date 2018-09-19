import matplotlib.pyplot as plt
import matplotlib.colors as colors
#from mpl_toolkits.axes_grid1 import make_axes_locatable
from astropy.wcs import WCS
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from mpl_toolkits.axes_grid1.colorbar import colorbar
from columndensity import columndensity
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
import numpy as np
def plot_fits(filetoplot,figure,contours=None, figure_coordinates=[0.1,0.1,0.8,0.8],cmap='jet',log=False,aspect=None):
    hdu = fits.open(filetoplot)[0]
    plt.rc('xtick', color='k', labelsize='medium', direction='in')
    plt.rc('ytick', color='k', labelsize='medium', direction='in')
    space=figure_coordinates[2]*7.5
    wcsfound=False
    try:
        wcs = WCS(hdu.header)
        fitsplt = figure.add_axes(figure_coordinates,projection=wcs)
        wcsfound=True
    except:
        fitsplt = figure.add_axes(figure_coordinates)
        xaxis=[hdu.header['CRVAL1']+(i-hdu.header['CRPIX1']+1)*(hdu.header['CDELT1']) for i in range(hdu.header['NAXIS1'])]
        yaxis=[hdu.header['CRVAL2']+(i-hdu.header['CRPIX2']+1)*(hdu.header['CDELT2']) for i in range(hdu.header['NAXIS2'])]
         
    # we have to check for nan blanks and set them to 0
    tmp = np.isnan(hdu.data)
    nonzero= np.array(np.where(tmp == True))
#    if not log:
#        print(nonzero,nonzero.shape[1],len(nonzero))
#        exit()

    if nonzero.shape[1] < 1:
        nonzero= np.array(hdu.data[np.where(hdu.data != 0)])
    else:
        nonzero= np.array(hdu.data[np.where(tmp == False)])
  
    #tr_fk5 = fitsplt.get_transform("fk5")
    #    figure.subplots_adjust(left=figure_coordinates[0],bottom=figure_coordinates[1],top=figure_coordinates[3],right=figure_coordinates[2])
    #nonzero= np.array(hdu.data[np.where(hdu.data != 0)])
  
    mini = np.min(nonzero)
    counter=0
    while np.array(np.where(nonzero < mini)[0]).shape[0] < nonzero.shape[0]/10.:
        mini = mini+abs(mini/10.)
        counter += 1
        if counter > 100:
            break
    counter=0
    maxi = np.max(nonzero)
    while np.array(np.where(nonzero > maxi)[0]).shape[0] < nonzero.shape[0]/200.:
        maxi = maxi-abs(maxi/50.)
        counter += 1
        if counter >100:
            break

   
    if log:
      #  try:
        print(mini)
        if mini < 0:
            print('Why are we not adjusting?')
            rows,cols=np.where(hdu.data > 1e-3)
            newarray=hdu.data[rows,cols]
            mini = np.min(newarray)
            print(mini)
            counter=0
            while np.array(np.where(newarray < mini)[0]).shape[0] < newarray.shape[0]/10.:
                mini = mini+abs(mini/10.)
                counter += 1
                if counter > 100:
                    break
        else:
            newarray=hdu.data
        print(mini)
       
        i = fitsplt.imshow(hdu.data, origin='lower',cmap=cmap,norm=colors.LogNorm(vmin=mini, vmax=maxi),aspect=aspect)
      #  except ValueError:
       #     i = fitsplt.imshow(hdu.data, origin='lower',cmap=cmap,vmin=mini, vmax=maxi,aspect=aspect) 
    else:
        i = fitsplt.imshow(hdu.data, origin='lower',cmap=cmap,vmin=mini, vmax=maxi,aspect=aspect)
    try: 
        cbar = figure.colorbar(i)
    except:
        pass
    if wcsfound:
        ra = fitsplt.coords[0]
        ra.set_major_formatter('hh:mm:ss')
        ra.set_ticks(number=round(space))
    else:
        plt.gca().set_xticks(range(len(xaxis))[0:-1:int(len(xaxis)/5)])
        plt.gca().set_yticks(range(len(yaxis))[0:-1:int(len(yaxis)/5)])
        plt.gca().set_xticklabels(['{:10.0f}'.format(i*3600) for i in xaxis[0:-1:int(len(xaxis)/5)]])
        plt.gca().set_yticklabels(['{:10.1f}'.format(i) for i in yaxis[0:-1:int(len(yaxis)/5)]])
    fitsplt.set_xlabel(hdu.header['CTYPE1'])
    fitsplt.set_ylabel(hdu.header['CTYPE2'])
    if log:
        levels=[1E20, 2E20, 4E20, 8E20,16E20,32E20]
        try:
            beam=[hdu.header['BMAJ']*3600.,hdu.header['BMIN']*3600.]
        except:
            beam=[1.,1.]
        levels = columndensity(levels,vwidth=1,beam=beam,ncolumn=True)/1000.
    else:
        if mini > 0:
            maxvrot=(maxi-mini)/2.
            velostep=int(((int(maxi-mini)-int(maxi-mini)/10.)/10.)*1.5)
            levels=[(i)*velostep+mini for i in range(15)]
            #            levels=[(i*((maxi-mini)/15))+mini for i in range(9)]
        else:
            
            sig1=np.std(hdu.data[0:10,-10:-1])
            sig2=np.std(hdu.data[0:10,0:10])
            sig3=np.std(hdu.data[-10:-1,-10:-1])
            sig4=np.std(hdu.data[-10:-1,0:10])
            sig = np.mean([sig1,sig2,sig3,sig4])
            if np.isnan(sig):
                sig=abs(mini)/3.
            while sig > maxi:
                sig =sig/2.
            levels=np.array([-2,-1,1,2,4,8,16,32,64,128])*1.5*sig
            print(levels)
#            levels=[-sig*3,-sig*1.5,sig*1.5,3*sig,9*sig]
#            for i in range(3,15):
#                if ((i**2)*((maxi+mini)/30.)) > 9*sig:
#                    levels.append(((i**2)*((maxi+mini)/30.)))
    contourdata = fitsplt.contour(hdu.data, levels, colors='k', origin='lower',linewidths=0.75)
    if contours:
       
        hdumod = fits.open(contours)[0]
        if mini < 0:
            contour = fitsplt.contour(hdumod.data, levels, colors='r', origin='lower')
        else:
            try:
                contour = fitsplt.contour(hdumod.data, levels, colors='w', origin='lower')
            except UserWarning:
                contour = fitsplt.contour(hdumod.data, levels/2., colors='w', origin='lower')
