from numarray import *

def int_profilev2(inputor, profile, pa=None, minor=None, header=None, xcenter=None, ycenter=None, range=None, axis=None, rotimage=None, noaverage=None):
   """
    NAME:
          INT_PROFILEV2
   
    PURPOSE:
          Program to make a radial or vertical profile from an image
   
    CATEGORY:
          Support
   
    CALLING SEQUENCE:
          INT_PROFILEV2,inputor,profile, PA = paor, /MINOR,
          HEADER = header, XCENTER = xcenter, YCENTER = ycenter,
          RANGE = range, AXIS = axis, ROTIMAGE = rotimage, /NOAVERAGE
   
   
    INPUTS:
          inputor = 2D array with the values of the image pixels
   
    OPTIONAL INPUTS:
          -
   
    KEYWORD PARAMETERS:
          /MINOR - extract the minor axis profile instead of the major axis
          PA = position angle of the galaxy be aware that that PA is
          defined from North so a galaxy aligned on the y-axis has a PA of 90
          HEADER = Image's header can be given to rotate the image
          around the images center or to give range in map units.
          XCENTER= rotation center on xaxis in pixels. Supersedes the
          header. Default is half the xaxis
          YCENTER= rotation center on yaxis in pixels. Supersedes the
          header default is half the yaxis
          RANGE = the range over which to collapse the map/cube without
          a header in pixels when header is given it assumed to be in
          map units if not given the whole axis will be collapsed, when
          2x2 array is assumed to be [xrange,yrange]
          AXIS = can be set to return a xaxis for your profile, if
          header is given it will be based on header otherwise it is pixels offset from the center.
          /NOAVERAGE - Set to not get an average profile but an
                       integrated one
          ROTIMAGE = If set the rotated image will be returned in this array
   
    OUTPUTS:
          profile = the extracted profile
   
    OPTIONAL OUTPUTS:
          -
   
    PROCEDURES CALLED:
          ROT(), SIZE(), STRUPCASE(), STRTRIM() STRCOMPRESS(),
          STR_SEP(), SXPAR(), TOTAL()
   
    EXAMPLE:
   
   
    MODIFICATION HISTORY:
          11-07-2018 P. Kamphuis; Added some check condintion to make
                                  sure the arrays can be created and the
                                  routine does not crash when range is
                                  smaller than a pixel.
          Written by P.Kamphuis 01-01-2015
   
    NOTE:
   
   """

   n_params = 2
   paor = pa
   _opt = (paor, minor, header, xcenter, ycenter, range, axis, rotimage, noaverage)
   def _ret():
      _optrv = zip(_opt, [paor, minor, header, xcenter, ycenter, range, axis, rotimage, noaverage])
      _rv = [inputor, profile]
      _rv += [_o[1] for _o in _optrv if _o[0] is not None]
      return tuple(_rv)
   
   # COMPILE_OPT IDL2
   input = inputor
   if array(paor, copy=0).nelements() > 0:   
      pa = paor
   else:   
      pa = 0
   arsize = size(input)
   if arsize[0] > 2:   
      moments(input, inputmap)
      input = dblarr(arsize[1], arsize[2])
      input = inputmap
      arsize = size(input)
   if (minor is not None):   
      pa = pa + 90
   
   if array(xcenter, copy=0).nelements() == 0:   
      if array(header, copy=0).nelements() == 0:   
         xcenter = array(arsize[1] / 2., copy=0).astype(Int32)
      else:   
         xcenter = sxpar(header, 'CRPIX1')
   if array(ycenter, copy=0).nelements() == 0:   
      if array(header, copy=0).nelements() == 0:   
         ycenter = array(arsize[2] / 2., copy=0).astype(Int32)
      else:   
         ycenter = sxpar(header, 'CRPIX2')
   
   if array(header, copy=0).nelements() != 0:   
      xcenterdeg = (xcenter - sxpar(header, 'CRPIX1')) * sxpar(header, 'CDELT1') + sxpar(header, 'CRVAL1')
      ycenterdeg = (ycenter - sxpar(header, 'CRPIX2')) * sxpar(header, 'CDELT2') + sxpar(header, 'CRVAL2')
   
   if pa != 90:   
      result = rot(input, pa - 90, 1.0, xcenter, ycenter, cubic=-0.5, pivot=True)
      input = result
   
   
   if array(header, copy=0).nelements() > 0:   
      buildaxii(header, xaxis, yaxis)
   pixpos = dblarr(2, 2)
   if array(range, copy=0).nelements() > 0:   
      rangesize = size(range)
      if array(header, copy=0).nelements() > 0:   
         if rangesize[0] == 1:   
            
            if (minor is not None):   
               pixpos[0,0] = 0
               pixpos[0,1] = array(input[:,0], copy=0).nelements() - 1
               if xaxis[0] - xaxis[1] < 0:   
                  tmp = where(ravel(xaxis > range[0] + xcenterdeg))[0]
                  pixpos[1,0] = tmp[0]
                  tmp = where(ravel(xaxis < range[1] + xcenterdeg))[0]
                  pixpos[1,1] = tmp[array(tmp, copy=0).nelements() - 1]
               else:   
                  tmp = where(ravel(xaxis < range[0] + xcenterdeg))[0]
                  pixpos[1,1] = tmp[0]
                  tmp = where(ravel(xaxis > range[1] + xcenterdeg))[0]
                  pixpos[1,0] = tmp[array(tmp, copy=0).nelements() - 1]
               
               
            else:   
               pixpos[0,0] = 0
               pixpos[0,1] = array(input[1,:], copy=0).nelements() - 1
               tmp = where(ravel(yaxis > range[0] + ycenterdeg))[0]
               pixpos[1,0] = tmp[0]
               tmp = where(ravel(yaxis < range[1] + ycenterdeg))[0]
               
               pixpos[1,1] = tmp[array(tmp, copy=0).nelements() - 1]
               
            
         if rangesize[0] == 2:   
            
            if (minor is not None):   
               if xaxis[0] - xaxis[1] < 0:   
                  tmp = where(ravel(xaxis > range[0,0] + xcenterdeg))[0]
                  pixpos[1,0] = tmp[0]
                  tmp = where(ravel(xaxis < range[0,1] + xcenterdeg))[0]
                  pixpos[1,1] = tmp[array(tmp, copy=0).nelements() - 1]
               else:   
                  tmp = where(ravel(xaxis < range[0,0] + xcenterdeg))[0]
                  pixpos[1,1] = tmp[0]
                  tmp = where(ravel(xaxis > range[0,1] + xcenterdeg))[0]
                  pixpos[1,0] = tmp[array(tmp, copy=0).nelements() - 1]
               tmp = where(ravel(yaxis > range[1,0] + ycenterdeg))[0]
               pixpos[0,0] = tmp[0]
               tmp = where(ravel(yaxis < range[1,1] + ycenterdeg))[0]
               pixpos[0,1] = tmp[array(tmp, copy=0).nelements() - 1]
            else:   
               if xaxis[0] - xaxis[1] < 0:   
                  tmp = where(ravel(xaxis > range[0,0] + xcenterdeg))[0]
                  pixpos[0,0] = tmp[0]
                  tmp = where(ravel(xaxis < range[0,1] + xcenterdeg))[0]
                  pixpos[0,1] = tmp[array(tmp, copy=0).nelements() - 1]
               else:   
                  tmp = where(ravel(xaxis < range[0,0] + xcenterdeg))[0]
                  pixpos[0,1] = tmp[0]
                  tmp = where(ravel(xaxis > range[0,1] + xcenterdeg))[0]
                  pixpos[0,0] = tmp[array(tmp, copy=0).nelements() - 1]
               tmp = where(ravel(yaxis > range[1,0] + ycenterdeg))[0]
               pixpos[1,0] = tmp[0]
               tmp = where(ravel(yaxis < range[1,1] + ycenterdeg))[0]
               pixpos[1,1] = tmp[array(tmp, copy=0).nelements() - 1]
            
      else:   
         if rangesize[0] == 1:   
            pixpos[0,0] = 0
            pixpos[0,1] = array(input[1,:], copy=0).nelements() - 1
            pixpos[1,0] = range[0]
            pixpos[1,1] = range[1]
         if rangesize[0] == 2:   
            if (minor is not None):   
               pixpos[1,0] = range[0,0]
               pixpos[1,1] = range[0,1]
               pixpos[0,0] = range[1,0]
               pixpos[0,1] = range[1,1]
            else:   
               pixpos[0,0] = range[0,0]
               pixpos[0,1] = range[0,1]
               pixpos[1,0] = range[1,0]
               pixpos[1,1] = range[1,1]
   else:   
      
      pixpos[:,0] = 0
      pixpos[0,1] = array(input[1,:], copy=0).nelements() - 1
      pixpos[1,1] = array(input[:,1], copy=0).nelements() - 1
   if pixpos[0,0] > pixpos[0,1]:   
      pixpos[0,:] = reverse(pixpos[0,:])
   if pixpos[1,0] > pixpos[1,1]:   
      pixpos[1,:] = reverse(pixpos[1,:])
   if array(pixpos[1,1] - pixpos[1,0], copy=0).astype(Int32) == 0:   
      pixpos[1,0] = pixpos[1,0] - 1
   if array(pixpos[0,1] + 1 - pixpos[0,0], copy=0).astype(Int32) == 0:   
      pixpos[0,0] = pixpos[0,0] - 1
   
   profile = dblarr(pixpos[0,1] + 1 - pixpos[0,0])
   clear = dblarr(array(input[0,:], copy=0).nelements(), array(input[:,0], copy=0).nelements())
   tmp2 = where(ravel(finite(input)))[0]
   clear[tmp2] = input[tmp2]
   tmp = dblarr(array(pixpos[1,1] - pixpos[1,0], copy=0).astype(Int32))
   
   for i in arange(pixpos[0,0], (pixpos[0,1])+(1)):
   
      tmp = clear[pixpos[1,0]:(pixpos[1,1])+1,i]
      tmp2 = where(ravel(tmp == 0))[0]
      if tmp2[0] == -1:   
         if (noaverage is not None):   
            profile[i - pixpos[0,0]] = total(clear[pixpos[1,0]:(pixpos[1,1])+1,i]) / (array(clear[pixpos[1,0]:(pixpos[1,1])+1,i], copy=0).nelements())
         else:   
            profile[i - pixpos[0,0]] = total(clear[pixpos[1,0]:(pixpos[1,1])+1,i]) / (array(clear[pixpos[1,0]:(pixpos[1,1])+1,i], copy=0).nelements())
      else:   
         if (noaverage is not None):   
            profile[i - pixpos[0,0]] = total(clear[pixpos[1,0]:(pixpos[1,1])+1,i])
         else:   
            profile[i - pixpos[0,0]] = total(clear[pixpos[1,0]:(pixpos[1,1])+1,i]) / (array(clear[pixpos[1,0]:(pixpos[1,1])+1,i], copy=0).nelements() - array(tmp2, copy=0).nelements())
   tmpindex = where(ravel(finite(profile) == 0.))[0]
   if tmpindex[0] != -1:   
      profile[where(ravel(finite(profile) == 0.))[0]] = 0.
   axis = dblarr(array(profile, copy=0).nelements())
   if array(header, copy=0).nelements() > 0:   
      if (minor is not None):   
         axis[:] = yaxis[pixpos[0,0]:(pixpos[0,1])+1] - ycenterdeg
      else:   
         axis[:] = xaxis[pixpos[0,0]:(pixpos[0,1])+1] - xcenterdeg
   else:   
      if (minor is not None):   
         axis = findgen(array(axis, copy=0).nelements()) + pixpos[0,0] - ycenter
      else:   
         axis = findgen(array(axis, copy=0).nelements()) + pixpos[0,0] - xcenter
   if array(rotimage, copy=0).nelements() > 0:   
      rotimage = dblarr(array(input[0,:], copy=0).nelements(), array(input[:,0], copy=0).nelements())
      rotimage = input
   
   
   return _ret()



