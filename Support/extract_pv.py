from numarray import *

def extract_pv(cube, header, pa, xv, center=None, xvheader=None):
   """
    NAME:
          EXTRACT_PV
   
    PURPOSE:
          Program to extract a PV diagram along a Position Angle from a
          Line emission cube. The strip width along which to extract the
          XV-diagram is the beam's major axis FWHM  if stated in
          the header else a single pixel is used.
   
    CATEGORY:
          Support
   
    CALLING SEQUENCE:
          EXTRACT_PV,Cube,header,pa,xv,CENTER=center,XVHEADER=new_header
   
   
    INPUTS:
          Cube = The array containing the Data Cube
          header = The header of the data cube
          pa = The position angle along which to extract this runs from
          north east wards
   
    OPTIONAL INPUTS:
          -
   
    KEYWORD PARAMETERS:
          CENTER = if set the rotation will happen around this center,
          if blank it will rotate around the center set in the header of
          the cube
          XVHEADER = A new header to write the xv array to a fits files
   
    OUTPUTS:
          xv = a 2-dimensional array with the XV - Diagram
   
    OPTIONAL OUTPUTS:
          -
   
    PROCEDURES CALLED:
          REVERSE(),ROT(),SXADDPAR,SXPAR(),SXDELPAR,TOTAL()
   
    EXAMPLE:
   
   
    MODIFICATION HISTORY:
          11-07-2018 P.Kamphuis; Replaced bitwise operator not with
                                 logical operator ~
          07-01-2016 P.Kamphuis; Replaced SUM commands with the proper
          TOTAL commands reducing the need for outside routines.
          05-01-2016 P.Kamphuis; Updated NAXIS1 in new header
          Written by P.Kamphuis 01-01-2015
   
    NOTE:
   
   """

   n_params = 4
   new_header = xvheader
   _opt = (center, new_header)
   def _ret():
      _optrv = zip(_opt, [center, new_header])
      _rv = [cube, header, pa, xv]
      _rv += [_o[1] for _o in _optrv if _o[0] is not None]
      return tuple(_rv)
   
   # COMPILE_OPT IDL2
   inheader = header
   if array(center, copy=0).nelements() == 0:   
      center = concatenate([sxpar(inheader, 'CRVAL1'), sxpar(inheader, 'CRVAL2')])
   if array(width, copy=0).nelements() == 0:   
      if logical_not((sxpar(inheader, 'BMAJ'))):   
         width = 0.
      else:   
         width = sxpar(inheader, 'BMAJ')
   xv = dblarr(array(cube[0,0,:], copy=0).nelements(), array(cube[:,0,0], copy=0).nelements())
   ypix = sxpar(inheader, 'CRPIX2') + (center[1] - sxpar(inheader, 'CRVAL2')) / sxpar(inheader, 'CDELT2')
   xpix = sxpar(inheader, 'CRPIX1') + (center[0] - sxpar(inheader, 'CRVAL1')) / sxpar(inheader, 'CDELT1') * cos(center[1] * _sys_dtor)
   xrange = concatenate([-xpix, sxpar(inheader, 'NAXIS1') - xpix])
   if xpix > sxpar(inheader, 'NAXIS1') - xpix - 1:   
      xsize = floor(sxpar(inheader, 'NAXIS1') - xpix - 1)
   else:   
      xsize = floor(xpix)
   xv = dblarr(array(2 * xsize, copy=0).astype(Int32), array(cube[:,0,0], copy=0).nelements())
   newcube = dblarr(array(cube[0,0,:], copy=0).nelements(), array(cube[0,:,0], copy=0).nelements(), array(cube[:,0,0], copy=0).nelements())
   for j in arange(0, (array(cube[:,0,0], copy=0).nelements() - 1)+(1)):
      newcube[j,:,:] = rot(cube[j,:,:], pa - 90., 1.0, xpix, ypix, missing=_sys_values.f_nan, cubic=-1, pivot=True) #For HI
   if width / (sxpar(inheader, 'CDELT2')) < 1:   
      xv[:,:] = newcube[:,round(ypix),array(xpix - xsize, copy=0).astype(Int32):(array(xpix + xsize - 1, copy=0).astype(Int32))+1]
   else:   
      centrpix = xpix
      istart = array(xpix - xsize, copy=0).astype(Int32)
      
      if istart < 0:   
         istart = 0
      if istart != 0:   
         centrpix = xpix - istart
      iend = array(xpix + xsize - 1, copy=0).astype(Int32)
      if iend - istart > array(xv[0,:], copy=0).nelements() - 1:   
         iend = array(xv[0,:], copy=0).nelements() - 1
      xv[:,0:(iend - istart)+1] = total(newcube[:,array(ypix - (width) / (2. * (sxpar(inheader, 'CDELT2'))), copy=0).astype(Int32):(array(ypix + (width) / (2. * (sxpar(inheader, 'CDELT2'))), copy=0).astype(Int32))+1,istart:(iend)+1], 2) / array(newcube[0,array(ypix - (width) / (2. * (sxpar(inheader, 'CDELT2'))), copy=0).astype(Int32):(array(ypix + (width) / (2. * (sxpar(inheader, 'CDELT2'))), copy=0).astype(Int32))+1,0], copy=0).nelements()
   if sxpar(inheader, 'CDELT1') < 0:   
      sxaddpar(inheader, 'CDELT1', abs(sxpar(inheader, 'CDELT1')))
      xv = reverse(xv)
   
   new_header = inheader
   sxaddpar(new_header, 'NAXIS', 2)
   sxaddpar(new_header, 'CDELT1', sxpar(inheader, 'CDELT1'))
   sxaddpar(new_header, 'CRVAL1', 0.)
   sxaddpar(new_header, 'CRPIX1', xsize + 1)
   sxaddpar(new_header, 'CTYPE1', 'ANGLE')
   sxaddpar(new_header, 'NAXIS1', array(2 * xsize, copy=0).astype(Int32))
   sxaddpar(new_header, 'CUNIT1', 'DEGREE', after='CRPIX1')
   sxaddpar(new_header, 'PA', pa, after='CTYPE2')
   sxaddpar(new_header, 'RA POS', center[0], after='PA')
   sxaddpar(new_header, 'DEC POS', center[1], after='RA POS')
   if width / (sxpar(inheader, 'CDELT2')) > 1:   
      sxaddpar(new_header, 'Strip Width', width, after='PA')
   sxaddpar(new_header, 'CDELT2', sxpar(inheader, 'CDELT3'))
   sxaddpar(new_header, 'CRVAL2', sxpar(inheader, 'CRVAL3'))
   sxaddpar(new_header, 'CRPIX2', sxpar(inheader, 'CRPIX3'))
   sxaddpar(new_header, 'CTYPE2', sxpar(inheader, 'CTYPE3'))
   sxaddpar(new_header, 'NAXIS2', sxpar(inheader, 'NAXIS3'))
   sxdelpar(new_header, concatenate(['NAXIS3', 'CDELT3', 'CRPIX3', 'CRVAL3', 'CTYPE3', 'LTYPE']))
   if sxpar(inheader, 'CUNIT3'):   
      sxaddpar(new_header, 'CUNIT2', sxpar(inheader, 'CUNIT3'), after='CTYPE2')
      sxdelpar(new_header, 'CUNIT3')
   if logical_not((sxpar(new_header, 'CUNIT2'))):   
      if sxpar(new_header, 'CDELT2') > 500.:   
         sxaddpar(new_header, 'CUNIT2', 'M/S', after='CTYPE2')
      else:   
         sxaddpar(new_header, 'CUNIT2', 'KM/S', after='CTYPE2')
   if strupcase(strtrim(sxpar(new_header, 'CUNIT2'), 2)) == 'M/S':   
      print linenumber() + 'EXTRACT_PV: We are converting to KM/S'
      sxaddpar(new_header, 'CDELT2', sxpar(new_header, 'CDELT2') / 1000.)
      sxaddpar(new_header, 'CRVAL2', sxpar(new_header, 'CRVAL2') / 1000.)
      sxaddpar(new_header, 'CUNIT2', 'KM/S')
   
   
   return _ret()







