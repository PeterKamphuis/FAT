from numarray import *

def momentsv2(cube, momentmap, header, map, blank_value=None, gdlidl=None):
   """
    NAME:
          MOMENTSV2
   
    PURPOSE:
          Program to calculate the moment zero map of a 3D fits cube  and update the header.
   
    CATEGORY:
          Support
   
    CALLING SEQUENCE:
          momentsv2,Cube,Momentmap,header,map
   
   
    INPUTS:
          Cube = A 3D-array containing the values of the cube
          header = The fits header of the data cube
          map = 1 or 0 indicating the requested moment
   
    OPTIONAL INPUTS:
          gdlidl = indicator of running gdl or idl. Is 1 when running GDL
   
    KEYWORD PARAMETERS:
          -
   
    OUTPUTS:
          momentmap = the calaculated moment map
          header = the header will be updated to reflect the reduction
          in dimension and proper units.
   
    OPTIONAL OUTPUTS:
          -
   
    PROCEDURES CALLED:
          STR_SEP(), STRTRIM(), STRCOMPRESS(), TOTAL(), SXADDPAR(),
          SXPAR(), SXDELPAR()
   
    EXAMPLE:
   
   
    MODIFICATION HISTORY:
          20-08-2017 P.Kamphuis; Rebin is broken on the mac 0.9.7 hence
                                 we have added a condition to rebin in
                                 loop when running GDL. This is only
                                 used for velocity fields so no need to
                                 add to intergrated moment map
          12-05-2017 P.Kamphuis; In the  maps the blanks should
                                 propagate.
          16-11-2016 P.Kamphuis; Now dealing with the (non-)presence of
                                 CUNIT3 properly
          01-06-2016 P.Kamphuis; Added a condition to check that datamax
                                 and datamin are finite.
          07-01-2016 P.Kamphuis; Replaced SUM commands with the proper
          TOTAL commands reducing the need for outside routines.
          Modified to deal with existing CUNIT3 without error message 12-08-2015 P. Kamphuis v2
          Modified to deal with missing CUNIT3  29-07-2015 P. Kamphuis v2
          Modified to use SUM which is much faster 25-05-2015 P. Kamphuis v2
          Written 15-07-2010 P.Kamphuis v1.0
   
    NOTE:
   
   """

   n_params = 4
   blanked = blank_value
   _opt = (blanked, gdlidl)
   def _ret():
      _optrv = zip(_opt, [blanked, gdlidl])
      _rv = [cube, momentmap, header, map]
      _rv += [_o[1] for _o in _optrv if _o[0] is not None]
      return tuple(_rv)
   
   if map == 0:   
      #In the case of moment 0 we want blanks to propagate through.
      #   blank=WHERE(FINITE(Cube) NE 1.)
      #   IF blank[0] NE -1 then Cube[blank]=0
      momentmap = zeros([array(cube[0,:,0], copy=0).nelements(), array(cube[0,0,:], copy=0).nelements()], Float32)
      momentmap[:,:] = total(cube, 3) * abs(sxpar(header, 'CDELT3'))
      if logical_not((sxpar(header, 'CUNIT3'))):   
         if sxpar(header, 'CDELT3') > 500.:   
            sxaddpar(header, 'CUNIT3', 'M/S')
         else:   
            sxaddpar(header, 'CUNIT3', 'KM/S')
      if strupcase(strtrim(sxpar(header, 'CUNIT3'), 2)) == 'M/S':   
         sxaddpar(header, 'CUNIT3', 'KM/S')
         momentmap = momentmap / 1000.
         print linenumber() + 'MOMENTSV2: We have converted the units to Jy/Beam x Km/s'
      sxaddpar(header, 'BUNIT', strtrim(strcompress(sxpar(header, 'BUNIT')), 2) + '.' + strtrim(strcompress(sxpar(header, 'CUNIT3')), 2))
      maxmap = max(momentmap, min=minmap)
      if finite(maxmap):   
         sxaddpar(header, 'DATAMAX', maxmap)
      else:   
         sxaddpar(header, 'DATAMAX', max(cube))
      if finite(minmap):   
         sxaddpar(header, 'DATAMIN', minmap)
      else:   
         sxaddpar(header, 'DATAMIN', 0.)
      
   if map == 1:   
      buildaxii(header, xaxis, yaxis, zaxis=zaxis)
      #   blank=WHERE(FINITE(Cube) NE 1.)
      #   IF blank[0] NE -1 then Cube[blank]=0
      if logical_not((sxpar(header, 'CUNIT3'))):   
         if sxpar(header, 'CDELT3') > 500.:   
            sxaddpar(header, 'CUNIT3', 'M/S')
         else:   
            sxaddpar(header, 'CUNIT3', 'KM/S')
      if strupcase(strtrim(sxpar(header, 'CUNIT3'), 2)) == 'M/S':   
         print linenumber() + 'MOMENTSV2: We are converting to KM/S'
         zaxis = zaxis / 1000.
      momentmap = zeros([array(cube[0,:,0], copy=0).nelements(), array(cube[0,0,:], copy=0).nelements()], Float32)
      
      if gdlidl:   
         #this is necessary for gdl as rebin is
         #broken in gdl on the mac. This seems
         #an issue of 0.9.7. We should check on ubuntu and
         #issu a ticket or some such.
         c = dblarr(array(cube[0,0,:], copy=0).nelements(), array(cube[0,:,0], copy=0).nelements(), array(cube[:,0,0], copy=0).nelements())
         for i in arange(0, (array(cube[:,0,0], copy=0).nelements() - 1)+(1)):
            c[i,:,:] = zaxis[i]
      else:   
         c = rebin(reform(zaxis, 1, 1, array(zaxis, copy=0).nelements()), array(cube[0,0,:], copy=0).nelements(), array(cube[0,:,0], copy=0).nelements(), array(cube[:,0,0], copy=0).nelements())
      momentmap = total(c * cube, 3) / total(cube, 3)
      maxmap = max(momentmap, min=minmap)
      if finite(maxmap):   
         sxaddpar(header, 'DATAMAX', maxmap)
      else:   
         sxaddpar(header, 'DATAMAX', zaxis[array(zaxis, copy=0).nelements() - 1])
      if finite(minmap):   
         sxaddpar(header, 'DATAMIN', minmap)
      else:   
         sxaddpar(header, 'DATAMIN', zaxis[0])
      sxaddpar(header, 'BUNIT', 'KM/S')
      c = 0
   if array(blanked, copy=0).nelements() != 0:   
      tmp = where(ravel(momentmap == 0.))[0]
      if tmp[0] != -1:   
         momentmap[tmp] = blanked
   sxdelpar(header, 'CUNIT3')
   sxdelpar(header, 'CTYPE3')
   sxdelpar(header, 'CRVAL3')
   sxdelpar(header, 'CDELT3')
   sxdelpar(header, 'NAXIS3')
   sxdelpar(header, 'CRPIX3')
   sxaddpar(header, 'NAXIS', 2)
   
   
   return _ret()



