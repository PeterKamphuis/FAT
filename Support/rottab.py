from numarray import *

def rottab(x, y, rotangle, xrotated=None, yrotated=None, center=None):
   """
    NAME:
          ROTTAB
   
    PURPOSE:
          Program to rotate tables consisting of x positions and ypositions;
   
    CATEGORY:
          Support
   
    CALLING SEQUENCE:
          ROTTAB,x,y,rotangle,xrotated=xrot,yrotated=yrot,center=center
   
    INPUTS:
            x = Array of x positions
            y = Array of y positions
     rotangle = Required angle of rotation
   
    OPTIONAL INPUTS:
           xrotated = array that contains the rotated values. If not
           present then x will contain the output.
           yrotated = array that contains the rotated values. If not
           present then y will contain the output
           center = Point of rotation, Default 0,0
   
    KEYWORD PARAMETERS:
   
    OUTPUTS:
   
    OPTIONAL OUTPUTS:
          -
   
    PROCEDURES CALLED:
         COS(), SIN()
   
    MODIFICATION HISTORY:
          Written 01-01-2009 P.Kamphuis v1.0
   
    NOTE:
   
   """

   n_params = 3
   xrot = xrotated
   yrot = yrotated
   _opt = (xrot, yrot, center)
   def _ret():
      _optrv = zip(_opt, [xrot, yrot, center])
      _rv = [x, y, rotangle]
      _rv += [_o[1] for _o in _optrv if _o[0] is not None]
      return tuple(_rv)
   
   # COMPILE_OPT IDL2
   if array(center, copy=0).nelements() >= 1:   
      tempx = dblarr(array(x, copy=0).nelements())
      tempx = dblarr(array(y, copy=0).nelements())
      tempx = x - center[0]
      tempy = y - center[1]
   else:   
      tempx = dblarr(array(x, copy=0).nelements())
      tempx = dblarr(array(y, copy=0).nelements())
      tempx = x
      tempy = y
   tempx2 = dblarr(array(x, copy=0).nelements())
   tempy2 = dblarr(array(y, copy=0).nelements())
   for i in arange(0, (array(tempx, copy=0).nelements() - 1)+(1)):
      tempx2[i] = (cos((rotangle * _sys_pi / 180.)) * tempx[i] - sin((rotangle * _sys_pi / 180.)) * tempy[i])
      tempy2[i] = (sin((rotangle * _sys_pi / 180.)) * tempx[i] + cos((rotangle * _sys_pi / 180.)) * tempy[i])
   if array(center, copy=0).nelements() >= 1:   
      tempx2[:] = tempx2[:] + center[0]
      tempy2[:] = tempy2[:] + center[1]
   if array(xrot, copy=0).nelements() < 1:   
      x = tempx2
   else:   
      xrot = dblarr(array(x, copy=0).nelements())
      xrot = tempx2
   if array(yrot, copy=0).nelements() < 1:   
      y = tempy2
   else:   
      yrot = dblarr(array(y, copy=0).nelements())
      yrot = tempy2
   
   return _ret()



