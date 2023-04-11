************************************************************************
**                                                                    **
**                   Fit a correlation function peak                  **
**                                                                    **
**                                                                    **
**                  Copyright (C) 1999 John L. Tonry                  **
**                                                                    **
************************************************************************

c Included in this file:
c
c subroutine peakfit
c function poly
c subroutine quarticpeak
c subroutine parabolapeak



************************************************************************
* subroutine peakfit -- Fit the correlation peak
************************************************************************

      subroutine peakfit(npt,lz1,lz2,cfn,center,height,width,a,error)
* Fit a quartic to a range between inflection points plus 20% or so
* unless there are too few points, in which case fall back to a parabola
      integer i
      parameter (MAXFIT=512)
      real cfn(0:npt-1)
      real scale(2)
      real abc(0:MAXFIT-1), ord(0:MAXFIT-1)
      real*8 a(5), aa(5), rms, ave, v, vl, vr, h, w, poly
      logical error, noroot
      parameter (EXTRA=0.2, ERRMAX=0.01)
C      poly4(z) = a(1) + z*(a(2) + z*(a(3) + z*(a(4) + z*a(5))))
C      poly2(z) = a(1) + z*(a(2) + z*a(3))
      mmod(i) = mod(i+npt, npt)

* Initialize variables
      error = .false.
      noroot = .false.
C Change to cmax=0 to avoid problems for z=0 correlations
C      cmax = cfn(1)
      cmax = 0
      localmax = 0

* First find the peak
      iz2 = npt
      if(lz1.lt.0) iz2 = npt + lz1
C      do 900 i=1, npt-2
      do 900 i=0, npt-1
         if(i.lt.lz1.or.(i.gt.lz2.and.i.lt.iz2)) goto 900
         if(cfn(i).gt.cmax) then
            cmax = cfn(i)
            imax = i
            localmax = 0
            if(i.eq.0) then
               if(cfn(i).ge.cfn(npt-1).and.cfn(i).ge.cfn(1)) localmax=1
            else if(i.eq.npt-1) then
               if(cfn(i).ge.cfn(i-1).and.cfn(i).ge.cfn(0))   localmax=1
            else
               if(cfn(i).ge.cfn(i-1).and.cfn(i).ge.cfn(i+1)) localmax=1
            end if
         end if
 900  continue
C      WRITE(6,*) 'PEAK FOUND AT', IMAX, CMAX

* Check for no maximum
      if(cmax.eq.0) then
         write(6,*) 'PEAKFIT: Correlation function is all zero!'
         stop
      end if

* No peak in z range => return immediately
      if(localmax.eq.0) then
         center = imax
         width = 0
         height = cfn(imax)
         return
      end if

* Find the inflection points
      nodat = 5*npt
      nsearch = npt / 4
      il = nodat
      dcl = 0
      ir = nodat
      dcr = 0
      do i=1, nsearch
         if(il.ne.nodat.and.ir.ne.nodat) goto 9011
C     WRITE(6,*) I, IL, IR, MMOD(IMAX+I), MMOD(IMAX-I)
         dc = cfn(mmod(imax+i-1)) - cfn(mmod(imax+i))
C     WRITE(6,*) cfn(mmod(imax+i)), cfn(mmod(imax+i-1)), dc
         if(dc.lt.dcr.and.ir.eq.nodat) then
            ir = imax + i - 1
         end if
         dcr = dc
         dc = cfn(mmod(imax-i+1)) - cfn(mmod(imax-i))
C     WRITE(6,*) cfn(mmod(imax-i)), cfn(mmod(imax-i+1)), dc
         if(dc.lt.dcl.and.il.eq.nodat) then
            il = imax - i + 1
         end if
         dcl = dc
C     WRITE(6,*) I, IL, IR, MMOD(IMAX+I), MMOD(IMAX-I)
      end do
 9011 continue

      if(ir-il.lt.6) then
C     WRITE(6,*) 'PEAKFIT: narrow peak', IMAX, IL, IR
         il = imax - 3
         ir = imax + 3
      end if
      ifit1 = il - EXTRA*(ir-il)
      ifit2 = ir + EXTRA*(ir-il)
C     WRITE(6,*) IFIT1, IL, IR, IFIT2

* Find the half peak points
      ilh = nodat
      irh = nodat
      do i=1, nsearch
         if(ilh.ne.nodat.and.irh.ne.nodat) goto 9021
         if(irh.eq.nodat.and.cfn(mmod(imax+i)).lt.cmax/2) then
            irh = imax + i
         end if
         if(ilh.eq.nodat.and.cfn(mmod(imax-i)).lt.cmax/2) then
            ilh = imax - i
         end if
      end do
 9021 continue
      ifit1 = min(ifit1,ilh)
      ifit2 = max(ifit2,irh)
C     WRITE(6,*) IFIT1, ILH, IRH, IFIT2

* Find the zero crossing points
      ilz = nodat
      irz = nodat
      do i=1, nsearch
         if(ilz.ne.nodat.and.irz.ne.nodat) goto 9031
         if(irz.eq.nodat.and.cfn(mmod(imax+i)).lt.0) then
            irz = imax + i
         end if
         if(ilz.eq.nodat.and.cfn(mmod(imax-i)).lt.0) then
            ilz = imax - i
         end if
      end do
 9031 continue
      ifit1 = max(ifit1,ilz)
      ifit2 = min(ifit2,irz)
C     WRITE(6,*) IFIT1, ILZ, IRZ, IFIT2

* Fit a quartic to the data
* The fit must not extend past the zero crossing, 
* nor quit before the half peak points
      lastf1 = 0
      lastf2 = 0
 1    nfit = ifit2 - ifit1 + 1
      if(ifit1.gt.ilh.or.ifit2.lt.irh) then
         write(6,*) 'PEAKFIT: fit quits before half peak points!'
         write(6,*) 'How did this happen? Center =', imax
         write(6,*) 'left: fit, half =',ifit1, ilh, 
     $        'right: fit, half =',ifit2, irh
         write(6,*) 'left: inflect, zero =',il, ilz, 
     $        'right: inflect, zero =',ir, irz
         write(6,'(12i7)') (i,i=imax-5,imax+5)
         write(6,'(12f7.3)') (cfn(i)/cmax,i=imax-5,imax+5)
         stop
      end if
      if(nfit.ge.8) then
         ncoeff = 5
      else
         ncoeff = 3
C     write(6,*) 'PEAKFIT: too few points to fit well'
C     write(6,'(12i7)') (i,i=imax-5,imax+5)
C     write(6,'(12f7.3)') (cfn(mmod(i))/cmax,i=imax-5,imax+5)
      end if
      if(nfit.gt.MAXFIT) then
         write(6,*) 'PEAKFIT: nfit > MAXFIT!'
         stop
      end if
      
      do i=0, nfit-1
         abc(i) = ifit1 + i
         ord(i) = cfn(mmod(ifit1+i))
C     WRITE(6,'(2I5,2F8.2)') I, IFIT1+I, ABC(I), ORD(I)
      end do
      
      scale(1) = ifit1
      scale(2) = ifit2
      call fitlpoly(nfit,abc,ord,scale,ncoeff,aa)
C     WRITE(6,'(A,1P5G12.4)') ' fit',(AA(I)/CMAX,I=1,5)
C     WRITE(6,*) IFIT1, IFIT2, NFIT
      call polyc(ncoeff,scale,aa,a)
      
* Check that the residuals are small enough
      rms = 0
      ave = 0
      do i=ifit1, ifit2
         ave = ave + (cfn(mmod(i)) - poly(abc(i-ifit1),ncoeff,a))
         rms = rms + (cfn(mmod(i)) - poly(abc(i-ifit1),ncoeff,a)) *
     $        (cfn(mmod(i)) - poly(abc(i-ifit1),ncoeff,a))
C     WRITE(6,'(2I5,3F8.2)') I, I-IFIT1, CFN(MMOD(I)), ABC(I-IFIT1),
C     $        POLY(ABC(I-IFIT1),NCOEFF,A)
      end do
      ave = ave / nfit
      rms = dsqrt(rms/nfit-ave*ave)
      if(rms/cmax.gt.ERRMAX) then
C     WRITE(6,'(1X,A,F8.4,A,I5)')
C     1	'PEAKFIT: large residuals, rms =', rms/cmax,'  N =',nfit
         if(noroot.or.(ifit1.eq.ilh.and.ifit2.eq.irh)) goto 10
         if(nfit.gt.16) then
            ifit1 = min(ifit1+2,ilh)
            ifit2 = max(ifit2-2,irh)
         else
            ifit1 = min(ifit1+1,ilh)
            ifit2 = max(ifit2-1,irh)
         end if
         goto 1
      end if

* Compute the center and width of the polynomial
 10   v = imax
      vl = ilh
      vr = irh
      if(ncoeff.eq.5) then
         call quarticpeak(a,v,vl,vr,h,w,noroot)
      else
         a(4) = 0d0
         a(5) = 0d0
         call parabolapeak(a,v,vl,vr,h,w,noroot)
      end if
C     WRITE(6,*) 'ROOTS: ',VL,VR,NOROOT
      if(lastf1.eq.ifit1.and.lastf2.eq.ifit2) then
C     write(6,*) 'PEAKFIT: Can''t seem to better', ifit1,ifit2
C     write(6,'(1x,6i8)') ilz,il,ilh,imax,irh,ir,irz
C     write(6,'(a,1p5g12.4)') ' FIT',(aa(i)/cmax,i=1,5)
         goto 999
      end if
      lastf1 = ifit1
      lastf2 = ifit2
      if(noroot) then
         ifit1 = max(ifit1-1,ilz)
         ifit2 = min(ifit2+1,irz)
         goto 1
      end if
      center = v
      width = w
      height = poly(center,ncoeff,a)
 
C     WRITE(6,*) V, W
      
 999  return
      end



************************************************************************
* function poly -- Evaluate polynomial at z
************************************************************************

      function poly(z,n,a)
      real*8 a(1), poly
      poly = 0
      do i=n, 1, -1
         poly = a(i) + z*poly 
      end do
      return
      end



************************************************************************
* subroutine quarticpeak -- Compute the center/half-max points of a quartic
************************************************************************

      subroutine quarticpeak(a,x0,xl,xr,h,w,noroot)
* x0 and xl, xr are provided with initial guesses

      implicit real*8 (a-h,o-z)
      real*8 a(5)
      logical noroot

      p(z) = a(1) + z*(a(2) + z*(a(3) + z*(a(4) + z*a(5))))
      d(z) = a(2) + z*(2*a(3) + z*(3*a(4) + z*4*a(5)))
      c(z) = 2*a(3) + z*(6*a(4) + z*12*a(5))
      acc = 1d-3

C      WRITE(6,*) X0, XL, XR

      noroot = .false.
* Get the center with Newton
      do i=1, 50
         x = x0 - d(x0)/c(x0)
         if(abs(x-x0).lt.acc) goto 1
         x0 = x
      end do
1     x0 = .5*(x+x0)
      h = p(x0)
* Now test to see whether the polynomial ever reaches half max
* and if it does compute where it is (fudge a0 and look for roots)
      a(1) = a(1) - 0.5d0*h
      do iroot=1, 2
         if(iroot.eq.1) then
            x1 = xl
         else
            x1 = xr
         end if
C	xguess = x1
         step = .1
C         do s=step, 3, step
         do i=1, 30
            s = real(i)*step
            if(p(x0+s*(x1-x0)).lt.0) then
               x = (s-step/2)*(x1-x0) + x0
               goto 3
            end if
         end do
* No sign change and we've probably got no root
C      write(6,*) 'QUARTICPEAK: poly never reaches half maximum'
         noroot = .true.
         x = x1
C     xinit = x
         goto 4
* Get the zero with Newton's method
 3       continue
C	xinit = x
         do i=1, 50
            x1 = x - p(x)/d(x)
            if(abs(x1-x).lt.acc) goto 2
            x = x1
         end do
 2       x = .5*(x1+x)
 4       if(iroot.eq.1) then
            xl = x
         else
            xr = x
         end if
C     write(6,3000) xguess,xinit,x,p(x)/h+0.5,p(xinit)/h+.05
C3000	format(1x,'QUARTIC: V0, V1, V, p(V), p(V1) =',3f10.1,2f8.4)
      end do
      a(1) = a(1) + 0.5d0*h
      w = xr - xl
      return
      end



************************************************************************
* subroutine parabolapeak -- Compute the center/half-max points of a parabola
************************************************************************

      subroutine parabolapeak(a,x0,xl,xr,h,w,noroot)
* Center x0, height h, width w, left & right half maxima xl, xr
      implicit real*8 (a-h,o-z)
      real*8 a(5)
      logical noroot
      noroot = .false.
      if(a(3).eq.0) then
          x0 = 0
          h = 1
          w = 1
          noroot = .true.
          return
      end if
      x0 = -a(2) / (2*a(3))
      h = a(1) - a(2)*a(2) / (4*a(3))
      w = -2*h / a(3)
      if(w.ge.0) then
          w = sqrt(w)
      else
          w = -sqrt(-w)
      end if
      xl = x0 - 0.5*w
      xr = x0 + 0.5*w
      return
      end

