************************************************************************
**                                                                    **
**                   Prepare a spectrum for an FFT                    **
**                                                                    **
**                                                                    **
**                 Copyright (C) 1999 John L. Tonry                   **
**                                                                    **
************************************************************************

c Included in this file:
c
c subroutine apodize
c subroutine rebin
c subroutine zerospec
c subroutine addspec
c subroutine meanzero
c subroutine splinedex
c subroutine overlap
c function xlog
c subroutine despace
c subroutine medfilt
c subroutine medwfilt



************************************************************************
* subroutine apodize -- Apodize the ends of a spectrum with a cosine bell
************************************************************************

      subroutine apodize(n,n1,n2,data,percent)
* The cosine bell starts rising from zero at n1 and falls to zero at n2
      real data(1)
      parameter (PI=3.14159265)

      nsquash = min(n*0.01*percent,real((n2-n1)/2))
      if(nsquash.lt.1) return
      do i=0, nsquash-1
         arg = PI * float(i)/(nsquash-1)
         factor = .5 * (1 - cos(arg))
         data(n1+i) = factor * data(n1+i)
         data(n2-i) = factor * data(n2-i)
      end do
      return
      end



************************************************************************
* subroutine rebin --  Bin a spectrum in log wavelength
************************************************************************

      subroutine rebin(nwave,wave,fsrc,nlog,w0,dwlog,fdest)

      real wave(nwave),fsrc(nwave),fdest(nlog)

      do i=1, nlog
         fdest(i) = 0
      end do

* Rebin each source pixel onto the destination array
      do l=1, nwave
* Give each source pixel boundaries half way between tabulated wavelengths
         if(l.eq.1) then
            s0 = 0.5*(3*wave(l)-wave(l+1))
            s1 = 0.5*(wave(l)+wave(l+1))
         else if(l.eq.nwave) then
            s0 = 0.5*(wave(l-1)+wave(l))
            s1 = 0.5*(3*wave(l)-wave(l-1))
         else
            s0 = 0.5*(wave(l-1)+wave(l))
            s1 = 0.5*(wave(l)+wave(l+1))
         end if
* What do these boundaries map to in log wavelength array space
         s0log = alog(s0/w0)/dwlog + 1
         s1log = alog(s1/w0)/dwlog + 1
* I *think* this is Fnu being tabulated
C         dnu = (s1-s0)*2.99793e10/(wave(l)*wave(l)*1e-8)
* Nope, Saurabh says Flambda
         dnu = (s1-s0)

* Run over rebinning loop
         do 10 i=int(s0log), int(s1log)
            if(i.lt.1 .or. i.gt.nlog) goto 10
            alen = min(s1log,real(i+1)) - max(s0log,real(i))
            flux = fsrc(l) * alen/(s1log-s0log) * dnu
            fdest(i) = fdest(i) + flux
 10      continue
      end do

      return
      end



************************************************************************
* subroutine zerospec -- Initialize flux array
************************************************************************

      subroutine zerospec(nw,flux)
      real flux(nw)
      do i=1, nw
         flux(i) = 0
      end do
      return
      end



************************************************************************
* subroutine addspec -- Add two flux arrays
************************************************************************

      subroutine addspec(nw,f1,f2)
      real f1(nw),f2(nw)
      do i=1, nw
         f1(i) = f1(i) + f2(i)
      end do
      return
      end



************************************************************************
* subroutine meanzero -- Divide flux array by spline
************************************************************************

      subroutine meanzero(n,l1,l2,y,ynorm,ioff,nknot,xknot,yknot)
* Array Y(N) has spline fitted between l1 and l2 and divided out
* Apply ioff as an offset in picking knots if > 0
      parameter (MAXKNOT=20)
      real xknot(1),yknot(1)
      logical err
      real y(n),ynorm(n)

* KNOTCHOICE = 1 means maximum; KNOTCHOICE = 2 means ave
      KNOTCHOICE = 2
* KNOTNUM is the normal number of knots over the N pixels
      KNOTNUM = 13

* Copy the array
      do i=1, n
         ynorm(i) = y(i)
      end do

* First find the range of non-zero data values, nuking edge pixels
      nedge = 1
      l1 = 1
      nuke = 0
 1    if(l1.lt.n .and. (y(l1).le.0 .or. nuke.lt.nedge)) then
         if(y(l1).gt.0) nuke = nuke + 1
         ynorm(l1) = 0
         l1 = l1 + 1
         goto 1
      end if
      
      l2 = n
      nuke = 0
 2    if(l2.gt.1 .and. (y(l2).le.0 .or. nuke.lt.nedge)) then
         if(y(l2).gt.0) nuke = nuke + 1
         ynorm(l2) = 0
         l2 = l2 - 1
         goto 2
      end if

      if(l2-l1 .lt. 3*KNOTNUM) then
         write(6,*) 'MEANZERO: This spectrum is zero!'
         return
      end if

* KNOTCHOICE = 1 means pick the maximum
      if(KNOTCHOICE.eq.1) then
* How many knots to use for the spline?
         nknot = max(5, (KNOTNUM*(l2-l1))/n)

* Choose knots for the spline
         top = 0
         do k=1, nknot
            if(k.eq.1) then
               i1 = l1
               i2 = l1
            else if(k.eq.nknot) then
               i1 = l2
               i2 = l2
            else
               i1 = (k-1.5)*(l2-l1+1)/float(nknot-1) + l1
               i2 = (k-0.5)*(l2-l1+1)/float(nknot-1) + l1 - 1
            end if
            biggie = -500
            do i=i1, i2
               if(y(i).gt.biggie) then
                  xknot(k) = i
                  yknot(k) = alog10(y(i))
                  biggie = y(i)
                  top = amax1(top, biggie)
               end if
            end do
C     WRITE(6,'(2I6,2F10.3)') I1, I2, XKNOT(K), YKNOT(K)
         end do

* KNOTCHOICE = 2 means use the average
      else if(KNOTCHOICE.eq.2) then
         nknot = 0
         kwidth = n / KNOTNUM
         nave = 0
         wave = 0
         fave = 0
         istart = 0
         if(ioff.gt.0) istart = mod(ioff,kwidth) - kwidth
C         WRITE(6,*) N, KWIDTH, ISTART
         do i=1, n
            if(i.gt.l1 .and. i.lt.l2) then
               nave = nave + 1
               wave = wave + i-0.5
               fave = fave + y(i)
            end if
            if(mod(i-istart,kwidth) .eq. 0) then
               if(nave.gt.0 .and. fave.gt.0) then
                  nknot = nknot + 1
                  xknot(nknot) = wave / nave
                  yknot(nknot) = alog10(fave/nave)
C               WRITE(6,*) NKNOT, XKNOT(NKNOT), YKNOT(NKNOT)
               end if
               nave = 0
               wave = 0
               fave = 0
            end if
         end do

      else
         write(6,*) 'MEANZERO: must make a choice for KNOTCHOICE'
         stop
      end if

* Calculate the spline
      call spline(xknot,yknot,nknot,err)
      if(err) then
         write(6,*) 'MEANZERO: Error determining spline'
         return
      end if

* Evaluate the spline at each point, divide and subtract 1
      do i=l1, l2
         call splineval(i-0.5,spl,xknot,yknot,nknot)
         ynorm(i) = y(i) / 10.0**spl - 1
C         ynorm(i) = (y(i) - 10**spl) / top
      end do

      return
      end



************************************************************************
* subroutine splinedex --  Fill an array with spline evaluations, 
*                          and dex the lot
************************************************************************

      subroutine splinedex(n,y,loff,nknot,xknot,yknot)
      real xknot(nknot),yknot(nknot),y(n)
      logical err
      call spline(xknot,yknot,nknot,err)
      do i=1, n
C         call splineval(i-0.5+loff, spl, xknot, yknot, nknot)
C Nah, hold to constant off of the region where the fit took place.
         call splineval(max(1,min(n,i+loff))-0.5,spl,xknot,yknot,nknot)
         y(i) = 10.0**spl
      end do
      return
      end



************************************************************************
* subroutine overlap --  Calculate overlap and trim shifted buffers
************************************************************************

      subroutine overlap(n,shift,lap,x0,x1,y0,y1,percent)
      real lap
      real y0(n),x0(n),x1(n),y1(n)
* Y has SHIFT to the right wrt X

* Find left edge of X
      lx1 = 1
 1    if(lx1.lt.n .and. x0(lx1).eq.0) then
         lx1 = lx1 + 1
         goto 1
      end if
      
* Find right edge of X
      lx2 = n
 2    if(lx2.gt.1 .and. x0(lx2).eq.0) then
         lx2 = lx2 - 1
         goto 2
      end if
      
* Find left edge of Y
      ly1 = 1
 3    if(ly1.lt.n .and. y0(ly1).eq.0) then
         ly1 = ly1 + 1
         goto 3
      end if
      
* Find right edge of Y
      ly2 = n
 4    if(ly2.gt.1 .and. y0(ly2).eq.0) then
         ly2 = ly2 - 1
         goto 4
      end if
      
      ishift = nint(shift)
* Desired edges of X and Y
      mx1 = max(lx1, ly1-ishift)
      mx2 = min(lx2, ly2-ishift)
      my1 = max(ly1, lx1+ishift)
      my2 = min(ly2, lx2+ishift)
      
      lap = mx2 - mx1 + 1
C      WRITE(6,'(9I8)') ISHIFT, LX1, LX2, MX1, MX2, LY1, LY2, MY1, MY2
* Copy and reapodize
      do i=1, n
         if(i.gt.mx1 .and. i.lt.mx2) then
            x1(i) = x0(i)
         else
            x1(i) = 0
         end if
         if(i.gt.my1 .and. i.lt.my2) then
            y1(i) = y0(i)
         else
            y1(i) = 0
         end if
      end do
      call apodize(n,mx1,mx2,x1,percent)
      call apodize(n,my1,my2,y1,percent)

      return
      end



************************************************************************
* function xlog -- Return alog10(x)
************************************************************************

      function xlog(x)
      if(x.lt.0) then
         xlog = -alog10(-x)
      else if(x.eq.0) then
         xlog = 0
      else
         xlog = alog10(x)
      end if
      return
      end



************************************************************************
* subroutine despace -- Remove blancks from a string
************************************************************************

      subroutine despace(s)
      character*(*) s
      l = len(s)
      do j=1, l
         if(s(1:1).eq.' ') then
            do i=2, l-1
               s(i-1:i-1) = s(i:i)
            end do
            s(l:l) = ' '
         else
            return
         end if
      end do
      return
      end



************************************************************************
* subroutine medfilt -- Replace data with a running median
************************************************************************

      subroutine medfilt(n,data,medlen,buf,tmp)
* Uses buffers buf and tmp
* Dumb version which just does a box-car median over available pixels
      real data(n),buf(n),tmp(1)
      do i=1, n
         buf(i) = data(i)
      end do
      medwidth = medlen / 2
      do k=1, n
         nmed = 0
         do i=max(1,k-medwidth), min(n,k+medwidth)
            nmed = nmed + 1
            tmp(nmed) = buf(i)
         end do
         data(k) = amedian(nmed, tmp)
      end do
      return
      end



************************************************************************
* subroutine medwfilt -- Replace data with a weighted running median
************************************************************************

      subroutine medwfilt(n,wave,data,fwmed,buf,tmp)
* Uses buffers buf and tmp
* Quasi-Gaussian weighted median filter over a width in wavelength
      parameter (MAXDUP=3)
      real wave(n),data(n),buf(n),tmp(n)
      real brkpt(MAXDUP)
* Stash a copy 
      do i=1, n
         buf(i) = data(i)
      end do

* Gaussian width of our FWHM
      sig = fwmed / 2.35

* Where are the count break points?
      do i=1, MAXDUP
         brkpt(i) = sig * sqrt(2*alog(MAXDUP/(i-0.5)))
      end do

      do k=1, n
         nmed = 0
         do 10 i = 1,n
            if(abs(wave(i)-wave(k)) .gt. brkpt(1)) goto 10
            do j=1, MAXDUP
               if(abs(wave(i)-wave(k)) .le. brkpt(j)) then
                  nmed = nmed + 1
                  if(nmed .gt. n) then
                     write(6,*) 'MEDWFILT: not enough array for this ',
     $                    'smoothing width', fwmed
                     return
                  end if
                  tmp(nmed) = buf(i)
               end if
            end do
 10      continue
         data(k) = amedian(nmed, tmp)
      end do
      return
      end
