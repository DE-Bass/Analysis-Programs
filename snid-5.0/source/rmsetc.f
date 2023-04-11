************************************************************************
**                                                                    **
**                Copy buffers and compute rms's, etc.                **
**                                                                    **
**                                                                    **
**                  Copyright (C) 1999 John L. Tonry                  **
**                                                                    **
************************************************************************

c Included in this file:
c
c subroutine rrvector
c subroutine rvscale
c subroutine cvscale
c subroutine ccvector
c subroutine rcvector
c subroutine crvector
c subroutine correlate
c function getrms
c subroutine aspart
c subroutine shiftit
c function rmsfilter
c subroutine filter
c subroutine phaseband



************************************************************************
* subroutine rrvector -- Copy a real buffer to another
************************************************************************

      subroutine rrvector(n,buf1,buf2)
      real buf1(1),buf2(1)
      do j=1, n
         buf2(j) = buf1(j)
      end do
      return
      end



************************************************************************
* subroutine rvscale -- Multiply a real buffer by a constant
************************************************************************

      subroutine rvscale(n,buf1,value)
      real buf1(1)
      do j=1, n
         buf1(j) = buf1(j) * value
      end do
      return
      end



************************************************************************
* subroutine cvscale -- Multiply a complex buffer by a constant
************************************************************************

      subroutine cvscale(n,x1,value)
      complex x1(1)
      do j=1, n
         x1(j) = x1(j) * value
      end do
      return
      end



************************************************************************
* subroutine ccvector -- Copy a complex buffer to another
************************************************************************

      subroutine ccvector(n,x1,x2)
      complex x1(1),x2(1)
      do j=1, n
         x2(j) = x1(j)
      end do
      return
      end



************************************************************************
* subroutine rcvector -- Copy a real buffer to a complex one
************************************************************************

      subroutine rcvector(n,buf,x)
      real buf(1)
      complex x(1)
      do j=n, 1, -1
         x(j) = cmplx(buf(j),0.0)
      end do
      return
      end



************************************************************************
* subroutine crvector -- Copy a complex buffer to a real one
************************************************************************

      subroutine crvector(n,x,buf)
      real buf(1)
      complex x(1)
      do j=1, n
         buf(j) = real(x(j))
      end do
      return
      end



************************************************************************
* subroutine correlate -- Compute correlation function
************************************************************************

      subroutine correlate(n,ft1,ft2,cfn)
      complex ft1(1),ft2(1),cfn(1)
      do j=1, n
         cfn(j) = ft1(j) * conjg(ft2(j))
      end do
      return
      end



************************************************************************
* function getrms -- Compute the rms of a complex buffer
************************************************************************

      function getrms(n,x)
* Assume a transform of a real function
      complex x(1)
      getrms = cabs(x(1))*cabs(x(1))+cabs(x(n/2+1))*cabs(x(n/2+1))
      do j=2, n/2
         getrms = getrms + 2 * cabs(x(j)) * cabs(x(j))
      end do
* Divide by N since it is the transform
      getrms = sqrt(getrms) / n
      return
      end



************************************************************************
* subroutine aspart -- Compute the (anti)symmetric parts of a complex buffer, 
*                      shifted by shift
************************************************************************

      subroutine aspart(n,k1,k2,k3,k4,shift,x,arms,srms)
* Assume a transform of a real function, and a cosine bell of k1,k2,k3,k4.
      complex x(n), phase
      parameter (pi = 3.14159265)

      arms = 0
      srms = 0
      do k=k1 ,k4
         angle = -2*pi*k*shift/n
         phase = cmplx(cos(angle),sin(angle))
         if(k.eq.0.or.k.eq.n/2) then
            f = 1
         else
            f = 2
         end if
         if(k.lt.k2) then
            arg = pi * float(k-k1)/float(k2-k1)
            f = f * .25 * (1-cos(arg)) * (1-cos(arg))
         else if(k.gt.k3) then
            arg = pi * float(k-k3)/float(k4-k3)
            f = f * .25 * (1+cos(arg)) * (1+cos(arg))
         end if
         arms = arms + f*aimag(phase*x(k+1))*aimag(phase*x(k+1))
         srms = srms + f*real(phase*x(k+1))*real(phase*x(k+1))
      end do
* Divide by N since it is the transform
      arms = sqrt(arms) / n
      srms = sqrt(srms) / n
      return
      end



************************************************************************
* subroutine shiftit -- Shift a complex buffer by -shift
************************************************************************

      subroutine shiftit(n,shift,x)
      complex x(n), phase
      parameter (pi = 3.14159265)
      do i=1, n
         if(i-1.le.n/2) then
            angle = -2*pi*(i-1)*shift/n
         else
            angle = -2*pi*(i-1-n)*shift/n
         endif
         phase = cmplx(cos(angle),sin(angle))
         x(i) = x(i) * phase
      end do
      return
      end



************************************************************************
* function rmsfilter -- Return the rms of a real function passed through 
*                       a bandpass filter
************************************************************************

      function rmsfilter(n,k1,k2,k3,k4,x)
* The bandpass is the square root of a cosine bell.
      complex*8 x(1)
      parameter (pi=3.14159)
      rmsfilter = 0
      do k=k1, k4
         if(k.eq.0.or.k.eq.n/2) then
            f = 1
         else
            f = 2
         end if
         if(k.lt.k2) then
            arg = pi * float(k-k1)/float(k2-k1)
            factor = .5 * (1-cos(arg))
         else if(k.gt.k3) then
            arg = pi * float(k-k3)/float(k4-k3)
            factor = .5 * (1+cos(arg))
         else
            factor = 1
         end if
         rmsfilter = rmsfilter + f * factor *
     1        (real(x(k+1))*real(x(k+1)) + aimag(x(k+1))*aimag(x(k+1)))
      end do
* Divide by N because it is the transform
      rmsfilter = sqrt(rmsfilter) / n
      return
      end



************************************************************************
* subroutine filter -- Filter a FT by multiplying by a cosine bell
************************************************************************

      subroutine filter(n,k1,k2,k3,k4,x)
      complex*8 x(1)
      parameter (pi = 3.14159265)
      do j=1, n
         if(j.lt.n/2+1) then
            number = j - 1
         else
            number = j - n - 1
         end if
         numa = abs(number)
         if(numa.lt.k1.or.numa.gt.k4) then
            factor = 0
         else if(numa.lt.k2) then
            arg = pi * float(numa-k1)/float(k2-k1)
            factor = .5 * (1-cos(arg))
         else if(numa.gt.k3) then
            arg = pi * float(numa-k3)/float(k4-k3)
            factor = .5 * (1+cos(arg))
         else
            factor = 1
         end if
         x(j) = x(j) * factor
      end do
      return
      end



************************************************************************
* subroutine phaseband -- Computes the phase of a complex function
************************************************************************

      subroutine phaseband(n,shift,x,phase)
* assumption that x(*) is within +/- pi of 2*pi*SHIFT*j/N
      real phase(1)
      complex x(1)
      parameter (pi = 3.14159265)
      do i=1, n/2
         angle = 2*pi*(i-1)*shift/n
         if(x(i).eq.0) then
            phase(i) = angle
         else
            phase(i) = atan2(aimag(x(i)),real(x(i)))
         end if
         npi = (phase(i)-angle) / pi
         if(npi.gt.0) then
            ncycle = (npi + 1) / 2
         else
            ncycle = (npi - 1) / 2
         end if
         phase(i) = phase(i) - 2*pi*ncycle - angle
         if(i.lt.n/2) phase(n-(i-1)) = phase(i) + angle
C         phase(i) = phase(i) - 2*pi*ncycle
C         if(i.lt.n/2) phase(n-(i-1)) = phase(i) - angle
C         phase(i) = phase(i) - 2*pi*ncycle
C         if(i.lt.n/2) phase(i+n/2) = phase(i) - angle
      end do
      return
      end

