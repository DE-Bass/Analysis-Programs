* Routine to fit a linear function
* Copyright (C) 1980 John L. Tonry
*
* 800905 v1.0 John Tonry 

      subroutine linearfit(npt, y, nx, x, param)
      parameter (maxx=20)
      real x(nx,npt), y(npt), param(nx)
C Y = data being fitted as a linear function of X(I,*)
C PARAM = returned parameters of fit:
C	Y(*) = PARAM(1)*X(1,*)
C	     + PARAM(2)*X(2,*)
C	     + PARAM(3)*X(3,*)   etc.
C
      real*8 v(maxx), am(maxx*maxx), par8(maxx)
      if(nx.gt.maxx) then
         write(6,*) 'Too many parameters requested:', nx, maxx
         return
      end if

      if(npt.lt.nx) then
         write(6,*) 'Too few points provided:', npt, nx
         return
      end if

      do 10 j = 1,nx
         v(j) = 0
         do 11 i = 1,nx
            am(i+(j-1)*nx) = 0
 11      continue
 10   continue

C Accumulate sums
      do 20 i = 1,npt
         do 22 j = 1,nx
            v(j) = v(j) + dble(y(i))*dble(x(j,i))
            do 23 k = 1,j
               am(k+(j-1)*nx) = am(k+(j-1)*nx) + 
     $              dble(x(j,i))*dble(x(k,i))
 23         continue
 22      continue
 20   continue

C Fill in symmetrical matrix
      do 30 l = 1,nx-1
         do 31 k = l+1,nx
            am(k+(l-1)*nx) = am(l+(k-1)*nx)
 31      continue
 30   continue

      call solve(nx, v, par8, am)

      do 40 i = 1,nx
         param(i) = sngl(par8(i))
 40   continue

      return
      end


* Routine to do a weighted fit a linear function
      subroutine wlinearfit(npt, y, w, nx, x, param)
      parameter (maxx=20)
      real x(nx,npt), w(npt), y(npt), param(nx)
C Y = data being fitted as a linear function of X(I,*)
C PARAM = returned parameters of fit:
C	Y(*) = PARAM(1)*X(1,*)
C	     + PARAM(2)*X(2,*)
C	     + PARAM(3)*X(3,*)   etc.
C
      real*8 v(maxx), am(maxx*maxx), par8(maxx)
      if(nx.gt.maxx) then
         write(6,*) 'Too many parameters requested:', nx, maxx
         return
      end if

      if(npt.lt.nx) then
         write(6,*) 'Too few points provided:', npt, nx
         return
      end if

      do 10 j = 1,nx
         v(j) = 0
         do 11 i = 1,nx
            am(i+(j-1)*nx) = 0
 11      continue
 10   continue

C Accumulate sums
      do 20 i = 1,npt
         do 22 j = 1,nx
            v(j) = v(j) + dble(y(i))*dble(x(j,i))*dble(w(i))
            do 23 k = 1,j
               am(k+(j-1)*nx) = am(k+(j-1)*nx) + 
     $              dble(x(j,i))*dble(x(k,i))*dble(w(i))
 23         continue
 22      continue
 20   continue

C Fill in symmetrical matrix
      do 30 l = 1,nx-1
         do 31 k = l+1,nx
            am(k+(l-1)*nx) = am(l+(k-1)*nx)
 31      continue
 30   continue

      call solve(nx, v, par8, am)

      do 40 i = 1,nx
         param(i) = sngl(par8(i))
 40   continue

      return
      end

      subroutine solve(n, y, x, a)
      implicit real*8 (a-h,o-z)
      real*8 x(n), y(n), a(n*n)

C Solve matrix by Gaussian elimination
      do 40 j = 1,n-1
         do 41 i = j+1,n
            r = a(i+(j-1)*n) / a(j+(j-1)*n)
            y(i) = y(i) - r*y(j)
            do 42 k = j+1,n
               a(i+(k-1)*n) = a(i+(k-1)*n) - r*a(j+(k-1)*n)
 42         continue
 41      continue
 40   continue

C Back substitute to solve for parameters
      x(n) = y(n) / a(n*n)
      do 50 i = 1,n-1
         j = n - i
         x(j) = y(j)
         do 51 k = j+1,n
            x(j) = x(j) - a(j+(k-1)*n)*x(k)
 51      continue
         x(j) = x(j) / a(j+(j-1)*n)
 50   continue

      return
      end
