* Various median functions
* Copyright (C) 1980 John L. Tonry

C Return the median of a buffer
      function amedian(n,x)
      real x(1)
      if(n.eq.0) then
         amedian = 0
         return
      else if(n.eq.1) then
         amedian = x(1)
         return
      end if
      do 910 k = 2,n
         xx = x(k)
         l = k - 1
 911     if(l.ge.1.and.x(l).gt.xx) then
            x(l+1) = x(l)
            l = l - 1
            goto 911
         end if
         x(l+1) = xx
 910  continue
      amedian = 0.5*(x((n+1)/2)+x((n+2)/2))
      return
      end

C Return the median of a buffer and sort the index along with it
      function amedidx(n,x,idx)
      real x(1)
      integer idx(1)
      if(n.eq.0) then
         amedidx = 0
         return
      else if(n.eq.1) then
         amedidx = x(1)
         return
      end if
      do 910 k = 2,n
         xx = x(k)
         isave = idx(k)
         l = k - 1
 911     if(l.ge.1.and.x(l).gt.xx) then
            x(l+1) = x(l)
            idx(l+1) = idx(l)
            l = l - 1
            goto 911
         end if
         x(l+1) = xx
         idx(l+1) = isave
 910  continue
      amedidx = 0.5*(x((n+1)/2)+x((n+2)/2))
      return
      end

      function median(n,buf)
* Get the median by counting the pixels, assume unsigned 16 bit
      integer*2 buf(n)
      integer count(0:65535)

      do 10 i = 0,65535
         count(i) = 0
 10   continue

      do 20 j = 1,n
         idx = buf(j)
 20   continue

      median = -1
      now = 0
 30   if(now.lt.n/2) then
         median = median + 1
         now = now + count(median)
         goto 30
      end if

      return
      end


