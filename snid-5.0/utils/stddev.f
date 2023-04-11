* Compute the standard deviation of an array
* Copyright (C) 2007 Stephane Blondin

      real function stddev(n,arr)

      implicit none

      integer i
      integer n
      real mean,var,arr(n)
      real sum,res,sumres

      if(n.lt.2) then
         write(6,*) 'stddev: need at least 2 elements to compute stddev'
         stop
      end if

c Compute SUM(arr) and mean(arr)
      sum=0.
      do i=1, n
        sum = sum + arr(i)
      end do
      mean = sum / n

c Compute mean residuals, sum thereof and sum thereof^2
      res    = 0.
      sumres = 0.
      var    = 0.
      do i=1, n
        res    = arr(i) - mean
        sumres = sumres + res
        var    = var + res*res
      end do

c Compute variance and stddev
      var = (var - sumres*sumres / n) / (n-1)
      stddev = sqrt(var)

      return
      end
