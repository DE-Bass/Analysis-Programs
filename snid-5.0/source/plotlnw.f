************************************************************************
**                                                                    **
**                              Plotlnw                               **
**                                                                    **
**                                                                    **
**                Plot .lnw files output from Logwave                 **
**                                                                    **
**                                                                    **
**                Copyright (C) 2007 Stephane Blondin                 **
**                                                                    **
************************************************************************
*
************************************************************************
**                                                                    **
**                         Contact information                        **
**                                                                    **
************************************************************************
*
* Stephane Blondin
* Harvard-Smithsonian Center for Astrophysics
* Optical and Infrared Astronomy Division
* 60 Garden Street, MS-20
* Cambridge, MA 02138
* USA
* Email: sblondin@cfa.harvard.edu
*
************************************************************************
**                                                                    **
**                GNU GENERAL PUBLIC LICENSE disclaimer               **
**                                                                    **
**                       (see also gpl-3.0.txt)                       **
**                                                                    **
************************************************************************
*
* This program is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program.  If not, see <http://www.gnu.org/licenses/>.
*
************************************************************************
**                                                                    **
**                          Revision history                          **
**                                                                    **
************************************************************************
*
* 20 Jul 2007 - created (Stephane Blondin)
*
************************************************************************


      program plotlnw

      implicit none

      include 'snid.inc'

c Further declarations
      integer i,j,k
      integer ittype,itstype
      integer ifilt
      integer ilhs,irhs
      integer idx0,idx1
      character*200 inbuf         
      character*100 ftemp,tmpfout 
      character*12 tname          
      character*10 ttype          
      integer nepoch
      integer nwx
      real w0x,w1x,ww
      real dta
      integer mostknots
      integer nk(MAXEPOCH)
      real fmean
      real xk(MAXKNOT,MAXEPOCH),yk(MAXKNOT,MAXEPOCH)
      real tmpep(MAXEPOCH),ep(MAXEPOCH)
      integer tf
      real tmp(MAXLOG),temp(MAXLOG,MAXTEMP)
      complex tft(MAXLOG)
      character*200 psfile
      integer ixmin,ixmax
      real yoff
      integer ierr,jps,pgopen

c For plotting
      integer pgbeg,ier
      real xmin,xmax,ymin,ymax
      character*100 tit,xtit,ytit
      character*10 cnepoch
      character*20 label
      real ylab

c Functions
      integer iargc,lnb,nlb

c Version number
      character*20 vnum
      data vnum /'1.0, 20 July 2007'/


      
c
c Display syntax on call
c
      if(iargc().lt.1) then
         write(6,*) 'Plotlnw (v'//vnum(:lnb(vnum))//')'
         write(6,*) 'Copyright (C) 2007 S. Blondin'
         write(6,*) ' '
         write(6,*) 'Usage: plotlnw [options] template.lnw'
         write(6,*) ' '
         write(6,*) 'options (default):'
         write(6,*) '    filter=0/1        apply bandpass filter (0)'   
         write(6,*) '    k{1,2,3,4}=       bandpass filter range',
     $        ' (1,4,85,102)'
         write(6,*) '    x{min,max}=       abscissa range (template)'
         write(6,*) '    yoff=             vertical offset between',
     $        ' spectra (1)'
         write(6,*) '    ps=0/1            create postscript output (0)'   
         write(6,*) ' '
         write(6,*) '[see file Howto.snid for more info]'
         stop
      end if



c
c Default values for parameters (override with options)
c

* Postscipt output?
      ips = 0

* Wavenumbers for filtering 
      k1 = 0
      k2 = 0
      k3 = 0
      k4 = 0

* Abscissa range
      ixmin = 0
      ixmax = 0

* Vertical offset between spectra
      yoff = 1.

c
c Have we got some options?
c
      do 10 i=1, iargc()         
         noption = i
         call getarg(noption,inbuf)
         if(index(inbuf, '=') .eq. 0) goto 20

         if(inbuf(:3) .eq. 'ps=') then
            read(inbuf(4:), '(i20)', err=10) ips

         else if(inbuf(:7) .eq. 'filter=') then
            read(inbuf(8:), '(i20)', err=10) ifilt

         else if(inbuf(:3) .eq. 'k1=') then
            read(inbuf(4:), '(i20)', err=10) k1

         else if(inbuf(:3) .eq. 'k2=') then
            read(inbuf(4:), '(i20)', err=10) k2

         else if(inbuf(:3) .eq. 'k3=') then
            read(inbuf(4:), '(i20)', err=10) k3

         else if(inbuf(:3) .eq. 'k4=') then
            read(inbuf(4:), '(i20)', err=10) k4

         else if(inbuf(:5) .eq. 'xmin=') then
            read(inbuf(6:), '(g20.0)', err=10) xmin
            ixmin = 1

         else if(inbuf(:5) .eq. 'xmax=') then
            read(inbuf(6:), '(g20.0)', err=10) xmax
            ixmax = 1

         else if(inbuf(:5) .eq. 'yoff=') then
            read(inbuf(6:), '(g20.0)', err=10) yoff

         else
            write(6,*) 'Illegal option: ',inbuf(:index(inbuf,"=")-1)
            stop
         end if
 10   continue
 20   noption = noption - 1



c
c Read input template file
c
      call getarg(noption+1,ftemp)
      call getroot(ftemp,tmpfout,froot)
      open(unit=1,file=ftemp,status='old')

* Read header information
      read(1,1000) nepoch,nwx,w0x,w1x,mostknots,tname,dta,ttype,ittype,
     $     itstype
 1000 format(2i5,2f10.2,i7,5x,a12,f7.2,2x,a10,2(i3))

* Wavenumbers for filtering 
      if(k1.eq.0) k1 = 1
      if(k2.eq.0) k2 = 4
      if(k3.eq.0) k3 = nwx/12
      if(k4.eq.0) k4 = nwx/10

* Read spline information
      read(1,*) i,(nk(j),fmean,j=1,nepoch)
      do i=1, mostknots
         read(1,*) k,(xk(i,j),yk(i,j),j=1,nepoch)
      end do

* Read flag and ages
      read(1,*) tf, (ep(i),i=1,nepoch)

* Read the normalized spectra
      do i=1, nwx
         read(1,*) ww, (tmpep(j),j=1,nepoch)
         wave(i) = ww
         do j=1, nepoch
            temp(i,j) = tmpep(j)
         end do
      end do

* Close file
      close(1)

* Apply bandpass filter?
      if(ifilt.gt.0) then
         do j=1, nepoch
c find non-zero ends
            ilhs = 0
            irhs = 0
            idx0 = 1
            idx1 = nwx
            do i=1, nwx
               if(ilhs.eq.0.and.temp(i,j).eq.0) then
                  idx0 = idx0 + 1
               else
                  ilhs = 1
               end if
               if(irhs.eq.0.and.temp(nwx+1-i,j).eq.0) then
                  idx1 = idx1 - 1
               else
                  irhs = 1
               end if
            end do
            do i=1, nwx
               tmp(i) = temp(i,j)
            end do
            call rcvector(nwx,tmp,tft)
            call four2(tft,nwx,1,+1,1,ierr)
            call cvscale(nwx,tft,float(nwx))
            call filter(nwx,k1,k2,k3,k4,tft)
            call four2(tft,nwx,1,-1,1,ierr)
            call cvscale(nwx,tft,1.0/float(nwx))
            call crvector(nwx,tft,tmp)           
            do k=idx0, idx1
               temp(k,j) = tmp(k)
            end do
         end do
      end if


c
c Plot
c
      
* Open graphics device
      jps = 0
      ier = pgbeg(0,'/XSERVE',1,1)
      if(ier.ne.1) then
         write(6,*) 'Could not open output graphics device for PGPLOT'
         stop
      endif
     
* Axis ranges
      if(ixmin.eq.0) xmin = w0x
      if(ixmax.eq.0) xmax = w1x
      ymin = 0.
      ymax = (nepoch+1)*yoff

* Plot box
      write(cnepoch,'(i3)') nepoch
      tit ='SNID template '//tname(:lnb(tname))//' ('//
     $     ttype(:lnb(ttype))//') ; '//
     $     cnepoch(nlb(cnepoch):lnb(cnepoch))//' spectra'
      xtit='Rest Wavelength [\\A]'
      ytit='Flattened Flux'
      if(ifilt.gt.0) tit=tit(:lnb(tit))//' - filtered'
 100  call pgslw(3)
      call pgenv(xmin,xmax,ymin,ymax,0,0)
      call pglabel(xtit,ytit,tit)

* Plot spectra and labels
      do j=1, nepoch
         write(label,*) nint(ep(j))
         label = label(nlb(label):lnb(label))
         if(nint(ep(j)).ge.0) label='+'//label(nlb(label):lnb(label))
         label = label(nlb(label):lnb(label))
         do i=1, nwx
            tmp(i) = temp(i,nepoch+1-j) + j*yoff
         end do
         call pgline(nwx,wave,tmp)
         ylab = (nepoch-j+1.35)/(nepoch+1.)
         call pgmtxt('RV',-3.0,ylab,0,label)
      end do

* Postscipt output?
      if(ips.gt.0) then
         if(jps.eq.0) then
            psfile = tmpfout(:lnb(tmpfout))//'.ps'
            istat = pgopen(psfile//'/ps')
            if(istat.le.0) then
               write(6,*) 'PGOPEN: Could not open PS device'
               stop
            endif
            jps = 1
            goto 100
         else 
            call pgclos
            write(6,'(2a)') ' Created PS file: ',psfile(:lnb(psfile))
         end if
      end if

* End PGPLOT
      call pgend






************************************************************************
**                                                                    **
**                       End of program plotlnw                       **
**                                                                    **
************************************************************************

 999  stop
      end

