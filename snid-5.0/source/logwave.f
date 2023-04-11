************************************************************************
**                                                                    **
**                            L o g w a v e                           **
**                                                                    **
**                                                                    **
**     Convert (w,f) data files to a common log wavelength scale      **
**                                                                    **
**                                                                    **
**    Copyright (C) 1999-2007 Stephane Blondin and John L. Tonry      **
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
* John L. Tonry
* Institute for Astronomy
* 2680 Woodlawn Dr.
* Honolulu, HI 96822
* USA
* Email: jt@ifa.hawaii.edu
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
* 31 Jul 2007, rev 2.0 - added various header info to match SNID v5.0
*     added age flag (SB, 31 Jul 2007)
*     call typeinfo (SB, 21 May 2007)
*     added options, implicit none (SB, 05 Mar 2007)
*     added type info (SB, 24 Apr 2006)
*
* Aug-Oct 1999 - created (John Tonry)
*
************************************************************************


      program logwave

      implicit none

      include 'snid.inc'

      integer i,j,k,l
      character*200 inbuf
      integer iep,maxep,nsn,nsntot,nfile,ntok,nwave
      character*100 fname(MAXTEMP)
      character*100 fspec
      character*200 line
      real zfactor,wmean
      integer meconvert,mostknots,nmean
      integer l1,l2
      real redshift(MAXTEMP)
      integer nepoch(MAXSN),iepoch(MAXSN),iab(MAXTEMP)
      real flog(MAXLOG,MAXEPOCH),wlog(MAXLOG+1)
      real fnorm(MAXLOG,MAXEPOCH)
      real fmean(MAXEPOCH)
      integer nknot(MAXEPOCH)
      real xknot(MAXKNOT,MAXEPOCH),yknot(MAXKNOT,MAXEPOCH)
      real wavemin(MAXTEMP),wavemax(MAXTEMP)
      real delta(MAXSN)
c Functions
      integer iargc,lnb,nlb

c Version number
      character*20 vnum
      data vnum /'2.0, 31 July 2007'/

* Print syntax on call
      if(iargc().lt.1) then
         write(6,*) 'Logwave (v'//vnum(:lnb(vnum))//')'
         write(6,*) 'Copyright (C) 1999-2007 S. Blondin and J. L. Tonry'
         write(6,*) ' '
         write(6,*) 'Usage: logwave [options] list'
         write(6,*) ' '
         write(6,*) 'list format:'
         write(6,*) '    Object Type Filename Age Age_flag Delta AB?',
     $        ' Wave_range [Redshift]'
         write(6,*) 'creates output files: Object.lnw'
         write(6,*) ' '
         write(6,*) 'options (default):'
         write(6,*) '    w0= w1=  wavelength range (2500,10000)'
         write(6,*) '    nw=      number of log bins (1024)'
         write(6,*) ' '
         write(6,*) '[see file Howto.snid for more info]'
         stop
      end if

* Get type info
      call typeinfo

* Default values for parameters (override with options)
      w0 = 2500
      w1 = 10000
      nw = 1024

* How much to apodize?
      percent = 5.0

* Have we got some options?
      do 1 i=1, iargc()-1
         noption = i
         call getarg(noption,inbuf)
         if(index(inbuf, '=') .eq. 0) goto 2
         if(inbuf(:3)      .eq. 'w0=') then
            read(inbuf(4:), '(g20.0)', err=1) w0
         else if(inbuf(:3) .eq. 'w1=') then
            read(inbuf(4:), '(g20.0)', err=1) w1
         else if(inbuf(:3) .eq. 'nw=') then
            read(inbuf(4:), '(i20)'  , err=1) nw
         end if
 1    continue
 2    noption = noption - 1

      if(nw.gt.MAXLOG) then
         write(6,*) 'ERROR! Allocate more space for MAXLOG'
         write(*,1000) nw,MAXLOG
 1000    format(' (nw,MAXLOG) = (',i5,',',i5,')')
         stop
      end if

* Set up the log wavelength array
      dwlog = alog(w1/w0) / nw
      do i=1, nw+1
         wlog(i) = w0 * exp(float(i-1)*dwlog)
      end do

* Read the list of data
      fspec = 'list'
      call getarg(noption+2,fspec)
      call getroot(fspec,fdata,froot)
      write(6,*) 'Reading data files from: ', fdata(:lnb(fdata))
      open(unit=1,file=fspec,status='old')
C      read(1,*)
      nsn = 0
      maxep = 0
      nfile = 0
      nsntot = 0

c
c Loop over files
c
      do 4 i=1, MAXTEMP
         read(1,'(a)',end=5) line

* Skip line if commented out (#)
         if(index(line(nlb(line):nlb(line)+1),'#').gt.0) goto 4

* Read in wavelength and redshift info
         call parser(line,MAXTOK,ntok,ltok,tok)
         read(tok(8),*) wavemin(nfile+1)
         read(tok(9),*) wavemax(nfile+1)
         redshift(nfile+1) = 0
         if(ntok.gt.9) read(tok(10),*) redshift(nfile+1)

* Skip this SN if there's no wavelength range
         if(wavemin(nfile+1).gt.w1.or.wavemax(nfile+1).lt.w0) then
            write(6,*) 'WARNING! No wavelength range for file',
     $           tok(3)(:lnb(tok(3)))
            write(6,*) 'w0/w1 = ',w0,w1
            write(6,*) 'wmin/wmax = ',wavemin(nfile+1),wavemax(nfile+1)
            goto 3
         end if
         nfile = nfile + 1

* Read in name, type, filename, age, AB?
         sname(nfile)= tok(1)
         stype(nfile)= tok(2)
         do it=1, NT
            do ist=2, NST
               if(index(typename(it,ist)(:lnb(typename(it,ist))),
     $              stype(nfile)(:lnb(stype(nfile)))).gt.0) then
                  itype(1,nfile) = it
                  itype(2,nfile) = ist
                  goto 3
               end if
            end do
         end do
         write(6,'(2a)') ' ERROR! Unknown type: ',stype(nfile)
         goto 4
 3       continue
         fname(nfile)= tok(3)
         read(tok(4),*) epoch(nfile)
         read(tok(5),*) tflag(nfile)
         read(tok(7),*) iab(nfile)
         if(nfile.eq.1 .or. 
     $        sname(nfile)(:lnb(sname(nfile))).ne.
     $        sname(nfile-1)(:lnb(sname(nfile)))) then
            nsn = nsn + 1
            iepoch(nsn) = nfile
         end if
         nepoch(nsn) = nepoch(nsn) + 1
         nsntot = nsntot + 1
         read(tok(6),*) delta(nsn)
         if(nepoch(nsn).gt.maxep) then
            maxep = nepoch(nsn)
            iep = nsn
         end if

 4    continue
 5    close(1)

      write(6,*) 'nfile = ',nfile,' ; nsn = ',nsn
      write(6,*) 'Max epochs =',maxep,' for ',sname(iepoch(iep))
      write(6,*) 'Total spectra = ',nsntot

      if(maxep.gt.MAXEPOCH) then
         write(6,*) 'ERROR! Too many epochs for MAXEPOCH =',MAXEPOCH
         stop
      end if

c
c Construct a data table for each SN
c
      do k=1, nsn

* Zero out the tables
         do i=1, MAXKNOT
            do j=1, MAXEPOCH
               xknot(i,j) = 1 ! =1 so that alog10(xknot(i,j))=0.0 
               yknot(i,j) = 0
            end do
         end do

         do i=1, nw
            do j=1, MAXEPOCH
               flog(i,j) = 0
            end do
         end do

* Visit each data file
         mostknots = 0
         do j=1, nepoch(k)
            fspec = fname(iepoch(k) + j - 1)
            meconvert  = iab(iepoch(k) + j - 1)
            wmin  = wavemin(iepoch(k) + j - 1)
            wmax  = wavemax(iepoch(k) + j - 1)
            zfactor = redshift(iepoch(k) + j - 1) + 1
            open(unit=1,file=fspec,status='old',err=666)

* Read each epoch
            nwave = 0
            do 6 l=1, MAXWAVE
               read(1,'(a)',end=7) line
               if(index(line(nlb(line):nlb(line)+1),'#').gt.0) goto 6
               call parser(line,MAXTOK,ntok,ltok,tok)
               read(tok(1),*) wave(l)
               read(tok(2),*) flux(l)
               nwave = nwave + 1
               if(meconvert.eq.1) flux(l) = 10.**(-0.4*flux(l))
               if(wave(l).lt.wmin .or. wave(l).gt.wmax) flux(l) = 0
               wave(l) = wave(l) / zfactor
 6          continue
 7          close(1)

* Rebin onto the common log wavelength scale
            call rebin(nwave,wave,flux,MAXLOG,w0,dwlog,flog(1,j))

* Remove the mean with a spline fit
            call meanzero(MAXLOG,l1,l2,flog(1,j),fnorm(1,j),-1,
     $           nknot(j),xknot(1,j),yknot(1,j))
            mostknots = max(mostknots,nknot(j))

* Apodize the end
            call apodize(maxlog,l1,l2,fnorm(1,j),percent)

* Find an appropriate scale factor for the non-normalized spectrum
            nmean = 0
            fmean(j) = 0
            do i=1, nw
               if(flog(i,j).gt.0) then
                  nmean = nmean + 1
                  fmean(j) = fmean(j) + flog(i,j)
               end if
            end do
            fmean(j) = fmean(j) / max(1,nmean)
            do i=1, nw
               flog(i,j) = flog(i,j) / fmean(j)
            end do
         end do                 ! [do j=1, nepoch(k)]

* Create an output file name
         l = lnb(sname(iepoch(k)))
         fspec = sname(iepoch(k))(:l)//'.lnw'
         open(unit=2,file=fspec,status='unknown')

* Write sizes, ranges, delta, and SN type (including indices) into first line
         write(2,1001) nepoch(k),nw,w0,w1,mostknots,
     $        sname(iepoch(k)),delta(iepoch(k)),stype(iepoch(k)),
     $        itype(1,iepoch(k)),itype(2,iepoch(k))
 1001    format(2i5,2f10.2,i7,5x,a12,f7.2,2x,a10,2(i3))

* Write nknot and knots into next lines: log(w),log(f/fmean)
         write(2,1002) mostknots,
     $        (nknot(j),alog10(fmean(j)),j=1,nepoch(k))
 1002    format(i7,300(i3,f14.5))
         do i=1, mostknots
            write(2,1003) i,(alog10(xknot(i,j)),
     $           yknot(i,j)-alog10(fmean(j)),j=1,nepoch(k))
         end do
 1003    format(i7,300(f8.4,f9.4))

* Write ages into header line (1st element is age flag)
         write(2,1004) tflag(k),(epoch(iepoch(k)+i),i=0,nepoch(k)-1)
 1004    format(i8,300f9.3)

* Write wavelength and log(flux) in subseqent lines
         do i=1, nw
            wmean = 0.5*(wlog(i)+wlog(i+1))
            write(2,1005) wmean,(fnorm(i,j),j=1,nepoch(k))
         end do
 1005    format(f8.2,300f9.3)
         close(2)

         write(*,1006) fspec(:lnb(fspec)),nepoch(k)
 1006    format(1x,a,' processed, with ',i3,' epochs')

      end do                    ! [do k=1, nsn]

      goto 999

 666  write(6,*) 'ERROR! Could not open', fspec






************************************************************************
**                                                                    **
**                       End of program logwave                       **
**                                                                    **
************************************************************************

 999  stop
      end
