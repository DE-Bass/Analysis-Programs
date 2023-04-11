************************************************************************
**                                                                    **
**                    Plotting routines for SNID                      **
**                                                                    **
**                                                                    **
**                 Copyright (C) 2007 Stephane Blondin                **
**                                                                    **
************************************************************************

c Included in this file:
c
c subroutine initplotvar
c subroutine setplotvar
c subroutine initbutt
c subroutine setbutt
c subroutine getlabels
c subroutine getdlabels
c subroutine getpdat
c subroutine plotflux
c subroutine plotflat
c subroutine plotxcor
c subroutine plotpeaks
c subroutine plotzt
c subroutine plotfrac
c subroutine plotfracall
c subroutine getpggm



************************************************************************
* subroutine initplotvar -- Initialize plot variables
************************************************************************

      subroutine initplotvar

      implicit none

      include 'snid.inc'

      integer i,j

      ips = 0
      ips2 = 0
      iascii = 0
      jflux = 0
      jflat = 0

      do i=1, 2
         jpeaks(i) = 0
         jfracall(i) = 0
      end do

      do i=1, MAXTEMP
         do j=1, 2
            jcompflux(j,i) = 0
            jcompflat(j,i) = 0
            jcompxcor(j,i) = 0
         end do
      end do

      do i=1, NST
         do j=1, 2
            jfrac(j,i) = 0
         end do
      end do

      return
      end

************************************************************************
* subroutine setplotvar -- Set plot variables
************************************************************************

      subroutine setplotvar(nb)

      implicit none

      include 'snid.inc'

      integer nb ! button number

* Initialize variables
c disable "peaks" button if forcez > -1
      if(nb.eq.30.and.forcez.gt.-1) return
      iflux = 0
      iflat = 0
      ixcor = 0
      ipeaks = 0
      if(forcez.gt.-1) ipeaks = 3
      izt = 0
      iia = 0
      iib = 0
      iic = 0
      iii = 0
      inotsn = 0
      iall = 0
      if(nb.eq.38) then
         iflux = 4
      else if(nb.eq.39) then
         iflat = 4
      else if(nb.eq.40) then
         ixcor = 4
      else if(nb.eq.30) then
         ipeaks = 4
      else if(nb.eq.42) then
         izt = 4
      else if(nb.eq.32) then
         iia = 4
      else if(nb.eq.33) then
         iib = 4
      else if(nb.eq.34) then
         iic = 4
      else if(nb.eq.44) then
         iii = 4
      else if(nb.eq.45) then
         inotsn = 4
      else if(nb.eq.46) then
         iall = 4
      end if      

      return
      end



************************************************************************
* subroutine initbutt -- Setup static plot/button regions
************************************************************************

      subroutine initbutt

      implicit none

      include '../button/button.inc'

      MAX_XBUTT=12    ! number of buttons in the x direction
      MAX_YBUTT=4     ! number of buttons in the y direction
      PGSCF_BUTT=1    ! PGPLOT font type for text in buttons
      PGSCH_BUTT=1.   ! PGPLOT font size for text in buttons
      YTEXT_BUTT=0.35 ! relative y-position of the text baseline in buttons
      X1VPORT=0.1     ! x-coordinate of the left hand edge of the plot region
      X2VPORT=0.95    ! x-coordinate of the right hand edge of the plot region
      Y1VPORT=0.1     ! y-coordinate of the bottom edge of the plot region
      Y2VPORT=0.70    ! y-coordinate of the top edge of the plot region
      X3VPORT=0.05    ! x-coordinate of the left hand edge of the button region
      X4VPORT=0.95    ! x-coordinate of the right hand edge of the button region
      Y3VPORT=0.80    ! y-coordinate of the bottom edge of the button region
      Y4VPORT=1.00    ! y-coordinate of the top edge of the button region
      
      call buttsxb(MAX_XBUTT)
      call buttsyb(MAX_YBUTT)

* Plot buttons
       call button(7,'S u p e r N o v a    '//
     $     'I D e n t i f i c a t i o n (SNID)',3)
      call button(12,'QUIT',0)
      call button(14,'Template Navigation',1)
      call button(16,'No.',0)
      call button(25,'<<',0)
      call button(26,'<',0)
      call button(27,'>',0)
      call button(28,'>>',0)
      call button(37,'Show:',1)
      call button(18,'Other',1)
      call button(20,'Fraction',1)
      call button(24,'Output',1)
    
      return
      end



************************************************************************
* subroutine setbutt -- Set dynamic buttons
************************************************************************

      subroutine setbutt
      
      implicit none

      include 'snid.inc'
     
      call button(38,'flux',iflux)
      call button(39,'flat',iflat)
      call button(40,'xcor',ixcor)
      call button(30,'peaks',ipeaks)
      call button(42,'z/t',izt)
      call button(32,'Ia',iia)
      call button(33,'Ib',iib)
      call button(34,'Ic',iic)
      call button(44,'II',iii)
      call button(45,'NotSN',inotsn)
      call button(46,'All',iall)
      call button(36,'PS',ips)
      call button(48,'ASCII',iascii)

      return
      end



************************************************************************
* subroutine getlabels -- Get static labels
************************************************************************

      subroutine getlabels

      implicit none

      include 'snid.inc'

      integer l,l0

c Functions
      integer lnb
      
      write(labels(1),'(f6.1)') rlapmin
      write(labels(2),'(f6.3)') zfilter
      write(labels(3),'(f6.3)') zuser    
      
      return
      end



************************************************************************
* subroutine getdlabels -- Get dynamic labels
************************************************************************

      subroutine getdlabels(itn,sid,type,t,z,zerr,r)

      implicit none

      include 'snid.inc'

      integer itn       ! template number
      character*12 sid  ! name of SN template
      character*10 type ! type of SN template
      real t,z,zerr,r   ! age, redshift, redshift error, r-value

c Functions
      integer lnb,nlb

      write(dlabels(1),*) itn
      dlabels(2) = sid
      dlabels(3) = type
      write(dlabels(4),*) nint(t)
      if(t.ge.0) dlabels(4)='+'//
     $     dlabels(4)(nlb(dlabels(4)):lnb(dlabels(4)))
      write(dlabels(5),'(f6.3)') z
      write(dlabels(6),'(f5.3)') zerr
      if(r.gt.MAXR) r = MAXR
      write(dlabels(7),'(f6.1)') r

* template number label for PS/ASCII output
      if(itn.lt.10) then
         write(dlabels(8)(1:3),'(a3)') '000'
         write(dlabels(8)(4:4),'(i1)') itn
      else if(itn.lt.100) then
         write(dlabels(8)(1:2),'(a2)') '00'
         write(dlabels(8)(3:4),'(i2)') itn
      else if(itn.lt.1000) then
         write(dlabels(8)(1:1),'(a1)') '0'
         write(dlabels(8)(2:4),'(i3)') itn
      else if(itn.lt.10000) then
         write(dlabels(8)(1:4),'(i4)') itn
      else
         write(dlabels(8)(1:4),'(a)') 'xxxx'         
      end if

      return
      end



************************************************************************
* subroutine getpdat -- Create pdat array for plotting
************************************************************************

      subroutine getpdat(pdat,fnorm,dft,temp,tft,cfn,cft,
     $        ctr,hgt,drms,trms,arms,srms,r,redshift,wid,z,
     $        nknod,xknod,yknod,nknot,xknot,yknot)

c See below for what data is stored for different values of pdat(1,X)
c
c X = 
c  1 : Object flux (mean)                  6 : Template flux (mean)
c  2 : Object flux                         7 : Template flux
c  3 : Object flux (filtered)              8 : Template flux (filtered)
c  4 : Object flux (trimmed)               9 : Template flux (trimmed)
c  5 : Object flux (trim+filt)            10 : Template flux (trim+filt)
c 15 : Object flux (mean), deredshifted   16 : Template flux (mean), redshifted
c 11 : Correlation function               13 : Phase band
c 12 : Correlation function (trimmed)     14 : Phase band (trimmed)
c 17 : Peak trials: lap
c 18 : Peak trials: r value
c 19 : Peak trials: redshift
c 20 : Peak trials: trim redshift

      implicit none

      include 'snid.inc'

      real pdat(MAXPPT,MAXPLOT)                   ! data array for plotting (cf. snid.pl)
      real fnorm(MAXLOG)                          ! normalized input spectrum
      complex dft(MAXLOG),tft(MAXLOG),cft(MAXLOG) ! data/template/correlation function fourier transforms
      real temp(MAXLOG)                           ! array containing the template spectrum [BEWARE! Different declaration from snid.f]
      real cfn(2*MAXLOG)                          ! (filtered) correlation function 
      real ctr,hgt,wid,r                          ! center,height,width and r-value of best peak
      real z                                      ! template redshift
      real drms,trms                              ! data/template RMS
      real arms,srms                              ! (anti)symmetric RMS of correlation function
      real redshift                               ! correlation redshift
      integer nknot                               ! number of knot points for the template spectrum [BEWARE! Different declaration from snid.f]
      real xknot(MAXKNOT),yknot(MAXKNOT)          ! knot point (w,f) coordinates for the template spectrum [BEWARE! Different declaration from snid.f]
      integer lbest                               ! log lambda bin number of best correlation peak
      real lap                                    ! data/template overlap at correlation redshift
      real ftrim(MAXLOG),ttrim(MAXLOG)            ! redshift-trimmed inpute/template spectra
     
* Compute the correlation function ab initio
      call docorr(fnorm,dft,temp,tft,cfn,cft,ctr,hgt,
     $     drms,trms,arms,srms,r,redshift,wid)

* Save the spectra
      call loadspec(fnorm,dft,temp,tft,pdat,2,3,7,8)

* Save the smoothed spline fits
      lbest = nint(alog(z+1)/dwlog)
      call splinedex(nw,pdat(1,1),0,nknod,xknod,yknod)
      call splinedex(nw,pdat(1,15),lbest,nknod,xknod,yknod)
      call splinedex(nw,pdat(1,6),0,nknot,xknot,yknot)
      call splinedex(nw,pdat(1,16),-lbest,nknot,xknot,yknot)

* Try again with the ends trimmed a bit
      ctr = lbest
      call overlap(nw,ctr,lap,temp,ttrim,fnorm,ftrim,percent)
      call rrvector(nw,cfn,pdat(1,11))
      call rvscale(nw,pdat(1,11),1/(drms*trms*nw))
      call phaseband(nw,ctr,cft,pdat(1,13))

* Compute the trimmed correlation function
      call docorr(ftrim,dft,ttrim,tft,cfn,cft,ctr,hgt,
     $     drms,trms,arms,srms,r,redshift,wid)

* Save the raw and smoothed spectra
      call loadspec(ftrim,dft,ttrim,tft,pdat,4,5,9,10)
      call rrvector(nw,cfn,pdat(1,12))
      call rvscale(nw,pdat(1,12),1/(drms*trms*nw))
      call phaseband(nw,ctr,cft,pdat(1,14))

      return
      end



************************************************************************
* subroutine plotflux -- Plot fluxed spectra
************************************************************************

      subroutine plotflux(itn,z,dfl,dmean,dmeanz,tfl)

      implicit none

      include 'snid.inc'
      include '../button/button.inc'

      integer i
      integer itn                                               ! template number
      real z                                                    ! template redshift
      real dfl(MAXPPT),dmean(MAXPPT),tfl(MAXPPT),dmeanz(MAXPPT) ! data flux and mean; template flux; de-redshifted data mean 
      character*200 psfile,asciifile                            ! name of file for PS/ASCII output
      real xmin,xmax,ymin,ymax,dy                               ! x- and y-ranges for plot
      real dflux(MAXWAVE),tflux(MAXWAVE)                        ! data/template flux with mean included
      real zwave(MAXWAVE)                                       ! redshifted wavelength axis for template
      character*100 tit,xtit,ytit                               ! (x-/y-)title for plot box

c Functions
      integer lnb,pgopen

* Erase previous plot      
      call rpgerasw(0.,1.,0.,Y3VPORT) 

* X-range
      xmin = wmin/(1.+z)
      xmax = wmax/(1.+z)

* Data/template flux normalized to data flux
      do i=1, nw                
         dflux(i) = dfl(i) + 1.
         dflux(i) = dflux(i) * dmean(i)
         tflux(i) = tfl(i) + 1.
         tflux(i) = tflux(i) * dmeanz(i)
      end do

* Y-range and wavelength axis
      ymin = 1e6                
      ymax = -1
      do i=1, nw
         wave(i) = w0 * exp(float(i-1)*dwlog)
         zwave(i)= wave(i) * (1.+z) 
         if(wave(i).ge.wmin.and.wave(i).le.wmax) then
            ymin = min(ymin,dflux(i))
            ymax = max(ymax,dflux(i))
         end if
      end do
      dy = ymax - ymin
      ymin = ymin - 0.025*dy
      ymax = ymax + 0.2*dy

* ASCII output
      if(iascii.gt.0) then
c input spectrum
         asciifile = froot(:lnb(froot))//'_snidflux.dat'
         if(jflux.eq.0) then
            open(unit=1,file=asciifile,status='replace')
            write(1,'(a)') '#wavelength[A] flux[arbitraty]'
            do i=1, nw
               if(wave(i).ge.wmin.and.wave(i).le.wmax) 
     $              write(1,*) wave(i),dflux(i)
            end do
            close(1)
            write(6,'(2a)') ' Created ASCII file: ',
     $           asciifile(:lnb(asciifile))
            jflux = 1
         end if
c template spectrum
         asciifile = froot(:lnb(froot))//'_comp'//
     $        dlabels(8)(:lnb(dlabels(8)))//
     $        '_snidflux.dat'
         if(jcompflux(1,itn).eq.0) then
            open(unit=1,file=asciifile,status='replace')
            write(1,'(a)') '# Template No.'//
     $           dlabels(1)(:lnb(dlabels(1)))//': '//
     $           dlabels(2)(:lnb(dlabels(2)))//' ('//
     $           dlabels(3)(:lnb(dlabels(3)))//'; '//
     $           dlabels(4)(:lnb(dlabels(4)))//')'//' ; z='//
     $           dlabels(5)(:lnb(dlabels(5)))//' +/- '//
     $           dlabels(6)(:lnb(dlabels(6)))
            write(1,'(a)') '#redshifted_wavelength[A] flux[arbitraty]'
            do i=1, nw
               if(zwave(i).ge.wmin.and.zwave(i).le.wmax) 
     $              write(1,*) zwave(i),tflux(i)
            end do
            close(1)
            write(6,'(2a)') ' Created ASCII file: ',
     $           asciifile(:lnb(asciifile))
            jcompflux(1,itn) = 1
         else
            write(6,'(2a)') ' Already created ASCII file: ',
     $           asciifile(:lnb(asciifile))            
         end if
      end if
      iascii = 0
      call button(48,'ASCII',iascii)

* Plot box; open PS device
      xtit='Observed Wavelength [\\A]'
      ytit='Relative Flux'
      tit='Rest Wavelength [\\A]'
      if(ips.gt.0) then
         psfile = froot(:lnb(froot))//'_comp'//
     $        dlabels(8)(:lnb(dlabels(8)))//
     $        '_snidflux.ps'
         if(jcompflux(2,itn).eq.0) then
            istat = pgopen(psfile//'/cps')
            if (istat.le.0) write(6,'(a)') 
     $           ' PGOPEN: Could not open PS device'
            call pgslw(3)
         else
            write(6,'(2a)') ' Already created PS file: ',
     $           psfile(:lnb(psfile))         
         end if
      end if
      call rpgenv(xmin,xmax,ymin,ymax,0,6)
      call rpgenv(wmin,wmax,ymin,ymax,0,3)
      call pglabel(xtit,ytit,tit)

* Plot spectra and labels
      call pgline(nw,wave,dflux)
      call pgmtxt('T',-1.75,.025,0,'Input: '//fdata)
      call pgsave
      call pgsci(2) 
      call pgmtxt('T',-3.25,.025,0,'No.'//
     $     dlabels(1)(:lnb(dlabels(1)))//': '//
     $     dlabels(2)(:lnb(dlabels(2)))//' ('//
     $     dlabels(3)(:lnb(dlabels(3)))//'; '//
     $     dlabels(4)(:lnb(dlabels(4)))//')'//' ; z='//
     $     dlabels(5)(:lnb(dlabels(5)))//'\\(2233)'//
     $     dlabels(6)(:lnb(dlabels(6))))
      call pgslw(5)          
      call pgline(nw,zwave,tflux)
      call pgunsa

* Close PS device
      if(ips.gt.0) then
         if(jcompflux(2,itn).eq.0) then
            call pgslw(1)
            call pgclos
            write(6,'(2a)') ' Created PS file: ',psfile(:lnb(psfile))
         end if
         ips = 0
         ips2 = 1
         jcompflux(2,itn) = 1
      end if
      
      return
      end



************************************************************************
* subroutine plotflat -- Plot flattened spectra
************************************************************************

      subroutine plotflat(itn,z,dfl,tfl)

      implicit none

      include 'snid.inc'
      include '../button/button.inc'
      
      integer i
      integer itn                    ! template number
      real z                         ! template redshift
      real dfl(MAXPPT),tfl(MAXPPT)   ! data and template flux
      character*200 psfile,asciifile ! name of file for PS/ASCII output
      real xmin,xmax,ymin,ymax,dy    ! x- and y-ranges for plot
      real zwave(MAXWAVE)            ! redshifted wavelength acis for template
      character*100 tit,xtit,ytit    ! (x-/y-)title for plot box

c Functions
      integer lnb,pgopen

* Erase previous plot      
      call rpgerasw(0.,1.,0.,Y3VPORT) 

* X-range
      xmin = wmin/(1.+z)
      xmax = wmax/(1.+z)

* Y-range and wavelength axis
      ymin = 1e6
      ymax = -1
      do i=1, nw
         wave(i) = w0 * exp(float(i-1)*dwlog)
         zwave(i) = wave(i) * (1.+z)
         if(wave(i).ge.wmin.and.wave(i).le.wmax) then
            ymin = min(ymin,dfl(i))
            ymax = max(ymax,dfl(i))
         end if
      end do
      dy = ymax - ymin
      ymin = ymin - 0.025*dy
      ymax = ymax + 0.2*dy

* ASCII output
      if(iascii.gt.0) then
c input spectrum
         asciifile = froot(:lnb(froot))//'_snidflat.dat'
         if(jflat.eq.0) then
            open(unit=1,file=asciifile,status='replace')
            write(1,'(a)') '#wavelength[A] flux[arbitraty]'
            do i=1, nw
               if(wave(i).ge.wmin.and.wave(i).le.wmax) 
     $              write(1,*) wave(i),dfl(i)
            end do
            close(1)
            write(6,'(2a)') ' Created ASCII file: ',
     $           asciifile(:lnb(asciifile))
            jflat = 1
         end if
c template spectrum
         asciifile = froot(:lnb(froot))//'_comp'//
     $        dlabels(8)(:lnb(dlabels(8)))//
     $        '_snidflat.dat'
         if(jcompflat(1,itn).eq.0) then
            open(unit=1,file=asciifile,status='replace')
            write(1,'(a)') '# Template No.'//
     $           dlabels(1)(:lnb(dlabels(1)))//': '//
     $           dlabels(2)(:lnb(dlabels(2)))//' ('//
     $           dlabels(3)(:lnb(dlabels(3)))//'; '//
     $           dlabels(4)(:lnb(dlabels(4)))//')'//' ; z='//
     $           dlabels(5)(:lnb(dlabels(5)))//' +/- '//
     $           dlabels(6)(:lnb(dlabels(6)))
            write(1,'(a)') '#redshifted_wavelength[A] flux[arbitraty]'
            do i=1, nw
               if(zwave(i).ge.wmin.and.zwave(i).le.wmax) 
     $              write(1,*) zwave(i),tfl(i)
            end do
            close(1)
            write(6,'(2a)') ' Created ASCII file: ',
     $           asciifile(:lnb(asciifile))
            jcompflat(1,itn) = 1
         else
            write(6,'(2a)') ' Already created ASCII file: ',
     $           asciifile(:lnb(asciifile))            
         end if
      end if
      iascii = 0
      call button(48,'ASCII',iascii)

* Plot box; open PS device
      xtit='Observed Wavelength [\\A]'
      ytit='Flattened Flux'
      tit='Rest Wavelength [\\A]'
      if(ips.gt.0) then
         psfile = froot(:lnb(froot))//'_comp'//
     $        dlabels(8)(:lnb(dlabels(8)))//
     $        '_snidflat.ps'
         if(jcompflat(2,itn).eq.0) then
            istat = pgopen(psfile//'/cps')
            if (istat.le.0) write(6,'(a)') 
     $           ' PGOPEN: Could not open PS device'
            call pgslw(3)
         else
            write(6,'(2a)') ' Already created PS file: ',
     $           psfile(:lnb(psfile))            
         end if
      end if
      call rpgenv(xmin,xmax,ymin,ymax,0,6)
      call rpgenv(wmin,wmax,ymin,ymax,0,3)
      call pglabel(xtit,ytit,tit)

* Plot spectra and labels
      call pgline(nw,wave,dfl)
      call pgmtxt('T',-1.75,.025,0,'Input: '//fdata)
      call pgsave
      call pgsci(2) 
      call pgmtxt('T',-3.25,.025,0,'No.'//
     $     dlabels(1)(:lnb(dlabels(1)))//': '//
     $     dlabels(2)(:lnb(dlabels(2)))//' ('//
     $     dlabels(3)(:lnb(dlabels(3)))//'; '//
     $     dlabels(4)(:lnb(dlabels(4)))//')'//' ; z='//
     $     dlabels(5)(:lnb(dlabels(5)))//'\\(2233)'//
     $     dlabels(6)(:lnb(dlabels(6))))
      call pgslw(5)          
      call pgline(nw,zwave,tfl)
      call pgunsa

* Close PS device
      if(ips.gt.0) then
         if(jcompflat(2,itn).eq.0) then
            call pgslw(1)
            call pgclos
            write(6,'(2a)') ' Created PS file: ',psfile(:lnb(psfile))
         end if
         ips = 0
         ips2 = 1
         jcompflat(2,itn) = 1
      end if

      return
      end



************************************************************************
* subroutine plotxcor -- Plot correlation functions
************************************************************************

      subroutine plotxcor(itn,npeak,zpeak,cfn,cfntrim)

      implicit none

      include 'snid.inc'
      include '../button/button.inc'
      
      integer i
      integer itn                               ! template number
      integer npeak                             ! total number of correlation peaks
      real zpeak(MAXPPT)                        ! z of peak
      real cfn(2*MAXLOG),cfntrim(2*MAXLOG)      ! correlation function and trimmed version
      real xcor(3*MAXLOG),xcortrim(3*MAXLOG)    ! same as cfn(trim) but extended as needed in redshift
      integer nz                                ! number of redshift bins between xmin and xmax
      character*200 psfile,asciifile,asciifile2 ! name of file for PS/ASCII output
      real xmin,xmax,ymin,ymax,dy               ! x- and y-ranges for plot
      real zaxis(MAXWAVE)                       ! redshifted axis
      character*100 tit,xtit,ytit               ! (x-/y-)title for plot box

c Functions
      integer lnb,pgopen

* Erase previous plot      
      call rpgerasw(0.,1.,0.,Y3VPORT) 

* X-range
      xmin = zmin - 2*EPSZ
      xmax = zmax + 2*EPSZ
      lz1 = nint(alog(xmin+1) / dwlog) 
      lz2 = nint(alog(xmax+1) / dwlog) 

* Y-range and redshift axis
      nz = 0
      ymin = 0.
      ymax = 1.
      do i=lz1, lz2
         nz = nz + 1
         if(i.le.0) then
            xcor(nz) = cfn(nw-mod(-i,nw))
            xcortrim(nz) = cfntrim(nw-mod(-i,nw))
         else if(i.le.nw) then
            xcor(nz) = cfn(i)
            xcortrim(nz) = cfntrim(i)            
         else
            xcor(nz) = cfn(mod(i,nw))
            xcortrim(nz) = cfntrim(mod(i,nw))
         end if
         zaxis(nz) = exp((i-1)*dwlog) - 1.
         if(zaxis(nz).ge.xmin.and.zaxis(nz).le.xmax) then
            ymin = min(ymin,xcor(nz))
            ymin = min(ymin,xcortrim(nz))
         end if
      end do
      dy = ymax - ymin
      ymin = ymin - 0.05*dy

* ASCII output
      if(iascii.gt.0) then
         asciifile = froot(:lnb(froot))//'_comp'//
     $        dlabels(8)(:lnb(dlabels(8)))//
     $        '_snidxcor.dat'
         asciifile2 = froot(:lnb(froot))//'_comp'//
     $        dlabels(8)(:lnb(dlabels(8)))//
     $        '_snidxcor_trim.dat'
         if(jcompxcor(1,itn).eq.0) then
c correlation function
            open(unit=1,file=asciifile,status='replace')
            write(1,'(a)') '# Template No.'//
     $           dlabels(1)(:lnb(dlabels(1)))//': '//
     $           froot(:lnb(froot))//' x '//
     $           dlabels(2)(:lnb(dlabels(2)))//' ('//
     $           dlabels(3)(:lnb(dlabels(3)))//'; '//
     $           dlabels(4)(:lnb(dlabels(4)))//')'
            write(1,'(a)') '#redshift correlation'
            do i=1, nz
               write(1,*) zaxis(i),xcor(i)
            end do
            close(1)
c trimmed correlation function
            open(unit=1,file=asciifile2,status='replace')
            write(1,'(a)') '# Template No.'//
     $           dlabels(1)(:lnb(dlabels(1)))//': '//
     $           froot(:lnb(froot))//' x '//
     $           dlabels(2)(:lnb(dlabels(2)))//' ('//
     $           dlabels(3)(:lnb(dlabels(3)))//'; '//
     $           dlabels(4)(:lnb(dlabels(4)))//')'
            write(1,'(a)') '#redshift correlation'
            do i=1, nz
               write(1,*) zaxis(i),xcortrim(i)
            end do
            close(1)
c
            write(6,'(2a,2x,a)') ' Created ASCII files: ',
     $           asciifile(:lnb(asciifile)),
     $           asciifile2(:lnb(asciifile2))
            jcompxcor(1,itn) = 1            
         else
            write(6,'(2a,2x,a)') ' Already created ASCII files: ',
     $           asciifile(:lnb(asciifile)),
     $           asciifile2(:lnb(asciifile2))
         end if
      end if
      iascii = 0
      call button(48,'ASCII',iascii)

* Plot box; open PS device
      xtit='Redshift'
      ytit='Normalized correlation amplitude'
      tit='Correlation function'
      if(ips.gt.0) then
         psfile = froot(:lnb(froot))//'_comp'//
     $        dlabels(8)(:lnb(dlabels(8)))//
     $        '_snidxcor.ps'
         if(jcompxcor(2,itn).eq.0) then
            istat = pgopen(psfile//'/cps')
            if (istat.le.0) write(6,'(a)') 
     $           ' PGOPEN: Could not open PS device'
            call pgslw(3)
         else
            write(6,'(2a)') ' Already created PS file: ',
     $           psfile(:lnb(psfile))            
         end if
      end if
      call rpgenv(xmin,xmax,ymin,ymax,0,0)
      call pglabel(xtit,ytit,tit)

* Plot correlation functions and labels
      call pgline(nz,zaxis,xcor)
      call pgmtxt('T',-1.75,.025,0,'No.'//
     $     dlabels(1)(:lnb(dlabels(1)))//': '//
     $     froot(:lnb(froot))//' x '//
     $     dlabels(2)(:lnb(dlabels(2)))//' ('//
     $     dlabels(3)(:lnb(dlabels(3)))//'; '//
     $     dlabels(4)(:lnb(dlabels(4)))//')')
      call pgsave
      call pgsci(2) 
      call pgmtxt('T',-3.25,.025,0,'Trimmed at z='//
     $     labels(3)(:lnb(labels(3)))//' ; r='//
     $     dlabels(7)(:lnb(dlabels(7))))
      call pgslw(3)
      call pgsls(2)          
      call pgline(nz,zaxis,xcortrim)
      call pgunsa            
      
* Close PS device
      if(ips.gt.0) then
         if(jcompxcor(2,itn).eq.0) then
            call pgslw(1)
            call pgclos
            write(6,'(2a)') ' Created PS file: ',psfile(:lnb(psfile))
         end if
         ips = 0
         ips2 = 1
         jcompxcor(2,itn) = 1
      end if

      return
      end



************************************************************************
* subroutine plotpeaks -- Plot correlation peaks
************************************************************************

      subroutine plotpeaks(npeak,zpeak,rpeak,lpeak)

      implicit none

      include 'snid.inc'
      include '../button/button.inc'

      integer i
      integer npeak                                  ! total number of correlation peaks
      real rpeak(MAXPPT),lpeak(MAXPPT),zpeak(MAXPPT) ! r,lap,z of peak
      real rlap(MAXPPT)                              ! rlap values for all spectra
      character*200 psfile,asciifile                 ! name of file for PS/ASCII output
      real xmin,xmax,ymin,ymax,dx,dy                 ! x- and y-ranges for plot
      real xline(2),yline(2)                         ! for plotting straight lines
      character*100 tit,xtit,ytit                    ! (x-/y-)title for plot box

c Functions
      integer lnb,pgopen
      
* Erase previous plot      
      call rpgerasw(0.,1.,0.,Y3VPORT) 

* X-range
      xmin = 1e6
      xmax = -1
      do i=1, npeak
         xmin = min(xmin,zpeak(i))
         if(xmin.lt.zmin) xmin = zmin         
         xmax = max(xmax,zpeak(i))
         if(xmax.gt.zmax) xmax = zmax
      end do
      dx = xmax - xmin
      if(dx.eq.0) dx = 1.
      xmin = xmin - 0.03*dx
      xmax = xmax + 0.03*dx
      
* Y-range
      ymin = 0.
      ymax = 0.
      do i=1, npeak
         rlap(i) = rpeak(i) * lpeak(i)
         ymax = max(ymax,rlap(i))
      end do
      dy = ymax - ymin
      if(dy.eq.0) dy = 1.
      ymax = real(nint(ymax)) + 0.15*dy
      
* ASCII output
      if(iascii.gt.0) then
         asciifile = froot(:lnb(froot))//'_snidpeaks.dat'
         if(jpeaks(1).eq.0) then
            open(unit=1,file=asciifile,status='replace')
            write(1,'(a)') '# Input: '//fdata(:lnb(fdata))
            write(1,'(a)') '#redshift rlap'
            do i=1, npeak
               if(zpeak(i).ge.xmin.and.zpeak(i).le.xmax) 
     $              write(1,*) zpeak(i),rlap(i)
            end do
            close(1)
            write(6,'(2a)') ' Created ASCII file: ',
     $           asciifile(:lnb(asciifile))
            jpeaks(1) = 1            
         else
            write(6,'(2a)') ' Already created ASCII file: ',
     $           asciifile(:lnb(asciifile))
         end if
      end if
      iascii = 0
      call button(48,'ASCII',iascii)

* Plot box; open PS device
      xtit='Redshift'
      ytit='rlap'
      tit='Initial correlation peaks for all template spectra'
      if(ips.gt.0) then
         psfile = froot(:lnb(froot))//'_snidpeaks.ps'
         if(jpeaks(2).eq.0) then
            istat = pgopen(psfile//'/cps')
            if (istat.le.0) write(6,'(a)') 
     $           ' PGOPEN: Could not open PS device'
            call pgslw(3)
         else
            write(6,'(2a)') ' Already created PS file: ',
     $           psfile(:lnb(psfile))            
         end if
      end if
      call rpgenv(xmin,xmax,ymin,ymax,0,0)
      call pglabel(xtit,ytit,tit)

* Plot correlation peaks and labels
      call pgpt(npeak,zpeak,rlap,0)
      call pgmtxt('T',-1.75,.025,0,'Input: '//fdata)
      call pgsave
      call pgsci(2)
      call pgsls(2)
      call pgslw(3)          
      xline(1) = xmin
      xline(2) = xmax
      yline(1) = 4.
      yline(2) = 4.
      call pgline(2,xline,yline)
      call pgunsa
      
* Close PS device
      if(ips.gt.0) then
         if(jpeaks(2).eq.0) then
            call pgslw(1)
            call pgclos
            write(6,'(2a)') ' Created PS file: ',psfile(:lnb(psfile))
         end if
         ips = 0
         ips2 = 1
         jpeaks(2) = 1
      end if

      return
      end



************************************************************************
* subroutine plotzt -- Plot redshift/age scatter plot
************************************************************************

      subroutine plotzt(ngood,nbad,zgood,tgood,zbad,tbad,zbest,tbest,
     $     zmed,agem)

      implicit none

      include 'snid.inc'
      include '../button/button.inc'

      integer i
      integer ngood,nbad                ! number of "good" (rlap>rlapmin AND |z-zmed|<zfilter) and "bad" (rlap>rlapmin AND |z-zmed|>zfilter) correlations
      real zgood(MAXTEMP),zbad(MAXTEMP) ! good and bad redshifts
      real tgood(MAXTEMP),tbad(MAXTEMP) ! good and bad ages
      real zbest,tbest                  ! redshift/age of best-match template
      real zmed,agem                    ! median redshift and age
      character*200 psfile,asciifile    ! name of file for PS/ASCII output
      real xmin,xmax,ymin,ymax,dx,dy    ! x- and y-ranges for plot
      real xline(2),yline(2)            ! for plotting straight lines
      character*100 tit,xtit,ytit       ! (x-/y-)title for plot box

c Functions
      integer lnb,pgopen
      
* Erase previous plot      
      call rpgerasw(0.,1.,0.,Y3VPORT) 

* X- and Y-ranges
      xmin = 1e6
      xmax = -1
      ymin = 1e6
      ymax = -1e6
      do i=1, ngood
         xmin = min(xmin,zgood(i))
         xmax = max(xmax,zgood(i))
         ymin = min(ymin,tgood(i))
         ymax = max(ymax,tgood(i))
      end do
      do i=1, nbad
         xmin = min(xmin,zbad(i))
         xmax = max(xmax,zbad(i))
         ymin = min(ymin,tbad(i))
         ymax = max(ymax,tbad(i))
      end do
      if(xmin.eq.1e6)  xmin = zmin
      if(xmax.eq.-1)   xmax = zmax
      if(ymin.eq.1e6)  ymin = agemin
      if(ymax.eq.-1e6) ymax = agemax
      dx = xmax - xmin
      if(dx.eq.0) dx = 1.
      dy = ymax - ymin
      if(dy.eq.0) dy = 1.
      xmin = xmin - 0.03*dx
      xmax = xmax + 0.03*dx
      ymin = ymin - 0.1*dy
      ymax = ymax + 0.15*dy

* ASCII output
      if(iascii.gt.0) then
         asciifile = froot(:lnb(froot))//'_snidzt.dat'
         if(jzt(1).eq.0) then
            open(unit=1,file=asciifile,status='replace')
            write(1,'(a)') '# Input: '//fdata(:lnb(fdata))
            if(ngood.gt.0) then 
               write(1,'(a,f6.3,2x,f6.1)') '# Median z/t: ',zmed,agem
            else
               write(1,'(a,f6.3,2x,f6.1)') '# Median z/t: ',zmed,-99.9
            end if
            write(1,'(a)') '# zfilter: '//labels(2)(:lnb(labels(2)))
            write(1,'(a)') '#redshift age grade'
            do i=1, ngood
               write(1,'(f7.4,2x,f6.1,2x,a4)') zgood(i),tgood(i),'good'
            end do
            do i=1, nbad
               write(1,'(f7.4,2x,f6.1,2x,a3)') zbad(i),tbad(i),'bad'
            end do
            close(1)
            write(6,'(2a)') ' Created ASCII file: ',
     $           asciifile(:lnb(asciifile))
            jzt(1) = 1            
         else
            write(6,'(2a)') ' Already created ASCII file: ',
     $           asciifile(:lnb(asciifile))
         end if
      end if
      iascii = 0
      call button(48,'ASCII',iascii)

* Plot box; open PS device
      xtit='Redshift'
      ytit='Age'
      tit='Redshift vs. Age for correlations with rlap \\(2244) '//
     $     labels(1)(:lnb(labels(1)))
      if(ips.gt.0) then
         psfile = froot(:lnb(froot))//'_snidzt.ps'
         if(jzt(2).eq.0) then
            istat = pgopen(psfile//'/cps')
            if (istat.le.0) write(6,'(a)') 
     $           ' PGOPEN: Could not open PS device'
            call pgslw(3)
         else
            write(6,'(2a)') ' Already created PS file: ',
     $           psfile(:lnb(psfile))            
         end if
      end if
      call rpgenv(xmin,xmax,ymin,ymax,0,0)
      call pglabel(xtit,ytit,tit)

* Plot z/t points and labels
      call pgpt(ngood,zgood,tgood,0)
      call pgmtxt('T',-1.75,.025,0,'Input: '//fdata)
      call pgsave
      xline(1) = xmin
      xline(2) = xmax
      yline(1) = agem
      yline(2) = agem
      call pgsls(2)
      if(ngood.gt.0) then 
         call pgline(2,xline,yline)
      end if
      xline(1) = zmed
      xline(2) = zmed
      yline(1) = ymin
      yline(2) = ymax
      call pgline(2,xline,yline)         
      if(ngood.gt.0) then 
         call pgsci(3)
         call pgsch(2.25)
         call pgpt1(zbest,tbest,17)
         call pgsch(1.0)
         call pgmtxt('T',-3.25,.75,0,char(17)//' best-match')
      end if
      call pgunsa
      call pgmtxt('T',-4.75,.75,0,'\\(0806) \\(0806) \\(0806) '//
     $     'median')
      if(nbad.gt.0) then
         call pgsci(2)
         call pgmtxt('T',-6.25,.75,0,'|z-'//
     $        labels(3)(:lnb(labels(3)))//'| \\(2244) '//
     $        labels(2)(:lnb(labels(2))))
         call pgpt(nbad,zbad,tbad,0)
      end if
      
* Close PS device
      if(ips.gt.0) then
         if(jzt(2).eq.0) then
            call pgslw(1)
            call pgclos
            write(6,'(2a)') ' Created PS file: ',psfile(:lnb(psfile))
         end if
         ips = 0
         ips2 = 1
         jzt(2) = 1
      end if

      return
      end



************************************************************************
* subroutine plotfrac -- Plot (sub)type fractions
************************************************************************

      subroutine plotfrac(it,nrl,rlarr,frl,frlerr)

      implicit none

      include 'snid.inc'
      include '../button/button.inc'      

      integer i
      integer nrl                        ! number of elements in rlap array (rlarr)
      real rlarr(MAXRLAP)                ! rlap array (rlaparr[0:rlapmax])
      real frl(MAXRLAP,NT,NST)           ! absolute fraction of templates with rlap >= rlap_i for each (sub)type
      real frlerr(MAXRLAP,NT,NST)        ! poisson error in absolute fraction of templates with rlap >= rlap_i for each (sub)type
      real x(MAXRLAP),y(MAXRLAP)         ! vectors for plotting
      real yerr1(MAXRLAP),yerr2(MAXRLAP) ! vectors for plotting
      integer nsub                       ! number of subtypes corresponding to type it
      character*200 psfile,asciifile     ! name of file for PS/ASCII output
      real xmin,xmax,ymin,ymax           ! x- and y-ranges for plot
      real xline(2),yline(2)             ! for plotting straight lines
      character*100 tit,xtit,ytit        ! (x-/y-)title for plot box
      integer pggmidx(NST)               ! PGPLOT graph marker indices

c Functions
      integer lnb,pgopen

* Erase previous plot      
      call rpgerasw(0.,1.,0.,Y3VPORT) 

* X- and Y-ranges
      xmin = 0.
      if(nrl.gt.MAXRLAP) nrl = MAXRLAP
      xmax = rlarr(nrl)
      if(xmax.eq.0) xmax=5.
      ymin = 0.
      ymax = 1.

* Number of subtypes for type index it
      nsub = 0
      do ist=1, NST
         if(index(typename(it,ist),'typename').gt.0) goto 10
         nsub = nsub + 1
      end do
 10   continue

* ASCII output
      if(iascii.gt.0) then
         asciifile = froot(:lnb(froot))//'_snidfrac_'//
     $        typename(it,1)(:lnb(typename(it,1)))//'.dat'
         if(jfrac(1,it).eq.0) then
            open(unit=1,file=asciifile,status='replace')
            write(1,'(a)') '# Input: '//fdata(:lnb(fdata))
            write(1,'(a)') '# Type: '//
     $           typename(it,1)(:lnb(typename(it,1)))       
            write(1,'(a)') '#rlap fraction fraction_error type'
            do ist=1, nsub
               do i=1, nrl
                  write(1,'(f5.1,2x,2(f7.4,2x),a10)') rlarr(i),
     $                 frl(i,it,ist),frlerr(i,it,ist),typename(it,ist)
               end do
            end do
            close(1)
            write(6,'(2a)') ' Created ASCII file: ',
     $           asciifile(:lnb(asciifile))
            jfrac(1,it) = 1            
         else
            write(6,'(2a)') ' Already created ASCII file: ',
     $           asciifile(:lnb(asciifile))
         end if
      end if
      iascii = 0
      call button(48,'ASCII',iascii)

* Plot box; open PS device
      xtit='rlap [rlapmin = '//labels(1)(:lnb(labels(1)))//']'
      ytit='Template fraction \\(2244) rlap'
      tit='Template fraction for type: '//
     $     typename(it,1)(:lnb(typename(it,1)))
      if(ips.gt.0) then
         psfile = froot(:lnb(froot))//'_snidfrac_'//
     $        typename(it,1)(:lnb(typename(it,1)))//'.ps'
         if(jfrac(2,it).eq.0) then
            istat = pgopen(psfile//'/cps')
            if (istat.le.0) write(6,'(a)') 
     $           ' PGOPEN: Could not open PS device'
            call pgslw(3)
         else
            write(6,'(2a)') ' Already created PS file: ',
     $           psfile(:lnb(psfile))            
         end if
      end if
      call rpgenv(xmin,xmax,ymin,ymax,0,0)
      call pglabel(xtit,ytit,tit)

* Plot fraction curves and labels
      call pgsave
      call getpggm(pggmidx)
      do ist=1, nsub
         do i=1, nrl
            x(i) = rlarr(i)
            y(i) = frl(i,it,ist)
            yerr1(i) = frl(i,it,ist)+frlerr(i,it,ist)
            yerr2(i) = frl(i,it,ist)-frlerr(i,it,ist)
         end do
         call pgsci(ist)
         call pgslw(3)
         call pgline(nrl,x,y)
         if(ist.gt.3) then
            call pgsch(2.25)
         else
            call pgsch(1.25)
         end if
         call pgpt(nrl,x,y,pggmidx(ist))
         call pgsch(1.0)
C         call pgerry(nrl,x,yerr1,yerr2,0.0)
         call pgmtxt('T',-2.-(ist-1)*1.5,.85,0,
     $        char(pggmidx(ist))//'  '//
     $        typename(it,ist)(:lnb(typename(it,ist))))
      end do
      xline(1) = rlapmin
      xline(2) = rlapmin
      yline(1) = ymin
      yline(2) = ymax
      call pgsci(1)
      call pgsls(4)
      call pgline(2,xline,yline)
      call pgunsa

* Close PS device
      if(ips.gt.0) then
         if(jfrac(2,it).eq.0) then
            call pgslw(1)
            call pgclos
            write(6,'(2a)') ' Created PS file: ',psfile(:lnb(psfile))
         end if
         ips = 0
         ips2 = 1
         jfrac(2,it) = 1
      end if

      return
      end



************************************************************************
* subroutine plotfracall -- Plot all type fractions
************************************************************************

      subroutine plotfracall(nrl,rlarr,frl,frlerr)

      implicit none

      include 'snid.inc'
      include '../button/button.inc'      

      integer i
      integer nrl                        ! number of elements in rlap array (rlarr)
      real rlarr(MAXRLAP)                ! rlap array (rlaparr[0:rlapmax])
      real frl(MAXRLAP,NT,NST)           ! absolute fraction of templates with rlap >= rlap_i for each (sub)type
      real frlerr(MAXRLAP,NT,NST)        ! poisson error in absolute fraction of templates with rlap >= rlap_i for each (sub)type
      real x(MAXRLAP),y(MAXRLAP)         ! vectors for plotting
      real yerr1(MAXRLAP),yerr2(MAXRLAP) ! vectors for plotting
      integer nsub                       ! number of subtypes corresponding to type it
      character*200 psfile,asciifile     ! name of file for PS/ASCII output
      real xmin,xmax,ymin,ymax           ! x- and y-ranges for plot
      real xline(2),yline(2)             ! for plotting straight lines
      character*100 tit,xtit,ytit        ! (x-/y-)title for plot box
      integer pggmidx(NST)               ! PGPLOT graph marker indices

c Functions
      integer lnb,pgopen

* Erase previous plot      
      call rpgerasw(0.,1.,0.,Y3VPORT) 

* X- and Y-ranges
      xmin = 0.
      if(nrl.gt.MAXRLAP) nrl = MAXRLAP
      xmax = rlarr(nrl)
      if(xmax.eq.0) xmax=5
      ymin = 0.
      ymax = 1.

* Number of subtypes for type it
      nsub = NT

* ASCII output
      if(iascii.gt.0) then
         asciifile = froot(:lnb(froot))//'_snidfrac_all.dat'
         if(jfracall(1).eq.0) then
            open(unit=1,file=asciifile,status='replace')
            write(1,'(a)') '# Input: '//fdata(:lnb(fdata))
            write(1,'(a)') '# Type: All'
            write(1,'(a)') '#rlap fraction fraction_error type'
            do it=1, nsub
               do i=1, nrl
                  write(1,'(f5.1,2x,2(f7.4,2x),a10)') rlarr(i),
     $                 frl(i,it,1),frlerr(i,it,1),typename(it,1)
               end do
            end do
            close(1)
            write(6,'(2a)') ' Created ASCII file: ',
     $           asciifile(:lnb(asciifile))
            jfracall(1) = 1            
         else
            write(6,'(2a)') ' Already created ASCII file: ',
     $           asciifile(:lnb(asciifile))
         end if
      end if
      iascii = 0
      call button(48,'ASCII',iascii)

* Plot box; open PS device
      xtit='rlap [rlapmin = '//labels(1)(:lnb(labels(1)))//']'
      ytit='Template fraction \\(2244) rlap'
      tit='Template fraction for all types'
      if(ips.gt.0) then
         psfile = froot(:lnb(froot))//'_snidfrac_all.ps'
         if(jfracall(2).eq.0) then
            istat = pgopen(psfile//'/cps')
            if (istat.le.0) write(6,'(a)') 
     $           ' PGOPEN: Could not open PS device'
            call pgslw(3)
         else
            write(6,'(2a)') ' Already created PS file: ',
     $           psfile(:lnb(psfile))            
         end if
      end if
      call rpgenv(xmin,xmax,ymin,ymax,0,0)
      call pglabel(xtit,ytit,tit)

* Plot fraction curves and labels
      call pgsave
      call getpggm(pggmidx)
      do it=1, nsub
         do i=1, nrl
            x(i) = rlarr(i)
            y(i) = frl(i,it,1)
            yerr1(i) = frl(i,it,1)+frlerr(i,it,1)
            yerr2(i) = frl(i,it,1)-frlerr(i,it,1)
         end do
         call pgsci(it)
         call pgslw(3)
         call pgline(nrl,x,y)
         if(it.gt.3) then
            call pgsch(2.25)
         else
            call pgsch(1.25)
         end if
         call pgpt(nrl,x,y,pggmidx(it))
         call pgsch(1.0)
C         call pgerry(nrl,x,yerr1,yerr2,0.0)
         call pgmtxt('T',-2.-(it-1)*1.5,.85,0,
     $        char(pggmidx(it))//'  '//
     $        typename(it,1)(:lnb(typename(it,1))))
      end do
      xline(1) = rlapmin
      xline(2) = rlapmin
      yline(1) = ymin
      yline(2) = ymax
      call pgsci(1)
      call pgsls(4)
      call pgline(2,xline,yline)
      call pgunsa

* Close PS device
      if(ips.gt.0) then
         if(jfracall(2).eq.0) then
            call pgslw(1)
            call pgclos
            write(6,'(2a)') ' Created PS file: ',psfile(:lnb(psfile))
         end if
         ips = 0
         ips2 = 1
         jfracall(2) = 1
      end if

      return
      end



************************************************************************
* subroutine getpggm -- Get PGPLOT graph markers
************************************************************************

      subroutine getpggm(pggmidx)

      implicit none

      include 'snid.inc'

      integer pggmidx(NST) ! PGPLOT graph marker indices
      
      pggmidx(1) = 0  ! empty square
      pggmidx(2) = 4  ! empty circle
      pggmidx(3) = 7  ! empty triangle
      pggmidx(4) = 16 ! filled square
      pggmidx(5) = 17 ! filled circle
      pggmidx(6) = 13 ! filled triangle
      
      return
      end
