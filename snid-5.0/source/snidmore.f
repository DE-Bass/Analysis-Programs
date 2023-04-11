************************************************************************
**                                                                    **
**                      Extra routines for SNID                       **
**                                                                    **
**                                                                    **
**    Copyright (C) 1999-2007  Stephane Blondin and John L. Tonry     ** 
**                                                                    **
************************************************************************

c Included in this file:
c
c subroutine gettempdir
c subroutine initialize
c subroutine getdata
c subroutine gettemp
c subroutine docorr
c subroutine allpeaks
c subroutine getzt
c subroutine wfout
c subroutine wfluxout
c subroutine wflatout
c subroutine wxcorout
c subroutine loadspec



************************************************************************
* subroutine gettempdir -- Get the template directory
************************************************************************

      subroutine gettempdir
      
      implicit none

      include 'snid.inc'

c IMPORTANT! Fortran 77 lines are limited to 72 characters. 
c If the line below extends beyond the line of asterisks, you will
c need to split the directory name into substrings, as follows:
c
c      tempdir='/somereallylongdirectorypath_part1/'//
c     $     'somereallylongdirectorypath_part2/'//
c     $     'snid-5.0/templates/'
c
c [there needs to be 6 spaces before "tempdir" and 5 spaces before
c each "$". The "//" symbols tell the compiler that the line
c continues on the next line.]
c
c ===>    NOTHING MUST EXTEND BEYOND THE LINE OF ASTERISKS BELOW!   <===
c
************************************************************************
      tempdir='./snid-5.0/templates/'
************************************************************************

      return
      end



************************************************************************
* subroutine initialize -- Initialize the variables
************************************************************************

      subroutine initialize
      
      implicit none

      include 'snid.inc'

      integer i,ii,j,l
      character*200 inbuf       ! string buffer to parse command-line input
      character*100 ftemp       ! template name
      character*10 dumuse       ! for {use,avoid}type index determination
      integer ierr              ! >1 if {use,avoid}type error
      integer ip1,ip2           ! indices used to build use/avoid lists
      integer itst              ! test level [test=level option]
      real z0,dz                ! z0 +/- dz used to compute [zmin,zmax] range
      real age0,dage            ! age0 +/- dage used to compute [agemin,agemax] range
      real delta0,ddelta        ! delta0 +/- ddelta used to compute [deltamin,deltamax] range
      integer nepoch            ! number of epochs for a given template
      character*200 line        ! length of string to parse w,f lines
      integer ntok              ! number of words in line 

c Functions
      integer iargc,lnb,nlb

c
c Display syntax on call
c
      if(iargc().lt.1) then
         write(6,*) 'S u p e r N o v a    I D e n t i f i c a t i o n',
     $        ' (SNID v'//version(:lnb(version))//')'
         write(6,*) 'Copyright (C) 1999-2007 S. Blondin and J. L. Tonry'
         write(6,*) ' '
         write(6,*) 'Usage: snid [options] spec1.dat[,spec2.dat,...]',
     $	      ' [templates.lnw]'
         write(6,*) ' '
         write(6,*) 'options (default):'
         write(6,*) '  INPUT SPECTRUM:'
         write(6,*) '    wmin=               lower wavelength limit',
     $        ' (data)'
         write(6,*) '    wmax=               upper wavelength limit',
     $        ' (data)'
         write(6,*) '    wmask=              file containing',
     $        ' wavelength mask (NA)'
         write(6,*) '    fwmed=              median filtering width',
     $        ' [Angstroms] (0.)'
         write(6,*) '    medlen=             median filtering width',
     $        ' [pixels] (1.)'
         write(6,*) '    emclip=             clip emission at this',
     $        ' redshift (NA)'
         write(6,*) '    skyclip=0/1         clip sky emission',
     $        ' lines (0)'
         write(6,*) '    emwid=              emission clip width in',
     $        ' Angstroms (40)'
         write(6,*) '    aband=0/1           keep A band (1)'
         write(6,*) '  TEMPLATE SPECTRA:'
         write(6,*) '    agemin=             lower age limit (-90)'
         write(6,*) '    agemax=             upper age limit (+1000)'
         write(6,*) '    [age= dage=         alternate age range]'
         write(6,*) '    deltamin=           lower delta limit (-9.99)'
         write(6,*) '    deltamax=           upper delta limit (99.9)'
         write(6,*) '    [delta= ddelta=     alternate delta range]'
         write(6,*) '    nminspec=           minimum number of spectra',
     $        ' per template (1)'
         write(6,*) '    use=temp1,...       only use these templates',
     $        ' (NA)'
         write(6,*) '    usetype=type1,...   only use these (sub)types',
     $        ' (NA)'
         write(6,*) '    avoid=temp1,...     avoid these templates',
     $        ' (NA)'
         write(6,*) '    avoidtype=type1,... avoid these (sub)types',
     $        ' (NA)'
         write(6,*) '  REDSHIFT/CORRELATION:'
         write(6,*) '    rlapmin=            minimum rlap listed',
     $        ' (5.0)'
         write(6,*) '    lapmin=             minimum overlap',
     $        ' allowed (0.4)'
         write(6,*) '    zmin=               lower redshift limit',
     $        ' (-0.01)'
         write(6,*) '    zmax=               upper redshift limit',
     $        ' (1.2)'
         write(6,*) '    [z= dz=             alternate z range]'
         write(6,*) '    forcez=             forced initial redshift',
     $        ' (NA)'
C         write(6,*) '    zfake=              extra redshift to apply',
C     $        ' [TEST] (0.)'
C         write(6,*) '    znsig=              pick initial redshift',
C     $        ' by sigma clipping [TEST] (0.)'
         write(6,*) '    zfilter=            redshift filter for final',
     $        ' median (0.02)'
         write(6,*) '    k{1,2,3,4}=         bandpass filter range',
     $        ' (1,4,85,102)'
         write(6,*) '  OTHER:'
         write(6,*) '    param=              get options from',
     $        ' specified parameter file (NA)'
         write(6,*) '    fout=(-1)-n         file output level (0)'
         write(6,*) '    fluxout=n           write top n fluxed',
     $        ' spectra to file (0)'
         write(6,*) '    flatout=n           write top n flattened',
     $        ' spectra to file (0)'
         write(6,*) '    xcorout=n           write top n correlation',
     $        ' functions to file (0)'
         write(6,*) '    tempdir=            template directory',
     $        ' (see snidmore.f)'
         write(6,*) '    plot=0/1            enable plotting (1)'
         write(6,*) '    iquery=0/1          allow user to',
     $        ' change initial redshift (1)'
         write(6,*) '    inter=0/1           interactive (1)'
         write(6,*) '    verbose=0-2         verbose level (2)'
         write(6,*) ' '
         write(6,*) '[see file Howto.snid for more info]'
         stop
      end if 

c
c Can't override these parameters with options
c

* constant for formal error calculation, zerr = zk * [wid/(1+rlap)]
      zk = 2.0  ! this may change...

* How much to apodize?
c do NOT play around with unless you know what you're doing!
      percent = 5.0



c
c Default values for parameters (override with options)
c

* Change to 1 if options read from parameter file
      iparam = 0

* Wavelength range
      wmin  = 99999.99 ! 2000
      wmax  = -9.99    ! 10000
      iwmin = 0
      iwmax = 0

* Wavelength mask
      iwmask = 0
      wmaskfile = 'wmaskfile'

* FWHM (angstroms) of weighted median filter of spectrum
      fwmed = 0.0

* Range of median filtering of raw spectrum (pixels)
      medlen = 1.

* Redshift to clip emission lines
      emclipz = -9.99

* Clip sky lines?
      skyclip = 0

* Width to clip emission lines (A)
      emx = 40.

* Do we include the A band?
      iaband = 1

* Restrictions on template age
      agemin = -90.
      agemax = +1000.
      age0 = -99.9
      dage = 0.0

* Restrictions on template decline rate
      deltamin = -9.99
      deltamax = 99.9
      delta0 = -9.99
      ddelta = 0.00

* Minimum number of epochs per template
      nminspec = 1

* Template name to use regardless
      nuse = 0

* Template type to use regardless
      nusetype = 0
      do it=1, NT
         do ist=1, NST
            usetype(it,ist) = 0
         end do
      end do

* Template name to avoid
      navoid = 0

* Template type to avoid
      navoidtype = 0
      do it=1, NT
         do ist=1, NST
            avoidtype(it,ist) = 0
         end do
      end do

* Minimum rlap for listing
      rlapmin = 5.0

* Minimum overlap of spectra
      lapmin = 0.40

* Restrictions on redshift search range
      zmin = -0.01            
      zmax = +1.2
      z0 = -9.99
      dz = 0.00
      forcez = -9.99

* Extra z to apply for testing purposes
      zfake = 0.

* Determine initial redshift based on good correlations whose rlap
* value falls at least znsig above the rlap distribution
      znsig = 0.

* Maximum difference between median and template redshift
      zfilter = 0.02

* Wavenumbers for filtering 
      k1 = 0
      k2 = 0
      k3 = 0
      k4 = 0

* File output
      fout     = 0
      fluxout  = 0
      flatout  = 0
      xcorout  = 0

* Enable plotting?
      iplot = 1               

* Query the user to change initial redshift
      iquery = 1              

* Interactive mode?
      inter = 1

* Verbose level
      verbose = 2


c
c Have we got some options?
c
      if(iargc().lt.2) then
         noption = 1
         goto 8
      end if

      do 9 i=1, iargc()
         noption = i
         call getarg(noption,inbuf)
         if(index(inbuf, '=') .eq. 0) goto 8

* Read options from parameter file?
         if(inbuf(:6) .eq. 'param=') then
            read(inbuf(7:), '(a200)', err=9) paramfile
            iparam = 1
            open(unit=1,file=paramfile(:lnb(paramfile)),status='old')

            do 12 j=1, MAXPARAM
               read(1,'(a)',end=11,err=12) line
               if(index(line(nlb(line):nlb(line)+1),'#').gt.0) goto 12
               call parser(line,MAXTOK,ntok,ltok,tok)

               if(tok(1)(:ltok(1)) .eq. 'iwmin') then
                  read(tok(2),'(i20)',err=12) iwmin

               else if(tok(1)(:ltok(1)) .eq. 'wmin') then
                  read(tok(2),'(g20.0)',err=12) wmin

               else if(tok(1)(:ltok(1)) .eq. 'iwmax') then
                  read(tok(2),'(i20)',err=12) iwmax

               else if(tok(1)(:ltok(1)) .eq. 'wmax') then
                  read(tok(2),'(g20.0)',err=12) wmax

               else if(tok(1)(:ltok(1)) .eq. 'wmask') then
                  read(tok(2),'(a200)',err=12) wmaskfile

               else if(tok(1)(:ltok(1)) .eq. 'fwmed') then
                  read(tok(2),'(g20.0)',err=12) fwmed

               else if(tok(1)(:ltok(1)) .eq. 'medlen') then
                  read(tok(2),'(g20.0)',err=12) medlen

               else if(tok(1)(:ltok(1)) .eq. 'emclip') then
                  read(tok(2),'(g20.0)',err=12) emclipz

               else if(tok(1)(:ltok(1)) .eq. 'skyclip') then
                  read(tok(2),'(i20)',err=12) skyclip

               else if(tok(1)(:ltok(1)) .eq. 'emwid') then
                  read(tok(2),'(g20.0)',err=12) emx

               else if(tok(1)(:ltok(1)) .eq. 'aband') then
                  read(tok(2),'(i20)',err=12) iaband

               else if(tok(1)(:ltok(1)) .eq. 'agemin') then
                  read(tok(2),'(g20.0)',err=12) agemin

               else if(tok(1)(:ltok(1)) .eq. 'agemax') then
                  read(tok(2),'(g20.0)',err=12) agemax
                  
               else if(tok(1)(:ltok(1)) .eq. 'age') then
                  read(tok(2),'(g20.0)',err=12) age0

               else if(tok(1)(:ltok(1)) .eq. 'dage') then
                  read(tok(2),'(g20.0)',err=12) dage

               else if(tok(1)(:ltok(1)) .eq. 'deltamin') then
                  read(tok(2),'(g20.0)',err=12) deltamin

               else if(tok(1)(:ltok(1)) .eq. 'deltamax') then
                  read(tok(2),'(g20.0)',err=12) deltamax
                  
               else if(tok(1)(:ltok(1)) .eq. 'delta') then
                  read(tok(2),'(g20.0)',err=12) delta0
                  
               else if(tok(1)(:ltok(1)) .eq. 'ddelta') then
                  read(tok(2),'(g20.0)',err=12) ddelta
                  
               else if(tok(1)(:ltok(1)) .eq. 'nminspec') then
                  read(tok(2),'(i20)',err=12) nminspec

               else if(tok(1)(:ltok(1)).eq.'use'.and.ltok(1).eq.3) then
                  read(tok(2),'(i20)',err=12) nuse
                  if(nuse.gt.0) then
                     read(tok(3), '(a200)', err=9) uses
                     ip1 = 1
                     do l=1, nuse
                        if(l.lt.nuse) then
                           ip2 = index(uses(ip1:), ',')
                           use(l) = uses(ip1:ip1+ip2-2)
                           ip1 = ip1 + ip2
                        else
                           use(l) = uses(ip1:)
                        end if
                     end do
                  end if

               else if(tok(1)(:ltok(1)) .eq. 'usetype') then
                  read(tok(2),'(i20)',err=12) nusetype
                  if(nusetype.gt.0) then
                     read(tok(3),'(a200)',err=12) usets
                     ip1 = 1
                     do l=1, nusetype
                        if(l.lt.nusetype) then
                           ip2 = index(usets(ip1:), ',')
                           dumuse = usets(ip1:ip1+ip2-2)
                           call typeidx(dumuse,usetype,ierr)
                           if(ierr.gt.0) then
                              write(6,*) 'ERROR! Forcing use of',
     $                             ' non-existing type: ',dumuse
                              stop
                           end if
                           ip1 = ip1 + ip2
                        else
                           dumuse = usets(ip1:)
                           call typeidx(dumuse,usetype,ierr)
                           if(ierr.gt.0) then
                              write(6,*) 'ERROR! Forcing use of',
     $                             ' non-existing type: ',dumuse
                              stop
                           end if
                        end if
                     end do
                  end if

               else if(tok(1)(:ltok(1)).eq.'avoid'.and.ltok(1).eq.5)then
                  read(tok(2),'(i20)',err=12) navoid
                  if(navoid.gt.0) then
                     read(tok(3), '(a200)', err=9) avoids
                     ip1 = 1
                     do l=1, navoid
                        if(l.lt.navoid) then
                           ip2 = index(avoids(ip1:), ',')
                           avoid(l) = avoids(ip1:ip1+ip2-2)
                           ip1 = ip1 + ip2
                        else
                           avoid(l) = avoids(ip1:)
                        end if
                     end do
                  end if

               else if(tok(1)(:ltok(1)) .eq. 'avoidtype') then
                  read(tok(2),'(i20)',err=12) navoidtype
                  if(navoidtype.gt.0) then
                     read(tok(3),'(a200)',err=12) avoidts
                     ip1 = 1
                     do l=1, navoidtype
                        if(l.lt.navoidtype) then
                           ip2 = index(avoidts(ip1:), ',')
                           dumuse = avoidts(ip1:ip1+ip2-2)
                           call typeidx(dumuse,avoidtype,ierr)
                           if(ierr.gt.0) then
                              write(6,*) 'ERROR! Attempt to avoid',
     $                             ' non-existing type: ',dumuse
                              stop
                           end if
                           ip1 = ip1 + ip2
                        else
                           dumuse = avoidts(ip1:)
                           call typeidx(dumuse,avoidtype,ierr)
                           if(ierr.gt.0) then
                              write(6,*) 'ERROR! Attempt to avoid',
     $                             ' non-existing type: ',dumuse
                              stop
                           end if
                        end if
                     end do
                  end if

               else if(tok(1)(:ltok(1)) .eq. 'rlapmin') then
                  read(tok(2),'(g20.0)',err=12) rlapmin

               else if(tok(1)(:ltok(1)) .eq. 'lapmin') then
                  read(tok(2),'(g20.0)',err=12) lapmin

               else if(tok(1)(:ltok(1)) .eq. 'zmin') then
                  read(tok(2),'(g20.0)',err=12) zmin

               else if(tok(1)(:ltok(1)) .eq. 'zmax') then
                  read(tok(2),'(g20.0)',err=12) zmax

               else if(tok(1)(:ltok(1)).eq.'z' .and. ltok(1).eq.1) then
                  read(tok(2),'(g20.0)',err=12) z0

               else if(tok(1)(:ltok(1)).eq.'dz') then
                  read(tok(2),'(g20.0)',err=12) dz

               else if(tok(1)(:ltok(1)).eq.'forcez') then
                  read(tok(2),'(g20.0)',err=12) forcez

               else if(tok(1)(:ltok(1)).eq.'zfake') then
                  read(tok(2),'(g20.0)',err=12) zfake

               else if(tok(1)(:ltok(1)).eq.'znsig') then
                  read(tok(2),'(g20.0)',err=12) znsig

               else if(tok(1)(:ltok(1)).eq.'zfilter') then
                  read(tok(2),'(g20.0)',err=12) zfilter

               else if(tok(1)(:ltok(1)).eq.'k1') then
                  read(tok(2),'(i20)',err=12) k1

               else if(tok(1)(:ltok(1)).eq.'k2') then
                  read(tok(2),'(i20)',err=12) k2

               else if(tok(1)(:ltok(1)).eq.'k3') then
                  read(tok(2),'(i20)',err=12) k3

               else if(tok(1)(:ltok(1)).eq.'k4') then
                  read(tok(2),'(i20)',err=12) k4

               else if(tok(1)(:ltok(1)).eq.'tempdir') then
                  read(tok(2),'(a200)',err=12) tempdir

               else if(tok(1)(:ltok(1)).eq.'fout') then
                  read(tok(2),'(i20)',err=12) fout

               else if(tok(1)(:ltok(1)).eq.'fluxout') then
                  read(tok(2),'(i20)',err=12) fluxout

               else if(tok(1)(:ltok(1)).eq.'flatout') then
                  read(tok(2),'(i20)',err=12) flatout

               else if(tok(1)(:ltok(1)).eq.'xcorout') then
                  read(tok(2),'(i20)',err=12) xcorout

               else if(tok(1)(:ltok(1)).eq.'plot') then
                  read(tok(2),'(i20)',err=12) iplot

               else if(tok(1)(:ltok(1)).eq.'iquery') then
                  read(tok(2),'(i20)',err=12) iquery

               else if(tok(1)(:ltok(1)).eq.'inter') then
                  read(tok(2),'(i20)',err=12) inter

               else if(tok(1)(:ltok(1)).eq.'verbose') then
                  read(tok(2),'(i20)',err=12) verbose

               else
                  write(6,*) 'ERROR! Illegal option: ',tok(1)(:ltok(1))
                  stop
               end if

 12         continue
 11         close(1)

* Other options
         else if(inbuf(:5) .eq. 'wmin=') then
            read(inbuf(6:), '(g20.0)', err=9) wmin
            iwmin = 1

         else if(inbuf(:5) .eq. 'wmax=') then
            read(inbuf(6:), '(g20.0)', err=9) wmax
            iwmax = 1

         else if(inbuf(:6) .eq. 'wmask=') then
            read(inbuf(7:), '(a200)', err=9) wmaskfile
            iwmask = 1

         else if(inbuf(:6) .eq. 'fwmed=') then
            read(inbuf(7:), '(g20.0)', err=9) fwmed

         else if(inbuf(:7) .eq. 'medlen=') then
            read(inbuf(8:), '(i20)', err=9) medlen

         else if(inbuf(:7) .eq. 'emclip=') then
            read(inbuf(8:), '(g20.0)', err=9) emclipz

         else if(inbuf(:8) .eq. 'skyclip=') then
            read(inbuf(9:), '(i20)', err=9) skyclip

         else if(inbuf(:6) .eq. 'emwid=') then
            read(inbuf(7:), '(g20.0)', err=9) emx

         else if(inbuf(:6) .eq. 'aband=') then
            read(inbuf(7:), '(i5)', err=9) iaband

         else if(inbuf(:7) .eq. 'agemin=') then
            read(inbuf(8:), '(g20.0)', err=9) agemin

         else if(inbuf(:7) .eq. 'agemax=') then
            read(inbuf(8:), '(g20.0)', err=9) agemax

         else if(inbuf(:4) .eq. 'age=') then
            read(inbuf(5:), '(g20.0)', err=9) age0

         else if(inbuf(:5) .eq. 'dage=') then
            read(inbuf(6:), '(g20.0)', err=9) dage

         else if(inbuf(:9) .eq. 'deltamin=') then
            read(inbuf(10:), '(g20.0)', err=9) deltamin

         else if(inbuf(:9) .eq. 'deltamax=') then
            read(inbuf(10:), '(g20.0)', err=9) deltamax

         else if(inbuf(:6) .eq. 'delta=') then
            read(inbuf(7:), '(g20.0)', err=9) delta0

         else if(inbuf(:7) .eq. 'ddelta=') then
            read(inbuf(8:), '(g20.0)', err=9) ddelta

         else if(inbuf(:9) .eq. 'nminspec=') then
            read(inbuf(10:), '(i20)', err=9) nminspec

         else if(inbuf(:4) .eq. 'use=') then
            read(inbuf(5:), '(a200)', err=9) uses
            ip1 = 1
            nuse = 0
            do ii=1, MAXUSE
               ip2 = index(uses(ip1:), ',')
               nuse = nuse + 1
               if(ip2.gt.0) then
                  use(nuse) = uses(ip1:ip1+ip2-2)
                  ip1 = ip1 + ip2
               else
                  use(nuse) = uses(ip1:)
                  goto 613
               end if
            end do
 613        continue

         else if(inbuf(:8) .eq. 'usetype=') then
            read(inbuf(9:), '(a200)', err=9) usets
            ip1 = 1
            nusetype = 0
            do ii=1, MAXUSE
               ip2 = index(usets(ip1:), ',')
               nusetype = nusetype + 1
               if(ip2.gt.0) then
                  dumuse = usets(ip1:ip1+ip2-2)
                  call typeidx(dumuse,usetype,ierr)
                  if(ierr.gt.0) then
                     write(6,*) 'ERROR! Forcing use of',
     $                    ' non-existing type: ',dumuse
                     stop
                  end if
                  ip1 = ip1 + ip2
               else
                  dumuse = usets(ip1:)
                  call typeidx(dumuse,usetype,ierr)
                  if(ierr.gt.0) then
                     write(6,*) 'ERROR! Forcing use of',
     $                    ' non-existing type: ',dumuse
                     stop
                  end if
                  goto 615
               end if
            end do
 615        continue

         else if(inbuf(:6) .eq. 'avoid=') then
            read(inbuf(7:), '(a200)', err=9) avoids
            ip1 = 1
            navoid = 0
            do ii=1, MAXUSE
               ip2 = index(avoids(ip1:), ',')
               navoid = navoid + 1
               if(ip2.gt.0) then
                  avoid(navoid) = avoids(ip1:ip1+ip2-2)
                  ip1 = ip1 + ip2
               else
                  avoid(navoid) = avoids(ip1:)
                  goto 609
               end if
            end do
 609        continue

         else if(inbuf(:10) .eq. 'avoidtype=') then
            read(inbuf(11:), '(a200)', err=9) avoidts
            ip1 = 1
            navoidtype = 0
            do ii=1, MAXUSE
               ip2 = index(avoidts(ip1:), ',')
               navoidtype = navoidtype + 1
               if(ip2.gt.0) then
                  dumuse = avoidts(ip1:ip1+ip2-2)
                  call typeidx(dumuse,avoidtype,ierr)
                  if(ierr.gt.0) then
                     write(6,*) 'ERROR! Attempt to avoid',
     $                    ' non-existing type: ',dumuse
                     stop
                  end if                        
                  ip1 = ip1 + ip2
               else
                  dumuse = avoidts(ip1:)
                  call typeidx(dumuse,avoidtype,ierr)
                  if(ierr.gt.0) then
                     write(6,*) 'ERROR! Attempt to avoid',
     $                    ' non-existing type: ',dumuse
                     stop
                  end if                        
                  goto 611
               end if
            end do
 611        continue

         else if(inbuf(:8) .eq. 'rlapmin=') then
            read(inbuf(9:), '(g20.0)', err=9) rlapmin

         else if(inbuf(:7) .eq. 'lapmin=') then
            read(inbuf(8:), '(g20.0)', err=9) lapmin

         else if(inbuf(:5) .eq. 'zmin=') then
            read(inbuf(6:), '(g20.0)', err=9) zmin

         else if(inbuf(:5) .eq. 'zmax=') then
            read(inbuf(6:), '(g20.0)', err=9) zmax

         else if(inbuf(:2)     .eq. 'z=') then
            read(inbuf(3:), '(g20.0)', err=9) z0

         else if(inbuf(:3) .eq. 'dz=') then
            read(inbuf(4:), '(g20.0)', err=9) dz

         else if(inbuf(:7) .eq. 'forcez=') then
            read(inbuf(8:), '(g20.0)', err=9) forcez

         else if(inbuf(:6) .eq. 'zfake=') then
            read(inbuf(7:), '(g20.0)', err=9) zfake

         else if(inbuf(:6) .eq. 'znsig=') then
            read(inbuf(7:), '(g20.0)', err=9) znsig

         else if(inbuf(:8) .eq. 'zfilter=') then
            read(inbuf(9:), '(g20.0)', err=9) zfilter

         else if(inbuf(:3) .eq. 'k1=') then
            read(inbuf(4:), '(i20)', err=9) k1

         else if(inbuf(:3) .eq. 'k2=') then
            read(inbuf(4:), '(i20)', err=9) k2

         else if(inbuf(:3) .eq. 'k3=') then
            read(inbuf(4:), '(i20)', err=9) k3

         else if(inbuf(:3) .eq. 'k4=') then
            read(inbuf(4:), '(i20)', err=9) k4

         else if(inbuf(:8) .eq. 'tempdir=') then
            read(inbuf(9:), '(a200)', err=9) tempdir

         else if(inbuf(:5) .eq.'fout=') then
            read(inbuf(6:), '(i20)', err=9) fout
            
         else if(inbuf(:8) .eq.'fluxout=') then
            read(inbuf(9:), '(i20)', err=9) fluxout
            
         else if(inbuf(:8) .eq.'flatout=') then
            read(inbuf(9:), '(i20)', err=9) flatout
            
         else if(inbuf(:8) .eq.'xcorout=') then
            read(inbuf(9:), '(i20)', err=9) xcorout
            
         else if(inbuf(:5) .eq. 'plot=') then
            read(inbuf(6:), '(i20)', err=9) iplot

         else if(inbuf(:7) .eq. 'iquery=') then
            read(inbuf(8:), '(i20)', err=9) iquery

         else if(inbuf(:6) .eq. 'inter=') then
            read(inbuf(7:), '(i20)', err=9) inter

         else if(inbuf(:8) .eq. 'verbose=') then
            read(inbuf(9:), '(i20)', err=9) verbose

         else
            write(6,*) 'Illegal option: ',inbuf(:index(inbuf,"=")-1)
            stop

         end if

 9    continue

 8    noption = noption - 1

c
c Check options
c

* wmin,wmax
      if((iwmin+iwmax.eq.2).and.wmax.lt.wmin) then
         write(6,1000) ' ERROR! wmin = ',wmin,' > wmax = ',wmax
         stop
      end if

* fwmed
      if(fwmed.lt.0) then
         write(6,*) 'ERROR! Negative median filtering width'
         write(6,'(a,f7.1)') ' fwmed [Angstroms] = ',fwmed
         stop
      end if

* medlen
      if(medlen.lt.0) then
         write(6,*) 'ERROR! Negative median filtering width'
         write(6,'(a,f7.1)') ' medlen [pixels] = ',medlen
         stop
      end if

* emclip
      if(emclipz.gt.-1.and.(emclipz.lt.zmin.or.emclipz.gt.zmax)) then
         write(6,*) 'ERROR! Clipping emission lines beyond',
     $        ' redshift bounds'
         if(emclipz.lt.zmin) then
            write(6,1001) ' emclipz = ',emclipz,' < zmin = ',zmin
         else
            write(6,1001) ' emclipz = ',emclipz,' > zmax = ',zmax
         end if
         stop
      end if

      if(emclipz.gt.-1.and.z0.gt.-1) then
         if(dz.gt.0) then
            if(emclipz.lt.z0-dz.or.emclipz.gt.z0+dz) then
               write(6,*) 'ERROR! Clipping emission lines beyond',
     $              ' redshift bounds'
               if(emclipz.lt.z0-dz) then
                  write(6,1001) ' emclipz = ',emclipz,' < zmin = ',z0-dz
               else
                  write(6,1001) ' emclipz = ',emclipz,' > zmax = ',z0+dz
               end if
               stop
            end if
         else
            if(emclipz.lt.z0-.1.or.emclipz.gt.z0+.1) then
               write(6,*) 'ERROR! Clipping emission lines beyond',
     $              ' redshift bounds'
               if(emclipz.lt.z0-.1) then
                  write(6,1001) ' emclipz = ',emclipz,' < zmin = ',z0-.1
               else
                  write(6,1001) ' emclipz = ',emclipz,' > zmax = ',z0+.1
               end if
               stop
            end if            
         end if
      end if

* emwid
      if(emx.lt.0) then
         write(6,*) 'ERROR! Negative emission clip width'
         write(6,'(a,f7.1)') ' emwid [Angstroms] = ',emx
         stop
      end if

* agemin,agemax
      if(agemax.lt.agemin) then
         write(6,1000) ' ERROR! agemin = ',agemin,' > agemax = ',agemax
         stop
      end if     

      if(age0.gt.-99.and.dage.lt.0) then
         write(6,'(a,f7.1)') ' ERROR! Negative dage: ',dage         
         stop
      end if

* deltamin,deltamax
      if(deltamax.lt.deltamin) then
         write(6,1000) ' ERROR! deltamin = ',deltamin,' > deltamax = ',
     $        deltamax
         stop
      end if     

      if(delta0.gt.-99.and.ddelta.lt.0) then
         write(6,'(a,f7.1)') ' ERROR! Negative ddelta: ',ddelta         
         stop
      end if

* nminspec
      if(nminspec.gt.MAXEPOCH) then
         write(6,*) 'ERROR! Cannot request more than MAXEPOCH spectra',
     $        ' per template'
         write(6,1002) ' nminspec = ',nminspec,' > MAXEPOCH = ',
     $        MAXEPOCH
         stop
      end if

* rlapmin
      if(rlapmin.lt.0) then
         write(6,'(a,f7.1)') ' ERROR! Negative rlapmin = ',rlapmin
         stop
      end if

* lapmin
      if(lapmin.lt.0) then
         write(6,'(a,f7.2)') ' ERROR! Negative lapmin = ',lapmin
         stop
      end if

* zmin,zmax
      if(zmax.lt.zmin) then
         write(6,1001) ' ERROR! zmin = ',zmin,' > zmax = ',zmax
         stop
      end if     

      if(z0.gt.-1.and.dz.lt.0) then
         write(6,'(a,f6.3)') ' ERROR! Negative dz: ',dz         
         stop
      end if

* forcez
      if(forcez.gt.-1.and.(forcez.lt.zmin.or.forcez.gt.zmax)) then
         write(6,*) 'ERROR! Forcing initial redshift beyond',
     $        ' redshift bounds'
         if(forcez.lt.zmin) then
            write(6,1001) ' forcez = ',forcez,' < zmin = ',zmin
         else
            write(6,1001) ' forcez = ',forcez,' > zmax = ',zmax
         end if
         stop
      end if

      if(forcez.gt.-1.and.z0.gt.-1) then
         if(dz.gt.0) then
            if(forcez.lt.z0-dz.or.forcez.gt.z0+dz) then
               write(6,*) 'ERROR! Clipping emission lines beyond',
     $              ' redshift bounds'
               if(forcez.lt.z0-dz) then
                  write(6,1001) ' forcez = ',forcez,' < zmin = ',z0-dz
               else
                  write(6,1001) ' forcez = ',forcez,' > zmax = ',z0+dz
               end if
               stop
            end if
         else
            if(forcez.lt.z0-.1.or.forcez.gt.z0+.1) then
               write(6,*) 'ERROR! Clipping emission lines beyond',
     $              ' redshift bounds'
               if(forcez.lt.z0-.1) then
                  write(6,1001) ' forcez = ',forcez,' < zmin = ',z0-.1
               else
                  write(6,1001) ' forcez = ',forcez,' > zmax = ',z0+.1
               end if
               stop
            end if            
         end if
      end if

* znsig
      if(znsig.lt.0) then
         write(6,'(a,f6.3)') ' ERROR! Negative sigma: znsig = ',znsig
         stop
      end if

* zfilter
      if(zfilter.le.0) then
         write(6,'(2a,f6.3)') ' ERROR! Negative or null redshift',
     $        ' filter: zfilter = ',zfilter
         stop
      end if

* k1,k2,k3,k4
c checked at a later stage (search for "check k1,k2,k3,k4")

 1000 format(2(a,f7.1))
 1001 format(2(a,f6.3))
 1002 format(2(a,i5))

c
c Write options to parameter file
c
      open(unit=1,file='snid.param',status='replace')
      write(1,'(a)') '# SNID parameter file -- EDIT WITH CAUTION!'
      write(1,'(a,f10.2)') 'wmin       ',wmin
      write(1,'(a,f10.2)') 'wmax       ',wmax
      write(1,'(2a)')      'wmask      ',wmaskfile(:lnb(wmaskfile))
      write(1,'(a,f10.2)') 'fwmed      ',fwmed
      write(1,'(a,f10.2)') 'medlen     ',medlen
      write(1,'(a,f10.4)') 'emclip     ',emclipz
      write(1,'(a,i10)')   'skyclip    ',skyclip
      write(1,'(a,f10.2)') 'emwid      ',emx
      write(1,'(a,i10)')   'aband      ',iaband
      write(1,'(a,f10.2)') 'agemin     ',agemin
      write(1,'(a,f10.2)') 'agemax     ',agemax
      write(1,'(a,f10.2)') 'age        ',age0
      write(1,'(a,f10.2)') 'dage       ',dage
      write(1,'(a,f10.2)') 'deltamin   ',deltamin
      write(1,'(a,f10.2)') 'deltamax   ',deltamax
      write(1,'(a,f10.2)') 'delta      ',delta0
      write(1,'(a,f10.2)') 'ddelta     ',ddelta
      write(1,'(a,i10)')   'nminspec   ',nminspec
      if(nuse.gt.0) then
         write(1,'(a,i10,2x,a)')
     $                     'use        ',nuse,uses(:lnb(uses))
      else
         write(1,'(a,i10)')'use        ',nuse
      endif
      if(nusetype.gt.0) then
         write(1,'(a,i10,2x,a)')
     $                     'usetype    ',nusetype,usets(:lnb(usets))
      else
         write(1,'(a,i10)')'usetype    ',nusetype
      endif
      if(navoid.gt.0) then
         write(1,'(a,i10,2x,a)')
     $                     'avoid      ',navoid,avoids(:lnb(avoids))
      else
         write(1,'(a,i10)')'avoid      ',navoid
      endif
      if(navoidtype.gt.0) then
         write(1,'(a,i10,2x,a)') 
     $                     'avoidtype  ',navoidtype,
     $        avoidts(:lnb(avoidts))
      else
         write(1,'(a,i10)')'avoidtype  ',navoidtype
      endif
      write(1,'(a,f10.2)') 'rlapmin    ',rlapmin
      write(1,'(a,f10.3)') 'lapmin     ',lapmin
      write(1,'(a,f10.4)') 'zmin       ',zmin
      write(1,'(a,f10.4)') 'zmax       ',zmax
      write(1,'(a,f10.4)') 'z          ',z0
      write(1,'(a,f10.4)') 'dz         ',dz
      write(1,'(a,f10.4)') 'forcez     ',forcez
      write(1,'(a,f10.4)') 'zfake      ',zfake
      write(1,'(a,f10.2)') 'znsig      ',znsig
      write(1,'(a,f10.4)') 'zfilter    ',zfilter
      write(1,'(a,i10)')   'k1         ',k1
      write(1,'(a,i10)')   'k2         ',k2
      write(1,'(a,i10)')   'k3         ',k3
      write(1,'(a,i10)')   'k4         ',k4
      write(1,'(2a)')      'tempdir    ',tempdir(:lnb(tempdir))
      write(1,'(a,i10)')   'fout       ',fout 
      write(1,'(a,i10)')   'fluxout    ',fluxout 
      write(1,'(a,i10)')   'flatout    ',flatout
      write(1,'(a,i10)')   'xcorout    ',xcorout 
      write(1,'(a,i10)')   'plot       ',iplot
      write(1,'(a,i10)')   'iquery     ',iquery
      write(1,'(a,i10)')   'inter      ',inter
      write(1,'(a,i10)')   'verbose    ',verbose    
      close(1)

* Can't be interactive if no terminal output!
      if(verbose.eq.0) inter=0

* Merge the constraints on redshift
      if(z0.gt.-1) then
         if(dz.gt.0) then
            zmin = z0 - dz
            zmax = z0 + dz
         else
            zmin = z0 - 0.1
            zmax = z0 + 0.1
         end if
      end if
C      z0 = 0.5*(zmin+zmax)
C      dz = 0.5*(zmax-zmin)

* Merge the constraints on age
      if(age0.gt.-99) then
         if(dage.gt.0) then
            agemin = age0 - dage
            agemax = age0 + dage
         else
            agemin = age0 - 3.
            agemax = age0 + 3.
         end if
      end if

* Merge the constraints on delta
      if(delta0.gt.-9) then
         if(ddelta.gt.-9) then
            deltamin = delta0 - ddelta
            deltamax = delta0 + ddelta
         else
            deltamin = delta0 - 0.2
            deltamax = delta0 + 0.2
         end if
      end if

      medlen = min(medlen,real(MAXLOG)-1)

* Sample the first template file to find out parameters
      if(iargc().lt.noption+2) then
         open(unit=1,file=tempdir(:lnb(tempdir))//'templist',
     $        status='old')
         read(1,'(a)') ftemp
         close(1)
         open(unit=1,file=tempdir(:lnb(tempdir))//ftemp,status='old')
      else
         call getarg(noption+2,inbuf)
         open(unit=1,file=inbuf,status='old')
      endif
      read(1,1101) nepoch, nw, w0, w1
 1101 format(2i5,2f10.2)
      if(nw.gt.MAXLOG) then
         write(6,1102) nw,MAXLOG
 1102    format(' (nw,MAXLOG) = (',i5,',',i5,')')
         write(6,*) 'Allocate more space for MAXLOG (in snid.inc)!'
         stop
      end if

* Set the log wavelength variables 
      dwlog = alog(w1/w0) / nw
      log2n = nint(alog(float(nw))/alog(2.))
      lz1 = nint(alog(zmin+1) / dwlog)   ! 0
      lz2 = nint(alog(zmax+1) / dwlog)   ! nw - 1
      close (1)

* Wavenumbers for filtering 
      if(k1.le.0.or.k1.gt.nw) k1 = 1
      if(k2.le.0.or.k2.gt.nw) k2 = 4
      if(k3.le.0.or.k3.gt.nw) k3 = nw/12
      if(k4.le.0.or.k4.gt.nw) k4 = nw/10

c check k1,k2,k3,k4
      if(.not.(k1.lt.k2.and.k2.lt.k3.and.k3.lt.k4)) then
         write(6,'(2a)') ' ERROR! Wavenumbers for filtering must be',
     $        ' in strict increasing order'
          write(6,'(a,4i4)') ' (k1,k2,k3,k4) = ',k1,k2,k3,k4
         stop
      end if


c
c Verbose output
c
      if(verbose.gt.0) then
         write(6,*) 'S u p e r N o v a    I D e n t i f i c a t i o n',
     $        ' (SNID v'//version(:lnb(version))//')'
         write(6,*) ' '
         if(iparam.gt.0) then
            write(6,'(2a)') ' Reading options from parameter file: ',
     $           paramfile(:lnb(paramfile))
         end if
         write(6,6112) zmin, zmax, rlapmin
 6112    format(' Searching in redshift range:',2f7.3,
     $        ' ; rlapmin = ',f5.2)
         write(6,6113) agemin, agemax
 6113    format(' Restricting to age range:',2f7.1)
         write(6,6114) deltamin, deltamax
 6114    format(' Restricting to delta range:',2f7.1)
         if(medlen.gt.1) then
            write(6,*) 'Median filtering with width ', medlen
         end if
         if(zfake.gt.0) then
            write(6,*) 'Fake redshift factor of ', zfake
         end if
         if(nuse.gt.0) then
            write(6,'(2a)') ' Forcing use of templates: ',
     $           uses(:lnb(uses))
         end if
         if(nusetype.gt.0) then
            write(6,'(2a)') ' Forcing use of spectral types: ',
     $           usets(:lnb(usets))
         end if
         if(navoid.gt.0) then
            write(6,'(2a)') ' Avoiding templates: ',
     $           avoids(:lnb(avoids))
         end if
         if(navoidtype.gt.0) then
            write(6,'(a,30a)') ' Avoiding spectral types: ',
     $           avoidts(:lnb(avoidts))
         end if
      end if

      return
      end



************************************************************************
* subroutine getdata -- Read in the input spectrum
************************************************************************

      subroutine getdata(fnorm,nknod,xknod,yknod)

      implicit none

      include 'snid.inc'

      integer i,l,l0                       ! loop indeces
      integer ip1,ip2                      ! indices used to build use/avoid lists
      character*200 fpdata                 ! input spectra, including path
      character*100 fthis                  ! one such input spectrum
      character*100 tmpfroot               ! temporary root name of input spectral file
      character*100 tmpfdata               ! temporary name of input spectral file
      character*200 line                   ! length of string to parse w,f lines
      integer ntok                         ! number of words in line 
      real flog(MAXLOG)                    ! input flux vector (binned in log lambda)
      integer nwave                        ! index for wave/flux vectors
      real fnorm(MAXLOG)                   ! normalized input flux vector
      real medbuf(MAXWAVE)                 ! buffer for median filtering
      real buf(MAXWAVE)                    ! buffer for median filtering
      real wminmax,fminmax                 ! wavelength/flux range of input spectrum
      real wmask0(MAXWAVE),wmask1(MAXWAVE) ! wavelength bounds of wavelength mask
      integer nwmask                       ! number of wavelength intervals in wavelength mask
      real faccum(MAXLOG)                  ! total input flux vector
      real aband0,aband1                   ! wavelength limits for atmospheric A-band
      real w                               ! wavelength of redshifted emission lines
      real z0                              ! z0 +/- dz used to compute [zmin,zmax] range
      integer izoff                        ! redshift offset (in log bins) used in spline fit
      integer l1,l2                        ! bounds (in log bins) for spline fit
      real fluxmean                        ! mean (spline-subtracted) flux vector
      integer nmean                        ! number of flux points in fluxmean

c Functions
      integer lnb,nlb

* Get filename and root
      call getarg(noption+1,fpdata)
      call getroot(fpdata,fdata,tmpfroot)

* Loop over different contributing data spectra
      call zerospec(nw,flog)
      ip1 = 1
 30   continue
      ip2 = index(fpdata(ip1:),',')
      if(ip2.eq.0) then
         fthis = fpdata(ip1:)
         call getroot(fthis,tmpfdata,tmpfroot)
         if(ip1.eq.1) then
            froot = tmpfroot            
         else
            froot = froot(:lnb(froot))//'_COMB_'//tmpfroot            
         end if
      else
         fthis = fpdata(ip1:ip1+ip2-2)
         call getroot(fthis,tmpfdata,tmpfroot)
         if(ip1.eq.1) then 
            froot = tmpfroot
         else
            froot = froot(:lnb(froot))//'_COMB_'//tmpfroot
         end if
      end if

* Verbose information
      if(verbose.gt.0) then
         if(ip1.eq.1) then
            write(6,'(2a)') ' Reading data file: ',
     $           tmpfdata(:lnb(tmpfdata))
         else
            write(6,'(2a)') ' Adding data file: ',
     $           tmpfdata(:lnb(tmpfdata))
         end if
      end if
      ip1 = ip1 + ip2

* Restrict wmin,wmax to range of input spectrum (if not set by the user)
      if(iwmin.eq.0.or.iwmax.eq.0) then
         open(unit=1,file=fthis,status='old')
         do l=1, MAXWAVE
            read(1,'(a)',end=33) line
            if(index(line(nlb(line):nlb(line)+1),'#').eq.0) then
               call parser(line,MAXTOK,ntok,ltok,tok)
               read(tok(1),*) wminmax
               read(tok(2),*) fminmax
               if(iwmin.eq.0) wmin=min(wmin,wminmax)
               if(iwmax.eq.0) wmax=max(wmax,wminmax)
            end if
         end do
 33      close(1)
      end if      
      if(verbose.gt.0.and.ip2.eq.0) write(6,'(a,2(f9.1))')
     $     ' Restricting to wavelength range: ',wmin,wmax


* Apply wavelength mask?
      nwmask = 0
      if(iwmask.gt.0) then
         if(verbose.gt.0.and.ip2.eq.0) write(6,'(2a)')
     $        ' Applying wavelength mask from file: ',
     $        wmaskfile(:lnb(wmaskfile))
         open(unit=2,file=wmaskfile(:lnb(wmaskfile)),status='old')
         do 20 i=1, MAXWAVE
            read(2,'(a)',end=21) line
            if(index(line(nlb(line):nlb(line)+1),'#').gt.0) goto 20
            call parser(line,MAXTOK,ntok,ltok,tok)
            nwmask = nwmask + 1
            read(tok(1),*) wmask0(nwmask)
            read(tok(2),*) wmask1(nwmask)
            if(wmask1(nwmask).lt.wmask0(nwmask)) then
               write(6,*) ' '
               write(6,*) 'WARNING! Bounds for wavelength mask must',
     $              'be in increasing order! '
               write(6,*) '(w0,w1) = ',wmask0(nwmask),wmask1(nwmask)
               write(6,*) 'Make sure to edit '//
     $              wmaskfile(:lnb(wmaskfile))
               write(6,*) ' '
               stop
            end if
 20      continue
 21      close(2)
      endif

      nwave = 0
      aband0 = 7575
      aband1 = 7675
      open(unit=1,file=fthis,status='old')

c
c Read in wavelength/flux values
c

* Verbose output
      if(skyclip.gt.0) then
         if(verbose.gt.0.and.ip2.eq.0)
     $        write(6,*) 'Clipping sky lines'
      end if
      
      if(emclipz.gt.zmin) then
         if(verbose.gt.0.and.ip2.eq.0)
     $        write(6,'(a,f5.3)') ' Clipping emission at z= ',emclipz
      end if

      do 4 l=1, MAXWAVE
         read(1,'(a)',end=5) line
         if(index(line(nlb(line):nlb(line)+1),'#').gt.0) goto 4
         call parser(line,MAXTOK,ntok,ltok,tok)
         read(tok(1),*) wave(nwave+1)
         read(tok(2),*) flux(nwave+1)

* Clip out the A band?
         if(iaband.eq.0.and.wave(nwave+1).ge.aband0.and.
     $        wave(nwave+1).le.aband1) goto 4

* Restrict wavelength range?
         if(wave(nwave+1).lt.wmin.or.wave(nwave+1).gt.wmax) goto 4

* Apply wavelength mask?
         if(iwmask.gt.0) then
            do 22 i=1, nwmask
               if(wave(nwave+1).lt.wmask0(i).or.
     $              wave(nwave+1).gt.wmask1(i)) then
                  goto 22
               else if(wave(nwave+1).ge.wmask0(i).and.
     $                 wave(nwave+1).le.wmask1(i)) then
                  goto 4
               end if
 22         continue
         endif

* Clip out emission lines (sky and galaxy)?
         if(skyclip.gt.0) then
            w = 5577.  ! OI sky line
            if(wave(nwave+1).ge.w-emx.and.wave(nwave+1).le.w+emx) goto 4
            w = 6300.2 ! OI sky line
            if(wave(nwave+1).ge.w-emx.and.wave(nwave+1).le.w+emx) goto 4
            w = 6364.  ! OI sky line
            if(wave(nwave+1).ge.w-emx.and.wave(nwave+1).le.w+emx) goto 4
         end if

         if(emclipz.gt.zmin) then
            w = 3727.3 * (1+emclipz)   ! OII
            if(wave(nwave+1).ge.w-emx.and.wave(nwave+1).le.w+emx) goto 4
            w = 4861.3 * (1+emclipz)   ! Hbeta
            if(wave(nwave+1).ge.w-emx.and.wave(nwave+1).le.w+emx) goto 4
            w = 4958.9 * (1+emclipz)   ! OIII
            if(wave(nwave+1).ge.w-emx.and.wave(nwave+1).le.w+emx) goto 4
            w = 5006.8 * (1+emclipz)   ! OIII
            if(wave(nwave+1).ge.w-emx.and.wave(nwave+1).le.w+emx) goto 4
            w = 6548.1 * (1+emclipz)   ! N1
            if(wave(nwave+1).ge.w-emx.and.wave(nwave+1).le.w+emx) goto 4
            w = 6562.8 * (1+emclipz)   ! Halpha
            if(wave(nwave+1).ge.w-emx.and.wave(nwave+1).le.w+emx) goto 4
            w = 6583.6 * (1+emclipz)   ! N2
            if(wave(nwave+1).ge.w-emx.and.wave(nwave+1).le.w+emx) goto 4
            w = 6716.4 * (1+emclipz)   ! S1
            if(wave(nwave+1).ge.w-emx.and.wave(nwave+1).le.w+emx) goto 4
            w = 6730.8 * (1+emclipz)   ! S2
            if(wave(nwave+1).ge.w-emx.and.wave(nwave+1).le.w+emx) goto 4
          end if

* Acceptable point
         nwave = nwave + 1

* Maybe to be buggered up by a test redshift!
         wave(nwave) = wave(nwave) * (1.+zfake)

 4    continue
 5    close(1)

* Rebin onto the common log wavelength scale
      if(fwmed.gt.0) then
         if(verbose.gt.0.and.ip2.eq.0)
     $        write(6,*) 'Median smoothing with fwmed',fwmed
         call medwfilt(nwave,wave,flux,fwmed,medbuf,buf)
      else if(medlen.gt.1) then
         call medfilt(nwave,flux,medlen,medbuf,buf)
      end if
      call rebin(nwave,wave,flux,nw,w0,dwlog,faccum)
      call addspec(nw,flog,faccum)

* Read in next spectrum if more than one input
      if(ip2.gt.0) goto 30

* Remove the mean with a spline fit
      izoff = -1
      z0 = 0.5*(zmin+zmax)
      if(z0.gt.-1) izoff = nint(alog(z0+1)/dwlog)
      call meanzero(nw,l1,l2,flog,fnorm,izoff,nknod,xknod,yknod)

* Remove the overall scale factor
      fluxmean = 0
      nmean = 0
      do i=1, nw
         if(flog(i).gt.0) then
            fluxmean = fluxmean + flog(i)
            nmean = nmean + 1
         end if
      end do
      fluxmean = alog10(fluxmean/nmean)
      do i=1, nknod
         yknod(i) = yknod(i) - fluxmean
      end do

* Apodize the end
      call apodize(nw,l1,l2,fnorm,percent)

      return
      end



************************************************************************
* subroutine gettemp -- Read in the template spectra
************************************************************************

      subroutine gettemp(ntemp,sname,stype,itype,delta,epoch,tflag,temp,
     $     nknot,xknot,yknot)

      implicit none
      
      include 'snid.inc'

      integer i,l,l0,ii,n,j,k
      integer ittype,itstype
      integer ktemp
      integer inouse
      integer nepoch
      integer nwx
      real w0x,w1x,ww
      real dta
      integer mostknots
      integer nk(MAXEPOCH)
      real fmean
      real tmp(MAXEPOCH),ep(MAXEPOCH)
      integer tf
      real xk(MAXKNOT,MAXEPOCH),yk(MAXKNOT,MAXEPOCH)
      real temp(MAXLOG,MAXTEMP)
      integer nknot(MAXTEMP)                     
      real xknot(MAXKNOT,MAXTEMP),yknot(MAXKNOT,MAXTEMP)
      real delta(MAXTEMP)
      character*100 ftemp,tmpfout
      character*12 tname
      character*10 ttype      
      integer ios

c Functions
      integer iargc,lnb

      if(verbose.gt.0) write(6,'(a,$)') ' Reading template files...'

c
c Return here until we've visited all the template files
c
      ntemp = 0
      ktemp = noption + 1
 1    continue

      ktemp = ktemp + 1

      ios = 0
      if(iargc().lt.noption+2) then
         open(unit=1,file=tempdir(:lnb(tempdir))//'templist',
     $        status='old')
         do i=1, ktemp-noption-1
            read(1,*,end=90,iostat=ios) ftemp
         enddo
 90      close(1)
         if(ios.lt.0) goto 999
         ftemp = tempdir(:lnb(tempdir))//ftemp
      else
         if(ktemp.gt.iargc()) goto 999
         call getarg(ktemp,ftemp)
      endif

* Skip template if avoided or not explicitely used     
      call getroot(ftemp,tmpfout,tname)

      if(nuse.gt.0) then
         inouse = 0
         do ii=1, nuse
            if (index(tname(:lnb(tname)),
     $           use(ii)(:lnb(use(ii)))).eq.0) inouse = inouse + 1
         end do
         if(inouse.eq.nuse) then
            close(1)
            goto 1
         end if
      end if

      if(navoid.gt.0) then
         inouse = 0
         do ii=1, navoid
            if(index(tname(:lnb(tname)),
     $           avoid(ii)(:lnb(avoid(ii)))).gt.0) inouse = inouse + 1  
         end do
         if(inouse.gt.0) then
            close(1)
            goto 1
         end if
      end if

* Open the template file
      open(unit=1,file=ftemp,status='old')

* Read header information
      read(1,1102) nepoch,nwx,w0x,w1x,mostknots,tname,dta,ttype,ittype,
     $     itstype
 1102 format(2i5,2f10.2,i7,5x,a12,f7.2,2x,a10,2(i3))
      
* Skip template if not of required type
      if(nusetype.gt.0) then
         if(usetype(ittype,itstype).eq.0) then
            close (1)
            goto 1
         end if
      end if

      if(navoidtype.gt.0) then
         if(avoidtype(ittype,itstype).eq.1) then
            close (1)
            goto 1
         end if
      end if

      if(nw.ne.nwx.or.w0.ne.w0x.or.w1.ne.w1x) then
         write(6,*) ' '
         write(6,*) 'GETTEMP: Sorry, cannot proceed.  This template'
         write(6,*) '(',ftemp(:lnb(ftemp)),') has a different'
         write(6,*) 'wavelength scale than the first (standard):'
         write(6,*) 'Standard:', nw, w0, w1
         write(6,*) 'Template:', nwx, w0x, w1x
      end if

      if(verbose.gt.1) then
         write(6,'(2a,$)') tname(:lnb(tname)),' '
         call flush(6)
      end if

* Read spline information
      read(1,*) i, (nk(j),fmean,j=1,nepoch)
      do i=1, mostknots
         read(1,*) k, (xk(i,j),yk(i,j),j=1,nepoch)
      end do

* Read age flag and ages
      read(1,*) tf, (ep(i),i=1,nepoch)

* Skip template if out of delta range (except for NotSN templates)
      if((dta.lt.deltamin.or.dta.gt.deltamax)
     $     .and.index(typename(ittype,1),'NotSN').eq.0) goto 1

      n = ntemp
      do 10 j=1, nepoch
* Skip template if out of age range (except for NotSN templates and tflag>=1)
         if((ep(j).lt.agemin.or.ep(j).gt.agemax).and.tf.ne.1
     $        .and.index(typename(ittype,1),'NotSN').eq.0) goto 10
         n = n + 1
         if(n.gt.MAXTEMP) then
            write(6,*) ' '
            write(6,*) 'ERROR! Too many templates!'
            write(6,*) 'Allocate more space for MAXTEMP (in snid.inc)'
            write(6,1003) n,MAXTEMP
 1003       format(' (n,MAXTEMP) = (',i5,',',i5,')')
            stop
         end if
         sname(n)  = tname
         stype(n)  = ttype
         itype(1,n)= ittype
         itype(2,n)= itstype
         epoch(n)  = ep(j)
         tflag(n)  = tf
         delta(n)  = dta
         nknot(n)  = nk(j)
         do 11 i=1, mostknots
            xknot(i,n) = 10.0 ** xk(i,j)
            yknot(i,n) = yk(i,j)
 11      continue
 10   continue

* Skip template if n < nminspec
      if(n.lt.nminspec) goto 1

* Read the normalized spectra
      do k=1, nw
         read(1,*) ww, (tmp(j),j=1,nepoch)
         n = ntemp
         do 21 j=1, nepoch
* Skip template if out of age range (except for NotSN templates and tflag=1)
            if((ep(j).lt.agemin.or.ep(j).gt.agemax).and.tf.ne.1
     $           .and.index(typename(ittype,1),'NotSN').eq.0) goto 21
            n = n + 1
            temp(k,n) = tmp(j)
 21      continue
      end do
      close(1)
      ntemp = n

      goto 1

 999  if(verbose.gt.0) then
         if(ntemp.gt.0) then
            write(6,*) 'done'
         else
            write(6,*) ' '
         end if
         write(6,*) 'Loaded ',ntemp,' spectra out of ',ktemp-noption-2,
     $        ' templates'
      end if
      
      if(ntemp.gt.0) then
* Issue error if attempting to use non-existing templates
         if(nuse.gt.0) then
            if(iargc().lt.noption+2) then
               do i=1, nuse
                  open(unit=1,file=tempdir(:lnb(tempdir))//'templist',
     $                 status='old')
                  inouse = 0
                  do j=1, ktemp-noption-2
                     read(1,*,end=100) ftemp
                     call getroot(ftemp,tmpfout,tname)
                     if(index(tname(:lnb(tname)),
     $                    use(i)(:lnb(use(i)))).eq.0) inouse = inouse+1
                  end do
 100              close(1)
                  if(inouse.eq.ktemp-noption-2) then
                     write(6,*) 'ERROR! Forcing use of',
     $                    ' non-existing template: ',
     $                    use(i)(:lnb(use(i)))
                     stop
                  end if
               end do
            end if
         end if
* Issue error if attempting to avoid non-existing templates
         if(navoid.gt.0) then
            if(iargc().lt.noption+2) then
               do i=1, navoid
                  open(unit=1,file=tempdir(:lnb(tempdir))//'templist',
     $                 status='old')
                  inouse = 0
                  do j=1, ktemp-noption-2
                     read(1,*,end=110) ftemp
                     call getroot(ftemp,tmpfout,tname)
                     if(index(tname(:lnb(tname)),
     $                    avoid(i)(:lnb(avoid(i)))).eq.0)inouse=inouse+1
                  end do
 110              close(1)
                  if(inouse.eq.ktemp-noption-2) then
                     write(6,*) 'ERROR! Attempt to avoid',
     $                    ' non-existing template: ',
     $                    avoid(i)(:lnb(avoid(i)))
                     stop
                  end if
               end do
            end if
         end if
* Stop here if no template spectra were loaded
      else
         stop
      end if

      return
      end



************************************************************************
* subroutine docorr -- Do the cross correlation
************************************************************************

      subroutine docorr(data,dft,temp,tft,cfn,cft,
     $     ctr,hgt,drms,trms,arms,srms,r,redshift,wid)

      implicit none

      include 'snid.inc'

      integer i
      real data(MAXLOG),temp(MAXLOG),cfn(MAXLOG)
      complex dft(MAXLOG),tft(MAXLOG),cft(MAXLOG)
      real*8 pc(5)
      integer ierr
      logical error
      real trms,drms,arms,srms
      real ctr,hgt,wid
      real r,redshift

c Functions
      real rmsfilter

* Initialize variables
      do i=1, 5
         pc(i) = 0.0d0
      end do

* Data transform
      call rcvector(nw,data,dft)
C 030915 JT fft2c is too slow, but need to scale by N to match old behavior
C      call fft2c(dft, log2n, +1)
      call four2(dft,nw,1,+1,1,ierr)
      call cvscale(nw,dft,float(nw))
      drms = rmsfilter(nw,k1,k2,k3,k4,dft)

* Template transform
      call rcvector(nw,temp,tft)
C 030915 JT fft2c is too slow, but need to scale by N to match old behavior
C      call fft2c(tft, log2n, +1)
      call four2(tft,nw,1,+1,1,ierr)
      call cvscale(nw,tft,float(nw))
      trms = rmsfilter(nw,k1,k2,k3,k4,tft)

* Cross-correlation between data and template transforms
      call correlate(nw,dft,tft,cft)
      call ccvector(nw,cft,cfn)
      call filter(nw,k1,k2,k3,k4,cfn)
C 030915 JT fft2c is too slow, but need to scale by N to match old behavior
C      call fft2c(cfn, log2n, -1)
      call four2(cfn,nw,1,-1,1,ierr)
      call cvscale(nw,cfn,1.0/float(nw))
      call crvector(nw,cfn,cfn)

* Fit correlation peak to determine center (ctr), height (hgt), and width (wid)
      call peakfit(nw,lz1,lz2,cfn,ctr,hgt,wid,pc,error)
      hgt = hgt / (drms*trms*nw)

* RMS of (anti)symmetric components of correlation function
      call aspart(nw,k1,k2,k3,k4,ctr,cft,arms,srms)
      arms = arms / (drms*trms*nw)
      srms = srms / (drms*trms*nw)

* Compute the correlation height-noise ratio (r-value)
C      r = 1000
C      if(arms.gt.0) r = hgt / (2*arms)
      if(arms.gt.0) then
         r = hgt / (2*arms)
      else if(arms.eq.0) then
         r = MAXR
      else
         write(6,*) 'ERROR! negative antisymmetric rms'
         write(6,*) 'DOCORR: arms = ', arms
         write(6,*) 'Setting r=0'
         r = 0.0
      end if

* Save away peak redshift and width
C      if(ctr.gt.0.9*nw) ctr = ctr - nw
      if(ctr.gt.(1.-2.*0.01*percent)*nw) ctr = ctr - nw
      redshift = exp(ctr*dwlog) - 1.
      wid      = exp(wid*dwlog) - 1.

      return
      end



************************************************************************
* subroutine allpeaks -- Look at the 10 highest correlation function peaks
************************************************************************

      subroutine allpeaks(cfn,scale,cft,npeak,peaks)

      implicit none

      include 'snid.inc'

      integer i,j,k
      real cfn(0:nw-1)
      complex cft(nw)
      integer ipeak(MAXPEAK)
      integer iz2
      integer npeak
      real cmax
      integer imax
      real ctr
      real arms,srms,scale
      real hgt,r,redshift

      iz2 = nw
      if(lz1.lt.0) iz2 = nw + lz1

      npeak = 0

* Loop over 10 highest correlation peaks
      do 10 k=1, 10
         cmax = 0
         imax = 0
C         do 900 i=1, nw-2
         do 900 i=0, nw-1
            if(i.lt.lz1.or.(i.gt.lz2.and.i.lt.iz2)) goto 900
            if(i.eq.0) then
               if(cfn(i).ge.cfn(nw-1).and.cfn(i).ge.cfn(1)) then
                  do j=1, npeak
                     if(i.eq.ipeak(j)) goto 900
                  end do
                  if(cfn(i).gt.cmax) then
                     cmax = cfn(i)
                     imax = i
                  end if
               end if
            else if (i.eq.nw-1) then
               if(cfn(i).ge.cfn(i-1).and.cfn(i).ge.cfn(0)) then
                  do j=1, npeak
                     if(i.eq.ipeak(j)) goto 900
                  end do
                  if(cfn(i).gt.cmax) then
                     cmax = cfn(i)
                     imax = i
                  end if
               end if                  
            else
               if(cfn(i).ge.cfn(i-1).and.cfn(i).ge.cfn(i+1)) then
                  do j=1, npeak
                     if(i.eq.ipeak(j)) goto 900
                  end do
                  if(cfn(i).gt.cmax) then
                     cmax = cfn(i)
                     imax = i
                  end if
               end if
            end if
 900     continue
         if(cmax.eq.0) goto 10
         npeak = npeak + 1
         ipeak(npeak) = imax
         ctr = imax
         call aspart(nw,k1,k2,k3,k4,ctr,cft,arms,srms)
         arms = arms / scale
         hgt = cmax / scale
         if(arms.gt.0) then
            r = hgt / (2*arms)
         else if(arms.eq.0) then
            r = MAXR
         else
            write(6,*) 'ERROR! negative antisymmetric rms'
            write(6,*) 'PEAKFIT: peak No. ',k,'arms = ',arms
            write(6,*) 'Setting r=0'
            r = 0.0
         end if
         redshift = exp(ctr*dwlog) - 1
         peaks(1,npeak) = redshift
         peaks(2,npeak) = r
         peaks(3,npeak) = imax
C         if(imax.gt.0.9*nw) peaks(3,npeak) = imax - nw
         if(imax.gt.(1.-2.*0.01*percent)*nw) peaks(3,npeak) = imax - nw

 10   continue

      return
      end



************************************************************************
* subroutine getzt -- Compute average/median redshift/age
************************************************************************

      subroutine getzt(nsave,save,sntype,idxa,idx,idxb,ngood,nbad,
     $     ncut,zave,zmed,zsdev,agea,agem,agesdev,ntype,zavetype,
     $     zmedtype,zsdevtype,ageatype,agemtype,agesdevtype,ztemp,
     $     ztemperr,agetemp)

      implicit none

      include 'snid.inc'

      integer k
      integer itemp                                              ! loop index for templates
      integer ibuf                                               ! loop index for number of templates in each type/subtype
      integer nsave                                              ! save away if rlap > rlapmin and lap>lapmin
      real save(9,MAXTEMP)                                       ! array of saved values for statistics and later analysis
      integer sntype(2,MAXTEMP)                                  ! type/subtype indeces
      integer idxa(MAXTEMP)                                      ! sorted index for all templates (based on rlap value)
      integer idx(MAXTEMP)                                       ! sorted index for good templates (based on rlap value and abs(z-zuser) < zfilter)
      integer idxb(MAXTEMP)                                      ! sorted index for bad templates (based on rlap value and abs(z-zuser) >= zfilter)
      integer ngood,nbad                                         ! number of "good" (rlap>rlapmin AND |z-zmed|<zfilter) and "bad" (rlap>rlapmin AND |z-zmed|>zfilter) correlations
      integer ncut,ncuta                                         ! number of correlations with rlap > rlapmin (ncuta for age flag = 0 only)
      integer ntype(NT,NST), ntypea(NT,NST)                      ! total number of SN in each type/subtype (ntypea for age flag = 0 only)
      real zmed,zave,zsdev                                       ! redshift median,weighted-mean,stddev for all good correlations
      real zmedtype(NT,NST),zavetype(NT,NST),zsdevtype(NT,NST)   ! ... and for each type/sub-type
      real ztemp(MAXTEMP,NT,NST),ztemperr(MAXTEMP,NT,NST)        ! redshift and formal error of ranked templates
      real agem,agea,agesdev                                     ! age median,weighted-mean,stddev for all good correlations
      real agemtype(NT,NST),ageatype(NT,NST),agesdevtype(NT,NST) ! ... and for each type/sub-type
      real agetemp(MAXTEMP,NT,NST)                               ! age of ranked templates
      integer idxtype(MAXTEMP,NT,NST)                            ! ranked SNe in each type/subtype

c Local variables
      real sumrl,sumrla                                   ! sum of rlap for all good correlations (sumrla for age flag = 0 only)
      real sumrltype(NT,NST),sumrltypea(NT,NST)           ! sum of rlap for each type/sub-type (sumrltypea for age flag = 0 only)
      integer nbuf,nadd                                   ! number of elements in buf(); nadd is to compute weighted median
      real buf(MAXPPT),buf1(MAXPPT)                       ! redshift/age buffers
      real buftype(NT,NST,MAXPPT),buf1type(NT,NST,MAXPPT) ! redshift/age buffers for each type/sub-type
      real tmpbuf(MAXPPT),tmpbuf1(MAXPPT)                 ! temporary buffers for type/sub-type statistics
      real hoser                                          ! dummy for median and rlap sorting

c Functions
      real amedian,amedidx
      real stddev

* Sort all templates by decreasing rlap
      do k=1, nsave
         buf(k) = -save(1,k)*save(2,k) ! to sort in decreasing order
         idxa(k) = k
      end do
      hoser = amedidx(nsave,buf,idxa)

* Sort only good templates (abs(z-zuser)<zfilter AND rlap>=rlapmin)
      ngood = 0
      nbad = 0
      do k=1, nsave
         if(save(1,idxa(k))*save(2,idxa(k)).ge.rlapmin) then
            if(abs(save(3,idxa(k))-zuser).lt.zfilter) then
               ngood = ngood + 1
               idx(ngood) = idxa(k)
            else
               nbad = nbad + 1
               idxb(nbad) = idxa(k)
            end if
         end if
      end do

* Initialize variables
      ncut = 0
      ncuta = 0
      zave = 0.
      zmed = 0.
      zsdev = 0.
      agea = 0.
      agem = 0.
      agesdev = 0.
      sumrl = 0.
      sumrla = 0.
      do it=1, NT
         do ist=1, NST
            ntype(it,ist)      = 0
            ntypea(it,ist)     = 0
            zavetype(it,ist)   = 0.
            zmedtype(it,ist)   = 0.
            zsdevtype(it,ist)  = 0.
            ageatype(it,ist)   = 0.
            agemtype(it,ist)   = 0.
            agesdevtype(it,ist)= 0.
            sumrltype(it,ist)  = 0.
            sumrltypea(it,ist) = 0.
            do itemp=1, nsave
               idxtype(itemp,it,ist) = 0
            end do
         end do
      end do

* Set up buffers for redshift/age according to a window in dz (zfilter)
      do 10 k=1, ngood

         if(save(1,idx(k))*save(2,idx(k)).ge.rlapmin) then
            ncut = ncut + 1
            buf(ncut) = save(3,idx(k))
            zave = zave + save(3,idx(k))*save(1,idx(k))*save(2,idx(k))
            sumrl = sumrl + save(1,idx(k))*save(2,idx(k))
* Age stats for age flag = 0 only
            if(nint(save(5,idx(k))).eq.0) then
               ncuta = ncuta + 1
               buf1(ncuta) = save(4,idx(k))
               agea = agea+save(4,idx(k))*save(1,idx(k))*save(2,idx(k))
               sumrla = sumrla + save(1,idx(k))*save(2,idx(k))
            end if

            do it=1, NT
               if(sntype(1,idx(k)).eq.it) then 
                  ntype(it,1) = ntype(it,1) + 1
                  buftype(it,1,ntype(it,1)) = save(3,idx(k))
                  zavetype(it,1) = zavetype(it,1) + 
     $                 save(3,idx(k))*save(1,idx(k))*save(2,idx(k))
                  sumrltype(it,1) = sumrltype(it,1) +
     $                 save(1,idx(k))*save(2,idx(k))
* Age stats for age flag = 0 only
                  if(nint(save(5,idx(k))).eq.0) then
                     ntypea(it,1) = ntypea(it,1) + 1                     
                     buf1type(it,1,ntypea(it,1)) = save(4,idx(k))
                     ageatype(it,1) = ageatype(it,1) + 
     $                    save(4,idx(k))*save(1,idx(k))*save(2,idx(k))
                     sumrltypea(it,1) = sumrltypea(it,1) +
     $                    save(1,idx(k))*save(2,idx(k))
                  end if
* idxtype(<rank number>,<type index>,<subtype index>)
                  idxtype(ntype(it,1),it,1) = k
                  ztemp(ntype(it,1),it,1)   = save(3,idx(k))
                  agetemp(ntype(it,1),it,1) = save(4,idx(k))
                  ztemperr(ntype(it,1),it,1)= save(9,idx(k))

                  do ist=2, NST
                     if(sntype(2,idx(k)).eq.ist) then
                        ntype(it,ist) = ntype(it,ist) + 1
                        buftype(it,ist,ntype(it,ist)) = save(3,idx(k))
                        zavetype(it,ist) = zavetype(it,ist) + 
     $                      save(3,idx(k))*save(1,idx(k))*save(2,idx(k))
                        sumrltype(it,ist) = sumrltype(it,ist) +
     $                      save(1,idx(k))*save(2,idx(k))
* Age stats for age flag = 0 only
                        if(nint(save(5,idx(k))).eq.0) then
                           ntypea(it,ist) = ntypea(it,ist) + 1
                           buf1type(it,ist,ntypea(it,ist)) = 
     $                          save(4,idx(k))
                           ageatype(it,ist) = ageatype(it,ist) + 
     $                          save(4,idx(k))*save(1,idx(k))*
     $                          save(2,idx(k))
                           sumrltypea(it,ist) = sumrltypea(it,ist) +
     $                          save(1,idx(k))*save(2,idx(k))
                        end if
* idxtype(<rank number>,<type index>,<subtype index>)
                        idxtype(ntype(it,ist),it,ist) = k
                        ztemp(ntype(it,ist),it,ist)   = save(3,idx(k))
                        agetemp(ntype(it,ist),it,ist) = save(4,idx(k))
                        ztemperr(ntype(it,ist),it,ist)= save(9,idx(k))
                        goto 10
                     end if
                  end do
               end if
            end do
         end if
 10   continue

* Compute median, rlap-weighted mean, and stddev redshift/age
      if(ncut.gt.0.and.sumrl.gt.0) then
         zmed = amedian(ncut,buf)
         zave = zave / sumrl
         if(ncut.gt.1) zsdev = stddev(ncut,buf) 
      end if

      if(ncuta.gt.0.and.sumrla.gt.0) then
         agem = amedian(ncuta,buf1)
         agea = agea / sumrla
         if(ncuta.gt.1) agesdev = stddev(ncuta,buf1)
      end if

* Same for each type/subtype
      do it=1, NT
         do ist=1, NST

            if(ntype(it,ist).gt.0.and.sumrltype(it,ist).gt.0) then
               do ibuf=1, ntype(it,ist)
                  tmpbuf(ibuf) = buftype(it,ist,ibuf)
               end do
               zmedtype(it,ist) = amedian(ntype(it,ist),tmpbuf)
               zavetype(it,ist) = zavetype(it,ist) / sumrltype(it,ist)
               if(ntype(it,ist).gt.1)  zsdevtype(it,ist) = 
     $              stddev(ntype(it,ist),tmpbuf)
            end if

            if(ntypea(it,ist).gt.0.and.sumrltypea(it,ist).gt.0) then
               do ibuf=1, ntypea(it,ist)
                  tmpbuf1(ibuf)= buf1type(it,ist,ibuf)
               end do
               agemtype(it,ist) = amedian(ntypea(it,ist),tmpbuf1)
               ageatype(it,ist) = ageatype(it,ist) / sumrltypea(it,ist)
               if(ntypea(it,ist).gt.1) agesdevtype(it,ist) = 
     $              stddev(ntypea(it,ist),tmpbuf1)
            end if

         end do
      end do

      return
      end



************************************************************************
* subroutine wfout -- Write generic output file
************************************************************************

      subroutine wfout(zinit,zmed,zsdev,agem,agesdev,
     $     ntype,fractype,slope,zmedtype,zsdevtype,agemtype,agesdevtype,
     $     nsave,ngood,nbad,save,idxa,sid,sntype)

      
      implicit none

      include 'snid.inc'

      integer i,l,l0
      real zinit                                ! initial redshift estimate (median)
      real zmed,zsdev                           ! redshift median,weighted-mean,stddev for all good correlations
      real agem,agesdev                         ! age median,weighted-mean,stddev for all good correlations
      integer ntype(NT,NST)                     ! total number of SN in each (sub)type
      real fractype(NT,NST)                     ! absolute fraction of SN in each (sub)type
      real slope(NT,NST)                        ! slope of fracrlap(rlaparr)
      real zmedtype(NT,NST),zsdevtype(NT,NST)   ! ... same for each (sub)type
      real agemtype(NT,NST),agesdevtype(NT,NST) ! ... same for each (sub)type
      integer nsave                             ! save away if lap>=lapmin
      integer ngood,nbad                        ! number of good and bad correlations
      real rlap                                 ! for displaying rlap value
      real save(9,MAXTEMP)                      ! array of saved values for statistics and later analysis
      integer idxa(MAXTEMP)                     ! sorted index for all templates (based on rlap value)
      character*12 sid(MAXTEMP)                 ! name of SN template
      integer sntype(2,MAXTEMP)                 ! (sub)type indeces      
      character*200 file                        ! name of file for output

c Functions
      integer lnb

* Open output file
      file = froot(:lnb(froot))//'_snid.output'
      open(unit=1,file=file,status='replace')
      
* Header
      write(1,'(a)') '### SNID output file ###'

* Input spectrum and options
c input
      write(1,'(a)') ' '
      write(1,'(a)') '### input spectrum and options ###'
      write(1,'(2a)') 
     $     '# Input spectrum                 : ',fdata(:lnb(fdata))
c paramfile
      if(iparam.gt.0) then
         write(1,'(2a)') 
     $     '# Parameter file                 : ',
     $        paramfile(:lnb(paramfile))
      else
         write(1,'(a)') 
     $     '# Parameter file                 : NA'           
      end if
c wmin,wmax
      write(1,'(a,2f10.2)') 
     $     '# Wavelength range               : ',wmin,wmax
c wmaskfile
      if(iwmask.gt.0) then
         write(1,'(2a)') 
     $     '# Mask file                      : ',wmaskfile
      else
         write(1,'(a)') 
     $     '# Mask file                      : NA'
      end if         
c fwmed
      if(fwmed.gt.0) then
         write(1,'(a,f10.2)')
     $     '# Median filtering width [A]     : ',fwmed
      else
         write(1,'(a)')
     $     '# Median filtering width [A]     : NA'
      end if
c medlen
      if(medlen.gt.0) then
         write(1,'(a,f10.2)')
     $     '# Median filtering width [pix]   : ',medlen
      else
         write(1,'(a)')
     $     '# Median filtering width [pix]   : NA'
      end if
c emclipz
      if(emclipz.gt.0) then
         write(1,'(a,f10.4)')
     $     '# Clip emission redshift         : ',emclipz
      else
         write(1,'(a)')
     $     '# Clip emission redshift         : NA'
      end if
c skyclip
      write(1,'(a,i1)') 
     $     '# Clip sky lines                 : ',skyclip
c emwid
      if(emclipz.gt.0.or.skyclip.gt.0) then
         write(1,'(a,f10.2)')
     $     '# Emission clip width [A]        : ',emclipz            
      else
         write(1,'(a)')
     $     '# Emission clip width [A]        : NA'
      end if
c aband
      write(1,'(a,i1)') 
     $     '# Keep A-band                    : ',iaband
c agemin,agemax
      write(1,'(a,2f10.2)') 
     $     '# Age range                      : ',agemin,agemax
c deltamin,deltamax
      write(1,'(a,2f10.2)') 
     $     '# Delta range                    : ',deltamin,deltamax
c nminspec
      write(1,'(a,i3)') 
     $     '# Min. number of spectra/template: ',nminspec
c use      
      if(nuse.gt.0) then
         write(1,'(2a)')
     $     '# Use only templates             : ',uses(:lnb(uses))
      else
         write(1,'(a)')
     $     '# Use only templates             : NA'
      end if
c usetype      
      if(nusetype.gt.0) then
         write(1,'(2a)')
     $     '# Use only types                 : ',usets(:lnb(usets))
      else
         write(1,'(a)')
     $     '# Use only types                 : NA'
      end if
c avoid      
      if(navoid.gt.0) then
         write(1,'(2a)')
     $     '# Avoid templates                : ',avoids(:lnb(avoids))
      else
         write(1,'(a)')
     $     '# Avoid templates                : NA'
      end if
c avoidtype      
      if(navoidtype.gt.0) then
         write(1,'(2a)')
     $     '# Avoid types                    : ',avoidts(:lnb(avoidts))
      else
         write(1,'(a)')
     $     '# Avoid types                    : NA'
      end if
c rlapmin
      write(1,'(a,f10.2)') 
     $     '# rlapmin                        : ',rlapmin      
c zmin,zmax
      write(1,'(a,2f10.4)') 
     $     '# Redshift range                 : ',zmin,zmax
c forcez
      if(forcez.gt.-1) then
         write(1,'(a,f10.4)') 
     $     '# Forced initial redshift        : ',forcez               
      else
         write(1,'(a)') 
     $     '# Forced initial redshift        : NA'
      end if
c zfake
      if(zfake.gt.0) then
         write(1,'(a,f10.4)') 
     $     '# Extra fake redshift            : ',zfake               
      else
         write(1,'(a)') 
     $     '# Extra fake redshift            : NA'
      end if
c znsig
      if(znsig.gt.0) then
         write(1,'(a,f10.2)') 
     $     '# Nsigma for initial redshift    : ',znsig               
      else
         write(1,'(a)') 
     $     '# Nsigma for initial redshift    : NA'
      end if
c zfilter
      write(1,'(a,f10.4)') 
     $     '# Redshift filter                : ',zfilter
c k1,k2,k3,k4
      write(1,'(a,4(i5))') 
     $     '# Wavenumbers k1,k2,k3,k4        : ',k1,k2,k3,k4
c tempdir
      write(1,'(2a)') 
     $     '# Template directory             : ',tempdir(:lnb(tempdir))

* Initial/user-input redshift
      write(1,'(a)') ' '
      write(1,'(a)') '### initial/user-input redshift ###'
      write(1,'(a,f10.4)') 'zinit   ',zinit
      write(1,'(a,f10.4)') 'zuser   ',zuser

* Overall median redshift/age
      write(1,'(a)') ' '
      write(1,'(a)') '### median redshift/age and error ###'
      write(1,'(a,2(f10.4))') 'zmed    ',zmed,zsdev
      write(1,'(a,2(f10.2))') 'agem    ',agem,agesdev

* Type fraction, redshift, and age
      write(1,'(a)') ' '
      write(1,'(a)') '### type fraction/redshift/age ###'
      write(1,'(2a)') '#type ntemp fraction slope',
     $     ' redshift redshift_error age age_error'
      do it=1, NT
         do ist=1, NST
            if(it.gt.1.and.ist.eq.1) write(1,'(a)') '###'
            if(index(typename(it,ist),'typename').eq.0) then
               write(1,'(a10,i7,f10.2,3(f10.4),2(f10.3))') 
     $              typename(it,ist),ntype(it,ist),fractype(it,ist),
     $              slope(it,ist),zmedtype(it,ist),zsdevtype(it,ist),
     $              agemtype(it,ist),agesdevtype(it,ist)
            end if
         enddo
      end do

* Ordered template listings
      write(1,'(a)') ' '
      write(1,'(a)') '### rlap-ordered template listings ###'
      write(1,'(a)') '#no. sn type lap rlap z zerr age age_flag grade'
      do i=1, nsave
         rlap = save(1,idxa(i))*save(2,idxa(i))
         if(rlap.gt.MAXR) rlap = MAXR
         if(rlap.ge.rlapmin) then
            if(abs(save(3,idxa(i))-zuser).lt.zfilter) then
               write(1,1000) i,sid(idxa(i)),
     $              typename(sntype(1,idxa(i)),sntype(2,idxa(i))),
     $              save(2,idxa(i)),rlap,
     $              save(3,idxa(i)),save(9,idxa(i)),save(4,idxa(i)),
     $              nint(save(5,idxa(i))),'   good'
            else
               write(1,1000) i,sid(idxa(i)),
     $              typename(sntype(1,idxa(i)),sntype(2,idxa(i))),
     $              save(2,idxa(i)),rlap,
     $              save(3,idxa(i)),save(9,idxa(i)),save(4,idxa(i)),
     $              nint(save(5,idxa(i))),'   bad'
            end if
         else
            if(i-1.eq.ngood+nbad) write(1,'(a)') 
     $           '#--- rlap cutoff'
            write(1,1000) i,sid(idxa(i)),
     $           typename(sntype(1,idxa(i)),sntype(2,idxa(i))),
     $           save(2,idxa(i)),rlap,
     $           save(3,idxa(i)),save(9,idxa(i)),save(4,idxa(i)),
     $           nint(save(5,idxa(i))),'   cut'               
         end if
      end do
 1000 format(i4,2x,a12,2x,a10,f7.4,f7.2,2(f9.4),f9.2,i3,a7)
      
* Close file
      close(1)
      if(verbose.gt.0) then
         write(6,*) ' '
         write(6,*) 'Created output file: '//froot(:lnb(froot))//
     $        '_snid.output'
      end if

      return
      end



************************************************************************
* subroutine wfluxout -- Write fluxed spectra to output file
************************************************************************

      subroutine wfluxout(itn,z,dfl,dmean,dmeanz,tfl)
      
      implicit none

      include 'snid.inc'

      integer i,l,l0
      integer itn                                               ! template number
      real z                                                    ! template redshift
      real dfl(MAXPPT),dmean(MAXPPT),tfl(MAXPPT),dmeanz(MAXPPT) ! data flux and mean; template flux; de-redshifted data mean 
      character*200 file                                        ! name of file for output
      real dflux(MAXWAVE),tflux(MAXWAVE)                        ! data/template flux with mean included
      real zwave(MAXWAVE)                                       ! redshifted wavelength axis for template

c Functions
      integer lnb

* Data/template flux normalized to data flux
      do i=1, nw                
         dflux(i) = dfl(i) + 1.
         dflux(i) = dflux(i) * dmean(i)
         tflux(i) = tfl(i) + 1.
         tflux(i) = tflux(i) * dmeanz(i)
      end do

* Wavelength axis
      do i=1, nw
         wave(i) = w0 * exp(float(i-1)*dwlog)
         zwave(i)= wave(i) * (1.+z) 
      end do

* File output

c input spectrum -- only once!
      if(itn.eq.1) then
         file = froot(:lnb(froot))//'_snidflux.dat'
         open(unit=1,file=file,status='replace')
         write(1,'(a)') '#wavelength[A] flux[arbitraty]'
         do i=1, nw
            if(wave(i).ge.wmin.and.wave(i).le.wmax) 
     $           write(1,*) wave(i),dflux(i)
         end do
         close(1)
         if(verbose.gt.0) write(6,*) 'Wrote fluxed input spectrum to',
     $        ' file: ',file(:lnb(file))
      end if

c template spectrum
      file = froot(:lnb(froot))//'_comp'//
     $     dlabels(8)(:lnb(dlabels(8)))//
     $     '_snidflux.dat'
      open(unit=1,file=file,status='replace')
      write(1,'(a)') '# Template No.'//
     $     dlabels(1)(:lnb(dlabels(1)))//': '//
     $     dlabels(2)(:lnb(dlabels(2)))//' ('//
     $     dlabels(3)(:lnb(dlabels(3)))//'; '//
     $     dlabels(4)(:lnb(dlabels(4)))//')'//' ; z='//
     $     dlabels(5)(:lnb(dlabels(5)))//' +/- '//
     $     dlabels(6)(:lnb(dlabels(6)))
      write(1,'(a)') '#redshifted_wavelength[A] flux[arbitraty]'
      do i=1, nw
         if(zwave(i).ge.wmin.and.zwave(i).le.wmax) 
     $        write(1,*) zwave(i),tflux(i)
      end do
      close(1)

      return
      end



************************************************************************
* subroutine wflatout -- Write flattened spectra to output file
************************************************************************

      subroutine wflatout(itn,z,dfl,tfl)
      
      implicit none

      include 'snid.inc'

      integer i,l,l0
      integer itn                        ! template number
      real z                             ! template redshift
      real dfl(MAXPPT),tfl(MAXPPT)       ! data and template flux
      character*200 file                 ! name of file for output
      real dflux(MAXWAVE),tflux(MAXWAVE) ! data/template flux with mean included
      real zwave(MAXWAVE)                ! redshifted wavelength axis for template

c Functions
      integer lnb

* Wavelength axis
      do i=1, nw
         wave(i) = w0 * exp(float(i-1)*dwlog)
         zwave(i)= wave(i) * (1.+z) 
      end do

* File output

c input spectrum -- only once!
      if(itn.eq.1) then
         file = froot(:lnb(froot))//'_snidflat.dat'
         open(unit=1,file=file,status='replace')
         write(1,'(a)') '#wavelength[A] flux[arbitraty]'
         do i=1, nw
            if(wave(i).ge.wmin.and.wave(i).le.wmax) 
     $           write(1,*) wave(i),dfl(i)
         end do
         close(1)
         if(verbose.gt.0) write(6,*) 'Wrote flattened input spectrum',
     $        ' to file: ',file(:lnb(file))
      end if

c template spectrum
      file = froot(:lnb(froot))//'_comp'//
     $     dlabels(8)(:lnb(dlabels(8)))//
     $     '_snidflat.dat'
      open(unit=1,file=file,status='replace')
      write(1,'(a)') '# Template No.'//
     $     dlabels(1)(:lnb(dlabels(1)))//': '//
     $     dlabels(2)(:lnb(dlabels(2)))//' ('//
     $     dlabels(3)(:lnb(dlabels(3)))//'; '//
     $     dlabels(4)(:lnb(dlabels(4)))//')'//' ; z='//
     $     dlabels(5)(:lnb(dlabels(5)))//' +/- '//
     $     dlabels(6)(:lnb(dlabels(6)))
      write(1,'(a)') '#redshifted_wavelength[A] flux[arbitraty]'
      do i=1, nw
         if(zwave(i).ge.wmin.and.zwave(i).le.wmax) 
     $        write(1,*) zwave(i),tfl(i)
      end do
      close(1)

      return
      end



************************************************************************
* subroutine wxcorout -- Write correlation functions to output file
************************************************************************

      subroutine wxcorout(itn,npeak,zpeak,cfn,cfntrim)
      
      implicit none

      include 'snid.inc'

      integer i,l,l0
      integer itn                          ! template number
      integer npeak                        ! total number of correlation peaks
      real zpeak(MAXPPT)                   ! z of peak
      real cfn(2*MAXLOG),cfntrim(2*MAXLOG) ! correlation function and trimmed version
      character*200 file,file2             ! name of file for output
      real zaxis(MAXWAVE)                  ! redshifted axis
      
c Functions
      integer lnb

* Redshift axis
      do i=1, nw
         zaxis(i) = exp(i*dwlog) - 1.
      end do

* File output
      file = froot(:lnb(froot))//'_comp'//
     $     dlabels(8)(:lnb(dlabels(8)))//'_snidxcor.dat'
      file2 = froot(:lnb(froot))//'_comp'//
     $     dlabels(8)(:lnb(dlabels(8)))//'_snidxcor_trim.dat'

c correlation function
      open(unit=1,file=file,status='replace')
      write(1,'(a)') '# Template No.'//
     $     dlabels(1)(:lnb(dlabels(1)))//': '//
     $     froot(:lnb(froot))//' x '//
     $     dlabels(2)(:lnb(dlabels(2)))//' ('//
     $     dlabels(3)(:lnb(dlabels(3)))//'; '//
     $     dlabels(4)(:lnb(dlabels(4)))//')'
      write(1,'(a)') '#redshift correlation'
      do i=1, nw
         if(zaxis(i).ge.zmin.and.zaxis(i).le.zmax) 
     $        write(1,*) zaxis(i),cfn(i)
      end do
      close(1)

c trimmed correlation function
      open(unit=1,file=file2,status='replace')
      write(1,'(a)') '# Template No.'//
     $     dlabels(1)(:lnb(dlabels(1)))//': '//
     $     froot(:lnb(froot))//' x '//
     $     dlabels(2)(:lnb(dlabels(2)))//' ('//
     $     dlabels(3)(:lnb(dlabels(3)))//'; '//
     $     dlabels(4)(:lnb(dlabels(4)))//')'
      write(1,'(a)') '#redshift correlation'
      do i=1, nw
         if(zaxis(i).ge.zmin.and.zaxis(i).le.zmax) 
     $        write(1,*) zaxis(i),cfntrim(i)
      end do
      close(1)

      return
      end



************************************************************************
* subroutine loadpsec -- Load up the plot buffer with data and template, 
*                        raw and smooth
************************************************************************

      subroutine loadspec(data,dft,temp,tft,pdat,idr,ids,itr,its)

      implicit none

      include 'snid.inc'

      real data(MAXLOG),temp(MAXLOG)
      complex dft(MAXLOG),tft(MAXLOG)
      real pdat(MAXPPT,1)
      complex ftbuf(MAXLOG)
      integer idr,itr,ids,its
      integer ierr

      call rrvector(nw,data,pdat(1,idr))
      call rrvector(nw,temp,pdat(1,itr))

      call ccvector(nw,dft,ftbuf)
      call filter(nw,k1,k2,k3,k4,ftbuf)
C 030915 JT fft2c is too slow, but need to scale by N to match old behavior
C      call fft2c(ftbuf,log2n,-1)
      call four2(ftbuf,nw,1,-1,1,ierr)
      call cvscale(nw,ftbuf,1.0/float(nw))
      call crvector(nw,ftbuf,pdat(1,ids))

      call ccvector(nw,tft,ftbuf)
      call filter(nw,k1,k2,k3,k4,ftbuf)
C 030915 JT fft2c is too slow, but need to scale by N to match old behavior
C      call fft2c(ftbuf,log2n,-1)
      call four2(ftbuf,nw,1,-1,1,ierr)
      call cvscale(nw,ftbuf,1.0/float(nw))
      call crvector(nw,ftbuf,pdat(1,its))

      return
      end
