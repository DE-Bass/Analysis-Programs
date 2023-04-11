************************************************************************
**                                                                    **
**                                                                    **
**           S u p e r N o v a    I D e n t i f i c a t i o n         **
**                                                                    **
**                                                                    **
**                              ( S N I D )                           **
**                                                                    **
**                                                                    **
**     Program to compare a data spectrum against template spectra    **
**                                                                    **
**                                                                    **
**     Copyright (C) 1999-2007  Stephane Blondin and John L. Tonry    ** 
**                                                                    **
**                                                                    **
**  v1.0-4.2 : Aug 1999 - Sep 2003 (John Tonry)                       **
**  v5.0     : Apr 2006 - Aug 2007 (Stephane Blondin)                 **
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
* 24 Aug 2007, rev 5.0 (SB) - New routines, type statistics,  new GUI
*     Major code clean-up
*     New user interface, more options
*     Output various type statistics
*     New GUI using PGPLOT and the BUTTON library (snidplot.f)     
* 24 Aug 2007 (SB)
*     Checked compilation with -fautomatic
*     Fixed bug in getdata: defined z0=0.5*(zmin+zmax)
*     Fixed bug in rebin(): maxwave -> nwave
* 15 Aug 2007 (SB)
*     Fixed small bug in peakfit which ignored z=0 peaks in trimmed xcor
* 24 Apr 2006 (SB)
*     Fixed small bug in avoid option
*     Removed nopeak option (not used)
*     Removed zspline variable in snidmore.f (not used)
*     Restricted default wmin,wmax to match range of input spectrum
*     Added {use,avoid}type options
*     Changed use option to mean "use only"
*     Added {d}delta,delta{min,max} options
*     Added {d}age options to conform to z,delta options
*     Added iquery option
*     Redshift/age for each (sub)type
*     Display formal redshift error
*     IMPLICIT NONE in logwave.f and snid*.f
*     changed variable snid-->sid (for "program snid" statement)
*
* 15 Sep 2003, rev 4.2 (JT) - 
*     Fixed a potential bug in rmsetc with CMPLX(), squashed a memory bug
*     Enhanced median filtering options "fwmed="
*     Changed fft2c() to four2()
*     Added emission line clipping "emclip=" and "emwid="
*     Changed "peaks" to show trim/peak and added "pub" plotting
*     Reworked how user's z is used to select trim
*               
* 25 Aug 2000, rev 4.1 (JT) - Altered control flow a lot:
*     Open correlate, peak search
*     Overlap, best z, stats
*     Select SN and plot
*               
* 21 Aug 2000, rev 4.01 (JT) - Major changes:
*     Modularized code so as to make it more readable
*     Added allpeaks subroutine and "peaks" mongo macro
*     Added some more options such as lapmin and nopeak
*     Cleaned up mongo macros to make displays more informative
*
* 15 Aug 2000, rev 4.0 (JT) - Add AGEMIN, AGEMAX, MEDLEN and a few other things
*
* 09 Nov 1999, rev 3.2.2 (JT) - make ZMIN and ZMAX really restrict peak search
*
* 08 Nov 1999, rev 3.2.1 (JT) - add use option
*
* 20 Oct 1999, rev 3.2 (JT) - add trim+smooth spectra
*
* 19 Oct 1999, rev 3.1 (JT) - remove plotting calc outside of main loop
*
* 16 Oct 1999, rev 3.0 (JT) - add spline info into templates
*
* 09 Oct 1999, rev 2.0 (JT) - add plotting enhancements into snid.{f,pl}
*
************************************************************************

************************************************************************
**                                                                    **
**                            Declarations                            **
**                                                                    **
************************************************************************

      program snid
     
      implicit none
     
      include 'snid.inc'
      include '../button/button.inc'

c Loop indices
      integer i,j,k
      integer itry  ! loop index for correlation peaks

c Spectral info
      character*12 sid(MAXTEMP) ! name of SN template
      integer sntype(2,MAXTEMP) ! (sub)type indeces
      real fnorm(MAXLOG)        ! normalized input spectrum
      real temp(MAXLOG,MAXTEMP) ! array containing all template spectra
      real delta(MAXTEMP)       ! delta value of SN template

c Spectral preparation and cross-correlation
      integer nknot(MAXTEMP)                             ! number of knot points for each template spectrum
      real xknot(MAXKNOT,MAXTEMP),yknot(MAXKNOT,MAXTEMP) ! knot point (w,f) coordinates for each template spectrum
      real ftrim(MAXLOG),ttrim(MAXLOG)                   ! redshift-trimmed inpute/template spectra
      real cfn(2*MAXLOG)                                 ! (filtered) correlation function 
      complex dft(MAXLOG),tft(MAXLOG),cft(MAXLOG)        ! data/template/correlation function fourier transforms
      integer npeak                                      ! total number of correlation peaks
      real rpeak(MAXPPT),lpeak(MAXPPT),zpeak(MAXPPT)     ! r,lap,z of peak
      real jpeak(MAXPPT)                                 ! index of peak
      real ctr,hgt,wid,r                                 ! center,height,width and r-value of best peak
      real drms,trms                                     ! data/template RMS
      real arms,srms                                     ! (anti)symmetric RMS of correlation function
      real redshift                                      ! correlation redshift
      integer npeaktry                                   ! choose best npeaktry peaks (out of 10)
      real ctx                                           ! center of best correlation peaks (to trim at correlation redshift)
      real lap                                           ! data/template overlap at correlation redshift

c General statistics
      integer nsave ! save away if lap>=lapmin
      integer ncut  ! number of correlations with rlap>=rlapmin

c Type statistics
      integer ntype(NT,NST)              ! total number of SN in each (sub)type
      integer idxtype(MAXTEMP,NT,NST)    ! ranked SNe in each (sub)type
      real fractype(NT,NST)              ! absolute fraction of SN in each (sub)type
      real rfractype(NT,NST)             ! relative fraction of SN in each (sub)type
      integer idxa(MAXTEMP)              ! sorted index for all templates (based on rlap value)
      integer idx(MAXTEMP)               ! sorted index for good templates (based on rlap value)
      integer idxb(MAXTEMP)              ! sorted index for bad templates (based on rlap value)
      integer idxt(NT)                   ! sorted index for types (based on absolute fraction) 
      integer idxst(NT*(NST-1))          ! sorted index for subtypes (based on absolute fraction) 
      integer idxstr(NT,NST-1)           ! sorted index for subtypes (based on relative fraction)
      integer idxts(NT)                  ! sorted index for types (based on slope) 
      integer idxsts(NT*(NST-1))         ! sorted index for subtypes (based on slope)
      integer idxnt(NT*(NST-1))          ! index of type name corresponding to each subtype (for absolute fraction)
      integer idxnst(NT*(NST-1))         ! index of subtype name (for absolute fraction)
      integer idxnts(NT*(NST-1))         ! index of type name corresponding to each subtype (for slope)
      integer idxnsts(NT*(NST-1))        ! index of subtype name (for slope)
      integer imbtemp                    ! number of multiple best templates (based on rlap value)
      integer imbtempr(NT,NST)           ! number of multiple best templates in each (sub)type (based on rlap value)
      integer imbt                       ! number of multiple best types (based on absolute fraction) 
      integer imbst                      ! number of multiple best subtypes (based on absolute fraction) 
      integer imbstr(NST)                ! number of multiple best subtypes (based on relative fraction) 
      integer imbts                      ! number of multiple best types (based on slope) 
      integer imbsts                     ! number of multiple best subtypes (based on slope) 
      integer nbt,nbst,nobst             ! number of best (sub)type based on n best templates
      integer samebt,samebst             ! number of best (sub)types
      integer rlapmax                    ! maximum integer rlap value
      real rlaparr(MAXRLAP)              ! rlap array (rlaparr[0:rlapmax])
      real fracrlap(MAXRLAP,NT,NST)      ! absolute fraction of templates with rlap >= rlap_i for each (sub)type
      real fracrlaperr(MAXRLAP,NT,NST)   ! poisson error in absolute fraction of templates with rlap >= rlap_i for each (sub)type
      real slope(NT,NST)                 ! slope of fracrlap(rlaparr)
      integer okt,okst                   ! secure (sub)type indicators
      integer bttemp,bsttemp             ! (sub)type of best-match templates
      integer btfrac,bstfrac,btstfrac    ! (sub)type with highest fraction
      integer btslope,bstslope,btstslope ! (sub)type with highest slope

c Redshift statistics
      real zinit                                               ! initial redshift estimate (median)
      real zmed,zave,zsdev                                     ! redshift median,weighted-mean,stddev for all good correlations
      real zmedtype(NT,NST),zavetype(NT,NST),zsdevtype(NT,NST) ! ... same for each (sub)type
      real ztemp(MAXTEMP,NT,NST),ztemperr(MAXTEMP,NT,NST)      ! redshift and formal error of ranked templates

c Age statistics
      real agem,agea,agesdev                                     ! age median,weighted-mean,stddev for all good correlations
      real agemtype(NT,NST),ageatype(NT,NST),agesdevtype(NT,NST) ! ... same for each (sub)type
      real agetemp(MAXTEMP,NT,NST)                               ! age of ranked templates

c For redshift/age statistics
      real save(9,MAXTEMP) ! array of saved values for statistics and later analysis
      real buf(MAXPPT)     ! redshift buffers
      integer nbuf,nadd    ! number of elements in buf(); nadd is to compute weighted median
      real buf1(MAXPPT)    ! rlap buffer
      integer nbuf1        ! number of elements in buf1()
      real sumbuf1,maxbuf1 ! sum of buf1(); max(buf1)
      real zmaxbuf1        ! redshift corresponding to maxbuf1
      real mbuf1,sbuf1     ! mean/sdev of buf1()

c Display stuff
      character*1 reply    ! interactive user input
      integer ngood,nbad   ! number of good and bad correlations
      integer ngoodt,nbadt ! number of good and bad correlations for ageflag=0
      real rlap            ! for displaying rlap value
      character*200 line   ! length of string to parse
      integer ntok         ! number of words in line 
      integer iexp         ! =1 if age from explosion is known but not that from maximum
      real tmax,texp       ! date of maximum and explosion

c File output
      integer nout,nfluxout,nflatout,nxcorout ! number of output files

c For plotting
      real xc,yc                        ! cursor coordinates
      character*1 ch                    ! character for mouse button
      integer itn                       ! template number
      integer nb                        ! button number
      real pdat(MAXPPT,MAXPLOT)         ! data array for plotting
      real zgood(MAXTEMP),zbad(MAXTEMP) ! good and bad redshifts
      real tgood(MAXTEMP),tbad(MAXTEMP) ! good and bad ages

c Functions
      real amedian,amedidx
      integer lnb,nlb
      real stddev

c IEEE flag
      integer iflag
      character*9 out
      integer ieee_flags

c Version number
      data version /'5.0, 24 Aug 2007'/






************************************************************************
**                                                                    **
**                              Preamble                              **
**                                                                    **
************************************************************************

* Get the template directory
      call gettempdir

* Get type info
      call typeinfo
      
* Interpret the options and initialize the control variables
      call initialize

* Read in the input spectrum
      call getdata(fnorm,nknod,xknod,yknod)

* Read in the template spectra
      call gettemp(ntemp,sname,stype,itype,delta,epoch,tflag,temp,
     $     nknot,xknot,yknot)
      
* Skip first round of correlations if forced redshift
      if(forcez.gt.-1) then
         if(verbose.gt.0) write(6,'(a,f7.3)') 
     $        ' Forcing initial redshift:',forcez
         zuser = forcez
         goto 30
      end if
      





************************************************************************
**                                                                    **
**    First run of cross-correlations -- Initial redshift estimate    **
**                                                                    **
************************************************************************

      if(verbose.gt.0) write(6,'(a,$)') 
     $     ' Searching all correlation peaks...'
      npeak = 0

c
c Loop over all template spectra
c
      do j=1, ntemp

* Compute the correlation function ab initio
         call docorr(fnorm,dft,temp(1,j),tft,cfn,cft,
     $        ctr,hgt,drms,trms,arms,srms,r,redshift,wid)

* Look at the 10 highest correlation function peaks
         call allpeaks(cfn,drms*trms*nw,cft,npeaktry,peaks)

c
c Loop over valid peaks
c
         do 10 itry=1, npeaktry
            ctx = peaks(3,itry)

* Trim input/template spectra at peak redshift
            call overlap(nw,ctx,lap,temp(1,j),ttrim,fnorm,ftrim,percent)

* Do the cross-correlation if the overlap is still sufficient
            if(lap*dwlog.lt.lapmin) goto 10
            call docorr(ftrim,dft,ttrim,tft,cfn,cft,
     $           ctr,hgt,drms,trms,arms,srms,r,redshift,wid)
            if(redshift.lt.zmin.or.redshift.gt.zmax) goto 10
            if(npeak.ge.MAXPPT) then
               write(6,*) ' '
               write(6,*) ' '
               write(6,*) 'WARNING! Buffer for peaks will be full!'
               write(6,*) 'Try using fewer templates or set MAXPPT to',
     $              'a larger value (in snid.inc). MAXPPT = ',MAXPPT
               write(6,*) ' '
               stop
            end if
            npeak = min(MAXPPT,npeak+1)
            pdat(npeak,17) = lap*dwlog
            pdat(npeak,18) = r
            pdat(npeak,19) = redshift
            pdat(npeak,20) = exp(ctx*dwlog) - 1
            jpeak(npeak) = j
            rpeak(npeak) = r
            lpeak(npeak) = lap*dwlog
            zpeak(npeak) = redshift
 10      continue

      end do

c
c Determine the most likely redshift (zinit)
c
      nbuf = 0
      do i=1, npeak
* Assemble a weighted median
         nadd = 0
         if(rpeak(i)*lpeak(i).gt.4) nadd = nadd + 1
         if(rpeak(i)*lpeak(i).gt.5) nadd = nadd + 2
         if(rpeak(i)*lpeak(i).gt.6) nadd = nadd + 2
         if(nbuf+nadd.gt.MAXPPT) then
            write(6,*) ' '
            write(6,*) ' '
            write(6,'(2a)') ' WARNING! Median buffer for z peaks is',
     $           ' full!'
            write(6,*) 'Try using fewer templates or set MAXPPT to',
     $           'a larger value (in snid.inc). MAXPPT = ',MAXPPT
            write(6,*) ' '            
            stop
         end if
         do k=1, nadd
            buf(nbuf+k) = zpeak(i)
         end do
         nbuf = nbuf + nadd
      end do
      zinit = amedian(nbuf,buf)
      if(verbose.gt.0) then
         write(6,*) 'done'
      end if

* Set zinit to highest peak if its rlap is znsig above the rlap distribution
* TEST!
      if(znsig.gt.0) then
         nbuf1 = 0
         sumbuf1 = 0
         maxbuf1 = 0
         zmaxbuf1 = 0
         mbuf1 = 0
         sbuf1 = 0
         do i=1, npeak
            if(rpeak(i)*lpeak(i).gt.4) then
               nbuf1 = nbuf1 + 1
               buf1(nbuf1) = rpeak(i)*lpeak(i)
               sumbuf1 = sumbuf1 + rpeak(i)*lpeak(i)
               if(buf1(nbuf1).gt.maxbuf1) then
                  maxbuf1 = buf1(nbuf1)
                  zmaxbuf1 = zpeak(i)
               end if
            end if
         end do
         if(nbuf1.gt.0) mbuf1 = sumbuf1 / real(nbuf1)
         if(nbuf1.gt.1) sbuf1 = stddev(nbuf1,buf1)
         if(maxbuf1.gt.mbuf1+znsig*sbuf1) then
            zinit = zmaxbuf1
            if(verbose.gt.0) then
               write(6,*) 'Maximum peak rlap is more than ',znsig,
     $              ' sigma above rlap distribution'
               write(6,*) 'Setting initial redshift to that peak...'
            end if
         end if
      end if
 
* Report initial redshift
      if(verbose.gt.0) then
         if(nbuf.gt.0) then
            write(6,'(a,f6.3)') ' Initial redshift estimate: z =', zinit
         else
            write(6,'(a,f6.3)') ' No peaks are good, setting z =', zinit
         end if
      end if


c
c Allow user to impose a different redshift
c
      zuser = zinit
      if (iquery.eq.0.or.inter.eq.0) goto 30

 20   continue
      write(6,'(2(a),$)') ' Do you want to enter a new redshift?',
     $     ' (y/n) [n]: '
      read(5,'(a)',err=20) reply
      if(lnb(reply).eq.0.or.reply.eq.'n') then
         goto 30
      else if(reply.eq.'y') then
 21      continue
         write(6,'(a,$)') ' Enter a new redshift: '
         read(5,*,err=21) zuser
      else
         goto 20
      end if






************************************************************************
**                                                                    **
**   Second run of cross-correlations -- Revised redshift estimate    **
**                                                                    **
************************************************************************

 30   continue

      if(verbose.gt.0) 
     $     write(6,'(a,f6.3)') ' Trimming spectra to match at z =',zuser

      nsave = 0
      ctx = nint(alog(zuser+1)/dwlog)

c
c Loop over all template spectra
c
      do 40 j=1, ntemp

* Trim input/template spectra at user-requested redshift
         call overlap(nw,ctx,lap,temp(1,j),ttrim,fnorm,ftrim,percent)

* Do the cross-correlation if the overlap is still sufficient
         if(lap*dwlog.lt.lapmin) goto 40

* Compute the trimmed correlation function
         call docorr(ftrim,dft,ttrim,tft,cfn,cft,
     $        ctr,hgt,drms,trms,arms,srms,r,redshift,wid)

c
c Save away for statistics and later analysis
c
         nsave = nsave + 1
 
* Template name, type, and subtype    
         sid(nsave)      = sname(j)
         sntype(1,nsave) = itype(1,j)
         sntype(2,nsave) = itype(2,j)

* All other useful parameters
         save(1,nsave) = r
         save(2,nsave) = lap*dwlog
         save(3,nsave) = redshift
         save(4,nsave) = epoch(j)
         save(5,nsave) = tflag(j)
         save(6,nsave) = hgt
         save(7,nsave) = delta(j)
         save(8,nsave) = j
         save(9,nsave) = zk*wid/(1.+r*lap*dwlog)

* Continue with the template iteration
 40   continue

c
c Prompt user for input if no templates had sufficient overlap
c
      if(nsave.eq.0) then
         write(6,'(2(a,f5.3))') ' No template meets lap >= ',lapmin,
     $        ' at redshift z = ',zuser
         if(inter.eq.0) then
            write(6,*) 'Interactive mode off. No output written.'
            stop
         end if
 50      continue
         write(6,'(a,$)') 
     $        ' Enter a new (1) redshift; (2) lapmin; or (q)uit: '
         read(5,'(a)',err=50) reply
         if(reply.eq.'1') then
 51         continue
            write(6,'(a,f6.3,a,$)') 
     $           ' Enter a new redshift (current value ',zuser,'): '
            read(5,*,err=51) zuser
            goto 30
         else if(reply.eq.'2') then
 52         continue
            write(6,'(a,f6.2,a,$)') 
     $           ' Enter a new lapmin (current value ',lapmin,'): ' 
            read(5,*,err=52) lapmin
            goto 30
         else if(reply.eq.'q') then
            goto 999
         else
            goto 50
         end if         
      end if

c
c Compute average/median redshift/age (also for each (sub)type) 
c
      call getzt(nsave,save,sntype,idxa,idx,idxb,ngood,nbad,
     $     ncut,zave,zmed,zsdev,agea,agem,agesdev,ntype,zavetype,
     $     zmedtype,zsdevtype,ageatype,agemtype,agesdevtype,ztemp,
     $     ztemperr,agetemp)

c    
c If nothing found, prompt user for some input
c
      if(ncut.eq.0) then
         write(6,'(2(a,f5.3),a,f5.1)') ' No template meets |z-',zuser,
     $        '| < ',zfilter,' and rlap >= ',rlapmin
         if(inter.eq.0) then
            write(6,*) 'Interactive mode off. No output written.'
            stop
         end if
 60      continue
         write(6,'(2a,$)')
     $        ' Enter a new (1) redshift; (2) zfilter; (3) rlapmin; ',
     $        '(4) lapmin; or (q)uit: '
         read(5,'(a)',err=60) reply
         if(reply.eq.'1') then
 61         continue
            write(6,'(a,f6.3,a,$)') 
     $           ' Enter a new redshift (current value ',zuser,'): '
            read(5,*,err=61) zuser
            goto 30
         else if(reply.eq.'2') then
 62         continue
            write(6,'(a,f6.3,a,$)') 
     $           ' Enter a new zfilter (current value ',zfilter,'): '
            read(5,*,err=62) zfilter
            goto 30
         else if(reply.eq.'3') then
 63         continue
            write(6,'(a,f6.1,a,$)') 
     $           ' Enter a new rlapmin (current value ',rlapmin,'): '
            read(5,*,err=63) rlapmin
            goto 30
         else if(reply.eq.'4') then
 64         continue
            write(6,'(a,f6.2,a,$)') 
     $           ' Enter a new lapmin (current value ',lapmin,'): '            
            read(5,*,err=64) lapmin
            goto 30
         else if(reply.eq.'q') then
            goto 999
         else
            goto 60
         end if
      end if






************************************************************************
**                                                                    **
**      Type statistics -- Absolute and relative fractions etc.       **
**                                                                    **
************************************************************************

* Determine best template(s) (based on rlap value)
      call getbtemp(ngood,idx,ntype,idxtype,save,imbtemp,imbtempr)
    
* Determine best (sub)type based on n best templates
      call getnbtemp(ngood,idx,sntype,nbt,nbst,nobst)

* Rank (sub)type based to absolute/relative fraction
      call rankfrac(ngood,ntype,fractype,rfractype,
     $     idxt,idxnt,idxnst,idxst,idxstr,imbt,imbst,imbstr)

* Rank (sub)type based on slope of absolute fraction with rlap [TEST]
      call rankslope(ngood,save,idx,sntype,rlapmax,rlaparr,fracrlap,
     $     fracrlaperr,slope,imbts,imbsts,idxnts,idxnsts,idxts,idxsts)






************************************************************************
**                                                                    **
**                      Write output to screen                        **
**                                                                    **
************************************************************************

c
c Best overall type and subtype
c
      bttemp = 0
      bsttemp = 0
      btfrac = 0
      btslope = 0
      btstfrac = 0
      bstfrac = 0
      btstslope = 0
      bstslope = 0
      okt = 0
      okst = 0

c
c Table header
c
      if(verbose.gt.0) then
         write(6,*) ' '
         write(6,*) ' '
         write(6,*) '                               SNID RESULTS'
         write(6,*) ' '
         write(6,1000) ' No.','Temp/Misc','Type','Subtype','Redshift',
     $        'Age'
         write(6,1001)
 1000    format(a4,2x,a9,5x,a4,8x,a7,6x,a8,7x,a5)
 1001    format('=====================================================',
     $        '=====================')
 1002    format('-----------------------------------------------------',
     $        '---------------------')
      end if
c
c Best-match templates (based on rlap value)
c
      if(verbose.gt.0) then
         write(6,'(a)') ' Best-match template(s):'
         do i=1, imbtemp
            write(6,1003) i,sid(idx(i)),
     $           typename(sntype(1,idx(i)),1),
     $           typename(sntype(1,idx(i)),sntype(2,idx(i))),
     $           save(3,idx(i)),save(9,idx(i)),
     $           save(4,idx(i))
         end do
 1003    format(i4,2x,a12,2x,2(a10,2x),f6.3,' (',f5.3,')',2x,f6.1)
      end if

      call notetemp(imbtemp,nbt,nbst,sntype,idx,bttemp,bsttemp,okt,okst,
     $     save,ngood,nbad,idxb)

c
c Best type (based on fraction of templates >= rlapmin)
c
      if(verbose.gt.0) then
         write(6,1002)
         write(6,'(a)') ' Best type(s):'
         write(6,'(a)') ' [fraction]'
         do it=1, imbt
            write(6,1004) it, fractype(idxt(it),1)*100,
     $           typename(idxt(it),1),
     $           zmedtype(idxt(it),1),zsdevtype(idxt(it),1),
     $           agemtype(idxt(it),1),agesdevtype(idxt(it),1)
         end do
 1004    format(i4,2x,f5.1,'%',8x,a10,2x,'---',9x,f6.3,' (',f5.3,')',2x,
     $        f6.1,' (',f5.1,')')
      end if

      call notetype(imbt,idxt,bttemp,btfrac,okt)

c
c Best type (based on slope of fraction vs. rlap >= rlapmin)
c
      if(verbose.gt.0) then
         write(6,'(a)') ' [slope]'
         do it=1, imbts
            write(6,10041) it, slope(idxts(it),1),
     $           typename(idxts(it),1),
     $           zmedtype(idxts(it),1),zsdevtype(idxts(it),1),
     $           agemtype(idxts(it),1),agesdevtype(idxts(it),1)
         end do
10041    format(i4,2x,f6.3,8x,a10,2x,'---',9x,f6.3,' (',f5.3,')',2x,
     $        f6.1,' (',f5.1,')')
      end if

      call notetypeslope(imbts,idxts,bttemp,btfrac,btslope,okt)

c
c Best subtype (based on fraction of templates >= rlapmin)
c
      if(verbose.gt.0) then
         write(6,1002)
         write(6,'(a)') ' Best subtype(s):'
         write(6,'(a)') ' [fraction]'
         do ist=1, imbst
            write(6,1005) ist, 
     $           fractype(idxnt(idxst(ist)),idxnst(idxst(ist)))*100,
     $           typename(idxnt(idxst(ist)),1),
     $           typename(idxnt(idxst(ist)),idxnst(idxst(ist))),
     $           zmedtype(idxnt(idxst(ist)),idxnst(idxst(ist))),
     $           zsdevtype(idxnt(idxst(ist)),idxnst(idxst(ist))),
     $           agemtype(idxnt(idxst(ist)),idxnst(idxst(ist))),
     $           agesdevtype(idxnt(idxst(ist)),idxnst(idxst(ist)))
         end do
 1005    format(i4,2x,f5.1,'%',8x,2(a10,2x),f6.3,' (',f5.3,')',2x,
     $        f6.1,' (',f5.1,')')
      end if
      
      call notesubt(imbst,idxt,idxnt,idxst,idxnst,btslope,btfrac,bttemp,
     $     bsttemp,btstfrac,bstfrac,okt,okst)

c
c Best subtype (based on slope of fraction vs. rlap >= rlapmin)
c
      if(verbose.gt.0) then
         write(6,'(a)') ' [slope]'
         do ist=1, imbsts
            write(6,10051) ist, 
     $           slope(idxnts(idxsts(ist)),idxnsts(idxsts(ist))),
     $           typename(idxnts(idxsts(ist)),1),
     $           typename(idxnts(idxsts(ist)),idxnsts(idxsts(ist))),
     $           zmedtype(idxnts(idxsts(ist)),idxnsts(idxsts(ist))),
     $           zsdevtype(idxnts(idxsts(ist)),idxnsts(idxsts(ist))),
     $           agemtype(idxnts(idxsts(ist)),idxnsts(idxsts(ist))),
     $           agesdevtype(idxnts(idxsts(ist)),idxnsts(idxsts(ist)))
         end do
10051    format(i4,2x,f6.3,8x,2(a10,2x),f6.3,' (',f5.3,')',2x,
     $        f6.1,' (',f5.1,')')
      end if

      call notesubtslope(imbsts,idxts,idxnts,idxsts,idxnsts,
     $     btslope,btfrac,bttemp,bsttemp,btstfrac,bstfrac,btstslope,
     $     bstslope,okt,okst)

c
c Favoured (sub)type
c
      if(verbose.gt.0) then      
         write(6,1001)
         write(6,*) ' '
         if(okt+okst.eq.8) then
            write(6,'(2a,2x,a)') ' OK! Favoured type and subtype are: ',
     $           typename(bttemp,1)(:lnb(typename(bttemp,1))),
     $           typename(bttemp,bsttemp)
     $           (:lnb(typename(bttemp,bsttemp)))
         else if(okt.eq.5) then
            write(6,'(3a)') ' NOTE: Favoured type is: ',
     $           typename(bttemp,1)(:lnb(typename(bttemp,1))),
     $           ', but there is no favoured subtype'
         else
            write(6,*) 'WARNING! No favoured type or subtype'
         end if
         if(nuse+nusetype+navoid+navoidtype.gt.0) then
            write(6,*) 'WARNING! You have used/avoided specific'//
     $           ' templates/types.'
            write(6,*) 'Re-run SNID with the full template suite'//
     $           ' for more reliable type statistics.'
         end if
      end if

c
c Get good/bad resdhift/age
c
      ngoodt = 0
      do k=1, ngood
         if(nint(save(5,idx(k))).eq.0) then
            ngoodt = ngoodt + 1
            zgood(ngoodt) = save(3,idx(k))
            tgood(ngoodt) = save(4,idx(k))
         end if
      end do

      nbadt = 0
      do k=1, nbad
         if(nint(save(5,idxb(k))).eq.0) then
            nbadt = nbadt + 1
            zbad(nbadt) = save(3,idxb(k))
            tbad(nbadt) = save(4,idxb(k))
         end if
      end do

c
c Pause for reflection
c      
      if(inter.gt.0) then
         write(6,*) ' '
 70      continue
         write(6,'(2a,$)') ' Hit <CR> to view template listings,',
     $        ' or (q)uit listings and move on: '
         read(5,'(a)',err=70) reply
         if(reply.eq.'q') then
            goto 300
         else if(lnb(reply).gt.0) then
            goto 70
         end if
      end if

c
c Template listings
c
      if(verbose.eq.0) goto 300
      write(6,*) ' '
      write(6,*) ' '
      write(6,*) '                         TEMPLATE LISTINGS'
      write(6,*) ' '
      write(6,2000)
 2000 format('   No.  Name          Type        lap    rlap    z',
     $     '    zerr    Age   NOTES')
      do i=1, ngood+nbad
         if(mod(i,MAXLIST).eq.0.and.inter.gt.0) then
            write(6,*) ' '
 90         continue
            write(6,'(2a,$)') ' Hit <CR> to view next set of templates',
     $           ' or (q)uit listings: '
            read(5,'(a)',err=90) reply
            if(reply.eq.'q') then
               goto 300
            else if(lnb(reply).gt.0) then
               goto 90
            end if
            write(6,*) ' '
            write(6,2000)
         end if
         rlap = save(1,idxa(i))*save(2,idxa(i))
         if(rlap.gt.MAXR) rlap = MAXR
c Good list
         if(abs(save(3,idxa(i))-zuser).lt.zfilter) then
            write(6,2001) i,sid(idxa(i)),
     $           typename(sntype(1,idxa(i)),sntype(2,idxa(i))),
     $           save(2,idxa(i)),rlap,
     $           save(3,idxa(i)),save(9,idxa(i)),save(4,idxa(i))            
 2001       format(2x,i4,2x,a12,2x,a10,1x,f6.3,f6.1,2(f7.3),f7.1,$)
* Age from explosion?
            iexp = 0
            open(unit=1,file=tempdir(:lnb(tempdir))//'texplist',
     $              status='old')
            do j=1, MAXTEMP
               read(1,'(a)',end=220) line
               if(index(line(nlb(line):nlb(line)+1),'#').eq.0) then
                  call parser(line,MAXTOK,ntok,ltok,tok)
                  if(index(sid(idxa(i)),tok(1)(:lnb(tok(1))))
     $                 .gt.0) then
                     read(tok(2),*) tmax
                     read(tok(3),*) texp
                     write(6,'(2x,a,f6.1,$)') 'Age from explosion =',
     $                    save(4,idxa(i))+tmax-texp
                     iexp = 1
                     goto 220
                  end if
               end if
            end do
 220        close(1)
* Age from first spectrum?
            if(nint(save(5,idxa(i))).eq.1) then
               if(iexp.eq.0) then
                  write(6,'(2x,a,$)') 'Age from 1st spectrum'
               else if(iexp.gt.0) then
                  write(6,'(a,$)') '; Age from 1st spectrum'
               end if
               open(unit=1,file=tempdir(:lnb(tempdir))//'tfirstlist',
     $              status='old')
               do j=1, MAXTEMP
                  read(1,'(a)',end=221) line
                  if(index(line(nlb(line):nlb(line)+1),'#').eq.0) then
                     call parser(line,MAXTOK,ntok,ltok,tok)
                     if(index(sid(idxa(i)),tok(1)(:lnb(tok(1))))
     $                    .gt.0) then
                        write(6,'(a,$)') 
     $                       ' ('//tok(2)(:lnb(tok(2)))//')'
                        goto 221
                     end if
                  end if
               end do
 221           close(1)
               write(6,*) ' '
            else if(nint(save(5,idxa(i))).eq.2) then
               if(index(typename(sntype(1,idxa(i)),sntype(2,idxa(i))),
     $              'M-star').gt.0) then
                  write(6,'(2x,a)') 'Age is MK type'
               else
                  write(6,'(2x,a)') 'Age N/A'
               end if
            else
               write(6,*) ' '
            end if
c Bad list
         else
            write(6,2002) 'X',i,sid(idxa(i)),
     $           typename(sntype(1,idxa(i)),sntype(2,idxa(i))),
     $           save(2,idxa(i)),rlap,
     $           save(3,idxa(i)),save(9,idxa(i)),save(4,idxa(i)),
     $           'WARNING! |z-',zuser,'| >',zfilter
 2002       format(a1,1x,i4,2x,a12,2x,a10,1x,f6.3,f6.1,2(f7.3),f7.1,2x,
     $           a12,f6.3,a4,f6.3,$)
* Age from explosion?
            open(unit=1,file=tempdir(:lnb(tempdir))//'texplist',
     $              status='old')
            do j=1, MAXTEMP
               read(1,'(a)',end=222) line
               if(index(line(nlb(line):nlb(line)+1),'#').eq.0) then
                  call parser(line,MAXTOK,ntok,ltok,tok)
                  if(index(sid(idxa(i)),tok(1)(:lnb(tok(1))))
     $                 .gt.0) then
                     read(tok(2),*) tmax
                     read(tok(3),*) texp
                     write(6,'(a,f6.1,$)') '; Age from explosion =',
     $                    save(4,idxa(i))+tmax-texp
                     goto 222
                  end if
               end if
            end do
 222        close(1)
* Age from first spectrum?
            if(nint(save(5,idxa(i))).eq.1) then
               write(6,'(a,$)') '; Age from 1st spectrum'
               open(unit=1,file=tempdir(:lnb(tempdir))//'tfirstlist',
     $              status='old')
               do j=1, MAXTEMP
                  read(1,'(a)',end=223) line
                  if(index(line(nlb(line):nlb(line)+1),'#').eq.0) then
                     call parser(line,MAXTOK,ntok,ltok,tok)
                     if(index(sid(idxa(i)),tok(1)(:lnb(tok(1))))
     $                    .gt.0) then
                        write(6,'(a,$)') 
     $                       ' ('//tok(2)(:lnb(tok(2)))//')'
                        goto 223
                     end if
                  end if
               end do
 223           close(1)
               write(6,*) ' '
            else if(nint(save(5,idxa(i))).eq.2) then
               if(index(typename(sntype(1,idxa(i)),sntype(2,idxa(i))),
     $              'M-star').gt.0) then
                  write(6,'(a)') '; Age is MK type'
               else
                  write(6,'(a)') '; Age N/A'
               end if
            else
               write(6,*) ' '
            end if
         end if
      end do






************************************************************************
**                                                                    **
**                            File output                             **
**                                                                    **
************************************************************************
 300  continue

c
c Generic output file
c
      if(fout.ge.0) call wfout(zinit,zmed,zsdev,agem,agesdev,
     $     ntype,fractype,slope,zmedtype,zsdevtype,agemtype,agesdevtype,
     $     nsave,ngood,nbad,save,idxa,sid,sntype)

c
c Other output
c
      nout = 0
      nfluxout = 0
      nflatout = 0
      nxcorout = 0

      if(fout.gt.0.or.fluxout.gt.0.or.flatout.gt.0.or.xcorout.gt.0) then

         nout = max(fout,fluxout)
         nout = max(nout,flatout)
         nout = max(nout,xcorout)
         if(nout.gt.nsave) then
            write(6,*) 'WARNING! Requesting too many output files'
            write(6,*) 'Setting fout = ',nsave
            nout = nsave
         end if
         
         nfluxout = max(fout,fluxout)
         nflatout = max(fout,flatout)
         nxcorout = max(fout,xcorout)
         if(nfluxout.gt.nsave) nfluxout = nsave
         if(nflatout.gt.nsave) nflatout = nsave
         if(nxcorout.gt.nsave) nxcorout = nsave

         do i=1, nout
            j = nint(save(8,idxa(i)))
            call getpdat(pdat,fnorm,dft,temp(1,j),tft,cfn,cft,
     $           ctr,hgt,drms,trms,arms,srms,r,redshift,wid,
     $           save(3,idxa(i)),nknod,xknod,yknod,nknot(j),
     $           xknot(1,j),yknot(1,j))
            call getdlabels(i,sid(idxa(i)),
     $           typename(sntype(1,idxa(i)),sntype(2,idxa(i))),
     $           save(4,idxa(i)),save(3,idxa(i)),save(9,idxa(i)),
     $           save(1,idxa(i)))

c fluxed spectra
            if(i.le.nfluxout) then
               call wfluxout(i,save(3,idxa(i)),pdat(1,2),pdat(1,1),
     $              pdat(1,15),pdat(1,8))
            end if

c flattened spectra
            if(i.le.nflatout) then
               call wflatout(i,save(3,idxa(i)),pdat(1,2),pdat(1,8))
            end if

c correlation function
            if(i.le.nxcorout) then
               call wxcorout(itn,npeak,pdat(1,19),pdat(1,11),pdat(1,12))
            end if

         end do

         if(verbose.gt.0) then
            if(nfluxout.gt.0) 
     $           write(6,*) 'Wrote top ',nfluxout,
     $           ' fluxed template spectra to file *_snidflux.dat'
            if(nflatout.gt.0) 
     $           write(6,*) 'Wrote top ',nflatout,
     $           ' flattened template spectra to file *_snidflat.dat'
            if(nxcorout.gt.0) 
     $           write(6,*) 'Wrote top ',nxcorout,
     $           ' correlation template functions to file',
     $           ' *_snidxcor.dat'
         end if

      end if






************************************************************************
**                                                                    **
**                       Interactive plotting                         **
**                                                                    **
************************************************************************
 100  continue
      if(iplot.eq.0) goto 999
      if(verbose.gt.0) then
         write(6,*) ' '
         write(6,'(a)') ' On to interactive plotting...'
      end if

c
c Initialize variables
c
      call initplotvar
      call setplotvar(0)

c
c Open graphics output
c
 110  continue
      call rpgbegok('/XWINDOW')

c
c Setup plot/button regions
c
      call initbutt
      call setbutt
      if(ips2.gt.0) goto 140

c
c Plot best-match template by default
c
      itn = 1
      nb = 38
      goto 121

c
c Main button handler
c
 120  continue
      call rpgband(0,0,0.,0.,xc,yc,ch)
      call ifbutton(xc,yc,nb)
 121  continue

c
c Determine template number or change what is being plotted
c
      if(nb.eq.0.or.nb.eq.7.or.nb.eq.14.or.nb.eq.37.or.nb.eq.18.or.
     $     nb.eq.20.or.nb.eq.24) then
         goto 120
      else if(nb.eq.12) then
         call button(12,'QUIT',5)
         goto 200
      else if(nb.eq.16) then
         call button(16,'No.',5)
 130     continue
         write(6,'(a,i4,a,$)') ' Enter template No. [1 -',
     $        ngood+nbad,']: '
         read(5,*,err=130) itn
         if(itn.eq.0.or.itn.gt.ngood+nbad) goto 130
         call button(16,'No.',0)
      else if(nb.eq.25) then
         itn = 1
      else if(nb.eq.26) then
         if(itn.gt.1) itn = itn - 1
      else if(nb.eq.27) then
         if(itn.lt.ngood+nbad) itn = itn + 1
      else if(nb.eq.28) then
         itn = ngood + nbad
      else if(nb.eq.36) then
         ips = 5
         call button(36,'PS',ips)
      else if(nb.eq.48) then
         iascii = 5
         call button(48,'ASCII',iascii)                  
      else
         call setplotvar(nb)
         call setbutt
      end if

c
c Plot spectra, xcor, peaks, or fractions
c
 140  continue
      ips2 = 0

* get static labels
      call getlabels

      if(iflux.gt.0.or.iflat.gt.0.or.ixcor.gt.0) then

* Re-run cross-correlation on template itn and get pdat array for plotting
         j = nint(save(8,idxa(itn)))
         call getpdat(pdat,fnorm,dft,temp(1,j),tft,cfn,cft,
     $        ctr,hgt,drms,trms,arms,srms,r,redshift,wid,
     $        save(3,idxa(itn)),nknod,xknod,yknod,nknot(j),
     $        xknot(1,j),yknot(1,j))

* get dynamic labels
         call getdlabels(itn,sid(idxa(itn)),
     $        typename(sntype(1,idxa(itn)),sntype(2,idxa(itn))),
     $        save(4,idxa(itn)),save(3,idxa(itn)),save(9,idxa(itn)),
     $        save(1,idxa(itn)))


c Plot fluxed spectra
         if(iflux.gt.0) then
            call plotflux(itn,save(3,idxa(itn)),pdat(1,2),pdat(1,1),
     $           pdat(1,15),pdat(1,8))
            if(ips2.gt.0) goto 110

c Plot flattened spectra
         else if(iflat.gt.0) then
            call plotflat(itn,save(3,idxa(itn)),pdat(1,2),pdat(1,8))
            if(ips2.gt.0) goto 110            

c Plot correlation functions
         else if(ixcor.gt.0) then
            call plotxcor(itn,npeak,pdat(1,19),pdat(1,11),pdat(1,12))
            if(ips2.gt.0) goto 110

         end if

c Plot correlation peaks
      else if(ipeaks.gt.0.and.forcez.lt.-1) then
         call plotpeaks(npeak,zpeak,rpeak,lpeak)
         if(ips2.gt.0) goto 110

c Plot redshift/age scatter plot
      else if(izt.gt.0) then
         call plotzt(ngoodt,nbadt,zgood,tgood,zbad,tbad,zgood(1),
     $        tgood(1),zmed,agem)
         if(ips2.gt.0) goto 110

c Plot (sub)type fractions
      else if(iia.gt.0) then      ! Ia
         it = 1
         call plotfrac(it,rlapmax+1,rlaparr,fracrlap,fracrlaperr)
         if(ips2.gt.0) goto 110
      else if(iib.gt.0) then      ! Ib
         it = 2
         call plotfrac(it,rlapmax+1,rlaparr,fracrlap,fracrlaperr)
         if(ips2.gt.0) goto 110
      else if(iic.gt.0) then      ! Ic
         it = 3 
         call plotfrac(it,rlapmax+1,rlaparr,fracrlap,fracrlaperr)
         if(ips2.gt.0) goto 110
      else if(iii.gt.0) then      ! II
         it = 4
         call plotfrac(it,rlapmax+1,rlaparr,fracrlap,fracrlaperr)
         if(ips2.gt.0) goto 110
      else if(inotsn.gt.0) then   ! NotSN
         it = 5
         call plotfrac(it,rlapmax+1,rlaparr,fracrlap,fracrlaperr)
         if(ips2.gt.0) goto 110
      else if(iall.gt.0) then     ! All
         call plotfracall(rlapmax+1,rlaparr,fracrlap,fracrlaperr)
         if(ips2.gt.0) goto 110

      end if
      goto 120      

 200  continue
      call pgend






************************************************************************
**                                                                    **
**                        End of program snid                         **
**                                                                    **
************************************************************************

 999  if(verbose.gt.0) then
         write(6,*) ' '
         write(6,'(a)') ' Thank you for using SNID! Goodbye.'
      end if

* clear exceptions (for Sun compilers)
c Uncomment the line below (call to ieee_flags) if you see:
c "Note: IEEE floating-point exception flags raised:"
C      iflag = ieee_flags('clear','exception','all',out)

      stop
      end
