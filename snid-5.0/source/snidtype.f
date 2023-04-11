************************************************************************
**                                                                    **
**                Type determination routines for SNID                **
**                                                                    **
**                                                                    **
**                Copyright (C) 2007 Stephane Blondin                 **
**                                                                    **
************************************************************************

c Included in this file:
c
c subroutine getbtemp
c subroutine getnbtemp
c subroutine rankfrac
c subroutine rankslope
c subroutine notetemp
c subroutine notetype
c subroutine notetypeslope
c subroutine notesubt
c subroutine notesubtslope



************************************************************************
* subroutine getbtemp -- Determine best template(s) based on rlap value
************************************************************************

      subroutine getbtemp(ngood,idx,ntype,idxtype,save,imbtemp,imbtempr)

      implicit none

      include 'snid.inc'

      integer i
      integer itemp                   ! loop index for templates
      integer ngood                   ! number of "good" (rlap>rlapmin AND |z-zmed|<zfilter) correlations
      integer idx(MAXTEMP)            ! sorted index for good templates (based on rlap value and abs(z-zuser) < zfilter)
      integer ntype(NT,NST)           ! total number of SN in each (sub)type
      integer idxtype(MAXTEMP,NT,NST) ! ranked SNe in each (sub)type
      real save(9,MAXTEMP)            ! array of saved values for statistics and later analysis
      integer imbtemp                 ! number of multiple best templates (based on rlap value)
      integer imbtempr(NT,NST)        ! number of multiple best templates in each (sub)type (based on rlap value)

* Initialize variables
      imbtemp = 0
      do it=1, NT
         do ist=1, NST
            imbtempr(it,ist) = 0
         end do
      end do
      if(ngood.eq.0) return
      
* Check for multiple best templates
      imbtemp = 1
      if(ngood.le.1) goto 10
      do itemp=2, ngood
         if(abs(save(1,idx(itemp))*save(2,idx(itemp)) - 
     $        save(1,idx(1))*save(2,idx(1))).lt.EPSRLAP) then
            imbtemp = imbtemp+1
         else
            goto 10
         end if
      end do
 10   continue

* Same within each (sub)type
      do it=1, NT
         do ist=1, NST            
            if(idxtype(1,it,ist).ne.0) then
               imbtempr(it,ist) = 1
               if(ntype(it,ist).gt.1) then
                  do itemp=2, ntype(it,ist)
                     if(abs(save(1,idx(idxtype(itemp,it,ist))) 
     $                    *save(2,idx(idxtype(itemp,it,ist))) - 
     $                    save(1,idx(idxtype(1,it,ist)))
     $                    *save(2,idx(idxtype(1,it,ist))))
     $                    .lt.EPSRLAP)then
                        imbtempr(it,ist) = imbtempr(it,ist)+1
                     else
                        goto 20
                     end if
                  end do
               end if
            end if
         end do
      end do
 20   continue

      return
      end



************************************************************************
* subroutine getnbtemp -- Determine (sub)type based on n best templates
************************************************************************

      subroutine getnbtemp(ngood,idx,sntype,nbt,nbst,nobst)

      implicit none

      include 'snid.inc'

      integer itemp             ! loop index for templates
      integer ngood             ! number of "good" (rlap>rlapmin AND |z-zmed|<zfilter) correlations
      integer idx(MAXTEMP)      ! sorted index for good templates (based on rlap value and abs(z-zuser) < zfilter)
      integer sntype(2,MAXTEMP) ! (sub)type indeces
      integer nbt,nbst,nobst    ! number of best (sub)type based on n best templates

      nbt = 0
      nbst = 0
      nobst = 0
      if(ngood.eq.0) return

      nbt = 1
      nbst = 1
      if(ngood.le.1) return
      do itemp=2, ngood
         if(sntype(1,idx(itemp)).eq.sntype(1,idx(1))) then
            nbt=nbt+1
            if(sntype(2,idx(itemp)).eq.sntype(2,idx(1)) .and.
     $           nobst.eq.0) then
               nbst=nbst+1 
            else
               nobst = 1
            end if
         else
            goto 10
         end if
      end do
 10   continue      

      return
      end



************************************************************************
* subroutine rankfrac -- Rank (sub)type based on absolute/relative fraction
************************************************************************

      subroutine rankfrac(ngood,ntype,fractype,rfractype,
     $     idxt,idxnt,idxnst,idxst,idxstr,imbt,imbst,imbstr)      

      implicit none

      include 'snid.inc'

      integer ngood              ! number of "good" (rlap>rlapmin AND |z-zmed|<zfilter) correlations
      integer ntype(NT,NST)      ! total number of SN in each (sub)type
      real fractype(NT,NST)      ! absolute fraction of SN in each (sub)type
      real rfractype(NT,NST)     ! relative fraction of SN in each (sub)type
      integer idxt(NT)           ! sorted index for types (based on absolute fraction) 
      integer idxnt(NT*(NST-1))  ! index of type name corresponding to each subtype (for absolute fraction)
      integer idxnst(NT*(NST-1)) ! index of subtype name (for absolute fraction)
      integer idxst(NT*(NST-1))  ! sorted index for subtypes (based on absolute fraction) 
      integer idxstr(NT,NST-1)   ! sorted index for subtypes (based on relative fraction)
      integer imbt               ! number of multiple best types (based on absolute fraction) 
      integer imbst              ! number of multiple best subtypes (based on absolute fraction) 
      integer imbstr(NST)        ! number of multiple best subtypes (based on relative fraction)     

c Local variables
      real tmpt(NT)            ! temporary array for ranking types (based on absolute fraction)
      real tmpst(NT*(NST-1))   ! temporary array for ranking subtypes (based on absolute fraction)
      integer ntmpst           ! number of elements in tmpst (dummy)  
      real tmpstr(NST-1)       ! temporary array for ranking subtypes (based on relative fraction within a given type)
      integer tmpidxstr(NST-1) ! temporary array holding values for idxstr
      real hoser               ! dummy for median and rlap sorting

c Functions
      real amedidx
      
* Initialize variables
      do it=1, NT
         do ist=1, NST
            if(it.eq.1) imbstr(ist) = 0
            fractype(it,ist) = 0
            rfractype(it,ist) = 0
         end do
      end do  
      ntmpst = 0
      imbt = 0
      imbst = 0
      if(ngood.eq.0) return

* Absolute/relative (sub)type fractions
      do it=1, NT
         do ist=1, NST
            fractype(it,ist) = (ntype(it,ist)/real(ngood))               
            if (ntype(it,1).gt.0)
     $      rfractype(it,ist)= (ntype(it,ist)/real(ntype(it,1)))
         end do
      end do

* Rank (sub)type based on absolute fraction 
* (check for multiple best (sub)types)
      do it=1, NT
         tmpt(it) = -fractype(it,1)
         idxt(it) = it
         do ist=2, NST
            ntmpst = ntmpst + 1
            tmpst(ntmpst)   = -fractype(it,ist)
            idxnt(ntmpst)   = it
            idxnst(ntmpst)  = ist
            idxst(ntmpst)   = ntmpst
         end do
      end do

      hoser = amedidx(NT,tmpt,idxt)
      hoser = amedidx(ntmpst,tmpst,idxst)

      if(tmpt(1).lt.0) then
         imbt = 1
         do it=2, NT
            if(abs(tmpt(it)-tmpt(1)).lt.EPSFRAC) imbt = imbt + 1
         end do
      end if

      if(tmpst(1).lt.0) then
         imbst = 1
         do ist=2, ntmpst
            if(abs(tmpst(ist)-tmpst(1)).lt.EPSFRAC) imbst = imbst + 1
         end do
      end if

* Rank subtype based on relative fraction within a given type 
* (check for multiple best types)
      do ist=1, NT
         do it=2, NST
            tmpstr(it-1) = -rfractype(ist,it)
            tmpidxstr(it-1) = it
         end do

         hoser = amedidx(NST-1,tmpstr,tmpidxstr)

         do it=1, NST-1 
            idxstr(ist,it) = tmpidxstr(it)
         end do

        if(tmpstr(1).lt.0) then
           imbstr(ist) = 1
           do it=2, NST
              if(abs(tmpstr(it)-tmpstr(1)).lt.EPSFRAC) 
     $             imbstr(ist) = imbstr(ist) + 1
           end do
        end if
      end do

      return
      end



************************************************************************
* subroutine rankslope -- Rank (sub)type based on slope of absolute fraction 
*                         with rlap
************************************************************************

      subroutine rankslope(ngood,save,idx,sntype,rlapmax,rlaparr,
     $     fracrlap,fracrlaperr,slope,imbts,imbsts,idxnts,idxnsts,
     $     idxts,idxsts)

      implicit none

      include 'snid.inc'

      integer i,j
      integer ngood                     ! number of "good" (rlap>rlapmin AND |z-zmed|<zfilter) correlations
      real save(9,MAXTEMP)              ! array of saved values for statistics and later analysis
      integer idx(MAXTEMP)              ! sorted index for good templates (based on rlap value and abs(z-zuser) < zfilter)
      integer sntype(2,MAXTEMP)         ! (sub)type indeces
      integer rlapmax                   ! maximum integer rlap value
      real rlaparr(MAXRLAP)             ! rlap array (rlaparr[0:rlapmax])
      real fracrlap(MAXRLAP,NT,NST)     ! absolute fraction of templates with rlap >= rlap_i for each (sub)type
      real fracrlaperr(MAXRLAP,NT,NST)  ! poisson error in absolute fraction of templates with rlap >= rlap_i for each (sub)type
      real slope(NT,NST)                ! slope of fracrlap(rlaparr)
      integer imbts                     ! number of multiple best types (based on slope) 
      integer imbsts                    ! number of multiple best subtypes (based on slope) 
      integer idxnts(NT*(NST-1))        ! index of type name corresponding to each subtype (for slope)
      integer idxnsts(NT*(NST-1))       ! index of subtype name (for slope)
      integer idxts(NT)                 ! sorted index for types (based on slope) 
      integer idxsts(NT*(NST-1))        ! sorted index for subtypes (based on slope)
      integer ncutrlap(MAXRLAP)         ! number of templates with rlap >= rlap_i
      integer ntyperlap(MAXRLAP,NT,NST) ! number of templates with rlap >= rlap_i for each (sub)type

c Local variables
      integer ntmpst                          ! number of elements in tmpst (dummy)  
      integer iposf                           ! to ensure at least 2 positive fractions for linear fit
      integer nrlap                           ! number of elements of [rlapmin,rlapmax]
      real x(2,MAXRLAP),y(MAXRLAP),w(MAXRLAP) ! {x,y} and weights for linearfit
      real param(2)                           ! fit parameters for linearfit subroutine
      real tmpts(NT)                          ! temporary array for ranking types (based on slope)
      real tmpsts(NT*(NST-1))                 ! temporary array for ranking subtypes (based on slope)
      real hoser                              ! dummy for median and rlap sorting

c Functions
      real amedidx      

* Create rlap array [0,rlapmax] and initialize variables
      do i=1, MAXRLAP
         rlaparr(i) = i-1
         ncutrlap(i) = 0
         do it=1, NT
            do ist=1, NST
               ntyperlap(i,it,ist) = 0
               fracrlap(i,it,ist) = 0.
               fracrlaperr(i,it,ist) = 1.
               if(i.eq.1) slope(it,ist) = 0.
            end do
         end do
      end do
      ntmpst = 0
      imbts = 0
      imbsts = 0

      rlapmax = 0
      if(ngood.eq.0) return
      rlapmax = int(save(1,idx(1))*save(2,idx(1)))
      if(rlapmax.ge.MAXRLAP) rlapmax = MAXRLAP-1

* Compute number of templates with rlap > rlap_i
      do 10 i=1, rlapmax+1
         do 20 j=1, ngood
            if(save(1,idx(j))*save(2,idx(j)).ge.i-1) then
               ncutrlap(i) = ncutrlap(i) + 1
            else
               goto 10
            end if
            do it=1, NT
               if(sntype(1,idx(j)).eq.it) then
                  ntyperlap(i,it,1) = ntyperlap(i,it,1) + 1
                  do ist=2, NST
                     if(sntype(2,idx(j)).eq.ist) then
                        ntyperlap(i,it,ist) = ntyperlap(i,it,ist) + 1
                        goto 20
                     end if
                  end do
               end if
            end do
 20      continue
 10   continue
      
* Compute absolute fractions
      do i=1, rlapmax+1
         if(ncutrlap(i).gt.0) then
            do it=1, NT
               do ist=1, NST
                  fracrlap(i,it,ist) = real(ntyperlap(i,it,ist)) / 
     $                 real(ncutrlap(i))
                  if(ntyperlap(i,it,ist).gt.0) then
                     fracrlaperr(i,it,ist) = 
     $                    sqrt(real(ntyperlap(i,it,ist))) / 
     $                    real(ncutrlap(i))
                  else
                     fracrlaperr(i,it,ist) = 1.
                  end if
               end do
            end do
         end if
      end do

* Determine slope of fraction vs. rlap
      if(rlapmax.eq.0) return
      nrlap = rlapmax-int(rlapmin)+1
      do it=1, NT
         do ist=1, NST
            iposf = 0
            do i=1, nrlap
               x(1,i) = 1.
               x(2,i) = rlaparr(int(rlapmin)+i)
               y(i) = fracrlap(int(rlapmin)+i,it,ist)
               if(fracrlaperr(int(rlapmin)+i,it,ist).gt.EPSFRAC) then
                  w(i) = 1./fracrlaperr(int(rlapmin)+i,it,ist)
               else
                  w(i) = 1.
               end if
               if(y(i).gt.EPSFRAC) iposf = iposf + 1
            end do
            if(iposf.gt.1) then 
               call wlinearfit(nrlap,y,w,2,x,param)
C               call linearfit(nrlap,y,2,x,param)
               slope(it,ist) = param(2)
            end if
         end do
      end do

* Rank (sub)type based on slope
* (check for multiple best (sub)types)
      do it=1, NT
         tmpts(it) = -slope(it,1)
         idxts(it) = it
         do ist=2, NST
            ntmpst = ntmpst + 1
            tmpsts(ntmpst)   = -slope(it,ist)
            idxnts(ntmpst)   = it
            idxnsts(ntmpst)  = ist
            idxsts(ntmpst)   = ntmpst
         end do
      end do
      hoser = amedidx(NT,tmpts,idxts)
      hoser = amedidx(ntmpst,tmpsts,idxsts)

      if(abs(tmpts(1)).gt.EPSSLOPE) then
         imbts = 1
         do it=2, NT
            if(abs(tmpts(it)-tmpts(1)).lt.EPSSLOPE) 
     $           imbts = imbts + 1
         end do
      end if

      if(abs(tmpsts(1)).gt.EPSSLOPE) then
         imbsts = 1
         do ist = 2, ntmpst
            if(abs(tmpsts(ist)-tmpsts(1)).lt.EPSSLOPE) 
     $           imbsts = imbsts + 1
         end do
      end if
      
C      if(verbose.gt.0) then
C         write(6,*) 'Best type according to slope: ', 
C     $        (typename(idxts(i),1),i=1,imbts)
C         write(6,*) 'Best subtype according to slope: ',
C     $        (typename(idxnts(idxsts(i)),idxnsts(idxsts(i))),
C     $        i=1,imbsts)
C      end if

      return
      end



************************************************************************
* subroutine notetemp -- Notes on best-match templates
************************************************************************

      subroutine notetemp(imbtemp,nbt,nbst,sntype,idx,bttemp,bsttemp,
     $     okt,okst,save,ngood,nbad,idxb)

      implicit none

      include 'snid.inc'

      integer imbtemp           ! number of multiple best templates (based on rlap value)
      integer nbt,nbst          ! number of best (sub)type based on n best templates
      integer sntype(2,MAXTEMP) ! (sub)type indeces
      integer idx(MAXTEMP)      ! sorted index for good templates (based on rlap value and abs(z-zuser) < zfilter)
      integer bttemp,bsttemp    ! (sub)type of best-match templates
      integer samebt,samebst    ! number of best (sub)types
      integer okt,okst          ! secure (sub)type indicators
      real save(9,MAXTEMP)      ! array of saved values for statistics and later analysis
      integer ngood,nbad        ! number of good and bad correlations
      integer idxb(MAXTEMP)     ! sorted index for bad templates (based on rlap value)

      bttemp = 0
      bsttemp = 0
      if(imbtemp.eq.1) then
         okt = okt + 1
         okst = okst + 1
         bttemp = sntype(1,idx(1))
         bsttemp = sntype(2,idx(1))
         if(verbose.gt.0) then
            write(6,*) 'NOTE: the top ',nbt,' (',nbst,')',
     $           ' templates have same type (subtype)'
         end if
      else if(imbtemp.gt.1) then
         samebt = 0
         do it=1, imbtemp
            if(sntype(1,idx(it)).eq.sntype(1,idx(1))) 
     $           samebt = samebt + 1
         end do
         if(samebt.eq.imbtemp) then
            okt = okt + 1
            samebst = 0
            do ist=1, imbtemp
               if(sntype(2,idx(ist)).eq.sntype(2,idx(1))) 
     $              samebst = samebst + 1
            end do
            if(samebst.eq.imbtemp) then
               okst = okst + 1
               bttemp = sntype(1,idx(1))
               bsttemp = sntype(2,idx(1))
               if(verbose.gt.0) then
                  write(6,*) 'NOTE: the top ',nbt,' (',nbst,')',
     $                 ' templates have same type/subtype' 
               end if
            else
               bttemp = sntype(1,idx(1))
               if(verbose.gt.0) then
                  write(6,*) 'WARNING! Best-match templates have'//
     $                 ' same type but different subtypes'
                  write(6,*) 'NOTE: the top ',nbt,
     $                 ' templates have same type'
               end if
            end if
         else if(verbose.gt.0) then
            write(6,*) 'WARNING! Best-match templates have'//
     $              ' different types'
         end if
      end if

* Beware of type if "bad" template has highest rlap
      if(nbad.gt.0.and.ngood.gt.0) then
         if(save(1,idxb(1))*save(2,idxb(1)).ge.
     $        save(1,idx(1))*save(2,idx(1))) then
            write(6,'(2a,f6.3)') ' WARNING! Template with highest rlap',
     $           ' has z =',save(3,idxb(1))
            if(sntype(1,idxb(1)).ne.sntype(1,idx(1))) then
               write(6,*) 'WARNING! Template with highest rlap has',
     $              ' different type: ',typename(sntype(1,idxb(1)),1)
               okt = 0
               okst = 0
            else if(sntype(1,idxb(1)).eq.sntype(1,idx(1)) .and.
     $              sntype(2,idxb(1)).ne.sntype(2,idx(1))) then
               write(6,*) 'WARNING! Template with highest rlap has',
     $              ' different subtype: ',
     $              typename(sntype(1,idxb(1)),sntype(2,idxb(1)))
               okst = 0
            end if
         end if
      end if

      return
      end



************************************************************************
* subroutine notetype -- Notes on best types
************************************************************************

      subroutine notetype(imbt,idxt,bttemp,btfrac,okt)

      implicit none

      include 'snid.inc'

      integer imbt     ! number of multiple best types (based on absolute fraction) 
      integer idxt(NT) ! sorted index for types (based on absolute fraction) 
      integer bttemp   ! type of best-match templates
      integer btfrac   ! type with highest fraction
      integer okt      ! secure type indicator      

      if(imbt.eq.1) then 
         btfrac = idxt(1)
         if(bttemp.ne.0.and.btfrac.eq.bttemp) okt = okt + 1
      else if(imbt.gt.1.and.verbose.gt.0) then
         write(6,*) 'WARNING! Multiple best types'
      end if      

      return
      end



************************************************************************
* subroutine notetypeslope -- Notes on best types [from slope]
************************************************************************

      subroutine notetypeslope(imbts,idxts,bttemp,btfrac,btslope,okt)

      implicit none

      include 'snid.inc'

      integer imbts     ! number of multiple best types (based on slope) 
      integer idxts(NT) ! sorted index for types (based on slope) 
      integer bttemp    ! type of best-match templates
      integer btfrac    ! type with highest fraction
      integer btslope   ! type with highest slope
      integer okt       ! secure type indicator      

      if(imbts.eq.0) then
         okt = okt + 1
         if(verbose.gt.0) then
            write(6,*) 'NOTE: all null slopes'
         end if
      else if(imbts.eq.1) then 
         btslope = idxts(1)
         if(bttemp.ne.0.and.btfrac.eq.bttemp.and.btslope.eq.btfrac) 
     $        okt = okt + 1
      else if(imbts.gt.1.and.verbose.gt.0) then
         write(6,*) 'WARNING! Multiple best types'
      end if      

      return
      end



************************************************************************
* subroutine notesubt -- Notes on best subtypes
************************************************************************

      subroutine notesubt(imbst,idxt,idxnt,idxst,idxnst,btslope,btfrac,
     $     bttemp,bsttemp,btstfrac,bstfrac,okt,okst)

      implicit none

      include 'snid.inc'

      integer i
      integer imbst                   ! number of multiple best subtypes (based on absolute fraction) 
      integer idxt(NT)                ! sorted index for types (based on absolute fraction) 
      integer idxnt(NT*(NST-1))       ! index of type name corresponding to each subtype (for absolute fraction)
      integer idxst(NT*(NST-1))       ! sorted index for subtypes (based on absolute fraction) 
      integer idxnst(NT*(NST-1))      ! index of subtype name (for absolute fraction)
      integer btslope                 ! type with highest slope
      integer btfrac,bstfrac,btstfrac ! (sub)type with highest fraction
      integer bttemp,bsttemp          ! (sub)type of best-match templates
      integer okt,okst                ! secure (sub)type indicators
      integer samebt                  ! number of best types

      if(imbst.eq.1) then

         btstfrac = idxnt(idxst(1))
         bstfrac = idxnst(idxst(1))

         if((btslope.eq.0.or.(btslope.gt.0.and.btstfrac.eq.btslope))
     $        .and.(btfrac.ne.0.and.btstfrac.eq.btfrac)
     $        .and.(bttemp.ne.0.and.btstfrac.eq.bttemp)) okt = okt + 1

         if((bttemp.ne.0.and.btstfrac.eq.bttemp)
     $        .and.(bsttemp.ne.0.and.bstfrac.eq.bsttemp)) okst = okst+1

      else if(imbst.gt.1) then

         samebt = 0
         do i=1, imbst
            if(idxnt(idxst(i)).eq.idxt(1)) samebt = samebt + 1
         end do

         if(samebt.eq.imbst) then
            btstfrac = idxnt(idxst(1))
            if((btslope.eq.0.or.(btslope.gt.0.and.btstfrac.eq.btslope))
     $           .and.(btfrac.ne.0.and.btstfrac.eq.btfrac)
     $           .and.(bttemp.ne.0.and.btstfrac.eq.bttemp)) okt = okt+1
            if(verbose.gt.0) then
               write(6,*) 'WARNING! Multiple best subtypes'//
     $              ' within same type'
            end if
         else if(verbose.gt.0) then
            write(6,*) 'WARNING! Multiple best subtypes of'//
     $           ' different types'
         end if

      end if

      return
      end



************************************************************************
* subroutine notesubtslope -- Notes on best subtypes [from slope]
************************************************************************

      subroutine notesubtslope(imbsts,idxts,idxnts,idxsts,idxnsts,
     $     btslope,btfrac,bttemp,bsttemp,btstfrac,bstfrac,btstslope,
     $     bstslope,okt,okst)

      implicit none

      include 'snid.inc'

      integer i
      integer imbsts                     ! number of multiple best subtypes (based on slope) 
      integer idxts(NT)                  ! sorted index for types (based on slope) 
      integer idxnts(NT*(NST-1))         ! index of type name corresponding to each subtype (for slope)
      integer idxsts(NT*(NST-1))         ! sorted index for subtypes (based on slope) 
      integer idxnsts(NT*(NST-1))        ! index of subtype name (for slope)
      integer btslope,bstslope,btstslope ! (sub)type with highest slope
      integer btfrac,bstfrac,btstfrac    ! (sub)type with highest fraction
      integer bttemp,bsttemp             ! (sub)type of best-match templates
      integer okt,okst                   ! secure (sub)type indicators
      integer samebt                     ! number of best types

      if(imbsts.eq.0) then
         okt = okt + 1
         okst = okst + 1
         if(verbose.gt.0) then
            write(6,*) 'NOTE: all null slopes'
         end if

      else if(imbsts.eq.1) then
         btstslope = idxnts(idxsts(1))
         bstslope  = idxnsts(idxsts(1))

         if((btslope.eq.0.or.(btslope.gt.0.and.btstslope.eq.btslope))
     $        .and.(btfrac.ne.0.and.btstslope.eq.btfrac)
     $        .and.(bttemp.ne.0.and.btstslope.eq.bttemp)
     $        .and.(btstfrac.ne.0.and.btstslope.eq.btstfrac)) okt=okt+1

         if((bttemp.ne.0.and.btstslope.eq.bttemp)
     $        .and.(bsttemp.ne.0.and.bstslope.eq.bsttemp)
     $        .and.(btstfrac.ne.0.and.btstslope.eq.btstfrac)
     $        .and.(bstfrac.ne.0.and.bstslope.eq.bstfrac)) okst = okst+1

      else if(imbsts.gt.1) then
         samebt = 0
         do i=1, imbsts
            if(idxnts(idxsts(i)).eq.idxts(1)) samebt = samebt + 1
         end do

         if(samebt.eq.imbsts) then
            btstslope = idxnts(idxsts(1))
            if((btslope.eq.0.or.(btslope.gt.0.and.btstslope.eq.btslope))
     $        .and.(btfrac.ne.0.and.btstslope.eq.btfrac)
     $        .and.(bttemp.ne.0.and.btstslope.eq.bttemp)
     $        .and.(btstfrac.ne.0.and.btstslope.eq.btstfrac)) okt=okt+1
            if(verbose.gt.0) then
               write(6,*) 'WARNING! Multiple best subtypes'//
     $              ' within same type'
            end if
         else if(verbose.gt.0) then
            write(6,*) 'WARNING! Multiple best subtypes of'//
     $           ' different types'
         end if

      end if

      return
      end
