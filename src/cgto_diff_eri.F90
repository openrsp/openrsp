! Copyright 2012      ?
!           2010-2012 Andreas J. Thorvaldsen
! This file is made available under the terms of the
! GNU Lesser General Public License version 3.

!> @file Contains module cgto_diff_eri

!> This module calculates differentiated electron repulsion integrals (diff_eri's)
!> over generally contracted solid-Harmonic Gaussian basis functions (cgto's).
!>
!>    geodiff_eri
!>          Geometry-differentiated 2e-integrals (d^n/dRk^n)
!>
module cgto_diff_eri

#ifndef OPENRSP_STANDALONE
  use dalton_ifc, only: boys_function
#endif /* OPENRSP_STANDALONE */
  use basis_set
  use eri_contractions
  ! use cgto_diff_eri_print
  implicit none
  public geodiff_eri
  public cgto
  real(8), parameter :: pi = acos(-1d0)
  private

contains


  !> Compute the geo'th order geometry differentiated eri's for basis
  !> quadrouple ABCD, and pass computed integrals to subroutine callback
  subroutine geodiff_eri(geo, A, B, C, D, ctrs)
    !> the four shells of basis functions
    type(cgto),    intent(in)    :: A, B, C, D
    !> order of geometry-differentiation
    integer,       intent(in)    :: geo
    !> vector of arguments to be passed to contraction routine
    type(ctr_arg)                :: ctrs(:)
    !------------------------------
    integer :: rec(4, A%mom + B%mom + C%mom + D%mom + geo  &
                    + (2*C%mom + 2*D%mom + geo)*(geo+1)/2)
    integer :: ma, mb, ep, lab, rab, mp, geoab, mt, i
    integer :: mc, md, eq, lcd, rcd, mq, geocd, j, k
    integer :: ht, hp, hab, gab, hq, hcd, gcd, lwork
    ! short-hands for dimensions
    ma = A%mom; mb = B%mom
    mc = C%mom; md = D%mom
    ep = size(A%exp) * size(B%exp)
    eq = size(C%exp) * size(D%exp)
    rab = size(A%ctr,1) * size(B%ctr,1)
    rcd = size(C%ctr,1) * size(D%ctr,1)
    lab = (2*ma+1) * (2*mb+1)
    lcd = (2*mc+1) * (2*md+1)
    mt = ma+mb+mc+md+geo
    ht = (mt+1)*(mt+2)*(mt+3)/6
    mp = ma+mb+geo
    hp = (mp+1)*(mp+2)*(mp+3)/6
    ! call inner twice. First call is to determine lwork,
    ! second to carry out integral calculation (and callback)
    lwork = 0
    call inner(.true.)   !dry run: measure up, dermine lwork
    call inner(.false.)  !actual run: calculate, contract

  contains

    subroutine inner(dryrun)
      logical, intent(in) :: dryrun
      real(8) :: work(eq*lwork)
      integer :: irec, lwork2, ixfer, n, n1, n2, n3, nab
      ! first calculate 1-electron 1-center Coulomb integrals on the 
      ! "P-Q" center, then transfer all moment to center P in work(:n)
      n = ht * ep !size
      if (dryrun) then
         lwork = merge(n,2*n,mp==0)
      else if (mp/=0) then
         call coulomb_PmQcdab(mt, A, B, C, D, work(eq*n+1:eq*n*2))
         call elxfer_Rqp_to_qPp(mt, A, B, C, D, work(eq*n+1:eq*n*2), work(:eq*n))
      else !if (mp==0) then
         call coulomb_PmQcdab(mt, A, B, C, D, work(:eq*n))
         call elxfer_Rqp_to_Qqp(mt, A, B, C, D, work(:eq*n))
      end if
      ! first electron's recursion, P => PAB
      if (dryrun .and. mp/=0) then
         rec(:,:mt) = rec_layout(ma, mb, geo, mc+md)
         hab = (ma+1)*(ma+2)/2 * (mb+1)*(mb+2)/2
         lwork = max(lwork, rec(3,mt) + rec(4,mt) &
                          + (ep-1) * max(0, hp*hab - ht))
      else if (mp/=0) then
         call recur_hpep_to_habep(A, B, geo, mc+md, &
                                  rec(:,:mt), eq, lwork, work, &
                                  C%exp, D%exp, ixfer)
      end if
      ! proceed with blocks of radial AB, angular AB, reorder AB,
      ! recursion Q => CD, radial CD, angular CD, reorder CD,
      ! callback for contraction
      irec = 0 !offset into rec
      mq = mt  !just for pre-increment of irec
      do geoab=0, geo
         geocd = geo - geoab
         irec = irec + mq !pre-increment
         mq = mc + md + geocd
         hq = (mq+1)*(mq+2)*(mq+3)/6
! if (dryrun) print *, 'rec1 lwork =', lwork
! if (dryrun) print *
         ! number of Hermite pairs on AB in this batch
         hab = sum((/((ma+geoab-i+1)*(ma+geoab-i+2)/2 &
                          * (mb+i+1)*(mb+i+2)/2, i=0,geoab)/))
         ! radial contraction for A and B
         n3 = hq*hab*ep
         n2 = hq*hab*(size(A%exp)+1) !scratch
         n1 = hq*rab*hab
         if (dryrun) then
            lwork = lwork + n3
            lwork = max(lwork, n1+n2+n3)
            if (mp==0) lwork = max(lwork, n1+n1)
         else if (mp/=0) then
            i = eq*ixfer
            call radial_habep_to_rabhab(hq*eq, A, B, geoab, hab, &
                                work(i+1:i+eq*n3), work(:eq*n1), &
                                work(i-eq*n2+1:i))
            ixfer = ixfer + n3
         else if (mp==0) then
            i = eq * max(n3+n2,n1)
            ! hqeqep instead resides at bottom. Do radial, then trivial angular
            call radial_habep_to_rabhab(hq*eq, A, B, 0, 1, &
                                work(:eq*n3), work(i+1:i+eq*n1), &
                                work(i-eq*n2+1:i))
            work(:eq*n1) = 1/(4*pi) * reshape(transpose(reshape( &
                           work(i+1:i+eq*n1),(/hq*eq,rab/))),(/rab*hq*eq/))
         end if
! if (dryrun) print *, 'rad1 lwork =', lwork
! if (dryrun) print *
         ! number of AB-diff pairs in this batch
         gab = sum((/((geoab-i+1)*(geoab-i+2)/2 &
                          * (i+1)*(i+2)/2, i=0, geoab)/))
         ! angular-geometric contraction for A and B, then final reorder
         n2 = hq*rab*(2*ma+1)*((geoab+1)*(geoab+2)/2) !scratch
         n3 = hq*rab*lab*gab
         if (mp/=0 .and. dryrun) then
            lwork = max(lwork, max(n1+n2+n3, n3+n3))
         else if (mp/=0) then
            i = eq * max(n1+n2,n3)
            call angular_hab_to_lagablb( &
                              ma, cart_to_spher_coef(ma), 0,   &
                              mb, cart_to_spher_coef(mb), 0,   &
                              geoab, hab, gab, hq*eq*rab,      &
                              work(:eq*n1), work(i+1:i+eq*n3), &
                              work(i-eq*n2+1:i))
            call reorder_rablagablb_to_gabralarblb(         &
                              size(A%ctr,1), (2*ma+1),      &
                              size(B%ctr,1), (2*mb+1), gab, &
                              hq*eq, work(i+1:i+eq*n3), work(:n3*eq))
         end if
! if (dryrun) print *, 'ang1 lwork =', lwork
! if (dryrun) print *
         ! leading dim for 2nd electron
         nab = gab*rab*lab
         ! recursion layout for this batch of Q => CD
         if (mq/=0 .and. dryrun) &
            rec(:,irec+1:irec+mq) = rec_layout(mc, md, geocd, 0)
         ! number of CD Hermite pairs in this batch
         hcd = sum((/((mc+geocd-i+1)*(mc+geocd-i+2)/2 &
                          * (md+i+1)*(md+i+2)/2, i=0,geocd)/))
         ! this batch of Q => CD recursion, then radial contraction
         n3 = hcd*eq
         n2 = hcd*(size(C%exp)+1) !scratch
         n1 = rcd*hcd
! if (mq/=0.and.dryrun) print *, 'rec2 lwork =', max(lwork, (nab*(rec(3,irec+mq) &
!                                    + rec(4,irec+mq) - hcd + n3)-1)/eq+1)
! if (mq/=0.and.dryrun) print *
         lwork2 = n1+n2+n3
         if (mq/=0) lwork2 = max(lwork2, rec(3,irec+mq) &
                                       + rec(4,irec+mq) + max(hq,hcd)*(eq-1))
         if (mq==0) lwork2 = max(lwork2, n1+n1)
         if (mq/=0 .and. .not.dryrun) then
            call recur_hpep_to_habep(C, D, geocd, 0, &
                                     rec(:,irec+1:irec+mq), nab, &
                                     lwork2, work(:nab*lwork2))
            i = nab * (lwork2-n3)
            call radial_habep_to_rabhab(nab, C, D, geocd, hcd, &
                                work(i+1:i+nab*n3), work(:nab*n1), &
                                work(i-nab*n2+1:i))
         else if (mq==0 .and. .not.dryrun) then
            i = nab * max(n3+n2,n1)
            ! gablaralbrbeq already resides at bottom, do radial
            call radial_habep_to_rabhab(nab, C, D, 0, 1, &
                                work(:nab*n3), work(i+1:i+nab*n1), &
                                work(i-nab*n2+1:i))
            ! then trivial angular, scale-copy
            work(:nab*n1) = 1/(4*pi) * work(i+1:i+nab*n1)
         end if
! if (dryrun) print *, 'rad2 lwork =', max(lwork, (nab*lwork2-1)/eq+1)
! if (dryrun) print *
         ! number of CD-diff pairs in this batch
         gcd = sum((/((geocd-i+1)*(geocd-i+2)/2 &
                          * (i+1)*(i+2)/2, i=0, geocd)/))
         ! angular-geometric contraction for CD, then final reorder
         n2 = rcd*(2*mc+1)*((geocd+1)*(geocd+2)/2) !scratch
         n3 = rcd*lcd*gcd
         if (mq/=0 .and. dryrun) then
            lwork2 = max(lwork2, max(n1+n2+n3, n3+n3))
         else if (mq/=0) then
            i = nab * max(n1+n2,n3)
            call angular_hab_to_lagablb( &
                              mc, cart_to_spher_coef(mc), 0, &
                              md, cart_to_spher_coef(md), 0, &
                              geocd, hcd, gcd, nab*rcd, work(:nab*n1), &
                              work(i+1:i+nab*n3), work(i-nab*n2+1:i))
            call reorder_rcdlcgcdld_to_gcdrclcrdld(    &
                              size(C%ctr,1), (2*mc+1), &
                              size(D%ctr,1), (2*md+1), &
                              gcd, gab, rab*lab,       &
                              work(i+1:i+nab*n3), work(:nab*n3))
         end if
         ! finally either update lwork or call callback routine 
         if (dryrun) then
            ! adjust lwork so (eq,lwork) extends past (nab,lwork2)
            lwork = max(lwork, (nab*lwork2-1)/eq+1)
         else
            ! call contraction routine
            call contract_with_density(A, B, C, D, geoab, geocd, gab*gcd, &
                                       work(:nab*n3), ctrs)
         end if
! if (dryrun) print *, 'ang2 lwork =', lwork
! if (dryrun) print *
      end do
    end subroutine

  end subroutine




  subroutine recur_hpep_to_habep(A, B, geo, mq, rec, &
                                 ldim, lwork, work, &
                                 ec, ed, ixfer)
    !> the two shells of basis functions
    type(cgto), intent(in) :: A, B
    !> order of geometry-differentiation
    integer,    intent(in) :: geo
    !> excess momentum on center P, to be electron-transferred to
    !> center Q in the result (ec, ed and ixfer must be present if mq/=0)
    integer,    intent(in) :: mq
    !> memory layout for recursion
    integer,    intent(in) :: rec(4, A%mom + B%mom + geo + mq)
    !> leading dimension
    integer,    intent(in) :: ldim
    !> length of work
    integer,    intent(in) :: lwork
    !> work array, with input integrals hPeP at the beginning
    real(8), intent(inout) :: work(ldim, lwork)
    !> exponents on centers C and D, used in electron-transfer
    real(8), optional, intent(in)  :: ec(:), ed(:)
    !> offset into work where electron-transferred results are placed
    integer, optional, intent(out) :: ixfer
    !-----------------------------------------------------------
    real(8) :: ooq(ldim) !only used when ixfer present
    integer :: ma, mb, mp, i, j, g, m
    integer :: ep, hp, ip, hab, q, lw, ix(geo+1)
    ! point hP offset to the end of hPeP (last pair a b)
    ma = A%mom
    mb = B%mom
    mp = ma + mb + geo + mq        !momentum on P in input
    hp = (mp+1)*(mp+2)*(mp+3)/6    !number of hermite components on P
    ep = size(A%exp) * size(B%exp) !number of exponent pairs
    ! set pointer to P to (past) the end of PhPe
    ip = hp * ep
    ! if present, prepare transfer factors and offsets
    if (present(ixfer)) then
       ! layout transfer offsets
       ix(geo+1) = lwork - ip
       do g = geo, 1, -1
          q = mq + geo - g
          ix(g) = ix(g+1) - (q+1)*(q+2)*(q+3)/6  &
                * sum((/((ma+g-m+1)*(ma+g-m+2)/2 &
                       * (mb+m+1)*(mb+m+2)/2, m=0,g)/)) * ep
       end do
       ! precalculate transfer factors
       ooq = (/(1/(ec + ed(i)), i=1, size(ed))/)
    else
       hab = sum((/((ma+geo-g+1)*(ma+geo-g+2)/2 &
                  * (mb+g+1)*(mb+g+2)/2, g=0,geo)/))
    end if
! call print_tensor((/ldim,hp,ep/), 'qPp', 'iai', work(:,:ip), &
!                   title = '###### initial')
    ! loop backwards over exponent pairs Aexp Bexp
    if (.not.present(ixfer)) lw = lwork - ip + hab
    do j = size(B%exp), 1,-1
       do i = size(A%exp), 1,-1
          ! pre-decrement P offset, and increment transfer offsets,
          ! then carry out recursion
          ip = ip - hp
          if (present(ixfer)) then
             ix(:) = ix(:) + hp
             call recur_P_to_AB(ma, mb, geo, mq,      &
                                A%cent - B%cent,      &
                                A%exp(i), B%exp(j),   &
                                rec, ldim, ix(geo+1), &
                                work(:,ip+1:ip+ix(geo+1)), &
                                ooq, ix)
          else
             lw = lw + hp - hab
             call recur_P_to_AB(ma, mb, geo, mq, &
                                A%cent - B%cent, &
                                A%exp(i), B%exp(j), rec, &
                                ldim, lw, work(:,ip+1:ip+lw))
          end if
       end do
    end do
    ! return offset of electron-transferred results
    if (present(ixfer)) ixfer = ix(1)
  end subroutine




  !> horizontal recursion of moment from product center P to cgto
  !> centers A and B
  subroutine recur_P_to_AB(ma, mb, geo, mq, AmB, ea, eb, rec, &
                           ldim, lwork, work, ooq, ixfer)
    !> moment on A, moment on B, order of geometry-
    !> differentiation, residual moment on P to be
    !> electron-transferred to Q, leading dim, length of work
    integer, intent(in)    :: ma, mb, geo, mq, ldim, lwork
    !> center displacement A-B, exponent on A, exponent on B
    real(8), intent(in)    :: AmB(3), ea, eb
    !> memory layout to be used in recursion
    integer, intent(in)    :: rec(4, ma + mb + geo + mq)
    !> recursion and result buffer
    real(8), intent(inout) :: work(ldim,lwork)
    !> electron-transfer factors
    real(8), intent(in),    optional :: ooq(ldim)
    !> offset into work at which to place electron-transferred results
    integer, intent(inout), optional :: ixfer(geo+1)
    !--------------------------------------------
    real(8) :: oo2p, AmP(3), BmP(3), bo2ap, ao2bp, poq(ldim)
    integer :: a, b, p, ab, q, minb, maxb, mina, maxa
    integer :: i, j, k, l, m, n, g
    integer :: ipab, iab, ia, iabb
    integer :: npab, nab, na, nabb, lpab, labb, np, npp, nq
    integer :: basoff, dx1(geo+1), dx3(geo+1), d1, d3
    logical :: odd
    ! calculate factors used in recursion
    oo2p  = 1/(2*(ea+eb))
    AmP   =  2 * eb * oo2p * AmB
    BmP   = -2 * ea * oo2p * AmB
    bo2ap = eb * oo2p / ea
    ao2bp = ea * oo2p / eb
    ! if 1st electron calculate electron-transfer factors p/q
    if (present(ixfer)) poq = (ea+eb)*ooq
    ! if 2nd electron (no ixfer) set base offset so that result
    ! of final recursion will end up at the end of work
    basoff = merge(0, lwork - rec(3,ma+mb+geo+mq) &
                            - rec(4,ma+mb+geo+mq), present(ixfer))
    ! if no momentum on A and B, also P00 blocks should be transferred
    if (present(ixfer) .and. ma==0 .and. mb==0) then
       q = mq+geo
       nq = (q+1)*(q+2)*(q+3)/6
       ixfer(1) = ixfer(1) - nq
       i = nq
       do p=q,0,-1
          np = (p+1)*(p+2)/2
          i = i - np
          call elxfer_QePh_to_QhQe(ldim, p, nq, 1, poq,  &
                                   work(:,i+1:i+np), &
                                   work(:,ixfer(1)+1 : ixfer(1)+nq))
       end do
    end if
    ! loop over descending momentum on P, increasing on AB
    do k=1, ma + mb + geo + mq
       p = ma + mb + geo + mq - k
       np = (p+1)*(p+2)/2
       npp = (p+2)*(p+3)/2
       ! momenta for this moment p
       a = min(ma, k)
       b = min(mb, max(0,k-ma))
       g = min(geo,max(0,k-ma-mb))
       ! initialize offsets and lengths
       iab  = (p+0)*(p+1)*(p+2)/6 !the zero block
       nab  = 1
       ia   = iab
       na   = 0
       ipab = iab + np !the zero block
       npab = 1
       iabb = rec(3,k) + basoff
       nabb = rec(4,k)
       ! loop over diagonals a+b
       do ab=1,a+b+g
          ! after ab=1, PAB should start on its head part
          if (ab==2) ipab = rec(3,k-1) + basoff
          if (ab==2) npab = rec(4,k-1)
          ! if head of PAB empty, start on tail
          if (npab == 0) ipab = rec(1,k-1) + basoff
          if (npab == 0) npab = rec(2,k-1)
          ! if head of ABB empty, start on tail
          if (nabb == 0) iabb = rec(1,k) + basoff
          if (nabb == 0) nabb = rec(2,k)
          ! next "anti-diagonal" of integrals
          minb = max(0, min(ab-a,(ab-a+b-g+1)/2))
          maxb = max(0, min(ab,  (ab-a+b+g)/2))
          mina = ab - maxb
          maxa = ab - minb
          ! sometimes edge-blocks in (P)AB should be skipped
          if (maxa-minb == a-b+g-1 .and. maxa>ma .and. minb/=0) then
             nab = nab - d1
             iab = iab + np * d1
             if (p<mq) npab = npab - d1
             if (p<mq) ipab = ipab + npp * d1 
          end if
          if (mina-maxb == a-b-g+1 .and. mina/=0 .and. maxb/=0) then
             nab = nab - d3
          end if
          ! replicate length of PAB, which will be adjusted by recur
          lpab = npab
          ! if PAB blocks are above ABB, make sure they don't overlap
          labb = merge((ipab-iabb)/np, nabb, ipab > iabb)
          call recur_P_to_AB_adia(ldim, p, ma, ab, minb, maxb,        &
                                  AmP, BmP, bo2ap, ao2bp, oo2p,       &
                                  lpab, work(:,ipab+1:ipab+npp*lpab), &
                                  nab,  work(:,iab +1:iab + np*nab ), &
                                  na,   work(:,ia  +1:ia  + np*na  ), &
                                  labb, work(:,iabb+1:iabb+ np*labb))
          ! reverse skipping of edge blocks
          if (maxa-minb == a-b+g-1 .and. maxa>ma .and. minb/=0) then
             nab = nab + d1
             iab = iab - np * d1
          end if
          if (mina-maxb == a-b-g+1 .and. mina/=0 .and. maxb/=0) then
             nab = nab + d3
             if (p<mq) lpab = lpab + d3
          end if
          ! update 'fringe' sizes
          d1 = (maxa+1)*(maxa+2)/2 * (minb+1)*(minb+2)/2
          d3 = (mina+1)*(mina+2)/2 * (maxb+1)*(maxb+2)/2
          ! if this is the first electron (ooq present), the blocks
          ! a>=ma, b>=mb should be electron-transferred (scaled by poq)
          if (present(ixfer) .and. ab>=ma+mb) then
             i = ab-ma-mb+1 !index in ixfer for this ab
             q = mq+geo-i+1 !total moment on q for this ab
             nq = (q+1)*(q+2)*(q+3)/6
             ! if first transfer of this order, pre-decrement
             ! transfer offset and initialize dx1 and dx3 to zero
             if (p>=mq .and. ab==a+b+g) then
                ixfer(i) = ixfer(i) - nq * labb
                dx1(i) = 0
                dx3(i) = 0
             else if (p>=mq) then
                if (maxa-minb == a-b+g) dx1(i) = dx1(i) + d1
                if (mina-maxb == a-b-g) dx3(i) = dx3(i) + d3
             end if
             ! electron transfer
             call elxfer_QePh_to_QhQe(ldim, p, nq, &
                             labb - dx1(i) - dx3(i), poq,         &
                             work(:,iabb + np * dx1(i) + 1        &
                                  : iabb + np * (labb - dx3(i))), &
                             work(:,ixfer(i)+1 : ixfer(i)         &
                                  + nq * (labb - dx1(i) - dx3(i))))
          end if
          ! 
          ! progress offsets and lengths
          ipab = ipab + npp * lpab
          npab = npab - lpab
          ia   = iab
          na   = nab
          iab  = iabb
          nab  = nabb
          iabb = iabb + np * labb
          nabb = nabb - labb
       end do
    end do
  end subroutine




  !> electron transfer from Hermite center P, the first electron's product
  !> center, to Q, the second electron's product center
  subroutine elxfer_QePh_to_QhQe(ne, mp, nQh, tdim, poq, PePh, QhQe)
    !> number of exponents ne, moment on P, no. Hermites on Q, trail dim
    integer, intent(in)    :: ne, mp, nQh, tdim
    !> electron transfer factors p/q = (a+b)/(c+d)
    real(8), intent(in)    :: poq(ne)
    !> input integrals with moment on P
    real(8), intent(in)    :: PePh(ne, (mp+1)*(mp+2)/2, tdim)
    !> output integrals with moment moved to Q, and transposed
    real(8), intent(inout) :: QhQe(nQh, ne, tdim)
    !--------------------------------------------
    integer :: i, j, k, n
    ! offset and size into QhQe
    i = (mp+0)*(mp+1)*(mp+2)/6
    n = size(PePh,2) !(mp+1)*(mp+2)/2
    do k=1,tdim
       do j=1,ne
          QhQe(i+1:i+n,j,k) = (-poq(j))**mp * PePh(j,:,k)
       end do
    end do
  end subroutine



  !> Carries out ABB = PAB + (B-P)AB + A + B
  !> ajt FIXME generalize to magnetic field
  subroutine recur_P_to_AB_adia(ldim, mp, Amom, mab, minb, maxb, &
                                AmP, BmP, bo2ap, ao2bp, oo2p,    &
                                lpaab, PPAAB, laab,  PAAB,       &
                                laa,   PAA,   laabb, PAABB)
    integer, intent(in)    :: Amom, minb, maxb, mp, mab, ldim
    real(8), intent(in)    :: AmP(3), Bmp(3), bo2ap, ao2bp, oo2p
    integer, intent(in)    :: laab, laa    !lengths, unchanged
    integer, intent(inout) :: lpaab, laabb !lengths, corrected on exit
    real(8), intent(in)    :: PPAAB(ldim, (mp+2)*(mp+3)/2, lpaab)
    real(8), intent(in)    :: PAAB( ldim, (mp+1)*(mp+2)/2, laab)
    real(8), intent(in)    :: PAA(  ldim, (mp+1)*(mp+2)/2, laa)
    real(8), intent(out)   :: PAABB(ldim, (mp+1)*(mp+2)/2, laabb)
    integer :: ma, mb, a0, b0, iaab, iaa, iab, iaabb
    integer :: np, na, nb, naab, naa, nab, naabb
    logical :: doa
    np = (mp+1)*(mp+2)/2
    iaab=0; iaa=0; iab=0; iaabb=0
    do mb = minb, maxb
       ma = mab-mb
       na = (ma+1)*(ma+2)/2
       nb = (mb+1)*(mb+2)/2
       naabb = na * nb
       doa = (mb==0 .or. ma > Amom)
       b0 = mb - merge(0,1,doa)
       a0 = ma - merge(1,0,doa)
       naab = (a0+1)*(a0+2)/2 * (b0+1)*(b0+2)/2
       naa  = (a0+1)*(a0+2)/2 * (b0+0)*(b0+1)/2
       nab  = (a0+0)*(a0+1)/2 * (b0+1)*(b0+2)/2
! print '(a,3(a1,i1)/)', '######## going to make ', 'P', mp, 'A', ma, 'B', mb
! print '(a,3f10.5/a,3f10.5)','A-P', AmP, 'B-P', BmP
! print '(3(a,f10.5/))', 'bo2ap ', bo2ap, 'ao2bp ', ao2bp, 'oo2p  ', oo2p
! call print_tensor((/ldim,(mp+2)*(mp+3)/2,(a0+1)*(a0+2)/2,(b0+1)*(b0+2)/2/), &
!                   'IPAB', 'iccc', PPAAB(:,:,iaab+1:iaab+naab), &
!                   title = '######## source top')
! call print_tensor((/ldim,(mp+1)*(mp+2)/2,(a0+1)*(a0+2)/2,(b0+1)*(b0+2)/2/), &
!                   'IPAB', 'iccc', PAAB(:,:,iaab+1:iaab+naab), &
!                   title = '######## source bottom')
       ! ajt In this call, 'ifort -check all' warns 'array temporary was created' (x3).
       !     The warning disappears if I join the two first dimensions, thus putting
       !     '(:,' instead of '(:,:,'. This must be a bug in ifort 12.1. If the same
       !     turns out for later versions, this should probably be worked around
       call recur_translate(ldim, mp, merge(1,na,doa),     &
                            merge(ma,mb,doa),              &
                            merge(nb, 1, doa .and. mb/=0), &
                            merge(AmP,BmP,doa),            &
                            PPAAB(:,:,iaab +1:iaab +naab), &
                            PAAB( :,:,iaab +1:iaab +naab), &
                            PAABB(:,:,iaabb+1:iaabb+naabb))
! if (doa .and. mb/=0) &
! call print_tensor((/ldim,(mp+1)*(mp+2)/2,(ma+0)*(ma+1)/2,(mb+0)*(mb+1)/2/), &
!                   'IPAB', 'iccc', PAA(:,:,iaa+1:iaa+naa), &
!                   title = '######## other A')
       if (doa .and. mb/=0) &
          call recur_commute_other_A(ldim * np, ma, mb, oo2p, &
                            PAA(  :,:,iaa  +1:iaa  +naa),     &
                            PAABB(:,:,iaabb+1:iaabb+naabb))
! if (.not.doa .and. mb/=1) &
! call print_tensor((/ldim,(mp+1)*(mp+2)/2,(ma+1)*(ma+2)/2,(mb-1)*(mb+0)/2/), &
!                   'IPAB', 'iccc', PAA(:,:,iaa+1:iaa+naa), &
!                   title = '######## same B')
       if (.not.doa .and. mb/=1) &
          call recur_commute_same(ldim * np * na, mb, 1, ao2bp, &
                            PAA(  :,:,iaa  +1:iaa  +naa),       &
                            PAABB(:,:,iaabb+1:iaabb+naabb))
! if (doa .and. ma/=1) &
! call print_tensor((/ldim,(mp+1)*(mp+2)/2,(ma-1)*(ma+0)/2,(mb+1)*(mb+2)/2/), &
!                   'IPAB', 'iccc', PAA(:,:,iaa+naa+1:iaa+naa+nab), &
!                   title = '######## same A')
       if (doa .and. ma/=1) &
          call recur_commute_same(ldim * np, ma, nb, bo2ap,   &
                            PAA(  :,:,iaa+naa+1:iaa+naa+nab), &
                            PAABB(:,:,iaabb  +1:iaabb+naabb))
! if (.not.doa .and. ma/=0) &
! call print_tensor((/ldim,(mp+1)*(mp+2)/2,(ma+0)*(ma+1)/2,(mb+0)*(mb+1)/2/), &
!                   'IPAB', 'iccc', PAA(:,:,iaa+naa+1:iaa+naa+nab), &
!                   title = '######## other B')
       if (.not.doa .and. ma/=0) &
          call recur_commute_other_B(ldim * np, ma, mb, oo2p, &
                            PAA(  :,:,iaa+naa+1:iaa+naa+nab), &
                            PAABB(:,:,iaabb  +1:iaabb+naabb))
! call print_tensor((/ldim,(mp+1)*(mp+2)/2,(ma+1)*(ma+2)/2,(mb+1)*(mb+2)/2/), &
!                   'IPAB', 'iccc', PAABB(:,:,iaabb+1:iaabb+naabb), &
!                   title = '######## result')
       ! pause if next block from same PPAAB PAAB PAA
       if (.not.(doa .and. mb /= maxb .and. ma <= Amom+1)) then
          iaab = iaab + naab 
          iaa  = iaa  + naa
       end if
       iaabb = iaabb + naabb
    end do
    ! move past final PAA block
    if (maxb /= 0) iaa = iaa + nab
    ! reduce to correct lengths
    lpaab = iaab
    laabb = iaabb
  end subroutine



  subroutine recur_translate(ldim, mp, mdim, ma, tdim, AmP, PPA, PA, PAA)
    integer, intent(in)  :: ldim, mp, mdim, ma, tdim
    real(8), intent(in)  :: AmP(3)
    real(8), intent(in)  :: PPA(ldim, (mp+2)*(mp+3)/2, mdim, &
                                      (ma+0)*(ma+1)/2, tdim)
    real(8), intent(in)  ::  PA(ldim, (mp+1)*(mp+2)/2, mdim, &
                                      (ma+0)*(ma+1)/2, tdim)
    real(8), intent(out) :: PAA(ldim, (mp+1)*(mp+2)/2, mdim, &
                                      (ma+1)*(ma+2)/2, tdim)
    integer :: i, j, k, n, m
    i = 0; j = 1
    do n = mp+1,1,-1
       ! Axxxx
       PAA(:,i+1:i+n,:,1,:) = PPA(:,j+0:j+n-1,:,1,:) &
                   - AmP(1) *  PA(:,i+1:i+n,  :,1,:)
       ! Axxxy ... Ayyyy
       PAA(:,i+1:i+n,:,2:ma+1,:) = PPA(:,j+1:j+n,:,1:ma,:) &
                        - AmP(2) *  PA(:,i+1:i+n,:,1:ma,:)
       i = i+n; j = j+n+1
    end do
    ! Axxxz ... Ayyyz; ...; Azzzz
    PAA(:,:,:,ma+2:,:) = PPA(:,mp+3:,:,:,:) - AmP(3) * PA(:,:,:,:,:)
  end subroutine



  subroutine recur_commute_same(ldim, ma, tdim, bo2ap, PA, PAAA)
    integer, intent(in)    :: ldim, ma, tdim
    real(8), intent(in)    :: bo2ap
    real(8), intent(in)    :: PA(  ldim, (ma-1)*(ma+0)/2, tdim)
    real(8), intent(inout) :: PAAA(ldim, (ma+1)*(ma+2)/2, tdim)
    integer :: i, j, n
    PAAA(:,1,:) = PAAA(:,1,:) - (ma-1) * bo2ap * PA(:,1,:)
    do i=1,ma-1
       PAAA(:,i+2,:) = PAAA(:,i+2,:) - i * bo2ap * PA(:,i,:)
    end do
    i = 0
    do n = ma-1,1,-1
       j = i + 2*ma+1
       PAAA(:,j+1:j+n,:) = PAAA(:,j+1:j+n,:) &
                         - (ma-n) * bo2ap * PA(:,i+1:i+n,:)
       i = i+n
    end do
  end subroutine



  subroutine recur_commute_other_A(ldim, ma, mb, oo2p, PAB, PAABB)
    integer, intent(in)    :: ldim, ma, mb
    real(8), intent(in)    :: oo2p
    real(8), intent(in)    :: PAB(  ldim,(ma+0)*(ma+1)/2, &
                                         (mb+0)*(mb+1)/2)
    real(8), intent(inout) :: PAABB(ldim,(ma+1)*(ma+2)/2, &
                                         (mb+1)*(mb+2)/2)
    integer :: i, j, n, m
    i = 1
    do n=mb,1,-1
       do m=1,n
          j = i+mb+1
          PAABB(:,1,j-n-1) = PAABB(:,1,j-n-1) &
                           +  (n-m+1) * oo2p * PAB(:,1,i)
          PAABB(:,2:ma+1,j-n) = PAABB(:,2:ma+1,j-n) &
                              +     m * oo2p * PAB(:,:ma,i)
          PAABB(:,ma+2:,j) = PAABB(:,ma+2:,j) &
                           + (mb-n+1) * oo2p * PAB(:,:,i)
          i = i+1
       end do
    end do
  end subroutine



  subroutine recur_commute_other_B(ldim, ma, mb, oo2p, PAB, PAABB)
    integer, intent(in)    :: ldim, ma, mb
    real(8), intent(in)    :: oo2p
    real(8), intent(in)    :: PAB(  ldim, (ma+0)*(ma+1)/2, &
                                          (mb+0)*(mb+1)/2)
    real(8), intent(inout) :: PAABB(ldim, (ma+1)*(ma+2)/2, &
                                          (mb+1)*(mb+2)/2)
    integer :: i, mxy, mx, my, j, n
    ! Bxxxx
    i = 0
    do mxy = ma, 1, -1
       do mx = mxy, 1, -1
          j = i+ma-mxy
          PAABB(:,j+1,1) = PAABB(:,j+1,1) &
               + mx*oo2p *   PAB(:,i+1,1)
          i = i + 1
       end do
    end do
    ! Bxxxy ... Byyyy
    i = 0
    do mxy = ma, 1, -1
       do my = 1, mxy
          j = i+ma-mxy+1
          PAABB(:,j+1,2:mb+1) = PAABB(:,j+1,2:mb+1) &
                    + my*oo2p *   PAB(:,i+1,1:mb)
          i = i + 1
       end do
    end do
    ! Bxxxz ... Byyyz; ...; Bzzzz
    i = 0
    do n = ma,1,-1
       j = i+ma+1
       PAABB(:,j+1:j+n,mb+2:) = PAABB(:,j+1:j+n,mb+2:) &
              + (ma+1-n)*oo2p *   PAB(:,i+1:i+n,:)
       i = i+n
    end do
  end subroutine



  !> Coulomb integrals over Hermite Gaussians with exponent expo,
  !> up to order/moment mom, centered at cent, i.e.
  !> derivatives wrt. cent up to order mom, divided by (2*expo)^n of:
  !>    integral [ exp(-expo*|r-cent|^2) / |r| dr ]
  subroutine hermite_coulomb(mom, expo, cent, erint, prefac)
    !> Gaussian exponent at center Cent
    real(8), intent(in)  :: expo
    !> center of basis functions
    real(8), intent(in)  :: cent(3)
    !> (highest) moment of basis
    integer, intent(in)  :: mom
    !> resulting Coulomb (repulsion) integrals
    real(8), intent(out) :: erint((mom+1)*(mom+2)*(mom+3)/6)
    !> prefactor, default 1
    real(8), optional, intent(in) :: prefac
    !-------------------------------------
    integer :: n, i, j, k, m, x, z, zz, xy
    real(8) :: er2, emer2, tmp
    er2 = expo * sum(cent**2)
    emer2 = exp(-er2)
#ifdef OPENRSP_STANDALONE
    print *, 'error: not part of standalone'
    stop 1
#else /* OPENRSP_STANDALONE */
    erint(1) = 2*pi / expo * boys_function(mom, er2, emer2)
#endif /* OPENRSP_STANDALONE */
    emer2    = 2*pi / expo * emer2
    if (present(prefac)) emer2    = emer2 * prefac
    if (present(prefac)) erint(1) = erint(1) * prefac
    ! first step: fill triangles with decending orders of the
    !             Boys function differentiated by x and y
    i = 1
    j = 2
    do m=1,mom
       do xy=m,1,-1
          erint(j:j+xy-2) = -cent(1) * erint(i:i+xy-2) &
                        - 1/(2*expo) * (/(x-1,x=xy,2,-1)/) &
                                     * erint(i+xy:i+2*xy-2)
          erint(j+xy-1) = -cent(1) * erint(i+xy-1)
          erint(j+xy)   = -cent(2) * erint(i+xy-1) &
               - 1/(2*expo)*(xy-1) * erint(i+2*xy-2)
          i = i + xy
          j = j + xy + 1
       end do
       erint(j) = (2*er2 * erint(i-1) + emer2) / (2*(mom-m)+1)
       j = j + 1
    end do
    ! time-consuming part: xyz-integrals from xy-integrals
    !ajt fixme Small/zero Cz will cause overflow/infinity here.
    !          Perhaps avoidable by using 1 instead of Cz in some of the
    !          early iterations of the loop
    ! ajt temporarily changed to prevent overflow for Cz==0
    !   do m=2,mom
    do m=1,mom
       i = 0
       j = 3
       n = 1
       ! ajt temporarily changed to prevent overflow for Cz==0
       !   do z=m,2,-1
       do z=m,1,-1
          ! ajt temporarily changed to prevent overflow for Cz==0
          !   erint(i+1:i+n) = erint(i+1:i+n) &
          !                  - erint(j+1:j+n) * (z-1) &
          !                  / (2*expo * cent(3)**2)
          erint(i+1:i+n) = -cent(3) * erint(i+1:i+n) &
                   - (z-1)/(2*expo) * erint(j+1:j+n)
          i = i + n
          j = j + n + m - z + 3
          n = n + m - z + 2
       end do
    end do
    ! multiply integrals by (-Cz)^mz while reordering
    j = mom*(mom+1)*(mom+2)/6
    do z=0,mom/2
       n = mom+1 - 2*z
       ! ajt temporarily changed to prevent overflow for Cz==0
       !   erint(j+1:j+n) = erint(j+1:j+n) * (-cent(3))**z
       i = j
       do zz=z+1,mom-z
          j = j + n
          i = i - (mom+1-zz)*(mom+2-zz)/2 - z
          n = n - 1
          ! swap an scale rows
          do k=1,n
             tmp = erint(j+k)
             ! ajt temporarily changed to prevent overflow for Cz==0
             !   erint(j+k) = erint(i+k) * (-cent(3))**zz
             !   erint(i+k) = tmp * (-cent(3))**z
             erint(j+k) = erint(i+k)
             erint(i+k) = tmp
          end do
       end do
       j = j - (mom+1-z  )*(mom+2-z  )/2 &
             - (mom-1-z*2)*(mom+0-z*2)/2 + 1
    end do
  end subroutine




  !> Determine circular layout for horizontal recursion
  function rec_layout(ma, mb, geo, mq) result(r)
    !> moment on center A
    integer,   intent(in) :: ma
    !> moment on center B
    integer,   intent(in) :: mb
    !> order of geometry differentiation wrt A and B
    integer,   intent(in) :: geo
    !> moment on next electron (mc+md if 1stel else zero)
    integer,   intent(in) :: mq
    !> buffer offsets and lengths for each iteration of recursion
    !> 1:offset tail, 2:length tail, 3:offset head, 4:length head
    integer :: r(4, ma + mb + geo + mq)
    !-----------------------------------------------
    integer :: p, i, a, b, g, ab, m, np, npp
    integer :: posp, headp, numpp, nump, maxb, minb
    integer :: shift, total, maxtail
    logical :: wrapp, wrap, odd
    !-----------------------------------------------
    ! first pass, reverse of recursion order. Split in heads 
    ! and tails, and determine shifts. This is the complicated part
    ! Loop over descending momentum on P, increasing on A and B
    headp = 0
    do p = -1, ma + mb + geo + mq - 2
       np = (p+1)*(p+2)/2    !leading dim in rec p
       npp = (p+2)*(p+3)/2   !leading dim in rec p+1 (previous)
       a = min(ma, ma + mb + geo + mq - p)
       b = min(mb, max(0, mb + geo + mq - p))
       g = min(geo,max(0, geo + mq - p))
       i = ma + mb + geo + mq - 1 - p !index in r(:) of p+1
       ! initialize shifts to minus huge, lengths to zero
       r(:,i) = (/-huge(1)/2,0,-huge(1)/2,0/)
       ! loop over anti-diagonals with same a+b
       posp  = 0
       wrap  = .false.
       wrapp = .false.
       do ab=1, min(i+1,ma+mb+geo+1)
          ! if reached end of head, wrap
          if (.not.wrap .and. posp >= headp) then
             posp  = 0
             wrap  = .true.
             wrapp = .true.
          end if
          ! determine length of this diagonal in p
          minb = max(0, min(ab-a, (ab-a+b-g+1)/2))
          maxb = max(0, min(ab,   (ab-a+b+g)/2))
          posp = posp + sum((/((ab-m+1)*(ab-m+2)/2 &
                             * (m+1)*(m+2)/2, m = minb, maxb)/))
          ! determine length of previous diagonal
          odd = (p<mq .and. mod(a+b+g-ab+2,2)==1)
          if (ab > a .and. (ab <= a+b-g+1 .or. odd)) minb = minb-1
          if (.not.(ab <= a-b-g+1 .or. odd)) maxb = maxb-1
          numpp = sum((/((ab-1-m+1)*(ab-1-m+2)/2 &
                       * (m+1)*(m+2)/2, m = minb, maxb)/))
          ! if prev with alignment will exceed this, wrap it
          if (.not.wrapp .and. max(r(3,i) + npp * r(4,i), &
                   np * posp) + npp * numpp > np * headp) &
             wrapp = .true.
          ! update alignment
          if (ab/=1 .and. ab /= ma + mb + geo + 1) then
             if (.not.wrapp .and. .not.wrap) &
                r(3,i) = max(r(3,i), np * posp - npp * r(4,i))
             if (wrapp .and. wrap) &
                r(1,i) = max(r(1,i), np * posp - npp * r(2,i))
          end if
          ! increment size of head or tail
          if (ab/=1) then
             if (.not.wrapp) r(4,i) = r(4,i) + numpp
             if (     wrapp) r(2,i) = r(2,i) + numpp
          end if
       end do
       ! after all-tail comes another all-tail. Make that all-head instead
       if (headp == 0) r(:,i) = (/r(3:4,i),r(1:2,i)/)
       headp = r(4,i)
    end do
!     ! print shifts, tails and heads
!     do i=1, ma + mb + geo + mq
!        print *, r(:,i)
!     end do
!     print *, '-----------------------------------------------------'
    ! second pass, reverse. Determine total size
    total = 0
    shift = 4
    do p=0, ma + mb + geo + mq - 1
       np = (p+1)*(p+2)/2
       a = min(ma, ma + mb + geo + mq - p)
       b = min(mb, max(0, mb + geo + mq - p))
       g = min(geo,max(0, geo + mq - p))
       i = ma + mb + geo + mq - p
       ! if no tail, only head
       if (r(2,i) == 0) then
          shift = max((p+2)*(p+3)*(p+4)/6, shift + r(3,i))
          total = max(total, shift + np * r(4,i))
       else
          ! if lowest block is ab=10/01, push above pp's ab=00
          m = merge(p+1, p, r(2,i) == 0 .or. r(4,i) == 0 &
                       .or.(r(4,i) <= 6 .and. b+g>a))
          shift = max((m+1)*(m+2)*(m+3)/6, shift + r(1,i))
          total = max(total, shift + np * r(2,i) + (np-p-1) * r(4,i+1))
       end if
!        print '(4(a,i2),a,i6)', 'a =', a, '  b =', b, &
!                              '  g =', g, '  p =', p, &
!                              '  total =', total
    end do
!     print *, '-----------------------------------------------------'
    ! third pass, forward. Layout heads and tails
    maxtail = total
    do p = ma + mb + geo + mq - 1, 0, -1
       np = (p+1)*(p+2)/2
       i = ma + mb + geo + mq - p
       if (r(2,i) == 0) then
          maxtail = total - np * r(4,i) - r(3,i)
          r(3,i) = total - np * r(4,i)
          r(1,i) = r(3,i)
       else
          r(3,i) = total - np * r(4,i)
          maxtail = min(maxtail, min(total - (np-p-1) &
                      * r(4,i+1), r(3,i)) - np * r(2,i))
          maxtail = maxtail - r(1,i)
          r(1,i) = maxtail + r(1,i)
       end if
!        print '(a,i2,4(a,i7))', 'p =', p, &
!                             '   it =', r(1,i), '   nt =', r(2,i), &
!                             '   ih =', r(3,i), '   nh =', r(4,i)
    end do
!     print *
!     print *
  end function




  !> Angular contraction from Hermite- to geometry-differentiated
  !> spherical (real solid-harmonic) bases
  subroutine angular_hab_to_lagablb(ma, Ylma, minga, &
                                    mb, Ylmb, mingb, &
                                    maxga, hab, gab, ldim, &
                                    habint, lagablb, scratch)
    !> angular momentum on first center A
    integer, intent(in)  :: ma
    !> Cartesian to spherical contraction coefficients for A
    real(8), intent(in)  :: Ylma(2*ma+1, (ma+1)*(ma+2)/2)
    !> lowest order of geometry diff. for A
    integer, intent(in)  :: minga
    !> angular momentum on second center B
    integer, intent(in)  :: mb
    !> Cartesian to spherical contraction coefficients for B
    real(8), intent(in)  :: Ylmb(2*mb+1, (mb+1)*(mb+2)/2)
    !> lowest order of geometry diff. for B
    integer, intent(in)  :: mingb
    !> highest order of geometry diff. for A, maxb = minb+maxa-mina
    integer, intent(in)  :: maxga
    !> total number of Hermite components AB
    integer, intent(in)  :: hab
    !> total number of geo.diff. components
    integer, intent(in)  :: gab
    !> leading dimension
    integer, intent(in)  :: ldim
    !> input Hermite integrals to contract
    real(8), intent(in)  :: habint(ldim, hab)
    !> output geo-diff'ed solid-harmonic integrals
    real(8), intent(out) :: lagablb(ldim, 2*ma+1, gab, 2*mb+1)
    !> scratch buffer for half-contracted integrals
    real(8) :: scratch(ldim, 2*ma+1, (maxga+1)*(maxga+2)/2)
    !------------------------------------------------------
    integer :: ga, gb, ay, az, by, bz, ihab
    integer :: g, gy, gz, i, y, z, n, l, igab
! integer :: na, nb
! i = 0
! do ga = maxga, minga, -1
! gb = mingb + (maxga-ga)
! na = (ma+ga+1) * (ma+ga+2) / 2
! nb = (mb+gb+1) * (mb+gb+2) / 2
! call print_tensor((/ldim,na,nb/), 'iAB', 'icc', &
!                   habint(:,i+1:i+na*nb), &
!                   title = '###### before angular')
! i = i + na*nb
! end do
    ! loop over hermite blocks AmaBmb
lagablb(:,:,:,:) = 0
    ihab = 0
    igab = 0
    do ga = maxga, minga, -1
       gb = mingb + (maxga-ga)
       n  = (ga+1)*(ga+2)/2
       ! loop over Hermite indices by,bz in block
       do bz=0,mb+gb
          do by=0,mb+gb-bz
             ! for this by,bz, contract Hermites ay,az to geo-diff'ed
             ! solid-harmonic basis la,ga in scratch(:,l,g)
scratch(:,:,:n) = 0
             do az=0,ma+ga
                do ay=0,ma+ga-az
                   ! pre-increment AB block index
                   ihab = ihab+1
                   ! loop over geodiff components gy,gz and undiffed Hermites
                   ! on A, y,z to which Hermite ay,az will contribute
                   gz = max(0,az-ma)
                   g = (2*ga+3-gz)*gz/2
                   z = az-gz
                   i = (2*ma+3-(z+1))*(z+1)/2
                   do gz=gz, min(az,ga)
                      gy = max(0,ga-gz-(ma+ga-ay-az)) !=ax
                      g = g+gy
                      y = ay-gy
                      i = i - (ma-z-y) !=x
                      do gy=gy, min(ga-gz,ay)
                         ! pre-increment geodiff index
                         g = g+1
                         ! add contribution to diffed-angular
                         do l=1,2*ma+1
                            scratch(:,l,g) = scratch(:,l,g) &
                                           + Ylma(l,i) * habint(:,ihab)
                         end do
                         ! post-decrement undiffed Hermite indices
                         i = i-1
                         y = y-1
                      end do
                      g = g + (ga-gz-gy+1) !=gx
                      i = i - (y+1)
                      z = z-1
                   end do
                end do
             end do
             ! this by,bz contracted for center A in scratch, now add
             ! all contributions to
             gz = max(0,bz-mb)
             g = (2*gb+3-gz)*gz/2
             z = bz-gz
             i = (2*mb+3-(z+1))*(z+1)/2
             do gz=gz, min(bz,gb)
                gy = max(0,gb-gz-(mb+gb-by-bz)) !=bx
                g = g+gy
                y = by-gy
                i = i - (mb-z-y) !=x
                do gy=gy, min(gb-gz,by)
                   ! pre-increment geodiff index
                   g = g+1
                   ! add contribution to diffed-angular
                   do l=1,2*mb+1
                      lagablb(:,:,igab+(g-1)*n+1:igab+g*n,l) &
                               = lagablb(:,:,igab+(g-1)*n+1:igab+g*n,l) &
                               + Ylmb(l,i) * scratch(:,:,:n)
                   end do
                   ! post-decrement undiffed Hermite indices
                   i = i-1
                   y = y-1
                end do
                g = g + (gb-gz-gy+1) !=gx
                i = i - (y+1)
                z = z-1
             end do
          end do
       end do
       igab = igab + n * (gb+1)*(gb+2)/2
    end do
! i = 0
! do ga = maxga, minga, -1
! gb = mingb + (maxga-ga)
! na = (ga+1) * (ga+2) / 2
! nb = (gb+1) * (gb+2) / 2
! call print_tensor((/ldim,2*ma+1,na,nb,2*mb+1/), 'iAABB', 'isccs', &
!                   lagablb(:,:,i+1:i+na*nb,:), &
!                   title = '###### after angular')
! i = i + na*nb
! end do

!             ! for this g, number of l's this xyz contributes to
!             n = (amom + 2 - az - mod(ay,2))/2
!             ! index of first contribution l
!             i = mod(amom-az,2) - merge(0, 2*n, mod(ay,2)==0) + 1
!             do l = i, i+2*(n-1), 2
  end subroutine




  !> Contract from primitive Hermite basis hABeP to contracted
  !> geometry-differentiated radial basis rABhAB (notice the transpose)
  subroutine radial_habep_to_rabhab(ldim, A, B, geo, hab, &
                                    habep, rabhab, scratch)
    !> common leading dimension of arays habep and rabhab
    integer,    intent(in)  :: ldim
    !> the two CGTO shells
    type(cgto), intent(in)  :: A, B
    !> order of geometry-differentiation
    integer,    intent(in)  :: geo
    !> second, third dimension of arrays. No. Hermite comp.
    integer,    intent(in)  :: hab
    !> integral array with ep = A%exp+B%exp primitive basis
    real(8),    intent(in)  :: habep(ldim, hab, size(A%exp), &
                                                size(B%exp))
    !> output integral array with rab in contracted radial basis
    real(8),    intent(out) :: rabhab(ldim, size(A%ctr,1), &
                                            size(B%ctr,1), hab)
    !> scratch needed during contraction
    real(8)                 :: scratch(ldim, hab, size(A%exp)+1)
    !-----------------------------------------------------------
    integer :: i, n, ea, eb, ra, rb, ga, gb
    real(8) :: fac
! integer :: na, nb
! i = 0
! do gb=0,geo
! na = (A%mom + geo-gb+1) * (A%mom + geo-gb+2) / 2
! nb = (B%mom + gb+1) * (B%mom + gb+2) / 2
! call print_tensor((/ldim,na,nb/), 'IAB', 'icc', &
!                   habep(:,i+1:i+na*nb,1,1), &
!                   title = '###### before radial ')
! i = i + na*nb
! end do
    do rb=1, size(B%ctr,1)
       ! contract this (rb) radial on B, at the same time as
       ! scaling the (geo+1) Hermite blocks (hab) by exponent powers
       ! (2ea)^ga (2eb)^gb, which come from geometry-differentiation
       do eb=1, size(B%exp)
          do ea=1, size(A%exp)
             i = 0 !offset into hAB
             do gb=0,geo
                ga = geo-gb
                ! number of Hermites
                n = (A%mom+ga+1)*(A%mom+ga+2)/2 * (B%mom+gb+1)*(B%mom+gb+2)/2
                ! contraction factor including geometry-differentiation
                fac = B%ctr(rb,eb) * (2*A%exp(ea))**ga * (2*B%exp(eb))**gb
                if (eb==1) then
                   scratch(:,i+1:i+n,ea) = fac * habep(:,i+1:i+n,ea,eb)
                else
                   scratch(:,i+1:i+n,ea) = scratch(:,i+1:i+n,ea) &
                                         + fac * habep(:,i+1:i+n,ea,eb)
                end if
                i = i + n
             end do
          end do
       end do
       ! for this rb, contract all ra
       do ra=1, size(A%ctr,1)
          ! contract one radial on A
          do ea=1, size(A%exp)
             if (ea==1) then !first scale
                scratch(:,:,size(A%exp)+1) = A%ctr(ra,ea) * scratch(:,:,ea)
             else !then axpy
                scratch(:,:,size(A%exp)+1) = scratch(:,:,size(A%exp)+1) &
                                           + A%ctr(ra,ea) * scratch(:,:,ea)
             end if
          end do
          ! transpose and store. FIXME: overlap with scaling/normalization
          rabhab(:,ra,rb,:) = scratch(:,:,size(A%exp)+1)
       end do
    end do
! i = 0
! do gb=0,geo
! na = (A%mom + geo-gb+1) * (A%mom + geo-gb+2) / 2
! nb = (B%mom + gb+1) * (B%mom + gb+2) / 2
! call print_tensor((/ldim,na,nb/), 'IAB','icc', &
!                   rABhAB(:,1,1,i+1:i+na*nb), &
!                   title = '###### after  radial ')
! i = i + na*nb
! end do
  end subroutine




  subroutine reorder_rablagablb_to_gabralarblb(ra, la, rb, lb, &
                           gab, ldim, rablagablb, gabralarblb)
    !> radial and angular dimensions of integral arrays
    integer, intent(in)  :: ra, la, rb, lb
    !> number of geometry-differentiated components
    integer, intent(in)  :: gab
    !> leading dimension, to become trailing
    integer, intent(in)  :: ldim
    !> input array of integrals to be transposed
    real(8), intent(in)  :: rablagablb(ldim, ra, rb, la, gab, lb)
    !> resulting reordered array of integrals
    real(8), intent(out) :: gabralarblb(gab, ra, la, rb, lb, ldim)
    !--------------------------------------------------------
    integer :: g, l, t
! call print_tensor((/ldim,la,gab,lb/), 'IAGB','isis', &
!                   rablagablb(:,1,1,:,:,:),           &
!                   title = '###### before reorder')
    do g=1,gab
       do l=1,la
          do t=1,ldim
             ! ajt FIXME post-scaling of la and lb to hide latency
             gabralarblb(g,:,l,:,:,t) = rablagablb(t,:,:,l,g,:)
          end do
       end do
    end do
! call print_tensor((/gab,la,lb,ldim/), 'GABI','issi', &
!                   gabralarblb(:,1,:,1,:,:),          &
!                   title = '###### after reorder')
  end subroutine




  subroutine reorder_rcdlcgcdld_to_gcdrclcrdld(rc, lc, rd, ld, &
                     gcd, ldim, mdim, rcdlcgcdld, gcdrclcrdld)
    !> radial and angular dimensions of integral arrays
    integer, intent(in)  :: rc, lc, rd, ld
    !> number of geometry-differentiated components
    integer, intent(in)  :: gcd
    !> leading and middle dimensions
    integer, intent(in)  :: ldim, mdim
    !> input array of integrals to be transposed
    real(8), intent(in)  :: rcdlcgcdld(ldim, mdim, &
                                       rc, rd, lc, gcd, ld)
    !> resulting reordered array of integrals
    real(8), intent(out) :: gcdrclcrdld(ldim, gcd, mdim, &
                                        rc, lc, rd, ld)
    !--------------------------------------------------------
    integer :: g, l, m
! call print_tensor((/ldim,mdim,lc,gcd,ld/), 'IICGD','iisis', &
!                   rcdlcgcdld(:,:,1,1,:,:,:),                &
!                   title = '###### before reorder')
    do g=1,gcd
       do l=1,lc
          ! ajt FIXME possible post-scaling of lc and ld to hide latency
          gcdrclcrdld(:,g,:,:,l,:,:) = rcdlcgcdld(:,:,:,:,l,g,:)
       end do
    end do
! call print_tensor((/ldim,gcd,mdim,lc,ld/), 'IGICD','iiiss', &
!                   gcdrclcrdld(:,:,:,1,:,1,:),               &
!                   title = '###### after reorder')
  end subroutine




  ! Coulomb integrals over the "P-Q" center Hermite basis
  subroutine coulomb_PmQcdab(mom, A, B, C, D, PmQcdab)
    integer,    intent(in)  :: mom
    type(cgto), intent(in)  :: A, B, C, D
    real(8),    intent(out) :: PmQcdab((mom+1)*(mom+2)*(mom+3)/6, &
                                       size(C%exp),size(D%exp),   &
                                       size(A%exp),size(B%exp))
    integer :: ea, eb, ec, ed
    real(8) :: pq, PmQ(3)
    do eb=1, size(B%exp)
       do ea=1, size(A%exp)
          do ed=1, size(D%exp)
             do ec=1, size(C%exp)
                pq  = (A%exp(ea)+B%exp(eb)) * (C%exp(ec)+D%exp(ed)) &
                    / (A%exp(ea)+B%exp(eb)  +  C%exp(ec)+D%exp(ed))
                PmQ = A%exp(ea) / (A%exp(ea) + B%exp(eb)) * A%cent &
                    + B%exp(eb) / (A%exp(ea) + B%exp(eb)) * B%cent &
                    - C%exp(ec) / (C%exp(ec) + D%exp(ed)) * C%cent &
                    - D%exp(ed) / (C%exp(ec) + D%exp(ed)) * D%cent
                call hermite_coulomb(mom, pq, PmQ, PmQcdab(:,ec,ed,ea,eb))
             end do
          end do
       end do
    end do
  end subroutine




  !> move all momentum from "difference center" R to 1st electron's
  !> product center P by scaling with (q/(p+q))**m
  subroutine elxfer_Rqp_to_qPp(mom, A, B, C, D, Rqp, qPp)
    integer,    intent(in)  :: mom
    type(cgto), intent(in)  :: A, B, C, D
    real(8),    intent(in)  :: Rqp((mom+1)*(mom+2)*(mom+3)/6, &
                                   size(C%exp), size(D%exp),  &
                                   size(A%exp), size(B%exp))
    real(8),    intent(out) :: qPp(size(C%exp), size(D%exp),  &
                                   (mom+1)*(mom+2)*(mom+3)/6, &
                                   size(A%exp), size(B%exp))
    integer :: i, n, m, ea, eb, ec, ed
    real(8) :: p, q, inv, fac
    do eb=1, size(B%exp)
       do ea=1, size(A%exp)
          p = A%exp(ea) + B%exp(eb)
          do ed=1, size(D%exp)
             do ec=1, size(C%exp)
                q = C%exp(ec) + D%exp(ed)
                inv = 1 / (p+q)
                fac = sqrt(pi * inv)**3 &
                    * exp(-A%exp(ea) * B%exp(eb) / p &
                         * sum((A%cent - B%cent)**2)  &
                         - C%exp(ec) * D%exp(ed) / q &
                         * sum((C%cent - D%cent)**2))
                i = 1
                n = 1
                do m=0, mom
                   qPp(ec,ed,i:i+n-1,ea,eb) = fac * Rqp(i:i+n-1,ec,ed,ea,eb)
                   fac = fac * q * inv
                   i = i + n
                   n = n + m + 2
                end do
             end do
          end do
       end do
    end do
  end subroutine




  !> move all momentum from "difference center" R to 2nd electron's
  !> product center Q by scaling with (-p/(p+q))**m (in-place)
  subroutine elxfer_Rqp_to_Qqp(mom, A, B, C, D, Rqp)
    integer,    intent(in)    :: mom
    type(cgto), intent(in)    :: A, B, C, D
    real(8),    intent(inout) :: Rqp((mom+1)*(mom+2)*(mom+3)/6, &
                                   size(C%exp), size(D%exp),  &
                                   size(A%exp), size(B%exp))
    integer :: i, n, m, ea, eb, ec, ed
    real(8) :: p, q, inv, fac
    do eb=1, size(B%exp)
       do ea=1, size(A%exp)
          p = A%exp(ea) + B%exp(eb)
          do ed=1, size(D%exp)
             do ec=1, size(C%exp)
                q = C%exp(ec) + D%exp(ed)
                inv = 1 / (p+q)
                fac = sqrt(pi * inv)**3 &
                    * exp(-A%exp(ea) * B%exp(eb) / p &
                         * sum((A%cent - B%cent)**2) &
                         - C%exp(ec) * D%exp(ed) / q &
                         * sum((C%cent - D%cent)**2))
                i = 1
                n = 1
                do m=0, mom
                   Rqp(i:i+n-1,ec,ed,ea,eb) = fac * Rqp(i:i+n-1,ec,ed,ea,eb)
                   fac = - fac * p * inv
                   i = i + n
                   n = n + m + 2
                end do
             end do
          end do
       end do
    end do
  end subroutine


end module
