! Copyright 2012 Andreas J. Thorvaldsen
! This file is made available under the terms of the
! GNU Lesser General Public License version 3.

!> @file Contains module eri_contractions

!> Electron-repulsion-integral contractions over sets of integrals 
!> (AB|CD) for one CGTO quadruple
module eri_contractions

  use basis_set,      only: cgto
  use matrix_backend, only: matrix

  implicit none

  public contract_with_density
  public ctr_arg
  public set_eri_contractions_xfac

  !> A 'contraction' is either of 'average' type, or of 'Fock' type,
  !> depending on whether pointer @param average is associated.
  !> A fock contraction adds two-electron contribution to *one*
  !> Fock matrix given in @param fock_or_dens, whose geometrical
  !> coordinate(s) is chosen by @param comp. @param ncor determines
  !> the dimensions of average (implicitly).
  type ctr_arg
     integer               :: geo, comp, ncor
     type(matrix), pointer :: dens, fock_or_dens
     real(8),      pointer :: average(:)
  end type
  
  private

! exchange factor, 1 for HF, 0 for non-hybrid DFT, between 0 and 1 for hybrids
  real(8), save :: xfac = 1.0d0
  
contains

  subroutine set_eri_contractions_xfac(f)
    real(8), intent(in) :: f
    xfac = f
  end subroutine

  !> 
  subroutine contract_with_density(A, B, C, D, gab, gcd, ng, eri, arg)
    type(cgto), intent(in) :: A, B, C, D
    integer,    intent(in) :: gab, gcd, ng
    real(8),    intent(in) :: eri(ng, size(A%ctr,1)*(2*A%mom+1), &
                                      size(B%ctr,1)*(2*B%mom+1), &
                                      size(C%ctr,1)*(2*C%mom+1), &
                                      size(D%ctr,1)*(2*D%mom+1))
    type(ctr_arg), intent(inout) :: arg(:)
    real(8) ctr(ng)
    integer iarg, i, j, k, l
    do iarg = 1, size(arg)
       ! average (any order)
       if (associated(arg(iarg)%average)) then
          call coul_exch_ave(arg(iarg)%dens%nrow, &
                             arg(iarg)%dens%elms, &
                             arg(iarg)%fock_or_dens%elms)
          call accumulate_ave(arg(iarg)%ncor, arg(iarg)%average)
       ! unperturbed Fock matrix
       else if (gab == 0 .and. gcd == 0) then
          call coul_exch_mat(1, arg(iarg)%dens%nrow, &
                                arg(iarg)%dens%elms, &
                                arg(iarg)%fock_or_dens%elms)
       ! first-order Fock matrix
       else if (gab + gcd == 1) then
          do i = 0, 1 !0: A or C, 1: B or D
             j = merge(merge(A%icent, B%icent, i==0), &
                       merge(C%icent, D%icent, i==0), gcd==0)
             k = arg(iarg)%comp - 3*(j-1)
             if (.not.(k>=1 .and. k<=3)) cycle !x, y or z
             call coul_exch_mat(k+3*i, arg(iarg)%dens%nrow, &
                                       arg(iarg)%dens%elms, &
                                       arg(iarg)%fock_or_dens%elms)
          end do
       ! second-order Fock matrix, both on same electron
       else if (max(gab,gcd) == 2 .and. min(gab,gcd) == 0) then
          do i = 0, 3 !0:AA/CC, 1:BA/DC, 2:AB/CD, 3:BB/DD
             j = merge(merge(B%icent, A%icent, btest(i,0)), &
                       merge(D%icent, C%icent, btest(i,0)), gab/=0)
             k = merge(merge(B%icent, A%icent, btest(i,1)), &
                       merge(D%icent, C%icent, btest(i,1)), gcd==0)
             l = arg(iarg)%comp - 3*(j-1) - arg(iarg)%ncor * 3*(k-1)
             if (l > 3 .and. l <= arg(iarg)%ncor) cycle
             if (l > 3) l = l - arg(iarg)%ncor + 3
             if (l > 6 .and. l <= arg(iarg)%ncor + 3) cycle
             if (l > 6) l = l - arg(iarg)%ncor + 3
             if (l <= 0 .or. l > 9) cycle
             ! if BA/DC or lower triangle of AA/BB/CC/DD, 'transpose' l
             if (i==1 .or. ((i==0 .or. i==3) .and. mod(l-1,3) > (l-1)/3)) &
                l = 1 + (l-1)/3 + 3 * mod(l-1,3)
             ! if AA, BB, CC or DD, go from 9 to 6 comps
             if (i==0 .or. i==3) &
                l = l - merge(merge(3,2,l>6),0,l>3)
             ! add offset
             l = l + merge(merge(15,6,i==3),0,i>=1)
             ! call
             call coul_exch_mat(l, arg(iarg)%dens%nrow, &
                                   arg(iarg)%dens%elms, &
                                   arg(iarg)%fock_or_dens%elms)
          end do
       ! second-order Fock matrix, one on each electron
       else if (gab == 1 .and. gcd == 1) then
          do i = 0, 7 !0:AC 1:BC 2:AD 3:BD 4:CA 5:CB 6:DA 7:DB
             j = merge(merge(B%icent, A%icent, btest(i,0)), &
                       merge(D%icent, C%icent, btest(i,1)), .not.btest(i,2))
             k = merge(merge(B%icent, A%icent, btest(i,0)), &
                       merge(D%icent, C%icent, btest(i,1)), btest(i,2))
             l = arg(iarg)%comp - 3*(j-1) - arg(iarg)%ncor * 3*(k-1)
             if (l > 3 .and. l <= arg(iarg)%ncor) cycle
             if (l > 3) l = l - arg(iarg)%ncor + 3
             if (l > 6 .and. l <= arg(iarg)%ncor + 3) cycle
             if (l > 6) l = l - arg(iarg)%ncor + 3
             if (l <= 0 .or. l > 9) cycle
             ! if CA/CB/DA/DB 'transpose' l
             if (btest(i,2)) l = 1 + (l-1)/3 + 3 * mod(l-1,3)
             ! leading dimension is 6, ie. AxAyAzBxByBz
             l = 1 + mod(l-1,3) + 6*((l-1)/3)
             ! add offset
             l = l + merge(3,0,btest(i,0)) + merge(18,0,btest(i,1))
             ! call
             call coul_exch_mat(l, arg(iarg)%dens%nrow, &
                                   arg(iarg)%dens%elms, &
                                   arg(iarg)%fock_or_dens%elms)
          end do
       else
          call quit('only averages and undiff, 1st and 2nd order Fock matrices implemented')
       end if
    end do

  contains


    subroutine coul_exch_ave(nbas, dens, densb)
      integer, intent(in) :: nbas
      real(8), intent(in) :: dens(nbas,nbas), densb(nbas,nbas)
      integer i, j, k, l, ii, jj, kk, ll
      ! contract with pair of density matrices
      ctr = 0
      l = D%ibas + 1
      do ll=1, size(eri,5)
         k = C%ibas + 1
         do kk=1, size(eri,4)
            j = B%ibas + 1
            do jj=1, size(eri,3)
               i = A%ibas + 1
               do ii=1, size(eri,2)
                  ctr(:) = ctr(:) +  (2*dens(i,j)*densb(k,l) &
                                 - xfac*dens(i,l)*densb(k,j)) * eri(:,ii,jj,kk,ll)
                  i = i + 2*A%mom + 1
                  if (i > A%ibas + size(eri,2)) i = i - size(eri,2) + 1
               end do
               j = j + 2*B%mom + 1
               if (j > B%ibas + size(eri,3)) j = j - size(eri,3) + 1
            end do
            k = k + 2*C%mom + 1
            if (k > C%ibas + size(eri,4)) k = k - size(eri,4) + 1
         end do
         l = l + 2*D%mom + 1
         if (l > D%ibas + size(eri,5)) l = l - size(eri,5) + 1
      end do
    end subroutine


    subroutine accumulate_ave(ncor, ave)
      integer, intent(in)    :: ncor
      real(8), intent(inout) :: ave(ncor**(gab+gcd))
      real(8) c33(3**(gab+gcd),2) !scratch for unpacking of indices to 3^m
      integer tab, iab, icd, na, nb, nc, nd, ga, gb, gc, gd, i2, i
      integer divab, divcd, offset
      ! the total first electron dimension of ctr
      tab = (gab+1)*(gab+2)*(gab+3)*(gab+4)*(gab+5)/120
      icd = 0
      do gd = 0, gcd
         gc = gcd-gd
         iab = 0
         ! factorial divisor for C and D
         divcd = product((/(i,i=1,gc)/)) * product((/(i,i=1,gd)/))
         do gb = 0, gab
            ga = gab-gb
            na = (ga+1)*(ga+2)/2
            nb = (gb+1)*(gb+2)/2
            nc = (gc+1)*(gc+2)/2
            nd = (gd+1)*(gd+2)/2
            ! factorial divisor for A and B
            divab = product((/(i,i=1,ga)/)) * product((/(i,i=1,gb)/))
            ! determine which c33(:,i2) to start in, in order to
            ! end up in c33(:,1)
            i2 = 1 + mod(count((/ga>1,gb>1,gc>1,gd>1/)), 2)
            ! copy-reshape-scale from ctr to c33
            call copy_scale_rows(iab, na*nb, tab, nc*nd, &
                                 divab * divcd,          &
                                 ctr(tab*icd+1 : tab*(icd+nc*nd)), &
                                 c33(:na*nb*nc*nd, i2))
            iab = iab + na*nb
            ! unpack indices of A, if needed
            if (ga > 1) then
               call unpack_multicart(1, ga, nb*nc*nd, &
                          c33(:   na*nb*nc*nd,   i2), &
                          c33(:3**ga*nb*nc*nd, 3-i2))
               na = 3**ga
               i2 = 3-i2
            end if
            ! unpack indices of B, if needed
            if (gb > 1) then
               call unpack_multicart(na, gb, nc*nd, &
                        c33(:na*   nb*nc*nd,   i2), &
                        c33(:na*3**gb*nc*nd, 3-i2))
               nb = 3**gb
               i2 = 3-i2
            end if
            ! unpack indices of C, if needed
            if (gc > 1) then
               call unpack_multicart(na*nb, gc, nd, &
                        c33(:na*nb*   nc*nd,   i2), &
                        c33(:na*nb*3**gc*nd, 3-i2))
               nc = 3**gc
               i2 = 3-i2
            end if
            ! unpack indices of D, if needed
            if (gd > 1) then
               call unpack_multicart(na*nb*nc, gd, 1, &
                                c33(:na*nb*nc*nd, 2), &
                                c33(:           , 1))
               nd = 3**gd
            end if
            ! compute offset into average
            offset = 0
            do i = 1, gd
               offset = 3*(D%icent-1) + ncor * offset
            end do
            do i = 1, gc
               offset = 3*(C%icent-1) + ncor * offset
            end do
            do i = 1, gb
               offset = 3*(B%icent-1) + ncor * offset
            end do
            do i = 1, ga
               offset = 3*(A%icent-1) + ncor * offset
            end do
            ! accumulate unpacked derivatives into ave
            call add_unpacked_to_ave(gab+gcd, c33(:,1), ncor, &
                                     offset, ave)
         end do
         ! increment offset
         icd = icd + (gc+1)*(gc+2)/2 * (gd+1)*(gd+2)/2
      end do
    end subroutine



    subroutine coul_exch_mat(comp, nbas, dens, fock)
      integer, intent(in)    :: comp, nbas
      real(8), intent(in)    :: dens(nbas,nbas)
      real(8), intent(inout) :: fock(nbas,nbas)
      integer i, j, k, l, ii, jj, kk, ll
      l = D%ibas + 1
      do ll=1, size(eri,5)
         k = C%ibas + 1
         do kk=1, size(eri,4)
            j = B%ibas + 1
            do jj=1, size(eri,3)
               i = A%ibas + 1
               do ii=1, size(eri,2)
                  ! Coulomb
                  fock(i,j) = fock(i,j) +    2*dens(k,l) * eri(comp,ii,jj,kk,ll)
                  ! Exchange
                  fock(i,l) = fock(i,l) - xfac*dens(k,j) * eri(comp,ii,jj,kk,ll)
                  i = i + 2*A%mom + 1
                  if (i > A%ibas + size(eri,2)) i = i - size(eri,2) + 1
               end do
               j = j + 2*B%mom + 1
               if (j > B%ibas + size(eri,3)) j = j - size(eri,3) + 1
            end do
            k = k + 2*C%mom + 1
            if (k > C%ibas + size(eri,4)) k = k - size(eri,4) + 1
         end do
         l = l + 2*D%mom + 1
         if (l > D%ibas + size(eri,5)) l = l - size(eri,5) + 1
      end do
    end subroutine

  end subroutine



  !> simple inline to copy some rows from matrix art over to arn,
  !> while dividing by integer div
  subroutine copy_scale_rows(i, n, t, n2, div, art, arn)
    integer, intent(in)    :: i, n, t, n2, div
    real(8), intent(in)    :: art(t,n2)
    real(8), intent(inout) :: arn(n,n2)
    arn(:,:) = 1d0/div * art(i+1:i+n,:)
  end subroutine


  !> unpack/unravel/stretch triangular indices for x^mx*y^my*z^mz (mx+my+mz=@param m)
  !> to full 3**@param m tensor. @param l and @param t are leading and trailing dimensions.
  subroutine unpack_multicart(l, m, t, p, r)
    integer, intent(in)  :: l, m, t !leading, momentum, trailing
    real(8), intent(in)  :: p(l,(m+1)*(m+2)/2,t)
    real(8), intent(out) :: r(l,3**m,t)
    integer ii, iy, iz, mm, my, mz, s
    iy=0; iz=0; mm=1; my=0; mz=0
    do ii = 1, 3**m
       r(:,ii,:) = p(:,mm,:)
       do s=0, m-1
          ! x -> y transition
          if (.not.btest(iy,s) .and. .not.btest(iz,s)) then
             iy = ibset(iy,s)
             mm = mm+1
             my = my+1
             exit
          ! y -> z transition
          else if (btest(iy,s)) then
             iy = ibclr(iy,s)
             iz = ibset(iz,s)
             mm = mm + (m-mz)
             my = my-1
             mz = mz+1
             exit
          ! z -> x transition (followed by 'carry')
          else
             iz = ibclr(iz,s)
             mm = mm - (m-mz) - 2
             mz = mz-1
          end if
       end do
    end do
  end subroutine


  !> add ctr(3,3,...) to ave(nc,nc,...) at effective offset off
  subroutine add_unpacked_to_ave(geo, ctr, nc, off, ave)
    integer, intent(in)    :: geo, nc, off
    real(8), intent(in)    :: ctr(3**geo)
    real(8), intent(inout) :: ave(nc**geo)
    integer i, y, z, j, s
    i=off+1; y=0; z=0; 
    do j = 1, 3**geo
       ave(i) = ave(i) + ctr(j)
       ! increment i
       do s = 0, geo-1
          ! x -> y transition
          if (.not.btest(y,s) .and. .not.btest(z,s)) then
             y = ibset(y,s)
             i = i + nc**s
             exit
          ! y -> z transition
          else if (btest(y,s)) then
             y = ibclr(y,s)
             z = ibset(z,s)
             i = i + nc**s
             exit
          ! z -> x transition (followed by 'carry')
          else
             z = ibclr(z,s)
             i = i - 2*nc**s
          end if
       end do
    end do
  end subroutine


end module
