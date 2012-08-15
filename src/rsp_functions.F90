! Copyright 2012      Magnus Ringholm
!           2012      Dan Jonsson
!           2009-2011 Radovan Bast
!           2009-2011 Andreas J. Thorvaldsen
! This file is made available under the terms of the
! GNU Lesser General Public License version 3.

!> @file Contains module rsp_functions

!> Front-end to OpenRsp. This module organizes, computes and prints
!> response function tensors.
module rsp_functions

  use matrix_defop
  use rsp_contribs
  use rsp_equations
  use interface_molecule
  use interface_rsp_solver
  use interface_xc

  implicit none

  public prop_test_gradient
  public prop_test_hessian
  public prop_test_cubicff
  public prop_test_quarticff
  public prop_test_diphes

  private

  !unit number for IO in this file
  !radovan: a unit nr this high may not work on all compilers and/or processors
  !         we should be careful going beyond 100
  integer, parameter :: iounit = 345645

contains


  !> Calculate the gradient: dE/dR  =  dh_nuc/dR  +  Tr dH/dR D
  !>            +  1/2 Tr dG/dR(D) D  +  dExc[D]/dR  -  Tr dS/dR DFD
  subroutine prop_test_gradient(ng, S, D, F)
    integer,       intent(in) :: ng
    type(matrix),  intent(in) :: S, D, F
    complex(8)      gra(ng), tmp(ng)
    type(matrix)    DFD
    type(rsp_field) geo(1)

    real(8) :: total_force(3)
    real(8) :: total_torque(3)
    integer :: iatom, ixyz, k

    real(8)             :: charges(ng/3), coords(ng)
    real(8)             :: masses(ng/3), totmass, origin(3)
    integer             :: isotopes(ng/3)

    geo(1) = rsp_field('GEO ', (0d0,0d0), 1, ng)
    ! first nuclear contribution
    tmp = 0.0d0
    call rsp_nucpot(geo, tmp)
    gra = tmp
    call print_tensor(shape(tmp), tmp, 'nucpot')
    ! overlap contribution
    DFD = D*F*D
    call rsp_ovlave(1, (/'GEO '/), (/1/), shape(tmp), DFD, &
                    tmp, (/(0d0,0d0)/), D)
    DFD = 0 !free
    gra = gra + tmp
    call print_tensor(shape(tmp), tmp, 'ovlave')
    ! 1-electron contribution
    call rsp_oneave(1, (/'GEO '/), (/1/), shape(tmp), D, tmp)
    gra = gra + tmp
    call print_tensor(shape(tmp), tmp, 'oneave')
    ! 2-electron contribution
    call rsp_twoave(1, (/'GEO '/), (/1/), shape(tmp), D, D, tmp)
    gra = gra + tmp/2
    call print_tensor(shape(tmp), tmp/2, 'twoave')
    ! Kohn-Sham exchange correlation average
    call rsp_xcave(pert='g', res=tmp, D=D)
    gra = gra + tmp
    call print_tensor(shape(tmp), tmp, 'xcave')

    ! print to screen
    call print_tensor(shape(gra), gra, 'gradient = Eg')

    total_force = 0.0d0
    k = 0
    do iatom = 1, ng/3
       do ixyz = 1, 3
          k = k + 1
          total_force(ixyz) = total_force(ixyz) - dreal(gra(k))
       end do
    end do
    write(*, '(a, 3e12.3)') 'total force = ', total_force

    do iatom = 1, ng/3
       charges(iatom)  = get_nuc_charge(iatom)
       isotopes(iatom) = get_nuc_isotope(iatom)
       coords((iatom-1)*3 + 1) = get_nuc_xyz(1, iatom)
       coords((iatom-1)*3 + 2) = get_nuc_xyz(2, iatom)
       coords((iatom-1)*3 + 3) = get_nuc_xyz(3, iatom)
    end do
    call vibmas(masses, totmass, isotopes, nint(charges), ng/3, coords, origin, 1)

    total_torque = 0.0d0
    k = 0
    do iatom = 1, ng/3
       total_torque(1) = total_torque(1) + (get_nuc_xyz(2, iatom) - origin(2))*dreal(gra(k + 3)) &
                                         - (get_nuc_xyz(3, iatom) - origin(3))*dreal(gra(k + 2))
       total_torque(2) = total_torque(2) + (get_nuc_xyz(3, iatom) - origin(3))*dreal(gra(k + 1)) &
                                         - (get_nuc_xyz(1, iatom) - origin(1))*dreal(gra(k + 3))
       total_torque(3) = total_torque(3) + (get_nuc_xyz(1, iatom) - origin(1))*dreal(gra(k + 2)) &
                                         - (get_nuc_xyz(2, iatom) - origin(2))*dreal(gra(k + 1))
       k = k + 3
    end do
    write(*, '(a, 3e12.3)') 'total torque = ', total_torque

    ! print to file
    open(unit=iounit, file='gradient', status='replace', action='write')
    call print_tensor(shape(gra), gra, unit=iounit)
    close(iounit)

  end subroutine

  subroutine get_fo_geo_perturbed_matrices(ng, S, D, F, Sg, Dg, Fg)

!   ----------------------------------------------------------------------------
    integer,       intent(in)    :: ng
    type(matrix),  intent(in)    :: S, D, F
    type(matrix),  intent(inout) :: Sg(ng), Dg(ng), Fg(ng)
!   ----------------------------------------------------------------------------
    type(matrix)                 :: X(1), FDSg(1)
    character(4)                 :: nof(0)
    integer                      :: noc(0)
    integer                      :: i
!   ----------------------------------------------------------------------------

!   calculate perturbed integrals
    call rsp_ovlint(S%nrow, 1, (/'GEO '/), (/1/), shape(Sg), Sg)
    call rsp_oneint(S%nrow, 1, (/'GEO '/), (/1/), shape(Fg), Fg)
    call rsp_twoint(S%nrow, 1, (/'GEO '/), (/1/), shape(Fg), D, Fg)
    call rsp_xcint(D=(/D/), Fg=Fg)

!   solve equations
    do i = 1, size(Dg)

       Dg(i) = -D*Sg(i)*D

       FDSg(1) = F*D*Sg(i) - Sg(i)*D*F
       FDSg(1) = FDSg(1)*D*S - S*D*FDSg(1) + Fg(i)

       call rsp_twoint(S%nrow, 0, nof, noc, noc, Dg(i), FDSg(1:1))
       call rsp_xcint(D=(/D, Dg(i)/), F=FDSg(1))

       FDSg(1) = FDSg(1)*D*S - S*D*FDSg(1)

       X(1) = mat_alloc_like(D)
       call rsp_mosolver_exec(FDSg(1), (/0d0/), X)
       FDSg(1) = 0

#ifdef PRG_DIRAC
       Dg(i) = Dg(i) + X(1)
#else
       X(1) = -2.0d0*X(1)
       Dg(i) = Dg(i) + X(1)*S*D - D*S*X(1)
#endif
       X(1) = 0

       call rsp_twoint(S%nrow, 0, nof, noc, noc, Dg(i), Fg(i:i))
       call rsp_xcint(D=(/D, Dg(i)/), F=Fg(i))

    end do

  end subroutine

  subroutine contract_hessian(ng, S, D, F, Sg, Dg, Fg)

!   ----------------------------------------------------------------------------
    integer,       intent(in) :: ng
    type(matrix),  intent(in) :: S, D, F
    type(matrix),  intent(in) :: Sg(ng), Dg(ng), Fg(ng)
!   ----------------------------------------------------------------------------
    type(rsp_field)           :: geo3(3)
    type(matrix)              :: DFD, DFDg
    complex(8)                :: hes(ng, ng), temp(ng, ng), xc_hes(ng, ng)
    integer                   :: i, j
    real(8)                   :: t
!   ----------------------------------------------------------------------------

    geo3(:3) = rsp_field('GEO ', (0d0,0d0), 1, ng)
    ! nuclear repulsion
    temp = 0
    call rsp_nucpot(geo3(:2), temp)
    hes = temp; call print_tensor(shape(temp), temp, 'hnuc_ab')
    ! overlap contribution
    DFD = D*F*D
    call rsp_ovlave(2, (/'GEO ','GEO '/), (/1,1/), shape(temp), DFD, &
                    temp); DFD=0
    hes = hes + temp; call print_tensor(shape(temp), temp, '-SabDFD')
    !
    do i = 1, size(Dg)
       DFDg = Dg(i)*F*D + D*Fg(i)*D + D*F*Dg(i)
       call rsp_ovlave(1, (/'GEO '/), (/1/), shape(temp(:,i)), DFDg, &
                       temp(:,i)); DFDg=0
    end do
    hes = hes + temp; call print_tensor(shape(temp), temp, '-SaDFDb')
    ! 1-electron contribution
    call rsp_oneave(2, (/'GEO ','GEO '/), (/1,1/), shape(temp), D, temp)
    hes = hes + temp; call print_tensor(shape(temp), temp, 'HabD')
    !
    do i = 1, size(Dg)
       call rsp_oneave(1, (/'GEO '/), (/1/), shape(temp(:,i)), Dg(i), temp(:,i))
    end do
    hes = hes + temp; call print_tensor(shape(temp), temp, 'HaDb')

#ifndef PRG_DIRAC
    ! 2-electron contribution
    call rsp_twoave(2, (/'GEO ','GEO '/), (/1,1/), shape(temp), D, D, temp)
    hes = hes + temp/2; call print_tensor(shape(temp), temp/2, 'Gab(D)D/2')
#endif

    do i = 1, size(Dg)
       call rsp_twoave(1, (/'GEO '/), (/1/), shape(temp(:,i)), D, Dg(i), temp(:,i))
    end do
    hes = hes + temp; call print_tensor(shape(temp), temp, 'Ga(D)Db')

    xc_hes = 0.0d0
    ! Exchange/correlation contribution
    call rsp_xcave(pert='gg', res=xc_hes, D=D, Dg=Dg)
    call print_tensor(shape(xc_hes), xc_hes, 'xc_hes')
    hes = hes + xc_hes

!   cheat: only symmetry unique elements are correctly calculated with xc
!          here copy to symmetry dependent to get correct, symmetric hessian
!          i (radovan) will later move this one lever down and hide
!          it in some xc backend
    do i = 1, size(Dg)
       do j = 1, i
          hes(j, i) = hes(i, j)
       end do
    end do

!   check sum rules
    do i = 1, size(Dg)
       t = 0.0d0
       do j = 1, size(Dg)
          t = t + real(hes(i, j))
       end do
       print *, 'hessian column sum', i, t
    end do

    ! print
    call print_tensor(shape(hes), hes, 'hessian = Egg')
    ! print formatted to file hessian
    open (unit=iounit, file='hessian', status='replace', action='write')
    call print_tensor(shape(hes), hes, unit=iounit)
    close (iounit)
!   !------------------------------------
  end subroutine

  subroutine contract_cubicff(ng, S, D, F, Sg, Dg, Fg)

!   ----------------------------------------------------------------------------
    integer,       intent(in) :: ng
    type(matrix),  intent(in) :: S, D, F
    type(matrix),  intent(in) :: Sg(ng), Dg(ng), Fg(ng)
!   ----------------------------------------------------------------------------
    type(rsp_field)           :: geo3(3)
    type(matrix)              :: DFD, DFDg, DFDgg, DSDgg, FDSgg
    type(matrix)              :: FgDS(ng), DgSD(ng)
    complex(8)                :: cub(ng, ng, ng), tmp(ng, ng, ng)
    complex(8)                :: xc_cub(ng, ng, ng)
    integer                   :: i, j, k
    character(4)              :: nof(0)
    integer                   :: noc(0)
!   ----------------------------------------------------------------------------


!   nuclear repulsion
!   =================

    geo3(:3) = rsp_field('GEO ', (0d0,0d0), 1, ng)
    tmp = 0.0d0
    call rsp_nucpot(geo3, tmp)
    cub = tmp


!   overlap contribution
!   ====================

    DFD = D*F*D
    tmp = 0.0d0
    call rsp_ovlave(3, (/'GEO ','GEO ','GEO '/), (/1,1,1/), &
                    shape(tmp), DFD, tmp)
    cub = cub + tmp
    DFD = 0
    do i = 1, size(Dg)
       DFDg = Dg(i)*F*D + D*Fg(i)*D + D*F*Dg(i)
       tmp = 0.0d0
       call rsp_ovlave(2, (/'GEO ','GEO '/), (/1,1/), shape(tmp(i,:,:)), &
                       DFDg, tmp(i,:,:))
       do j = 1, size(Dg)
          do k = 1, size(Dg)
             cub(i, j, k) = cub(i, j, k) + tmp(i, j, k)
             cub(j, i, k) = cub(j, i, k) + tmp(i, j, k)
             cub(j, k, i) = cub(j, k, i) + tmp(i, j, k)
          end do
       end do
       DFDg = 0
       do j = 1, i
          DFDgg = Dg(j)*Fg(i)*D + Dg(j)*F*Dg(i) + D*Fg(j)*Dg(i) &
                + Dg(i)*Fg(j)*D + Dg(i)*F*Dg(j) + D*Fg(i)*Dg(j)
          tmp = 0.0d0
          call rsp_ovlave(1, (/'GEO '/), (/1/), shape(tmp(:,i,j)), &
                          DFDgg, tmp(:,i,j))
          do k = 1, size(Dg)
             cub(k, i, j) = cub(k, i, j) + tmp(k, i, j)
             if (i /= j) then
                cub(k, j, i) = cub(k, j, i) + tmp(k, i, j)
             end if
          end do
          DFDgg = 0
       end do
    end do


!   1-electron contribution
!   =======================

    tmp = 0.0d0
    call rsp_oneave(3, (/'GEO ','GEO ','GEO '/), (/1,1,1/), shape(tmp), D, tmp)
    cub = cub + tmp
    do i = 1, size(Dg)
       tmp = 0.0d0
       call rsp_oneave(2, (/'GEO ','GEO '/), (/1,1/), shape(tmp(i,:,:)), &
                       Dg(i), tmp(i,:,:))
       do j = 1, size(Dg)
          do k = 1, size(Dg)
             cub(i, j, k) = cub(i, j, k) + tmp(i, j, k)
             cub(j, i, k) = cub(j, i, k) + tmp(i, j, k)
             cub(j, k, i) = cub(j, k, i) + tmp(i, j, k)
          end do
       end do
    end do


!   2-electron contribution
!   =======================

    tmp = 0.0d0
    call rsp_twoave(3, (/'GEO ','GEO ','GEO '/), (/1,1,1/), shape(tmp), D, D, tmp)
    cub = cub + tmp/2
    do i = 1, size(Dg)
       tmp = 0.0d0
       call rsp_twoave(2, (/'GEO ','GEO '/), (/1,1/), shape(tmp(i,:,:)), &
                       D, Dg(i), tmp(i,:,:))
       do j = 1, size(Dg)
          do k = 1, size(Dg)
             cub(i, j, k) = cub(i, j, k) + tmp(i, j, k)
             cub(j, i, k) = cub(j, i, k) + tmp(i, j, k)
             cub(j, k, i) = cub(j, k, i) + tmp(i, j, k)
          end do
       end do
       do j = 1, i
          tmp = 0.0d0
          call rsp_twoave(1, (/'GEO '/), (/1/), shape(tmp(i,j,:)), &
                          Dg(i), Dg(j), tmp(i,j,:))
          do k = 1, size(Dg)
             cub(i, j, k) = cub(i, j, k) + tmp(i, j, k)
             cub(i, k, j) = cub(i, k, j) + tmp(i, j, k)
             cub(k, i, j) = cub(k, i, j) + tmp(i, j, k)
             if (i /= j) then
                cub(j, i, k) = cub(j, i, k) + tmp(i, j, k)
                cub(j, k, i) = cub(j, k, i) + tmp(i, j, k)
                cub(k, j, i) = cub(k, j, i) + tmp(i, j, k)
             end if
          end do
       end do
    end do


!   xc contribution
!   ===============

    xc_cub = 0.0d0
    call rsp_xcave(pert='ggg', res=xc_cub, D=D, Dg=Dg)
    call print_tensor(shape(xc_cub), xc_cub, 'xc_cub')
    cub = cub + xc_cub


!   idempotency multiplier contribution
!   ===================================

    do i = 1, size(Dg)
       FgDS(i) = Fg(i)*D*S + S*D*Fg(i) - Fg(i) - F*D*Sg(i) - Sg(i)*D*F
    end do
    do k = 1, size(Dg)
       do j = 1, size(Dg)
          DSDgg = Dg(k)*(Sg(j)*D + S*Dg(j)) + D*Sg(k)*Dg(j) &
                + Dg(j)*(Sg(k)*D + S*Dg(k)) + D*Sg(j)*Dg(k)
          do i = 1, size(Dg)
             tmp(i,j,k) = -tr(FgDS(i),DSDgg)
          end do; DSDgg=0
       end do
    end do; FgDS=0;
    cub = cub + tmp


!   self-consistency multiplier contribution
!   ========================================

    do i = 1, size(Dg)
       DgSD(i) = Dg(i)*S*D - D*S*Dg(i)
    end do
    do k = 1, size(Dg)
       do j = 1, size(Dg)
          FDSgg = (Fg(k)*Dg(j)+Fg(j)*Dg(k))*S - (Sg(k)*Dg(j)+Sg(j)*Dg(k))*F &
                + (Fg(k)*D+F*Dg(k))*Sg(j) - (Sg(k)*D+S*Dg(k))*Fg(j) &
                + (Fg(j)*D+F*Dg(j))*Sg(k) - (Sg(j)*D+S*Dg(j))*Fg(k)
          do i = 1, size(Dg)
             tmp(i,j,k) = -tr(DgSD(i),FDSgg)
          end do; FDSgg=0
       end do
    end do; DgSD=0;
    cub = cub + tmp


    ! print
    call print_tensor(shape(cub), cub, 'cubicff = Eggg')
    ! print formatted to file cubicff
    open (unit=iounit, file='cubicff', status='replace', action='write')
    write (iounit,*)
    do i = 1, size(cub,3)
       call print_tensor(shape(cub(:,:,i)), cub(:,:,i), unit=iounit)
    end do
    close (iounit)

  end subroutine

  subroutine prop_test_hessian(ng, S, D, F)

    integer,       intent(in) :: ng
    type(matrix),  intent(in) :: S, D, F

    type(matrix)              :: Sg(ng), Dg(ng), Fg(ng)

    call print_nuclear_masses(ng)
    call prop_test_gradient(ng, S, D, F)
    call get_fo_geo_perturbed_matrices(ng, S, D, F, Sg, Dg, Fg)
    call contract_hessian(ng, S, D, F, Sg, Dg, Fg)

    Sg = 0
    Dg = 0
    Fg = 0

  end subroutine

  subroutine prop_test_cubicff(ng, S, D, F)

    integer,       intent(in) :: ng
    type(matrix),  intent(in) :: S, D, F

    type(matrix)              :: Sg(ng), Dg(ng), Fg(ng)

    call print_nuclear_masses(ng)
    call prop_test_gradient(ng, S, D, F)
    call get_fo_geo_perturbed_matrices(ng, S, D, F, Sg, Dg, Fg)
    call contract_hessian(ng, S, D, F, Sg, Dg, Fg)
    call contract_cubicff(ng, S, D, F, Sg, Dg, Fg)

    Sg = 0
    Dg = 0
    Fg = 0

  end subroutine

  subroutine print_nuclear_masses(ng)

    integer, intent(in) :: ng

    real(8)             :: charges(ng/3), coords(ng)
    real(8)             :: masses(ng/3), totmass, origin(3)
    integer             :: isotopes(ng/3)
    integer             :: i

    do i = 1, ng/3
       charges(i)  = get_nuc_charge(i)
       isotopes(i) = get_nuc_isotope(i)
       coords((i-1)*3 + 1) = get_nuc_xyz(1, i)
       coords((i-1)*3 + 2) = get_nuc_xyz(2, i)
       coords((i-1)*3 + 3) = get_nuc_xyz(3, i)
    end do

    call VIBMAS(masses, totmass, isotopes, nint(charges), ng/3, coords, origin, 1)
    open (unit=iounit, file='masses', status='replace', action='write')
    call print_tensor(shape(masses), masses*(1d0,0d0), unit=iounit)
    close (iounit)

  end subroutine




! MR: Working on quartic force field in this routine
! The first routines are ajt routines to get the gradient, Hessian
! and cubic force field. Then comes the new code that is intended to get
! the quartic force field.

  subroutine prop_test_quarticff(ng, S, D, F)
    integer,       intent(in) :: ng
    type(matrix),  intent(in) :: S, D, F
    type(rsp_field) geo4(4)
    type(matrix) DFD, Sg(ng), Dg(ng), Fg(ng), DSDg, FDSg(1), DFDg, A, Ag(ng)
    type(matrix) X(1), FgDS(ng), DgSD(ng), FDSgg, DSDgg, DFDgg
    type(matrix) Dgg(ng,ng), Sgg(ng,ng), Fgg(ng,ng) ! NOT INITIALIZED
    type(matrix) DFDggg2p, DSDggg2p, FDSggg2p, RHS(1), Dh
    complex(8)   gra(ng), tm1(ng), hes(ng,ng), tm2(ng,ng)
    complex(8)   cub(ng,ng,ng), tm3(ng,ng,ng)
    complex(8)   qua(ng,ng,ng,ng), tmp(ng,ng,ng,ng)
    complex(8)   xc_qua(ng,ng,ng,ng)
    character(4) nof(0) !no-field
    integer      i, j, k, l, m, noc(0) !no-comp

    geo4(1) = rsp_field('GEO ', (0d0,0d0), 1, ng)
    geo4(2) = rsp_field('GEO ', (0d0,0d0), 1, ng)
    geo4(3) = rsp_field('GEO ', (0d0,0d0), 1, ng)
    geo4(4) = rsp_field('GEO ', (0d0,0d0), 1, ng)

    call print_nuclear_masses(ng)
    call prop_test_gradient(ng, S, D, F)
    call get_fo_geo_perturbed_matrices(ng, S, D, F, Sg, Dg, Fg)
    call contract_hessian(ng, S, D, F, Sg, Dg, Fg)
    call contract_cubicff(ng, S, D, F, Sg, Dg, Fg)


! MR: BEGIN QUARTIC FORCE FIELD
! USES 1,2 RULE


! CALCULATE Sgg

    call rsp_ovlint(S%nrow, 2, (/'GEO ','GEO '/), (/1,1/), shape(Sgg), Sgg)


!   calculate Fgg
!   =============

    call rsp_oneint(S%nrow, 2, (/'GEO ','GEO '/), (/1,1/), shape(Fgg), Fgg)
    call rsp_twoint(S%nrow, 2, (/'GEO ','GEO '/), (/1,1/), shape(Fgg), D, Fgg)
    call rsp_xcint(D=(/D/), Fgg=Fgg)

!   auxiliary matrices to avoid
!   calculating same integrals twice
    A = tiny(0.0d0)*D
    do j = 1, size(Dg)
       Ag(j) = tiny(0.0d0)*D
    end do

    do i = 1, size(Dg)
       do j = 1, size(Dg)
          Ag(j)%elms_alpha = 0.0d0
       end do
       call rsp_twoint(S%nrow, 1, (/'GEO '/), (/1/), shape(Ag), Dg(i), Ag)
       call rsp_xcint(D=(/D, Dg(i)/), Fg=Ag)
       do j = 1, size(Dg)
          Fgg(i, j) = Fgg(i, j) + Ag(j)
          Fgg(j, i) = Fgg(j, i) + Ag(j)
       end do
       do j = 1, i
          A%elms_alpha = 0.0d0
          call rsp_xcint(D=(/D, Dg(i), Dg(j)/), F=A)
          Fgg(i, j) = Fgg(i, j) + A
          if (i /= j) then
             Fgg(j, i) = Fgg(j, i) + A
          end if
       end do
    end do

    A = 0
    do j = 1, size(Dg)
       Ag(j) = 0
    end do


! CALCULATE Dgg
! DJ triangular loop
    do i = 1, size(Sg)
       do j = i, size(Sg)
          Dgg(i,j) = D*Sg(i)*Dg(j) + D*Sgg(i,j)*D + D*Sg(j)*Dg(i) + &
               Dg(i)*S*Dg(j) + Dg(i)*Sg(j)*D + Dg(j)*S*Dg(i) + Dg(j)*Sg(i)*D
          Dgg(i,j) = Dgg(i,j) - D*S*Dgg(i,j) - Dgg(i,j)*S*D
          call rsp_twoint(S%nrow, 0, (nof), (noc), shape(Fgg(i,j)), Dgg(i,j), Fgg(i,j))
          call rsp_xcint(D=(/D, Dgg(i, j)/), F=Fgg(i, j))
          RHS(1) = F*D*Sgg(i,j) + F*Dg(i)*Sg(j) + F*Dgg(i,j)*S + F*Dg(j)*Sg(i)  &
                 + Fg(i)*D*Sg(j) + Fg(i)*Dg(j)*S + Fgg(i,j)*D*S + Fg(j)*D*Sg(i) &
                 + Fg(j)*Dg(i)*S - S*D*Fgg(i,j) - S*Dg(i)*Fg(j) - S*Dgg(i,j)*F  &
                 - S*Dg(j)*Fg(i) - Sg(i)*D*Fg(j) - Sg(i)*Dg(j)*F - Sgg(i,j)*D*F &
                 - Sg(j)*D*Fg(i) - Sg(j)*Dg(i)*F
          X(1) = 0*RHS(1)
          call rsp_mosolver_exec(RHS(1), (/0d0/), X)
          X(1)=-2d0*X(1)
          RHS(1)=0
          Dh = X(1)*S*D
          Dh = Dh - D*S*X(1)
          X(1)=0
          Dgg(i,j) = Dgg(i,j) + Dh
          call rsp_twoint(S%nrow, 0, (nof), (noc), shape(Fgg(i,j)), Dh, Fgg(i,j))
          call rsp_xcint(D=(/D, Dh/), F=Fgg(i, j))
       end do
    end do
! Symmetrize
    do i = 1, size(Sg)
       do j = i+1, size(Sg)
          Dgg(j,i) = Dgg(i,j)
          Fgg(j,i) = Fgg(i,j)
       end do
    end do


! BEGIN ADDING CONTRIBUTIONS TO QUARTIC FORCE FIELD TENSOR

! A, B, C, D PERTURBATIONS IN SAME BLOCK - NO "OTHER" BLOCKS TREATING
! ONE PERTURBATION SEPARATELY


! NUCLEAR POTENTIAL BLOCK
! MR (01/2012: THE NUCPOT CALL SEEMS TO RETURN THE CORRECT RESULTS)

    tmp = 0
    call rsp_nucpot(geo4, tmp)
    qua = tmp; call print_tensor(shape(tmp), tmp, 'Hnuc_abcd')


! OVL BLOCK: THIS IS THE BLOCK CONTAINING ALL PULAY-TYPE (-SW) CONTRIBUTIONS
! MR (01/2012: THE OVLAVE CALLS SEEM TO RETURN THE CORRECT RESULTS)

    DFD = D*F*D
    call rsp_ovlave(4, (/'GEO ','GEO ','GEO ','GEO '/), (/1,1,1,1/), &
                    shape(tmp), DFD, tmp); DFD=0
    qua = qua + tmp; call print_tensor(shape(tmp), tmp, '-SabcdDFD')


    do i = 1, size(Dg)
       tmp = 0.0d0
       DFDg = Dg(i)*F*D + D*Fg(i)*D + D*F*Dg(i) ! Wa
       call rsp_ovlave(3, (/'GEO ','GEO ','GEO '/), (/1,1,1/), shape(tmp(i,:,:,:)), &
                       DFDg, tmp(i,:,:,:))
       DFDg=0
       do j = 1, size(Dg)
          do k = 1, size(Dg)
             do l = 1, size(Dg)
                qua(i, j, k, l) = qua(i, j, k, l) + tmp(i, j, k, l)
                qua(j, i, k, l) = qua(j, i, k, l) + tmp(i, j, k, l)
                qua(j, k, i, l) = qua(j, k, i, l) + tmp(i, j, k, l)
                qua(j, k, l, i) = qua(j, k, l, i) + tmp(i, j, k, l)
             end do
          end do
       end do
    end do

!   note 1:
!   here it can look like 3 terms are missing
!   but this is fine (2n + 1 rule)
!   they are omitted because the rule
!   used is the "1, 2" rule, i.e., contributions with contractions with perturbed
!   density matrices involving perturbation "a" will always be omitted if the order
!   of the perturbed density matrix is more than 1.
    do j = 1, size(Dg)
       do i = 1, j
          tmp = 0.0d0
          DFDgg = Dgg(i,j)*F*D + Dg(i)*Fg(j)*D + Dg(i)*F*Dg(j) + &
                  D*F*Dgg(i,j) + D*Fg(j)*Dg(i) + Dg(j)*F*Dg(i) + &
                  Dg(j)*Fg(i)*D + D*Fgg(i,j)*D + D*Fg(i)*Dg(j) ! Wcd
          call rsp_ovlave(2, (/'GEO ','GEO '/), (/1,1/), shape(tmp(:,:,i,j)), &
                          DFDgg, tmp(:,:,i,j)); DFDgg=0
          do k = 1, size(Dg)
             do l = 1, size(Dg)
                qua(k, l, i, j) = qua(k, l, i, j) + tmp(k, l, i, j)
                qua(k, i, l, j) = qua(k, i, l, j) + tmp(k, l, i, j)
                qua(k, i, j, l) = qua(k, i, j, l) + tmp(k, l, i, j)
                if (i /= j) then
                   qua(k, l, j, i) = qua(k, l, j, i) + tmp(k, l, i, j)
                   qua(k, j, l, i) = qua(k, j, l, i) + tmp(k, l, i, j)
                   qua(k, j, i, l) = qua(k, j, i, l) + tmp(k, l, i, j)
                end if
             end do
          end do
       end do
    end do

! W2'-TYPE CALL
    do k = 1, size(Dg)
      do j = 1, size(Dg)
        do i = 1, size(Dg)
            DFDggg2p = Dgg(i,j)*(Fg(k)*D+F*Dg(k))+(D*Fg(k)+Dg(k)*F)*Dgg(i,j) + &
                      Dgg(i,k)*(Fg(j)*D+F*Dg(j))+(D*Fg(j)+Dg(j)*F)*Dgg(i,k) + &
                      Dgg(j,k)*(Fg(i)*D+F*Dg(i))+(D*Fg(i)+Dg(i)*F)*Dgg(j,k) + &
                      Dg(i)*(Fgg(j,k)*D+Fg(j)*Dg(k))+(D*Fgg(j,k)+Dg(k)*Fg(j))*Dg(i) + &
                      Dg(j)*(Fgg(i,k)*D+Fg(k)*Dg(i))+(D*Fgg(i,k)+Dg(i)*Fg(k))*Dg(j) + &
                      Dg(k)*(Fgg(i,j)*D+Fg(i)*Dg(j))+(D*Fgg(i,j)+Dg(j)*Fg(i))*Dg(k)
            call rsp_ovlave(1, (/'GEO '/), (/1/), shape(tmp(:,i,j,k)), &
                            DFDggg2p, tmp(:,i,j,k)); DFDggg2p=0
        end do
      end do
    end do
    qua = qua + tmp; call print_tensor(shape(tmp), tmp, '-SaDFDbcd2p')


! ONEEL BLOCK: THIS IS THE BLOCK CONTAINING ALL "ONE-ELECTRON AVERAGE" CONTRIBUTIONS
! MR (01/2012: THE ONEAVE CALLS SEEM TO RETURN THE CORRECT RESULTS)

! E0-TYPE CALLS


    tmp = 0.0d0
    call rsp_oneave(4, (/'GEO ','GEO ','GEO ','GEO '/), (/1,1,1,1/), shape(tmp), D, tmp)
    qua = qua + tmp

! E1-TYPE CALLS
! Dg CALLS
    do i = 1, size(Dg)
       tmp = 0.0d0
       call rsp_oneave(3, (/'GEO ','GEO ','GEO '/), (/1,1,1/), shape(tmp(i,:,:,:)), &
                       Dg(i), tmp(i,:,:,:))
       do j = 1, size(Dg)
          do k = 1, size(Dg)
             do l = 1, size(Dg)
                qua(i, j, k, l) = qua(i, j, k, l) + tmp(i, j, k, l)
                qua(j, i, k, l) = qua(j, i, k, l) + tmp(i, j, k, l)
                qua(j, k, i, l) = qua(j, k, i, l) + tmp(i, j, k, l)
                qua(j, k, l, i) = qua(j, k, l, i) + tmp(i, j, k, l)
             end do
          end do
       end do
    end do

!   (see note 1)
    do j = 1, size(Dg)
       do i = 1, j
          tmp = 0.0d0
          call rsp_oneave(2, (/'GEO ','GEO '/), (/1,1/), shape(tmp(:,:,i,j)), &
                          Dgg(i,j), tmp(:,:,i,j))
          do k = 1, size(Dg)
             do l = 1, size(Dg)
                qua(k, l, i, j) = qua(k, l, i, j) + tmp(k, l, i, j)
                qua(k, i, l, j) = qua(k, i, l, j) + tmp(k, l, i, j)
                qua(k, i, j, l) = qua(k, i, j, l) + tmp(k, l, i, j)
                if (i /= j) then
                   qua(k, l, j, i) = qua(k, l, j, i) + tmp(k, l, i, j)
                   qua(k, j, l, i) = qua(k, j, l, i) + tmp(k, l, i, j)
                   qua(k, j, i, l) = qua(k, j, i, l) + tmp(k, l, i, j)
                end if
             end do
          end do
       end do
    end do


! TWOEL BLOCK: THIS IS THE BLOCK CONTAING ALL THE "TWO-ELECTRON AVERAGE" CONTRIBUTIONS
! MR (01/2012: THE TWOAVE CALLS SEEM TO RETURN THE CORRECT RESULTS)

! E0-TYPE CALLS
    call rsp_twoave(4, (/'GEO ','GEO ','GEO ','GEO '/), (/1,1,1,1/), shape(tmp), D, D, tmp)
    qua = qua + tmp/2; call print_tensor(shape(tmp), tmp/2, 'Gabcd(D)D/2')

! E1-TYPE CALLS
! D-Dg CALLS
    do i = 1, size(Dg)
       tmp = 0.0d0
       call rsp_twoave(3, (/'GEO ','GEO ','GEO '/), (/1,1,1/), shape(tmp(i,:,:,:)), &
                       D, Dg(i), tmp(i,:,:,:))
       do j = 1, size(Dg)
          do k = 1, size(Dg)
             do l = 1, size(Dg)
                qua(i, j, k, l) = qua(i, j, k, l) + tmp(i, j, k, l)
                qua(j, i, k, l) = qua(j, i, k, l) + tmp(i, j, k, l)
                qua(j, k, i, l) = qua(j, k, i, l) + tmp(i, j, k, l)
                qua(j, k, l, i) = qua(j, k, l, i) + tmp(i, j, k, l)
             end do
          end do
       end do
    end do

! D-Dgg CALLS
!   (see note 1)
    do j = 1, size(Dg)
       do i = 1, j
          tmp = 0.0d0
          call rsp_twoave(2, (/'GEO ','GEO '/), (/1,1/), shape(tmp(:,:,i,j)), &
                          D, Dgg(i,j), tmp(:,:,i,j))
          do k = 1, size(Dg)
             do l = 1, size(Dg)
                qua(k, l, i, j) = qua(k, l, i, j) + tmp(k, l, i, j)
                qua(k, i, l, j) = qua(k, i, l, j) + tmp(k, l, i, j)
                qua(k, i, j, l) = qua(k, i, j, l) + tmp(k, l, i, j)
                if (i /= j) then
                   qua(k, l, j, i) = qua(k, l, j, i) + tmp(k, l, i, j)
                   qua(k, j, l, i) = qua(k, j, l, i) + tmp(k, l, i, j)
                   qua(k, j, i, l) = qua(k, j, i, l) + tmp(k, l, i, j)
                end if
             end do
          end do
       end do
    end do


! E2-TYPE CALLS
! Dg-Dg CALLS
    do i = 1, size(Dg)
       do j = 1, i
          tmp = 0.0d0
          call rsp_twoave(2, (/'GEO ','GEO '/), (/1,1/), shape(tmp(i,j,:,:)), &
                          Dg(i), Dg(j), tmp(i,j,:,:))
          do k = 1, size(Dg)
             do l = 1, size(Dg)
                qua(i, j, k, l) = qua(i, j, k, l) + tmp(i, j, k, l)
                qua(i, k, j, l) = qua(i, k, j, l) + tmp(i, j, k, l)
                qua(i, k, l, j) = qua(i, k, l, j) + tmp(i, j, k, l)
                qua(k, i, j, l) = qua(k, i, j, l) + tmp(i, j, k, l)
                qua(k, i, l, j) = qua(k, i, l, j) + tmp(i, j, k, l)
                qua(k, l, i, j) = qua(k, l, i, j) + tmp(i, j, k, l)
                if (i /= j) then
                   qua(j, i, k, l) = qua(j, i, k, l) + tmp(i, j, k, l)
                   qua(j, k, i, l) = qua(j, k, i, l) + tmp(i, j, k, l)
                   qua(j, k, l, i) = qua(j, k, l, i) + tmp(i, j, k, l)
                   qua(k, j, i, l) = qua(k, j, i, l) + tmp(i, j, k, l)
                   qua(k, j, l, i) = qua(k, j, l, i) + tmp(i, j, k, l)
                   qua(k, l, j, i) = qua(k, l, j, i) + tmp(i, j, k, l)
                end if
             end do
          end do
       end do
    end do

    do k = 1, size(Dg)
      do j = 1, k
        do i = 1, size(Dg)
           tmp = 0.0d0
           call rsp_twoave(1, (/'GEO '/), (/1/), shape(tmp(i,j,k,:)), &
                         Dg(i), Dgg(j,k), tmp(i,j,k,:))
           do l = 1, size(Dg)
              qua(i, j, k, l) = qua(i, j, k, l) + tmp(i, j, k, l)
              qua(i, j, l, k) = qua(i, j, l, k) + tmp(i, j, k, l)
              qua(i, l, j, k) = qua(i, l, j, k) + tmp(i, j, k, l)
              qua(l, i, j, k) = qua(l, i, j, k) + tmp(i, j, k, l)
              qua(l, k, i, j) = qua(l, k, i, j) + tmp(i, j, k, l)
              qua(l, k, j, i) = qua(l, k, j, i) + tmp(i, j, k, l)
              if (j /= k) then
                 qua(i, k, j, l) = qua(i, k, j, l) + tmp(i, j, k, l)
                 qua(i, k, l, j) = qua(i, k, l, j) + tmp(i, j, k, l)
                 qua(i, l, k, j) = qua(i, l, k, j) + tmp(i, j, k, l)
                 qua(l, i, k, j) = qua(l, i, k, j) + tmp(i, j, k, l)
                 qua(l, j, i, k) = qua(l, j, i, k) + tmp(i, j, k, l)
                 qua(l, j, k, i) = qua(l, j, k, i) + tmp(i, j, k, l)
              end if
           end do
        end do
      end do
    end do


!   xc contribution to the energy derivative
!   ========================================

    xc_qua = 0.0d0
    call rsp_xcave(pert='gggg', res=xc_qua, D=D, Dg=Dg, Dgg=Dgg)
    call print_tensor(shape(xc_qua), xc_qua, 'xc_qua')
    qua = qua + xc_qua


! IDEMPOTENCY BLOCK
    ! 'Zeta g'
    do i = 1, size(Dg)
       FgDS(i) = Fg(i)*D*S + S*D*Fg(i) - Fg(i) - F*D*Sg(i) - Sg(i)*D*F
    end do

    ! 'Zeta g' WITH DSDggg2p
    do k = 1, size(Dg)
      do j = 1, size(Dg)
        do i = 1, size(Dg)
          DSDggg2p = Dgg(i,j)*(Sg(k)*D + S*Dg(k))+(D*Sg(k)+Dg(k)*S)*Dgg(i,j) + &
                     Dgg(i,k)*(Sg(j)*D + S*Dg(j))+(D*Sg(j)+Dg(j)*S)*Dgg(i,k) + &
                     Dgg(j,k)*(Sg(i)*D + S*Dg(i))+(D*Sg(i)+Dg(i)*S)*Dgg(j,k) + &
                     Dg(i)*(Sgg(j,k)*D+Sg(j)*Dg(k))+(D*Sgg(j,k)+Dg(k)*Sg(j))*Dg(i) + &
                     Dg(j)*(Sgg(i,k)*D+Sg(k)*Dg(i))+(D*Sgg(i,k)+Dg(i)*Sg(k))*Dg(j) + &
                     Dg(k)*(Sgg(i,j)*D+Sg(i)*Dg(j))+(D*Sgg(i,j)+Dg(j)*Sg(i))*Dg(k)

            do m = 1, size(Dg)
              tmp(m,i,j,k) = -tr(FgDS(m),DSDggg2p)
            end do; DSDggg2p=0

        end do
      end do
    end do
 FgDS=0

    qua = qua + tmp; call print_tensor(shape(tmp), tmp, '-FaDS DSDbcd2p')

! SELF-CONSISTENCY BLOCK
    ! 'Lambda g'
    do i = 1, size(Dg)
       DgSD(i) = Dg(i)*S*D - D*S*Dg(i)
    end do

    do k = 1, size(Dg)
      do j = 1, size(Dg)
        do i = 1, size(Dg)
          FDSggg2p = Fgg(i,j)*(Dg(k)*S+D*Sg(k))-(S*Dg(k)+Sg(k)*D)*Fgg(i,j) + &
                     Fgg(i,k)*(Dg(j)*S+D*Sg(j))-(S*Dg(j)+Sg(j)*D)*Fgg(i,k) + &
                     Fgg(j,k)*(Dg(i)*S+D*Sg(i))-(S*Dg(i)+Sg(i)*D)*Fgg(j,k) + &
                     Fg(i)*(Dgg(j,k)*S+Dg(j)*Sg(k)+Dg(k)*Sg(j)+D*Sgg(j,k)) - &
                     (S*Dgg(j,k)+Sg(k)*Dg(j)+Sg(j)*Dg(k)+Sgg(j,k)*D)*Fg(i) + &
                     Fg(j)*(Dgg(i,k)*S+Dg(i)*Sg(k)+Dg(k)*Sg(i)+D*Sgg(i,k)) - &
                     (S*Dgg(i,k)+Sg(k)*Dg(i)+Sg(i)*Dg(k)+Sgg(i,k)*D)*Fg(j) + &
                     Fg(k)*(Dgg(i,j)*S+Dg(j)*Sg(i)+Dg(i)*Sg(j)+D*Sgg(i,j)) - &
                     (S*Dgg(i,j)+Sg(i)*Dg(j)+Sg(j)*Dg(i)+Sgg(i,j)*D)*Fg(k) + &
                     F*(Dgg(i,j)*Sg(k)+Dgg(i,k)*Sg(j)+Dgg(j,k)*Sg(i)) - &
                     (Sg(k)*Dgg(i,j)+Sg(j)*Dgg(i,k)+Sg(i)*Dgg(j,k))*F + &
                     F*(Dg(i)*Sgg(j,k)+Dg(j)*Sgg(i,k)+Dg(k)*Sgg(i,j)) - &
                     (Sgg(j,k)*Dg(i)+Sgg(i,k)*Dg(j)+Sgg(i,j)*Dg(k))*F
          do m = 1, size(Dg)
             tmp(m,i,j,k) = -tr(DgSD(m),FDSggg2p)
          end do

        end do
      end do
    end do
    FDSggg2p=0
    DgSD=0

    qua = qua + tmp; call print_tensor(shape(tmp), tmp, '-DaSD FDSbcd2p')

    call print_tensor(shape(qua), qua, 'quarticff = Egggg')

! MR: NOTE INDICES OF OUTPUT (EVEN THOUGH THEY DON'T REALLY MATTER A LOT HERE IN Egggg)
    open (unit=iounit, file='quarticff', status='replace', action='write')
    write (iounit,*)
    do j = 1, size(qua,3)
      do i = 1, size(qua,4)
         call print_tensor(shape(qua(:,:,j,i)), qua(:,:,j,i), unit=iounit)
      end do
    end do
    close (iounit)

  end subroutine

  subroutine prop_test_diphes(ng, S, D, F)

    integer,       intent(in) :: ng
    type(matrix),  intent(in) :: S, D, F
    type(rsp_field) geo2_el(3)
    type(matrix) DFD, Sg(ng), Dg(ng), Fg(ng)
    type(matrix) Df(3), Ff(3), FDSf(1), DFDf
    type(matrix) X(1), FgDS(ng), DgSD(ng)
    type(matrix) FDSgf, DSDgf, DFDgf
    complex(8)   dip(3), pol(3,3), gra(ng), tm1(ng)
    complex(8)   Eggf(ng,ng,3), tmp(ng,ng,3)
    character(4) nof(0) !no-field
    integer      i, j, k, noc(0) !no-comp
    geo2_el(:2) = rsp_field('GEO ', (0d0,0d0), 1, ng)
    geo2_el(3)  = rsp_field('EL  ', (0d0,0d0), 1, 3)
    !-------------------------------------------------
    dip = 0
    call rsp_nucpot(geo2_el(3:3), dip)
    call rsp_oneave(1, (/'EL  '/), (/1/), shape(dip), D, dip)
    call print_tensor(shape(dip), -dip, 'dipole moment')
    ! gradient, first nuclear contribution
    tm1 = 0
    call rsp_nucpot(geo2_el(:1), tm1)
    gra = tm1; call print_tensor(shape(tm1), tm1, 'hnuc_g')
    ! overlap contribution
    DFD = D*F*D
    call rsp_ovlave(1, (/'GEO '/), (/1/), shape(tm1), DFD, &
                    tm1, (/(0d0,0d0)/), D); DFD=0
    gra = gra + tm1; call print_tensor(shape(tm1), tm1, '-SgDFD')
    ! 1-electron contribution
    call rsp_oneave(1, (/'GEO '/), (/1/), shape(tm1), D, tm1)
    gra = gra + tm1; call print_tensor(shape(tm1), tm1, 'HgD')
    ! 2-electron contribution
    call rsp_twoave(1, (/'GEO '/), (/1/), shape(tm1), D, D, tm1)
    gra = gra + tm1/2; call print_tensor(shape(tm1), tm1/2, 'Gg(D)D/2')
    ! Kohn-Sham exchange correlation average
!   call rsp_xcave(mol, 1, (/'GEO '/), (/1/), shape(tm1), 1, (/D/), tm1)
    gra = gra + tm1; call print_tensor(shape(tm1), tm1, 'Exc_g(D)')
    ! print
    call print_tensor(shape(gra), gra, 'gradient = Eg')
    ! print formatted to file gradient
    open (unit=iounit, file='gradient', status='replace', action='write')
    call print_tensor(shape(gra), gra, unit=iounit)
    close (iounit)
    ! test dipole gradient average
    call rsp_oneave(2, (/'GEO ','EL  '/), (/1,1/), shape(tmp(:,1,:)), &
                       D, tmp(:,1,:))
    call print_tensor(shape(tmp(:,1,:)), tmp(:,1,:), 'dpgave form rsp_oneave')

! density independent contribution to Ff
    call rsp_oneint(S%nrow, 1, (/'EL  '/), (/1/), shape(Ff), Ff)
! solve equations for Df and construct Ff
    pol = 0
    do i = 1, size(Df)
       Df(i) = 0d0*D
       FDSf(1) = Ff(i)*D*S
       FDSf(1) = FDSf(1) - S*D*Ff(i)
       X(1) = 0*FDSf(1)
       call mat_ensure_alloc(X(1))
       call rsp_mosolver_exec(FDSf(1), (/0d0/), X)
       X(1)=-2d0*X(1); FDSf(1)=0
       Df(i) = Df(i) + X(1)*S*D
       Df(i) = Df(i) - D*S*X(1); X(1)=0
       ! Df contribution to Ff
       call rsp_twoint(S%nrow, 0, nof, noc, noc, Df(i), Ff(i:i))
       !polarzability
       call rsp_oneave(1, (/'EL  '/), (/1/), shape(pol(:,i)), Df(i), pol(:,i))
    end do
    call print_tensor(shape(pol), -pol, 'polarizability')



    call get_fo_geo_perturbed_matrices(ng, S, D, F, Sg, Dg, Fg)
    call contract_hessian(ng, S, D, F, Sg, Dg, Fg)


    !------------------------------------
    ! Dipole Hessian contributions start here
    !
    ! overlap contribution
    DFD = D*F*D
    !
    do i = 1, size(Df)
       DFDf = Df(i)*F*D + D*Ff(i)*D + D*F*Df(i)
       call rsp_ovlave(2, (/'GEO ','GEO '/), (/1,1/), shape(tmp(:,:,i)), &
                       DFDf, tmp(:,:,i))
       DFDf=0
    end do
    Eggf = tmp; call print_tensor(shape(tmp), tmp, '-SabDFDf')
    !
    do j = 1, size(Df)
       do i = 1, size(Dg)
          DFDgf = Df(j)*Fg(i)*D + Df(j)*F*Dg(i) + D*Ff(j)*Dg(i) &
                + Dg(i)*Ff(j)*D + Dg(i)*F*Df(j) + D*Fg(i)*Df(j)
          call rsp_ovlave(1, (/'GEO '/), (/1/), shape(tmp(:,i,j)), &
                          DFDgf, tmp(:,i,j))
          DFDgf=0
       end do
    end do
    Eggf = Eggf + tmp; call print_tensor(shape(tmp), tmp, '-SaDFDbf1')

    ! 1-electron contribution
    ! fixme: finite field calculation in rsp_oneave
    call rsp_oneave(3, (/'GEO ','GEO ','EL  '/), (/1,1,1/), shape(tmp), D, tmp)
    Eggf = Eggf + tmp; call print_tensor(shape(tmp), tmp, 'HabfD')
    do i = 1, size(Df)
       call rsp_oneave(2, (/'GEO ','GEO '/), (/1,1/), shape(tmp(:,:,i)), &
                       Df(i), tmp(:,:,i))
    end do
    Eggf = Eggf + tmp; call print_tensor(shape(tmp), tmp, 'HabDf')
    do i = 1, size(Dg)
       tmp = 0.0d0
       call rsp_oneave(2, (/'GEO ','EL  '/), (/1,1/), shape(tmp(i,:,:)), &
                       Dg(i), tmp(i,:,:))
       do j = 1, size(Dg)
          do k = 1, size(Df)
             Eggf(i, j, k) = Eggf(i, j, k) + tmp(i, j, k)
             Eggf(j, i, k) = Eggf(j, i, k) + tmp(i, j, k)
          end do
       end do
    end do

    ! 2-electron contribution
    do i = 1, size(Df)
       call rsp_twoave(2, (/'GEO ','GEO '/), (/1,1/), shape(tmp(:,:,i)), &
                       D, Df(i), tmp(:,:,i))
    end do
    Eggf = Eggf + tmp; call print_tensor(shape(tmp), tmp, 'Gab(D)Df')
    do j = 1, size(Df)
       do i = 1, size(Dg)
          tmp = 0.0d0
          call rsp_twoave(1, (/'GEO '/), (/1/), shape(tmp(:,i,j)), &
                          Dg(i), Df(j), tmp(:,i,j))
          do k = 1, size(Dg)
             Eggf(k, i, j) = Eggf(k, i, j) + tmp(k, i, j)
             Eggf(i, k, j) = Eggf(i, k, j) + tmp(k, i, j)
          end do
       end do
    end do

    ! idempotency multiplier contribution
    do i = 1, size(Dg)
       FgDS(i) = Fg(i)*D*S + S*D*Fg(i) - Fg(i) - F*D*Sg(i) - Sg(i)*D*F
    end do
    do k = 1, size(Df)
       do j = 1, size(Dg)
          DSDgf = Df(k)*(Sg(j)*D + S*Dg(j))+ (Dg(j)*S + D*Sg(j))*Df(k)
          do i = 1, size(Dg)
             tmp(i,j,k) = -tr(FgDS(i),DSDgf)
          end do; DSDgf=0
       end do
    end do; FgDS=0;
    Eggf = Eggf + tmp; call print_tensor(shape(tmp), tmp, '-FaDS DSDbf1')
    ! self-consistency multiplier contribution
    do i = 1, size(Dg)
       DgSD(i) = Dg(i)*S*D - D*S*Dg(i)
    end do
    do k = 1, size(Df)
       do j = 1, size(Dg)
          FDSgf = (Ff(k)*Dg(j)+Fg(j)*Df(k))*S - Sg(j)*(D*Ff(k) + Df(k)*F) &
                + (Ff(k)*D+F*Df(k))*Sg(j) - S*(Dg(j)*Ff(k) + Df(k)*Fg(j))
          do i = 1, size(Dg)
             tmp(i,j,k) = -tr(DgSD(i),FDSgf)
          end do; FDSgf=0
       end do
    end do; DgSD=0;
    Eggf = Eggf + tmp; call print_tensor(shape(tmp), tmp, '-DaSD FDSbf1')
    ! print
    call print_tensor(shape(Eggf), Eggf, 'diphes = Eggf')
    ! print formatted to file diphes
    open (unit=iounit, file='diphes', status='replace', action='write')
    write (iounit,*)
    do i = 1, size(Eggf,3)
       !call print_tensor(shape(Eggf(:,:,i)), Eggf(:,:,i))
       call print_tensor(shape(Eggf(:,:,i)), Eggf(:,:,i), unit=iounit)
    end do
    close (iounit)
    ! free
    Sg=0; Dg=0; Fg=0
  end subroutine


  !> ajt This is used for printing response function tensors
  !> TODO: Add comp. lables argument. Make space between blocks of rows
  subroutine print_tensor(dims, tensor, title, freqs, colwidth, unit)
    integer,                intent(in) :: dims(:)
    complex(8),             intent(in) :: tensor(*) !(product(dims))
    character(*), optional, intent(in) :: title
    integer,      optional, intent(in) :: colwidth, unit
    complex(8),   optional, intent(in) :: freqs(:)
    integer     :: uni, siz, dec, cwid
    character*9 :: fmt
    uni = 6
    if (present(unit)) uni = unit
    if (present(title)) then
       if (present(freqs)) then
          fmt = '(a,nf6.3)'
          write (fmt(4:4),'(i1)') size(freqs)
          write (uni,fmt) title // ' w:', dreal(freqs)
       else
          write (uni,'(a)') title
       end if
    end if
    siz = product(dims)
    dec = max(2, 1+ceiling(log10(maxval(abs(tensor(1:siz))))))
    !dec = max(dec,ceiling(log10( maxval(tensor(1:siz)))))
    cwid = 15
    if (present(colwidth)) cwid = max(colwidth,dec+2)
    dec = cwid - dec - 1
    fmt = '(fww.dd)'
    write (fmt(3:7),'(i2,a1,i2)') cwid, '.', dec
    if (any(dimag(tensor(1:siz)) /= 0)) write (uni,'(a)') '    (real part)'
    call subr(dims, dreal(tensor(1:siz)), uni, cwid, fmt)
    if (any(dimag(tensor(1:siz)) /= 0)) then
       write (uni,'(a)') '    (imaginary part)'
       call subr(dims, dimag(tensor(1:siz)), uni, cwid, fmt)
    end if
    write (uni,'()') !final blank line
  contains
    subroutine subr(dims, tnz, unit, cwid, fmt)
      integer,     intent(in) :: dims(:), unit, cwid
      real(8),     intent(in) :: tnz(*)
      character*8, intent(in) :: fmt
      integer :: i, j, l
      character*(dims(1)*(cwid+1)) :: line
      do i = 1, product(dims), dims(1)
         l = 1
         do j = 1, dims(1)
            line(l:l) = ' '
            write (line(l+1:l+cwid), fmt) tnz(i+j-1)
            l = l + 1 + cwid
         end do
         write (unit,'(a)') line
      end do
    end subroutine
  end subroutine


end module
