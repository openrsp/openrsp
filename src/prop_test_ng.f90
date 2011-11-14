!> @file
!> Contains module prop_test

!> ajt/radovan: Response-related testing routines and some calculations 
module prop_test_ng

  use matrix_defop_ng
  use rsp_contribs_ng
  use rsp_equations_ng

  implicit none

  public prop_test_gradient
  public prop_test_hessian
  public prop_test_cubicff
  public prop_test_diphes

  private

  !physical constants for conversion to non-atomic units
  real(8), parameter :: cm1 = 1/219474.631371d0, & !1 centimeter-to-minus-one in au
                        cvl = 137.03599907d0,    & !c-the-velocity-of-light in au
                        nm  = 10/0.52917706d0,   & !1 nanometer in au
                        pi  = 3.14159265358979323846D0 !acos(-1.d0)

  !field component lables for printing
  character(2) :: fc(3) = (/'Fx','Fy','Fz'/), &
                  bc(3) = (/'Bx','By','Bz'/)

  ! unit number for IO in this file
  integer, parameter :: iounit = 345645

contains


  !> Calculate the gradient: dE/dR  =  dh_nuc/dR  +  Tr dH/dR D
  !>            +  1/2 Tr dG/dR(D) D  +  dExc[D]/dR  -  Tr dS/dR DFD
  subroutine prop_test_gradient(mol, S, D, F)
    type(rsp_cfg), intent(in) :: mol
    type(matrix),  intent(in) :: S, D, F
    complex(8), dimension(3*mol%natom) :: gra, tmp
    type(matrix) :: DFD
    ! first nuclear contribution
    call rsp_nucpot(mol, 1, (/'GEO '/), (/(0d0,0d0)/), (/1/), &
                    shape(tmp), tmp)
    gra = tmp
    call print_tensor(shape(tmp), tmp, 'nucpot')
    ! overlap contribution
    DFD = D*F*D
    call rsp_ovlave(mol, 1, (/'GEO '/), (/1/), shape(tmp), DFD, &
                    tmp, (/(0d0,0d0)/), D)
    DFD = 0 !free
    gra = gra + tmp
    call print_tensor(shape(tmp), tmp, 'ovlave')
    ! 1-electron contribution
    call rsp_oneave(mol, 1, (/'GEO '/), (/1/), shape(tmp), D, tmp)
    gra = gra + tmp
    call print_tensor(shape(tmp), tmp, 'oneave')
    ! 2-electron contribution
    call rsp_twoave(mol, 1, (/'GEO '/), (/1/), shape(tmp), D, D, tmp)
    gra = gra + tmp/2
    call print_tensor(shape(tmp), tmp/2, 'twoave')
    ! Kohn-Sham exchange correlation average
    call rsp_excave(mol, 1, (/'GEO '/), (/1/), shape(tmp), 1, (/D/), tmp)
    gra = gra + tmp
    call print_tensor(shape(tmp), tmp, 'excave')
    ! print
    call print_tensor((/3*mol%natom/), gra, 'gradient = Eg')
  end subroutine


  subroutine prop_test_hessian(mol, S, D, F)
    use dalton_ifc_ng, only: rsp_mosolver_exec
    type(rsp_cfg), intent(in) :: mol
    type(matrix),  intent(in) :: S, D, F
    integer, parameter :: j = 1
    complex(8)    hes(3*mol%natom,1) !3*mol%natom)
    complex(8)    tmp(3*mol%natom,1) !3*mol%natom)
    type(matrix)  DFD, DFDg(1)
    type(matrix)  Sg(1), Dg(1), Fg(1), DSDg(1), FDSg(1), X(1)
    character(4)  nof(0)
    complex(8)    w, freq(1)
    integer       i, noc(0)
    w = (0d0,0d0)
!     call mat_print(D, label='D', width=15)
!     call mat_print(F, label='F', width=15)
    ! 1-electron contribution
    call rsp_oneint(mol, 1, (/'GEO '/), (/j/), shape(Fg), Fg)
    ! overlap contribution
    call rsp_ovlint(mol, 1, (/'GEO '/), (/j/), shape(Sg), Sg, (/w/), Fg)
    ! 2-electron contribution
    call rsp_twoint(mol, 1, (/'GEO '/), (/j/), shape(Fg), D, Fg)
    ! Kohn-Sham exchange correlation average
    call rsp_excint(mol, 1, (/'GEO '/), (/j/), shape(Fg), 1, (/D/), Fg)
    ! print
!     call mat_print(Sg(1), label='Sgx')
    ! solve equations manually, print
    do i = 1, size(Dg)
       FDSg(i) = (F-w/2*S)*D*Sg(i)
       FDSg(i) = FDSg(i) - Sg(i)*D*(F+w/2*S)
       Dg(i)   = -D*Sg(i)*D
       X(1)    = FDSg(i)
       FDSg(i) = X(1)*D*S
       FDSg(i) = FDSg(i) - S*D*X(1)
       FDSg(i) = FDSg(i) + Fg(i)
       call rsp_twoint(mol, 0, nof, noc, noc, Dg(i), FDSg(i:i))
       X(1)    = FDSg(i)
       FDSg(i) = X(1)*D*S
       FDSg(i) = FDSg(i) - S*D*X(1)
       !call mat_init(X(1), Dg(i), zero=.true.)
       call rsp_mosolver_exec(FDSg(i:i), (/dreal(w)/), X); X(1)=-2d0*X(1)
       ! ajt FIXME if I put these on the same line, I get something else
       Dg(i) = Dg(i) + X(1)*S*D
       Dg(i) = Dg(i) - D*S*X(1); X(1)=0
       call rsp_twoint(mol, 0, nof, noc, noc, Dg(i), Fg(i:i))
    end do; FDSg=0; X=0
!     call mat_print(Dg(1), label='Dgx', width=15)
!     call mat_print(Fg(1), label='Fgx', width=15)
    ! verify equations, print
!     do i = 1, size(Dg)
!        DSDg(i) = Dg(i)*S*D + D*Sg(i)*D + D*S*Dg(i) - Dg(i)
!        call mat_print(DSDg(i), label='DSDg'); DSDg=0
!        FDSg(i) = (F-w/2*S)*(Dg(i)*S + D*Sg(i))
!        FDSg(i) = FDSg(i) - (Sg(i)*D + S*Dg(i))*(F+w/2*S)
!        FDSg(i) = FDSg(i) + Fg(i)*D*S
!        FDSg(i) = FDSg(i) - S*D*Fg(i)
!        call mat_print(FDSg(i), label='FDSg'); FDSg=0
!     end do
    ! contract
    call rsp_nucpot(mol, 2, (/'GEO ','GEO '/), (0d0,0d0)*(/0,0/), &
                    (/1,j/), shape(tmp), tmp)
    hes = tmp
    call print_tensor(shape(tmp), tmp, 'nucpot')
    ! overlap contribution
    DFD = D*F*D
    call rsp_ovlave(mol, 2, (/'GEO ','GEO '/), (/1,j/), shape(tmp), DFD, &
                    tmp); DFD=0 !, (0d0,0d0)*(/0,0/), D)
    hes = hes + tmp; call print_tensor(shape(tmp), tmp, '-SabDFD')
    DFDg(1) = Dg(1)*(F+w/2*S)*D + D*Fg(1)*D + D*(F-w/2*S)*Dg(1)
    call rsp_ovlave(mol, 1, (/'GEO '/), (/1/), shape(tmp(:,1)), DFDg(1), &
                    tmp); DFDg(1)=0 !, (0d0,0d0)*(/0,0/), D)
    hes = hes + tmp; call print_tensor(shape(tmp), tmp, '-SaDFDb')
    ! 1-electron contribution
    call rsp_oneave(mol, 2, (/'GEO ','GEO '/), (/1,j/), shape(tmp), D, tmp)
    hes = hes + tmp; call print_tensor(shape(tmp), tmp, 'HabD')
    call rsp_oneave(mol, 1, (/'GEO '/), (/1/), shape(tmp(:,1)), Dg(1), tmp)
    hes = hes + tmp; call print_tensor(shape(tmp), tmp, 'HaDb')
    ! 2-electron contribution
    call rsp_twoave(mol, 2, (/'GEO ','GEO '/), (/1,j/), shape(tmp), D, D, tmp)
    hes = hes + tmp/2; call print_tensor(shape(tmp), tmp/2, 'Gab(D)D/2')
    call rsp_twoave(mol, 1, (/'GEO '/), (/1/), shape(tmp(:,1)), D, Dg(1), tmp)
    hes = hes + tmp; call print_tensor(shape(tmp), tmp, 'Ga(D)Db')
    ! print
    call print_tensor(shape(hes), hes, 'hessian = Egg')
    ! free
    Sg=0; Dg=0; Fg=0
  end subroutine


  subroutine prop_test_cubicff(mol, ng, S, D, F)
    use matrix_genop_ng, only: mat_alloc, matrix_genop_debug
    use dalton_ifc_ng,   only: rsp_mosolver_exec
    type(rsp_cfg), intent(in) :: mol
    integer,       intent(in) :: ng
    type(matrix),  intent(in) :: S, D, F
    type(matrix) DFD, Sg(ng), Dg(ng), Fg(ng), DSDg, FDSg(1), DFDg
    type(matrix) X(1), FgDS(ng), DgSD(ng), FDSgg, DSDgg, DFDgg
    complex(8)   gra(ng), tm1(ng), hes(ng,ng), tm2(ng,ng)
    complex(8)   cub(ng,ng,ng), tmp(ng,ng,ng)
    real(8)      masses(mol%natom), totmass, origin(3)
    integer      isotop(mol%natom)
    character(4) nof(0) !no-field
    integer      i, j, k, noc(0) !no-comp
    ! print nuclear masses, for use when diagonalizing Hessian
    isotop(:) = 1
    call VIBMAS(masses, totmass, isotop, nint(mol%charge), mol%natom, mol%coord, origin, 1)
    open (unit=iounit, file='masses', status='replace', action='write')
    call print_tensor(shape(masses), masses*(1d0,0d0), unit=iounit)
    close (iounit)
    !-------------------------------------
    ! gradient, first nuclear contribution
    call rsp_nucpot(mol, 1, (/'GEO '/), (/(0d0,0d0)/), (/1/), &
                    shape(tm1), tm1)
    gra = tm1; call print_tensor(shape(tm1), tm1, 'hnuc_g')
    ! overlap contribution
    DFD = D*F*D
    call rsp_ovlave(mol, 1, (/'GEO '/), (/1/), shape(tm1), DFD, &
                    tm1, (/(0d0,0d0)/), D); DFD=0
    gra = gra + tm1; call print_tensor(shape(tm1), tm1, '-SgDFD')
    ! 1-electron contribution
    call rsp_oneave(mol, 1, (/'GEO '/), (/1/), shape(tm1), D, tm1)
    gra = gra + tm1; call print_tensor(shape(tm1), tm1, 'HgD')
    ! 2-electron contribution
    call rsp_twoave(mol, 1, (/'GEO '/), (/1/), shape(tm1), D, D, tm1)
    gra = gra + tm1/2; call print_tensor(shape(tm1), tm1/2, 'Gg(D)D/2')
    ! Kohn-Sham exchange correlation average
    call rsp_excave(mol, 1, (/'GEO '/), (/1/), shape(tm1), 1, (/D/), tm1)
    gra = gra + tm1; call print_tensor(shape(tm1), tm1, 'Exc_g(D)')
    ! print
    call print_tensor(shape(gra), gra, 'gradient = Eg')
    ! print formatted to file gradient
    open (unit=iounit, file='gradient', status='replace', action='write')
    call print_tensor(shape(gra), gra, unit=iounit)
    close (iounit)
    !------------------------------
    ! calculate perturbed integrals
    do i = 1, size(Sg)
       call rsp_ovlint(mol, 1, (/'GEO '/), (/1/), shape(Sg), Sg)
       call rsp_oneint(mol, 1, (/'GEO '/), (/1/), shape(Fg), Fg)
       call rsp_twoint(mol, 1, (/'GEO '/), (/1/), shape(Fg), D, Fg)
    end do
    ! solve equations
    do i = 1, size(Dg)
       FDSg(1) = F*D*Sg(i) - Sg(i)*D*F
       Dg(i) = -D*Sg(i)*D
       FDSg(1) = FDSg(1)*D*S - S*D*FDSg(1) + Fg(i)
       call rsp_twoint(mol, 0, nof, noc, noc, Dg(i), FDSg(1:1))
       FDSg(1) = FDSg(1)*D*S - S*D*FDSg(1)
       X(1) = 0*FDSg(1)
       call mat_alloc(X(1))
       call rsp_mosolver_exec(FDSg(1), (/0d0/), X)
       X(1)=-2*X(1); FDSg(1)=0
       Dg(i) = Dg(i) + X(1)*S*D - D*S*X(1); X(1)=0
       call rsp_twoint(mol, 0, nof, noc, noc, Dg(i), Fg(i:i))
    end do
    !------------------------------------
    ! contract Hessian, nuclear repulsion
    call rsp_nucpot(mol, 2, (/'GEO ','GEO '/), (0d0,0d0)*(/0,0/), &
                    (/1,1/), shape(tm2), tm2)
    hes = tm2; call print_tensor(shape(tm2), tm2, 'hnuc_ab')
    ! overlap contribution
    DFD = D*F*D
    call rsp_ovlave(mol, 2, (/'GEO ','GEO '/), (/1,1/), shape(tm2), DFD, &
                    tm2); DFD=0 !, (0d0,0d0)*(/0,0/), D)
    hes = hes + tm2; call print_tensor(shape(tm2), tm2, '-SabDFD')
    !
    do i = 1, size(Dg)
       DFDg = Dg(i)*F*D + D*Fg(i)*D + D*F*Dg(i)
       call rsp_ovlave(mol, 1, (/'GEO '/), (/1/), shape(tm2(:,i)), DFDg, &
                       tm2(:,i)); DFDg=0 !, (0d0,0d0)*(/0,0/), D)
    end do
    hes = hes + tm2; call print_tensor(shape(tm2), tm2, '-SaDFDb')
    ! 1-electron contribution
    call rsp_oneave(mol, 2, (/'GEO ','GEO '/), (/1,1/), shape(tm2), D, tm2)
    hes = hes + tm2; call print_tensor(shape(tm2), tm2, 'HabD')
    !
    do i = 1, size(Dg)
       call rsp_oneave(mol, 1, (/'GEO '/), (/1/), shape(tm2(:,i)), Dg(i), tm2(:,i))
    end do
    hes = hes + tm2; call print_tensor(shape(tm2), tm2, 'HaDb')
    ! 2-electron contribution
    call rsp_twoave(mol, 2, (/'GEO ','GEO '/), (/1,1/), shape(tm2), D, D, tm2)
    hes = hes + tm2/2; call print_tensor(shape(tm2), tm2/2, 'Gab(D)D/2')
    !
    do i = 1, size(Dg)
       call rsp_twoave(mol, 1, (/'GEO '/), (/1/), shape(tm2(:,i)), D, Dg(i), tm2(:,i))
    end do
    hes = hes + tm2; call print_tensor(shape(tm2), tm2, 'Ga(D)Db')
    ! print
    call print_tensor(shape(hes), hes, 'hessian = Egg')
    ! print formatted to file hessian
    open (unit=iounit, file='hessian', status='replace', action='write')
    call print_tensor(shape(hes), hes, unit=iounit)
    close (iounit)
    !------------------------------------
    ! contract cubicff, nuclear repulsion
    call rsp_nucpot(mol, 3, (/'GEO ','GEO ','GEO '/), (0d0,0d0)*(/0,0,0/), &
                    (/1,1,1/), shape(tmp), tmp)
    cub = tmp; call print_tensor(shape(tmp), tmp, 'Hnuc_abc')
    ! overlap contribution
    DFD = D*F*D
    call rsp_ovlave(mol, 3, (/'GEO ','GEO ','GEO '/), (/1,1,1/), &
                    shape(tmp), DFD, tmp); DFD=0 !, (0d0,0d0)*(/0,0/), D)
    cub = cub + tmp; call print_tensor(shape(tmp), tmp, '-SabcDFD')
    !
    do i = 1, size(Dg)
       DFDg = Dg(i)*F*D + D*Fg(i)*D + D*F*Dg(i)
       call rsp_ovlave(mol, 2, (/'GEO ','GEO '/), (/1,1/), shape(tmp(:,:,i)), &
                       DFDg, tmp(:,:,i)); DFDg=0 !, (0d0,0d0)*(/0,0/), D)
    end do
    cub = cub + tmp; call print_tensor(shape(tmp), tmp, '-SabDFDc')
    !
    do i = 1, size(Dg)
       DFDg = Dg(i)*F*D + D*Fg(i)*D + D*F*Dg(i)
       call rsp_ovlave(mol, 2, (/'GEO ','GEO '/), (/1,1/), shape(tmp(:,i,:)), &
                       DFDg, tmp(:,i,:)); DFDg=0 !, (0d0,0d0)*(/0,0/), D)
    end do
    cub = cub + tmp; call print_tensor(shape(tmp), tmp, '-SacDFDb')
    !
    do j = 1, size(Dg)
       do i = 1, size(Dg)
          DFDgg = Dg(j)*Fg(i)*D + Dg(j)*F*Dg(i) + D*Fg(j)*Dg(i) &
                + Dg(i)*Fg(j)*D + Dg(i)*F*Dg(j) + D*Fg(i)*Dg(j)
          call rsp_ovlave(mol, 1, (/'GEO '/), (/1/), shape(tmp(:,i,j)), &
                          DFDgg, tmp(:,i,j)); DFDgg=0 !, (0d0,0d0)*(/0/), D)
       end do
    end do
    cub = cub + tmp; call print_tensor(shape(tmp), tmp, '-SaDFDbc1')
    ! 1-electron contribution
    call rsp_oneave(mol, 3, (/'GEO ','GEO ','GEO '/), (/1,1,1/), shape(tmp), D, tmp)
    cub = cub + tmp; call print_tensor(shape(tmp), tmp, 'HabcD')
    !
    do i = 1, size(Dg)
       call rsp_oneave(mol, 2, (/'GEO ','GEO '/), (/1,1/), shape(tmp(:,:,i)), &
                       Dg(i), tmp(:,:,i))
    end do
    cub = cub + tmp; call print_tensor(shape(tmp), tmp, 'HabDc')
    !
    do i = 1, size(Dg)
       call rsp_oneave(mol, 2, (/'GEO ','GEO '/), (/1,1/), shape(tmp(:,i,:)), &
                       Dg(i), tmp(:,i,:))
    end do
    cub = cub + tmp; call print_tensor(shape(tmp), tmp, 'HacDb')
    ! 2-electron contribution
    call rsp_twoave(mol, 3, (/'GEO ','GEO ','GEO '/), (/1,1,1/), shape(tmp), D, D, tmp)
    cub = cub + tmp/2; call print_tensor(shape(tmp), tmp/2, 'Gabc(D)D/2')
    !
    do i = 1, size(Dg)
       call rsp_twoave(mol, 2, (/'GEO ','GEO '/), (/1,1/), shape(tmp(:,:,i)), &
                       D, Dg(i), tmp(:,:,i))
    end do
    cub = cub + tmp; call print_tensor(shape(tmp), tmp, 'Gab(D)Dc')
    !
    do i = 1, size(Dg)
       call rsp_twoave(mol, 2, (/'GEO ','GEO '/), (/1,1/), shape(tmp(:,i,:)), &
                       D, Dg(i), tmp(:,i,:))
    end do
    cub = cub + tmp; call print_tensor(shape(tmp), tmp, 'Gac(D)Db')
    !
    do j = 1, size(Dg)
       do i = 1, size(Dg)
          call rsp_twoave(mol, 1, (/'GEO '/), (/1/), shape(tmp(:,i,j)), &
                          Dg(i), Dg(j), tmp(:,i,j))
       end do
    end do
    cub = cub + tmp; call print_tensor(shape(tmp), tmp, 'Ga(Db)Dc')
    ! idempotency multiplier contribution
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
    cub = cub + tmp; call print_tensor(shape(tmp), tmp, '-FaDS DSDbc1')
    ! self-consistency multiplier contribution
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
    cub = cub + tmp; call print_tensor(shape(tmp), tmp, '-DaSD FDSbc1')
    ! other overlap contribution
    do i = 1, size(Dg)
       DFDg = Dg(i)*F*D + D*Fg(i)*D + D*F*Dg(i)
       call rsp_ovlave(mol, 2, (/'GEO ','GEO '/), (/1,1/), shape(tmp(i,:,:)), &
                       DFDg, tmp(i,:,:)); DFDg=0 !, (0d0,0d0)*(/0,0/), D)
    end do
    cub = cub + tmp; call print_tensor(shape(tmp), tmp, '-SbcDFDa')
    ! other 1-electron contribution
    do i = 1, size(Dg)
       call rsp_oneave(mol, 2, (/'GEO ','GEO '/), (/1,1/), shape(tmp(i,:,:)), &
                       Dg(i), tmp(i,:,:))
    end do
    cub = cub + tmp; call print_tensor(shape(tmp), tmp, 'HbcDa')
    !
    do i = 1, size(Dg)
       call rsp_twoave(mol, 2, (/'GEO ','GEO '/), (/1,1/), shape(tmp(i,:,:)), &
                       D, Dg(i), tmp(i,:,:))
    end do
    cub = cub + tmp; call print_tensor(shape(tmp), tmp, 'Gbc(D)Da')
    !
    do j = 1, size(Dg)
       do i = 1, size(Dg)
          call rsp_twoave(mol, 1, (/'GEO '/), (/1/), shape(tmp(i,:,j)), &
                          Dg(i), Dg(j), tmp(i,:,j))
       end do
    end do
    cub = cub + tmp; call print_tensor(shape(tmp), tmp, 'Gb(Da)Dc')
    !
    do j = 1, size(Dg)
       do i = 1, size(Dg)
          call rsp_twoave(mol, 1, (/'GEO '/), (/1/), shape(tmp(i,j,:)), &
                          Dg(i), Dg(j), tmp(i,j,:))
       end do
    end do
    cub = cub + tmp; call print_tensor(shape(tmp), tmp, 'Gc(Da)Db')
    ! print
    call print_tensor(shape(cub), cub, 'cubicff = Eggg')
    ! print formatted to file cubicff
    open (unit=iounit, file='cubicff', status='replace', action='write')
    write (iounit,*)
    do i = 1, size(cub,3)
       call print_tensor(shape(cub(:,:,i)), cub(:,:,i), unit=iounit)
    end do
    close (iounit)
    ! free
    Sg=0; Dg=0; Fg=0
  end subroutine


  subroutine prop_test_diphes(mol, ng, S, D, F)
    use matrix_genop_ng, only: mat_alloc
    use dalton_ifc_ng,   only: rsp_mosolver_exec
    type(rsp_cfg), intent(in) :: mol
    integer,       intent(in) :: ng
    type(matrix),  intent(in) :: S, D, F
    type(matrix) DFD, Sg(ng), Dg(ng), Fg(ng), DSDg, FDSg(1), DFDg
    type(matrix) Df(3), Ff(3), DSDf, FDSf(1), DFDf
    type(matrix) X(1), FgDS(ng), DgSD(ng), FDSgg, DSDgg, DFDgg
    type(matrix) FfDS(3), DfSD(3), FDSgf, DSDgf, DFDgf
    complex(8)   dip(3), pol(3,3), gra(ng), tm1(ng), hes(ng,ng), tm2(ng,ng)
    complex(8)   Eggf(ng,ng,3), tmp(ng,ng,3)
    character(4) nof(0) !no-field
    integer      i, j, k, noc(0) !no-comp
    !-------------------------------------
    dip = 0
    call rsp_oneave(mol, 1, (/'EL  '/), (/1/), shape(dip), D, dip)
    call print_tensor(shape(dip), -dip, 'dipole moment')
    ! gradient, first nuclear contribution
    call rsp_nucpot(mol, 1, (/'GEO '/), (/(0d0,0d0)/), (/1/), &
                    shape(tm1), tm1)
    gra = tm1; call print_tensor(shape(tm1), tm1, 'hnuc_g')
    ! overlap contribution
    DFD = D*F*D
    call rsp_ovlave(mol, 1, (/'GEO '/), (/1/), shape(tm1), DFD, &
                    tm1, (/(0d0,0d0)/), D); DFD=0
    gra = gra + tm1; call print_tensor(shape(tm1), tm1, '-SgDFD')
    ! 1-electron contribution
    call rsp_oneave(mol, 1, (/'GEO '/), (/1/), shape(tm1), D, tm1)
    gra = gra + tm1; call print_tensor(shape(tm1), tm1, 'HgD')
    ! 2-electron contribution
    call rsp_twoave(mol, 1, (/'GEO '/), (/1/), shape(tm1), D, D, tm1)
    gra = gra + tm1/2; call print_tensor(shape(tm1), tm1/2, 'Gg(D)D/2')
    ! Kohn-Sham exchange correlation average
    call rsp_excave(mol, 1, (/'GEO '/), (/1/), shape(tm1), 1, (/D/), tm1)
    gra = gra + tm1; call print_tensor(shape(tm1), tm1, 'Exc_g(D)')
    ! print
    call print_tensor(shape(gra), gra, 'gradient = Eg')
    ! print formatted to file gradient
    open (unit=iounit, file='gradient', status='replace', action='write')
    call print_tensor(shape(gra), gra, unit=iounit)
    close (iounit)
    ! test dipole gradient average
    call rsp_oneave(mol, 2, (/'GEO ','EL  '/), (/1,1/), shape(tmp(:,1,:)), &
                       D, tmp(:,1,:))
    call print_tensor(shape(tmp(:,1,:)), tmp(:,1,:), 'dpgave form rsp_oneave')
    !
    !------------------------------
    ! calculate perturbed integrals
    !do i = 1, size(Sg)
       call rsp_ovlint(mol, 1, (/'GEO '/), (/1/), shape(Sg), Sg)
       call rsp_oneint(mol, 1, (/'GEO '/), (/1/), shape(Fg), Fg)
       call rsp_twoint(mol, 1, (/'GEO '/), (/1/), shape(Fg), D, Fg)
    !end do
! density independent contribution to Ff
    call rsp_oneint(mol, 1, (/'EL  '/), (/1/), shape(Ff), Ff)
! solve equations for Df and construct Ff
    pol = 0
    do i = 1, size(Df)
       Df(i) = 0d0*D
       FDSf(1) = Ff(i)*D*S
       FDSf(1) = FDSf(1) - S*D*Ff(i)
       X(1) = 0*FDSf(1)
       call mat_alloc(X(1))
       call rsp_mosolver_exec(FDSf(1), (/0d0/), X)
       X(1)=-2d0*X(1); FDSf(1)=0
       Df(i) = Df(i) + X(1)*S*D
       Df(i) = Df(i) - D*S*X(1); X(1)=0
       ! Df contribution to Ff
       call rsp_twoint(mol, 0, nof, noc, noc, Df(i), Ff(i:i))
       !polarzability
       call rsp_oneave(mol, 1, (/'EL  '/), (/1/), shape(pol(:,i)), Df(i), pol(:,i))
    end do
    call print_tensor(shape(pol), -pol, 'polarizability')
    ! solve equations for Dg
    do i = 1, size(Dg)
       FDSg(1) = F*D*Sg(i)
       FDSg(1) = FDSg(1) - Sg(i)*D*F
       Dg(i) = -D*Sg(i)*D
       !ajt Sorry, doesn't work atm, while aliasing is only partially implemented
       !ajt X(1) = FDSg(1)
       !ajt FDSg(1) = X(1)*D*S
       !ajt FDSg(1) = FDSg(1) - S*D*X(1)
       !ajt FDSg(1) = FDSg(1) + Fg(i)
       FDSg(1) = FDSg(1)*D*S - S*D*FDSg(1) + Fg(i)
       !/ajt
       call rsp_twoint(mol, 0, nof, noc, noc, Dg(i), FDSg(1:1))
       !ajt Sorry, doesn't work atm, while aliasing is only partially implemented
       !ajt X(1) = FDSg(1)
       !ajt FDSg(1) = X(1)*D*S
       !ajt FDSg(1) = FDSg(1) - S*D*X(1)
       FDSg(1) = FDSg(1)*D*S - S*D*FDSg(1)
       !/ajt
       X(1) = 0*FDSg(1); call mat_alloc(X(1))
       call rsp_mosolver_exec(FDSg(1), (/0d0/), X); X(1)=-2d0*X(1); FDSg(1)=0
       ! ajt FIXME if I put these on the same line, I get something else
       Dg(i) = Dg(i) + X(1)*S*D
       Dg(i) = Dg(i) - D*S*X(1); X(1)=0
       call rsp_twoint(mol, 0, nof, noc, noc, Dg(i), Fg(i:i))
    end do
    !------------------------------------
    ! contract Hessian, nuclear repulsion
    call rsp_nucpot(mol, 2, (/'GEO ','GEO '/), (0d0,0d0)*(/0,0/), &
                    (/1,1/), shape(tm2), tm2)
    hes = tm2; call print_tensor(shape(tm2), tm2, 'hnuc_ab')
    ! overlap contribution
    DFD = D*F*D
    call rsp_ovlave(mol, 2, (/'GEO ','GEO '/), (/1,1/), shape(tm2), DFD, &
                    tm2); DFD=0 !, (0d0,0d0)*(/0,0/), D)
    hes = hes + tm2; call print_tensor(shape(tm2), tm2, '-SabDFD')
    !
    do i = 1, size(Dg)
       DFDg = Dg(i)*F*D + D*Fg(i)*D + D*F*Dg(i)
       call rsp_ovlave(mol, 1, (/'GEO '/), (/1/), shape(tm2(:,i)), DFDg, &
                       tm2(:,i)); DFDg=0 !, (0d0,0d0)*(/0,0/), D)
    end do
    hes = hes + tm2; call print_tensor(shape(tm2), tm2, '-SaDFDb')
    ! 1-electron contribution
    call rsp_oneave(mol, 2, (/'GEO ','GEO '/), (/1,1/), shape(tm2), D, tm2)
    hes = hes + tm2; call print_tensor(shape(tm2), tm2, 'HabD')
    !
    do i = 1, size(Dg)
       call rsp_oneave(mol, 1, (/'GEO '/), (/1/), shape(tm2(:,i)), Dg(i), tm2(:,i))
    end do
    hes = hes + tm2; call print_tensor(shape(tm2), tm2, 'HaDb')
    ! 2-electron contribution
    call rsp_twoave(mol, 2, (/'GEO ','GEO '/), (/1,1/), shape(tm2), D, D, tm2)
    hes = hes + tm2/2; call print_tensor(shape(tm2), tm2/2, 'Gab(D)D/2')
    !
    do i = 1, size(Dg)
       call rsp_twoave(mol, 1, (/'GEO '/), (/1/), shape(tm2(:,i)), D, Dg(i), tm2(:,i))
    end do
    hes = hes + tm2; call print_tensor(shape(tm2), tm2, 'Ga(D)Db')
    ! print
    call print_tensor(shape(hes), hes, 'hessian = Egg')
    ! print formatted to file hessian
    open (unit=iounit, file='hessian', status='replace', action='write')
    call print_tensor(shape(hes), hes, unit=iounit)
    close (iounit)
    !------------------------------------
    ! Dipole Hessian contributions start here
    !
    ! overlap contribution
    DFD = D*F*D
    !
    do i = 1, size(Df)
       DFDf = Df(i)*F*D + D*Ff(i)*D + D*F*Df(i)
       call rsp_ovlave(mol, 2, (/'GEO ','GEO '/), (/1,1/), shape(tmp(:,:,i)), &
                       DFDf, tmp(:,:,i))
       DFDf=0   !, (0d0,0d0)*(/0,0/), D)
    end do
    Eggf = tmp; call print_tensor(shape(tmp), tmp, '-SabDFDf')
    !
    do j = 1, size(Df)
       do i = 1, size(Dg)
          DFDgf = Df(j)*Fg(i)*D + Df(j)*F*Dg(i) + D*Ff(j)*Dg(i) &
                + Dg(i)*Ff(j)*D + Dg(i)*F*Df(j) + D*Fg(i)*Df(j)
          call rsp_ovlave(mol, 1, (/'GEO '/), (/1/), shape(tmp(:,i,j)), &
                          DFDgf, tmp(:,i,j))
          DFDgf=0 !, (0d0,0d0)*(/0/), D)
       end do
    end do
    Eggf = Eggf + tmp; call print_tensor(shape(tmp), tmp, '-SaDFDbf1')
    ! 1-electron contribution 
    ! fixme: finite field calculation in rsp_oneave
    call rsp_oneave(mol, 3, (/'GEO ','GEO ','EL  '/), (/1,1,1/), shape(tmp), D, tmp)
    Eggf = Eggf + tmp; call print_tensor(shape(tmp), tmp, 'HabfD')
    !call print_tensor(shape(tmp), tmp, 'HabfD')
    !
    do i = 1, size(Df)
       call rsp_oneave(mol, 2, (/'GEO ','GEO '/), (/1,1/), shape(tmp(:,:,i)), &
                       Df(i), tmp(:,:,i))
    end do
    Eggf = Eggf + tmp; call print_tensor(shape(tmp), tmp, 'HabDf')
    !
    do i = 1, size(Dg)
    ! fixme: take care of this case in rsp_oneave
       call rsp_oneave(mol, 2, (/'GEO ','EL  '/), (/1,1/), shape(tmp(:,i,:)), &
                       Dg(i), tmp(:,i,:))
    end do
    Eggf = Eggf + tmp; call print_tensor(shape(tmp), tmp, 'HafDb')
    ! other 1-electron contribution
    do i = 1, size(Dg)
       call rsp_oneave(mol, 2, (/'GEO ','EL  '/), (/1,1/), shape(tmp(i,:,:)), &
                       Dg(i), tmp(i,:,:))
    end do
    Eggf = Eggf + tmp; call print_tensor(shape(tmp), tmp, 'HbfDa')
    ! 2-electron contribution
    do i = 1, size(Df)
       call rsp_twoave(mol, 2, (/'GEO ','GEO '/), (/1,1/), shape(tmp(:,:,i)), &
                       D, Df(i), tmp(:,:,i))
    end do
    Eggf = Eggf + tmp; call print_tensor(shape(tmp), tmp, 'Gab(D)Df')
    !
    do j = 1, size(Df)
       do i = 1, size(Dg)
          call rsp_twoave(mol, 1, (/'GEO '/), (/1/), shape(tmp(:,i,j)), &
                          Dg(i), Df(j), tmp(:,i,j))
       end do
    end do
    Eggf = Eggf + tmp; call print_tensor(shape(tmp), tmp, 'Ga(Db)Df')
    !
    do j = 1, size(Df)
       do i = 1, size(Dg)
          call rsp_twoave(mol, 1, (/'GEO '/), (/1/), shape(tmp(i,:,j)), &
                          Dg(i), Df(j), tmp(i,:,j))
       end do
    end do
    Eggf = Eggf + tmp; call print_tensor(shape(tmp), tmp, 'Gb(Da)Df')
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
    integer     :: uni, siz, dec, cwid, i
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
