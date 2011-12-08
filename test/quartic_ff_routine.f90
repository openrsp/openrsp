
! MR: Working on quartic force field in this routine

  subroutine prop_test_quarticff(mol, ng, S, D, F)
    use dalton_ifc_ng, only: rsp_mosolver_exec
    type(rsp_cfg), intent(in) :: mol
    integer,       intent(in) :: ng
    type(matrix),  intent(in) :: S, D, F
    type(matrix) DFD, Sg(ng), Dg(ng), Fg(ng), DSDg, FDSg(1), DFDg
    type(matrix) X(1), FgDS(ng), DgSD(ng), FDSgg, DSDgg, DFDgg
    type(matrix) Dgg(ng,ng), Sgg(ng,ng), Fg(ng,ng) ! NOT INITIALIZED
    complex(8)   gra(ng), tm1(ng), hes(ng,ng), tm2(ng,ng)
    !complex(8)   cub(ng,ng,ng), tmp(ng,ng,ng)
    complex(8)   cub(ng,ng,ng), tm3(ng,ng,ng) ! NOT ALL INITIALIZED
    complex(8) 	 qua(ng,ng,ng,ng), tmp(ng,ng,ng,ng) ! NOT ALL INITIALIZED
    character(4) nof(0) !no-field
    integer      i, j, k, m, noc(0) !no-comp
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
       FDSg(1) = F*D*Sg(i)
       FDSg(1) = FDSg(1) - Sg(i)*D*F
       Dg(i) = -D*Sg(i)*D
       X(1) = FDSg(1)
       FDSg(1) = X(1)*D*S
       FDSg(1) = FDSg(1) - S*D*X(1)
       FDSg(1) = FDSg(1) + Fg(i)
       call rsp_twoint(mol, 0, nof, noc, noc, Dg(i), FDSg(1:1))
       X(1) = FDSg(1)
       FDSg(1) = X(1)*D*S
       FDSg(1) = FDSg(1) - S*D*X(1)
       call mat_init(X(1), Dg(i), zero=.true.)
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
! MR: BEGIN QUARTIC FORCE FIELD
! USES 1,2 RULE
! NOTE: THIS ROUTINE HAS NOT BEEN TESTED YET (11/2011)
! NOT ALL NEW VARIABLES HAVE BEEN DECLARED YET
! RECORD HERE: WHAT NEW VARIABLES ARE INTRODUCED?
! DFDggg2p, DSDggg2p, FDSggg2p
!
! RECORD HERE: WHAT NEW VARIABLES ARE NEEDED 
! AND MUST BE CALCULATED BEFORE THIS?
! Dgg, Fgg, Sgg

! Sg, Fg AND Dg COPIED FROM CUBIC FF ROUTINE ABOVE

    do i = 1, size(Sg)
! MR: WHAT ABOUT XC?
       call rsp_ovlint(mol, 1, (/'GEO '/), (/1/), shape(Sg), Sg)
       call rsp_oneint(mol, 1, (/'GEO '/), (/1/), shape(Fg), Fg)
       call rsp_twoint(mol, 1, (/'GEO '/), (/1/), shape(Fg), D, Fg)
    end do
    ! solve equations
    do i = 1, size(Dg)
       FDSg(1) = F*D*Sg(i)
       FDSg(1) = FDSg(1) - Sg(i)*D*F
       Dg(i) = -D*Sg(i)*D
       X(1) = FDSg(1)
       FDSg(1) = X(1)*D*S
       FDSg(1) = FDSg(1) - S*D*X(1)
       FDSg(1) = FDSg(1) + Fg(i)
       call rsp_twoint(mol, 0, nof, noc, noc, Dg(i), FDSg(1:1))
       X(1) = FDSg(1)
       FDSg(1) = X(1)*D*S
       FDSg(1) = FDSg(1) - S*D*X(1)
       call mat_init(X(1), Dg(i), zero=.true.)
       call rsp_mosolver_exec(FDSg(1), (/0d0/), X); X(1)=-2d0*X(1); FDSg(1)=0
       ! ajt FIXME if I put these on the same line, I get something else
       Dg(i) = Dg(i) + X(1)*S*D
       Dg(i) = Dg(i) - D*S*X(1); X(1)=0
       call rsp_twoint(mol, 0, nof, noc, noc, Dg(i), Fg(i:i))
    end do

! CALCULATE Sgg
! MR: IS THIS CORRECT?
    call rsp_ovlint(mol, 2, (/'GEO ','GEO '/), (/1,1/), shape(Sgg), Sgg)

! CALCULATE Fgg
! MR: WHAT ABOUT XC?
! MR: IS THIS CORRECT?
    call rsp_oneint(mol, 2, (/'GEO ','GEO '/), (/1,1/), shape(Fgg), Fgg)
    call rsp_twoint(mol, 2, (/'GEO ','GEO '/), (/1,1/), shape(Fgg), D, Fgg)  
    do i = 1, size(Dg)
      call rsp_twoint(mol, 1, (/'GEO '/), (/1/), shape(Fgg(:,i)), Dg(i), Fgg(:,i))
    end do



! REMAINS:

! CALCULATE Dgg




! BEGIN ADDING CONTRIBUTIONS TO QUARTIC FORCE FIELD TENSOR

! A, B, C, D PERTURBATIONS IN SAME BLOCK - NO "OTHER" BLOCKS

! NUCPOT BLOCK: DONE: YES

    call rsp_nucpot(mol, 4, (/'GEO ','GEO ','GEO ','GEO '/), (0d0,0d0)*(/0,0,0,0/), &
                    (/1,1,1,1/), shape(tmp), tmp)

    qua = tmp; call print_tensor(shape(tmp), tmp, 'Hnuc_abcd')
! OVL BLOCK: DONE: YES

    DFD = D*F*D
    call rsp_ovlave(mol, 4, (/'GEO ','GEO ','GEO ','GEO '/), (/1,1,1,1/), &
                    shape(tmp), DFD, tmp); DFD=0 !, (0d0,0d0)*(/0,0/), D)
    qua = qua + tmp; call print_tensor(shape(tmp), tmp, '-SabcdDFD')

! CONSIDER DOING ALL OF THE NEXT FOUR INSIDE THE SAME LOOP TO SAVE MATRIX MULTIPLICATIONS
! S3W1-TYPE CALLS
    do i = 1, size(Dg)
       DFDg = Dg(i)*F*D + D*Fg(i)*D + D*F*Dg(i) ! Wd
       call rsp_ovlave(mol, 3, (/'GEO ','GEO ','GEO '/), (/1,1,1/), shape(tmp(:,:,:,i)), &
                       DFDg, tmp(:,:,:,i)); DFDg=0 !, (0d0,0d0)*(/0,0/), D)
    end do
    qua = qua + tmp; call print_tensor(shape(tmp), tmp, '-SabcDFDd')
    do i = 1, size(Dg)
       DFDg = Dg(i)*F*D + D*Fg(i)*D + D*F*Dg(i) ! Wc
       call rsp_ovlave(mol, 3, (/'GEO ','GEO ','GEO '/), (/1,1,1/), shape(tmp(:,:,i,:)), &
                       DFDg, tmp(:,:,i,:)); DFDg=0 !, (0d0,0d0)*(/0,0/), D)
    end do
    qua = qua + tmp; call print_tensor(shape(tmp), tmp, '-SabdDFDc')
    do i = 1, size(Dg)
       DFDg = Dg(i)*F*D + D*Fg(i)*D + D*F*Dg(i) ! Wb
       call rsp_ovlave(mol, 3, (/'GEO ','GEO ','GEO '/), (/1,1,1/), shape(tmp(:,i,:,:)), &
                       DFDg, tmp(:,i,:,:)); DFDg=0 !, (0d0,0d0)*(/0,0/), D)
    end do
    qua = qua + tmp; call print_tensor(shape(tmp), tmp, '-SabdDFDb')
    do i = 1, size(Dg)
       DFDg = Dg(i)*F*D + D*Fg(i)*D + D*F*Dg(i) ! Wa
       call rsp_ovlave(mol, 3, (/'GEO ','GEO ','GEO '/), (/1,1,1/), shape(tmp(i,:,:,:)), &
                       DFDg, tmp(i,:,:,:)); DFDg=0 !, (0d0,0d0)*(/0,0/), D)
    end do
    qua = qua + tmp; call print_tensor(shape(tmp), tmp, '-SabdDFDa')
! CONSIDER DOING ALL OF THE NEXT THREE INSIDE THE SAME LOOP TO SAVE MATRIX MULTIPLICATIONS
! S2W2-TYPE CALLS
    do j = 1, size(Dg)
       do i = 1, size(Dg) ! UNSURE OF DOUBLE INDICES BELOW
          DFDgg = Dgg(i,j)*F*D + Dg(i)*Fg(j)*D + Dg(i)*F*Dg(j) + &
		  D*F*Dgg(i,j) + D*Fg(j)*Dg(i) + Dg(j)*F*Dg(i) + &
		  Dg(j)*Fg(i)*D + D*Fgg(i,j)*D + D*Fg(i)*Fg(j) ! Wcd
          call rsp_ovlave(mol, 2, (/'GEO ','GEO '/), (/1,1/), shape(tmp(:,:,i,j)), &
                          DFDgg, tmp(:,:,i,j)); DFDgg=0 !, (0d0,0d0)*(/0/), D)
       end do
    end do
    qua = qua + tmp; call print_tensor(shape(tmp), tmp, '-SabDFDcd')

    do j = 1, size(Dg)
       do i = 1, size(Dg) ! UNSURE OF DOUBLE INDICES BELOW
          DFDgg = Dgg(i,j)*F*D + Dg(i)*Fg(j)*D + Dg(i)*F*Dg(j) + &
		  D*F*Dgg(i,j) + D*Fg(j)*Dg(i) + Dg(j)*F*Dg(i) + &
		  Dg(j)*Fg(i)*D + D*Fgg(i,j)*D + D*Fg(i)*Fg(j) ! Wbc
          call rsp_ovlave(mol, 2, (/'GEO ','GEO '/), (/1,1/), shape(tmp(:,i,:,j)), &
                          DFDgg, tmp(:,i,:,j)); DFDgg=0 !, (0d0,0d0)*(/0/), D)
       end do
    end do
    qua = qua + tmp; call print_tensor(shape(tmp), tmp, '-SacDFDbd')



    do j = 1, size(Dg)
       do i = 1, size(Dg) ! UNSURE OF DOUBLE INDICES BELOW
          DFDgg = Dgg(i,j)*F*D + Dg(i)*Fg(j)*D + Dg(i)*F*Dg(j) + &
		  D*F*Dgg(i,j) + D*Fg(j)*Dg(i) + Dg(j)*F*Dg(i) + &
		  Dg(j)*Fg(i)*D + D*Fgg(i,j)*D + D*Fg(i)*Fg(j) ! Wbc
          call rsp_ovlave(mol, 2, (/'GEO ','GEO '/), (/1,1/), shape(tmp(:,i,j,:)), &
                          DFDgg, tmp(:,i,j,:)); DFDgg=0 !, (0d0,0d0)*(/0/), D)
       end do
    end do
    qua = qua + tmp; call print_tensor(shape(tmp), tmp, '-SadDFDbc')


! W2'-TYPE CALL
! MAYBE LOOK AT AGAIN

    do k = 1, size(Dg)
      do j = 1, size(Dg)
	do i = 1, size(Dg)
	    DFDggg2p = Dgg(i,j)*(Fg(k)*D+F*Dg(k))+(D*Fg(k)+Dg(k)*F)*Dgg(i,j) + &
		      Dgg(i,k)*(Fg(j)*D+F*Dg(j))+(D*Fg(j)+Dg(j)*F)*Dgg(i,k) + &
		      Dgg(j,k)*(Fg(i)*D+F*Dg(i))+(D*Fg(i)+Dg(i)*F)*Dgg(j,k) + &
		      Dg(i)*(Fgg(j,k)*D+Fg(j)*Dg(k))+(D*Fgg(j,k)+Dg(k)*Fg(j))*Dg(i) + &
		      Dg(j)*(Fgg(i,k)*D+Fg(k)*Dg(i))+(D*Fgg(i,k)+Dg(i)*Fg(k))*Dg(j) + &
		      Dg(k)*(Fgg(i,j)*D+Fg(i)*Dg(j))+(D*Fgg(i,j)+Dg(j)*Fg(i))*Dg(k)
	    call rsp_ovlave(mol, 1, (/'GEO '/), (/1/), shape(tmp(:,i,j,k)), &
			    DFDggg2p, tmp(:,i,j,k)); DFDggg2p=0 !, (0d0,0d0)*(/0/), D)
	end do
      end do
    end do
    qua = qua + tmp; call print_tensor(shape(tmp), tmp, '-SaDFDbcd2p')





! ONEEL BLOCK: DONE: YES
! E0-TYPE CALLS
    call rsp_oneave(mol, 4, (/'GEO ','GEO ','GEO ','GEO '/), (/1,1,1,1/), shape(tmp), D, tmp)
    qua = qua + tmp; call print_tensor(shape(tmp), tmp, 'HabcdD')

! E1-TYPE CALLS
! Dg CALLS
    do i = 1, size(Dg)
       call rsp_oneave(mol, 3, (/'GEO ','GEO ','GEO '/), (/1,1,1/), shape(tmp(:,:,:,i)), &
                       Dg(i), tmp(:,:,:,i))
    end do
    qua = qua + tmp; call print_tensor(shape(tmp), tmp, 'HabcDd')

    do i = 1, size(Dg)
       call rsp_oneave(mol, 3, (/'GEO ','GEO ','GEO '/), (/1,1,1/), shape(tmp(:,:,i,:)), &
                       Dg(i), tmp(:,:,i,:))
    end do
    qua = qua + tmp; call print_tensor(shape(tmp), tmp, 'HabdDc')

    do i = 1, size(Dg)
       call rsp_oneave(mol, 3, (/'GEO ','GEO ','GEO '/), (/1,1,1/), shape(tmp(:,i,:,:)), &
                       Dg(i), tmp(:,i,:,:))
    end do
    qua = qua + tmp; call print_tensor(shape(tmp), tmp, 'HacdDb')

    do i = 1, size(Dg)
       call rsp_oneave(mol, 3, (/'GEO ','GEO ','GEO '/), (/1,1,1/), shape(tmp(i,:,:,:)), &
                       Dg(i), tmp(i,:,:,:))
    end do
    qua = qua + tmp; call print_tensor(shape(tmp), tmp, 'HbcdDa')


! Dgg CALLS
    do i = 1, size(Dg)
      do j = 1, size(Dg)

	 call rsp_oneave(mol, 2, (/'GEO ','GEO '/), (/1,1/), shape(tmp(:,:,i,j)), &
                       Dgg(i,j), tmp(:,:,i,j))

      end do
    end do
    qua = qua + tmp; call print_tensor(shape(tmp), tmp, 'HabDcd')


    do i = 1, size(Dg)
      do j = 1, size(Dg)

	 call rsp_oneave(mol, 2, (/'GEO ','GEO '/), (/1,1/), shape(tmp(:,i,:,j)), &
                       Dgg(i,j), tmp(:,i,:,j))

      end do
    end do
    qua = qua + tmp; call print_tensor(shape(tmp), tmp, 'HacDbd')


    do i = 1, size(Dg)
      do j = 1, size(Dg)

	 call rsp_oneave(mol, 2, (/'GEO ','GEO '/), (/1,1/), shape(tmp(:,i,j,:)), &
                       Dgg(i,j), tmp(:,i,j,:))

      end do
    end do
    qua = qua + tmp; call print_tensor(shape(tmp), tmp, 'HadDbc')



! TWOEL BLOCK: DONE: YES

! E0-TYPE CALLS

    call rsp_twoave(mol, 4, (/'GEO ','GEO ','GEO ','GEO '/), (/1,1,1,1/), shape(tmp), D, D, tmp)
    qua = qua + tmp/2; call print_tensor(shape(tmp), tmp/2, 'Gabcd(D)D/2')


! E1-TYPE CALLS
! D-Dg CALLS
    do i = 1, size(Dg)
       call rsp_twoave(mol, 3, (/'GEO ','GEO ','GEO '/), (/1,1,1/), shape(tmp(:,:,:,i)), &
                       D, Dg(i), tmp(:,:,:,i))
    end do
    qua = qua + tmp; call print_tensor(shape(tmp), tmp, 'Gabc(D)Dd')

    do i = 1, size(Dg)
       call rsp_twoave(mol, 3, (/'GEO ','GEO ','GEO '/), (/1,1,1/), shape(tmp(:,:,i,:)), &
                       D, Dg(i), tmp(:,:,i,:))
    end do
    qua = qua + tmp; call print_tensor(shape(tmp), tmp, 'Gabd(D)Dc')

    do i = 1, size(Dg)
       call rsp_twoave(mol, 3, (/'GEO ','GEO ','GEO '/), (/1,1,1/), shape(tmp(:,i,:,:)), &
                       D, Dg(i), tmp(:,i,:,:))
    end do
    qua = qua + tmp; call print_tensor(shape(tmp), tmp, 'Gacd(D)Db')

    do i = 1, size(Dg)
       call rsp_twoave(mol, 3, (/'GEO ','GEO ','GEO '/), (/1,1,1/), shape(tmp(i,:,:,:)), &
                       D, Dg(i), tmp(i,:,:,:))
    end do
    qua = qua + tmp; call print_tensor(shape(tmp), tmp, 'Gbcd(D)Da')

! D-Dgg CALLS
    do j = 1, size(Dg)
       do i = 1, size(Dg)
          call rsp_twoave(mol, 2, (/'GEO ','GEO '/), (/1,1/), shape(tmp(:,:,i,j)), &
                          D, Dgg(i,j), tmp(:,:,i,j))
       end do
    end do
    qua = qua + tmp; call print_tensor(shape(tmp), tmp, 'Gab(D)Dcd')


    do j = 1, size(Dg)
       do i = 1, size(Dg)
          call rsp_twoave(mol, 2, (/'GEO ','GEO '/), (/1,1/), shape(tmp(:,i,:,j)), &
                          D, Dgg(i,j), tmp(:,i,:,j))
       end do
    end do
    qua = qua + tmp; call print_tensor(shape(tmp), tmp, 'Gac(D)Dbd')


    do j = 1, size(Dg)
       do i = 1, size(Dg)
          call rsp_twoave(mol, 2, (/'GEO ','GEO '/), (/1,1/), shape(tmp(:,i,j,:)), &
                          D, Dgg(i,j), tmp(:,i,j,:))
       end do
    end do
    qua = qua + tmp; call print_tensor(shape(tmp), tmp, 'Gad(D)Dbc')


! E2-TYPE CALLS

! Dg-Dg CALLS
    do j = 1, size(Dg)
       do i = 1, size(Dg)
          call rsp_twoave(mol, 2, (/'GEO ','GEO '/), (/1,1/), shape(tmp(:,:,i,j)), &
                          Dg(i), Dg(j), tmp(:,:,i,j))
       end do
    end do
    qua = qua + tmp; call print_tensor(shape(tmp), tmp, 'Gab(Dc)Dd')


    do j = 1, size(Dg)
       do i = 1, size(Dg)
          call rsp_twoave(mol, 2, (/'GEO ','GEO '/), (/1,1/), shape(tmp(:,i,:,j)), &
                          Dg(i), Dg(j), tmp(:,i,:,j))
       end do
    end do
    qua = qua + tmp; call print_tensor(shape(tmp), tmp, 'Gac(Db)Dd')


    do j = 1, size(Dg)
       do i = 1, size(Dg)
          call rsp_twoave(mol, 2, (/'GEO ','GEO '/), (/1,1/), shape(tmp(:,i,j,:)), &
                          Dg(i), Dg(j), tmp(:,i,j,:))
       end do
    end do
    qua = qua + tmp; call print_tensor(shape(tmp), tmp, 'Gad(Db)Dc')

! MR: SHOULD THIS BE HERE? PROBABLY: MAYBE LOOK AGAIN
    do j = 1, size(Dg)
       do i = 1, size(Dg)
          call rsp_twoave(mol, 2, (/'GEO ','GEO '/), (/1,1/), shape(tmp(i,:,:,j)), &
                          Dg(i), Dg(j), tmp(i,:,:,j))
       end do
    end do
    qua = qua + tmp; call print_tensor(shape(tmp), tmp, 'Gbc(Da)Dd')


    do j = 1, size(Dg)
       do i = 1, size(Dg)
          call rsp_twoave(mol, 2, (/'GEO ','GEO '/), (/1,1/), shape(tmp(i,:,j,:)), &
                          Dg(i), Dg(j), tmp(i,:,j,:))
       end do
    end do
    qua = qua + tmp; call print_tensor(shape(tmp), tmp, 'Gbd(Da)Dc')


    do j = 1, size(Dg)
       do i = 1, size(Dg)
          call rsp_twoave(mol, 2, (/'GEO ','GEO '/), (/1,1/), shape(tmp(i,j,:,:)), &
                          Dg(i), Dg(j), tmp(i,j,:,:))
       end do
    end do
    qua = qua + tmp; call print_tensor(shape(tmp), tmp, 'Gcd(Da)Db')

! Dgg-Dg CALLS

    do k = 1, size(Dg)
      do j = 1, size(Dg)
	do i = 1, size(Dg)
	    call rsp_twoave(mol, 1, (/'GEO '/), (/1/), shape(tmp(:,i,j,k)), &
                          Dgg(i,j), Dg(k), tmp(:,i,j,k))
	end do
      end do
    end do
    qua = qua + tmp; call print_tensor(shape(tmp), tmp, 'Ga(Dbc)Dd')


! MR: SOME TROUBLE WITH INDICES
! MR: ARE THE TMP INDICES AND Dgg, Dg INDICES CORRECT?

    do k = 1, size(Dg)
      do j = 1, size(Dg)
	do i = 1, size(Dg)
	    call rsp_twoave(mol, 1, (/'GEO '/), (/1/), shape(tmp(:,i,k,j)), &
                          Dgg(i,j), Dg(k), tmp(:,i,k,j))
	end do
      end do
    end do
    qua = qua + tmp; call print_tensor(shape(tmp), tmp, 'Ga(Dbd)Dc')

! MR: SOME CONTRIBUTIONS LACKING? WHY NOT Ga(Dcd)Db? SOLVED: Ga(Dcd)Db = Ga(Db)Dcd



! Dg-Dgg CALLS
! MR: SOME TROUBLE WITH INDICES?
! MR: THE : IS PLACED RIGHT
! MR: WHAT ABOUT THE TMP INDICES? THEY SEEM OK 
! MR: WHAT ABOUT THE Dg and Dgg INDICES? THEY SEEM OK
! MR: IN FACT, EVERYTHING LOOKS OK BUT KEEP COMMENT IN CASE THERE IS SOME ERROR
    do k = 1, size(Dg)
      do j = 1, size(Dg)
	do i = 1, size(Dg)
	    call rsp_twoave(mol, 1, (/'GEO '/), (/1/), shape(tmp(:,i,j,k)), &
                          Dg(i), Dgg(j,k), tmp(:,i,j,k))
	end do
      end do
    end do
    qua = qua + tmp; call print_tensor(shape(tmp), tmp, 'Ga(Db)Dcd')

    do k = 1, size(Dg)
      do j = 1, size(Dg)
	do i = 1, size(Dg)
	    call rsp_twoave(mol, 1, (/'GEO '/), (/1/), shape(tmp(i,:,j,k)), &
                          Dg(i), Dgg(j,k), tmp(i,:,j,k))
	end do
      end do
    end do
    qua = qua + tmp; call print_tensor(shape(tmp), tmp, 'Gb(Da)Dcd')

    do k = 1, size(Dg)
      do j = 1, size(Dg)
	do i = 1, size(Dg)
	    call rsp_twoave(mol, 1, (/'GEO '/), (/1/), shape(tmp(i,j,:,k)), &
                          Dg(i), Dgg(j,k), tmp(i,j,:,k))
	end do
      end do
    end do
    qua = qua + tmp; call print_tensor(shape(tmp), tmp, 'Gc(Da)Dbd')

    do k = 1, size(Dg)
      do j = 1, size(Dg)
	do i = 1, size(Dg)
	    call rsp_twoave(mol, 1, (/'GEO '/), (/1/), shape(tmp(i,j,k,:)), &
                          Dg(i), Dgg(j,k), tmp(i,j,k,:))
	end do
      end do
    end do
    qua = qua + tmp; call print_tensor(shape(tmp), tmp, 'Gd(Da)Dbc')





! XC BLOCK: DONE: YES

! MR: IS THERE ALWAYS A MULTIPLICATION WITH AN UNPERTURBED D FIRST?
! MR: ASSUME SO FOR NOW

! E0-TYPE CALLS

    call rsp_excave(mol, 4, (/'GEO ','GEO ','GEO ','GEO '/), (/1,1,1,1/), &
	shape(tmp(:,:,:,:)), 1, (/D/), tmp(:,:,:,:))
    qua = qua + tmp; call print_tensor(shape(tmp), tmp, 'Exc_abcd(D)')


! E1-TYPE CALLS
! D-Dg CALLS

    do i = 1, size(Dg)
      call rsp_excave(mol, 3, (/'GEO ','GEO ','GEO '/), (/1,1,1/), &
	  shape(tmp(:,:,:,i)), 2, (/D, Dg(i)/), tmp(:,:,:,i))
    end do
    qua = qua + tmp; call print_tensor(shape(tmp), tmp, 'Exc_abc(D)Dd')

    do i = 1, size(Dg)
      call rsp_excave(mol, 3, (/'GEO ','GEO ','GEO '/), (/1,1,1/), &
	  shape(tmp(:,:,i,:)), 2, (/D, Dg(i)/), tmp(:,:,i,:))
    end do
    qua = qua + tmp; call print_tensor(shape(tmp), tmp, 'Exc_abd(D)Dc')

    do i = 1, size(Dg)
      call rsp_excave(mol, 3, (/'GEO ','GEO ','GEO '/), (/1,1,1/), &
	  shape(tmp(:,i,:,:)), 2, (/D, Dg(i)/), tmp(:,i,:,:))
    end do
    qua = qua + tmp; call print_tensor(shape(tmp), tmp, 'Exc_acd(D)Db')

    do i = 1, size(Dg)
      call rsp_excave(mol, 3, (/'GEO ','GEO ','GEO '/), (/1,1,1/), &
	  shape(tmp(i,:,:,:)), 2, (/D, Dg(i)/), tmp(i,:,:,:))
    end do
    qua = qua + tmp; call print_tensor(shape(tmp), tmp, 'Exc_bcd(D)Da')

! D-Dgg CALLS

    do j = 1, size(Dg)
      do i = 1, size(Dg)
	call rsp_excave(mol, 2, (/'GEO ','GEO '/), (/1,1/), &
	    shape(tmp(:,:,i,j)), 2, (/D, Dgg(i,j)/), tmp(:,:,i,j))
      end do
    end do
    qua = qua + tmp; call print_tensor(shape(tmp), tmp, 'Exc_ab(D)Dcd')

    do j = 1, size(Dg)
      do i = 1, size(Dg)
	call rsp_excave(mol, 2, (/'GEO ','GEO '/), (/1,1/), &
	    shape(tmp(:,i,:,j)), 2, (/D, Dgg(i,j)/), tmp(:,i,:,j))
      end do
    end do
    qua = qua + tmp; call print_tensor(shape(tmp), tmp, 'Exc_ac(D)Dbd')

    do j = 1, size(Dg)
      do i = 1, size(Dg)
	call rsp_excave(mol, 2, (/'GEO ','GEO '/), (/1,1/), &
	    shape(tmp(:,i,j,:)), 2, (/D, Dgg(i,j)/), tmp(:,i,j,:))
      end do
    end do
    qua = qua + tmp; call print_tensor(shape(tmp), tmp, 'Exc_ad(D)Dbc')


! E2-TYPE CALLS
! D-Dg-Dg CALLS

    do j = 1, size(Dg)
      do i = 1, size(Dg)
	call rsp_excave(mol, 2, (/'GEO ','GEO '/), (/1,1/), &
	    shape(tmp(:,:,i,j)), 3, (/D, Dg(i), Dg(j)/), tmp(:,:,i,j))
      end do
    end do
    qua = qua + tmp; call print_tensor(shape(tmp), tmp, 'Exc_ab(D)DcDd')


    do j = 1, size(Dg)
      do i = 1, size(Dg)
	call rsp_excave(mol, 2, (/'GEO ','GEO '/), (/1,1/), &
	    shape(tmp(:,i,:,j)), 3, (/D, Dg(i), Dg(j)/), tmp(:,i,:,j))
      end do
    end do
    qua = qua + tmp; call print_tensor(shape(tmp), tmp, 'Exc_ac(D)DbDd')

    do j = 1, size(Dg)
      do i = 1, size(Dg)
	call rsp_excave(mol, 2, (/'GEO ','GEO '/), (/1,1/), &
	    shape(tmp(:,i,j,:)), 3, (/D, Dg(i), Dg(j)/), tmp(:,i,j,:))
      end do
    end do
    qua = qua + tmp; call print_tensor(shape(tmp), tmp, 'Exc_ad(D)DbDc')

    do j = 1, size(Dg)
      do i = 1, size(Dg)
	call rsp_excave(mol, 2, (/'GEO ','GEO '/), (/1,1/), &
	    shape(tmp(i,:,:,j)), 3, (/D, Dg(i), Dg(j)/), tmp(i,:,:,j))
      end do
    end do
    qua = qua + tmp; call print_tensor(shape(tmp), tmp, 'Exc_bc(D)DaDd')

    do j = 1, size(Dg)
      do i = 1, size(Dg)
	call rsp_excave(mol, 2, (/'GEO ','GEO '/), (/1,1/), &
	    shape(tmp(i,:,j,:)), 3, (/D, Dg(i), Dg(j)/), tmp(i,:,j,:))
      end do
    end do
    qua = qua + tmp; call print_tensor(shape(tmp), tmp, 'Exc_bd(D)DaDc')

    do j = 1, size(Dg)
      do i = 1, size(Dg)
	call rsp_excave(mol, 2, (/'GEO ','GEO '/), (/1,1/), &
	    shape(tmp(i,j,:,:)), 3, (/D, Dg(i), Dg(j)/), tmp(i,j,:,:))
      end do
    end do
    qua = qua + tmp; call print_tensor(shape(tmp), tmp, 'Exc_cd(D)DaDb')


! D-Dgg-Dg CALLS

    do k = 1, size(Dg)
      do j = 1, size(Dg)
	do i = 1, size(Dg)
	  call rsp_excave(mol, 1, (/'GEO '/), (/1/), &
	      shape(tmp(:,i,j,k)), 3, (/D, Dgg(i,j), Dg(k)/), tmp(:,i,j,k))
	end do
      end do
    end do
    qua = qua + tmp; call print_tensor(shape(tmp), tmp, 'Exc_a(D)DbcDd')


    do k = 1, size(Dg)
      do j = 1, size(Dg)
	do i = 1, size(Dg)
	  call rsp_excave(mol, 1, (/'GEO '/), (/1/), &
	      shape(tmp(:,i,k,j)), 3, (/D, Dgg(i,j), Dg(k)/), tmp(:,i,k,j))
	end do
      end do
    end do
    qua = qua + tmp; call print_tensor(shape(tmp), tmp, 'Exc_a(D)DbdDc')


! D-Dg-Dgg CALLS
! MAYBE LOOK AT AGAIN

    do k = 1, size(Dg)
      do j = 1, size(Dg)
	do i = 1, size(Dg)
	  call rsp_excave(mol, 1, (/'GEO '/), (/1/), &
	      shape(tmp(:,i,j,k)), 3, (/D, Dg(i), Dgg(j,k)/), tmp(:,i,j,k))
	end do
      end do
    end do
    qua = qua + tmp; call print_tensor(shape(tmp), tmp, 'Exc_a(D)DbDcd')


    do k = 1, size(Dg)
      do j = 1, size(Dg)
	do i = 1, size(Dg)
	  call rsp_excave(mol, 1, (/'GEO '/), (/1/), &
	      shape(tmp(i,:,j,k)), 3, (/D, Dg(i), Dgg(j,k)/), tmp(i,:,j,k))
	end do
      end do
    end do
    qua = qua + tmp; call print_tensor(shape(tmp), tmp, 'Exc_b(D)DaDcd')


    do k = 1, size(Dg)
      do j = 1, size(Dg)
	do i = 1, size(Dg)
	  call rsp_excave(mol, 1, (/'GEO '/), (/1/), &
	      shape(tmp(i,j,:,k)), 3, (/D, Dg(i), Dgg(j,k)/), tmp(i,j,:,k))
	end do
      end do
    end do
    qua = qua + tmp; call print_tensor(shape(tmp), tmp, 'Exc_c(D)DaDbd')

    do k = 1, size(Dg)
      do j = 1, size(Dg)
	do i = 1, size(Dg)
	  call rsp_excave(mol, 1, (/'GEO '/), (/1/), &
	      shape(tmp(i,j,k,:)), 3, (/D, Dg(i), Dgg(j,k)/), tmp(i,j,k,:))
	end do
      end do
    end do
    qua = qua + tmp; call print_tensor(shape(tmp), tmp, 'Exc_d(D)DaDbc')


! E3-TYPE CALLS

! D-Dg-Dg-Dg CALLS

    do k = 1, size(Dg)
      do j = 1, size(Dg)
	do i = 1, size(Dg)
	  call rsp_excave(mol, 1, (/'GEO '/), (/1/), &
	      shape(tmp(:,i,j,k)), 4, (/D, Dg(i), Dg(j), Dg(k)/), tmp(:,i,j,k))
	end do
      end do
    end do
    qua = qua + tmp; call print_tensor(shape(tmp), tmp, 'Exc_a(D)DbDcDd')

    do k = 1, size(Dg)
      do j = 1, size(Dg)
	do i = 1, size(Dg)
	  call rsp_excave(mol, 1, (/'GEO '/), (/1/), &
	      shape(tmp(i,:,j,k)), 4, (/D, Dg(i), Dg(j), Dg(k)/), tmp(i,:,j,k))
	end do
      end do
    end do
    qua = qua + tmp; call print_tensor(shape(tmp), tmp, 'Exc_b(D)DaDcDd')

    do k = 1, size(Dg)
      do j = 1, size(Dg)
	do i = 1, size(Dg)
	  call rsp_excave(mol, 1, (/'GEO '/), (/1/), &
	      shape(tmp(i,j,:,k)), 4, (/D, Dg(i), Dg(j), Dg(k)/), tmp(i,j,:,k))
	end do
      end do
    end do
    qua = qua + tmp; call print_tensor(shape(tmp), tmp, 'Exc_c(D)DaDbDd')

    do k = 1, size(Dg)
      do j = 1, size(Dg)
	do i = 1, size(Dg)
	  call rsp_excave(mol, 1, (/'GEO '/), (/1/), &
	      shape(tmp(i,j,k,:)), 4, (/D, Dg(i), Dg(j), Dg(k)/), tmp(i,j,k,:))
	end do
      end do
    end do
    qua = qua + tmp; call print_tensor(shape(tmp), tmp, 'Exc_d(D)DaDbDc')

! D-Dg-Dgg-Dg CALLS
! MAYBE LOOK AT AGAIN

! WILL THE EMPTY ARRAYS IN THE rsp_excave CALLS CAUSE PROBLEMS?

    do m = 1, size(Dg)
      do k = 1, size(Dg)
	do j = 1, size(Dg)
	  do i = 1, size(Dg)
	    call rsp_excave(mol, 0, (//), (//), &
		shape(tmp(i,j,k,m)), 4, (/D, Dg(i), Dgg(j,k), Dg(m)/), tmp(i,j,k,m))
	  end do
	end do
      end do
    end do
    qua = qua + tmp; call print_tensor(shape(tmp), tmp, 'Exc(D)DaDbcDd')


    do m = 1, size(Dg)
      do k = 1, size(Dg)
	do j = 1, size(Dg)
	  do i = 1, size(Dg)
	    call rsp_excave(mol, 0, (//), (//), &
		shape(tmp(i,j,m,k)), 4, (/D, Dg(i), Dgg(j,k), Dg(m)/), tmp(i,j,m,k))
	  end do
	end do
      end do
    end do
    qua = qua + tmp; call print_tensor(shape(tmp), tmp, 'Exc(D)DaDbdDc')

! D-Dg-Dg-Dgg CALLS

    do m = 1, size(Dg)
      do k = 1, size(Dg)
	do j = 1, size(Dg)
	  do i = 1, size(Dg)
	    call rsp_excave(mol, 0, (//), (//), &
		shape(tmp(i,j,k,m)), 4, (/D, Dg(i), Dg(j), Dgg(k,m)/), tmp(i,j,k,m))
	  end do
	end do
      end do
    end do
    qua = qua + tmp; call print_tensor(shape(tmp), tmp, 'Exc(D)DaDbDcd')

! E4-TYPE CALLS

! D-Dg-Dg-Dg-Dg CALLS

    do m = 1, size(Dg)
      do k = 1, size(Dg)
	do j = 1, size(Dg)
	  do i = 1, size(Dg)
	    call rsp_excave(mol, 0, (//), (//), &
		shape(tmp(i,j,k,m)), 5, (/D, Dg(i), Dg(j), Dg(k), Dg(m)/), tmp(i,j,k,m))
	  end do
	end do
      end do
    end do
    qua = qua + tmp; call print_tensor(shape(tmp), tmp, 'Exc(D)DaDbDcDd')


! IDEMPOTENCY BLOCK: DONE: YES

    ! 'Zeta g'
    do i = 1, size(Dg)
       FgDS(i) = Fg(i)*D*S + S*D*Fg(i) - Fg(i) - F*D*Sg(i) - Sg(i)*D*F
    end do

    ! 'Zeta g' WITH DSDggg2p
    do k = 1, size(Dg)
      do j = 1, size(Dg)
	do i = 1, size(Dg)
	  ! MR: TAKEN FROM DFDggg2p AND ADAPTED (SWITCHED FROM F TO S)
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
    end do; FgDS=0;
! MR: WHY END WITH SEMICOLON ON PREVIOUS LINE?
    qua = qua + tmp; call print_tensor(shape(tmp), tmp, '-FaDS DSDbcd2p')

! SELF-CONSISTENCY BLOCK: DONE: YES

    ! 'Lambda g'
    do i = 1, size(Dg)
       DgSD(i) = Dg(i)*S*D - D*S*Dg(i)
    end do


    
    do k = 1, size(Dg)
      do j = 1, size(Dg)
	do i = 1, size(Dg)
	  ! MR: MAYBE LOOK AT AGAIN
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
          end do; FDSggg2p=0
       end do
    end do; DgSD=0;
    qua = qua + tmp; call print_tensor(shape(tmp), tmp, '-DaSD FDSbcd2p')

! ALL CONTRIBUTIONS ARE NOW ADDED TO QUARTIC FORCE FIELD TENSOR


! PRINT OUTPUT: DONE: YES

! MR: NOTE INDICES OF OUTPUT (EVEN THOUGH THEY DON'T REALLY MATTER A LOT HERE IN Egggg)
    open (unit=iounit, file='quarticff', status='replace', action='write')
    write (iounit,*)
    do j = 1, size(qua,3)
      do i = 1, size(qua,4)
	call print_tensor(shape(cub(:,:,j,i)), cub(:,:,j,i), unit=iounit)
      end do
    end do
    close (iounit)





  end subroutine
