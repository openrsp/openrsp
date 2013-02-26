! commented out by Bin Gao, 2012-11-05
! radovan: probably because it made problems with MPI?
#define GRCONT_NOT_AVAILABLE

module interface_2el

   use matrix_defop
   use matrix_lowlevel, only: mat_init
   use interface_molecule
   use interface_basis
   use interface_scf
 
   implicit none
 
   public rsp_twoint
   public rsp_twoave
 
   private

contains

  !> Average 2-electron integrals perturbed by fields f over the
  !> product of density matrces D1 and D2
  subroutine rsp_twoave(nf, f, c, nc, D1, D2, propsize, ave)

    use eri_contractions, only: ctr_arg
    use eri_basis_loops,  only: unopt_geodiff_loop
    use interface_interest

    !> number of fields
    integer,              intent(in)  :: nf, propsize
    !> field labels in std order
    character(4),         intent(in)  :: f(nf)
    !> first and number of- components in each field
    integer,              intent(in)  :: c(nf), nc(nf)
    !> density matrix to average over
    type(matrix), target, intent(in)  :: D1, D2
    !> output average
    complex(8),           intent(out) :: ave(propsize)
    !----------------------------------------------
    real(8), allocatable              :: real_ave(:)
    real(8), pointer :: tmp(:,:,:,:) !scratch
    real(8), pointer :: tmp_5(:,:,:,:,:) !scratch
    real(8), pointer :: tmp_6(:,:,:,:,:,:) !scratch
    type(matrix)  A(1) !scratch matrices
    type(ctr_arg) arg(1)
    real(8)       r
    integer       h, i, j, k, l, m, n, p, q, ncor

    if (any(f == 'EL  ')) then
       ave = 0.0d0
       return
    end if

    if (nf==0) then

       ! contract second density to Fock matrix, then trace with first
       A(1) = 0*D1
       call mat_ensure_alloc(A(1), only_alloc=.true.)
       call interface_scf_get_g(D2, A(1)) !Coulomb and exchange
       ave(1) = trace(A(1),D1)





    else if (nf==1 .and. f(1)=='GEO ') then

#ifdef PRG_DALTON
       ncor = 3 * get_nr_atoms()
       allocate(tmp(ncor,1,1,1))
       tmp = 0.0
#ifdef GRCONT_NOT_AVAILABLE
       arg(1) = ctr_arg(1, -huge(1), ncor, D1, D2, &
                        rank_one_pointer(ncor, tmp(:,1,1,1)))
       call unopt_geodiff_loop(basis_large, &
                               basis_small, &
                               arg)
       tmp = 2.0d0*tmp
#else
       n = D1%nrow
       f77_memory(     :n*n)   = reshape(D1%elms,(/n*n/))
       f77_memory(n*n+1:n*n*2) = reshape(D2%elms,(/n*n/))
       call GRCONT(f77_memory(n*n*2+1:), size(f77_memory)-n*n*2, &
                   tmp(:,1,1,1), ncor, .true., .false., &
                   1, 0, .true., .false., f77_memory(:n*n*2), 2)
#endif
       ave(:nc(1)) = tmp(c(1):c(1)+nc(1)-1,1,1,1)
       deallocate(tmp)
#endif /* ifdef PRG_DALTON */

#ifdef PRG_DIRAC
       allocate(real_ave(size(ave)))
       real_ave = 0.0
       call interest_mpi_wake_up()
       call interest_get_int(D1%nrow, D1%elms, D2%elms, 1, 0, size(real_ave), real_ave)
       do i = 1, size(ave)
          ave(i) = 2.0d0*real_ave(i)
       end do
       deallocate(real_ave)
#endif /* ifdef PRG_DIRAC */






    else if (nf==2 .and. all(f==(/'GEO ','GEO '/))) then

       ncor = 3 * get_nr_atoms()
#ifdef PRG_DALTON
       allocate(tmp(ncor,ncor,1,1))
       tmp = 0.0
#ifdef GRCONT_NOT_AVAILABLE
       arg(1) = ctr_arg(2, -huge(1), ncor, D1, D2, &
                        rank_one_pointer(ncor**2, tmp(:,:,1,1)))
       call unopt_geodiff_loop(basis_large, &
                               basis_small, &
                               arg)
       ! symmetrize
       do j = 1, ncor
          do i = 1, j
             r = tmp(i, j, 1, 1) + tmp(j, i, 1, 1)
             tmp(i, j, 1, 1) = 2.0d0*r
             tmp(j, i, 1, 1) = 2.0d0*r
          end do
       end do
#else
       n = D1%nrow
       f77_memory(     :n*n)   = reshape(D1%elms,(/n*n/))
       f77_memory(n*n+1:n*n*2) = reshape(D2%elms,(/n*n/))
       call GRCONT(f77_memory(n*n*2+1:), size(f77_memory)-n*n*2, &
                   tmp(:,:,1,1), ncor**2, .true., .false., &
                   2, 0, .true., .false., f77_memory(:n*n*2), 2)
#endif
       h = 0
       do i = 1, ncor
          do j = i, ncor
             h = h + 1
             ave(h) = tmp(i,j,1,1)
          end do
       end do

       deallocate(tmp)
#endif  /* ifdef PRG_DALTON */

#ifdef PRG_DIRAC
       allocate(real_ave(ncor*ncor))
       real_ave = 0.0
       call interest_mpi_wake_up()
       call interest_get_int(D1%nrow, D1%elms, D2%elms, 2, 0, size(real_ave), real_ave)
       p = 0
       q = 0
       do i = 1, ncor
          do j = 1, ncor
             q = q + 1
             if (j >= i) then
                p = p + 1
                ave(p) = 2.0d0*real_ave(q)
             end if
          end do
       end do
       deallocate(real_ave)
#endif /* ifdef PRG_DIRAC */




    else if (nf==3 .and. all(f==(/'GEO ','GEO ','GEO '/))) then

#ifdef PRG_DIRAC
       print *, 'error: twoave contribution not programmed'
       stop 1
#endif
       ! contract FULL cubic in tmp, unsymmetrized divided by six
       ncor = 3 * get_nr_atoms()
       allocate(tmp(ncor,ncor,ncor,1))
       tmp = 0.0
       arg(1) = ctr_arg(3, -huge(1), ncor, D1, D2, &
                        rank_one_pointer(ncor**3, tmp(:,:,:,1)))
       call unopt_geodiff_loop(basis_large, &
                               basis_small, &
                               arg)
       ! symmetrize
       do k = 1, ncor
          do j = 1, k
             do i = 1, j
                r = tmp(i,j,k,1) + tmp(i,k,j,1) + tmp(k,i,j,1) &
                  + tmp(k,j,i,1) + tmp(j,k,i,1) + tmp(j,i,k,1)
                tmp(i,j,k,1) = r
             end do
          end do
       end do

       h = 0
       do i = 1, ncor
          do j = i, ncor
             do k = j, ncor
                h = h + 1
                ! MR: UNSURE ABOUT ORIGINS OF FACTOR 2 IN NEXT LINE
                ave(h) = 2 * tmp(i,j,k,1)
             end do
          end do
       end do

       deallocate(tmp)






    else if (nf==4 .and. all(f==(/'GEO ','GEO ','GEO ','GEO '/))) then
 
#ifdef PRG_DIRAC
       print *, 'error: twoave contribution not programmed'
       stop 1
#endif
      ncor = 3 * get_nr_atoms()
       allocate(tmp(ncor,ncor,ncor,ncor))
       tmp = 0.0
       ! contract FULL quartic in tmp, unsymmetrized divided by 24
       arg(1) = ctr_arg(4, -huge(1), ncor, D1, D2, &
                        rank_one_pointer(ncor**4, tmp))
       call unopt_geodiff_loop(basis_large, &
                               basis_small, &
                               arg)
       ! symmetrize
       do l = 1, ncor
          do k = 1, l
             do j = 1, k
                do i = 1, j
                   r = tmp(i,j,k,l) + tmp(i,k,j,l) + tmp(k,i,j,l) &
                     + tmp(k,j,i,l) + tmp(j,k,i,l) + tmp(j,i,k,l) &
                     + tmp(i,j,l,k) + tmp(i,k,l,j) + tmp(k,i,l,j) &
                     + tmp(k,j,l,i) + tmp(j,k,l,i) + tmp(j,i,l,k) &
                     + tmp(i,l,j,k) + tmp(i,l,k,j) + tmp(k,l,i,j) &
                     + tmp(k,l,j,i) + tmp(j,l,k,i) + tmp(j,l,i,k) &
                     + tmp(l,i,j,k) + tmp(l,i,k,j) + tmp(l,k,i,j) &
                     + tmp(l,k,j,i) + tmp(l,j,k,i) + tmp(l,j,i,k)
                   tmp(i,j,k,l) = r
                end do
             end do
          end do
       end do

       h = 0
       do i = 1, ncor
          do j = i, ncor
             do k = j, ncor
                do m = k, ncor
                   h = h + 1
                   ! MR: UNSURE ABOUT ORIGINS OF FACTOR 2 IN NEXT LINE
                   ave(h) = 2 * tmp(i,j,k,m)
                end do
             end do
          end do
       end do

       deallocate(tmp)






    else if (nf==5 .and. all(f==(/'GEO ','GEO ','GEO ','GEO ', 'GEO '/))) then

#ifdef PRG_DIRAC
       print *, 'error: twoave contribution not programmed'
       stop 1
#endif
       ncor = 3 * get_nr_atoms()
       allocate(tmp_5(ncor,ncor,ncor,ncor,ncor))
       tmp_5 = 0.0
       ! contract FULL quartic in tmp, unsymmetrized divided by 24
       arg(1) = ctr_arg(5, -huge(1), ncor, D1, D2, &
                        rank_one_pointer(ncor**5, tmp_5))
       call unopt_geodiff_loop(basis_large, &
                               basis_small, &
                               arg)
       ! symmetrize
       do m = 1, ncor
          do l = 1, m
             do k = 1, l
                do j = 1, k
                   do i = 1, j

                      r = tmp_5(i,j,k,l,m) + tmp_5(i,j,k,m,l) + tmp_5(i,j,l,k,m) + &
                          tmp_5(i,j,l,m,k) + tmp_5(i,j,m,k,l) + tmp_5(i,j,m,l,k) + &
                          tmp_5(i,k,j,l,m) + tmp_5(i,k,j,m,l) + tmp_5(i,k,l,j,m) + &
                          tmp_5(i,k,l,m,j) + tmp_5(i,k,m,j,l) + tmp_5(i,k,m,l,j) + &
                          tmp_5(i,l,j,k,m) + tmp_5(i,l,j,m,k) + tmp_5(i,l,k,j,m) + &
                          tmp_5(i,l,k,m,j) + tmp_5(i,l,m,j,k) + tmp_5(i,l,m,k,j) + &
                          tmp_5(i,m,j,k,l) + tmp_5(i,m,j,l,k) + tmp_5(i,m,k,j,l) + &
                          tmp_5(i,m,k,l,j) + tmp_5(i,m,l,j,k) + tmp_5(i,m,l,k,j) + &
                          tmp_5(j,i,k,l,m) + tmp_5(j,i,k,m,l) + tmp_5(j,i,l,k,m) + &
                          tmp_5(j,i,l,m,k) + tmp_5(j,i,m,k,l) + tmp_5(j,i,m,l,k) + &
                          tmp_5(j,k,i,l,m) + tmp_5(j,k,i,m,l) + tmp_5(j,k,l,i,m) + &
                          tmp_5(j,k,l,m,i) + tmp_5(j,k,m,i,l) + tmp_5(j,k,m,l,i) + &
                          tmp_5(j,l,i,k,m) + tmp_5(j,l,i,m,k) + tmp_5(j,l,k,i,m) + &
                          tmp_5(j,l,k,m,i) + tmp_5(j,l,m,i,k) + tmp_5(j,l,m,k,i) + &
                          tmp_5(j,m,i,k,l) + tmp_5(j,m,i,l,k) + tmp_5(j,m,k,i,l) + &
                          tmp_5(j,m,k,l,i) + tmp_5(j,m,l,i,k) + tmp_5(j,m,l,k,i) + &
                          tmp_5(k,i,j,l,m) + tmp_5(k,i,j,m,l) + tmp_5(k,i,l,j,m) + &
                          tmp_5(k,i,l,m,j) + tmp_5(k,i,m,j,l) + tmp_5(k,i,m,l,j) + &
                          tmp_5(k,j,i,l,m) + tmp_5(k,j,i,m,l) + tmp_5(k,j,l,i,m) + &
                          tmp_5(k,j,l,m,i) + tmp_5(k,j,m,i,l) + tmp_5(k,j,m,l,i) + &
                          tmp_5(k,l,i,j,m) + tmp_5(k,l,i,m,j) + tmp_5(k,l,j,i,m) + &
                          tmp_5(k,l,j,m,i) + tmp_5(k,l,m,i,j) + tmp_5(k,l,m,j,i) + &
                          tmp_5(k,m,i,j,l) + tmp_5(k,m,i,l,j) + tmp_5(k,m,j,i,l) + &
                          tmp_5(k,m,j,l,i) + tmp_5(k,m,l,i,j) + tmp_5(k,m,l,j,i) + &
                          tmp_5(l,i,j,k,m) + tmp_5(l,i,j,m,k) + tmp_5(l,i,k,j,m) + &
                          tmp_5(l,i,k,m,j) + tmp_5(l,i,m,j,k) + tmp_5(l,i,m,k,j) + &
                          tmp_5(l,j,i,k,m) + tmp_5(l,j,i,m,k) + tmp_5(l,j,k,i,m) + &
                          tmp_5(l,j,k,m,i) + tmp_5(l,j,m,i,k) + tmp_5(l,j,m,k,i) + &
                          tmp_5(l,k,i,j,m) + tmp_5(l,k,i,m,j) + tmp_5(l,k,j,i,m) + &
                          tmp_5(l,k,j,m,i) + tmp_5(l,k,m,i,j) + tmp_5(l,k,m,j,i) + &
                          tmp_5(l,m,i,j,k) + tmp_5(l,m,i,k,j) + tmp_5(l,m,j,i,k) + &
                          tmp_5(l,m,j,k,i) + tmp_5(l,m,k,i,j) + tmp_5(l,m,k,j,i) + &
                          tmp_5(m,i,j,k,l) + tmp_5(m,i,j,l,k) + tmp_5(m,i,k,j,l) + &
                          tmp_5(m,i,k,l,j) + tmp_5(m,i,l,j,k) + tmp_5(m,i,l,k,j) + &
                          tmp_5(m,j,i,k,l) + tmp_5(m,j,i,l,k) + tmp_5(m,j,k,i,l) + &
                          tmp_5(m,j,k,l,i) + tmp_5(m,j,l,i,k) + tmp_5(m,j,l,k,i) + &
                          tmp_5(m,k,i,j,l) + tmp_5(m,k,i,l,j) + tmp_5(m,k,j,i,l) + &
                          tmp_5(m,k,j,l,i) + tmp_5(m,k,l,i,j) + tmp_5(m,k,l,j,i) + &
                          tmp_5(m,l,i,j,k) + tmp_5(m,l,i,k,j) + tmp_5(m,l,j,i,k) + &
                          tmp_5(m,l,j,k,i) + tmp_5(m,l,k,i,j) + tmp_5(m,l,k,j,i) 
                          tmp_5(i,j,k,l,m) = r

                   end do
                end do
             end do
          end do
       end do

       h = 0
       do i = 1, ncor
          do j = i, ncor
             do k = j, ncor
                do l = k, ncor
                   do m = l, ncor
                      h = h + 1
                      ! MR: UNSURE ABOUT ORIGINS OF FACTOR 2 IN NEXT LINE
                      ave(h) = 2 * tmp_5(i,j,k,l,m)
                   end do
                end do
             end do
          end do
       end do

       deallocate(tmp_5)






    else if (nf==6 .and. all(f==(/'GEO ','GEO ','GEO ','GEO ', 'GEO ', 'GEO '/))) then

#ifdef PRG_DIRAC
       print *, 'error: twoave contribution not programmed'
       stop 1
#endif
       ncor = 3 * get_nr_atoms()
       allocate(tmp_6(ncor,ncor,ncor,ncor,ncor,ncor))
       tmp_6 = 0.0
       ! contract FULL quartic in tmp, unsymmetrized divided by 24
       arg(1) = ctr_arg(6, -huge(1), ncor, D1, D2, &
                        rank_one_pointer(ncor**6, tmp_6))
       call unopt_geodiff_loop(basis_large, &
                               basis_small, &
                               arg)
       ! symmetrize
       do n = 1, ncor
          do m = 1, n
             do l = 1, m
                do k = 1, l
                   do j = 1, k
                      do i = 1, j

                         r = tmp_6(i,j,k,l,m,n) + tmp_6(i,j,k,l,n,m) + tmp_6(i,j,k,m,l,n) + &
                             tmp_6(i,j,k,m,n,l) + tmp_6(i,j,k,n,l,m) + tmp_6(i,j,k,n,m,l) + &
                             tmp_6(i,j,l,k,m,n) + tmp_6(i,j,l,k,n,m) + tmp_6(i,j,l,m,k,n) + &
                             tmp_6(i,j,l,m,n,k) + tmp_6(i,j,l,n,k,m) + tmp_6(i,j,l,n,m,k) + &
                             tmp_6(i,j,m,k,l,n) + tmp_6(i,j,m,k,n,l) + tmp_6(i,j,m,l,k,n) + &
                             tmp_6(i,j,m,l,n,k) + tmp_6(i,j,m,n,k,l) + tmp_6(i,j,m,n,l,k) + &
                             tmp_6(i,j,n,k,l,m) + tmp_6(i,j,n,k,m,l) + tmp_6(i,j,n,l,k,m) + &
                             tmp_6(i,j,n,l,m,k) + tmp_6(i,j,n,m,k,l) + tmp_6(i,j,n,m,l,k) + &
                             tmp_6(i,k,j,l,m,n) + tmp_6(i,k,j,l,n,m) + tmp_6(i,k,j,m,l,n) + &
                             tmp_6(i,k,j,m,n,l) + tmp_6(i,k,j,n,l,m) + tmp_6(i,k,j,n,m,l) + &
                             tmp_6(i,k,l,j,m,n) + tmp_6(i,k,l,j,n,m) + tmp_6(i,k,l,m,j,n) + &
                             tmp_6(i,k,l,m,n,j) + tmp_6(i,k,l,n,j,m) + tmp_6(i,k,l,n,m,j) + &
                             tmp_6(i,k,m,j,l,n) + tmp_6(i,k,m,j,n,l) + tmp_6(i,k,m,l,j,n) + &
                             tmp_6(i,k,m,l,n,j) + tmp_6(i,k,m,n,j,l) + tmp_6(i,k,m,n,l,j) + &
                             tmp_6(i,k,n,j,l,m) + tmp_6(i,k,n,j,m,l) + tmp_6(i,k,n,l,j,m) + &
                             tmp_6(i,k,n,l,m,j) + tmp_6(i,k,n,m,j,l) + tmp_6(i,k,n,m,l,j) + &
                             tmp_6(i,l,j,k,m,n) + tmp_6(i,l,j,k,n,m) + tmp_6(i,l,j,m,k,n) + &
                             tmp_6(i,l,j,m,n,k) + tmp_6(i,l,j,n,k,m) + tmp_6(i,l,j,n,m,k) + &
                             tmp_6(i,l,k,j,m,n) + tmp_6(i,l,k,j,n,m) + tmp_6(i,l,k,m,j,n) + &
                             tmp_6(i,l,k,m,n,j) + tmp_6(i,l,k,n,j,m) + tmp_6(i,l,k,n,m,j) + &
                             tmp_6(i,l,m,j,k,n) + tmp_6(i,l,m,j,n,k) + tmp_6(i,l,m,k,j,n) + &
                             tmp_6(i,l,m,k,n,j) + tmp_6(i,l,m,n,j,k) + tmp_6(i,l,m,n,k,j) + &
                             tmp_6(i,l,n,j,k,m) + tmp_6(i,l,n,j,m,k) + tmp_6(i,l,n,k,j,m) + &
                             tmp_6(i,l,n,k,m,j) + tmp_6(i,l,n,m,j,k) + tmp_6(i,l,n,m,k,j) + &
                             tmp_6(i,m,j,k,l,n) + tmp_6(i,m,j,k,n,l) + tmp_6(i,m,j,l,k,n) + &
                             tmp_6(i,m,j,l,n,k) + tmp_6(i,m,j,n,k,l) + tmp_6(i,m,j,n,l,k) + &
                             tmp_6(i,m,k,j,l,n) + tmp_6(i,m,k,j,n,l) + tmp_6(i,m,k,l,j,n) + &
                             tmp_6(i,m,k,l,n,j) + tmp_6(i,m,k,n,j,l) + tmp_6(i,m,k,n,l,j) + &
                             tmp_6(i,m,l,j,k,n) + tmp_6(i,m,l,j,n,k) + tmp_6(i,m,l,k,j,n) + &
                             tmp_6(i,m,l,k,n,j) + tmp_6(i,m,l,n,j,k) + tmp_6(i,m,l,n,k,j) + &
                             tmp_6(i,m,n,j,k,l) + tmp_6(i,m,n,j,l,k) + tmp_6(i,m,n,k,j,l) + &
                             tmp_6(i,m,n,k,l,j) + tmp_6(i,m,n,l,j,k) + tmp_6(i,m,n,l,k,j) + &
                             tmp_6(i,n,j,k,l,m) + tmp_6(i,n,j,k,m,l) + tmp_6(i,n,j,l,k,m) + &
                             tmp_6(i,n,j,l,m,k) + tmp_6(i,n,j,m,k,l) + tmp_6(i,n,j,m,l,k) + &
                             tmp_6(i,n,k,j,l,m) + tmp_6(i,n,k,j,m,l) + tmp_6(i,n,k,l,j,m) + &
                             tmp_6(i,n,k,l,m,j) + tmp_6(i,n,k,m,j,l) + tmp_6(i,n,k,m,l,j) + &
                             tmp_6(i,n,l,j,k,m) + tmp_6(i,n,l,j,m,k) + tmp_6(i,n,l,k,j,m) + &
                             tmp_6(i,n,l,k,m,j) + tmp_6(i,n,l,m,j,k) + tmp_6(i,n,l,m,k,j) + &
                             tmp_6(i,n,m,j,k,l) + tmp_6(i,n,m,j,l,k) + tmp_6(i,n,m,k,j,l) + &
                             tmp_6(i,n,m,k,l,j) + tmp_6(i,n,m,l,j,k) + tmp_6(i,n,m,l,k,j) + &
                             tmp_6(j,i,k,l,m,n) + tmp_6(j,i,k,l,n,m) + tmp_6(j,i,k,m,l,n) + &
                             tmp_6(j,i,k,m,n,l) + tmp_6(j,i,k,n,l,m) + tmp_6(j,i,k,n,m,l) + &
                             tmp_6(j,i,l,k,m,n) + tmp_6(j,i,l,k,n,m) + tmp_6(j,i,l,m,k,n) + &
                             tmp_6(j,i,l,m,n,k) + tmp_6(j,i,l,n,k,m) + tmp_6(j,i,l,n,m,k) + &
                             tmp_6(j,i,m,k,l,n) + tmp_6(j,i,m,k,n,l) + tmp_6(j,i,m,l,k,n) + &
                             tmp_6(j,i,m,l,n,k) + tmp_6(j,i,m,n,k,l) + tmp_6(j,i,m,n,l,k) + &
                             tmp_6(j,i,n,k,l,m) + tmp_6(j,i,n,k,m,l) + tmp_6(j,i,n,l,k,m) + &
                             tmp_6(j,i,n,l,m,k) + tmp_6(j,i,n,m,k,l) + tmp_6(j,i,n,m,l,k) + &
                             tmp_6(j,k,i,l,m,n) + tmp_6(j,k,i,l,n,m) + tmp_6(j,k,i,m,l,n) + &
                             tmp_6(j,k,i,m,n,l) + tmp_6(j,k,i,n,l,m) + tmp_6(j,k,i,n,m,l) + &
                             tmp_6(j,k,l,i,m,n) + tmp_6(j,k,l,i,n,m) + tmp_6(j,k,l,m,i,n) + &
                             tmp_6(j,k,l,m,n,i) + tmp_6(j,k,l,n,i,m) + tmp_6(j,k,l,n,m,i) + &
                             tmp_6(j,k,m,i,l,n) + tmp_6(j,k,m,i,n,l) + tmp_6(j,k,m,l,i,n) + &
                             tmp_6(j,k,m,l,n,i) + tmp_6(j,k,m,n,i,l) + tmp_6(j,k,m,n,l,i) + &
                             tmp_6(j,k,n,i,l,m) + tmp_6(j,k,n,i,m,l) + tmp_6(j,k,n,l,i,m) + &
                             tmp_6(j,k,n,l,m,i) + tmp_6(j,k,n,m,i,l) + tmp_6(j,k,n,m,l,i) + &
                             tmp_6(j,l,i,k,m,n) + tmp_6(j,l,i,k,n,m) + tmp_6(j,l,i,m,k,n) + &
                             tmp_6(j,l,i,m,n,k) + tmp_6(j,l,i,n,k,m) + tmp_6(j,l,i,n,m,k) + &
                             tmp_6(j,l,k,i,m,n) + tmp_6(j,l,k,i,n,m) + tmp_6(j,l,k,m,i,n) + &
                             tmp_6(j,l,k,m,n,i) + tmp_6(j,l,k,n,i,m) + tmp_6(j,l,k,n,m,i) + &
                             tmp_6(j,l,m,i,k,n) + tmp_6(j,l,m,i,n,k) + tmp_6(j,l,m,k,i,n) + &
                             tmp_6(j,l,m,k,n,i) + tmp_6(j,l,m,n,i,k) + tmp_6(j,l,m,n,k,i) + &
                             tmp_6(j,l,n,i,k,m) + tmp_6(j,l,n,i,m,k) + tmp_6(j,l,n,k,i,m) + &
                             tmp_6(j,l,n,k,m,i) + tmp_6(j,l,n,m,i,k) + tmp_6(j,l,n,m,k,i) + &
                             tmp_6(j,m,i,k,l,n) + tmp_6(j,m,i,k,n,l) + tmp_6(j,m,i,l,k,n) + &
                             tmp_6(j,m,i,l,n,k) + tmp_6(j,m,i,n,k,l) + tmp_6(j,m,i,n,l,k) + &
                             tmp_6(j,m,k,i,l,n) + tmp_6(j,m,k,i,n,l) + tmp_6(j,m,k,l,i,n) + &
                             tmp_6(j,m,k,l,n,i) + tmp_6(j,m,k,n,i,l) + tmp_6(j,m,k,n,l,i) + &
                             tmp_6(j,m,l,i,k,n) + tmp_6(j,m,l,i,n,k) + tmp_6(j,m,l,k,i,n) + &
                             tmp_6(j,m,l,k,n,i) + tmp_6(j,m,l,n,i,k) + tmp_6(j,m,l,n,k,i) + &
                             tmp_6(j,m,n,i,k,l) + tmp_6(j,m,n,i,l,k) + tmp_6(j,m,n,k,i,l) + &
                             tmp_6(j,m,n,k,l,i) + tmp_6(j,m,n,l,i,k) + tmp_6(j,m,n,l,k,i) + &
                             tmp_6(j,n,i,k,l,m) + tmp_6(j,n,i,k,m,l) + tmp_6(j,n,i,l,k,m) + &
                             tmp_6(j,n,i,l,m,k) + tmp_6(j,n,i,m,k,l) + tmp_6(j,n,i,m,l,k) + &
                             tmp_6(j,n,k,i,l,m) + tmp_6(j,n,k,i,m,l) + tmp_6(j,n,k,l,i,m) + &
                             tmp_6(j,n,k,l,m,i) + tmp_6(j,n,k,m,i,l) + tmp_6(j,n,k,m,l,i) + &
                             tmp_6(j,n,l,i,k,m) + tmp_6(j,n,l,i,m,k) + tmp_6(j,n,l,k,i,m) + &
                             tmp_6(j,n,l,k,m,i) + tmp_6(j,n,l,m,i,k) + tmp_6(j,n,l,m,k,i) + &
                             tmp_6(j,n,m,i,k,l) + tmp_6(j,n,m,i,l,k) + tmp_6(j,n,m,k,i,l) + &
                             tmp_6(j,n,m,k,l,i) + tmp_6(j,n,m,l,i,k) + tmp_6(j,n,m,l,k,i) + &
                             tmp_6(k,i,j,l,m,n) + tmp_6(k,i,j,l,n,m) + tmp_6(k,i,j,m,l,n) + &
                             tmp_6(k,i,j,m,n,l) + tmp_6(k,i,j,n,l,m) + tmp_6(k,i,j,n,m,l) + &
                             tmp_6(k,i,l,j,m,n) + tmp_6(k,i,l,j,n,m) + tmp_6(k,i,l,m,j,n) + &
                             tmp_6(k,i,l,m,n,j) + tmp_6(k,i,l,n,j,m) + tmp_6(k,i,l,n,m,j) + &
                             tmp_6(k,i,m,j,l,n) + tmp_6(k,i,m,j,n,l) + tmp_6(k,i,m,l,j,n) + &
                             tmp_6(k,i,m,l,n,j) + tmp_6(k,i,m,n,j,l) + tmp_6(k,i,m,n,l,j) + &
                             tmp_6(k,i,n,j,l,m) + tmp_6(k,i,n,j,m,l) + tmp_6(k,i,n,l,j,m) + &
                             tmp_6(k,i,n,l,m,j) + tmp_6(k,i,n,m,j,l) + tmp_6(k,i,n,m,l,j) + &
                             tmp_6(k,j,i,l,m,n) + tmp_6(k,j,i,l,n,m) + tmp_6(k,j,i,m,l,n) + &
                             tmp_6(k,j,i,m,n,l) + tmp_6(k,j,i,n,l,m) + tmp_6(k,j,i,n,m,l) + &
                             tmp_6(k,j,l,i,m,n) + tmp_6(k,j,l,i,n,m) + tmp_6(k,j,l,m,i,n) + &
                             tmp_6(k,j,l,m,n,i) + tmp_6(k,j,l,n,i,m) + tmp_6(k,j,l,n,m,i) + &
                             tmp_6(k,j,m,i,l,n) + tmp_6(k,j,m,i,n,l) + tmp_6(k,j,m,l,i,n) + &
                             tmp_6(k,j,m,l,n,i) + tmp_6(k,j,m,n,i,l) + tmp_6(k,j,m,n,l,i) + &
                             tmp_6(k,j,n,i,l,m) + tmp_6(k,j,n,i,m,l) + tmp_6(k,j,n,l,i,m) + &
                             tmp_6(k,j,n,l,m,i) + tmp_6(k,j,n,m,i,l) + tmp_6(k,j,n,m,l,i) + &
                             tmp_6(k,l,i,j,m,n) + tmp_6(k,l,i,j,n,m) + tmp_6(k,l,i,m,j,n) + &
                             tmp_6(k,l,i,m,n,j) + tmp_6(k,l,i,n,j,m) + tmp_6(k,l,i,n,m,j) + &
                             tmp_6(k,l,j,i,m,n) + tmp_6(k,l,j,i,n,m) + tmp_6(k,l,j,m,i,n) + &
                             tmp_6(k,l,j,m,n,i) + tmp_6(k,l,j,n,i,m) + tmp_6(k,l,j,n,m,i) + &
                             tmp_6(k,l,m,i,j,n) + tmp_6(k,l,m,i,n,j) + tmp_6(k,l,m,j,i,n) + &
                             tmp_6(k,l,m,j,n,i) + tmp_6(k,l,m,n,i,j) + tmp_6(k,l,m,n,j,i) + &
                             tmp_6(k,l,n,i,j,m) + tmp_6(k,l,n,i,m,j) + tmp_6(k,l,n,j,i,m) + &
                             tmp_6(k,l,n,j,m,i) + tmp_6(k,l,n,m,i,j) + tmp_6(k,l,n,m,j,i) + &
                             tmp_6(k,m,i,j,l,n) + tmp_6(k,m,i,j,n,l) + tmp_6(k,m,i,l,j,n) + &
                             tmp_6(k,m,i,l,n,j) + tmp_6(k,m,i,n,j,l) + tmp_6(k,m,i,n,l,j) + &
                             tmp_6(k,m,j,i,l,n) + tmp_6(k,m,j,i,n,l) + tmp_6(k,m,j,l,i,n) + &
                             tmp_6(k,m,j,l,n,i) + tmp_6(k,m,j,n,i,l) + tmp_6(k,m,j,n,l,i) + &
                             tmp_6(k,m,l,i,j,n) + tmp_6(k,m,l,i,n,j) + tmp_6(k,m,l,j,i,n) + &
                             tmp_6(k,m,l,j,n,i) + tmp_6(k,m,l,n,i,j) + tmp_6(k,m,l,n,j,i) + &
                             tmp_6(k,m,n,i,j,l) + tmp_6(k,m,n,i,l,j) + tmp_6(k,m,n,j,i,l) + &
                             tmp_6(k,m,n,j,l,i) + tmp_6(k,m,n,l,i,j) + tmp_6(k,m,n,l,j,i) + &
                             tmp_6(k,n,i,j,l,m) + tmp_6(k,n,i,j,m,l) + tmp_6(k,n,i,l,j,m) + &
                             tmp_6(k,n,i,l,m,j) + tmp_6(k,n,i,m,j,l) + tmp_6(k,n,i,m,l,j) + &
                             tmp_6(k,n,j,i,l,m) + tmp_6(k,n,j,i,m,l) + tmp_6(k,n,j,l,i,m) + &
                             tmp_6(k,n,j,l,m,i) + tmp_6(k,n,j,m,i,l) + tmp_6(k,n,j,m,l,i) + &
                             tmp_6(k,n,l,i,j,m) + tmp_6(k,n,l,i,m,j) + tmp_6(k,n,l,j,i,m) + &
                             tmp_6(k,n,l,j,m,i) + tmp_6(k,n,l,m,i,j) + tmp_6(k,n,l,m,j,i) + &
                             tmp_6(k,n,m,i,j,l) + tmp_6(k,n,m,i,l,j) + tmp_6(k,n,m,j,i,l) + &
                             tmp_6(k,n,m,j,l,i) + tmp_6(k,n,m,l,i,j) + tmp_6(k,n,m,l,j,i) + &
                             tmp_6(l,i,j,k,m,n) + tmp_6(l,i,j,k,n,m) + tmp_6(l,i,j,m,k,n) + &
                             tmp_6(l,i,j,m,n,k) + tmp_6(l,i,j,n,k,m) + tmp_6(l,i,j,n,m,k) + &
                             tmp_6(l,i,k,j,m,n) + tmp_6(l,i,k,j,n,m) + tmp_6(l,i,k,m,j,n) + &
                             tmp_6(l,i,k,m,n,j) + tmp_6(l,i,k,n,j,m) + tmp_6(l,i,k,n,m,j) + &
                             tmp_6(l,i,m,j,k,n) + tmp_6(l,i,m,j,n,k) + tmp_6(l,i,m,k,j,n) + &
                             tmp_6(l,i,m,k,n,j) + tmp_6(l,i,m,n,j,k) + tmp_6(l,i,m,n,k,j) + &
                             tmp_6(l,i,n,j,k,m) + tmp_6(l,i,n,j,m,k) + tmp_6(l,i,n,k,j,m) + &
                             tmp_6(l,i,n,k,m,j) + tmp_6(l,i,n,m,j,k) + tmp_6(l,i,n,m,k,j) + &
                             tmp_6(l,j,i,k,m,n) + tmp_6(l,j,i,k,n,m) + tmp_6(l,j,i,m,k,n) + &
                             tmp_6(l,j,i,m,n,k) + tmp_6(l,j,i,n,k,m) + tmp_6(l,j,i,n,m,k) + &
                             tmp_6(l,j,k,i,m,n) + tmp_6(l,j,k,i,n,m) + tmp_6(l,j,k,m,i,n) + &
                             tmp_6(l,j,k,m,n,i) + tmp_6(l,j,k,n,i,m) + tmp_6(l,j,k,n,m,i) + &
                             tmp_6(l,j,m,i,k,n) + tmp_6(l,j,m,i,n,k) + tmp_6(l,j,m,k,i,n) + &
                             tmp_6(l,j,m,k,n,i) + tmp_6(l,j,m,n,i,k) + tmp_6(l,j,m,n,k,i) + &
                             tmp_6(l,j,n,i,k,m) + tmp_6(l,j,n,i,m,k) + tmp_6(l,j,n,k,i,m) + &
                             tmp_6(l,j,n,k,m,i) + tmp_6(l,j,n,m,i,k) + tmp_6(l,j,n,m,k,i) + &
                             tmp_6(l,k,i,j,m,n) + tmp_6(l,k,i,j,n,m) + tmp_6(l,k,i,m,j,n) + &
                             tmp_6(l,k,i,m,n,j) + tmp_6(l,k,i,n,j,m) + tmp_6(l,k,i,n,m,j) + &
                             tmp_6(l,k,j,i,m,n) + tmp_6(l,k,j,i,n,m) + tmp_6(l,k,j,m,i,n) + &
                             tmp_6(l,k,j,m,n,i) + tmp_6(l,k,j,n,i,m) + tmp_6(l,k,j,n,m,i) + &
                             tmp_6(l,k,m,i,j,n) + tmp_6(l,k,m,i,n,j) + tmp_6(l,k,m,j,i,n) + &
                             tmp_6(l,k,m,j,n,i) + tmp_6(l,k,m,n,i,j) + tmp_6(l,k,m,n,j,i) + &
                             tmp_6(l,k,n,i,j,m) + tmp_6(l,k,n,i,m,j) + tmp_6(l,k,n,j,i,m) + &
                             tmp_6(l,k,n,j,m,i) + tmp_6(l,k,n,m,i,j) + tmp_6(l,k,n,m,j,i) + &
                             tmp_6(l,m,i,j,k,n) + tmp_6(l,m,i,j,n,k) + tmp_6(l,m,i,k,j,n) + &
                             tmp_6(l,m,i,k,n,j) + tmp_6(l,m,i,n,j,k) + tmp_6(l,m,i,n,k,j) + &
                             tmp_6(l,m,j,i,k,n) + tmp_6(l,m,j,i,n,k) + tmp_6(l,m,j,k,i,n) + &
                             tmp_6(l,m,j,k,n,i) + tmp_6(l,m,j,n,i,k) + tmp_6(l,m,j,n,k,i) + &
                             tmp_6(l,m,k,i,j,n) + tmp_6(l,m,k,i,n,j) + tmp_6(l,m,k,j,i,n) + &
                             tmp_6(l,m,k,j,n,i) + tmp_6(l,m,k,n,i,j) + tmp_6(l,m,k,n,j,i) + &
                             tmp_6(l,m,n,i,j,k) + tmp_6(l,m,n,i,k,j) + tmp_6(l,m,n,j,i,k) + &
                             tmp_6(l,m,n,j,k,i) + tmp_6(l,m,n,k,i,j) + tmp_6(l,m,n,k,j,i) + &
                             tmp_6(l,n,i,j,k,m) + tmp_6(l,n,i,j,m,k) + tmp_6(l,n,i,k,j,m) + &
                             tmp_6(l,n,i,k,m,j) + tmp_6(l,n,i,m,j,k) + tmp_6(l,n,i,m,k,j) + &
                             tmp_6(l,n,j,i,k,m) + tmp_6(l,n,j,i,m,k) + tmp_6(l,n,j,k,i,m) + &
                             tmp_6(l,n,j,k,m,i) + tmp_6(l,n,j,m,i,k) + tmp_6(l,n,j,m,k,i) + &
                             tmp_6(l,n,k,i,j,m) + tmp_6(l,n,k,i,m,j) + tmp_6(l,n,k,j,i,m) + &
                             tmp_6(l,n,k,j,m,i) + tmp_6(l,n,k,m,i,j) + tmp_6(l,n,k,m,j,i) + &
                             tmp_6(l,n,m,i,j,k) + tmp_6(l,n,m,i,k,j) + tmp_6(l,n,m,j,i,k) + &
                             tmp_6(l,n,m,j,k,i) + tmp_6(l,n,m,k,i,j) + tmp_6(l,n,m,k,j,i) + &
                             tmp_6(m,i,j,k,l,n) + tmp_6(m,i,j,k,n,l) + tmp_6(m,i,j,l,k,n) + &
                             tmp_6(m,i,j,l,n,k) + tmp_6(m,i,j,n,k,l) + tmp_6(m,i,j,n,l,k) + &
                             tmp_6(m,i,k,j,l,n) + tmp_6(m,i,k,j,n,l) + tmp_6(m,i,k,l,j,n) + &
                             tmp_6(m,i,k,l,n,j) + tmp_6(m,i,k,n,j,l) + tmp_6(m,i,k,n,l,j) + &
                             tmp_6(m,i,l,j,k,n) + tmp_6(m,i,l,j,n,k) + tmp_6(m,i,l,k,j,n) + &
                             tmp_6(m,i,l,k,n,j) + tmp_6(m,i,l,n,j,k) + tmp_6(m,i,l,n,k,j) + &
                             tmp_6(m,i,n,j,k,l) + tmp_6(m,i,n,j,l,k) + tmp_6(m,i,n,k,j,l) + &
                             tmp_6(m,i,n,k,l,j) + tmp_6(m,i,n,l,j,k) + tmp_6(m,i,n,l,k,j) + &
                             tmp_6(m,j,i,k,l,n) + tmp_6(m,j,i,k,n,l) + tmp_6(m,j,i,l,k,n) + &
                             tmp_6(m,j,i,l,n,k) + tmp_6(m,j,i,n,k,l) + tmp_6(m,j,i,n,l,k) + &
                             tmp_6(m,j,k,i,l,n) + tmp_6(m,j,k,i,n,l) + tmp_6(m,j,k,l,i,n) + &
                             tmp_6(m,j,k,l,n,i) + tmp_6(m,j,k,n,i,l) + tmp_6(m,j,k,n,l,i) + &
                             tmp_6(m,j,l,i,k,n) + tmp_6(m,j,l,i,n,k) + tmp_6(m,j,l,k,i,n) + &
                             tmp_6(m,j,l,k,n,i) + tmp_6(m,j,l,n,i,k) + tmp_6(m,j,l,n,k,i) + &
                             tmp_6(m,j,n,i,k,l) + tmp_6(m,j,n,i,l,k) + tmp_6(m,j,n,k,i,l) + &
                             tmp_6(m,j,n,k,l,i) + tmp_6(m,j,n,l,i,k) + tmp_6(m,j,n,l,k,i) + &
                             tmp_6(m,k,i,j,l,n) + tmp_6(m,k,i,j,n,l) + tmp_6(m,k,i,l,j,n) + &
                             tmp_6(m,k,i,l,n,j) + tmp_6(m,k,i,n,j,l) + tmp_6(m,k,i,n,l,j) + &
                             tmp_6(m,k,j,i,l,n) + tmp_6(m,k,j,i,n,l) + tmp_6(m,k,j,l,i,n) + &
                             tmp_6(m,k,j,l,n,i) + tmp_6(m,k,j,n,i,l) + tmp_6(m,k,j,n,l,i) + &
                             tmp_6(m,k,l,i,j,n) + tmp_6(m,k,l,i,n,j) + tmp_6(m,k,l,j,i,n) + &
                             tmp_6(m,k,l,j,n,i) + tmp_6(m,k,l,n,i,j) + tmp_6(m,k,l,n,j,i) + &
                             tmp_6(m,k,n,i,j,l) + tmp_6(m,k,n,i,l,j) + tmp_6(m,k,n,j,i,l) + &
                             tmp_6(m,k,n,j,l,i) + tmp_6(m,k,n,l,i,j) + tmp_6(m,k,n,l,j,i) + &
                             tmp_6(m,l,i,j,k,n) + tmp_6(m,l,i,j,n,k) + tmp_6(m,l,i,k,j,n) + &
                             tmp_6(m,l,i,k,n,j) + tmp_6(m,l,i,n,j,k) + tmp_6(m,l,i,n,k,j) + &
                             tmp_6(m,l,j,i,k,n) + tmp_6(m,l,j,i,n,k) + tmp_6(m,l,j,k,i,n) + &
                             tmp_6(m,l,j,k,n,i) + tmp_6(m,l,j,n,i,k) + tmp_6(m,l,j,n,k,i) + &
                             tmp_6(m,l,k,i,j,n) + tmp_6(m,l,k,i,n,j) + tmp_6(m,l,k,j,i,n) + &
                             tmp_6(m,l,k,j,n,i) + tmp_6(m,l,k,n,i,j) + tmp_6(m,l,k,n,j,i) + &
                             tmp_6(m,l,n,i,j,k) + tmp_6(m,l,n,i,k,j) + tmp_6(m,l,n,j,i,k) + &
                             tmp_6(m,l,n,j,k,i) + tmp_6(m,l,n,k,i,j) + tmp_6(m,l,n,k,j,i) + &
                             tmp_6(m,n,i,j,k,l) + tmp_6(m,n,i,j,l,k) + tmp_6(m,n,i,k,j,l) + &
                             tmp_6(m,n,i,k,l,j) + tmp_6(m,n,i,l,j,k) + tmp_6(m,n,i,l,k,j) + &
                             tmp_6(m,n,j,i,k,l) + tmp_6(m,n,j,i,l,k) + tmp_6(m,n,j,k,i,l) + &
                             tmp_6(m,n,j,k,l,i) + tmp_6(m,n,j,l,i,k) + tmp_6(m,n,j,l,k,i) + &
                             tmp_6(m,n,k,i,j,l) + tmp_6(m,n,k,i,l,j) + tmp_6(m,n,k,j,i,l) + &
                             tmp_6(m,n,k,j,l,i) + tmp_6(m,n,k,l,i,j) + tmp_6(m,n,k,l,j,i) + &
                             tmp_6(m,n,l,i,j,k) + tmp_6(m,n,l,i,k,j) + tmp_6(m,n,l,j,i,k) + &
                             tmp_6(m,n,l,j,k,i) + tmp_6(m,n,l,k,i,j) + tmp_6(m,n,l,k,j,i) + &
                             tmp_6(n,i,j,k,l,m) + tmp_6(n,i,j,k,m,l) + tmp_6(n,i,j,l,k,m) + &
                             tmp_6(n,i,j,l,m,k) + tmp_6(n,i,j,m,k,l) + tmp_6(n,i,j,m,l,k) + &
                             tmp_6(n,i,k,j,l,m) + tmp_6(n,i,k,j,m,l) + tmp_6(n,i,k,l,j,m) + &
                             tmp_6(n,i,k,l,m,j) + tmp_6(n,i,k,m,j,l) + tmp_6(n,i,k,m,l,j) + &
                             tmp_6(n,i,l,j,k,m) + tmp_6(n,i,l,j,m,k) + tmp_6(n,i,l,k,j,m) + &
                             tmp_6(n,i,l,k,m,j) + tmp_6(n,i,l,m,j,k) + tmp_6(n,i,l,m,k,j) + &
                             tmp_6(n,i,m,j,k,l) + tmp_6(n,i,m,j,l,k) + tmp_6(n,i,m,k,j,l) + &
                             tmp_6(n,i,m,k,l,j) + tmp_6(n,i,m,l,j,k) + tmp_6(n,i,m,l,k,j) + &
                             tmp_6(n,j,i,k,l,m) + tmp_6(n,j,i,k,m,l) + tmp_6(n,j,i,l,k,m) + &
                             tmp_6(n,j,i,l,m,k) + tmp_6(n,j,i,m,k,l) + tmp_6(n,j,i,m,l,k) + &
                             tmp_6(n,j,k,i,l,m) + tmp_6(n,j,k,i,m,l) + tmp_6(n,j,k,l,i,m) + &
                             tmp_6(n,j,k,l,m,i) + tmp_6(n,j,k,m,i,l) + tmp_6(n,j,k,m,l,i) + &
                             tmp_6(n,j,l,i,k,m) + tmp_6(n,j,l,i,m,k) + tmp_6(n,j,l,k,i,m) + &
                             tmp_6(n,j,l,k,m,i) + tmp_6(n,j,l,m,i,k) + tmp_6(n,j,l,m,k,i) + &
                             tmp_6(n,j,m,i,k,l) + tmp_6(n,j,m,i,l,k) + tmp_6(n,j,m,k,i,l) + &
                             tmp_6(n,j,m,k,l,i) + tmp_6(n,j,m,l,i,k) + tmp_6(n,j,m,l,k,i) + &
                             tmp_6(n,k,i,j,l,m) + tmp_6(n,k,i,j,m,l) + tmp_6(n,k,i,l,j,m) + &
                             tmp_6(n,k,i,l,m,j) + tmp_6(n,k,i,m,j,l) + tmp_6(n,k,i,m,l,j) + &
                             tmp_6(n,k,j,i,l,m) + tmp_6(n,k,j,i,m,l) + tmp_6(n,k,j,l,i,m) + &
                             tmp_6(n,k,j,l,m,i) + tmp_6(n,k,j,m,i,l) + tmp_6(n,k,j,m,l,i) + &
                             tmp_6(n,k,l,i,j,m) + tmp_6(n,k,l,i,m,j) + tmp_6(n,k,l,j,i,m) + &
                             tmp_6(n,k,l,j,m,i) + tmp_6(n,k,l,m,i,j) + tmp_6(n,k,l,m,j,i) + &
                             tmp_6(n,k,m,i,j,l) + tmp_6(n,k,m,i,l,j) + tmp_6(n,k,m,j,i,l) + &
                             tmp_6(n,k,m,j,l,i) + tmp_6(n,k,m,l,i,j) + tmp_6(n,k,m,l,j,i) + &
                             tmp_6(n,l,i,j,k,m) + tmp_6(n,l,i,j,m,k) + tmp_6(n,l,i,k,j,m) + &
                             tmp_6(n,l,i,k,m,j) + tmp_6(n,l,i,m,j,k) + tmp_6(n,l,i,m,k,j) + &
                             tmp_6(n,l,j,i,k,m) + tmp_6(n,l,j,i,m,k) + tmp_6(n,l,j,k,i,m) + &
                             tmp_6(n,l,j,k,m,i) + tmp_6(n,l,j,m,i,k) + tmp_6(n,l,j,m,k,i) + &
                             tmp_6(n,l,k,i,j,m) + tmp_6(n,l,k,i,m,j) + tmp_6(n,l,k,j,i,m) + &
                             tmp_6(n,l,k,j,m,i) + tmp_6(n,l,k,m,i,j) + tmp_6(n,l,k,m,j,i) + &
                             tmp_6(n,l,m,i,j,k) + tmp_6(n,l,m,i,k,j) + tmp_6(n,l,m,j,i,k) + &
                             tmp_6(n,l,m,j,k,i) + tmp_6(n,l,m,k,i,j) + tmp_6(n,l,m,k,j,i) + &
                             tmp_6(n,m,i,j,k,l) + tmp_6(n,m,i,j,l,k) + tmp_6(n,m,i,k,j,l) + &
                             tmp_6(n,m,i,k,l,j) + tmp_6(n,m,i,l,j,k) + tmp_6(n,m,i,l,k,j) + &
                             tmp_6(n,m,j,i,k,l) + tmp_6(n,m,j,i,l,k) + tmp_6(n,m,j,k,i,l) + &
                             tmp_6(n,m,j,k,l,i) + tmp_6(n,m,j,l,i,k) + tmp_6(n,m,j,l,k,i) + &
                             tmp_6(n,m,k,i,j,l) + tmp_6(n,m,k,i,l,j) + tmp_6(n,m,k,j,i,l) + &
                             tmp_6(n,m,k,j,l,i) + tmp_6(n,m,k,l,i,j) + tmp_6(n,m,k,l,j,i) + &
                             tmp_6(n,m,l,i,j,k) + tmp_6(n,m,l,i,k,j) + tmp_6(n,m,l,j,i,k) + &
                             tmp_6(n,m,l,j,k,i) + tmp_6(n,m,l,k,i,j) + tmp_6(n,m,l,k,j,i)
                             
                         tmp_6(i,j,k,l,m,n) = r

                      end do
                   end do
                end do
             end do
          end do
       end do

       h = 0
       do i = 1, ncor
          do j = i, ncor
             do k = j, ncor
                do l = k, ncor
                   do m = l, ncor
                      do n = m, ncor
                         h = h + 1
                         ! MR: UNSURE ABOUT ORIGINS OF FACTOR 2 IN NEXT LINE
                         ave(h) = 2 * tmp_6(i,j,k,l,m,n)
                      end do
                   end do
                end do
             end do
          end do
       end do

       deallocate(tmp_6)


    else
       print *, 'rsp_twoave error: Contribution not implemented or in wrong order - ', &
                (' ' // f(i), i=1,nf)
       call quit('rsp_twoave error: Contribution not implemented or in wrong order')
    end if

  end subroutine

  !> Contract 2-electron integrals perturbed by fields 'f' with density
  !> matrix 'dens', and add to Fock matrices 'fock' Average 2-electron integrals perturbed by fields f over the
  !> product of density matrces D1 and D2
  subroutine rsp_twoint(nr_ao, nf, f, c, nc, dens, propsize, fock)
    use eri_contractions, only: ctr_arg
    use eri_basis_loops,  only: unopt_geodiff_loop
    use interface_interest
    use rsp_indices_and_addressing, only: make_triangulated_indices
    !> number of fields
    integer,              intent(in)    :: nf, propsize
    !> field labels in std order
    character(4),         intent(in)    :: f(nf)
    !> first and number of- components in each field
    integer,              intent(in)    :: c(nf), nc(nf)
    !> density matrix to average over
    type(matrix), target, intent(in)    :: dens
    !> Fock matrix to which two-electron contribution is ADDED
    type(matrix), target, intent(inout) :: fock(propsize)
    type(matrix), allocatable, target, dimension(:) :: tmpfock
    !--------------------------------------------------
    real(8), pointer :: null_ptr(:) !because null() isn't f90
    real(8)  :: fdistep = 2d0**(-25)
    integer       h, hincr, i, j, k, n, ij, ijk, incr, ncor, nr_ao
    integer       lwork
    type(ctr_arg) arg(1)
    type(matrix)  A !scratch
    real(8)          :: dummy(1)
    integer          :: idummy(1)
    real(8), allocatable :: f77_memory(:)
    real(8), allocatable :: temp_fmat(:)
    real(8), allocatable :: temp_dmat(:)
    !--------------------------------------------------
    integer, allocatable, dimension(:,:) :: indices

    if (any(f == 'EL  ')) then
       ! nothing to add
       return
    end if

    if (nf == 1 .and. f(1) == 'MAG ') then
       !fixme wild guess
       lwork = 100000000
       allocate(f77_memory(lwork))

       allocate(temp_dmat(nr_ao*nr_ao))
       call dcopy(nr_ao*nr_ao, dens%elms, 1, temp_dmat, 1)

       allocate(temp_fmat(3*nr_ao*nr_ao))
       temp_fmat = 0.0d0
       call grcont(f77_memory, lwork, temp_fmat, (3*nr_ao*nr_ao), .false., .true., 1, &
                   0, .false., .true., temp_dmat, 1)
       do i = 1, 3
          call dcopy(nr_ao*nr_ao, temp_fmat((i-1)*nr_ao*nr_ao + 1), 1, fock(i)%elms, 1)
       end do
       deallocate(f77_memory)
       deallocate(temp_fmat)
       deallocate(temp_dmat)
    end if

#ifdef PRG_DALTON
! Begin new general geo code


       nullify(null_ptr)

       if (nf==0) then
       call mat_init(A, fock(1)%nrow, fock(1)%ncol)
          A = 0*dens
          call mat_ensure_alloc(A)
! write(*,*) 'A', A%elms
! write(*,*) 'dens', dens%elms

          call interface_scf_get_g(dens, A)

! write(*,*) 'A 2', A%elms
          fock(1) = fock(1) + A
          A = 0

       else

          if ((nf == count(f == 'GEO '))) then

             ncor = 3 * get_nr_atoms()
             allocate(indices(propsize, nf))
             call make_triangulated_indices(1, (/1, nf, nc(1)/), propsize, indices)

             do i = 1, propsize

                if (iszero(fock(i))) then
                   call mat_ensure_alloc(fock(i))
                   fock(i)%elms = 0 !ajt FIXME use mat_axpy
                end if
             
                k = 1

                do j = 1, nf

                   k = k + (indices(i,j) - 1) * ( nc(1)**(nf - j ) )

                end do

                   arg(1) = ctr_arg(nf, k, ncor, dens, fock(i), null_ptr)
                   call unopt_geodiff_loop(basis_large, basis_small, arg)

             end do

             deallocate(indices)

          end if

       end if


! End new general geo code

! MaR: I don't understand all the program-specific ifdefs below, so I kept this
! non-general code for "not-Dalton" runs

#else
       if (nf==0) then
          A = 0*dens
          call mat_ensure_alloc(A)
          call interface_scf_get_g(dens, A)
          fock(1) = fock(1) + A
          A = 0
       else if (nf==1 .and. f(1)=='GEO ') then
          n = nr_ao
          do i = 0, nc(1)-1
             if (iszero(fock(i+1))) then
                call mat_ensure_alloc(fock(i+1))
             end if

             call interest_mpi_wake_up()
             call interest_get_int(dens%nrow, dens%elms, fock(i+1)%elms, 1, i+1, 0, dummy)
   
          end do

       else
          print *, 'error in rsp_twoint: not implemented or in wrong order - ', &
                  (' ' // f(i), i=1,nf)
          call quit('error in rsp_twoint: not implemented or in wrong order')
       end if
#endif

  end subroutine

  function rank_one_pointer(siz, arr) result(ptr)
     integer,         intent(in) :: siz
     real(8), target, intent(in) :: arr(siz)
     real(8), pointer            :: ptr(:)
     ptr => arr
  end function

end module
