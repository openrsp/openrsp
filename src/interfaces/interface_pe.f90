module interface_pe

    use matrix_defop
    use matrix_lowlevel

    implicit none

    private

    public :: pe_rsp

contains

subroutine pe_rsp(nr_ao, nf, f, c, nc, dens, propsize, fock)

    use pe_variables, only: peqm
    use polarizable_embedding, only: pe_master

    !> number of fields
    integer, intent(in) :: nf, propsize, nr_ao
    !> field labels in std order
    character(4), intent(in) :: f(nf)
    !> first and number of- components in each field
    integer, intent(in) :: c(nf), nc(nf)
    !> density matrix
    type(matrix), target, intent(in) :: dens
    !> Fock matrix to which the PE contribution is ADDED
    type(matrix), target, intent(inout) :: fock(propsize)

    integer :: i, lwork
    real(8), dimension(:), allocatable :: f77_memory
    real(8), dimension(:,:), allocatable :: fmat_full, dmat_full
    real(8), dimension(:), allocatable :: fmat_packed, dmat_packed

    logical :: debug
    debug = .false.

    if (.not. peqm) return

    if (any(f /= 'EL  ')) then
        stop 'ERROR: PE-OpenRSP not implemented for other than EL.'
    end if

    !fixme wild guess
    lwork = 100000000
    allocate(f77_memory(lwork))
    allocate(dmat_packed(nr_ao * (nr_ao + 1) / 2), dmat_full(nr_ao,nr_ao))

    dmat_full(:,:) = dens%elms(:,:,1)
!    call dgefsp(nr_ao, dmat_full, dmat_packed)
    call dsitsp(nr_ao, dmat_full, dmat_packed)
    if (debug) then
       call dsptge(nr_ao, dmat_packed, dmat_full)

       do i = 1, nr_ao
          print*, 'i is', i
          print *,dens%elms(i,:,1)
          print*, ' '
          print *,dmat_full(i,:)
       enddo
    endif
    deallocate(dmat_full)

    allocate(fmat_packed(nr_ao * (nr_ao + 1) / 2))
    fmat_packed = 0.0d0

    call pe_master(runtype='response', denmats=dmat_packed,&
                  & fckmats=fmat_packed, nmats=1, dalwrk=f77_memory)

    deallocate(dmat_packed, f77_memory)

    allocate(fmat_full(nr_ao, nr_ao))
    fmat_full = 0.0d0
    call dsptge(nr_ao, fmat_packed, fmat_full)
    fock(propsize)%elms(:,:,1) = fock(propsize)%elms(:,:,1) + fmat_full(:,:)
    deallocate(fmat_full)
    deallocate(fmat_packed)

end subroutine pe_rsp

end module interface_pe
