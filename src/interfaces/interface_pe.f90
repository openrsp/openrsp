module interface_pe

    use matrix_defop
    use matrix_lowlevel

    implicit none

    private

    public :: pe_rsp

contains

subroutine pe_rsp(dens, fock, nr_ao, propsize)

!    use pe_variables, only: peqm
    use polarizable_embedding, only: pe_master

    !> number of fields
    integer, intent(in) :: propsize, nr_ao
    !> density matrix
    type(matrix), target, intent(in) :: dens
    !> Fock matrix to which the PE contribution is ADDED
    type(matrix), target, intent(inout) :: fock(propsize)

    integer :: i, lwork
    real(8), dimension(:), allocatable :: f77_memory
    real(8), dimension(:,:), allocatable :: fmat_full, dmat_full
    real(8), dimension(:), allocatable :: fmat_packed, dmat_packed

    !fixme wild guess
    lwork = 100000000
    allocate(f77_memory(lwork))
    allocate(dmat_packed(nr_ao * (nr_ao + 1) / 2), dmat_full(nr_ao,nr_ao))

    dmat_full(:,:) = dens%elms(:,:,1)
    call dgefsp(nr_ao, dmat_full, dmat_packed)
    deallocate(dmat_full)

    allocate(fmat_packed(nr_ao * (nr_ao + 1) / 2))

    fmat_packed = 0.0d0

    call pe_master(runtype='response', denmats=dmat_packed,&
                  & fckmats=fmat_packed, nmats=1, dalwrk=f77_memory)

    deallocate(dmat_packed, f77_memory)
    allocate(fmat_full(nr_ao, nr_ao))
    fmat_full = 0.0d0

    call dsptge(nr_ao, fmat_packed, fmat_full)

    do i = 1, propsize
       fock(i)%elms(:,:,1) = fock(i)%elms(:,:,1) + 2.0d0*fmat_full(:,:)
    enddo

    deallocate(fmat_full)
    deallocate(fmat_packed)

end subroutine pe_rsp

end module interface_pe
