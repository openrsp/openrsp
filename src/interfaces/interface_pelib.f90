module interface_pelib

    use matrix_defop
    use matrix_lowlevel
    use polarizable_embedding, only: pe_master

    implicit none

    private

    public :: pe_response_operator, pe_add_full_operator

contains

subroutine pe_add_full_operator(dmat, fmat, energy)

    real(8), dimension(:), intent(in) :: dmat
    real(8), dimension(:), intent(out) :: fmat
    real(8), dimension(:), intent(out) :: energy

    integer :: nbas
    real(8), dimension(1) :: temp
    real(8), dimension(:), allocatable :: packed_fmat, packed_dmat

    fmat = 0.0d0
    energy = 0.0d0
    nbas = int(sqrt(real(size(dmat))))
    allocate(packed_fmat(nbas * (nbas + 1) / 2))
    allocate(packed_dmat(nbas * (nbas + 1) / 2))
    packed_fmat = 0.0d0
    packed_dmat = 0.0d0
    call dgefsp(nbas, dmat, packed_dmat)
    call pe_master('fock', packed_dmat, packed_fmat, 1, energy, temp)
    call dsptge(nbas, packed_fmat, fmat)
    deallocate(packed_fmat, packed_dmat)

end subroutine pe_add_full_operator

subroutine pe_response_operator(dmats, fmats, nbas, ndens)

    integer, intent(in) :: nbas, ndens
    real(8), dimension(:), intent(in) :: dmats
    real(8), dimension(:), intent(out) :: fmats

    integer :: i, j, k, l, m
    integer :: nnbas
    real(8), dimension(1) :: temp
    real(8), dimension(:), allocatable :: packed_fmats, packed_dmats

    fmats = 0.0d0
    return
    nnbas = nbas * (nbas + 1) / 2
    allocate(packed_fmats(ndens * nnbas))
    allocate(packed_dmats(ndens * nnbas))
    packed_fmats = 0.0d0
    packed_dmats = 0.d00
    do i = 1, ndens
        j = 1 + (i - 1) * nbas * nbas
        k = i * nbas*nbas
        l = 1 + (i - 1) * nnbas
        m = i * nnbas
        call dgefsp(nbas, dmats(j), packed_dmats(l))
    end do
    call pe_master(runtype='response', denmats=packed_dmats,&
                  & fckmats=packed_fmats, nmats=ndens, dalwrk=temp)
    do i = 1, ndens
        j = 1 + (i - 1) * nnbas
        k = i * nnbas
        l = 1 + (i - 1) * nbas * nbas
        m = i * nbas*nbas
        call dsptge(nbas, packed_fmats(j), fmats(l))
    end do

end subroutine pe_response_operator

end module interface_pelib
