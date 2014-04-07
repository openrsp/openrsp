module interface_pelib

    use polarizable_embedding, only: pe_master

    implicit none

    private

    public :: pe_response_operator, pe_add_full_operator

contains

subroutine pe_add_full_operator(dmat, fmat, energy)

    real(8), dimension(:,:), intent(in) :: dmat
    real(8), dimension(:,:), intent(out) :: fmat
    real(8), dimension(:), intent(out) :: energy

    integer :: nbas
    integer :: i, j, k
    real(8), dimension(:), allocatable :: packed_fmat, packed_dmat

    fmat = 0.0d0
    energy = 0.0d0
    nbas = int(sqrt(real(size(dmat))))
    allocate(packed_fmat(nbas*(nbas+1)/2))
    allocate(packed_dmat(nbas*(nbas+1)/2))
    packed_fmat = 0.0d0
    packed_dmat = 0.0d0
    k = 1
    do i = 1, nbas
        do j = 1, i
            if (i == j) then
                packed_dmat(k) = dmat(i,i)
            else
                packed_dmat(k) = dmat(j,i) + dmat(i,j)
            end if
            k = k + 1
        end do
    end do
    call pe_master(runtype='fock', denmats=packed_dmat, fckmats=packed_fmat,&
                  & nmats=1, energies=energy)
    deallocate(packed_dmat)
    k = 1
    do i = 1, nbas
        do j = 1, i
            if (i == j) then
                fmat(i,i) = packed_fmat(k)
            else
                fmat(j,i) = packed_fmat(k)
                fmat(i,j) = packed_fmat(k)
            end if
            k = k + 1
        end do
    end do
    deallocate(packed_fmat)

end subroutine pe_add_full_operator

subroutine pe_response_operator(dmat, fmat, nbas, ndens)

    integer, intent(in) :: nbas, ndens
    real(8), dimension(:,:), intent(in) :: dmat
    real(8), dimension(:,:), intent(out) :: fmat

    integer :: i, j, k, l, m
    integer :: nnbas
    real(8), dimension(:), allocatable :: packed_fmat, packed_dmat

    if (ndens > 1) then
        stop 'ERROR: pe_response_operator can only use one density matrix'
    end if

    fmat = 0.0d0
    nnbas = nbas * (nbas + 1) / 2
    allocate(packed_fmat(ndens*nnbas))
    allocate(packed_dmat(ndens*nnbas))
    packed_fmat = 0.0d0
    packed_dmat = 0.d00
    k = 1
    do i = 1, nbas
        do j = 1, i
            if (i == j) then
                packed_dmat(k) = dmat(i,i)
            else
                packed_dmat(k) = dmat(j,i) + dmat(i,j)
            end if
            k = k + 1
        end do
    end do
!    do i = 1, ndens
!        j = 1 + (i - 1) * nbas * nbas
!        k = i * nbas*nbas
!        l = 1 + (i - 1) * nnbas
!        m = i * nnbas
!        call dgefsp(nbas, dmats(j), packed_dmats(l))
!    end do
    call pe_master(runtype='response', denmats=packed_dmat,&
                  & fckmats=packed_fmat, nmats=ndens)
    deallocate(packed_dmat)
    k = 1
    do i = 1, nbas
        do j = 1, i
            if (i == j) then
                fmat(i,i) = packed_fmat(k)
            else
                fmat(j,i) = packed_fmat(k)
                fmat(i,j) = packed_fmat(k)
            end if
            k = k + 1
        end do
    end do
    deallocate(packed_fmat)
!    do i = 1, ndens
!        j = 1 + (i - 1) * nnbas
!        k = i * nnbas
!        l = 1 + (i - 1) * nbas * nbas
!        m = i * nbas*nbas
!        call dsptge(nbas, packed_fmats(j), fmats(l))
!    end do

end subroutine pe_response_operator

end module interface_pelib
