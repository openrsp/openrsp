module nuc_contributions

   use interface_molecule
   use interface_io

   implicit none

   public nuclear_potential
   public gradient_nuc
   public hessian_nuc
   public cubicff_nuc
   public quarticff_nuc
   public nucrep_deriv
   public dipnuc_ifc
   public qdrnuc_ifc
   public dpgnuc_ifc
   public aatnuc_ifc

   private

contains

  !> (derivatives of) nuclear repulsion and nuclei--field interaction
  subroutine nuclear_potential(derv, ncor, field, ncomp, nucpot)
    !> order of geometry derivative requested
    integer,      intent(in)  :: derv
    !> number of geometry coordinates
    integer,      intent(in)  :: ncor
    !> one external potential label, or 'NONE'
    character(4), intent(in)  :: field
    !> number of components of the field (must be all, ie. 1/3/6)
    integer,      intent(in)  :: ncomp
    !> output tensor
    real(8),      intent(out) :: nucpot(ncor**derv * ncomp)
    !-------------------------------------------------------
    integer i, na
    na = get_nr_atoms()
    if (ncor /= 3*na) &
       call quit('rsp_backend nuclear_potential error: expected ncor = 3*natom')
    if ((field == 'NONE' .and. ncomp /= 1) .or. &
        (field == 'EL  ' .and. ncomp /= 3)) &
       call quit('rsp_backend nuclear_potential error: expected ncomp = 1 or 3')
    ! switch
    if (derv == 0 .and. field == 'NONE') then
       call quit('rsp_ error: unperturbed nuc.rep. not implemented')
    else if (derv == 0 .and. field == 'EL  ') then
       call DIPNUC_ifc(nucpot)
       nucpot = -nucpot !dip = -dE/dF
    else if (derv == 1 .and. field == 'NONE') then
       call gradient_nuc(na, nucpot)
    else if (derv == 2 .and. field == 'NONE') then
       call hessian_nuc(na, nucpot)
    else if (derv == 3 .and. field == 'NONE') then
       call cubicff_nuc(na, nucpot)
    else if (derv == 4 .and. field == 'NONE') then
       call quarticff_nuc(na, nucpot)
    else
       print *,'rsp_backend nuclear_potential error: not implented, derv =', &
               derv, '  field = ', field
       call quit('rsp_backend nuclear_potential error: not implented')
    end if
  end subroutine


   subroutine gradient_nuc(nr_atoms, tensor)

      integer, intent(in)  :: nr_atoms
      real(8), intent(out) :: tensor(3, nr_atoms)

      real(8)              :: r(3)
      real(8)              :: t(3)
      integer              :: i, j, k

      tensor = 0.0d0

      ! loop over pairs i<j of nuclei
      do j = 2, nr_atoms
         do i = 1, j-1
            r(1) = get_nuc_xyz(1, j) - get_nuc_xyz(1, i)
            r(2) = get_nuc_xyz(2, j) - get_nuc_xyz(2, i)
            r(3) = get_nuc_xyz(3, j) - get_nuc_xyz(3, i)

            ! apply scale factor
            t = r*get_nuc_charge(i)*get_nuc_charge(j)/(sum(r**2.0d0)**(3.0d0/2.0d0))

            tensor(:, i) = tensor(:, i) + t
            tensor(:, j) = tensor(:, j) - t
         end do
      end do

   end subroutine

   subroutine hessian_nuc(nr_atoms, tensor)

      integer, intent(in)  :: nr_atoms
      real(8), intent(out) :: tensor(3, nr_atoms, 3, nr_atoms)

      real(8)              :: r(3)
      real(8)              :: t(3, 3)
      real(8)              :: s
      integer              :: i, j, k

      tensor = 0.0d0

      ! loop over pairs i<j of nuclei
      do j = 2, nr_atoms
         do i = 1, j-1
            ! construct triple tensor prod with traces removed
            r(1) = get_nuc_xyz(1, j) - get_nuc_xyz(1, i)
            r(2) = get_nuc_xyz(2, j) - get_nuc_xyz(2, i)
            r(3) = get_nuc_xyz(3, j) - get_nuc_xyz(3, i)

            t = reshape((/(r(:)*r(k), k=1,3)/), (/3, 3/))
            s = (t(1, 1) + t(2, 2) + t(3, 3))/3.0d0
            do k = 1, 3
               t(k, k) = t(k, k) - s
            end do

            ! apply scale factor
            t = t*3.0d0*get_nuc_charge(i)*get_nuc_charge(j)/(sum(r**2.0d0)**(5.0d0/2.0d0))

            tensor(:, i, :, i) = tensor(:, i, :, i) + t !ii
            tensor(:, i, :, j) =                    - t !ij
            tensor(:, j, :, i) =                    - t !ji
            tensor(:, j, :, j) = tensor(:, j, :, j) + t !jj
         end do
      end do

   end subroutine


  !> Nuclear repulsion contribution to cubic force field, where \param na
  !> is the number of atoms, \param chg the charges of the nuclei, and
  !> \param cor the coordinates of the nuclei
  subroutine cubicff_nuc(na, cub)
     integer, intent(in)  :: na
     real(8), intent(out) :: cub(3,na,3,na,3,na)
     real(8) r(3), rrr(3,3,3), trace(3)
     integer i, j, k, l
     cub = 0 !start from zero
     ! loop over pairs i<j of nuclei
     do j = 2, na
        do i = 1, j-1
           ! construct triple tensor prod with traces removed
           r(1) = get_nuc_xyz(1, j) - get_nuc_xyz(1, i)
           r(2) = get_nuc_xyz(2, j) - get_nuc_xyz(2, i)
           r(3) = get_nuc_xyz(3, j) - get_nuc_xyz(3, i)
           rrr = reshape((/((r(:)*r(k)*r(l),k=1,3),l=1,3)/),(/3,3,3/))
           trace = rrr(:,1,1) + rrr(:,2,2) + rrr(:,3,3)
           do k = 1, 3
              rrr(:,k,k) = rrr(:,k,k) - trace/5
              rrr(k,:,k) = rrr(k,:,k) - trace/5
              rrr(k,k,:) = rrr(k,k,:) - trace/5
           end do
           ! apply scale factor 15 Qi Qj / r^7
           rrr = rrr * (15 * get_nuc_charge(i) * get_nuc_charge(j) / sum(r**2)**(7/2d0))
           cub(:,i,:,i,:,i) = cub(:,i,:,i,:,i) + rrr !iii
           cub(:,i,:,i,:,j) = -rrr !iij
           cub(:,i,:,j,:,i) = -rrr !iji
           cub(:,j,:,i,:,i) = -rrr !jii
           cub(:,i,:,j,:,j) =  rrr !ijj
           cub(:,j,:,i,:,j) =  rrr !jij
           cub(:,j,:,j,:,i) =  rrr !jji
           cub(:,j,:,j,:,j) = cub(:,j,:,j,:,j) - rrr !jjj
        end do
     end do
  end subroutine


  !> Nuclear repulsion contribution to quartic force field, where \param na
  !> is the number of atoms, \param chg the charges of the nuclei, and
  !> \param cor the coordinates of the nuclei
  !> Note: Dimensions 7 and 8 of 'qua' have been joined, to comply with
  !> fortran standard (max 7 dims)
  subroutine quarticff_nuc(na, qua)
     integer, intent(in)  :: na
     real(8), intent(out) :: qua(3,na,3,na,3,na,3*na)
     real(8) r(3), rrrr(3,3,3,3), trace(3,3), trace2
     integer i, j, k, l, m
     qua = 0 !start from zero
     ! loop over pairs i<j of nuclei
     do j = 2, na
        do i = 1, j-1
           ! construct quadruple tensor product with traces removed
           r(1) = get_nuc_xyz(1, j) - get_nuc_xyz(1, i)
           r(2) = get_nuc_xyz(2, j) - get_nuc_xyz(2, i)
           r(3) = get_nuc_xyz(3, j) - get_nuc_xyz(3, i)
           rrrr = reshape((/(((r(:)*r(k)*r(l)*r(m), &
                          k=1,3),l=1,3),m=1,3)/),(/3,3,3,3/))
           trace  = rrrr(:,:,1,1) + rrrr(:,:,2,2) + rrrr(:,:,3,3)
           trace2 = trace(1,1) + trace(2,2) + trace(3,3)
           do k = 1, 3
              trace(k,k) = trace(k,k) - trace2/10
           end do
           do k = 1, 3
              rrrr(:,:,k,k) = rrrr(:,:,k,k) - trace/7
              rrrr(:,k,:,k) = rrrr(:,k,:,k) - trace/7
              rrrr(k,:,:,k) = rrrr(k,:,:,k) - trace/7
              rrrr(:,k,k,:) = rrrr(:,k,k,:) - trace/7
              rrrr(k,:,k,:) = rrrr(k,:,k,:) - trace/7
              rrrr(k,k,:,:) = rrrr(k,k,:,:) - trace/7
           end do
           ! apply scale factor 105 Qi Qj / r^9
           rrrr = rrrr * (105 * get_nuc_charge(i) * get_nuc_charge(j) / sum(r**2)**(9/2d0))
           qua(:,i,:,i,:,i,3*i-2:3*i) = qua(:,i,:,i,:,i,3*i-2:3*i) + rrrr !iiii
           qua(:,i,:,i,:,i,3*j-2:3*j) = -rrrr !iiij
           qua(:,i,:,i,:,j,3*i-2:3*i) = -rrrr !iiji
           qua(:,i,:,j,:,i,3*i-2:3*i) = -rrrr !ijii
           qua(:,j,:,i,:,i,3*i-2:3*i) = -rrrr !jiii
           qua(:,i,:,i,:,j,3*j-2:3*j) =  rrrr !iijj
           qua(:,i,:,j,:,i,3*j-2:3*j) =  rrrr !ijij
           qua(:,j,:,i,:,i,3*j-2:3*j) =  rrrr !jiij
           qua(:,i,:,j,:,j,3*i-2:3*i) =  rrrr !ijji
           qua(:,j,:,i,:,j,3*i-2:3*i) =  rrrr !jiji
           qua(:,j,:,j,:,i,3*i-2:3*i) =  rrrr !jjii
           qua(:,i,:,j,:,j,3*j-2:3*j) = -rrrr !ijjj
           qua(:,j,:,i,:,j,3*j-2:3*j) = -rrrr !jijj
           qua(:,j,:,j,:,i,3*j-2:3*j) = -rrrr !jjij
           qua(:,j,:,j,:,j,3*i-2:3*i) = -rrrr !jjji
           qua(:,j,:,j,:,j,3*j-2:3*j) = qua(:,j,:,j,:,j,3*j-2:3*j) + rrrr !jjjj
        end do
     end do
  end subroutine



  ! derivatives of nuclear repulsion energy
  subroutine nucrep_deriv(na, charges, coords, deriv, nucrep)
    integer, intent(in)  :: na, deriv, charges(na)
    real(8), intent(in)  :: coords(3,na)
    real(8), intent(out) :: nucrep((3*na)**deriv)
    !--------------------------------------------
    real(8) displ(3), power(3**deriv)
    integer nuci, nucj, k, l, m, ncor_pow(deriv)
    ! powers of ncor=3*na, for navigating through nucrep
    ncor_pow(:) = (/((3*na)**k, k=0,deriv-1)/)
    ! start from zero
    nucrep(:) = 0
    ! loop over pairs i<j of nuclei
    do nucj = 2, na
       do nuci = 1, nucj-1
          ! internuclear displacement
          displ(:) = coords(:,nucj) - coords(:,nuci)
          ! first unperturbed repulsion energy
          power(1) = charges(nuci) * charges(nucj) / sqrt(sum(displ**2))
          ! construct deriv'th tensor power of displacement vector,
          ! carrying the prefactor QiQj (2k-1)!! / r^(2k+1)
          do k = 1, deriv
             l = 3**(k-1)
             power(2*l+1:3*l) = power(:l) * displ(3) * (2*k-1) / sum(displ**2)
             power(  l+1:2*l) = power(:l) * displ(2) * (2*k-1) / sum(displ**2)
             power(     :  l) = power(:l) * displ(1) * (2*k-1) / sum(displ**2)
          end do
          ! remove traces from displacement power
          if (deriv >= 2) call remove_traces(deriv, power)
          ! add traceless displacement power at the 2^deriv different
          ! places in the nucrep derivative tensor
          call place_power_tensor
       end do
    end do

  contains

    recursive subroutine remove_traces(m,a)
      integer, intent(in)    :: m
      real(8), intent(inout) :: a(3**(m-2),3,3)
      real(8) b(3**(m-2))
      integer i, j
      b(:) = a(:,1,1) + a(:,2,2) + a(:,3,3)
      if (m >= 4) call remove_traces(m-2,b)
      do j = 2, m
         do i = 1, j-1
            call subtract_from_diag( &
                         (deriv+m-1) * (1+(deriv-m)/2), &
                         3**(i-1), 3**(j-i-1), 3**(m-j), b, a)
         end do
      end do
    end subroutine

    subroutine subtract_from_diag(denom, l, m, t, b, a)
      integer, intent(in)    :: denom, l, m, t
      real(8), intent(in)    :: b(l,m,t)
      real(8), intent(inout) :: a(l,3,m,3,t)
      a(:,1,:,1,:) = a(:,1,:,1,:) - b(:,:,:)/denom
      a(:,2,:,2,:) = a(:,2,:,2,:) - b(:,:,:)/denom
      a(:,3,:,3,:) = a(:,3,:,3,:) - b(:,:,:)/denom
    end subroutine

    subroutine place_power_tensor
      integer offset, i, j, y, z, k, s
      logical neg
      offset = (3*(nuci-1)) * sum(ncor_pow)
      neg = .false.
      ! loop over all index permutations iiiii jiiii ... ijjjj jjjjj
      do k = 0, 2**deriv - 1
         ! loop over Cartesian indices xxxxx yxxxx ... yzzzz zzzzz
         i = offset + 1
         y = 0
         z = 0
         do j = 1, 3**deriv
            ! add or set
            if (k==0 .or. k == 2**deriv-1) then
               nucrep(i) = nucrep(i) + merge(-1,1,neg) * power(j)
            else
               nucrep(i) = merge(-1,1,neg) * power(j)
            end if
            ! increment index i
            do s = 0, deriv-1
               ! x -> y transition
               if (.not.btest(y,s) .and. .not.btest(z,s)) then
                  y = ibset(y,s)
                  i = i + ncor_pow(s+1)
                  exit
               ! y -> z transition
               else if (btest(y,s)) then
                  y = ibclr(y,s)
                  z = ibset(z,s)
                  i = i + ncor_pow(s+1)
                  exit
               ! z -> x transition (followed by 'carry')
               else
                  z = ibclr(z,s)
                  i = i - 2*ncor_pow(s+1)
               end if
            end do
         end do
         ! increment offset ofs for next index permutation
         do s = 0, deriv-1
            neg = .not.neg
            if (.not.btest(k,s)) then
               offset = offset + 3*(nucj-nuci) * ncor_pow(s+1)
               exit
            else
               offset = offset - 3*(nucj-nuci) * ncor_pow(s+1)
            end if
         end do
      end do
    end subroutine

  end subroutine


  !> \brief gets the nuclear contribution to electric dipole moment
  !> \author Bin Gao
  !> \date 2009-12-10
  !> \return DN contains the nuclear contribution to electric dipole moment
  !> \note modified on linsca/linears/VIBCTL_interface.F
  subroutine DIPNUC_ifc( DN )
    ! uses MXCOOR
#include "mxcent.h"
    ! uses DIPMN and DDIPN
#include "dipole.h"
    real(8), intent(out) :: DN(3)
    real(8), parameter :: zero = 0.0D+00
    call DIPNUC( (/zero/), (/zero/), 0, .true. )
    DN = DIPMN
  end subroutine


  !> \brief gets the nuclear contribution to quadrupole moments
  !> \author Bin Gao
  !> \date 2009-12-10
  !> \return Q contains the nuclear contribution to quadrupole moments
  !> \note modified on linsca/linears/VIBCTL_interface.F
  subroutine QDRNUC_ifc( Q )
#include "mxcent.h"
#include "nuclei.h"
#ifndef PRG_DIRAC
#include "quadru.h"
#endif
    real(8), intent(out) :: Q(6)
    real(8), parameter :: zero = 0.0D+00
#ifdef PRG_DIRAC
    print *, 'implement qdrnuc_ifc in dirac'
    stop 1
#else
    call NUCQDR( CORD(:,1:NUCDEP), (/zero,zero,zero/), get_print_unit(), 0 )
    Q = (/QDRNUC(1,1),QDRNUC(1:2,2),QDRNUC(1:3,3)/)
#endif
  end subroutine


  !> \brief gets the nuclear contribution to electric dipole gradient
  !> \author Bin Gao
  !> \date 2009-12-10
  !> \param na is the number of atoms
  !> \return DGN contains the nuclear contributions to dipole gradient
  !> \note modified on linsca/linears/VIBCTL_interface.F
  subroutine DPGNUC_ifc( na, DGN )
    ! uses MXCOOR
#include "mxcent.h"
    ! uses DDIPN
#include "dipole.h"
    integer, intent(in) :: na
    real(8), intent(out) :: DGN( 3*na, 3 )
    real(8), parameter :: zero = 0.0D+00
    call DIPNUC( (/zero/), (/zero/), 0, .true. )
    DGN = transpose( DDIPN( :, 1:3*na ) )
  end subroutine


  !> Nuclear contribution to the atomic axial tenaor (AAT),
  !> needed for vibrational circular dichroism (VCD)
  !> In the quasienergy formalism, the AAT is:
  !> d^3E/dR(-w)dB(w)dw |w=0
  subroutine AATNUC_ifc( na, AATN )
#include "mxcent.h"
#include "aatens.h"
    integer, intent(in) :: na
    real(8), intent(out) :: AATN( 3, 3*na )
    real(8) :: CSTRA( 3*na, 3*na ), SCTRA( 3*na, 3*na )
#ifdef PRG_DIRAC
    print *, 'fix nucaat call'
    stop 1
#else
    call NUCAAT( CSTRA, SCTRA, 0 )
         AATN(:,:) = AATNUC( :, :3*na )
#endif
  end subroutine

end module
