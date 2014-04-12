module interface_nuclear

   use interface_molecule
   use interface_io

   implicit none

   public nuclear_potential
   public gradient_nuc
   public hessian_nuc
   public aatnuc_ifc
   public get_nuc_dipole_gradient
#ifndef VAR_LSDALTON
   public dipnuc_ifc
   public qdrnuc_ifc
#endif
   private

contains

  !> (derivatives of) nuclear repulsion and nuclei--field interaction
  subroutine nuclear_potential(deriv, ncor, field, ncomp, nucpot)
    !> order of geometry derivative requested
    integer,      intent(in)  :: deriv
    !> number of geometry coordinates
    integer,      intent(in)  :: ncor
    !> external potential labels, or 'NONE'
    character(4), intent(in)  :: field(2)
    !> number of components of the field (must be all, ie. 1/3/6)
    integer,      intent(in)  :: ncomp
    !> output tensor
    real(8),      intent(out) :: nucpot(ncor**deriv * ncomp)
    !-------------------------------------------------------
    integer i, j, na, charges(ncor/3)
    real(8) coords(ncor)
    ! early return if deriv higher than moment (zero)
    if ((field(1) == 'EL  ' .and. deriv > 1) .or. &
        (field(1) == 'ELGR' .and. deriv > 2)) then
        nucpot(:) = 0
        return
    end if
    ! nonzero, so fetch charges and coordinates
    na = ncor/3
    charges(:) = (/(get_nuc_charge(i), i=1,na)/)
    coords(:)  = (/((get_nuc_xyz(i,j), i=1,3), j=1,na)/)
    ! ---- switch by field(1), field(2) and deriv
    ! all orders of nuclear repulsion
    if (field(1) == 'NONE') then
       call nucrep_deriv(ncor/3, charges, coords, deriv, nucpot)
    ! minus dipole moment
    else if (field(1) == 'EL  ' .and. deriv == 0) then
       call DIPNUC_ifc(nucpot)
       nucpot = -nucpot !dip = -dE/dF
    ! minus dipole moment gradient
    else if (field(1) == 'EL  ' .and. deriv == 1) then
       call get_nuc_dipole_gradient(na, nucpot)
       nucpot = -nucpot !d/dR dip = -d2E/dR/dF
    ! minus quadrupole moment
    else if (field(1) == 'ELGR' .and. deriv == 0) then
       call QDRNUC_ifc(nucpot)
       nucpot = -nucpot !qua = -dE/dG
    ! atomic axial tensor
    else if (field(1) == 'MAG ' .and. field(2) == 'NONE' .and. deriv == 1) then
       call AATNUC_ifc(na, nucpot)
       ! change sign, transpose
       nucpot(:) = (/-nucpot(1::3), -nucpot(2::3), -nucpot(3::3)/)
    else
       print *,'rsp_backend nuclear_potential error: not implented, deriv =', &
               deriv, '  field = ', field
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

#ifndef VAR_LSDALTON
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
#endif


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


   subroutine get_nuc_dipole_gradient(nr_atoms, dipole_gradient)

      integer, intent(in)  :: nr_atoms
      real(8), intent(out) :: dipole_gradient(3*nr_atoms, 3)

      integer              :: iatom
      integer              :: ixyz
      integer              :: n

      dipole_gradient = 0.0d0
      n = 0
      do iatom = 1, nr_atoms
         do ixyz = 1, 3
            n = n + 1
            dipole_gradient(n, ixyz) = get_nuc_charge(iatom)
         end do
      end do

   end subroutine


   !> Nuclear contribution to the atomic axial tenaor (AAT),
   !> needed for vibrational circular dichroism (VCD)
   !> In the quasienergy formalism, the AAT is:
   !> d^3E/dR(-w)dB(w)dw |w=0
   subroutine aatnuc_ifc(na, tensor)

      integer, intent(in)  :: na
      real(8), intent(out) :: tensor(3, 3*na)

      integer :: k, l, ia, ib
      integer :: lc(3, 3)
      real(8) :: flc(3, 3)

      lc = 0
      lc(1, 2) = 3
      lc(2, 1) = 3
      lc(1, 3) = 2
      lc(3, 1) = 2
      lc(2, 3) = 1
      lc(3, 2) = 1

      flc = 0.0d0
      flc(1, 2) =  1.0d0
      flc(2, 1) = -1.0d0
      flc(1, 3) = -1.0d0
      flc(3, 1) =  1.0d0
      flc(2, 3) =  1.0d0
      flc(3, 2) = -1.0d0

      tensor = 0.0d0
      l = 0
      do k = 1, get_num_atoms()
         do ia = 1, 3
            l = l + 1
            do ib = 1, 3
               tensor(ib, l) = get_nuc_charge(k)*get_nuc_xyz(lc(ia, ib), k)*flc(ia, ib)
            end do
         end do
      end do
      tensor = 0.25d0*tensor

   end subroutine

end module
