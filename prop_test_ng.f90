!> @file
!> Contains module prop_test

!> ajt/radovan: Response-related testing routines and some calculations 
module prop_test_ng

  use matrix_defop_ng
  use rsp_contribs_ng
  use rsp_equations_ng

  implicit none

  public prop_test_gradient_ng

  private

  !physical constants for conversion to non-atomic units
  real(8), parameter :: cm1 = 1/219474.631371d0, & !1 centimeter-to-minus-one in au
                        cvl = 137.03599907d0,    & !c-the-velocity-of-light in au
                        nm  = 10/0.52917706d0,   & !1 nanometer in au
                        pi  = 3.14159265358979323846D0 !acos(-1.d0)

  !field component lables for printing
  character(2) :: fc(3) = (/'Fx','Fy','Fz'/), &
                  bc(3) = (/'Bx','By','Bz'/)

contains


  !> Calculate the gradient: dE/dR = dh_nuc/dR + Tr dH/dR D + 1/2 Tr dG/dR(D) D
  !>                               + dExc/dR - Tr dS/dR DFD
  subroutine prop_test_gradient_ng(mol, S, D, F)
    type(rsp_cfg), intent(in) :: mol
    type(matrix),  intent(in) :: S, D, F
    complex(8), dimension(3*mol%natom) :: gra, tmp
    type(matrix) :: DFD
    ! first nuclear contribution
    call rsp_nucpot(mol, 1, (/'GEO'/), (/(0d0,0d0)/), (/1/), &
                    shape(tmp), tmp)
    gra = tmp
    ! overlap contribution
    DFD = D*F*D
    call rsp_ovlave(mol, 1, (/'GEO'/), (/(0d0,0d0)/), (/1/), &
                    shape(tmp), D, DFD, tmp)
    DFD = 0 !free
    gra = gra + tmp
    ! 1-electron contribution
    call rsp_oneave(mol, 1, (/'GEO'/), (/1/), shape(tmp), D, tmp)
    gra = gra + tmp
    ! 2-electron contribution
    call rsp_twoave(mol, 1, (/'GEO'/), (/1/), shape(tmp), D, D, tmp)
    ! Kohn-Sham exchange correlation average
    call rsp_ksmave(mol, 1, (/'GEO'/), (/1/), shape(tmp), 1, (/D/), tmp)
    gra = gra + tmp/2
    ! print
    call print_tensor(shape(gra), gra, 'gradient = Eg')
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
