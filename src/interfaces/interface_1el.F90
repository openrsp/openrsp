module interface_1el

   use matrix_defop
   use interface_molecule
   use interface_basis
   use interface_f77_memory
   use interface_pcm
   use interface_io

   implicit none

   public interface_1el_ovlave
   public interface_1el_oneave
   public interface_1el_ovlint
   public interface_1el_oneint

   public oneint_ave
   public get1in_ave_ifc
   public onedrv_ave_ifc
   public di_read_operator_int

   private

contains

   !> \brief host program routine to get the average f-perturbed overlap integrals
   !>        with perturbed density D and energy-weighted density DFD
   subroutine interface_1el_ovlave(nf, f, c, nc, DFD, ave, w, D)
     ! Gen1Int interface
#ifdef VAR_LINSCA
     use gen1int_host
#else
     use gen1int_api
#endif
     !> number of fields
     integer,       intent(in)  :: nf
     !> field labels in std order
     character(4),  intent(in)  :: f(nf)
     !> first and number of- components in each field
     integer,       intent(in)  :: c(nf), nc(nf)
     !> energy-weighted density matrix
     type(matrix),  intent(in)  :: DFD
     !> output average
     complex(8),    intent(out) :: ave(product(nc))
     !> field frequencies corresponding to each field
     complex(8),    intent(in), optional  :: w(nf)
     !> density matrix to contract half-differentiated overlap against
     type(matrix),  intent(in), optional  :: D
     !----------------------------------------------
     type(matrix) A(2)
     real(8), parameter   :: fdistep = 2d0**(-25)
     real(8), allocatable :: tmp(:,:,:)
     real(8), allocatable :: fdi(:,:)
     integer i
     integer order_geo                      !order of total geometric derivatives
     integer num_atom                       !number of atoms
     integer num_coord                      !number of atomic coordinates
     integer num_geom                       !number of total geometric derivatives
     integer num_expt                       !number of all expectation values
     real(8), allocatable :: val_expt(:,:)  !expectation values, real numbers
     integer ierr                           !error information
     if (present(w) .and. .not.present(D)) &
        call quit("error in interface_1el_ovlave: frequencies 'w' and density 'D' " &
               // 'must both be present or both absent')
 
   if (any(f=='EL  ')) then
      ave = 0.0
   else
 
     ! gets the order of total geometric derivatives
     order_geo = count(f=='GEO ')
     if (order_geo/=nf) &
       call quit("interface_1el_ovlave>> only geometric derivatives implemented!")
     ! sets the number of total geometric derivatives
     num_atom = get_nr_atoms()
     num_coord = 3*num_atom
     num_geom = num_coord**order_geo
     ! allocates memory for expectation values
     num_expt = num_geom
     allocate(val_expt(num_expt,1), stat=ierr)
     if (ierr/=0) call quit("interface_1el_ovlave>> failed to allocate val_expt!")
     val_expt = 0.0
 !FIXME changes to call Gen1Int
 !FIXME \sum_{j+k=n} (-\sum_{j}w_{j}+\sum_{k}w_{k})/2 S(>)^{j}(<)^{k}
 ! When it comes to w, I think it's a better idea to ask the integral program for S>> and S<>, and multiply the resulting average or integral by w afterwards, rather than to send w into the integral program.
     ! with field frequencies
     if (present(w)) then
       if (order_geo==1) then
         ! allocate matrices for integrals
         A(1) = mat_alloc_like(DFD)
         if (present(w)) then
            A(2) = mat_alloc_like(DFD)
         end if
         ! loop over nuclear coordinates
         do i = 0, nc(1)-1
            ! (half-) perturbed overlap -i/2 Tg into A(1), Sg in A(2)
            if (present(w)) then !w=0 means no -i/2 Tg contribution
               call di_read_operator_int('SQHDR' // prefix_zeros(c(1)+i,3), A(1))
               A(1) = -A(1) !SQHDR is really -dS>/dg
               A(2) = (-w(1)/2) * (A(1) + trps(A(1)))
               A(1) = A(1) + trps(A(1)) !=1DOVL
            else
               call di_read_operator_int('1DOVL' // prefix_zeros(c(1)+i,3), A(1))
               A(1) = -A(1) !1DOVL is really -dS/dg
            end if
            ave(1+i) = -tr(A(1),DFD)
            if (present(w)) ave(1+i) = ave(1+i) + tr(A(2),D)
         end do
         A(1:2) = 0 !deallocate
       else
         call quit('interface_1el_ovlave>> GEO(>1) with freqencies not implemented!')
       end if
     else
       ! calculates the expectaion values of overlap matrix
       call gen1int_host_get_expt(NON_LAO, INT_OVERLAP,       &
                                  0,                          &  !multipole moments
                                  0, 0, 0,                    &  !magnetic derivatives
                                  0, 0, 0,                    &  !derivatives w.r.t. total RAM
                                  0, 0,                       &  !partial geometric derivatives
                                  min(2,order_geo,num_atom),  &  !total geometric derivatives
                                  order_geo,                  &
                                  0, (/0/),                   &
                                  REDUNDANT_GEO,              &
                                  .false., .false., .false.,  &  !not implemented yet
                                  1, (/DFD/), num_expt,       &  !expectation values
                                  val_expt, .false.,          &
                                  get_print_unit(), 5)
       val_expt = -val_expt
       ! assigns the output average
       if (order_geo==0) then
         ave = val_expt(:,1)
       else
         call gen1int_reorder(num_coord=num_coord, num_field=nf, &
                              first_comp=c, num_comp=nc,         &
                              order_mom=0, order_geo=order_geo,  &
                              val_expect=val_expt(:,1), rsp_expect=ave)
       end if
     end if
     ! frees space
     deallocate(val_expt)
 
   end if
   end subroutine

   !> \brief host program routine to get the average 1-electron integrals perturbed by fields f
   !>        with the (perturbed) density matrix D
   subroutine interface_1el_oneave(nf, f, c, nc, D, ave)
     ! Gen1Int interface
#ifdef VAR_LINSCA
     use gen1int_host
#else
     use gen1int_api
#endif
     !> number of fields
     integer,       intent(in)  :: nf
     !> field labels in std order
     character(4),  intent(in)  :: f(nf)
     !> first and number of- components in each field
     integer,       intent(in)  :: c(nf), nc(nf)
     !> density matrix to average over
     type(matrix),  intent(in)  :: D
     !> output average
     complex(8),    intent(out) :: ave(product(nc))
     !----------------------------------------------
     integer order_mom                      !order of Cartesian multipole moments
     integer num_mom                        !number of Cartesian multipole moments
     integer order_geo                      !order of total geometric derivatives
     integer num_atom                       !number of atoms
     integer num_coord                      !number of atomic coordinates
     integer num_geom                       !number of total geometric derivatives
     integer num_expt                       !number of all expectation values
     real(8), allocatable :: val_expt(:,:)  !expectation values, real numbers
     integer ierr                           !error information
     ! gets the order of Cartesian multipole moments
     order_mom = count(f=='EL  ')
 
   if (order_mom > 1) then
      ave = 0.0
   else
 
     ! gets the order of total geometric derivatives
     order_geo = count(f=='GEO ')
     if (order_mom+order_geo/=nf) &
       call quit("interface_1el_oneave>> only electric and geometric perturbations implemented!")
     ! sets the number of operators and derivatives
     num_mom = (order_mom+1)*(order_mom+2)/2
     num_atom = get_nr_atoms()
     num_coord = 3*num_atom
     num_geom = num_coord**order_geo
     ! allocates memory for expectation values
     num_expt = num_mom*num_geom
     allocate(val_expt(num_expt,1), stat=ierr)
     if (ierr/=0) call quit("interface_1el_oneave>> failed to allocate val_expt!")
     val_expt = 0.0
     ! electric perturbations
     if (order_mom/=0) then
       call gen1int_host_get_expt(NON_LAO, INT_CART_MULTIPOLE, &
                                  order_mom,                   &  !multipole moments
                                  0, 0, 0,                     &  !magnetic derivatives
                                  0, 0, 0,                     &  !derivatives w.r.t. total RAM
                                  0, 0,                        &  !partial geometric derivatives
                                  min(3,order_geo,num_atom),   &  !total geometric derivatives
                                  order_geo,                   &
                                  0, (/0/),                    &
                                  REDUNDANT_GEO,               &
                                  .false., .false., .false.,   &  !not implemented yet
                                  1, (/D/), num_expt,          &  !expectation values
                                  val_expt, .false.,           &
                                  get_print_unit(), 5)
     ! only geometric perturbations
     else
       call gen1int_host_get_expt(NON_LAO, INT_ONE_HAMIL,    &
                                  0,                         &  !multipole moments
                                  0, 0, 0,                   &  !magnetic derivatives
                                  0, 0, 0,                   &  !derivatives w.r.t. total RAM
                                  0, 0,                      &  !partial geometric derivatives
                                  min(3,order_geo,num_atom), &  !total geometric derivatives
                                  order_geo,                 &
                                  0, (/0/),                  &
                                  REDUNDANT_GEO,             &
                                  .false., .false., .false., &  !not implemented yet
                                  1, (/D/), num_expt,        &  !expectation values
                                  val_expt, .false.,         &
                                  get_print_unit(), 5)
     end if
     ! assigns the output average
     call gen1int_reorder(num_coord=num_coord, num_field=nf,        &
                          first_comp=c, num_comp=nc,                &
                          order_mom=order_mom, order_geo=order_geo, &
                          val_expect=val_expt(:,1), rsp_expect=ave)
     ! frees space
     deallocate(val_expt)
 
   end if
   end subroutine

   !> \brief host program routine to compute differentiated overlap matrices, and optionally
   !>        add half-differentiated overlap contribution to Fock matrices
   subroutine interface_1el_ovlint(nr_ao, nf, f, c, nc, ovl, w, fock)
     ! Gen1Int interface
#ifdef VAR_LINSCA
     use gen1int_host
#else
     use gen1int_api
#endif
     integer, intent(in)          :: nr_ao
     !> number of fields
     integer,       intent(in)    :: nf
     !> field labels in std order
     character(4),  intent(in)    :: f(nf)
     !> first and number of- components in each field
     integer,       intent(in)    :: c(nf), nc(nf)
     !> resulting overlap integral matrices (incoming content deleted)
     type(matrix),  intent(inout) :: ovl(product(nc))
     !> frequencies of each field
     complex(8),    intent(in),    optional :: w(nf)
     !> Fock matrices to which the half-differentiated overlap
     !> contribution is ADDED
     type(matrix),  intent(inout), optional :: fock(product(nc))
     !------------------------------------------------
     integer      i
     integer order_geo  !order of total geometric derivatives
     integer num_atom   !number of atoms
     integer num_coord  !number of atomic coordinates
     integer num_geom   !number of total geometric derivatives
     integer num_ints   !number of integral matrices
     if (present(w) .and. .not.present(fock))                                      &
        call quit("error in interface_1el_ovlint: frequencies 'w' and Fock matrix 'fock' "// &
                  "must both be present or both absent")
 
   if (any(f=='EL  ')) then
 
      do i = 1, product(nc)
         call mat_init(ovl(i), nrow=nr_ao, ncol=nr_ao, closed_shell=.true.)
      end do
 
   else
 
     ! gets the order of total geometric derivatives
     order_geo = count(f=='GEO ')
     if (order_geo/=nf) &
       call quit("interface_1el_ovlint>> only geometric derivatives implemented!")
     ! sets the number of total geometric derivatives
     num_atom = get_nr_atoms()
     num_coord = 3*num_atom
     num_geom = num_coord**order_geo
     ! sets the number of integral matrices
     num_ints = num_geom
     if (num_ints/=size(ovl)) &
       call quit("interface_1el_ovlint>> returning specific components not implemented!")
 !FIXME changes to call Gen1Int
     ! with field frequencies
     if (present(w)) then
       if (order_geo==1) then
         ! loop over nuclear coordinates
         do i = 0, nc(1)-1
            ! allocate, if needed
            if (.not.isdef(ovl(1+i))) then
               call mat_init(ovl(1+i), nrow=nr_ao, ncol=nr_ao, closed_shell=.true.)
            end if
            ! overlap into ovl, half-perturbed overlap -i/2 Tg added to fock
            if (present(w)) then
               call di_read_operator_int('SQHDR' // prefix_zeros(c(1)+i,3), ovl(1+i))
               ovl(1+i)  = -ovl(1+i) !SQHDR is really -dS>/dg
               fock(1+i) = fock(1+i) - w(1)/2 * ovl(1+i)
               fock(1+i) = fock(1+i) + w(1)/2 * trps(ovl(1+i))
               ovl(1+i)  = ovl(1+i)  + trps(ovl(1+i)) !=dS/dg=-1DOVL
            else
               call di_read_operator_int('1DOVL' // prefix_zeros(c(1)+i,3), ovl(1+i))
               ovl(1+i) = -ovl(1+i) !1DOVL is really -dS/dg
            end if
         end do
 !FIXME to Andreas: do we need higher order geometric derivatives of overlap integrals with frequencies?
       else
         call quit('interface_1el_ovlint>> GEO(>1) with freqencies not implemented!')
       end if
     else
       ! allocates matrices
       do i = 1, num_ints
         if (.not.isdef(ovl(i))) then
           call mat_init(ovl(i), nrow=nr_ao, ncol=nr_ao, closed_shell=.true.)
         end if
       end do
       ! calculates the overlap matrix
       call gen1int_host_get_int(NON_LAO, INT_OVERLAP,       &
                                 0,                          &  !multipole moments
                                 0, 0, 0,                    &  !magnetic derivatives
                                 0, 0, 0,                    &  !derivatives w.r.t. total RAM
                                 0, 0,                       &  !partial geometric derivatives
                                 min(2,order_geo,num_atom),  &  !total geometric derivatives
                                 order_geo,                  &
                                 0, (/0/),                   &
                                 REDUNDANT_GEO,              &
                                 .false., .false., .false.,  &  !not implemented yet
                                 num_ints, ovl, .false.,     &  !integral matrices
                                 get_print_unit(), 5)
     end if
 
   end if
   end subroutine

   subroutine interface_1el_oneint(nr_ao, nf, f, c, nc, oneint)
     ! Gen1Int interface
#ifdef VAR_LINSCA
     use gen1int_host
#else
     use gen1int_api
#endif
     integer, intent(in)          :: nr_ao
     !> number of fields
     integer,       intent(in)    :: nf
     !> field labels in std order
     character(4),  intent(in)    :: f(nf)
     !> first and number of- components in each field
     integer,       intent(in)    :: c(nf), nc(nf)
     !> output perturbed integrals
     type(matrix),  intent(inout) :: oneint(product(nc))
     !--------------------------------------------------
     integer order_mom  !order of Cartesian multipole moments
     integer num_mom    !number of Cartesian multipole moments
     integer order_geo  !order of total geometric derivatives
     integer num_atom   !number of atoms
     integer num_coord  !number of atomic coordinates
     integer num_geom   !number of total geometric derivatives
     integer num_ints   !number of all integral matrices
     integer imat       !incremental recorder over matrices
     integer :: i
     type(matrix) :: A
 
   if (count(f=='EL  ') > 1) then
      call mat_init(A, nrow=nr_ao, ncol=nr_ao, closed_shell=.true.)
      do i = 1, product(nc)
         if (iszero(oneint(i))) then
            call mat_ensure_alloc(oneint(i))
            oneint(i)%elms = oneint(i)%elms + A%elms
         else
            oneint(i)%elms = oneint(i)%elms + A%elms
         end if
      end do
 
   else
 
     ! gets the order of Cartesian multipole moments
     order_mom = count(f=='EL  ')
     ! gets the order of total geometric derivatives
     order_geo = count(f=='GEO ')
     if (order_mom+order_geo/=nf) &
       call quit("interface_1el_oneint>> only electric and geometric perturbations implemented!")
     ! sets the number of operators and derivatives
     num_mom = (order_mom+1)*(order_mom+2)/2
     num_atom = get_nr_atoms()
     num_coord = 3*num_atom
     num_geom = num_coord**order_geo
     num_ints = num_mom*num_geom
     if (num_ints/=size(oneint)) &
       call quit("interface_1el_oneint>> returning specific components is not implemented!")
 !FIXME: it is better that we use unique components for higher order
     if (order_mom>1) &
       call quit("interface_1el_oneint>> only the first Cartesian multipole moments implemented!")
     do imat = 1, size(oneint)
       if (.not.isdef(oneint(imat))) then
         call mat_init(oneint(imat), nrow=nr_ao, ncol=nr_ao, closed_shell=.true.)
       end if
     end do
     ! electric perturbations
     if (order_mom/=0) then
       call gen1int_host_get_int(NON_LAO, INT_CART_MULTIPOLE, &
                                 order_mom,                   &  !multipole moments
                                 0, 0, 0,                     &  !magnetic derivatives
                                 0, 0, 0,                     &  !derivatives w.r.t. total RAM
                                 0, 0,                        &  !partial geometric derivatives
                                 min(3,order_geo,num_atom),   &  !total geometric derivatives
                                 order_geo,                   &
                                 0, (/0/),                    &
                                 REDUNDANT_GEO,               &
                                 .false., .false., .false.,   &  !not implemented yet
                                 num_ints, oneint, .false.,   &  !integral matrices
                                 get_print_unit(), 5)
 
     ! only geometric perturbations
     else
       call gen1int_host_get_int(NON_LAO, INT_ONE_HAMIL,    &
                                 0,                         &  !multipole moments
                                 0, 0, 0,                   &  !magnetic derivatives
                                 0, 0, 0,                   &  !derivatives w.r.t. total RAM
                                 0, 0,                      &  !partial geometric derivatives
                                 min(3,order_geo,num_atom), &  !total geometric derivatives
                                 order_geo,                 &
                                 0, (/0/),                  &
                                 REDUNDANT_GEO,             &
                                 .false., .false., .false., &  !not implemented yet
                                 num_ints, oneint, .false., &  !integral matrices
                                 get_print_unit(), 5)
     end if
 
   end if
   end subroutine

  !> \brief reorders and assigns the expectation values and/or integral matrices from
  !>        Gen1int to OpenRsp
  !> \author Bin Gao
  !> \date 2012-04-25
  !> \param num_coord is the number of atomic coordinates
  !> \param num_field is the number of fields
  !> \param first_comp constains the first component in each field ("GEO", "EL")
  !> \param num_comp contains the number of components in each field ("GEO", "EL")
  !> \param order_mom is the order of Cartesian multipole moments
  !> \param order_geo is the order of total geometric derivatives
  !> \param val_expect contains the redudant expectation values from Gen1Int
  !> \param val_ints contains the redudant integral matrices from Gen1Int
  !> \return rsp_expect contains the redudant expectation values for OpenRsp
  !> \return rsp_ints contains the redudant integral matrices for OpenRsp
  !> \note all the expectation values and integral matrices are in the order ("EL", "GEO")
  !>       in memory; this routine is implemented temporarily since it will cost much memory
  !>       for assigning integral matrices, it would be better that the integral
  !>       code returns the required integral matrices by itself
  subroutine gen1int_reorder(num_coord, num_field, first_comp, num_comp, &
                             order_mom, order_geo, val_expect, val_ints, &
                             rsp_expect, rsp_ints)
    integer, intent(in) :: num_coord
    integer, intent(in) :: num_field
    integer, intent(in) :: first_comp(num_field)
    integer, intent(in) :: num_comp(num_field)
    integer, intent(in) :: order_mom
    integer, intent(in) :: order_geo
    real(8), optional, intent(in)  :: val_expect(:)
    type(matrix), optional, intent(in) :: val_ints(:)
    complex(8), optional, intent(out) :: rsp_expect(product(num_comp))
    type(matrix), optional, intent(inout) :: rsp_ints(product(num_comp))
    logical assign_expect                   !if assigning expectation values
    logical assign_ints                     !if assigning integral matrices
    integer size_rsp_comp                   !size of all components for OpenRsp
    integer, allocatable :: offset_geom(:)  !offset of total geometric derivatives
    integer, allocatable :: last_comp(:)    !last component in each field
    integer, allocatable :: idx_comp(:)     !indices of specific components
    integer addr_rsp                        !address of specific components in the results for OpenRsp
    integer addr_gen1int                    !address of specific components in the results from Gen1Int
    integer icomp                           !incremental recorder over components
    integer ymom, zmom                      !orders of yz components of Cartesian multipole moments
    logical not_found                       !indicates the next specific components are not found
    integer ierr                            !error information
    ! sets what jobs are going to do
    assign_expect = present(val_expect) .and. present(rsp_expect)
    assign_ints = present(val_ints) .and. present(rsp_ints)
    ! sets the size of all components in the fields
    if (assign_expect) then
      size_rsp_comp = size(rsp_expect)
      if (assign_ints) then
        if (size_rsp_comp/=size(rsp_ints)) &
          call quit("gen1int_reorder>> requires the same size of expectation "// &
                    "values and integral matrices!")
      end if
    else if (assign_ints) then
      size_rsp_comp = size(rsp_ints)
    else
      return
    end if
    ! sets the offset of total geometric derivatives if there are
    if (order_geo>0) then
      allocate(offset_geom(order_geo), stat=ierr)
      if (ierr/=0) call quit("gen1int_reorder>> failed to allocate offset_geom!")
      offset_geom(order_geo) = (order_mom+1)*(order_mom+2)/2  !Cartesian multipole moments come first in memory
      do icomp = order_geo-1, 1, -1
        offset_geom(icomp) = num_coord*offset_geom(icomp+1)
      end do
    end if
    ! sets the last component in each field
    allocate(last_comp(num_field), stat=ierr)
    if (ierr/=0) call quit("gen1int_reorder>> failed to allocate last_comp!")
    do icomp = 1, num_field
      last_comp(icomp) = first_comp(icomp)+num_comp(icomp)-1
    end do
    ! indices of specific components
    allocate(idx_comp(num_field), stat=ierr)
    if (ierr/=0) call quit("gen1int_reorder>> failed to allocate idx_comp!")
    idx_comp = first_comp
    ! gets the orders of yz components of Cartesian multipole moments
    ymom = 0
    zmom = 0
    do icomp = num_field, order_geo+1, -1
      select case(idx_comp(icomp))
      case(1)
      case(2)
        ymom = ymom+1
      case(3)
        zmom = zmom+1
      case default
        call quit("gen1int_reorder>> unknown components for Cartesian multipole moments!")
      end select
    end do
    ! the address of x^{l}y^{m}z^{n} with l+m+n=order_mom is 1+m+(2*order_mom+3-n)*n/2
    addr_gen1int = 1+ymom+(2*order_mom+3-zmom)*zmom/2
    ! gets the address of total geometric derivatives in the results from Gen1Int
    do icomp = order_geo, 1, -1
      addr_gen1int = addr_gen1int+(idx_comp(icomp)-1)*offset_geom(icomp)
    end do
    ! assigns the expectation values
    if (assign_expect) rsp_expect(1) = val_expect(addr_gen1int)
    ! assigns the integral matrices
    if (assign_ints) rsp_ints(1) = val_ints(addr_gen1int)
    ! sets other components
    do addr_rsp = 2, size_rsp_comp
      ! generates new specific components
      not_found = .true.
      do icomp = num_field, 1, -1
        ! we walk to the next component in this field, the new specific components are found
        if (idx_comp(icomp)<last_comp(icomp)) then
          idx_comp(icomp) = idx_comp(icomp)+1
          not_found = .false.
          exit
        ! we have walked to the last component of this field, back to the first component
        else
          idx_comp(icomp) = first_comp(icomp)
        end if
      end do
      ! we could not find new specific components
      if (not_found) call quit("gen1int_reorder>> can not find next components!")
      ! gets the orders of yz components of Cartesian multipole moments
      ymom = 0
      zmom = 0
    do icomp = num_field, order_geo+1, -1
        select case(idx_comp(icomp))
        case(1)
        case(2)
          ymom = ymom+1
        case(3)
          zmom = zmom+1
        case default
          call quit("gen1int_reorder>> unknown components for Cartesian multipole moments!")
        end select
      end do
      ! the address of x^{l}y^{m}z^{n} with l+m+n=order_mom is 1+m+(2*order_mom+3-n)*n/2
      addr_gen1int = 1+ymom+(2*order_mom+3-zmom)*zmom/2
      ! gets the address of total geometric derivatives in the results from Gen1Int
      do icomp = order_geo, 1, -1
        addr_gen1int = addr_gen1int+(idx_comp(icomp)-1)*offset_geom(icomp)
      end do
      ! assigns the expectation values
      if (assign_expect) rsp_expect(addr_rsp) = val_expect(addr_gen1int)
      ! assigns the integral matrices
      if (assign_ints) rsp_ints(addr_rsp) = val_ints(addr_gen1int)
    end do
    ! frees space
    if (allocated(offset_geom)) deallocate(offset_geom)
    deallocate(last_comp)
    deallocate(idx_comp)
  end subroutine

   subroutine oneint_ave(nr_atoms, what, D, DFD, R)
      integer,           intent(in)  :: nr_atoms
      character(*),      intent(in)  :: what
      type(matrix),      intent(in)  :: D, DFD
      real(8),           intent(out) :: R(:)
#ifndef LSDALTON_ONLY
#include "mxcent.h"
#include "taymol.h"
#endif
      real(8), pointer :: wrk(:)
      integer          :: lwrk
#ifdef LSDALTON_ONLY
      call quit('Cannot run oneint_ave, only new integral code is compiled',-1)
#else
      call save_D_and_DFD_for_ABACUS(.false., D, DFD)
      lwrk = 50*D%nrow**2+10000*D%nrow+50000000
      call f77_memory_select(work_len=lwrk, work=wrk)
      HESMOL(:3*nr_atoms,:3*nr_atoms) = 0
      ! SUBROUTINE ONEDRV(WORK,LWORK,IPRINT,PROPTY,MAXDIF,DIFINT,NODC,
      ! &                  NODV,DIFDIP,HFONLY,NCLONE)
      call ONEDRV(wrk,lwrk,0,.true.,len(what),.true.,.true., &
                  .true.,.false.,.true.,.false.)
      ! HESMOL will contain either HESSKE+HESSNA or HESFS2 or their sum
      R(1:9*nr_atoms**2) = reshape(HESMOL(:3*nr_atoms,:3*nr_atoms), (/9*nr_atoms**2/))
      call f77_memory_deselect(work_len=lwrk, work=wrk)
#endif
   end subroutine

   subroutine save_D_and_DFD_for_ABACUS(anti, D, DFD)
   !ajt aug09 Same as set_dsofso in fock-eval.f90, but pertrubed
   !          D and DFD is input instead of unperturbed D and F
      logical,      intent(in) :: anti !whether the integrals D and DFD
         ! are to be contracted with are symmetric or anti-symmetric
      type(matrix), intent(in) :: D, DFD
         ! perturbed (or un-) density and
         ! 'generalized Fock matrix' (energy-weighted density matrix)
      real(8)   Dtri(D%nrow*(D%nrow+1)/2)
      real(8) DFDtri(D%nrow*(D%nrow+1)/2)
      if (iszero(D)) then
         Dtri = 0
      else
         if (.not.anti) call DGEFSP(D%nrow, D%elms, Dtri)
         if (     anti) call DGETAP(D%nrow, D%elms, Dtri)
         ! scale elms by 4 if anti, 2 if symm
         Dtri = Dtri * merge(4,2,anti)
      end if
      if (iszero(DFD)) then
         DFDtri = 0
      else
         if (.not.anti) call DGEFSP(D%nrow, DFD%elms, DFDtri)
         if (     anti) call DGETAP(D%nrow, DFD%elms, DFDtri)
         ! scale elms by 4 if anti, 2 if symm
         DFDtri = DFDtri * merge(4,2,anti)
      end if
      ! write fo files
#if defined(LSDALTON_ONLY) || defined(PRG_DIRAC)
      call quit('Cannot call write_dsofso, only new integral code is compiled',-1)
#else
      call write_dsofso(Dtri,DFDtri)
#endif
   end subroutine

  subroutine ONEDRV_ave_ifc(fld, siz, ave, D, DFD)
    character(4),           intent(in)  :: fld(:)
    integer,                intent(in)  :: siz
    real(8),                intent(out) :: ave(siz)
    type(matrix), optional, intent(in)  :: D, DFD
    !--------------------------------------------
#include "mxcent.h"
#include "taymol.h"

    real(8), allocatable :: Dtri(:)
    real(8), allocatable :: DFDtri(:)

    logical anti
    integer nc

    allocate(Dtri(  get_nr_ao()*(get_nr_ao()+1)/2))
    allocate(DFDtri(get_nr_ao()*(get_nr_ao()+1)/2))

    ! create triangularly packed matrices from D, DFD
    anti = (mod(count(fld=='MAG '),2) == 1)
    if (.not.present(D)) then
       Dtri = 0
    else if (anti) then
       call DGETAP(D%nrow, D%elms, Dtri)
    else !symm
       call DGEFSP(D%nrow, D%elms, Dtri)
    end if
    if (.not.present(DFD)) then
       DFDtri = 0
    else if (anti) then
       call DGETAP(DFD%nrow, DFD%elms, DFDtri)
    else !symm
       call DGEFSP(DFD%nrow, DFD%elms, DFDtri)
    end if
    ! write to files
#ifdef PRG_DIRAC
    print *, 'fix WRITE_DSOFSO'
    stop 1
#else
    call WRITE_DSOFSO(Dtri, DFDtri)
#endif
    nc = 3 * get_nr_atoms()
    HESMOL(:nc,:nc) = 0
    !  SUBROUTINE ONEDRV(WORK,LWORK,IPRINT,PROPTY,MAXDIF,
    ! &                  DIFINT,NODC,NODV,DIFDIP,DIFQDP,
    ! &                  HFONLY,NCLONE,PCM)
    call ONEDRV(f77_memory, size(f77_memory), 5, .true., size(fld), &
                .true., .true., .true., .false., .false., &
                .true., .false., .false.)
    ! HESMOL will contain either HESSKE+HESSNA or HESFS2 or their sum
    if (size(fld)==1) &
       call quit('error in ONEDRV_ave_ifc: not implemented')
    if (size(fld)==2) &
       ave = reshape(2*HESMOL(:nc,:nc), (/nc*nc/)) !factor 2 for total dens
    deallocate(Dtri)
    deallocate(DFDtri)
  end subroutine

  !> Call GET1IN in ABACUS
  subroutine GET1IN_ave_ifc(fld, siz, ave, D)
    character(4),           intent(in)  :: fld(:)
    integer,                intent(in)  :: siz
    real(8),                intent(out) :: ave(siz)
    type(matrix),           intent(in)  :: D
    !--------------------------------------------
    real(8)                             :: dummy(2)
#include "mxcent.h"
    real(8), allocatable :: Dtri(:)
    character(8), dimension(9*MXCENT) :: labint
    integer, dimension(9*MXCENT) :: intrep, intadr
    integer ncomp
    allocate(Dtri(  get_nr_ao()*(get_nr_ao()+1)/2))
    !ajt Dzero is dangerous! Should have been izero. =0 safe
    intrep = 0 !call dzero(intrep,9*MXCENT)
    intadr = 0 !call dzero(intadr,9*MXCENT)
    ! create triangularly packed matrix from D
    call DGEFSP(D%nrow, D%elms, Dtri)
    ncomp = 0
!      SUBROUTINE GET1IN(SINTMA,WORD,NCOMP,WORK,LWORK,LABINT,INTREP,
!     &                  INTADR,MPQUAD,TOFILE,KPATOM,TRIMAT,EXPVAL,
!     &                  EXP1VL,DENMAT,NPRINT)
#ifdef PRG_DIRAC
    print *, 'fix get1in call'
    stop 1
#else
    call GET1IN(dummy,'DPLGRA ',ncomp,f77_memory,size(f77_memory),labint,intrep, &
                       intadr,0,.false.,0,.false.,ave, &
                       .true.,Dtri,0)
#endif
    deallocate(Dtri)
  end subroutine


  !> \brief gets property integrals from file AOPROPER
  !> \author Bin Gao
  !> \date 2009-12-08
  !> \param prop_lab is the label of integral
  !> \param init_prop indicates if initialize the integral matrix
  !> \return prop_int contains the integral matrix
  subroutine di_read_operator_int( prop_lab, prop_int )
    character*(8), intent(in) :: prop_lab
    type(matrix), intent(inout) :: prop_int
    ! uses NBAST, NNBAST, NNBASX, N2BASX
#include "inforb.h"
    ! IO units, use LUPROP for file AOPROPER
#include "inftap.h"
    ! external DALTON function finding the corresponding label
    logical FNDLB2
    ! information when calling subroutine FNDLB2
    character*8 RTNLBL(2)
    ! dummy stuff
    integer IDUMMY
    logical anti
    ! one-electron Hamiltonian
    if ( prop_lab == 'ONEHAMIL' ) then
      call QUIT( 'Not implemented!' )
#ifdef PRG_DIRAC
    print *, 'fix rdonel call'
    stop 1
#else
      anti = .false. !radovan: was undefined, setting it to false
      call RDONEL( 'ONEHAMIL', ANTI, f77_memory(get_f77_memory_next()), NNBASX )
#endif
      call DSPTSI( NBAST, f77_memory(get_f77_memory_next()), prop_int%elms )
    else
      ! closes file AOPROPER first
      if ( LUPROP > 0 ) call GPCLOSE( LUPROP, 'KEEP' )
      call GPOPEN( LUPROP, 'AOPROPER', 'OLD', ' ', 'UNFORMATTED', IDUMMY, .false. )
      ! finds the label
      if ( FNDLB2( prop_lab, RTNLBL, LUPROP ) ) then
        ! square matrix
        if ( RTNLBL(2) == 'SQUARE' ) then
          call READT( LUPROP, N2BASX, prop_int%elms )
        ! symmetric matrix
        else if ( RTNLBL(2) == 'SYMMETRI' ) then
          call READT( LUPROP, NNBASX, f77_memory(get_f77_memory_next()) )
          call DSPTSI( NBAST, f77_memory(get_f77_memory_next()), prop_int%elms )
        ! anti-symmetric matrix
        else if ( RTNLBL(2) == 'ANTISYMM' ) then
          call READT( LUPROP, NNBASX, f77_memory(get_f77_memory_next()) )
          call DAPTGE( NBAST, f77_memory(get_f77_memory_next()), prop_int%elms )
        else
          call QUIT( 'Error: No symmetry label on AOPROPER!' )
        end if
        ! closes file AOPROPER
        call GPCLOSE( LUPROP, 'KEEP' )
      else
        call GPCLOSE( LUPROP, 'KEEP' )
        call QUIT( 'Integrals with label '''//prop_lab// &
                   ''' not found on file ''AOPROPER''!' )
      end if
    end if
  end subroutine

  function prefix_zeros(n, l)
    integer, intent(in) :: n, l !number, length
    character(l)        :: prefix_zeros !resulting n in ascii
    character(1), parameter :: char0to9(0:9) &
          = (/'0','1','2','3','4','5','6','7','8','9'/)
    integer :: i, k
    k = n
    do i = l, 1, -1
       prefix_zeros(i:i) = char0to9(mod(k,10))
       k = k / 10
    end do
    if (k /= 0) call quit('prefix_zeros error: Argument integer does not fit ' &
                       // 'in the specified number of ASCII caracters')
  end function

end module
