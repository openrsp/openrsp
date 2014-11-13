subroutine get_xc_contrib(pert_spec, num_mat_in_tuple, num_components, dens_collection, contrib)

undecided type/structure, intent(in)  :: pert_spec ! Specification of perturbations for property
integer, intent(in)  :: num_mat_in_tuple ! Number of (perturbed) density matrices in tuple, e.g. (/D, D(i), D(j), D(ij)/) gives num_mat_in_tuple = 4
integer, intent(in)  :: num_components ! Number of components to get result for
undecided type/structure, intent(in) :: pert_spec_tuple ! Specification of which perturbations the density matrices in tuple are differentiated with respect to
matrix, pointer, dimension(num_mat_in_tuple, num_components), intent(in)  :: dens_collection ! Collection of (perturbed) density matrixes needed by XC code
complex, dimension(num_components), intent(out)  :: contrib ! Resulting contribution
