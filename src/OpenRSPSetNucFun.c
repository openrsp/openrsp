subroutine get_nucpot_contrib(pert_spec, contrib)

undecided type/structure, intent(in)  :: pert_spec ! Specification of perturbations
complex, dimension(num_uniq_p_tuple_components), intent(out) :: contrib ! Resulting contribution: num_uniq_p_tuple_components is size of contribution
