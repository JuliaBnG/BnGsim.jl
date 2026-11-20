# The Breeding side of BnGsim (docs/src/upstream-contract.md). BnGsim
# imports Breeding only here.
module BreedingAdapter

using Breeding: Breeding, MME, mme, RandomEffect, solve_mme, incidence, TM1997, feasible, rocs, pev_candidates,
    reml, MultiTrait, mt_stack, mt_blocks, mt_rinv

# The upstream package (for its version).
const PACKAGE = Breeding

export MME, mme, RandomEffect, solve_mme, incidence, TM1997, feasible, rocs, pev_candidates, reml, MultiTrait,
    mt_stack, mt_blocks, mt_rinv

end
