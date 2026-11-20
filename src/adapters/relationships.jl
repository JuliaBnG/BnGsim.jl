# The RelationshipMatrices side of BnGsim (docs/src/upstream-contract.md).
# BnGsim imports RelationshipMatrices only here. `founder_segments` is
# imported so that BnGsim can add its method for a `Population`.
module RelationshipsAdapter

using RelationshipMatrices: RelationshipMatrices, nrm, nrm_diag, ainv, grm, hinv, irm, Meioses,
    HaplotypeSegments, segment_inbreeding
import RelationshipMatrices: founder_segments

# The upstream package (for its version).
const PACKAGE = RelationshipMatrices

export nrm, nrm_diag, ainv, grm, hinv, irm, Meioses, HaplotypeSegments, founder_segments, segment_inbreeding

end
