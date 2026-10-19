#pragma once

#include <vector>

#include "nbody/system.hpp"

namespace nbody {

// Direct-summation gravitational accelerations, O(N^2).
//
// Each pair (i, j) is visited once and both accelerations are updated from
// the same displacement vector, which halves the work and makes the
// pairwise forces exactly antisymmetric, so total momentum is conserved to
// round-off regardless of the integrator.
//
// With softening eps the potential is Plummer-softened:
//   a_i = G * sum_j m_j (r_j - r_i) / (|r_j - r_i|^2 + eps^2)^(3/2)
void compute_accelerations(const System& s, std::vector<Vec3>& acc);

}  // namespace nbody
