#include "nbody/forces.hpp"

#include <cmath>

namespace nbody {

void compute_accelerations(const System& s, std::vector<Vec3>& acc) {
    const std::size_t n = s.size();
    acc.assign(n, Vec3{});
    const double eps2 = s.softening * s.softening;

    for (std::size_t i = 0; i < n; ++i) {
        const Vec3 ri = s.pos[i];
        Vec3 ai;  // accumulate locally; one write to acc[i] at the end of the row
        for (std::size_t j = i + 1; j < n; ++j) {
            const Vec3 d = s.pos[j] - ri;
            const double r2 = dot(d, d) + eps2;
            const double inv_r = 1.0 / std::sqrt(r2);
            const double inv_r3 = inv_r * inv_r * inv_r;
            // Same scalar for both bodies -> forces are exactly equal and opposite.
            const Vec3 f = (s.G * inv_r3) * d;
            ai += s.mass[j] * f;
            acc[j] -= s.mass[i] * f;
        }
        acc[i] += ai;
    }
}

}  // namespace nbody
