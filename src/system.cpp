#include "nbody/system.hpp"

#include <cmath>
#include <utility>

namespace nbody {

void System::add_body(double m, const Vec3& r, const Vec3& v, std::string label) {
    mass.push_back(m);
    pos.push_back(r);
    vel.push_back(v);
    name.push_back(std::move(label));
}

void System::to_com_frame() {
    double m_tot = 0.0;
    Vec3 r_cm, v_cm;
    for (std::size_t i = 0; i < size(); ++i) {
        m_tot += mass[i];
        r_cm += mass[i] * pos[i];
        v_cm += mass[i] * vel[i];
    }
    if (m_tot <= 0.0) return;
    r_cm *= 1.0 / m_tot;
    v_cm *= 1.0 / m_tot;
    for (std::size_t i = 0; i < size(); ++i) {
        pos[i] -= r_cm;
        vel[i] -= v_cm;
    }
}

double kinetic_energy(const System& s) {
    double ke = 0.0;
    for (std::size_t i = 0; i < s.size(); ++i) ke += 0.5 * s.mass[i] * dot(s.vel[i], s.vel[i]);
    return ke;
}

double potential_energy(const System& s) {
    // Must match the softened force law in forces.cpp, otherwise "energy
    // conservation" is measured against the wrong Hamiltonian.
    const double eps2 = s.softening * s.softening;
    double pe = 0.0;
    for (std::size_t i = 0; i < s.size(); ++i) {
        for (std::size_t j = i + 1; j < s.size(); ++j) {
            const Vec3 d = s.pos[j] - s.pos[i];
            pe -= s.G * s.mass[i] * s.mass[j] / std::sqrt(dot(d, d) + eps2);
        }
    }
    return pe;
}

double total_energy(const System& s) { return kinetic_energy(s) + potential_energy(s); }

Vec3 total_momentum(const System& s) {
    Vec3 p;
    for (std::size_t i = 0; i < s.size(); ++i) p += s.mass[i] * s.vel[i];
    return p;
}

Vec3 total_angular_momentum(const System& s) {
    Vec3 L;
    for (std::size_t i = 0; i < s.size(); ++i) L += s.mass[i] * cross(s.pos[i], s.vel[i]);
    return L;
}

Vec3 center_of_mass(const System& s) {
    double m_tot = 0.0;
    Vec3 r;
    for (std::size_t i = 0; i < s.size(); ++i) {
        m_tot += s.mass[i];
        r += s.mass[i] * s.pos[i];
    }
    return m_tot > 0.0 ? r * (1.0 / m_tot) : r;
}

}  // namespace nbody
