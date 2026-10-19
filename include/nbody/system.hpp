#pragma once

#include <string>
#include <vector>

#include "nbody/vec3.hpp"

namespace nbody {

// State of an N-body system stored as structure-of-arrays.
//
// Units are chosen by the caller; the gravitational constant G is a field
// so the same code runs in SI, astronomical units, or dimensionless units.
struct System {
    std::vector<double> mass;
    std::vector<Vec3> pos;
    std::vector<Vec3> vel;
    std::vector<std::string> name;

    double G = 1.0;          // gravitational constant in the chosen units
    double softening = 0.0;  // Plummer softening length epsilon (0 = exact Newtonian)
    double time = 0.0;

    std::size_t size() const { return mass.size(); }

    void add_body(double m, const Vec3& r, const Vec3& v, std::string label = "");

    // Shift to the centre-of-mass frame so total momentum is zero and the
    // system does not drift across the domain.
    void to_com_frame();
};

// Diagnostics. All are exact functions of the current state.
double kinetic_energy(const System& s);
double potential_energy(const System& s);
double total_energy(const System& s);
Vec3 total_momentum(const System& s);
Vec3 total_angular_momentum(const System& s);
Vec3 center_of_mass(const System& s);

}  // namespace nbody
