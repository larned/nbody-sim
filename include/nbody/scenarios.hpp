#pragma once

#include <string>
#include <vector>

#include "nbody/system.hpp"

namespace nbody {

// Built-in initial conditions. Each returns a System already shifted to the
// centre-of-mass frame.

// Two-body Kepler orbit: a light body (mass m_light) around a heavy one
// (mass 1), semi-major axis a, eccentricity e, G = 1. Starts at periapsis.
// The period is 2*pi*sqrt(a^3 / (G*(m1+m2))).
System kepler(double eccentricity = 0.0, double semi_major_axis = 1.0, double m_light = 1e-3);

// The Chenciner-Montgomery figure-eight: three equal masses chasing each
// other along a single figure-eight curve (G = 1, m = 1). Period ~6.3259.
System figure_eight();

// Sun plus the eight planets on circular, coplanar orbits at their mean
// distances. Units: AU, years, solar masses (G = 4*pi^2).
System solar_system();

// N bodies drawn from a Plummer sphere (dimensionless units, total mass 1,
// G = 1). Deterministic for a given seed. Used for performance benchmarks.
System plummer(std::size_t n, unsigned seed = 42, double softening = 0.01);

// Look up a scenario by name: "kepler", "figure8", "solar", "plummer".
// For "plummer", n sets the number of bodies.
System make_scenario(const std::string& name, std::size_t n = 256);

std::vector<std::string> scenario_names();

// Analytic Kepler orbit for testing: given the initial state produced by
// kepler(e, a, m), return the exact relative position (light - heavy) at
// time t, by solving Kepler's equation with Newton's method.
Vec3 kepler_exact_relative_position(double eccentricity, double semi_major_axis, double total_mass,
                                    double t);

// Same, for the relative velocity.
Vec3 kepler_exact_relative_velocity(double eccentricity, double semi_major_axis, double total_mass,
                                    double t);

}  // namespace nbody
