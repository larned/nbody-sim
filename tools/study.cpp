// nbody_study: generates the data behind the plots in the README.
//
//   nbody_study convergence OUT.csv   error after one Kepler period vs dt
//   nbody_study longterm    OUT.csv   energy error over many Kepler periods
//   nbody_study scaling     OUT.csv   force-evaluation time vs number of bodies

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "nbody/forces.hpp"
#include "nbody/integrators.hpp"
#include "nbody/scenarios.hpp"
#include "nbody/system.hpp"

namespace {

using namespace nbody;

constexpr double kPi = 3.14159265358979323846;
constexpr double kEcc = 0.5;
constexpr double kMassLight = 1e-3;

double kepler_period() { return 2.0 * kPi / std::sqrt(1.0 + kMassLight); }

Vec3 relative_position(const System& s) { return s.pos[1] - s.pos[0]; }

// Integrate 0.6 of an orbit rather than a whole one: after exactly one period
// symplectic Euler's boundary half-kicks cancel and it looks second order.
constexpr double kStudyFraction = 0.6;

int convergence(const std::string& path) {
    std::ofstream out(path);
    out << "integrator,order,dt,steps,force_evals,position_error,phase_space_error,energy_error\n"
        << std::setprecision(12);
    const double T = kepler_period();
    for (const auto& name : integrator_names()) {
        // Steps per period from 2^6 to 2^15.
        for (int p = 6; p <= 15; ++p) {
            const long steps = 1L << p;
            const double dt = kStudyFraction * T / static_cast<double>(steps);
            System s = kepler(kEcc, 1.0, kMassLight);
            auto integ = make_integrator(name);
            const double E0 = total_energy(s);
            for (long k = 0; k < steps; ++k) integ->step(s, dt);
            const double m_tot = 1.0 + kMassLight;
            const Vec3 r = kepler_exact_relative_position(kEcc, 1.0, m_tot, s.time);
            const Vec3 v = kepler_exact_relative_velocity(kEcc, 1.0, m_tot, s.time);
            const double pos_err = norm(relative_position(s) - r);
            const double ps_err = pos_err + norm((s.vel[1] - s.vel[0]) - v);
            const double e_err = std::abs((total_energy(s) - E0) / E0);
            out << name << ',' << integ->order() << ',' << dt << ',' << steps << ','
                << steps * integ->force_evals_per_step() << ',' << pos_err << ',' << ps_err << ','
                << e_err << '\n';
        }
    }
    std::cout << "wrote " << path << '\n';
    return 0;
}

int longterm(const std::string& path) {
    std::ofstream out(path);
    out << "integrator,periods,max_rel_energy_error\n" << std::setprecision(10);
    const double T = kepler_period();
    const int periods = 1000;
    // Equal cost per period for all methods: 600 force evaluations per orbit.
    const int evals_per_period = 600;
    for (const auto& name : integrator_names()) {
        auto integ = make_integrator(name);
        const long steps_per_period = evals_per_period / integ->force_evals_per_step();
        const double dt = T / static_cast<double>(steps_per_period);
        System s = kepler(kEcc, 1.0, kMassLight);
        const double E0 = total_energy(s);
        // Record the largest error seen within each orbit (the envelope). Sampling
        // only at orbit boundaries is misleading: the energy error oscillates within
        // an orbit, and as the orbit precesses the boundary sample slides along that
        // oscillation, which makes a bounded error look like slow growth.
        for (int p = 1; p <= periods; ++p) {
            double worst = 0.0;
            for (long k = 0; k < steps_per_period; ++k) {
                integ->step(s, dt);
                worst = std::max(worst, std::abs((total_energy(s) - E0) / E0));
            }
            // Log-spaced samples keep the file small while covering the whole run.
            const bool sample = p <= 10 || (p <= 100 && p % 5 == 0) || p % 25 == 0;
            if (sample) out << name << ',' << p << ',' << std::max(worst, 1e-17) << '\n';
        }
    }
    std::cout << "wrote " << path << '\n';
    return 0;
}

int scaling(const std::string& path) {
    std::ofstream out(path);
    out << "n,seconds_per_eval,interactions_per_second\n" << std::setprecision(8);
    std::vector<Vec3> acc;
    for (std::size_t n : {64, 128, 256, 512, 1024, 2048, 4096, 8192}) {
        System s = plummer(n);
        compute_accelerations(s, acc);  // warm-up
        // Repeat until at least ~0.3 s has elapsed so small N is measured accurately.
        long reps = 0;
        const auto start = std::chrono::steady_clock::now();
        double elapsed = 0.0;
        while (elapsed < 0.3) {
            compute_accelerations(s, acc);
            ++reps;
            elapsed =
                std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
        }
        const double per_eval = elapsed / static_cast<double>(reps);
        const double pairs = 0.5 * static_cast<double>(n) * static_cast<double>(n - 1);
        out << n << ',' << per_eval << ',' << pairs / per_eval << '\n';
        std::cout << "  n=" << std::setw(5) << n << "  " << std::scientific << std::setprecision(3)
                  << per_eval << " s/eval\n"
                  << std::defaultfloat;
    }
    std::cout << "wrote " << path << '\n';
    return 0;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc != 3) {
        std::cerr << "usage: nbody_study {convergence|longterm|scaling} OUT.csv\n";
        return 2;
    }
    const std::string mode = argv[1];
    if (mode == "convergence") return convergence(argv[2]);
    if (mode == "longterm") return longterm(argv[2]);
    if (mode == "scaling") return scaling(argv[2]);
    std::cerr << "unknown mode: " << mode << '\n';
    return 2;
}
