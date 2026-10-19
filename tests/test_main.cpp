// Unit and physics tests. Uses a tiny self-contained harness so the project
// builds with no third-party dependencies.

#include <cmath>
#include <functional>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "nbody/forces.hpp"
#include "nbody/integrators.hpp"
#include "nbody/scenarios.hpp"
#include "nbody/system.hpp"

using namespace nbody;

// ------------------------------------------------------------- harness ----

namespace {

struct TestCase {
    std::string name;
    std::function<void()> fn;
};

std::vector<TestCase>& registry() {
    static std::vector<TestCase> tests;
    return tests;
}

struct Register {
    Register(std::string name, std::function<void()> fn) {
        registry().push_back({std::move(name), std::move(fn)});
    }
};

struct Failure : std::runtime_error {
    using std::runtime_error::runtime_error;
};

#define TEST(name)                                     \
    void name();                                       \
    static const Register register_##name(#name, name); \
    void name()

#define CHECK(cond)                                                                       \
    do {                                                                                  \
        if (!(cond)) {                                                                    \
            std::ostringstream os_;                                                       \
            os_ << __FILE__ << ':' << __LINE__ << ": CHECK(" #cond ") failed";            \
            throw Failure(os_.str());                                                     \
        }                                                                                 \
    } while (0)

#define CHECK_LT(a, b)                                                                    \
    do {                                                                                  \
        const double a_ = (a), b_ = (b);                                                  \
        if (!(a_ < b_)) {                                                                 \
            std::ostringstream os_;                                                       \
            os_ << __FILE__ << ':' << __LINE__ << ": expected " #a " < " #b ", got " << a_ \
                << " vs " << b_;                                                          \
            throw Failure(os_.str());                                                     \
        }                                                                                 \
    } while (0)

constexpr double kPi = 3.14159265358979323846;
constexpr double kEcc = 0.5;
constexpr double kMassLight = 1e-3;
const double kPeriod = 2.0 * kPi / std::sqrt(1.0 + kMassLight);

double rel_energy_error(const System& s, double E0) { return std::abs((total_energy(s) - E0) / E0); }

// Fraction of an orbit to integrate for convergence studies. Deliberately not
// a whole period: symplectic Euler equals leapfrog plus a half-kick at each
// end, and after exactly one period those half-kicks cancel, which would make
// a first-order method look second order.
constexpr double kStudyFraction = 0.6;

// Phase-space error |dr| + |dv| of the relative orbit at t = kStudyFraction * T.
double kepler_error(const std::string& integrator, long steps) {
    System s = kepler(kEcc, 1.0, kMassLight);
    auto integ = make_integrator(integrator);
    const double dt = kStudyFraction * kPeriod / static_cast<double>(steps);
    for (long k = 0; k < steps; ++k) integ->step(s, dt);
    const double m_tot = 1.0 + kMassLight;
    const Vec3 r = kepler_exact_relative_position(kEcc, 1.0, m_tot, s.time);
    const Vec3 v = kepler_exact_relative_velocity(kEcc, 1.0, m_tot, s.time);
    return norm((s.pos[1] - s.pos[0]) - r) + norm((s.vel[1] - s.vel[0]) - v);
}

// Least-squares slope of log(error) against log(dt).
double observed_order(const std::string& integrator, int p_min, int p_max) {
    double sx = 0, sy = 0, sxx = 0, sxy = 0;
    int n = 0;
    for (int p = p_min; p <= p_max; ++p) {
        const long steps = 1L << p;
        const double x = std::log(kStudyFraction * kPeriod / static_cast<double>(steps));
        const double y = std::log(kepler_error(integrator, steps));
        sx += x; sy += y; sxx += x * x; sxy += x * y; ++n;
    }
    return (n * sxy - sx * sy) / (n * sxx - sx * sx);
}

}  // namespace

// --------------------------------------------------------------- forces ----

TEST(two_body_force_matches_newton) {
    System s;
    s.G = 2.5;
    s.add_body(3.0, {0, 0, 0}, {});
    s.add_body(5.0, {2, 0, 0}, {});
    std::vector<Vec3> a;
    compute_accelerations(s, a);
    // |a_0| = G m_1 / r^2 toward body 1; |a_1| = G m_0 / r^2 toward body 0.
    CHECK_LT(std::abs(a[0].x - 2.5 * 5.0 / 4.0), 1e-14);
    CHECK_LT(std::abs(a[1].x + 2.5 * 3.0 / 4.0), 1e-14);
    CHECK_LT(std::abs(a[0].y) + std::abs(a[0].z), 1e-15);
}

TEST(forces_obey_newtons_third_law) {
    System s = plummer(200);
    std::vector<Vec3> a;
    compute_accelerations(s, a);
    Vec3 net;
    double scale = 0.0;
    for (std::size_t i = 0; i < s.size(); ++i) {
        net += s.mass[i] * a[i];
        scale += s.mass[i] * norm(a[i]);
    }
    CHECK_LT(norm(net) / scale, 1e-13);
}

TEST(softening_limits_close_encounter_force) {
    System s;
    s.softening = 0.1;
    s.add_body(1.0, {0, 0, 0}, {});
    s.add_body(1.0, {1e-8, 0, 0}, {});
    std::vector<Vec3> a;
    compute_accelerations(s, a);
    // Softened force peaks at r = eps/sqrt(2) with value 2/(3*sqrt(3)) / eps^2.
    const double bound = 2.0 / (3.0 * std::sqrt(3.0)) / (0.1 * 0.1);
    CHECK(std::isfinite(a[0].x));
    CHECK_LT(norm(a[0]), bound + 1e-12);
}

// -------------------------------------------------------- conservation ----

TEST(all_integrators_conserve_linear_momentum) {
    for (const auto& name : integrator_names()) {
        System s = plummer(64, 7);
        auto integ = make_integrator(name);
        for (int k = 0; k < 200; ++k) integ->step(s, 1e-3);
        // Started in the CoM frame, so total momentum should stay at round-off.
        CHECK_LT(norm(total_momentum(s)), 1e-13);
    }
}

TEST(leapfrog_methods_conserve_angular_momentum_to_roundoff) {
    // Kicks change v parallel to a central-force direction and drifts move r
    // parallel to v, so each substep preserves r x v exactly for a two-body problem.
    for (const std::string name : {"symplectic-euler", "verlet", "yoshida4"}) {
        System s = kepler(kEcc, 1.0, kMassLight);
        const double L0 = norm(total_angular_momentum(s));
        auto integ = make_integrator(name);
        for (int k = 0; k < 20000; ++k) integ->step(s, kPeriod / 500.0);
        CHECK_LT(std::abs(norm(total_angular_momentum(s)) - L0) / L0, 1e-11);
    }
}

TEST(symplectic_energy_error_is_bounded) {
    // Over 200 orbits the Verlet energy error oscillates but does not grow.
    System s = kepler(kEcc, 1.0, kMassLight);
    VelocityVerlet integ;
    const double E0 = total_energy(s);
    const double dt = kPeriod / 400.0;
    double max_first = 0.0, max_last = 0.0;
    for (int p = 0; p < 200; ++p) {
        for (int k = 0; k < 400; ++k) {
            integ.step(s, dt);
            const double err = rel_energy_error(s, E0);
            if (p < 20) max_first = std::max(max_first, err);
            if (p >= 180) max_last = std::max(max_last, err);
        }
    }
    CHECK_LT(max_last, 1.5 * max_first);
}

TEST(non_symplectic_energy_error_grows) {
    // RK4 at a coarse step drifts secularly: the error after 200 orbits is far
    // larger than after 10. This is the failure mode symplectic methods avoid.
    System s = kepler(kEcc, 1.0, kMassLight);
    RK4 integ;
    const double E0 = total_energy(s);
    const double dt = kPeriod / 100.0;
    double err_10 = 0.0;
    for (int p = 1; p <= 200; ++p) {
        for (int k = 0; k < 100; ++k) integ.step(s, dt);
        if (p == 10) err_10 = rel_energy_error(s, E0);
    }
    CHECK_LT(10.0 * err_10, rel_energy_error(s, E0));
}

// ---------------------------------------------------- order of accuracy ----

TEST(observed_convergence_orders_match_theory) {
    struct Case {
        const char* name;
        double expected;
    };
    for (const Case c : {Case{"euler", 1}, Case{"symplectic-euler", 1}, Case{"verlet", 2},
                         Case{"rk4", 4}, Case{"yoshida4", 4}}) {
        // Euler needs many steps to leave the pre-asymptotic regime on an eccentric orbit.
        const bool first_order = c.expected == 1;
        const double p = first_order ? observed_order(c.name, 12, 15) : observed_order(c.name, 7, 10);
        if (std::abs(p - c.expected) > 0.25) {
            std::ostringstream os;
            os << c.name << ": observed order " << p << ", expected " << c.expected;
            throw Failure(os.str());
        }
    }
}

// ------------------------------------------------------------ symmetry ----

TEST(verlet_is_time_reversible) {
    System s = figure_eight();
    const System start = s;
    VelocityVerlet fwd;
    for (int k = 0; k < 1000; ++k) fwd.step(s, 1e-3);
    for (auto& v : s.vel) v = -v;
    VelocityVerlet back;
    for (int k = 0; k < 1000; ++k) back.step(s, 1e-3);
    double err = 0.0;
    for (std::size_t i = 0; i < s.size(); ++i) err += norm(s.pos[i] - start.pos[i]);
    CHECK_LT(err, 1e-11);
}

TEST(figure_eight_returns_after_one_period) {
    System s = figure_eight();
    const System start = s;
    Yoshida4 integ;
    const double T = 6.32591398;
    const long steps = 20000;
    for (long k = 0; k < steps; ++k) integ.step(s, T / steps);
    double err = 0.0;
    for (std::size_t i = 0; i < s.size(); ++i) err += norm(s.pos[i] - start.pos[i]);
    // Limited by the 8-digit published initial conditions, not by the integrator.
    CHECK_LT(err, 1e-5);
}

// ------------------------------------------------------------- utility ----

TEST(kepler_exact_solution_is_periodic) {
    const Vec3 r0 = kepler_exact_relative_position(0.7, 2.0, 1.0, 0.0);
    const double T = 2.0 * kPi * std::sqrt(8.0);
    const Vec3 r1 = kepler_exact_relative_position(0.7, 2.0, 1.0, T);
    CHECK_LT(norm(r0 - Vec3{0.6, 0, 0}), 1e-14);
    CHECK_LT(norm(r1 - r0), 1e-12);

    // Velocity should equal the derivative of position (central difference).
    const double h = 1e-5, t = 1.234;
    const Vec3 fd = (kepler_exact_relative_position(0.7, 2.0, 1.0, t + h) -
                     kepler_exact_relative_position(0.7, 2.0, 1.0, t - h)) * (0.5 / h);
    CHECK_LT(norm(fd - kepler_exact_relative_velocity(0.7, 2.0, 1.0, t)), 1e-8);
}

TEST(factory_rejects_unknown_names) {
    bool threw = false;
    try {
        make_integrator("leapfrog-but-wrong");
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    CHECK(threw);
    for (const auto& name : integrator_names()) CHECK(make_integrator(name)->name() == name);
}

TEST(solar_system_earth_returns_after_one_year) {
    System s = solar_system();
    std::size_t earth = 0;
    while (s.name[earth] != "Earth") ++earth;
    const Vec3 start = s.pos[earth] - s.pos[0];
    Yoshida4 integ;
    for (int k = 0; k < 3650; ++k) integ.step(s, 1.0 / 3650.0);
    // Planetary perturbations and the Sun's reflex motion keep this from
    // being exact; a tolerance of 1e-3 AU (~150,000 km) is still a strong check.
    CHECK_LT(norm((s.pos[earth] - s.pos[0]) - start), 1e-3);
}

// ---------------------------------------------------------------- main ----

int main() {
    int failed = 0;
    for (const auto& t : registry()) {
        try {
            t.fn();
            std::cout << "[ PASS ] " << t.name << '\n';
        } catch (const std::exception& e) {
            ++failed;
            std::cout << "[ FAIL ] " << t.name << "\n         " << e.what() << '\n';
        }
    }
    std::cout << '\n' << registry().size() - failed << '/' << registry().size() << " tests passed\n";
    return failed == 0 ? 0 : 1;
}
