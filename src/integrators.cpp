#include "nbody/integrators.hpp"

#include <cmath>
#include <stdexcept>

#include "nbody/forces.hpp"

namespace nbody {

namespace {

void drift(System& s, double h) {
    for (std::size_t i = 0; i < s.size(); ++i) s.pos[i] += h * s.vel[i];
}

void kick(System& s, const std::vector<Vec3>& acc, double h) {
    for (std::size_t i = 0; i < s.size(); ++i) s.vel[i] += h * acc[i];
}

}  // namespace

// ---------------------------------------------------------------- Euler ----

void ExplicitEuler::step(System& s, double dt) {
    compute_accelerations(s, acc_);
    drift(s, dt);         // uses the old velocity
    kick(s, acc_, dt);    // uses the old acceleration
    s.time += dt;
}

// ------------------------------------------------------ Symplectic Euler ----

void SymplecticEuler::step(System& s, double dt) {
    compute_accelerations(s, acc_);
    kick(s, acc_, dt);    // v_{n+1} = v_n + dt * a(r_n)
    drift(s, dt);         // r_{n+1} = r_n + dt * v_{n+1}
    s.time += dt;
}

// -------------------------------------------------------- Velocity Verlet ----

void VelocityVerlet::step(System& s, double dt) {
    // The acceleration at the end of one step is the acceleration at the
    // start of the next, so reuse it unless the system changed underneath us.
    const bool cache_valid = cached_for_ == &s && cached_time_ == s.time && acc_.size() == s.size();
    if (!cache_valid) compute_accelerations(s, acc_);

    kick(s, acc_, 0.5 * dt);
    drift(s, dt);
    compute_accelerations(s, acc_);
    kick(s, acc_, 0.5 * dt);

    s.time += dt;
    cached_for_ = &s;
    cached_time_ = s.time;
}

// ------------------------------------------------------------------- RK4 ----

void RK4::step(System& s, double dt) {
    const std::size_t n = s.size();
    r0_ = s.pos;
    v0_ = s.vel;

    auto stage = [&](std::vector<Vec3>& kr, std::vector<Vec3>& kv) {
        kr = s.vel;                     // dr/dt = v
        compute_accelerations(s, kv);   // dv/dt = a(r)
    };
    auto set_state = [&](const std::vector<Vec3>& kr, const std::vector<Vec3>& kv, double h) {
        for (std::size_t i = 0; i < n; ++i) {
            s.pos[i] = r0_[i] + h * kr[i];
            s.vel[i] = v0_[i] + h * kv[i];
        }
    };

    stage(k1r_, k1v_);
    set_state(k1r_, k1v_, 0.5 * dt);
    stage(k2r_, k2v_);
    set_state(k2r_, k2v_, 0.5 * dt);
    stage(k3r_, k3v_);
    set_state(k3r_, k3v_, dt);
    stage(k4r_, k4v_);

    const double w = dt / 6.0;
    for (std::size_t i = 0; i < n; ++i) {
        s.pos[i] = r0_[i] + w * (k1r_[i] + 2.0 * k2r_[i] + 2.0 * k3r_[i] + k4r_[i]);
        s.vel[i] = v0_[i] + w * (k1v_[i] + 2.0 * k2v_[i] + 2.0 * k3v_[i] + k4v_[i]);
    }
    s.time += dt;
}

// -------------------------------------------------------------- Yoshida4 ----

void Yoshida4::step(System& s, double dt) {
    // Triple-jump composition of leapfrog (Yoshida 1990):
    //   w1 = 1 / (2 - 2^(1/3)),  w0 = -2^(1/3) / (2 - 2^(1/3))
    // Written as drift-kick-drift so the three substeps share drifts:
    //   D(c1) K(d1) D(c2) K(d2) D(c3) K(d3) D(c4)
    static const double cbrt2 = std::cbrt(2.0);
    static const double w1 = 1.0 / (2.0 - cbrt2);
    static const double w0 = -cbrt2 / (2.0 - cbrt2);
    static const double c1 = 0.5 * w1, c2 = 0.5 * (w0 + w1);
    static const double d1 = w1, d2 = w0;

    drift(s, c1 * dt);
    compute_accelerations(s, acc_);
    kick(s, acc_, d1 * dt);
    drift(s, c2 * dt);
    compute_accelerations(s, acc_);
    kick(s, acc_, d2 * dt);
    drift(s, c2 * dt);
    compute_accelerations(s, acc_);
    kick(s, acc_, d1 * dt);
    drift(s, c1 * dt);
    s.time += dt;
}

// --------------------------------------------------------------- Factory ----

std::unique_ptr<Integrator> make_integrator(const std::string& name) {
    if (name == "euler") return std::make_unique<ExplicitEuler>();
    if (name == "symplectic-euler") return std::make_unique<SymplecticEuler>();
    if (name == "verlet") return std::make_unique<VelocityVerlet>();
    if (name == "rk4") return std::make_unique<RK4>();
    if (name == "yoshida4") return std::make_unique<Yoshida4>();
    throw std::invalid_argument("unknown integrator: " + name);
}

std::vector<std::string> integrator_names() {
    return {"euler", "symplectic-euler", "verlet", "rk4", "yoshida4"};
}

}  // namespace nbody
