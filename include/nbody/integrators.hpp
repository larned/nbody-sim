#pragma once

#include <memory>
#include <string>
#include <vector>

#include "nbody/system.hpp"

namespace nbody {

// Interface for a fixed-step time integrator.
//
// Implementations own their scratch buffers so that step() performs no heap
// allocation after the first call.
class Integrator {
public:
    virtual ~Integrator() = default;

    // Advance the system by one step of size dt.
    virtual void step(System& s, double dt) = 0;

    virtual std::string name() const = 0;

    // Global order of accuracy p: error ~ C * dt^p.
    virtual int order() const = 0;

    // Symplectic methods conserve a modified Hamiltonian, so energy error
    // stays bounded instead of drifting.
    virtual bool symplectic() const = 0;

    // Force evaluations per step; used to compare methods at equal cost.
    virtual int force_evals_per_step() const = 0;
};

// Explicit (forward) Euler. First order, not symplectic. Included as the
// baseline every other method is compared against.
class ExplicitEuler final : public Integrator {
public:
    void step(System& s, double dt) override;
    std::string name() const override { return "euler"; }
    int order() const override { return 1; }
    bool symplectic() const override { return false; }
    int force_evals_per_step() const override { return 1; }

private:
    std::vector<Vec3> acc_;
};

// Semi-implicit (symplectic) Euler: kick then drift. First order, but
// symplectic, so its energy error does not grow secularly.
class SymplecticEuler final : public Integrator {
public:
    void step(System& s, double dt) override;
    std::string name() const override { return "symplectic-euler"; }
    int order() const override { return 1; }
    bool symplectic() const override { return true; }
    int force_evals_per_step() const override { return 1; }

private:
    std::vector<Vec3> acc_;
};

// Velocity Verlet (kick-drift-kick leapfrog). Second order, symplectic,
// time-reversible. The workhorse of molecular dynamics and N-body codes.
// Caches the acceleration at the end of each step so it costs one force
// evaluation per step.
class VelocityVerlet final : public Integrator {
public:
    void step(System& s, double dt) override;
    std::string name() const override { return "verlet"; }
    int order() const override { return 2; }
    bool symplectic() const override { return true; }
    int force_evals_per_step() const override { return 1; }

private:
    std::vector<Vec3> acc_;
    const System* cached_for_ = nullptr;
    double cached_time_ = 0.0;
};

// Classical fourth-order Runge-Kutta on the first-order system
// (r, v)' = (v, a(r)). Very accurate per step but not symplectic, so energy
// drifts over long integrations.
class RK4 final : public Integrator {
public:
    void step(System& s, double dt) override;
    std::string name() const override { return "rk4"; }
    int order() const override { return 4; }
    bool symplectic() const override { return false; }
    int force_evals_per_step() const override { return 4; }

private:
    std::vector<Vec3> r0_, v0_, k1v_, k2v_, k3v_, k4v_, k1r_, k2r_, k3r_, k4r_;
};

// Yoshida (1990) fourth-order symplectic integrator: three leapfrog
// substeps with weights w1, w0, w1 chosen so the leading error terms cancel.
// Fourth order and symplectic at three force evaluations per step.
class Yoshida4 final : public Integrator {
public:
    void step(System& s, double dt) override;
    std::string name() const override { return "yoshida4"; }
    int order() const override { return 4; }
    bool symplectic() const override { return true; }
    int force_evals_per_step() const override { return 3; }

private:
    std::vector<Vec3> acc_;
};

// Factory: "euler", "symplectic-euler", "verlet", "rk4", "yoshida4".
// Throws std::invalid_argument for an unknown name.
std::unique_ptr<Integrator> make_integrator(const std::string& name);

std::vector<std::string> integrator_names();

}  // namespace nbody
