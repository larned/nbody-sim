#include "nbody/scenarios.hpp"

#include <cmath>
#include <random>
#include <stdexcept>

namespace nbody {

namespace {

constexpr double kPi = 3.14159265358979323846;

// Uniform double in (0, 1) built from raw mt19937 output. std::uniform_real_distribution
// is implementation-defined, so using it would make benchmarks differ across compilers.
double uniform01(std::mt19937& gen) { return (static_cast<double>(gen()) + 0.5) / 4294967296.0; }

}  // namespace

System kepler(double e, double a, double m_light) {
    if (e < 0.0 || e >= 1.0) throw std::invalid_argument("kepler: eccentricity must be in [0, 1)");
    System s;
    s.G = 1.0;
    const double M = 1.0;
    const double m_tot = M + m_light;

    // Relative orbit (light minus heavy) starting at periapsis on +x, moving +y.
    const double r_p = a * (1.0 - e);
    const double v_p = std::sqrt(s.G * m_tot * (1.0 + e) / (a * (1.0 - e)));
    const Vec3 r_rel{r_p, 0.0, 0.0};
    const Vec3 v_rel{0.0, v_p, 0.0};

    // Split the relative coordinates by mass so the centre of mass is fixed at the origin.
    s.add_body(M, -(m_light / m_tot) * r_rel, -(m_light / m_tot) * v_rel, "primary");
    s.add_body(m_light, (M / m_tot) * r_rel, (M / m_tot) * v_rel, "secondary");
    s.to_com_frame();
    return s;
}

namespace {

double mean_motion(double a, double m_tot) { return std::sqrt(m_tot / (a * a * a)); }  // G = 1

// Solve Kepler's equation E - e sin E = M for the eccentric anomaly with Newton's method.
double eccentric_anomaly(double e, double a, double m_tot, double t) {
    const double mean_anomaly = std::fmod(mean_motion(a, m_tot) * t, 2.0 * kPi);
    double E = e < 0.8 ? mean_anomaly : kPi;
    for (int it = 0; it < 50; ++it) {
        const double dE = (E - e * std::sin(E) - mean_anomaly) / (1.0 - e * std::cos(E));
        E -= dE;
        if (std::abs(dE) < 1e-15) break;
    }
    return E;
}

}  // namespace

Vec3 kepler_exact_relative_position(double e, double a, double m_tot, double t) {
    const double E = eccentric_anomaly(e, a, m_tot, t);
    return {a * (std::cos(E) - e), a * std::sqrt(1.0 - e * e) * std::sin(E), 0.0};
}

Vec3 kepler_exact_relative_velocity(double e, double a, double m_tot, double t) {
    const double E = eccentric_anomaly(e, a, m_tot, t);
    const double Edot = mean_motion(a, m_tot) / (1.0 - e * std::cos(E));
    return {-a * std::sin(E) * Edot, a * std::sqrt(1.0 - e * e) * std::cos(E) * Edot, 0.0};
}

System figure_eight() {
    // Initial conditions from Chenciner & Montgomery (2000), as computed by Simo.
    System s;
    s.G = 1.0;
    const Vec3 r1{0.97000436, -0.24308753, 0.0};
    const Vec3 v3{-0.93240737, -0.86473146, 0.0};
    s.add_body(1.0, r1, -0.5 * v3, "A");
    s.add_body(1.0, -r1, -0.5 * v3, "B");
    s.add_body(1.0, Vec3{}, v3, "C");
    s.to_com_frame();
    return s;
}

System solar_system() {
    System s;
    s.G = 4.0 * kPi * kPi;  // AU^3 / (M_sun * yr^2)

    struct Planet {
        const char* name;
        double mass;  // solar masses
        double a;     // AU
    };
    static const Planet planets[] = {
        {"Mercury", 1.660e-7, 0.387}, {"Venus", 2.447e-6, 0.723},  {"Earth", 3.003e-6, 1.000},
        {"Mars", 3.227e-7, 1.524},    {"Jupiter", 9.545e-4, 5.203}, {"Saturn", 2.858e-4, 9.537},
        {"Uranus", 4.366e-5, 19.19},  {"Neptune", 5.151e-5, 30.07},
    };

    s.add_body(1.0, Vec3{}, Vec3{}, "Sun");
    // Spread starting longitudes by the golden angle so no two planets start aligned.
    const double golden = kPi * (3.0 - std::sqrt(5.0));
    int k = 0;
    for (const auto& p : planets) {
        const double theta = golden * k++;
        const double v = std::sqrt(s.G * (1.0 + p.mass) / p.a);
        s.add_body(p.mass, Vec3{p.a * std::cos(theta), p.a * std::sin(theta), 0.0},
                   Vec3{-v * std::sin(theta), v * std::cos(theta), 0.0}, p.name);
    }
    s.to_com_frame();
    return s;
}

System plummer(std::size_t n, unsigned seed, double softening) {
    // Sampling procedure from Aarseth, Henon & Wielen (1974).
    System s;
    s.G = 1.0;
    s.softening = softening;
    std::mt19937 gen(seed);
    const double m = 1.0 / static_cast<double>(n);

    auto random_direction = [&]() {
        const double z = 2.0 * uniform01(gen) - 1.0;
        const double phi = 2.0 * kPi * uniform01(gen);
        const double rho = std::sqrt(1.0 - z * z);
        return Vec3{rho * std::cos(phi), rho * std::sin(phi), z};
    };

    while (s.size() < n) {
        // Radius from the inverse cumulative mass profile; discard the far tail.
        const double X = uniform01(gen);
        const double r = 1.0 / std::sqrt(std::pow(X, -2.0 / 3.0) - 1.0);
        if (r > 10.0) continue;

        // Speed via von Neumann rejection on g(q) = q^2 (1 - q^2)^(7/2), q = v / v_esc.
        double q = 0.0;
        while (true) {
            q = uniform01(gen);
            const double g = q * q * std::pow(1.0 - q * q, 3.5);
            if (0.1 * uniform01(gen) < g) break;
        }
        const double v_esc = std::sqrt(2.0) * std::pow(1.0 + r * r, -0.25);
        s.add_body(m, r * random_direction(), q * v_esc * random_direction());
    }
    s.to_com_frame();
    return s;
}

System make_scenario(const std::string& name, std::size_t n) {
    if (name == "kepler") return kepler(0.5);
    if (name == "figure8") return figure_eight();
    if (name == "solar") return solar_system();
    if (name == "plummer") return plummer(n);
    throw std::invalid_argument("unknown scenario: " + name);
}

std::vector<std::string> scenario_names() { return {"kepler", "figure8", "solar", "plummer"}; }

}  // namespace nbody
