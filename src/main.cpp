// nbody: command-line driver.
//
// Runs one scenario with one integrator and writes two CSV files:
//   trajectory  t, body, x, y, z           (every --output-every steps)
//   diagnostics t, energy, rel_energy_error, momentum, angular_momentum

#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>

#include "nbody/integrators.hpp"
#include "nbody/scenarios.hpp"
#include "nbody/system.hpp"

namespace {

struct Options {
    std::string scenario = "kepler";
    std::string integrator = "verlet";
    double dt = 0.01;
    double t_end = 10.0;
    long output_every = 10;
    std::size_t n = 256;
    double softening = -1.0;  // < 0 means "use the scenario default"
    std::string trajectory_path = "trajectory.csv";
    std::string diagnostics_path = "diagnostics.csv";
    bool quiet = false;
};

void print_usage() {
    std::cout <<
        R"(Usage: nbody [options]

  --scenario NAME      kepler | figure8 | solar | plummer      (default: kepler)
  --integrator NAME    euler | symplectic-euler | verlet | rk4 | yoshida4
                                                                (default: verlet)
  --dt X               time step                               (default: 0.01)
  --t-end X            total simulated time                    (default: 10)
  --output-every K     write trajectory/diagnostics every K steps (default: 10)
  --n N                number of bodies for plummer            (default: 256)
  --softening X        override Plummer softening length
  --trajectory PATH    trajectory CSV path ("" to skip)        (default: trajectory.csv)
  --diagnostics PATH   diagnostics CSV path ("" to skip)       (default: diagnostics.csv)
  --quiet              suppress the summary
  --help               show this message

Units: kepler, figure8 and plummer use G = 1; solar uses AU, years, solar masses.
)";
}

Options parse(int argc, char** argv) {
    Options o;
    std::map<std::string, std::string*> strings{{"--scenario", &o.scenario},
                                                {"--integrator", &o.integrator},
                                                {"--trajectory", &o.trajectory_path},
                                                {"--diagnostics", &o.diagnostics_path}};
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--help" || arg == "-h") {
            print_usage();
            std::exit(0);
        }
        if (arg == "--quiet") {
            o.quiet = true;
            continue;
        }
        if (i + 1 >= argc) throw std::invalid_argument("missing value for " + arg);
        const std::string val = argv[++i];
        if (auto it = strings.find(arg); it != strings.end()) *it->second = val;
        else if (arg == "--dt") o.dt = std::stod(val);
        else if (arg == "--t-end") o.t_end = std::stod(val);
        else if (arg == "--output-every") o.output_every = std::stol(val);
        else if (arg == "--n") o.n = std::stoul(val);
        else if (arg == "--softening") o.softening = std::stod(val);
        else throw std::invalid_argument("unknown option: " + arg);
    }
    if (o.dt <= 0.0) throw std::invalid_argument("--dt must be positive");
    if (o.t_end <= 0.0) throw std::invalid_argument("--t-end must be positive");
    if (o.output_every < 1) throw std::invalid_argument("--output-every must be >= 1");
    return o;
}

}  // namespace

int main(int argc, char** argv) {
    using namespace nbody;
    Options opt;
    try {
        opt = parse(argc, argv);
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << "\n\n";
        print_usage();
        return 2;
    }

    try {
        System s = make_scenario(opt.scenario, opt.n);
        if (opt.softening >= 0.0) s.softening = opt.softening;
        auto integ = make_integrator(opt.integrator);

        std::ofstream traj, diag;
        if (!opt.trajectory_path.empty()) {
            traj.open(opt.trajectory_path);
            if (!traj) throw std::runtime_error("cannot open " + opt.trajectory_path);
            traj << "t,body,x,y,z\n" << std::setprecision(10);
        }
        if (!opt.diagnostics_path.empty()) {
            diag.open(opt.diagnostics_path);
            if (!diag) throw std::runtime_error("cannot open " + opt.diagnostics_path);
            diag << "t,energy,rel_energy_error,momentum,angular_momentum\n" << std::setprecision(12);
        }

        const double E0 = total_energy(s);
        double max_rel_err = 0.0;

        auto record = [&]() {
            const double E = total_energy(s);
            const double rel = std::abs((E - E0) / E0);
            if (rel > max_rel_err) max_rel_err = rel;
            if (traj.is_open()) {
                for (std::size_t i = 0; i < s.size(); ++i) {
                    traj << s.time << ',' << i << ',' << s.pos[i].x << ',' << s.pos[i].y << ','
                         << s.pos[i].z << '\n';
                }
            }
            if (diag.is_open()) {
                diag << s.time << ',' << E << ',' << rel << ',' << norm(total_momentum(s)) << ','
                     << norm(total_angular_momentum(s)) << '\n';
            }
        };

        const long steps = std::lround(opt.t_end / opt.dt);
        record();
        const auto start = std::chrono::steady_clock::now();
        for (long k = 1; k <= steps; ++k) {
            integ->step(s, opt.dt);
            if (k % opt.output_every == 0 || k == steps) record();
        }
        const double secs =
            std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();

        if (!opt.quiet) {
            std::cout << "scenario        " << opt.scenario << " (" << s.size() << " bodies)\n"
                      << "integrator      " << integ->name() << " (order " << integ->order()
                      << (integ->symplectic() ? ", symplectic" : "") << ")\n"
                      << "steps           " << steps << " x dt=" << opt.dt << "\n"
                      << std::scientific << std::setprecision(3)
                      << "max |dE/E0|     " << max_rel_err << "\n"
                      << std::fixed << std::setprecision(3)
                      << "wall time       " << secs << " s\n";
        }
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << '\n';
        return 1;
    }
    return 0;
}
