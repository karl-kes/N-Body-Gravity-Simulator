// Integration-step scaling benchmark: serial vs OpenMP across a range of N.
// Times full Yoshida integration steps (3 force evaluations + drift/kick
// updates). Sweeps both the direct O(N^2) kernel and Barnes-Hut O(N log N)
// so the crossover is visible in a single run.
//
// Build:  cmake --build build --target benchmark
// Run:    ./build/benchmark
//         ./build/benchmark --max-n 65536 --trials 5
//         ./build/benchmark --force bh --theta 0.3
//         ./build/benchmark --force direct

#include "../src/Particle/Particle.hpp"
#include "../src/Force/Force.hpp"
#include "../src/Force/BarnesHut.hpp"
#include "../src/Integrator/Integrator.hpp"
#include "../src/Config.hpp"

#include <iostream>
#include <iomanip>
#include <cmath>
#include <chrono>
#include <vector>
#include <memory>
#include <random>
#include <string>
#include <algorithm>
#include <sstream>

#include <omp.h>

enum class ForceKind { Direct, BarnesHut };

static char const* force_kind_label( ForceKind k ) {
    return k == ForceKind::Direct ? "direct" : "bh";
}

// Populate N bodies with random state. Enforces a minimum pairwise separation
// of `min_sep` meters to avoid near-singular force evaluations that would
// produce non-representative dynamics during the timed integration.
static void populate_random( Particles &p, std::size_t const N, double const min_sep = 1e9,
                             unsigned const seed = 42 ) {
    std::mt19937_64 rng{ seed };
    std::uniform_real_distribution<double> pos_dist{ -1e12, 1e12 };
    std::uniform_real_distribution<double> vel_dist{ -1e4, 1e4 };
    std::uniform_real_distribution<double> mass_dist{ 1e20, 1e30 };

    for ( std::size_t i{}; i < N; ++i ) {
        // Rejection-sample positions until minimum separation is satisfied.
        // For the densities used here, rejection is rare (< 1% of attempts).
        bool accepted{ false };
        while ( !accepted ) {
            double const px{ pos_dist( rng ) };
            double const py{ pos_dist( rng ) };
            double const pz{ pos_dist( rng ) };

            accepted = true;
            for ( std::size_t j{}; j < i; ++j ) {
                double const dx{ px - p.pos_x()[j] };
                double const dy{ py - p.pos_y()[j] };
                double const dz{ pz - p.pos_z()[j] };
                if ( dx*dx + dy*dy + dz*dz < min_sep * min_sep ) {
                    accepted = false;
                    break;
                }
            }

            if ( accepted ) {
                p.pos_x()[i] = px;
                p.pos_y()[i] = py;
                p.pos_z()[i] = pz;
            }
        }

        p.vel_x()[i] = vel_dist( rng );
        p.vel_y()[i] = vel_dist( rng );
        p.vel_z()[i] = vel_dist( rng );
        p.acc_x()[i] = 0.0;
        p.acc_y()[i] = 0.0;
        p.acc_z()[i] = 0.0;
        p.mass()[i]  = mass_dist( rng );
    }
}

struct BenchResult {
    std::size_t N;
    std::size_t steps;

    bool has_direct;
    double direct_serial_ms;
    double direct_omp_ms;
    double direct_serial_gflops;
    double direct_omp_gflops;
    double direct_speedup;

    bool has_bh;
    double bh_serial_ms;
    double bh_omp_ms;
    double bh_speedup;

    // Cross-method ratio at OMP setting; only meaningful when both run.
    double bh_over_direct_omp;
};

static std::unique_ptr<Force> make_force( ForceKind kind, double theta ) {
    if ( kind == ForceKind::BarnesHut ) {
        return std::make_unique<Gravity_BarnesHut>( theta );
    }
    return std::make_unique<Gravity>();
}

// Caller must set the OMP thread count before calling.
static double run_timed( std::size_t const N, std::size_t const steps,
                         ForceKind const kind, double const theta ) {
    Particles p{ N };
    populate_random( p, N );

    std::vector<std::unique_ptr<Force>> forces;
    forces.push_back( make_force( kind, theta ) );
    Yoshida integ{ 900.0 };

    // Warmup: 3 steps to stabilize thread pool and caches.
    for ( std::size_t s{}; s < 3; ++s ) {
        integ.integrate( p, forces );
    }

    auto const start{ std::chrono::high_resolution_clock::now() };

    for ( std::size_t s{}; s < steps; ++s ) {
        integ.integrate( p, forces );
    }

    auto const end{ std::chrono::high_resolution_clock::now() };

    return std::chrono::duration<double, std::milli>( end - start ).count();
}

static double run_median( std::size_t const N, std::size_t const steps,
                          int const num_threads, std::size_t const trials,
                          ForceKind const kind, double const theta ) {
    omp_set_num_threads( num_threads );

    std::vector<double> timings;
    timings.reserve( trials );

    for ( std::size_t t{}; t < trials; ++t ) {
        timings.push_back( run_timed( N, steps, kind, theta ) );
    }

    std::sort( timings.begin(), timings.end() );
    return timings[trials / 2];
}

// 3 force evaluations x N x N pairwise interactions x ~27 FLOPs per pair
// (sub, mul, add for dx/dy/dz, R_sq, 1/sqrt, mul chain, mask, accumulate)
// plus drift/kick updates: ~84N FLOPs per step. The BH kernel is data
// dependent and is not modeled here.
static double direct_flops_per_step( std::size_t const N ) {
    return 3.0 * N * N * 27.0 + 84.0 * N;
}

// Calibrate step count so the chosen kernel's serial runtime is approximately
// `target_ms`. The same step count is then reused by both kernels at this N,
// so wall times are directly comparable.
static std::size_t calibrate_steps( std::size_t const N, double const target_ms,
                                    ForceKind const kind, double const theta ) {
    omp_set_num_threads( 1 );
    double const calibration_ms{ run_timed( N, 5, kind, theta ) };
    double const ms_per_step{ calibration_ms / 5.0 };
    std::size_t const estimated{ static_cast<std::size_t>( target_ms / ms_per_step ) };
    return std::clamp( estimated, static_cast<std::size_t>( 3 ), static_cast<std::size_t>( 500 ) );
}

int main( int argc, char* argv[] ) {
    std::size_t max_n{ 65536 };
    std::size_t num_trials{ 3 };
    double target_ms{ 2000.0 };
    double theta{ 0.5 };
    std::string mode{ "both" };
    int const max_threads{ omp_get_max_threads() };
    int omp_threads{ max_threads };

    for ( int i{ 1 }; i < argc; ++i ) {
        std::string const arg{ argv[i] };
        if ( arg == "--max-n" && i + 1 < argc ) { max_n = std::stoull( argv[++i] ); }
        else if ( arg == "--trials" && i + 1 < argc ) { num_trials = std::stoull( argv[++i] ); }
        else if ( arg == "--target-ms" && i + 1 < argc ) { target_ms = std::stod( argv[++i] ); }
        else if ( arg == "--threads" && i + 1 < argc ) { omp_threads = std::stoi( argv[++i] ); }
        else if ( arg == "--theta" && i + 1 < argc ) { theta = std::stod( argv[++i] ); }
        else if ( arg == "--force" && i + 1 < argc ) { mode = argv[++i]; }
        else if ( arg == "-h" || arg == "--help" ) {
            std::cout << "Usage: benchmark [--max-n N] [--trials N] [--target-ms MS]\n"
                      << "                 [--threads N] [--force {direct|bh|both}] [--theta T]\n"
                      << "  --max-n N               Maximum N for sweep (default: 65536)\n"
                      << "  --trials N              Trials per config, reports median (default: 3)\n"
                      << "  --target-ms MS          Target serial runtime per trial in ms (default: 2000)\n"
                      << "  --threads N             OMP thread count for parallel runs (default: max)\n"
                      << "  --force MODE            'direct', 'bh', or 'both' (default: both)\n"
                      << "  --theta T               BH opening angle (default: 0.5)\n";
            return 0;
        }
    }

    if ( mode != "direct" && mode != "bh" && mode != "both" ) {
        std::cerr << "error: --force must be 'direct', 'bh', or 'both'\n";
        return 1;
    }

    std::cout << "\n<--- N-Body Scaling Benchmark --->\n"
              << "  Integrator:        Yoshida 4th-order (dt = 900 s)\n"
              << "  Force mode:        " << mode << "\n"
              << "  Theta (BH):        " << theta << "\n"
              << "  OMP threads:       " << omp_threads << "\n"
              << "  Trials/config:     " << num_trials << " (median)\n"
              << "  Target serial:     " << std::fixed << std::setprecision( 0 ) << target_ms << " ms per trial\n"
              << "  OMP threshold:     N >= " << config::OMP_THRESHOLD << "\n";
    std::cout << "\n";

    // Start at OMP_THRESHOLD (rounded up to next power of 2) so the parallel
    // code path is always active.
    std::vector<std::size_t> N_values;
    std::size_t n_start{ 1 };
    while ( n_start < config::OMP_THRESHOLD ) n_start *= 2;
    for ( std::size_t n{ n_start }; n <= max_n; n *= 2 ) {
        N_values.push_back( n );
    }

    std::vector<BenchResult> results;

    std::cout << std::left
              << std::setw( 8 )  << "N"
              << std::setw( 8 )  << "Steps"
              << std::setw( 14 ) << "Direct(ms)"
              << std::setw( 14 ) << "Dir-OMP(ms)"
              << std::setw( 11 ) << "Dir-Spd"
              << std::setw( 14 ) << "BH(ms)"
              << std::setw( 14 ) << "BH-OMP(ms)"
              << std::setw( 11 ) << "BH-Spd"
              << std::setw( 13 ) << "Winner(OMP)"
              << "\n";
    std::cout << std::string( 107, '=' ) << "\n";

    for ( std::size_t const N : N_values ) {
        bool const run_direct{ mode == "direct" || mode == "both" };
        bool const run_bh{ mode == "bh" || mode == "both" };

        ForceKind const calib_kind{ run_bh ? ForceKind::BarnesHut : ForceKind::Direct };
        std::size_t const steps{ calibrate_steps( N, target_ms, calib_kind, theta ) };

        BenchResult r{};
        r.N = N;
        r.steps = steps;
        r.has_direct = run_direct;
        r.has_bh = run_bh;

        if ( run_direct ) {
            double const ser{ run_median( N, steps, 1, num_trials, ForceKind::Direct, theta ) };
            double const omp{ run_median( N, steps, omp_threads, num_trials, ForceKind::Direct, theta ) };
            double const tot_flops{ direct_flops_per_step( N ) * steps };
            r.direct_serial_ms = ser;
            r.direct_omp_ms = omp;
            r.direct_serial_gflops = tot_flops / ( ser * 1e6 );
            r.direct_omp_gflops = tot_flops / ( omp * 1e6 );
            r.direct_speedup = ser / omp;
        }

        if ( run_bh ) {
            double const ser{ run_median( N, steps, 1, num_trials, ForceKind::BarnesHut, theta ) };
            double const omp{ run_median( N, steps, omp_threads, num_trials, ForceKind::BarnesHut, theta ) };
            r.bh_serial_ms = ser;
            r.bh_omp_ms = omp;
            r.bh_speedup = ser / omp;
        }

        if ( run_direct && run_bh ) {
            // > 1: BH wins; < 1: Direct wins.
            r.bh_over_direct_omp = r.direct_omp_ms / r.bh_omp_ms;
        }

        results.push_back( r );

        auto fmt_speedup = []( double spd ) -> std::string {
            std::ostringstream oss;
            oss << std::fixed << std::setprecision( 2 ) << spd << "x";
            return oss.str();
        };

        std::cout << std::left << std::fixed
                  << std::setw( 8 ) << N
                  << std::setw( 8 ) << steps;
        if ( run_direct ) {
            std::cout << std::setw( 14 ) << std::setprecision( 1 ) << r.direct_serial_ms
                      << std::setw( 14 ) << std::setprecision( 1 ) << r.direct_omp_ms
                      << std::setw( 11 ) << fmt_speedup( r.direct_speedup );
        } else {
            std::cout << std::setw( 14 ) << "--"
                      << std::setw( 14 ) << "--"
                      << std::setw( 11 ) << "--";
        }
        if ( run_bh ) {
            std::cout << std::setw( 14 ) << std::setprecision( 1 ) << r.bh_serial_ms
                      << std::setw( 14 ) << std::setprecision( 1 ) << r.bh_omp_ms
                      << std::setw( 11 ) << fmt_speedup( r.bh_speedup );
        } else {
            std::cout << std::setw( 14 ) << "--"
                      << std::setw( 14 ) << "--"
                      << std::setw( 11 ) << "--";
        }
        if ( run_direct && run_bh ) {
            std::ostringstream winner;
            if ( r.bh_over_direct_omp >= 1.0 ) {
                winner << "BH " << std::fixed << std::setprecision( 2 ) << r.bh_over_direct_omp << "x";
            } else {
                winner << "Dir " << std::fixed << std::setprecision( 2 ) << ( 1.0 / r.bh_over_direct_omp ) << "x";
            }
            std::cout << std::setw( 13 ) << winner.str();
        } else {
            std::cout << std::setw( 13 ) << "--";
        }
        std::cout << "\n" << std::flush;
    }

    std::cout << std::string( 107, '=' ) << "\n\n";

    omp_set_num_threads( max_threads );

    // CSV output for plotting. GFLOP/s reported only for the direct kernel.
    std::cout << "CSV (for plotting):\n";
    std::cout << "N,steps,method,serial_ms,omp_ms,speedup,gflops_serial,gflops_omp\n";
    for ( auto const &r : results ) {
        if ( r.has_direct ) {
            std::cout << r.N << "," << r.steps << "," << force_kind_label( ForceKind::Direct ) << ","
                      << std::fixed << std::setprecision( 2 )
                      << r.direct_serial_ms << "," << r.direct_omp_ms << ","
                      << std::setprecision( 3 ) << r.direct_speedup << ","
                      << r.direct_serial_gflops << "," << r.direct_omp_gflops << "\n";
        }
        if ( r.has_bh ) {
            std::cout << r.N << "," << r.steps << "," << force_kind_label( ForceKind::BarnesHut ) << ","
                      << std::fixed << std::setprecision( 2 )
                      << r.bh_serial_ms << "," << r.bh_omp_ms << ","
                      << std::setprecision( 3 ) << r.bh_speedup << ","
                      << ",\n";  // empty gflops fields; BH FLOP model is data-dependent
        }
    }

    return 0;
}
