#include "Force/Force.hpp"
#include "Force/BarnesHut.hpp"
#include "Integrator/Integrator.hpp"
#include "Particle/Particle.hpp"
#include "Config.hpp"

#include <iostream>
#include <iomanip>
#include <sstream>
#include <chrono>
#include <vector>
#include <memory>
#include <random>
#include <algorithm>
#include <functional>
#include <string>
#include <string_view>
#include <cstddef>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

constexpr double      DT{ 900.0 };
constexpr double      THETA{ 0.5 };
constexpr std::size_t DEFAULT_TRIALS{ 3 };
constexpr double      TARGET_MS{ 2000.0 };

// Threshold above which single-threaded measurements are skipped.
// At large N, serial direct takes minutes per trial. OMP results are what
// matter for the portfolio, so skipping serial above this saves wall time
// without losing useful data.
constexpr std::size_t SERIAL_MAX_N{ 4096 };


void populate_random( Particles &p, unsigned const seed ) {
    std::mt19937_64 rng{ seed };
    std::uniform_real_distribution<double> pos_dist{ -1e11, 1e11 };
    std::uniform_real_distribution<double> vel_dist{ -1e3, 1e3 };
    std::uniform_real_distribution<double> mass_dist{ 1e23, 1e26 };

    std::size_t const N{ p.num_particles() };
    for ( std::size_t i{}; i < N; ++i ) {
        p.pos_x()[i] = pos_dist( rng );
        p.pos_y()[i] = pos_dist( rng );
        p.pos_z()[i] = pos_dist( rng );
        p.vel_x()[i] = vel_dist( rng );
        p.vel_y()[i] = vel_dist( rng );
        p.vel_z()[i] = vel_dist( rng );
        p.mass()[i]  = mass_dist( rng );
        p.acc_x()[i] = 0.0;
        p.acc_y()[i] = 0.0;
        p.acc_z()[i] = 0.0;
    }
}


using ForceFactory = std::function<std::unique_ptr<Force>()>;


// Run one trial: integrate `steps` steps from a fresh random IC.
double time_trial( std::size_t const N, std::size_t const steps,
                   ForceFactory const &make_force ) {
    Particles p{ N };
    populate_random( p, 1234 );

    Yoshida integ{ DT };
    std::vector<std::unique_ptr<Force>> forces;
    forces.push_back( make_force() );

    auto const start{ std::chrono::high_resolution_clock::now() };
    for ( std::size_t s{}; s < steps; ++s ) {
        integ.integrate( p, forces );
    }
    auto const end{ std::chrono::high_resolution_clock::now() };

    return std::chrono::duration<double, std::milli>( end - start ).count();
}


// Median of `trials` runs.
double median_of( std::size_t const N, std::size_t const steps,
                  std::size_t const trials,
                  ForceFactory const &make_force ) {
    std::vector<double> times;
    times.reserve( trials );
    for ( std::size_t t{}; t < trials; ++t ) {
        times.push_back( time_trial( N, steps, make_force ) );
    }
    std::sort( times.begin(), times.end() );
    return times[ trials / 2 ];
}


// Calibrate step count so a single trial hits ~TARGET_MS. Probe is direct OMP
// since it runs on every row (skipping serial above SERIAL_MAX_N doesn't skip
// calibration).
std::size_t calibrate_steps( std::size_t const N ) {
    constexpr std::size_t probe_steps{ 4 };
    double const probe_ms{ time_trial( N, probe_steps,
                                       [&] { return std::make_unique<Gravity>(); } ) };
    if ( probe_ms <= 0.0 ) return 1000;
    double const factor{ TARGET_MS / probe_ms };
    return std::max<std::size_t>( 3, static_cast<std::size_t>( probe_steps * factor ) );
}


std::string fmt( double const v ) {
    std::stringstream s;
    s << std::fixed << std::setprecision( 1 ) << v;
    return s.str();
}


std::string fmt_speedup( double const v ) {
    std::stringstream s;
    s << std::fixed << std::setprecision( 2 ) << v << "x";
    return s.str();
}

}  // namespace


int main( int argc, char* argv[] ) {
    std::size_t max_n{ 65536 };
    std::size_t trials{ DEFAULT_TRIALS };

    for ( int i{ 1 }; i < argc; ++i ) {
        std::string_view const arg{ argv[i] };
        if      ( arg == "--max-n"  && i + 1 < argc ) { max_n  = std::stoul( argv[++i] ); }
        else if ( arg == "--trials" && i + 1 < argc ) { trials = std::stoul( argv[++i] ); }
        else if ( arg == "-h" || arg == "--help" ) {
            std::cout << "Usage: benchmark [--max-n N] [--trials N]\n"
                      << "  --max-n N    Largest N to test (default: 65536)\n"
                      << "  --trials N   Trials per config (default: 3, takes median)\n"
                      << "\nSingle-threaded measurements are skipped for N > "
                      << SERIAL_MAX_N << ".\n";
            return 0;
        }
    }

    int omp_threads{ 1 };
    #ifdef _OPENMP
        omp_threads = omp_get_max_threads();
    #endif

    std::cout << "<--- N-Body Scaling Benchmark --->\n"
              << "  Integrator:        Yoshida 4th-order (dt = " << DT << " s)\n"
              << "  Force mode:        both\n"
              << "  Theta (BH):        " << THETA << "\n"
              << "  OMP threads:       " << omp_threads << "\n"
              << "  Trials/config:     " << trials << " (median)\n"
              << "  Target serial:     " << static_cast<int>( TARGET_MS ) << " ms per trial\n"
              << "  OMP threshold:     N >= " << config::OMP_THRESHOLD << "\n"
              << "  Skip serial above: N = " << SERIAL_MAX_N << "\n\n";

    std::cout << std::left
              << std::setw( 8 )  << "N"
              << std::setw( 8 )  << "Steps"
              << std::setw( 14 ) << "Direct(ms)"
              << std::setw( 14 ) << "Dir-OMP(ms)"
              << std::setw( 11 ) << "Dir-Spd"
              << std::setw( 14 ) << "BH(ms)"
              << std::setw( 14 ) << "BH-OMP(ms)"
              << std::setw( 11 ) << "BH-Spd"
              << "Winner(OMP)"
              << "\n";
    std::cout << std::string( 107, '=' ) << "\n";

    for ( std::size_t N{ 512 }; N <= max_n; N *= 2 ) {
        std::size_t const steps{ calibrate_steps( N ) };
        bool const skip_serial{ N > SERIAL_MAX_N };

        ForceFactory make_direct{ [&] { return std::make_unique<Gravity>(); } };
        ForceFactory make_bh{     [&] { return std::make_unique<Gravity_BarnesHut>( THETA ); } };

        // Serial runs (when requested).
        double dir_ms{}, bh_ms{};
        if ( !skip_serial ) {
            #ifdef _OPENMP
                omp_set_num_threads( 1 );
            #endif
            dir_ms = median_of( N, steps, trials, make_direct );
            bh_ms  = median_of( N, steps, trials, make_bh );
        }

        // OMP runs (always).
        #ifdef _OPENMP
            omp_set_num_threads( omp_threads );
        #endif
        double const dir_omp_ms{ median_of( N, steps, trials, make_direct ) };
        double const bh_omp_ms{  median_of( N, steps, trials, make_bh ) };

        // Format columns. Serial-skip rows print "skip" instead of a number.
        std::string const dir_col{ skip_serial ? "skip" : fmt( dir_ms ) };
        std::string const bh_col{  skip_serial ? "skip" : fmt( bh_ms ) };
        std::string const dir_spd_col{ skip_serial ? "--" : fmt_speedup( dir_ms / dir_omp_ms ) };
        std::string const bh_spd_col{  skip_serial ? "--" : fmt_speedup( bh_ms  / bh_omp_ms ) };
        std::string const winner_col{ "BH " + fmt_speedup( dir_omp_ms / bh_omp_ms ) };

        std::cout << std::left
                  << std::setw( 8 )  << N
                  << std::setw( 8 )  << steps
                  << std::setw( 14 ) << dir_col
                  << std::setw( 14 ) << fmt( dir_omp_ms )
                  << std::setw( 11 ) << dir_spd_col
                  << std::setw( 14 ) << bh_col
                  << std::setw( 14 ) << fmt( bh_omp_ms )
                  << std::setw( 11 ) << bh_spd_col
                  << winner_col
                  << "\n";
    }

    return 0;
}