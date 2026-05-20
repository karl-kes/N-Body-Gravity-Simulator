#include "Force/Force.hpp"
#include "Force/BarnesHut.hpp"
#include "Integrator/Integrator.hpp"
#include "Particle/Particle.hpp"
#include "Simulation/Simulation.hpp"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdint>
#include <cstddef>
#include <cmath>
#include <cstring>
#include <string>
#include <string_view>
#include <vector>

// Binary IC format (all little-endian, double = IEEE 754):
//
//   Header (40 bytes):
//     uint64_t N           # particle count
//     double   dt          # integration timestep (s)
//     double   total_time  # integration duration (s)
//     double   output_dt   # output cadence (s)
//     double   eps         # softening length (m)
//
//   Names (32 * N bytes):
//     char[32] per particle, null-padded
//
//   Particle data (56 * N bytes):
//     double mass, x, y, z, vx, vy, vz   per particle (all SI)

int main( int argc, char* argv[] ) {
    std::string_view force_kind{ "direct" };
    double theta{ 0.5 };
    std::string ic_file{ "tests/sim_ic.bin" };
    std::string out_file{ "tests/sim_output.bin" };

    for ( int i{ 1 }; i < argc; ++i ) {
        std::string_view const arg{ argv[i] };
        if      ( arg == "--force" && i + 1 < argc ) { force_kind = argv[++i]; }
        else if ( arg == "--theta" && i + 1 < argc ) { theta = std::stod( argv[++i] ); }
        else if ( arg == "--ic"    && i + 1 < argc ) { ic_file = argv[++i]; }
        else if ( arg == "--out"   && i + 1 < argc ) { out_file = argv[++i]; }
        else if ( arg == "-h" || arg == "--help" ) {
            std::cout << "Usage: main [--ic FILE] [--out FILE] [--force {direct|bh}] [--theta T]\n"
                      << "  --ic FILE       Binary IC file (default: tests/sim_ic.bin)\n"
                      << "  --out FILE      Output binary path (default: tests/sim_output.bin)\n"
                      << "  --force direct  Direct O(N^2) summation (default)\n"
                      << "  --force bh      Barnes-Hut O(N log N) approximation\n"
                      << "  --theta T       Opening angle for BH (default 0.5)\n";
            return 0;
        }
    }

    if ( force_kind != "direct" && force_kind != "bh" ) {
        std::cerr << "error: --force must be 'direct' or 'bh'\n";
        return 1;
    }

    std::ifstream f{ ic_file, std::ios::binary };
    if ( !f ) {
        std::cerr << "error: failed to open IC file: " << ic_file << "\n"
                  << "       run jpl_compare.py fetch (solar) or generate_galaxies.py (galaxy) first\n";
        return 1;
    }

    // Header.
    uint64_t N{};
    double dt{}, total_time{}, output_dt{}, eps{};
    f.read( reinterpret_cast<char*>( &N ),          sizeof( N ) );
    f.read( reinterpret_cast<char*>( &dt ),         sizeof( dt ) );
    f.read( reinterpret_cast<char*>( &total_time ), sizeof( total_time ) );
    f.read( reinterpret_cast<char*>( &output_dt ),  sizeof( output_dt ) );
    f.read( reinterpret_cast<char*>( &eps ),        sizeof( eps ) );

    if ( !f || N == 0 || dt <= 0.0 || total_time <= 0.0 || output_dt <= 0.0 ) {
        std::cerr << "error: malformed IC header\n";
        return 1;
    }

    std::size_t const num_bodies{ static_cast<std::size_t>( N ) };
    std::size_t const total_steps{ static_cast<std::size_t>( total_time / dt ) };
    std::size_t const output_interval{ static_cast<std::size_t>( output_dt / dt ) };

    if ( output_interval == 0 ) {
        std::cerr << "error: output_dt (" << output_dt << "s) is smaller than dt (" << dt << "s)\n";
        return 1;
    }

    // Names block: 32 bytes per particle, null-padded.
    std::vector<std::string> names{};
    names.reserve( num_bodies );
    {
        std::vector<char> buf( 32 );
        for ( std::size_t i{}; i < num_bodies; ++i ) {
            f.read( buf.data(), 32 );
            auto const end{ std::find( buf.begin(), buf.end(), '\0' ) };
            names.emplace_back( buf.data(), static_cast<std::size_t>( end - buf.begin() ) );
        }
    }

    Simulation sim{
        num_bodies, total_steps, output_interval, dt, eps,
        std::move( names ), std::move( out_file )
    };

    // Particle data: 7 doubles per particle (mass, x, y, z, vx, vy, vz).
    // Single bulk read amortizes syscall overhead at large N.
    constexpr std::size_t doubles_per_p{ 7 };
    std::vector<double> data( num_bodies * doubles_per_p );
    f.read( reinterpret_cast<char*>( data.data() ),
            static_cast<std::streamsize>( data.size() * sizeof( double ) ) );

    for ( std::size_t i{}; i < num_bodies; ++i ) {
        std::size_t const o{ i * doubles_per_p };
        sim.particles().mass()[i]  = data[o + 0];
        sim.particles().pos_x()[i] = data[o + 1];
        sim.particles().pos_y()[i] = data[o + 2];
        sim.particles().pos_z()[i] = data[o + 3];
        sim.particles().vel_x()[i] = data[o + 4];
        sim.particles().vel_y()[i] = data[o + 5];
        sim.particles().vel_z()[i] = data[o + 6];
        sim.particles().acc_x()[i] = 0.0;
        sim.particles().acc_y()[i] = 0.0;
        sim.particles().acc_z()[i] = 0.0;
    }

    if ( force_kind == "bh" ) {
        sim.add_force( std::make_unique<Gravity_BarnesHut>( theta, 8, eps ) );
    } else {
        sim.add_force( std::make_unique<Gravity>( eps ) );
    }
    sim.set_integrator( std::make_unique<Yoshida>( dt ) );
    sim.run();

    return 0;
}