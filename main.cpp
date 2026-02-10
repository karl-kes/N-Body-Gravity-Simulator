#include "Classes/Force/Force.hpp"
#include "Classes/Integrator/Integrator.hpp"
#include "Classes/Particle/Particle.hpp"
#include "Classes/Simulation/Simulation.hpp"
#include "Body.hpp"

#include <iostream>
#include <iomanip>
#include <cstddef>
#include <cmath>
#include <string>

int main() {
    /*
        To compile and run:

        python jpl_compare.py fetch --moons --start 1950-01-01
        g++ -std=c++17 -O3 -march=native -fopenmp *.cpp Classes/Force/*.cpp Classes/Integrator/*.cpp Classes/Particle/*.cpp Classes/Simulation/*.cpp -o main.exe
        ./main.exe
        python jpl_compare.py compare
        python visualize.py

        Config lives in Config.hpp:
            dt           = integration timestep (seconds)
            num_years    = simulation duration
            output_hours = output interval in hours
    */

    // Determine number of bodies:
    static constexpr std::size_t num_bodies{ sizeof(bodies) / sizeof(bodies[0]) };

    // Initialize Simulation:
    Simulation sim{ num_bodies, config::total_steps, config::output_interval };
    sim.add_force( std::make_unique<Gravity>() );
    sim.set_integrator( std::make_unique<Yoshida>( config::dt ) );
    initialize_bodies( sim.particles(), num_bodies );

    // Run simulation:
    sim.initial_output();
    sim.run();

    return 0;
}