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
        g++ -std=c++17 -O3 -march=native *.cpp Classes/Force/*.cpp Classes/Integrator/*.cpp Classes/Particle/*.cpp Classes/Simulation/*.cpp -o main.exe
        ./main.exe
    */

    // Determine number of bodies:
    static constexpr std::size_t num_bodies{ sizeof(bodies) / sizeof(bodies[0]) };

    // Initialize and prepare simulation:
    Simulation sim{ num_bodies, constant::total_steps, constant::output_interval };
    sim.add_force( std::make_unique<Gravity>() );
    sim.set_integrator( std::make_unique<Yoshida>( constant::dt ) );
    initialize_bodies( sim.particles(), num_bodies );

    // Run simulation:
    sim.initial_output();
    sim.run();
    sim.final_output( bodies );

    return 0;
}