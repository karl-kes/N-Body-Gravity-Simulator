#include "Force.hpp"
#include "Integrator.hpp"
#include "Particle.hpp"
#include "Simulation.hpp"

#include <iostream>

int main() {
    /*
        To compile and run:
        g++ -std=c++17 -O3 -march=native *.cpp -o main.exe
        ./main.exe
    */
    
    // Physical constants
    constexpr double M_sun{ 1.989e30 };      // kg
    constexpr double M_earth{ 5.972e24 };    // kg
    constexpr double AU{ 1.496e11 };         // m
    constexpr double v_earth{ 29780.0 };     // m/s (orbital velocity)

    // Simulation parameters
    constexpr double dt{ 3600.0 };           // 1 hour timestep
    constexpr std::size_t steps_per_year{ 8766 };
    constexpr std::size_t num_years{ 10 };
    constexpr std::size_t total_steps{ steps_per_year * num_years };
    constexpr std::size_t output_interval{ steps_per_year };

    std::cout << "=== Earth-Sun Two-Body Simulation ===" << std::endl;
    std::cout << "Integrator: Velocity Verlet (2nd order)" << std::endl;
    std::cout << "Timestep: " << dt << " s" << std::endl;
    std::cout << "Duration: " << num_years << " years (" << total_steps << " steps)" << std::endl;
    std::cout << std::endl;

    // Create simulation with 2 particles
    Simulation sim{ 2, total_steps, output_interval };

    // Add gravity
    sim.add_force( std::make_unique<Gravity>() );

    // Set integrator
    sim.set_integrator( std::make_unique<Velocity_Verlet>( dt ) );

    // Initialize Sun (at origin, stationary)
    sim.particles().mass()[0] = M_sun;
    sim.particles().pos_x()[0] = 0.0;
    sim.particles().pos_y()[0] = 0.0;
    sim.particles().pos_z()[0] = 0.0;
    sim.particles().vel_x()[0] = 0.0;
    sim.particles().vel_y()[0] = 0.0;
    sim.particles().vel_z()[0] = 0.0;

    // Initialize Earth (1 AU from Sun, circular orbit)
    sim.particles().mass()[1] = M_earth;
    sim.particles().pos_x()[1] = AU;
    sim.particles().pos_y()[1] = 0.0;
    sim.particles().pos_z()[1] = 0.0;
    sim.particles().vel_x()[1] = 0.0;
    sim.particles().vel_y()[1] = v_earth;
    sim.particles().vel_z()[1] = 0.0;

    // Run simulation
    sim.run();

    // Print final position
    double final_dist = std::sqrt( 
        sim.particles().pos_x( 1 ) * sim.particles().pos_x( 1 ) +
        sim.particles().pos_y( 1 ) * sim.particles().pos_y( 1 ) +
        sim.particles().pos_z( 1 ) * sim.particles().pos_z( 1 )
    );
    std::cout << "Final Earth distance: " << final_dist / AU << " AU" << std::endl;

    // Test Yoshida integrator
    std::cout << "\n=== Testing Yoshida 4th Order ===" << std::endl;
    std::cout << "Timestep: " << dt * 10 << " s (10x larger)" << std::endl;

    Simulation sim2{ 2, total_steps / 10, output_interval / 10 };
    sim2.add_force( std::make_unique<Gravity>() );
    sim2.set_integrator( std::make_unique<Yoshida>( dt * 10 ) );

    // Same initial conditions
    sim2.particles().mass()[0] = M_sun;
    sim2.particles().pos_x()[0] = 0.0;
    sim2.particles().pos_y()[0] = 0.0;
    sim2.particles().pos_z()[0] = 0.0;
    sim2.particles().vel_x()[0] = 0.0;
    sim2.particles().vel_y()[0] = 0.0;
    sim2.particles().vel_z()[0] = 0.0;

    sim2.particles().mass()[1] = M_earth;
    sim2.particles().pos_x()[1] = AU;
    sim2.particles().pos_y()[1] = 0.0;
    sim2.particles().pos_z()[1] = 0.0;
    sim2.particles().vel_x()[1] = 0.0;
    sim2.particles().vel_y()[1] = v_earth;
    sim2.particles().vel_z()[1] = 0.0;

    sim2.run();

    double final_dist2 = std::sqrt( 
        sim2.particles().pos_x( 1 ) * sim2.particles().pos_x( 1 ) +
        sim2.particles().pos_y( 1 ) * sim2.particles().pos_y( 1 ) +
        sim2.particles().pos_z( 1 ) * sim2.particles().pos_z( 1 )
    );
    std::cout << "\nFinal Earth distance: " << final_dist2 / AU << " AU" << std::endl;

    return 0;
}