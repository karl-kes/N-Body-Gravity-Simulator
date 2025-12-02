#include "Simulation.hpp"

// Constructor:
Simulation::Simulation( std::vector<Body> &new_bodies,
                        std::string new_file_name,
                        double new_dt):
bodies_{ new_bodies },
file_name{ new_file_name },
dt_{ new_dt } {
    // Empty constructor.
}

// Getters:
const Body &Simulation::get_body( std::size_t idx ) const { 
    return bodies_[idx];
}
const double &Simulation::get_dt() const {
    return dt_;
}
const int &Simulation::get_steps() const {
    return num_steps_;
}
const int &Simulation::get_outputs() const {
    return num_outputs_;
}

// Setters:
void Simulation::set_dt( double const &new_dt ) {
    dt_ = new_dt;
}
void Simulation::set_steps( double const &new_steps ) {
    num_steps_ = new_steps;
}
void Simulation::set_outputs( double const &new_outputs ) {
    num_outputs_ = new_outputs;
}

// Helpers:
double Simulation::calculate_total_energy() const {
    double total_energy{ 0.0 };

    for ( std::size_t idx = 0; idx < bodies_.size(); ++idx ) {
        for ( size_t jdx = idx + 1; jdx < bodies_.size(); ++jdx ) {
            Vec_3D R{ get_body(idx).get_pos() - get_body(jdx).get_pos() };
            double dist{ R.norm() + EPSILON };

            // Potential Energy: V = GMm / R
            total_energy -= G * get_body(idx).get_mass() * get_body(jdx).get_mass() / dist;
        }
        // Kinetic Energy: T = 1/2 mv^2
        total_energy += 0.5 * get_body(idx).get_mass() * get_body(idx).get_vel().norm_squared();
    }
    return total_energy;
}

void Simulation::load_csv_bodies() {
    std::ifstream file( file_name );
    if ( !file.is_open() ) {
        std::cout << "Error opening " << file_name << std::endl;
        return;
    }

    std::string line{};
    while ( std::getline( file, line ) ) {
        if ( line.empty() || line[0] == '#' ) continue;


        std::stringstream string_stream( line );
        std::string segment{};
        std::vector<double> values{};

        while ( std::getline( string_stream, segment, ',' ) ) {
            values.push_back( std::stod( segment ) );
        }

        if ( values.size() == 7 ) {
            
            bodies_.emplace_back( Vec_3D{ values[0], values[1], values[2] }, // Positions
                                    Vec_3D{ values[3], values[4], values[5] }, // Velocities
                                    values[6]                                  // Mass
                                );
        }
    }
    file.close();
}

// Simulation:
void Simulation::configure_sim() {
    double temp_val{};

    std::cout << "<--- N-Body Simulation --->" << std::endl;

    std::cout << "Enter number of steps: ";
    std::cin >> temp_val;
    set_steps( temp_val );

    std::cout << "Enter number of outputs: ";
    std::cin >> temp_val;
    set_outputs( temp_val );

    std::cout << "\nStarting N-Body Simulation..." << std::endl;
}
void Simulation::run_simulation() {
    std::ofstream out_file( "trajectories.csv" );
    out_file << "step,body_id,x,y,z\n";

    // Keep values in scientific notation to 3 sig figs.
    std::cout << std::scientific << std::setprecision( 3 );
    double initial_energy{ calculate_total_energy() };
    double max_energy_drift{};

    auto start_time{ std::chrono::high_resolution_clock::now() };

    #pragma omp parallel
    {
        for( int current_step = 0; current_step < get_steps(); ++current_step ) {
            // Calculates new acceleration.
            #pragma omp for
            for ( std::size_t idx = 0; idx < bodies_.size(); ++idx ) {
                bodies_[idx].calculate_new_acc( bodies_, idx );
            }

            // Updates position and velocity for all bodies.
            #pragma omp for
            for ( std::size_t idx = 0; idx < bodies_.size(); ++idx ) {
                bodies_[idx].update( get_dt() );
            }

            // Outputs information in specific number of outputs.
            int output_interval{ get_steps() / get_outputs() };

            #pragma omp single
            {
                if ( current_step % output_interval == 0 || current_step == get_steps() - 1 ) {
                    // Calculates maximum energy drift in system.
                    double current_energy{ calculate_total_energy() };
                    double energy_drift_percent{ 100.0 * std::abs( current_energy - initial_energy ) / std::abs( initial_energy ) };
                    max_energy_drift = std::max( max_energy_drift, energy_drift_percent );

                    // Outputs the current position for all bodies.
                    for ( std::size_t idx = 0; idx < bodies_.size(); ++idx ) {
                        const Vec_3D &curr_body_pos = bodies_[idx].get_pos();

                        out_file << current_step << "," 
                                << idx << ","
                                << curr_body_pos.get_x() << ","
                                << curr_body_pos.get_y() << ","
                                << curr_body_pos.get_z() << "\n";
                    }
                }

                if ( current_step % 10 == 0 || current_step == get_steps() - 1 ) {
                    float progress = ( ( current_step + 1 ) * 100.0 / get_steps() );

                    std::cout << "Progress: " << std::fixed << std::setprecision( 1 ) 
                                << progress << "%\r" << std::flush;
                }
            }
        }
    }
    out_file.close();
    auto time_elapsed{ ( std::chrono::high_resolution_clock::now() - start_time ).count() * CONVERT_TO_SEC };

    std::cout << std::endl << std::fixed << std::scientific << std::setprecision( 4 );
    std::cout << "\nMax Energy Drift: " << max_energy_drift << "%." << std::endl;
    std::cout << "Time elapsed: " << time_elapsed << " seconds." << std::endl;
    std::cout << "\n<--- End of Simulation --->" << std::endl;
}