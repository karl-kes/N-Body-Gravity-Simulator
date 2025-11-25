#include <iostream>
#include <cassert>
#include <cmath>
#include <iomanip>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <chrono>
#include <omp.h>

static constexpr double G{ 6.67430e-11 };
static constexpr double EPSILON{ 1.0e5 };
static constexpr double CONVERT_TO_KMS{ 1e-3 };
static constexpr double CONVERT_TO_KM{ 1e-3 };
static constexpr double CONVERT_TO_SEC{ 1.0e-9 };

class Vec_3D {
private:
    double x_;
    double y_;
    double z_;

public:
    Vec_3D ( double new_x = 0, double new_y = 0, double new_z = 0 ):
    x_{ new_x },
    y_{ new_y },
    z_{ new_z } {
        // Empty constructor
    }

    // Getters:
    const double &get_x() const {
        return x_;
    }
    const double &get_y() const {
        return y_;
    }
    const double &get_z() const {
        return z_;
    }

    // Setters:
    void set_x( double new_x ) {
        x_ = new_x;
    }
    void set_y( double new_y ) {
        y_ = new_y;
    }
    void set_z( double new_z ) {
        z_ = new_z;
    }

    // Helper functions:
    double norm() const {
        return std::sqrt( get_x()*get_x() + get_y()*get_y() + get_z()*get_z() );
    }
    double norm_squared() const {
        return get_x()*get_x() + get_y()*get_y() + get_z()*get_z();
    }

    // Operator overloads:
    Vec_3D operator+( Vec_3D const &other_vec ) const {
        return { get_x() + other_vec.get_x(), get_y() + other_vec.get_y(), get_z() + other_vec.get_z() };
    }
    Vec_3D operator-( Vec_3D const &other_vec ) const {
        return { get_x() - other_vec.get_x(), get_y() - other_vec.get_y(), get_z() - other_vec.get_z() };
    }
    Vec_3D operator*( double const &constant ) const {
        return { get_x()*constant, get_y()*constant, get_z()*constant };
    }
    Vec_3D &operator+=( Vec_3D const &other_vec ) {
        set_x( get_x() + other_vec.get_x() );
        set_y( get_y() + other_vec.get_y() );
        set_z( get_z() + other_vec.get_z() );
        return *this;
    }
};

class Body {
private:
    Vec_3D pos_;
    Vec_3D vel_;
    Vec_3D acc_;
    Vec_3D old_acc_;
    double mass_;

public:
    Body ( Vec_3D new_pos, Vec_3D new_vel, double new_mass ):
    pos_( new_pos ),
    vel_( new_vel ), 
    acc_{ 0, 0, 0 }, 
    old_acc_{ 0, 0, 0 }, 
    mass_( new_mass ) {
        // Empty constructor.
    }
    
    // Calculates new acceleration based on forces from other bodies.
    void calculate_new_acc( std::vector<Body> const &other_bodies, std::size_t const &self_idx ) {
        set_old_acc( acc_ );
        Vec_3D total_force{};

        for ( std::size_t idx = 0; idx < other_bodies.size(); ++idx ) {
            if ( idx == self_idx ) continue;
            
            Vec_3D R{ other_bodies[idx].get_pos() - get_pos() };
            double dist_squared{ R.norm_squared() + EPSILON * EPSILON };
            double force_mag{ ( G * mass_ * other_bodies[idx].get_mass() ) / ( dist_squared ) };

            double dist{ R.norm() };
            total_force += R * ( force_mag / dist );
        }
        set_acc( total_force * ( 1.0 / get_mass() ) );
    }

    // Updates body.
    void update( double const &dt ) {
        set_pos( get_pos() + ( get_vel() * dt + get_acc() * ( 0.5 * dt * dt ) ) );
        set_vel( get_vel() + ( ( get_acc() + get_old_acc() ) * ( 0.5 * dt ) ) );
    }

    // Body getters:
    const Vec_3D &get_pos() const { return pos_; }
    const Vec_3D &get_vel() const { return vel_; }
    const Vec_3D &get_acc() const { return acc_; }
    const Vec_3D &get_old_acc() const { return old_acc_; }
    double get_mass() const { return mass_; }

    // Body setters:
    void set_pos( Vec_3D new_pos ) {
        pos_ = new_pos;
    }
    void set_vel( Vec_3D new_vel ) {
        vel_ = new_vel;
    }
    void set_acc( Vec_3D new_acc ) {
        acc_ = new_acc;
    }
    void set_old_acc( Vec_3D new_old_acc ) {
        old_acc_ = new_old_acc;
    }
};

// Calculates total energy of all bodies.
double calculate_total_energy( std::vector<Body> const &bodies ) {
    double T{ 0.0 };
    double V{ 0.0 };

    for ( std::size_t i = 0; i < bodies.size(); ++i ) {
        for ( size_t j = i + 1; j < bodies.size(); ++j ) {
            Vec_3D R{ bodies[i].get_pos() - bodies[j].get_pos() };
            double dist{ R.norm() + EPSILON };

            V -= G * bodies[i].get_mass() * bodies[j].get_mass() / dist;
        }
        T += 0.5 * bodies[i].get_mass() * bodies[i].get_vel().norm_squared();
    }
    
    return T + V;
}

// Load bodies.
void load_csv_bodies( std::string const &file_name, std::vector<Body> &bodies ) {
    std::ifstream file( file_name );
    if ( !file.is_open() ) {
        std::cout << "Error opening " << file_name << ".csv." << std::endl;
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
            bodies.emplace_back( Vec_3D{ values[0], values[1], values[2] }, // Positions
                                 Vec_3D{ values[3], values[4], values[5] }, // Velocities
                                 values[6]                                  // Mass
                                );
        }
    }
    file.close();
}

int main() {
    double dt{ 1000 };
    int steps{};
    int num_outputs{};

    std::vector<Body> bodies{};
    load_csv_bodies( "bodies.csv", bodies );

    std::cout << "<--- N-Body Simulation --->" << std::endl;

    std::cout << "Enter number of steps: ";
    std::cin >> steps;
    std::cout << "Enter number of outputs: ";
    std::cin >> num_outputs;

    std::cout << "\nStarting N-Body Simulation..." << std::endl;
    std::ofstream out_file( "trajectories.csv" );
    out_file << "step,body_id,x,y,z\n";

    // Keep values in scientific notation to 3 sig figs.
    std::cout << std::scientific << std::setprecision( 3 );
    double initial_energy{ calculate_total_energy( bodies ) };
    double max_energy_drift{};

    auto start_time{ std::chrono::high_resolution_clock::now() };

    #pragma omp parallel
    {
        for( int current_step = 0; current_step < steps; ++current_step ) {
            // Calculates new acceleration.
            #pragma omp for
            for ( std::size_t idx = 0; idx < bodies.size(); ++idx ) {
                bodies[idx].calculate_new_acc( bodies, idx );
            }

            // Updates position and velocity for all bodies.
            #pragma omp for
            for ( std::size_t idx = 0; idx < bodies.size(); ++idx ) {
                bodies[idx].update( dt );
            }

            // Outputs information in specific number of outputs.
            int output_interval{ steps / num_outputs };

            #pragma omp single
            {
                if ( current_step % output_interval == 0 || current_step == steps - 1 ) {
                    // Calculates maximum energy drift in system.
                    double current_energy{ calculate_total_energy( bodies ) };
                    double energy_drift_percent{ 100.0 * std::abs( current_energy - initial_energy ) / std::abs( initial_energy ) };
                    max_energy_drift = std::max( max_energy_drift, energy_drift_percent );

                    // Outputs the current position for all bodies.
                    for ( std::size_t idx = 0; idx < bodies.size(); ++idx ) {
                        const Vec_3D &curr_body_pos = bodies[idx].get_pos();

                        out_file << current_step << "," 
                                 << idx << ","
                                 << curr_body_pos.get_x() << ","
                                 << curr_body_pos.get_y() << ","
                                 << curr_body_pos.get_z() << "\n";
                    }
                }

                if ( current_step % 10 == 0 || current_step == steps - 1 ) {
                    float progress = ( ( current_step + 1 ) * 100.0 / steps );

                    std::cout << "Progress: " << std::fixed << std::setprecision(1) 
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

    return 0;
}