#pragma once

#include "../../Constants.hpp"
#include "../Body/Body.hpp"
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <chrono>
#include <omp.h>

class Simulation {
private:
    std::vector<Body> bodies_;
    std::string file_name;
    double dt_;
    int num_steps_;
    int num_outputs_;

public:
    Simulation( std::vector<Body> &new_bodies, std::string new_file_name = "bodies.csv", double new_dt = 500 );

    // Getters:
    const Body &get_body( std::size_t idx ) const;
    const double &get_dt() const;
    const int &get_steps() const;
    const int &get_outputs() const;

    // Setters:
    void set_dt( double const &new_dt );
    void set_steps( double const &new_steps );
    void set_outputs( double const &new_outputs );

    // Helpers:
    double calculate_total_energy() const;

    void load_csv_bodies();

    // Simulation:
    void configure_sim();
    void run_simulation();
};