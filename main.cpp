#include "Vec_3D.hpp"
#include "Constants.hpp"
#include "Body.hpp"
#include "Simulation.hpp"

int main() {
    // Use:
    // g++ main.cpp Vec_3D.cpp Body.cpp Simulation.cpp -o main.exe
    // ./main.exe
    // to compile and run.
    
    std::vector<Body> bodies{};
    Simulation Simulation{ bodies, "bodies.csv" };

    Simulation.load_csv_bodies();
    Simulation.configure_sim();
    Simulation.run_simulation();

    return 0;
}