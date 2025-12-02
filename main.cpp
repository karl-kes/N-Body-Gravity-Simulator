#include "Classes/Vec_3D/Vec_3D.hpp"
#include "Classes/Body/Body.hpp"
#include "Classes/Simulation/Simulation.hpp"
#include "Constants.hpp"

int main() {
    // Use (No Parallel):
    // g++ main.cpp Classes/Vec_3D/Vec_3D.cpp Classes/Body/Body.cpp Classes/Simulation/Simulation.cpp -o main.exe
    // OR (Parallel):
    // g++ main.cpp Classes/Vec_3D/Vec_3D.cpp Classes/Body/Body.cpp Classes/Simulation/Simulation.cpp -o main.exe -fopenmp
    // ./main.exe
    // to compile and run.
    
    std::vector<Body> bodies{};
    Simulation Simulation{ bodies, "bodies.csv" };

    Simulation.load_csv_bodies();
    Simulation.configure_sim();
    Simulation.run_simulation();

    return 0;
}