#include "Classes/Vec_3D/Vec_3D.hpp"
#include "Classes/Body/Body.hpp"
#include "Classes/Simulation/Simulation.hpp"
#include "Constants.hpp"

int main() {
    /* 
        To compile and run.

        For No Parallel ( <1000 bodies ):
        g++ -std=c++17 main.cpp Classes/Vec_3D/Vec_3D.cpp Classes/Body/Body.cpp Classes/Simulation/Simulation.cpp -o main.exe

        For Parallel ( >1000 bodies ):
        g++ -std=c++17 main.cpp Classes/Vec_3D/Vec_3D.cpp Classes/Body/Body.cpp Classes/Simulation/Simulation.cpp -o main.exe -fopenmp

        ./main.exe
    */
    
    std::vector<Body> bodies{};
    Simulation sim{ bodies, "bodies.csv", 1000 };

    sim.load_csv_bodies();
    sim.configure_sim();
    sim.run_simulation();

    return 0;
}