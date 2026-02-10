#pragma once

#include <cstddef>

namespace config {
    // PHYSICAL CONSTANTS - do not edit these:

    inline static constexpr double G{ 6.6743e-11 };
    inline static constexpr double EPS{ 1e-9 };
    inline static constexpr double C{ 299792458.0 }; // m/s
    inline static constexpr double C_SQ{ C * C };
    inline static constexpr double AU{ 1.496e11 }; // m
    inline static constexpr double KM_TO_M{ 1e3 };

    // TIME CONSTANTS - do not edit these:

    inline static constexpr double SECONDS_PER_HOUR{ 3600.0 };
    inline static constexpr double SECONDS_PER_DAY{ 86400.0 };
    inline static constexpr double DAYS_PER_YEAR{ 365.25 };
    inline static constexpr double SECONDS_PER_YEAR{ SECONDS_PER_DAY * DAYS_PER_YEAR };

    // USER CONFIG - change these, everything else adapts:

    inline static constexpr double dt{ 3600.0 };             // Integration timestep (seconds)
    inline static constexpr std::size_t num_years{ 249 };    // Simulation duration (years)
    inline static constexpr std::size_t output_hours{ 487 }; // Output every N hours (must match JPL --step)

    // DERIVED - do not edit these:

    inline static constexpr std::size_t steps_per_year{ static_cast<std::size_t>( SECONDS_PER_YEAR / dt ) };
    inline static constexpr std::size_t total_steps{ steps_per_year * num_years };
    inline static constexpr std::size_t output_interval{ static_cast<std::size_t>( output_hours * SECONDS_PER_HOUR / dt ) };

    inline static constexpr double OMP_THRESHOLD{ 500.0 };

    // COMPILE-TIME CHECK - output_interval must divide evenly:
    static_assert (
        static_cast<std::size_t>( output_hours * SECONDS_PER_HOUR ) % static_cast<std::size_t>( dt ) == 0,
        "output_hours * 3600 must be divisible by dt"
    );
}