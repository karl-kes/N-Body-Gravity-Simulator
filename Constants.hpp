#pragma once

#include <cstddef>

namespace constant {
    // Physical Constants:
    inline static constexpr double G{ 6.6743e-11 };
    inline static constexpr double EPS{ 1e-9 };
    inline static constexpr double AU{ 1.496e11 }; // m
    inline static constexpr double KM_TO_M{ 1e3 }; // km to m
\
    inline static constexpr double SECONDS_PER_HOUR{ 3600.0 };
    inline static constexpr double SECONDS_PER_DAY{ 86400.0 };
    inline static constexpr double DAYS_PER_YEAR{ 365.25 };
    inline static constexpr double SECONDS_PER_YEAR{ SECONDS_PER_DAY * DAYS_PER_YEAR };

    // Parameters:
    inline static constexpr double dt{ SECONDS_PER_HOUR };
    inline static constexpr std::size_t steps_per_year{ static_cast<std::size_t>( SECONDS_PER_YEAR / dt ) };
    inline static constexpr std::size_t num_years{ 10 };
    inline static constexpr std::size_t total_steps{ steps_per_year * num_years };
    inline static constexpr std::size_t output_interval{ steps_per_year };
}