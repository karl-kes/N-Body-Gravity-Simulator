#pragma once

#include <cstddef>

namespace constant {
    // Physical Constants:
    inline static constexpr double G{ 6.6743e-11 };
    inline static constexpr double EPS{ 1e-9 };
    inline static constexpr double AU{ 1.496e11 }; // m
    inline static constexpr double KM_TO_M{ 1e3 }; // km to m

    // Parameters:
    inline static constexpr double dt{ 3600 }; // s
    inline static constexpr std::size_t steps_per_year{ 8766 };
    inline static constexpr std::size_t num_years{ 1 };
    inline static constexpr std::size_t total_steps{ steps_per_year * num_years };
    inline static constexpr std::size_t output_interval{ steps_per_year };
}