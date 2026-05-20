#pragma once

#include <cstddef>

namespace config {
    // PHYSICAL CONSTANTS:

    inline static constexpr double G{ 6.6743e-11 };
    inline static constexpr double AU{ 1.496e11 };
    inline static constexpr double KM_TO_M{ 1e3 };

    // TIME UNIT CONVERSIONS (used by IC generators and validators):

    inline static constexpr double SECONDS_PER_HOUR{ 3600.0 };
    inline static constexpr double SECONDS_PER_DAY{ 86400.0 };
    inline static constexpr double DAYS_PER_YEAR{ 365.25 };
    inline static constexpr double SECONDS_PER_YEAR{ SECONDS_PER_DAY * DAYS_PER_YEAR };

    // PARALLELISM:

    inline static constexpr std::size_t OMP_THRESHOLD{ 350 };
}