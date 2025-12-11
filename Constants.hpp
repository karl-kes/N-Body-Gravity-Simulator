#pragma once

static constexpr double G{ 6.67430e-11 };
static constexpr double EPSILON{ 1.0e3 };
static constexpr double MAX_INTERACTION_DIST{ 1e13 };
static constexpr double MAX_INTERACTION_DIST_SQ{ MAX_INTERACTION_DIST * MAX_INTERACTION_DIST };
static constexpr double CONVERT_TO_KMS{ 1e-3 };
static constexpr double CONVERT_TO_KM{ 1e-3 };
static constexpr double CONVERT_TO_SEC{ 1.0e-9 };
static constexpr double SEC_TO_DAY{ 1.0 / 86400 };
static constexpr double SEC_TO_YEAR{ 1.0 / ( 3.154e+7 ) };