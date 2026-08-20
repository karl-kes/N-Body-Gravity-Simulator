#pragma once

#include <xpu/config.hpp>

template <typename C, typename T> [[nodiscard]] CUDA_CALLABLE
inline constexpr C to(T arg) noexcept {
  return static_cast<C>(arg)
}

template <typename C, typename T> [[nodiscard]] CUDA_CALLABLE
inline constexpr C enum_sum(T a, T b) noexcept {
  return static_cast<C>(a) + static_cast<C>(b);
}

enum class Axis : std::size_t { x, y, z, count };

inline constexpr auto x_axis{to<std::size_t>(Axis::x)};
inline constexpr auto y_axis{to<std::size_t>(Axis::y)};
inline constexpr auto z_axis{to<std::size_t>(Axis::z)};
inline constexpr auto axis_count{to<std::size_t>(Axis::count)};
inline constexpr auto scalar_count{1uz};