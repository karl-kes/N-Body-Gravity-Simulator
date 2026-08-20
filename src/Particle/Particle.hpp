#pragma once

#include <xpu/soa.hpp>

#include "../util/enum.hpp"
#include "../util/floating_point.hpp"

class Particles {
private:
  enum class ArrayIndex : std::size_t {
    mass = 0uz,
    pos = enum_sum<std::size_t>(mass, scalar_count),
    vel = enum_sum<std::size_t>(pos, axis_count),
    acc = enum_sum<std::size_t>(vel, axis_count),
    count = enum_sum<std::size_t>(acc, axis_count)
  };

  xpu::soa<fp_t, to<std::size_t>(ArrayIndex::count)> data_;

public:
  explicit Particles(const std::size_t num_particles)
    : data_{num_particles}
  { }

  [[nodiscard]] CUDA_CALLABLE inline
  xpu::soa_view<fp_t, scalar_count> mass() {
    return data_.view<scalar_count, to<std::size_t>(ArrayIndex::mass)>();
  }

  [[nodiscard]] CUDA_CALLABLE inline
  xpu::soa_view<fp_t, axis_count> pos() {
    return data_.view<axis_count, to<std::size_t>(ArrayIndex::pos)>();
  }

  [[nodiscard]] CUDA_CALLABLE inline
  xpu::soa_view<fp_t, axis_count> vel() {
    return data_.view<axis_count, to<std::size_t>(ArrayIndex::vel)>();
  }

  [[nodiscard]] CUDA_CALLABLE inline
  xpu::soa_view<fp_t, axis_count> acc() {
    return data_.view<axis_count, to<std::size_t>(ArrayIndex::acc)>();
  }
};
