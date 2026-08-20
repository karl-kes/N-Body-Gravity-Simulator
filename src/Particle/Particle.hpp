#pragma once

#include <xpu/soa.hpp>

#include "../util/enum.hpp"
#include "../util/floating_point.hpp"

class Particles {
private:
  enum class ViewIndex : std::size_t {
    mass = enum_sum<std::size_t>(0uz, 0uz),
    pos = enum_sum<std::size_t>(mass, scalar_count),
    vel = enum_sum<std::size_t>(pos, axis_count),
    acc = enum_sum<std::size_t>(vel, axis_count),
    count = enum_sum<std::size_t>(acc, axis_count)
  };

  xpu::soa<fp_t, to<std::size_t>(ViewIndex::count)> data_;

public:
  explicit Particles(const std::size_t num_particles)
    : data_{num_particles}
  { }

  struct ParticleView {
    xpu::soa_view<fp_t, scalar_count> mass;
    xpu::soa_view<fp_t, axis_count> pos;
    xpu::soa_view<fp_t, axis_count> vel;
    xpu::soa_view<fp_t, axis_count> acc;
    std::size_t count;
  };

  [[nodiscard]] inline
  ParticleView view() {
    return {
      .mass = mass(),
      .pos = pos(),
      .vel = vel(),
      .acc = acc(),
      .count = mass().count()
    };
  }

  [[nodiscard]] inline
  xpu::soa_view<fp_t, scalar_count> mass() {
    return data_.view<scalar_count, to<std::size_t>(ViewIndex::mass)>();
  }

  [[nodiscard]] inline
  xpu::soa_view<fp_t, axis_count> pos() {
    return data_.view<axis_count, to<std::size_t>(ViewIndex::pos)>();
  }

  [[nodiscard]] inline
  xpu::soa_view<fp_t, axis_count> vel() {
    return data_.view<axis_count, to<std::size_t>(ViewIndex::vel)>();
  }

  [[nodiscard]] inline
  xpu::soa_view<fp_t, axis_count> acc() {
    return data_.view<axis_count, to<std::size_t>(ViewIndex::acc)>();
  }
};
