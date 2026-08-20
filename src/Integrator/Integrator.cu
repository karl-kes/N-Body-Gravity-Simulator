#include "Integrator.hpp"

#include "../util/enum.hpp"
#include "../util/floating_point.hpp"

#include <xpu/soa.hpp>
#include <xpu/memory.hpp>
#include <omp.h>

void Yoshida::integrate(
  Particles::ParticleView particles,
  std::vector<std::unique_ptr<Force>> const &forces
) const {

  auto calculate_pos = [this, &particles](fp_t c) {
    const auto c_dt{c * dt()};

    for (auto axis{x_axis}; axis < axis_count; ++axis) {
      #pragma omp simd
      for (auto i = 0uz; i < particles.count; ++i) {
        particles.pos[axis][i] += c_dt * particles.vel[axis][i];
      }
    }
  };

  auto apply_force = [&forces, &particles]() {
    for (auto axis{x_axis}; axis < axis_count; ++axis) {
      xpu::zero_n(particles.acc[axis], particles.count);
    }

    for (const auto &force : forces) {
      force->apply(particles);
    };
  };

  auto calculate_vel = [this, &particles](const auto d) {
    const auto d_dt{d * dt()};

    #pragma omp simd
    for (auto axis{x_axis}; axis < axis_count; ++axis) {
      for (auto i = 0uz; i < particles.count; ++i) {
        particles.vel[axis][i] += d_dt * particles.acc[axis][i];
      }
    }
  };

  calculate_pos(c_1());
  apply_force();
  calculate_vel(d_1());

  calculate_pos(c_2());
  apply_force();
  calculate_vel(d_2());

  calculate_pos(c_3());
  apply_force();
  calculate_vel(d_3());

  calculate_pos(c_4());
}
