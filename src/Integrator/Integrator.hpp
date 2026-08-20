#pragma once

#include "../Particle/Particle.hpp"
#include "../Force/Force.hpp"

#include <vector>
#include <cstddef>
#include <memory>
#include <string>

class Integrator {
protected:
  fp_t dt_;
  std::string name_;
  
public:
  Integrator(fp_t dt, const std::string &name) 
    : dt_{dt}
    , name_{name}
  { }
  
  virtual ~Integrator() = default;
  virtual void integrate(
    Particles::ParticleView particle,
    std::vector<std::unique_ptr<Force>> const &force
  ) const = 0;
  [[nodiscard]] fp_t dt() const { return dt_; }
  [[nodiscard]] std::string const &name() const { return name_; }
};

class Yoshida : public Integrator {
private:
  static constexpr auto cbrt_2_{1.2599210498948732_fp};
  static constexpr auto w_0_{-cbrt_2_ / (2.0_fp - cbrt_2_)};
  static constexpr auto w_1_{1.0_fp / (2.0_fp - cbrt_2_)};

  static constexpr auto c_1_{w_1_ / 2.0_fp};
  static constexpr auto c_2_{(w_0_ + w_1_) / 2.0_fp};
  static constexpr auto c_3_{c_2_};
  static constexpr auto c_4_{c_1_};

  static constexpr auto d_1_{w_1_};
  static constexpr auto d_2_{w_0_};
  static constexpr auto d_3_{w_1_};

public:
  Yoshida(const fp_t dt = 900.0_fp)
    : Integrator{dt, "Yoshida"}
  { }

  void integrate(
    Particles::ParticleView particle,
    const std::vector<std::unique_ptr<Force>>& force
  ) const override;

  [[nodiscard]] static constexpr fp_t cbrt_2() { return cbrt_2_; }
  [[nodiscard]] static constexpr fp_t w_0() { return w_0_; }
  [[nodiscard]] static constexpr fp_t w_1() { return w_1_; }

  [[nodiscard]] static constexpr fp_t c_1() { return c_1_; }
  [[nodiscard]] static constexpr fp_t c_2() { return c_2_; }
  [[nodiscard]] static constexpr fp_t c_3() { return c_3_; }
  [[nodiscard]] static constexpr fp_t c_4() { return c_4_; }

  [[nodiscard]] static constexpr fp_t d_1() { return d_1_; }
  [[nodiscard]] static constexpr fp_t d_2() { return d_2_; }
  [[nodiscard]] static constexpr fp_t d_3() { return d_3_; }
};