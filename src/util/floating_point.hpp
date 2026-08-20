#pragma once

#if defined(FP32)
  using fp_t = float;
  inline constexpr auto tolerance{static_cast<fp_t>(1e-5)};
#else
  using fp_t = double;
  inline constexpr auto tolerance{static_cast<fp_t>(1e-10)};
#endif

constexpr fp_t operator""_fp(long double value) {
  return static_cast<fp_t>(value);
}

constexpr fp_t operator""_fp(unsigned long long value) {
  return static_cast<fp_t>(value);
}
