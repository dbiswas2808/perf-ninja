#include "solution.hpp"

uint16_t checksum(const Blob &blob) {
  std::array<uint16_t, N / 16> sums;
  sums.fill(0);
  for (int ii = 0; ii < N / 16; ++ii) {
    const auto offsetIdx = ii * 16;
    for (int jj = 0; jj < 16; ++jj) {
      sums[jj] += blob[offsetIdx + jj];
      sums[jj] += sums[jj] < blob[offsetIdx + jj];  // add carry
    }
  }

  uint16_t acc = 0;
  for (auto sum : sums) {
    acc += sum;
    acc += acc < sum;  // add carry
  }

  return acc;
}
