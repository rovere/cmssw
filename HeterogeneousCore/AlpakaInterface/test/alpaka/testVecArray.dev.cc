#include <cstdio>
#include <random>
#include <numeric>

#include <alpaka/alpaka.hpp>

#define CATCH_CONFIG_MAIN
#include <catch.hpp>

#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaInterface/interface/memory.h"
#include "HeterogeneousCore/AlpakaInterface/interface/workdivision.h"
#include "HeterogeneousCore/AlpakaInterface/interface/VecArray.h"

// each test binary is built for a single Alpaka backend
using namespace ALPAKA_ACCELERATOR_NAMESPACE;

static constexpr auto s_tag = "[" ALPAKA_TYPE_ALIAS_NAME(alpakaTestKernel) "]";

struct VecArrayFillKernel {
  template <typename TAcc, typename T>
  ALPAKA_FN_ACC void operator()(TAcc const& acc, T* __restrict__ out, size_t size) const {
    for (auto index : cms::alpakatools::elements_with_stride(acc, size)) {
      out->push_back(acc, index);
    }
  }
};

TEST_CASE("Standard checks of " ALPAKA_TYPE_ALIAS_NAME(alpakaTestKernel), s_tag) {
  SECTION("VecArrayFillKernel") {
    // get the list of devices on the current platform
    auto const& devices = cms::alpakatools::devices<Platform>();
    if (devices.empty()) {
      std::cout << "No devices available on the platform " << EDM_STRINGIZE(ALPAKA_ACCELERATOR_NAMESPACE)
                << ", the test will be skipped.\n";
      return;
    }

    // buffer size
    constexpr size_t size = 1024 * 1024;

    // allocate input and output host buffers in pinned memory accessible by the Platform devices
    auto vecArr_h = cms::alpakatools::make_host_buffer<cms::alpakatools::VecArray<unsigned int, size>>();

    // initialize original index locations
    std::vector<size_t> idx(size);
    std::iota(idx.begin(), idx.end(), 0);

    // run the test on each device
    for (auto const& device : devices) {
      std::cout << "Test VecArray filling on " << alpaka::getName(device) << '\n';
      auto queue = Queue(device);

      // allocate buffer on the device
      auto vecArr_d = cms::alpakatools::make_device_buffer<cms::alpakatools::VecArray<unsigned int, size>>(queue);

      // zero memory on the device
      alpaka::memset(queue, vecArr_d, 0x00);

      // launch the 1-dimensional kernel with scalar size
      auto div = cms::alpakatools::make_workdiv<Acc1D>(64, 64);
      alpaka::exec<Acc1D>(queue, div, VecArrayFillKernel{}, vecArr_d.data(), size);

      // copy the results from the device to the host
      alpaka::memcpy(queue, vecArr_h, vecArr_d);

      // wait for all the operations to complete
      alpaka::wait(queue);

      std::stable_sort(idx.begin(), idx.end(), [&vecArr_h](size_t i1, size_t i2) {
        auto const& v = *(vecArr_h.data());
        return v[i1] < v[i2];
      });

      // check the results
      auto const& v = *(vecArr_h.data());
      for (size_t i = 0; i < size; ++i) {
        REQUIRE(v[idx[i]] == i);
      }
    }
  }
}
