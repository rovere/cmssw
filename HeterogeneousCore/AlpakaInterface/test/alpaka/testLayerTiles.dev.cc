#include <cstdio>
#include <random>
#include <numeric>

#include <alpaka/alpaka.hpp>

#define CATCH_CONFIG_MAIN
#include <catch.hpp>

#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaInterface/interface/memory.h"
#include "HeterogeneousCore/AlpakaInterface/interface/workdivision.h"
#include "HeterogeneousCore/AlpakaInterface/interface/LayerTiles.h"

#include "DataFormats/Math/interface/constexpr_cmath.h"

// each test binary is built for a single Alpaka backend
using namespace ALPAKA_ACCELERATOR_NAMESPACE;

static constexpr auto s_tag = "[" ALPAKA_TYPE_ALIAS_NAME(alpakaTestKernel) "]";

struct TilesConstants {
  static constexpr float tileSize = 1.f;
  static constexpr float minDim1 = 0.f;
  static constexpr float maxDim1 = 64.f;
  static constexpr float minDim2 = 0.f;
  static constexpr float maxDim2 = 64.f;
  static constexpr int nColumns = reco::ceil((maxDim1 - minDim1) / tileSize);
  static constexpr int nRows = reco::ceil((maxDim2 - minDim2) / tileSize);
  static constexpr float invDim1BinSize = nColumns/(maxDim1 - minDim1);
  static constexpr float invDim2BinSize = nRows/(maxDim2 - minDim2);
  static constexpr int nTiles = nColumns * nRows;
};

struct LayerTilesFillKernel {
  template <typename TAcc, typename T>
  ALPAKA_FN_ACC void operator()(TAcc const& acc, T* __restrict__ out, size_t size) const {
    auto blockIdx  = alpaka::getIdx<alpaka::Grid, alpaka::Blocks>(acc);
    auto elements = alpaka::getWorkDiv<alpaka::Block, alpaka::Elems>(acc);
    for (auto index : cms::alpakatools::elements_with_stride(acc, size)) {
      auto x = T::type::tileSize * blockIdx[0u] + 0.5f;
      auto y = T::type::tileSize * (index % elements[0u]) + 0.5f ;
      out->fill(acc, x, y, index);
    }
  }
};


TEST_CASE("Standard checks of " ALPAKA_TYPE_ALIAS_NAME(alpakaTestKernel), s_tag) {
  SECTION("LayerTilesFillKernel") {
    // get the list of devices on the current platform
    auto const& devices = cms::alpakatools::devices<Platform>();
    if (devices.empty()) {
      std::cout << "No devices available on the platform " << EDM_STRINGIZE(ALPAKA_ACCELERATOR_NAMESPACE)
                << ", the test will be skipped.\n";
      return;
    }

    int constexpr maxEntriesPerBin = 20;
    using Tile = cms::alpakatools::LayerTiles<TilesConstants, maxEntriesPerBin>;

    // allocate input and output host buffers in pinned memory accessible by the Platform devices
    auto layerTiles_h = cms::alpakatools::make_host_buffer<Tile>();

    // run the test on each device
    for (auto const& device : devices) {
      std::cout << "Test VecArray filling on " << alpaka::getName(device) << '\n';
      auto queue = Queue(device);

      // allocate buffer on the device
      auto layerTiles_d = cms::alpakatools::make_device_buffer<Tile>(queue);

      // zero memory on the device
      alpaka::memset(queue, layerTiles_d, 0x00);

      // launch the 1-dimensional kernel with scalar size
      constexpr size_t size = 4096; // 64*64
      auto div = cms::alpakatools::make_workdiv<Acc1D>(64, 64);
      alpaka::exec<Acc1D>(queue, div, LayerTilesFillKernel{}, layerTiles_d.data(), size);

      // copy the results from the device to the host
      alpaka::memcpy(queue, layerTiles_h, layerTiles_d);

      // wait for all the operations to complete
      alpaka::wait(queue);

      // check the results
      auto const& t = *(layerTiles_h.data());
      for (int i = 0; i < 64; ++i) {
        for (int j = 0; j < 64; ++j) {
          auto b = t.getGlobalBin(i + 0.5f, j);
          REQUIRE(t[b][0] == i*64+j);
        }
      }
    }
  }
}
