#ifndef HeterogeneousCore_AlpakaInterface_interface_LayerTiles_h
#define HeterogeneousCore_AlpakaInterface_interface_LayerTiles_h

#include <memory>
#include <cmath>
#include <algorithm>
#include <cstdint>

#include "HeterogeneousCore/AlpakaInterface/interface/VecArray.h"

#if !defined(ALPAKA_ACC_GPU_CUDA_ENABLED) && !defined(ALPAKA_ACC_GPU_HIP_ENABLED)
struct int4 {
  int x, y, z, w;
};
#endif

namespace cms::alpakatools {
  // The type T is used to pass the number of bins in each dimension and the
  // allowed ranges spanned. Anchillary quantitied, like the inverse of the bin
  // width should also be provided. Code will not compile if any such
  // information is missing.
  template<typename T, int maxEntriesPerBin>
    class LayerTiles {
      using alpakaVect = VecArray<int, maxEntriesPerBin>;
      public:
        typedef T type;

        template <typename TAcc>
          ALPAKA_FN_ACC inline constexpr void fill(TAcc& acc, const std::vector<float>& dim1, const std::vector<float>& dim2) {
            assert(dim1.size() == dim2.size());
            auto cellsSize = dim1.size();
            for (unsigned int i = 0; i < cellsSize; ++i) {
              layerTiles_[getGlobalBin(dim1[i], dim2[i])].push_back(acc, i);
            }
          }

        template <typename TAcc>
          ALPAKA_FN_ACC inline constexpr void fill(TAcc& acc, float dim1, float dim2, int i) {
            layerTiles_[getGlobalBin(dim1, dim2)].push_back(acc, i);
          }

        ALPAKA_FN_HOST_ACC inline constexpr int getDim1Bin(float dim1) const {
          int dim1Bin = (dim1 - T::minDim1) * T::invDim1BinSize;
          /**
          dim1Bin = (dim1Bin < T::nColumns ? dim1Bin : T::nColumns - 1);
          bool dim1BinPositive = dim1Bin > 0;
          dim1Bin = dim1BinPositive*dim1Bin;
          */
          dim1Bin = std::clamp(dim1Bin, 0, T::nColumns - 1);
          return dim1Bin;
        }

        ALPAKA_FN_HOST_ACC inline constexpr int getDim2Bin(float dim2) const {
          int dim2Bin = (dim2 - T::minDim2) * T::invDim2BinSize;
          /**
          dim2Bin = (dim2Bin < T::nRows ? dim2Bin : T::nRows - 1);
          bool dim2BinPositive = dim2Bin > 0;
          dim2Bin = dim2BinPositive*dim2Bin;
          */
          dim2Bin = std::clamp(dim2Bin, 0, T::nRows - 1);
          return dim2Bin;
        }

        ALPAKA_FN_HOST_ACC inline constexpr int getGlobalBin(float dim1, float dim2) const {
          return getDim1Bin(dim1) + getDim2Bin(dim2) * T::nColumns;
        }

        ALPAKA_FN_HOST_ACC inline constexpr int getGlobalBinByBin(int dim1Bin, int dim2Bin) const {
          return dim1Bin + dim2Bin * T::nColumns;
        }

        ALPAKA_FN_HOST_ACC inline constexpr int4 searchBox(float dim1Min, float dim1Max, float dim2Min, float dim2Max) {
          return int4{getDim1Bin(dim1Min), getDim1Bin(dim1Max), getDim2Bin(dim2Min), getDim2Bin(dim2Max)};
        }

        ALPAKA_FN_HOST_ACC inline constexpr void clear() {
          for (auto& t : layerTiles_)
            t.reset();
        }

        ALPAKA_FN_HOST_ACC inline constexpr void clear(int i) {
          layerTiles_[i].reset();
        }

        ALPAKA_FN_HOST_ACC inline constexpr auto size() {
          return T::nTiles;
        }



        ALPAKA_FN_HOST_ACC inline constexpr const alpakaVect& operator[](int globalBinId) const { return layerTiles_[globalBinId]; }

      private:
        VecArray<VecArray<int, maxEntriesPerBin>, T::nTiles> layerTiles_;
    };

}  // end namespace cms::alpakatools

#endif // HeterogeneousCore_AlpakaInterface_interface_LayerTiles_h
