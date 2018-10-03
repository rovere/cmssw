#include <Eigen/Core>
#include <Eigen/Eigenvalues>

#include "HeterogeneousCore/CUDAUtilities/interface/cudaCheck.h"
#include "PixelTrackReconstructionGPU.h"

using namespace Eigen;

__global__ void
KernelFastFitAllHits(float *hits_and_covariances,
    int hits_in_fit,
    int cumulative_size,
    float B,
    Rfit::helix_fit *results,
    Rfit::Matrix3xNd *hits,
    Rfit::Matrix3Nd *hits_cov,
    Rfit::circle_fit *circle_fit,
    Vector4d *fast_fit,
    Rfit::line_fit *line_fit)
{
  // Reshape Eigen components from hits_and_covariances, using proper thread and block indices
  // Perform the fit
  // Store the results in the proper vector, using again correct indices

  // Loop for hits_in_fit times:
  //   first 3 are the points
  //   the rest is the covariance matrix, 3x3
  int start = (blockIdx.x * blockDim.x + threadIdx.x) * hits_in_fit * 12;
  int helix_start = (blockIdx.x * blockDim.x + threadIdx.x);
  if (start >= cumulative_size) {
    return;
  }

#ifdef GPU_DEBUG
  printf("BlockDim.x: %d, BlockIdx.x: %d, threadIdx.x: %d, start: %d, cumulative_size: %d\n",
      blockDim.x, blockIdx.x, threadIdx.x, start, cumulative_size);
#endif

  hits[helix_start].resize(3, hits_in_fit);
  hits_cov[helix_start].resize(3 * hits_in_fit, 3 * hits_in_fit);

  // Prepare data structure (stack)
  for (unsigned int i = 0; i < hits_in_fit; ++i) {
    hits[helix_start].col(i) << hits_and_covariances[start],
        hits_and_covariances[start + 1], hits_and_covariances[start + 2];
    start += 3;

    for (auto j = 0; j < 3; ++j) {
      for (auto l = 0; l < 3; ++l) {
        hits_cov[helix_start](i + j * hits_in_fit, i + l * hits_in_fit) =
            hits_and_covariances[start];
        start++;
      }
    }
  }

  Rfit::Fast_fit(hits[helix_start], fast_fit[helix_start]);
}

__global__
void kernelComputeCircleWeights(int hits_in_fit, int cumulative_size, float B,
    Rfit::Matrix3xNd *hits,
    Rfit::Matrix3Nd *hits_cov,
    Vector4d *fast_fit,
    Rfit::covariancesForCircle * cov)
{
  int start = (blockIdx.x * blockDim.x + threadIdx.x) * hits_in_fit * 12;
  int helix_start = (blockIdx.x * blockDim.x + threadIdx.x);
  if (start >= cumulative_size) {
    return;
  }

#ifdef GPU_DEBUG
  printf("blockDim.x: %d, blockIdx.x: %d, threadIdx.x: %d, helix_start: %d, cumulative_size: %d\n",
         blockDim.x, blockIdx.x, threadIdx.x, helix_start, cumulative_size);
#endif
  auto n = hits[helix_start].cols();

  Rfit::VectorNd rad = (hits[helix_start].block(0, 0, 2, n).colwise().norm());

  cov[helix_start].V.resize(2 * n, 2 * n);
  cov[helix_start].V = hits_cov[helix_start].block(0, 0, 2 * n, 2 * n);
  Rfit::printIt(&cov[helix_start].V, "kernelComputeCircleWeights: ");
  Rfit::ComputeCircleWeights(hits[helix_start].block(0, 0, 2, n),
      fast_fit[helix_start], rad, B, cov[helix_start]);
}

__global__
void kernelCircleFitComputeCircleParametersAndErrors(
    int hits_in_fit,
    int cumulative_size,
    Rfit::Matrix3xNd *hits,
    Rfit::circle_fit *circle_fit,
    Rfit::covariancesForCircle * cov)
{
  int start = (blockIdx.x * blockDim.x + threadIdx.x) * hits_in_fit * 12;
  int helix_start = (blockIdx.x * blockDim.x + threadIdx.x);
  if (start >= cumulative_size) {
    return;
  }

#ifdef GPU_DEBUG
  printf("blockDim.x: %d, blockIdx.x: %d, threadIdx.x: %d, helix_start: %d, cumulative_size: %d\n",
         blockDim.x, blockIdx.x, threadIdx.x, helix_start, foundNtuplets->size());
#endif
  auto n = hits[helix_start].cols();

  Rfit::ComputeCircleParametersAndErrors(hits[helix_start].block(0, 0, 2, n),
                                         cov[helix_start],
                                         circle_fit[helix_start]);

}

__global__ void
__launch_bounds__(256, 2)
KernelCircleFitAllHits(int hits_in_fit,
                       int cumulative_size, float B, Rfit::helix_fit *results,
                       Rfit::Matrix3xNd *hits, Rfit::Matrix3Nd *hits_cov,
                       Rfit::circle_fit *circle_fit, Vector4d *fast_fit,
                       Rfit::line_fit *line_fit, Rfit::covariancesForCircle * cov)// Rfit::ArrayNd * vcs, Rfit::MatrixNd * C, Rfit::MatrixNd * D)
{
  // Reshape Eigen components from hits_and_covariances, using proper thread and block indices
  // Perform the fit
  // Store the results in the proper vector, using again correct indices

  // Loop for hits_in_fit times:
  //   first 3 are the points
  //   the rest is the covariance matrix, 3x3
  int start = (blockIdx.x * blockDim.x + threadIdx.x) * hits_in_fit * 12;
  int helix_start = (blockIdx.x * blockDim.x + threadIdx.x);
  if (start >= cumulative_size) {
    return;
  }

#ifdef GPU_DEBUG
    printf("BlockDim.x: %d, BlockIdx.x: %d, threadIdx.x: %d, start: %d, "
           "cumulative_size: %d\n",
           blockDim.x, blockIdx.x, threadIdx.x, start, cumulative_size);
#endif
  u_int n = hits[helix_start].cols();

  Rfit::VectorNd rad = (hits[helix_start].block(0, 0, 2, n).colwise().norm());

  Rfit::Circle_fit(hits[helix_start].block(0, 0, 2, n),
                   hits_cov[helix_start].block(0, 0, 2 * n, 2 * n),
                   fast_fit[helix_start], rad, B,
                   cov[helix_start],
                   circle_fit[helix_start], true);

#ifdef GPU_DEBUG
    printf("KernelCircleFitAllHits circle.par(0): %d %f\n", helix_start,
           circle_fit[helix_start].par(0));
    printf("KernelCircleFitAllHits circle.par(1): %d %f\n", helix_start,
           circle_fit[helix_start].par(1));
    printf("KernelCircleFitAllHits circle.par(2): %d %f\n", helix_start,
           circle_fit[helix_start].par(2));
#endif

}

__global__ void
KernelLineFitAllHits(int hits_in_fit,
                     int cumulative_size, float B, Rfit::helix_fit *results,
                     Rfit::Matrix3xNd *hits, Rfit::Matrix3Nd *hits_cov,
                     Rfit::circle_fit *circle_fit, Vector4d *fast_fit,
                     Rfit::line_fit *line_fit)
{
  // Reshape Eigen components from hits_and_covariances, using proper thread and block indices
  // Perform the fit
  // Store the results in the proper vector, using again correct indices

  // Loop for hits_in_fit times:
  //   first 3 are the points
  //   the rest is the covariance matrix, 3x3
  int start = (blockIdx.x * blockDim.x + threadIdx.x) * hits_in_fit * 12;
  int helix_start = (blockIdx.x * blockDim.x + threadIdx.x);
  if (start >= cumulative_size) {
    return;
  }

#ifdef GPU_DEBUG

  printf("BlockDim.x: %d, BlockIdx.x: %d, threadIdx.x: %d, start: %d, "
      "cumulative_size: %d\n",
      blockDim.x, blockIdx.x, threadIdx.x, start, cumulative_size);
#endif

  Rfit::Line_fit(hits[helix_start], hits_cov[helix_start],
                 circle_fit[helix_start], fast_fit[helix_start], B, line_fit[helix_start], true);

  transformCircleParamsToTrackRepresentation(circle_fit[helix_start], B, true);

  // Grab helix_fit from the proper location in the output vector
  Rfit::helix_fit &helix = results[helix_start];
  helix.par << circle_fit[helix_start].par, line_fit[helix_start].par;

  // TODO: pass properly error booleans

  helix.cov = MatrixXd::Zero(5, 5);
  helix.cov.block(0, 0, 3, 3) = circle_fit[helix_start].cov;
  helix.cov.block(3, 3, 2, 2) = line_fit[helix_start].cov;

  helix.q = circle_fit[helix_start].q;
  helix.chi2_circle = circle_fit[helix_start].chi2;
  helix.chi2_line = line_fit[helix_start].chi2;

#ifdef GPU_DEBUG

    printf("KernelLineFitAllHits line.par(0): %d %f\n", helix_start,
           circle_fit[helix_start].par(0));
    printf("KernelLineFitAllHits line.par(1): %d %f\n", helix_start,
           line_fit[helix_start].par(1));
#endif
}

void PixelTrackReconstructionGPU::launchKernelFit(
    float *hits_and_covariancesGPU, int cumulative_size, int hits_in_fit,
    float B, Rfit::helix_fit *results)
{
  const dim3 threads_per_block(256, 1);
  int num_blocks = cumulative_size / (hits_in_fit * 12) / threads_per_block.x + 1;
  auto numberOfSeeds = cumulative_size / (hits_in_fit * 12);

  Rfit::Matrix3xNd *hitsGPU;
  cudaCheck(cudaMalloc(&hitsGPU, numberOfSeeds * sizeof(Rfit::Matrix3xNd(3, 4))));
  cudaCheck(cudaMemset(hitsGPU, 0x00, numberOfSeeds * sizeof(Rfit::Matrix3xNd(3, 4))));

  Rfit::Matrix3Nd *hits_covGPU = nullptr;
  cudaCheck(cudaMalloc(&hits_covGPU, numberOfSeeds * sizeof(Rfit::Matrix3Nd(12, 12))));
  cudaCheck(cudaMemset(hits_covGPU, 0x00, numberOfSeeds * sizeof(Rfit::Matrix3Nd(12, 12))));

  Vector4d *fast_fit_resultsGPU = nullptr;
  cudaCheck(cudaMalloc(&fast_fit_resultsGPU, numberOfSeeds * sizeof(Vector4d)));
  cudaCheck(cudaMemset(fast_fit_resultsGPU, 0x00, numberOfSeeds * sizeof(Vector4d)));

  Rfit::circle_fit *circle_fit_resultsGPU = nullptr;
  cudaCheck(cudaMalloc(&circle_fit_resultsGPU, numberOfSeeds * sizeof(Rfit::circle_fit)));
  cudaCheck(cudaMemset(circle_fit_resultsGPU, 0x00, numberOfSeeds * sizeof(Rfit::circle_fit)));

  Rfit::line_fit *line_fit_resultsGPU = nullptr;
  cudaCheck(cudaMalloc(&line_fit_resultsGPU, numberOfSeeds * sizeof(Rfit::line_fit)));
  cudaCheck(cudaMemset(line_fit_resultsGPU, 0x00, numberOfSeeds * sizeof(Rfit::line_fit)));

  Rfit::covariancesForCircle * cov_for_circle = nullptr;
  cudaCheck(cudaMalloc(&cov_for_circle, sizeof(Rfit::covariancesForCircle)*numberOfSeeds));
  cudaCheck(cudaMemset(cov_for_circle, 0x00, sizeof(Rfit::covariancesForCircle)*numberOfSeeds));

  KernelFastFitAllHits<<<num_blocks, threads_per_block>>>(
      hits_and_covariancesGPU, hits_in_fit, cumulative_size, B, results,
      hitsGPU, hits_covGPU, circle_fit_resultsGPU, fast_fit_resultsGPU,
      line_fit_resultsGPU);
  cudaCheck(cudaGetLastError());

  kernelComputeCircleWeights<<<num_blocks, threads_per_block>>>(
      hits_in_fit, cumulative_size, B,
      hitsGPU, hits_covGPU, fast_fit_resultsGPU, cov_for_circle);
  cudaCheck(cudaGetLastError());

  KernelCircleFitAllHits<<<num_blocks, threads_per_block>>>(
      hits_in_fit, cumulative_size, B, results,
      hitsGPU, hits_covGPU, circle_fit_resultsGPU, fast_fit_resultsGPU,
      line_fit_resultsGPU, cov_for_circle);
  cudaCheck(cudaGetLastError());

  kernelCircleFitComputeCircleParametersAndErrors<<<num_blocks, threads_per_block>>>(
      hits_in_fit, cumulative_size, hitsGPU, circle_fit_resultsGPU, cov_for_circle);
  cudaCheck(cudaGetLastError());

  KernelLineFitAllHits<<<num_blocks, threads_per_block>>>(
      hits_in_fit, cumulative_size, B, results,
      hitsGPU, hits_covGPU, circle_fit_resultsGPU, fast_fit_resultsGPU,
      line_fit_resultsGPU);
  cudaCheck(cudaGetLastError());

  cudaFree(hitsGPU);
  cudaFree(hits_covGPU);
  cudaFree(fast_fit_resultsGPU);
  cudaFree(circle_fit_resultsGPU);
  cudaFree(line_fit_resultsGPU);
  cudaFree(cov_for_circle);
}
