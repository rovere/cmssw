#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "TrackstersPCA.h"
#include "TPrincipal.h"

#include <iostream>

void ticl::assignPCAtoTracksters(std::vector<Trackster> & tracksters,
    const std::vector<reco::CaloCluster> &layerClusters, double z_limit_em, bool energyWeight) {
  TPrincipal pca(3, "N");
  LogDebug("TrackstersPCA") << "-------" << std::endl;
  for (auto &trackster : tracksters) {
    pca.Clear();
    trackster.raw_energy = 0.;
    trackster.raw_em_energy = 0.;
    trackster.raw_pt = 0.;
    trackster.raw_em_pt = 0.;
    for (size_t i = 0; i < trackster.vertices.size(); ++i) {
      auto fraction = 1.f / trackster.vertex_multiplicity[i];
      trackster.raw_energy += layerClusters[trackster.vertices[i]].energy() * fraction;
      if (std::abs(layerClusters[trackster.vertices[i]].z()) <= z_limit_em)
        trackster.raw_em_energy += layerClusters[trackster.vertices[i]].energy() * fraction;
      if (!energyWeight) {
        double point[3];
        pca.AddRow(vertexToArray(layerClusters[trackster.vertices[i]], point));
      }
    }

    float weights_sum = 0.f;
    if (energyWeight) {
      for (size_t i = 0; i < trackster.vertices.size(); ++i) {
        float weight = 1.f;
        auto fraction = 1.f / trackster.vertex_multiplicity[i];
        if (trackster.raw_energy)
          weight = (trackster.vertices.size() / trackster.raw_energy) * layerClusters[trackster.vertices[i]].energy() * fraction;
        weights_sum += weight;
        double point[3];
        pca.AddRow(vertexToArray(layerClusters[trackster.vertices[i]], point, weight));
      }
    }

    pca.MakePrincipals();
    const auto & barycenter = *(pca.GetMeanValues());
    double sigmasPCA[3] = {0};
    for (size_t i = 0; i < trackster.vertices.size(); ++i) {
      double point[3];
      double p[3];
      pca.X2P(vertexToArray(layerClusters[trackster.vertices[i]], point), p);
      for (size_t j=0; j<3; ++j) {
        sigmasPCA[j] += (p[j] - barycenter[j]) * (p[j] - barycenter[j]);
      }
    }

    // Add trackster attributes
    trackster.barycenter = ticl::Trackster::Vector(barycenter[0],
        barycenter[1],
        barycenter[2]);
    for (size_t i=0; i<3; ++i) {
      sigmasPCA[i] = std::sqrt(sigmasPCA[i]/trackster.vertices.size());
      trackster.sigmas[i] = (float)(*(pca.GetSigmas()))[i];
      trackster.sigmasPCA[i] = sigmasPCA[i];
      trackster.eigenvalues[i] = (float)(*(pca.GetEigenValues()))[i];
      trackster.eigenvectors[i] = ticl::Trackster::Vector((*(pca.GetEigenVectors()))[0][i],
        (*(pca.GetEigenVectors()))[1][i],
        (*(pca.GetEigenVectors()))[2][i] );
    }
    if (trackster.eigenvectors[0].z() * trackster.barycenter.z() < 0.0) {
      trackster.eigenvectors[0] = -ticl::Trackster::Vector((*(pca.GetEigenVectors()))[0][0],
          (*(pca.GetEigenVectors()))[1][0],
          (*(pca.GetEigenVectors()))[2][0] );
    }
    auto norm = std::sqrt(trackster.eigenvectors[0].Unit().perp2());
    trackster.raw_pt = norm * trackster.raw_energy;
    trackster.raw_em_pt = norm * trackster.raw_em_energy;
    const auto & mean = *(pca.GetMeanValues());
    const auto & eigenvectors = *(pca.GetEigenVectors());
    const auto & eigenvalues = *(pca.GetEigenValues());
    const auto & sigmas = *(pca.GetSigmas());

    LogDebug("TrackstersPCA") << "Trackster characteristics: " << std::endl;
    LogDebug("TrackstersPCA") << "Size: " << trackster.vertices.size() << std::endl;
    LogDebug("TrackstersPCA") << "Mean: " << mean[0] << ", " << mean[1] << ", " << mean[2] << std::endl;
    LogDebug("TrackstersPCA") << "Weights sum:" << weights_sum << std::endl;
    LogDebug("TrackstersPCA") << "EigenValues: " << eigenvalues[0] << ", " << eigenvalues[1] << ", " << eigenvalues[2]  << std::endl;
    LogDebug("TrackstersPCA") << "EigeVectors 1: " << eigenvectors(0, 0) << ", " << eigenvectors(1, 0) << ", " << eigenvectors(2, 0) <<std::endl;
    LogDebug("TrackstersPCA") << "EigeVectors 2: " << eigenvectors(0, 1) << ", " << eigenvectors(1, 1) << ", " << eigenvectors(2, 1) <<std::endl;
    LogDebug("TrackstersPCA") << "EigeVectors 3: " << eigenvectors(0, 2) << ", " << eigenvectors(1, 2) << ", " << eigenvectors(2, 2) <<std::endl;
    LogDebug("TrackstersPCA") << "Original sigmas: " << sigmas[0] << ", " << sigmas[1] << ", " << sigmas[2] << std::endl;
    LogDebug("TrackstersPCA") << "SigmasPCA: " << sigmasPCA[0] << ", " << sigmasPCA[1] << ", " << sigmasPCA[2] << std::endl;

    std::cout << "\nTrackster characteristics: " << std::endl;
    std::cout << "Size: " << trackster.vertices.size() << std::endl;
    std::cout << "Energy: " << trackster.raw_energy << std::endl;
    std::cout << "Mean: " << mean[0] << ", " << mean[1] << ", " << mean[2] << std::endl;
    std::cout << "Weights sum:" << weights_sum << std::endl;
    std::cout << "EigenValues: " << eigenvalues[0] << ", " << eigenvalues[1] << ", " << eigenvalues[2]  << std::endl;
    std::cout << "EigeVectors 1: " << eigenvectors(0, 0) << ", " << eigenvectors(1, 0) << ", " << eigenvectors(2, 0) <<std::endl;
    std::cout << "EigeVectors 2: " << eigenvectors(0, 1) << ", " << eigenvectors(1, 1) << ", " << eigenvectors(2, 1) <<std::endl;
    std::cout << "EigeVectors 3: " << eigenvectors(0, 2) << ", " << eigenvectors(1, 2) << ", " << eigenvectors(2, 2) <<std::endl;
    std::cout << "Original sigmas: " << sigmas[0] << ", " << sigmas[1] << ", " << sigmas[2] << std::endl;
    std::cout << "SigmasPCA: " << sigmasPCA[0] << ", " << sigmasPCA[1] << ", " << sigmasPCA[2] << std::endl;
  }
}


double* ticl::vertexToArray(const reco::CaloCluster cluster, double* point, const float weight) {
  point[0] = weight * cluster.x();
  point[1] = weight * cluster.y();
  point[2] = weight * cluster.z();
  return point;
}
