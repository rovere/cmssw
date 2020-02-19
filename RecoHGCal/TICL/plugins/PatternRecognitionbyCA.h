// Author: Felice Pantaleo,Marco Rovere - felice.pantaleo@cern.ch, marco.rovere@cern.ch
// Date: 09/2018

#ifndef __RecoHGCal_TICL_PRbyCA_H__
#define __RecoHGCal_TICL_PRbyCA_H__
#include <memory>  // unique_ptr
#include "RecoHGCal/TICL/plugins/PatternRecognitionAlgoBase.h"
#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"
#include "PhysicsTools/TensorFlow/interface/TensorFlow.h"
#include "TMVA/Factory.h"
#include "TMVA/Reader.h"
#include "CommonTools/MVAUtils/interface/TMVAZipReader.h"


class HGCGraph;

namespace ticl {
  class PatternRecognitionbyCA final : public PatternRecognitionAlgoBase {
  public:
    PatternRecognitionbyCA(const edm::ParameterSet& conf, const CacheBase* cache);
    ~PatternRecognitionbyCA() override;

    void makeTracksters(const PatternRecognitionAlgoBase::Inputs& input, std::vector<Trackster>& result) override;
    void TracksterEmVsHadMVAReader(std::vector<Trackster> & tracksters, TMVA::Reader* reader_);
    void energyRegressionAndID(const std::vector<reco::CaloCluster>& layerClusters, std::vector<Trackster>& result);

  private:
    const std::unique_ptr<HGCGraph> theGraph_;
    const bool out_in_dfs_;
    const unsigned int max_out_in_hops_;
    const float min_cos_theta_;
    const float min_cos_pointing_;
    const int missing_layers_;
    const int min_clusters_per_ntuplet_;
    const float max_delta_time_;
    const std::string eidInputName_;
    const std::string eidOutputNameEnergy_;
    const std::string eidOutputNameId_;
    const float eidMinClusterEnergy_;
    const int eidNLayers_;
    const int eidNClusters_;
    const edm::FileInPath bdtweights_;
    
    TMVA::Reader* reader_;
    float ts_energy_, ts_x_, ts_y_, ts_z_, ts_pcaeigval0_, ts_pcasig0_, ts_pcaeigval1_, ts_pcasig1_, ts_pcaeigval2_, ts_pcasig2_;

    hgcal::RecHitTools rhtools_;
    tensorflow::Session* eidSession_;

    static const int eidNFeatures_ = 3;
  };
}  // namespace ticl
#endif
