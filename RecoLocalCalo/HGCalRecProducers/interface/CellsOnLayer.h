struct CellsOnLayer {
  std::vector<DetId> detid;
  std::vector<float> dim1;
  std::vector<float> dim2;

  std::vector<float> weight;
  std::vector<float> rho;

  std::vector<float> delta;
  std::vector<int> nearestHigher;
  std::vector<int> clusterIndex;
  std::vector<float> sigmaNoise;
  std::vector<std::vector<int>> followers;
  std::vector<bool> isSeed;

  void clear() {
    detid.clear();
    dim1.clear();
    dim2.clear();
    weight.clear();
    rho.clear();
    delta.clear();
    nearestHigher.clear();
    clusterIndex.clear();
    sigmaNoise.clear();
    followers.clear();
    isSeed.clear();
  }

  void shrink_to_fit() {
    detid.shrink_to_fit();
    dim1.shrink_to_fit();
    dim2.shrink_to_fit();
    weight.shrink_to_fit();
    rho.shrink_to_fit();
    delta.shrink_to_fit();
    nearestHigher.shrink_to_fit();
    clusterIndex.shrink_to_fit();
    sigmaNoise.shrink_to_fit();
    followers.shrink_to_fit();
    isSeed.shrink_to_fit();
  }
};

struct CellsOnLayers {
  std::vector<DetId> detid;
  std::vector<float> dim1;
  std::vector<float> dim2;
  std::vector<int> layer;

  std::vector<float> weight;
  std::vector<float> rho;

  std::vector<float> delta;
  std::vector<int> nearestHigher;
  std::vector<int> clusterIndex;
  std::vector<float> sigmaNoise;
  std::vector<std::vector<int>> followers;
  std::vector<uint8_t> isSeed;

  void clear() {
    detid.clear();
    dim1.clear();
    dim2.clear();
    layer.clear();
    weight.clear();
    rho.clear();
    delta.clear();
    nearestHigher.clear();
    clusterIndex.clear();
    sigmaNoise.clear();
    followers.clear();
    isSeed.clear();
  }

  void shrink_to_fit() {
    detid.shrink_to_fit();
    dim1.shrink_to_fit();
    dim2.shrink_to_fit();
    layer.shrink_to_fit();
    weight.shrink_to_fit();
    rho.shrink_to_fit();
    delta.shrink_to_fit();
    nearestHigher.shrink_to_fit();
    clusterIndex.shrink_to_fit();
    sigmaNoise.shrink_to_fit();
    followers.shrink_to_fit();
    isSeed.shrink_to_fit();
  }
};