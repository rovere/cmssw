#ifndef CkfBaseTrajectoryFilter_H
#define CkfBaseTrajectoryFilter_H

#include "TrackingTools/TrajectoryFiltering/interface/TrajectoryFilter.h"

#include "TrackingTools/TrajectoryFiltering/interface/ChargeSignificanceTrajectoryFilter.h"
#include "TrackingTools/TrajectoryFiltering/interface/MaxConsecLostHitsTrajectoryFilter.h"
#include "TrackingTools/TrajectoryFiltering/interface/MaxHitsTrajectoryFilter.h"
#include "TrackingTools/TrajectoryFiltering/interface/MaxLostHitsTrajectoryFilter.h"
#include "TrackingTools/TrajectoryFiltering/interface/MinHitsTrajectoryFilter.h"
#include "TrackingTools/TrajectoryFiltering/interface/MinPtTrajectoryFilter.h"
#include "TrackingTools/TrajectoryFiltering/interface/LostHitsFractionTrajectoryFilter.h"
#include "TrackingTools/TrajectoryFiltering/interface/LooperTrajectoryFilter.h"


class CkfBaseTrajectoryFilter : public TrajectoryFilter {
public:

  explicit CkfBaseTrajectoryFilter( const edm::ParameterSet & pset, edm::ConsumesCollector& iC):
    //define the filters by default in the BaseCkfTrajectoryBuilder
    theChargeSignificanceTrajectoryFilter(new ChargeSignificanceTrajectoryFilter(pset, iC)),
    theMaxConsecLostHitsTrajectoryFilter(new MaxConsecLostHitsTrajectoryFilter(pset, iC)),
    theMaxHitsTrajectoryFilter(new MaxHitsTrajectoryFilter(pset, iC)),
    theMaxLostHitsTrajectoryFilter(new MaxLostHitsTrajectoryFilter(pset, iC)),
    theLostHitsFractionTrajectoryFilter(new LostHitsFractionTrajectoryFilter(pset, iC)),
    theMinHitsTrajectoryFilter(new MinHitsTrajectoryFilter(pset, iC)),
    theMinPtTrajectoryFilter(new MinPtTrajectoryFilter(pset, iC)),
    theLooperTrajectoryFilter(new LooperTrajectoryFilter(pset, iC))
  {}

  void setEvent(const edm::Event& iEvent, const edm::EventSetup& iSetup) override {
    theChargeSignificanceTrajectoryFilter->setEvent(iEvent, iSetup);
    theMaxLostHitsTrajectoryFilter->setEvent(iEvent, iSetup);
    theMaxConsecLostHitsTrajectoryFilter->setEvent(iEvent, iSetup);
    theMinPtTrajectoryFilter->setEvent(iEvent, iSetup);
    theMaxHitsTrajectoryFilter->setEvent(iEvent, iSetup);
    theMinHitsTrajectoryFilter->setEvent(iEvent, iSetup);
    theLostHitsFractionTrajectoryFilter->setEvent(iEvent, iSetup);
    theLooperTrajectoryFilter->setEvent(iEvent, iSetup);
  }

  virtual bool qualityFilter( const Trajectory& traj) const {return QF<Trajectory>(traj);}
  virtual bool qualityFilter( const TempTrajectory& traj) const {return QF<TempTrajectory>(traj);}
 
  virtual bool toBeContinued( Trajectory& traj) const {return TBC<Trajectory>(traj);}
  virtual bool toBeContinued( TempTrajectory& traj) const {return TBC<TempTrajectory>(traj);}

  virtual  std::string name() const { return "CkfBaseTrajectoryFilter";}
  
protected:

  template <class T> bool QF(const T& traj) const{
    if (!theChargeSignificanceTrajectoryFilter->qualityFilter(traj)) return false;            
    if (!theMinHitsTrajectoryFilter->qualityFilter(traj)) return false;
    if (!theMinPtTrajectoryFilter->qualityFilter(traj)) return false;
    if (!theLooperTrajectoryFilter->qualityFilter(traj)) return false;
    return true;}

  template <class T> bool TBC(T& traj) const{
    if (!theMaxHitsTrajectoryFilter->toBeContinued(traj)) {
      LogDebug("CkfBaseTrajectoryFilter") << "Killing PR due to maxHits" << std::endl;
      return false;     
    }
    if (!theMaxLostHitsTrajectoryFilter->toBeContinued(traj)) {
      LogDebug("CkfBaseTrajectoryFilter") << "Killing PR due to maxLostHits" << std::endl;
      return false;
    }
    if (!theMaxConsecLostHitsTrajectoryFilter->toBeContinued(traj)) {
      LogDebug("CkfBaseTrajectoryFilter") << "Killing PR due to maxConsecLostHits" << std::endl;
      return false;
    }
    if (!theLostHitsFractionTrajectoryFilter->toBeContinued(traj)) {
      LogDebug("CkfBaseTrajectoryFilter") << "Killing PR due to LostHitsFrac" << std::endl;
      return false;
    }
    if (!theMinPtTrajectoryFilter->toBeContinued(traj)) {
      LogDebug("CkfBaseTrajectoryFilter") << "Killing PR due to minPt" << std::endl;
      return false;     
    }
    if (!theChargeSignificanceTrajectoryFilter->toBeContinued(traj)) {
      LogDebug("CkfBaseTrajectoryFilter") << "Killing PR due to chargeSignificance" << std::endl;
      return false;
    }
    if (!theLooperTrajectoryFilter->toBeContinued(traj)) {
      LogDebug("CkfBaseTrajectoryFilter") << "Killing PR due to Looper" << std::endl;
      return false;
    }
    LogDebug("CkfBaseTrajectoryFilter") << "Trajectory is good." << std::endl;
    return true;
  }

  

  std::unique_ptr<ChargeSignificanceTrajectoryFilter> theChargeSignificanceTrajectoryFilter;
  std::unique_ptr<MaxConsecLostHitsTrajectoryFilter> theMaxConsecLostHitsTrajectoryFilter;
  std::unique_ptr<MaxHitsTrajectoryFilter> theMaxHitsTrajectoryFilter;
  std::unique_ptr<MaxLostHitsTrajectoryFilter> theMaxLostHitsTrajectoryFilter;
  std::unique_ptr<LostHitsFractionTrajectoryFilter> theLostHitsFractionTrajectoryFilter;
  std::unique_ptr<MinHitsTrajectoryFilter> theMinHitsTrajectoryFilter;
  std::unique_ptr<MinPtTrajectoryFilter> theMinPtTrajectoryFilter;
  std::unique_ptr<LooperTrajectoryFilter> theLooperTrajectoryFilter;
};

#endif
