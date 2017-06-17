#ifndef __PV_HE_H_

#define __PV_HE_H_

#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/EDConsumerBase.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Common/interface/Handle.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include <TTree.h>
#include <string>
#include <vector>
#include <map>
#include "baseTree.h"
using namespace std;
using namespace edm;


class PVSelector : public baseTree{

 public:
  PVSelector(std::string name, TTree* tree, bool debug, const edm::ParameterSet& cfg, edm::ConsumesCollector&& icc);
  ~PVSelector();
  void Fill(const edm::Event& iEvent);
  void SetBranches();
  void Clear();

 private:
  PVSelector(){};
//  edm::InputTag _pileupInfoSrc;
//  edm::InputTag _beamSpot;
//  edm::InputTag _vertexInputTag;
  edm::EDGetTokenT<reco::VertexCollection> _vertexInputTag;
  edm::EDGetTokenT<reco::BeamSpot> _beamSpot;
  edm::EDGetTokenT<std::vector< PileupSummaryInfo > > _pileupInfoSrc;

  // variables which will become branches
  int nBestVtx;
  int nObservedOutOfTimePUVertices, nObservedInTimePUVertices, nObservedMinus1BXPUVertices, nObservedPlus1BXPUVertices;
  vector<float> pvertex_x, pvertex_y, pvertex_z;
  vector<float> pvertex_xError, pvertex_yError, pvertex_zError;
  vector<float> beamSpot_x0, beamSpot_y0, beamSpot_z0, beamSpot_xWidth, beamSpot_yWidth;
  float nTruePUInteractions;

  // Primary vertex cuts
  double _Pvtx_ndof_min;
  double _Pvtx_vtx_max;
  double _Pvtx_vtxdxy_max;
  bool   _is_data;
  bool   _super_TNT;
};

#endif

