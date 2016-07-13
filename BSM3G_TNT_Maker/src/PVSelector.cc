#include "NtupleMaker/BSM3G_TNT_Maker/interface/PVSelector.h"

PVSelector::PVSelector(std::string name, TTree* tree, bool debug, const pset& iConfig, edm::ConsumesCollector&& iCC):baseTree(name,tree,debug) {
  if(debug) std::cout << "BSM3G TNT Maker: In the PVSelector Constructor --> getting parameters & calling SetBranches()." << std::endl;
//  _vertexInputTag               = iConfig.getParameter<edm::InputTag>("vertices");
  _vertexInputTag               = iCC.consumes<reco::VertexCollection >(iConfig.getParameter<edm::InputTag>("vertices"));
//  _beamSpot        		= iConfig.getParameter<edm::InputTag>("beamSpot");
  _beamSpot                     = iCC.consumes<reco::BeamSpot >(iConfig.getParameter<edm::InputTag>("beamSpot"));
//  _pileupInfoSrc 		= iConfig.getParameter<edm::InputTag>("pileupInfo");
  _pileupInfoSrc                = iCC.consumes<std::vector< PileupSummaryInfo > >(iConfig.getParameter<edm::InputTag>("pileupInfo"));
  _Pvtx_ndof_min   		= iConfig.getParameter<double>("Pvtx_ndof_min");
  _Pvtx_vtx_max    		= iConfig.getParameter<double>("Pvtx_vtx_max");
  _Pvtx_vtxdxy_max 		= iConfig.getParameter<double>("Pvtx_vtxdxy_max");
  _is_data         		= iConfig.getParameter<bool>("is_data");
  _super_TNT       		= iConfig.getParameter<bool>("super_TNT");
  SetBranches();
}

PVSelector::~PVSelector() {
  delete tree_;
}

void PVSelector::Fill(const edm::Event& iEvent) {

  // clear the vectors before analyzing each event
  Clear(); 

  // grab the user defined vertex collection
  edm::Handle<reco::VertexCollection> vtx;
//  iEvent.getByLabel(_vertexInputTag,vtx);
  iEvent.getByToken(_vertexInputTag,vtx);

  if(debug_) std::cout << "     PVSelector: Cleared the vectors, grabbed the vertex collection handle, and looping over PVs." << std::endl;
  
  nBestVtx = 0; 
  for(reco::VertexCollection::const_iterator vtxIt = vtx->begin(); vtxIt!= vtx->end(); ++vtxIt) {
    if((vtxIt->isValid()) && !(vtxIt->isFake())) {
      if(vtxIt->ndof() < _Pvtx_ndof_min) continue; 
      if(abs(vtxIt->z()) >= _Pvtx_vtx_max) continue;
      if(sqrt((vtxIt->x()*vtxIt->x()) + (vtxIt->y()*vtxIt->y())) >= _Pvtx_vtxdxy_max) continue; 
      nBestVtx++;
      pvertex_x.push_back(vtxIt->x());
      pvertex_y.push_back(vtxIt->y());
      pvertex_z.push_back(vtxIt->z());
      pvertex_xError.push_back(vtxIt->xError());
      pvertex_yError.push_back(vtxIt->yError());
      pvertex_zError.push_back(vtxIt->zError());
    }
  }
  
  std::vector<PileupSummaryInfo>::const_iterator PVI;
  nObservedOutOfTimePUVertices = 0;
  nObservedInTimePUVertices    = 0;
  nObservedMinus1BXPUVertices  = 0;
  nObservedPlus1BXPUVertices   = 0;
  nTruePUInteractions          = 0;
 
  if(!_is_data) {
    Handle<std::vector< PileupSummaryInfo > > PupInfo;
//    iEvent.getByLabel(_pileupInfoSrc, PupInfo); 
    iEvent.getByToken(_pileupInfoSrc, PupInfo); 

    if(debug_) std::cout << "     PVSelector: Grabbed the Pileup collection handle, looping over bunch crossings, and extracting in-time and out-of-time PU info." << std::endl;

    for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) { //loop over pileup info
      if(PVI->getBunchCrossing() == 0) { nObservedInTimePUVertices += PVI->getPU_NumInteractions(); } // number of observed in-time pileup interactions
      if(PVI->getBunchCrossing() == 0) { nTruePUInteractions = PVI->getTrueNumInteractions(); } // number of true in-time pileup interactions
      if(!_super_TNT) {
        if(PVI->getBunchCrossing()== -1) { nObservedMinus1BXPUVertices += PVI->getPU_NumInteractions(); } // number of observed pileup interactions in BX -1
        if(PVI->getBunchCrossing()==  1) { nObservedPlus1BXPUVertices += PVI->getPU_NumInteractions(); } // number of observed pileup interactions in BX +1
        if(abs(PVI->getBunchCrossing()) >= 2) { nObservedOutOfTimePUVertices += PVI->getPU_NumInteractions(); } // number of observed out-of-time pileup interactions
      }
    
      if(debug_) std::cout << "          PU info: bunch crossing, observed interactions, true interactions = " << PVI->getBunchCrossing() << ", " << PVI->getPU_NumInteractions() << ", " << PVI->getTrueNumInteractions() << std::endl;
    }
  }  

  if(debug_) std::cout << "     PVSelector: Grabbed the beamspot handle, and extracting beamspot position and width" << std::endl;

  reco::BeamSpot beamSpot;
  edm::Handle<reco::BeamSpot> beamSpotHandle;
//  iEvent.getByLabel(_beamSpot, beamSpotHandle);
  iEvent.getByToken(_beamSpot, beamSpotHandle);

  if(beamSpotHandle.isValid()) {
    beamSpot = *beamSpotHandle;
    if(debug_) std::cout << "          Beamspot info: x, y, z = " << beamSpot.x0() << ", " << beamSpot.y0() << ", " << beamSpot.z0() << std::endl;
  } else {
    std::cout << "     PVSelector: NO BEAMSPOT AVAILABLE FROM EVENTSETUP" << std::endl;
  }

  beamSpot_x0.push_back(beamSpot.x0());
  beamSpot_y0.push_back(beamSpot.y0());
  beamSpot_z0.push_back(beamSpot.z0());
  beamSpot_xWidth.push_back(beamSpot.BeamWidthX());
  beamSpot_yWidth.push_back(beamSpot.BeamWidthY());

}

void PVSelector::SetBranches(){
  if(debug_) std::cout << "     PVSelector: Setting branches by calling AddBranch of baseTree." << std::endl;

  AddBranch(&nObservedInTimePUVertices, "nObservedInTimePUVertices");
  AddBranch(&nTruePUInteractions,       "nTruePUInteractions");
  if(!_super_TNT){ 
    AddBranch(&nObservedOutOfTimePUVertices, "nObservedOutOfTimePUVertices");
    AddBranch(&nObservedPlus1BXPUVertices,   "nObservedPlus1BXPUVertices");
    AddBranch(&nObservedMinus1BXPUVertices,  "nObservedMinus1BXPUVertices");
    AddBranch(&nBestVtx        ,"bestVertices"); 
    AddBranch(&pvertex_x       ,"pvertex_x");
    AddBranch(&pvertex_y       ,"pvertex_y");
    AddBranch(&pvertex_z       ,"pvertex_z");
    AddBranch(&pvertex_xError  ,"pvertex_xError");
    AddBranch(&pvertex_yError  ,"pvertex_yError");
    AddBranch(&pvertex_zError  ,"pvertex_zError");
    AddBranch(&beamSpot_x0     ,"beamSpot_x0");
    AddBranch(&beamSpot_y0     ,"beamSpot_y0");
    AddBranch(&beamSpot_z0     ,"beamSpot_z0");
    AddBranch(&beamSpot_xWidth ,"beamSpot_xWidth");
    AddBranch(&beamSpot_yWidth ,"beamSpot_yWidth");
  }

  if(debug_) std::cout << "     PVSelector: Finished setting branches." << std::endl;
}

void PVSelector::Clear(){
  pvertex_x.clear();
  pvertex_y.clear();
  pvertex_z.clear(); 
  pvertex_xError.clear();
  pvertex_yError.clear();
  pvertex_zError.clear();
  beamSpot_x0.clear();
  beamSpot_y0.clear();
  beamSpot_z0.clear();
  beamSpot_xWidth.clear();
  beamSpot_yWidth.clear();
}
