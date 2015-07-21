#include "NtupleMaker/BSM3G_TNT_Maker/interface/PVSelector.h"

PVSelector::PVSelector(std::string name, TTree* tree, bool debug, const pset& iConfig):baseTree(name,tree,debug){
  if(debug) std::cout<<"in PVSelector constructor"<<std::endl;
  if(debug) std::cout<<"in pileup constructor: calling SetBrances()"<<std::endl;

  _Pvtx_ndof_min   = iConfig.getParameter<double>("Pvtx_ndof_min");
  _Pvtx_vtx_max    = iConfig.getParameter<double>("Pvtx_vtx_max");
  _Pvtx_vtxdxy_max = iConfig.getParameter<double>("Pvtx_vtxdxy_max");
  _is_data         = iConfig.getParameter<bool>("is_data");
  _super_TNT       = iConfig.getParameter<bool>("super_TNT");
  _beamSpot        = iConfig.getParameter<edm::InputTag>("beamSpot");
  SetBranches();
}

PVSelector::~PVSelector(){
  delete tree_;
}

void PVSelector::Fill(const edm::Event& iEvent){
   
  Clear(); 

  if(debug_)    std::cout<<"getting pileup info"<<std::endl;
  
  edm::Handle<reco::VertexCollection> vtx;
  iEvent.getByLabel("offlineSlimmedPrimaryVertices",vtx);
  
  nBestVtx = 0; 
  for( reco::VertexCollection::const_iterator vtxIt = vtx->begin(); vtxIt!= vtx->end(); ++vtxIt){
    
    if((vtxIt->isValid()) && !(vtxIt->isFake())) {
      
      if(vtxIt->ndof() < _Pvtx_ndof_min) continue; 
      if(abs(vtxIt->z()) >= _Pvtx_vtx_max) continue;
      double vtxdxy = sqrt((vtxIt->x()*vtxIt->x()) + (vtxIt->y()*vtxIt->y() ));
      if(vtxdxy >=  _Pvtx_vtxdxy_max) continue; 
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
  ootnpuVertices   = 0;
  npuVertices      = 0;
  npuVerticesm1    = 0;
  npuVerticesp1    = 0;
  trueInteractions = 0.;
 
  if (!_is_data){
    Handle<std::vector< PileupSummaryInfo > >  PupInfo;
    iEvent.getByLabel(std::string("addPileupInfo"), PupInfo); 
    for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
      if(PVI->getBunchCrossing()==  0){ npuVertices   += PVI->getPU_NumInteractions();  }
      if(!_super_TNT){
        if(PVI->getBunchCrossing()== -1){ npuVerticesm1 += PVI->getPU_NumInteractions(); }
        if(PVI->getBunchCrossing()==  1){ npuVerticesp1 += PVI->getPU_NumInteractions(); }
        if(abs(PVI->getBunchCrossing()) >= 2){
          ootnpuVertices += PVI->getPU_NumInteractions();
        }
      }
      if(PVI->getBunchCrossing() == 0)trueInteractions = PVI->getTrueNumInteractions();
    
      if(debug_)std::cout << " Pileup Information: bunchXing, nvtx,true nvtx: " << PVI->getBunchCrossing() << " " << PVI->getPU_NumInteractions()<< " "<< PVI->getTrueNumInteractions()<< std::endl;
    
    }//loop over pileup info
  }  
  if(debug_)    std::cout<<"got pileup info"<<std::endl;


  reco::BeamSpot beamSpot;
  edm::Handle<reco::BeamSpot> beamSpotHandle;
  iEvent.getByLabel(_beamSpot, beamSpotHandle);

  if ( beamSpotHandle.isValid() )
   {
     beamSpot = *beamSpotHandle;

   } else {
    edm::LogInfo("MyAnalyzer")
      << "No beam spot available from EventSetup \n";
   }

  beamSpot_x0.push_back(beamSpot.x0());
  beamSpot_y0.push_back(beamSpot.y0());
  beamSpot_z0.push_back(beamSpot.z0());

}

void PVSelector::SetBranches(){
  if(debug_)    std::cout<<"setting branches: calling AddBranch of baseTree"<<std::endl;
  AddBranch(&npuVertices     ,"npuVertices");
  AddBranch(&trueInteractions,"trueInteractions");

  if(!_super_TNT){ 
    AddBranch(&ootnpuVertices  ,"ootnpuVertices");
    AddBranch(&npuVerticesp1   ,"npuVerticesp1");
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
  }
  if(debug_)    std::cout<<"set branches"<<std::endl;
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
}
