// 
// Authors:  Andres Florez: Universidad de los Andes, Colombia. 
// kaur amandeepkalsi: Panjab University, India. 
//
#ifndef __PHOTON_H_                                                                                                                                  
#define __PHOTON_H_


#include <memory>
#include <iostream>
#include <cmath>
#include <vector>
#include <TBranch.h>
#include <TTree.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <string>
#include <map>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <TRandom3.h>
#include <TBranch.h>                                                                    
#include <TClonesArray.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "Math/VectorUtil.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/TFileDirectory.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "baseTree.h"

using namespace std;
using namespace pat;
using namespace edm;
//
// class declaration
//

class PhotonSelector : public  baseTree{

public:

  PhotonSelector(std::string name, TTree* tree, bool debug, const edm::ParameterSet& cfg);
  ~PhotonSelector();
  void Fill(const edm::Event& iEvent);
  void SetBranches();
  void Clear();
 private:
  PhotonSelector(){};
  
  // ----------member data ---------------------------
  
  vector<float> Photon_pt , Photon_eta, Photon_phi, Photon_energy;
  vector<float> Photon_et, Photon_HoverE, Photon_phoR9, Photon_SigmaIEtaIEta;
  vector<float> Photon_SigmaIPhiIPhi, Photon_PFChIso, Photon_PFPhoIso, Photon_PFNeuIso;
  vector<int>   Photon_EleVeto, Photon_hasPixelSeed;
 
  // confit variables
  edm::InputTag _PhotonToken;
  double _Photon_pt_min;
  double _Photon_eta_max;
};

#endif
