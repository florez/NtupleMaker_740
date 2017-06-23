//
// Authors:  Andres Florez: Universidad de los Andes, Colombia.
// kaur amandeepkalsi: Panjab University, India.
//
#ifndef __PHOTON_H_
#define __PHOTON_H_


#include <vector>
#include <string>
#include <map>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/EDConsumerBase.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "baseTree.h"

using namespace std;
using namespace pat;
using namespace edm;
//
// class declaration
//

class PhotonSelector : public  baseTree{

public:

  PhotonSelector(std::string name, TTree* tree, bool debug, const edm::ParameterSet& cfg, edm::ConsumesCollector&& icc);
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
//  edm::InputTag _PhotonToken;
  edm::EDGetTokenT<edm::View<pat::Photon> > _PhotonToken;
  double _Photon_pt_min;
  double _Photon_eta_max;
};

#endif
