// 
//  Authors:  Andres Florez: Universidad de los Andes, Colombia. 
//  kaur amandeepkalsi: Panjab University, India. 
//

#ifndef __TAU_MU_H_                                                                                                                                  
#define __TAU_MU_H_

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
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/TFileDirectory.h"

#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerEventWithRefs.h"
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"
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

class TauSelector : public baseTree{

public:

  TauSelector(std::string name, TTree* tree, bool debug, const edm::ParameterSet& cfg);
  ~TauSelector();
  void Fill(const edm::Event& iEvent);
  void SetBranches();
  void Clear();

private:
  TauSelector(){};
  
  vector <double> Tau_eta, Tau_phi, Tau_pt, Tau_energy, Tau_charge, Tau_chargedIsoPtSum, Tau_neutralIsoPtSum, Tau_puCorrPtSum ;
  vector<double>  Tau_leadChargedCandPt, Tau_leadChargedCandCharge, Tau_leadChargedCandEta, Tau_leadChargedCandPhi,Tau_nProngs;
  
  vector <int> Tau_decayModeFinding, Tau_decayModeFindingOldDMs, Tau_decayModeFindingNewDMs; 
  vector <int> Tau_byLooseCombinedIsolationDeltaBetaCorr,Tau_byLooseCombinedIsolationDeltaBetaCorr3Hits;
  vector <int> Tau_byMediumCombinedIsolationDeltaBetaCorr,Tau_byMediumCombinedIsolationDeltaBetaCorr3Hits; 
  vector <int> Tau_byTightCombinedIsolationDeltaBetaCorr,Tau_byTightCombinedIsolationDeltaBetaCorr3Hits; 
  vector <int> Tau_byLooseIsolationMVA3newDMwLT,Tau_byLooseIsolationMVA3newDMwoLT,Tau_byMediumIsolationMVA3newDMwLT; 
  vector <int> Tau_byMediumIsolationMVA3newDMwoLT, Tau_byLooseIsolationMva3oldDMwLT, Tau_byLooseIsolationMVA3oldDMwoLT; 
  vector <int> Tau_byMediumIsolationMva3oldDMwLT, Tau_byMediumIsolationMVA3oldDMwoLT, Tau_byTightIsolationMVA3newDMwLT;  
  vector <int> Tau_byTightIsolationMVA3newDMwoLT, Tau_byTightIsolationMva3oldDMwLT, Tau_byTightIsolationMVA3oldDMwoLT; 
  vector <int> Tau_againstMuonLoose2, Tau_againstMuonLoose3, Tau_againstMuonTight2, Tau_againstMuonTight3; 
  vector <int> Tau_againstElectronMVALooseMVA5, Tau_againstElectronMVAMediumMVA5, Tau_againstElectronMVATightMVA5; 
  vector <int> Tau_byVLooseCombinedIsolationDeltaBetaCorr; 

  // config variables
  edm::InputTag tauToken_;
  edm::InputTag _vertexInputTag;

  double _Tau_pt_min;
  double _Tau_eta_max;
  bool   _super_TNT;
};

#endif

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
