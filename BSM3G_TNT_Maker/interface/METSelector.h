// Authors:  Alfredo Gurrola (Vanderbilt University)
// Andres Florez: Universidad de los Andes, Colombia.
// kaur amandeepkalsi: Panjab University, India.

#ifndef __MET_HE_H_

#define __MET_HE_H_

#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/EDConsumerBase.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include <TTree.h>
#include <string>
#include <vector>
#include <map>
#include "baseTree.h"
using namespace std;
using namespace edm;


class METSelector : public baseTree{

 public:
  METSelector(std::string name, TTree* tree, bool debug, const edm::ParameterSet& cfg, edm::ConsumesCollector&& icc);
  ~METSelector();
  void Fill(const edm::Event& iEvent);
  void SetBranches();
  void Clear();

 private:
  METSelector(){};

//  edm::InputTag metToken_;
//  edm::InputTag puppi_metToken_;
//  edm::InputTag metNoHFToken_;
  edm::EDGetTokenT<pat::METCollection> metToken_;
  edm::EDGetTokenT<pat::METCollection> puppi_metToken_;
  edm::EDGetTokenT<pat::METCollection> metNoHFToken_;
  edm::EDGetTokenT<ROOT::Math::SMatrix<double,2,2,ROOT::Math::MatRepSym<double,2> > > metCovToken_;

  //variables which would become branches

  double Met_type1PF_pt, Met_type1PF_px, Met_type1PF_py, Met_type1PF_pz, Met_type1PF_phi, Met_type1PF_sumEt, Gen_Met, Met_type1PF_shiftedPtUp, Met_type1PF_shiftedPtDown,Met_type1PF_UnclEnshiftedPtUp,Met_type1PF_UnclEnshiftedPtDown,Met_type1PF_UnclEnshiftedPhiUp , Met_type1PF_UnclEnshiftedPhiDown;

  double Met_puppi_pt, Met_puppi_px, Met_puppi_py, Met_puppi_pz, Met_puppi_phi, Met_puppi_sumEt, Met_puppi_shiftedPtUp, Met_puppi_shiftedPtDown;
  double Met_NoHF_pt, Met_NoHF_px, Met_NoHF_py, Met_NoHF_pz, Met_NoHF_phi, Met_NoHF_sumEt;
  double Met_type1PF_cov00, Met_type1PF_cov01, Met_type1PF_cov10, Met_type1PF_cov11;
  bool _is_data;
  bool _super_TNT;
};

#endif

