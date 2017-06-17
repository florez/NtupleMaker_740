// Authors:  Alfredo Gurrola (Vanderbilt University)
// Andres Florez: Universidad de los Andes, Colombia.
// kaur amandeepkalsi: Panjab University, India.

#ifndef __GEN_PAR_H_

#define __GEN_PAR_H_

#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/EDConsumerBase.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include <Math/VectorUtil.h>
#include <TTree.h>
#include <TMath.h>
#include <string>
#include <vector>
#include <map>
#include <stdlib.h>
#include "baseTree.h"

using namespace std;
using namespace edm;


class GenParticleSelector : public baseTree{

 public:
  GenParticleSelector(std::string name, TTree* tree, bool debug, const edm::ParameterSet& cfg, edm::ConsumesCollector&& icc);
  ~GenParticleSelector();
  void Fill(const edm::Event& iEvent);
  void SetBranches();
  void Clear();

 private:
  GenParticleSelector(){};

  bool   _is_data;
  edm::EDGetTokenT<reco::GenParticleCollection> genToken_;

 //variables which would become branches
 vector <double> Gen_pt, Gen_eta, Gen_phi, Gen_status, Gen_pdg_id, Gen_motherpdg_id;
 vector <double> Gen_energy, Gen_vx, Gen_vy, Gen_vz;
 vector <double> Gen_charge, Gen_numDaught, Gen_numMother;
 vector <int> Gen_BmotherIndices, Gen_BdaughtIndices, Gen_BmotherIndex;
 double weighttoppt;
 double Gen_HT;
};

#endif

