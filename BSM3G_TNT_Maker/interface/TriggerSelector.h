//
//  Authors:  Andres Florez: Universidad de los Andes, Colombia.
//  kaur amandeepkalsi: Panjab University, India.
//

#ifndef __TRIGGER_H_

#define __TRIGGER_H_

#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/EDConsumerBase.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "baseTree.h"

using namespace std;
using namespace edm;

class TriggerSelector : public baseTree{

 public:
  TriggerSelector(std::string name, TTree* tree, bool debug, const edm::ParameterSet& cfg, edm::ConsumesCollector&& icc, std::vector<string>& _triggernames);
  ~TriggerSelector();
  virtual void startTrigger (edm::EventSetup const& , edm::Run const &);
  void Fill(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  void SetBranches();
  void Clear();

 private:
  TriggerSelector(){};
  std::vector<string>* triggernames;
  vector <int> Trigger_decision;
  //vector <string> Trigger_names;
  vector <int> Trigger_prescale;
  edm::EDGetTokenT<edm::TriggerResults> triggerResultsTag_;
  edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescales_;
};

#endif

