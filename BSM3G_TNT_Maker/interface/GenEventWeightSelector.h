// Authors:  Alfredo Gurrola (Vanderbilt University)
// Andres Florez: Universidad de los Andes, Colombia.
// kaur amandeepkalsi: Panjab University, India.

#ifndef __GEN_EVENTWEIGHT_H_

#define __GEN_EVENTWEIGHT_H_

#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/EDConsumerBase.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include <TTree.h>
#include <string>
#include <vector>
#include <map>
#include "baseTree.h"
using namespace std;
using namespace edm;


class GenEventWeightSelector : public baseTree{

 public:
  GenEventWeightSelector(std::string name, TTree* tree, bool debug, const edm::ParameterSet& cfg, edm::ConsumesCollector&& icc);
  ~GenEventWeightSelector();
  void Fill(const edm::Event& iEvent);
  void SetBranches();
  void Clear();

 private:
  GenEventWeightSelector(){};

  bool   _is_data;
  edm::EDGetTokenT<GenEventInfoProduct> generatorToken_;

 //variables which would become branches
 double weightevt;
};

#endif

