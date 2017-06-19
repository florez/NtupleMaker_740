// Authors:  Alfredo Gurrola (Vanderbilt University)
// Andres Florez: Universidad de los Andes, Colombia.
// kaur amandeepkalsi: Panjab University, India.

#ifndef __RUN_INFO_H_

#define __RUN_INFO_H_

#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include <string>
#include <vector>
#include <map>
#include "baseTree.h"
using namespace std;
using namespace edm;


class RunInfoSelector : public baseTree{

 public:
  RunInfoSelector(std::string name, TTree* tree, bool debug, const edm::ParameterSet& cfg, edm::ConsumesCollector&& icc);
  ~RunInfoSelector();
  void Fill(const edm::Event& iEvent);
  void SetBranches();
  void Clear();

 private:
  RunInfoSelector(){};

 //variables which would become branches
 int runNumber, eventNumber, lumiBlock;

};

#endif

