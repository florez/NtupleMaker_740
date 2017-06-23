// -*- C++ -*-
//
// Package:    MiniAOD/BSM3G_TNT_Maker
// Class:      BSM3G_TNT_Maker
//
//
// Authors:  Alfredo Gurrola (Vanderbilt University)
// Andres Florez: Universidad de los Andes, Colombia.
// kaur amandeepkalsi: Panjab University, India.
//

#ifndef  TREE_MAKER_H
#define  TREE_MAKER_H


// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "NtupleMaker/BSM3G_TNT_Maker/interface/TriggerSelector.h"
#include "NtupleMaker/BSM3G_TNT_Maker/interface/MuonSelector.h"
#include "NtupleMaker/BSM3G_TNT_Maker/interface/TauSelector.h"
#include "NtupleMaker/BSM3G_TNT_Maker/interface/JetSelector.h"
#include "NtupleMaker/BSM3G_TNT_Maker/interface/PVSelector.h"
#include "NtupleMaker/BSM3G_TNT_Maker/interface/GenParticleSelector.h"
#include "NtupleMaker/BSM3G_TNT_Maker/interface/GenEventWeightSelector.h"
#include "NtupleMaker/BSM3G_TNT_Maker/interface/RunInfoSelector.h"
#include "NtupleMaker/BSM3G_TNT_Maker/interface/METSelector.h"
#include "NtupleMaker/BSM3G_TNT_Maker/interface/ElectronPatSelector.h"
#include "NtupleMaker/BSM3G_TNT_Maker/interface/PhotonSelector.h"
#include "NtupleMaker/BSM3G_TNT_Maker/interface/PDFSelector.h"
#include "baseTree.h"
#include <TBranch.h>
#include <TTree.h>
#include <TFile.h>
#include <vector>

#ifdef __MAKECINT__
#pragma link C++ class std::vector<std::vector<int> >+;
#pragma link C++ class std::vector<std::vector<std::string> >+;
#pragma link C++ class std::vector<std::vector<TString> >+;
#pragma link C++ class std::vector<std::vector<float> >+;
#pragma link C++ class std::vector<std::vector<double> >+;
#pragma link C++ class std::vector<std::vector<bool> >+;
#pragma extra_include "std::vector";
#endif



//
// class declaration
//

class BSM3G_TNT_Maker : public edm::EDAnalyzer {
public:
  explicit BSM3G_TNT_Maker(const edm::ParameterSet&);
  ~BSM3G_TNT_Maker();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  virtual void beginJob() override;
  virtual void beginRun(edm::Run const &, edm::EventSetup const&);
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;

  // ----------member data ---------------------------

  TFile* file;
  TTree* tree_;
  TTree* storetree_;

  std::vector<string> pdf_name;
  std::vector<string> triggernames;

  const size_t MaxN;
  bool debug_;

  bool _filltriggerinfo;
  bool _fillmuoninfo;
  bool _fillelectronpatinfo;
  bool _filltauinfo;
  bool _filljetinfo;
  bool _fillgeninfo;
  bool _fillruninfo;
  bool _fillgenweightinfo;
  bool _fillPVinfo;
  bool _fillMETinfo;
  bool _fillphotoninfo;
  bool _fillPDFinfo;

  TriggerSelector        *trselector;
  MuonSelector           *muselector;
  TauSelector            *tauselector;
  JetSelector            *jetselector;
  GenParticleSelector    *genselector;
  GenEventWeightSelector *genweightselector;
  RunInfoSelector        *runinfoselector;
  PVSelector             *pvselector;
  METSelector            *metselector;
  ElectronPatSelector    *elpatselector;
  PhotonSelector         *photonselector;
  PDFSelector            *pdfselector;
};

#endif


