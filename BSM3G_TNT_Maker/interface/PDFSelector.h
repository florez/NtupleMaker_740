// Authors:  Alfredo Gurrola (Vanderbilt University)
// Andres Florez: Universidad de los Andes, Colombia.
// kaur amandeepkalsi: Panjab University, India.

#ifndef __PDF_H_

#define __PDF_H_

#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/EDConsumerBase.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include <string>
#include <vector>
#include <map>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include "baseTree.h"
using namespace std;
using namespace edm;


class Range
{
public:
  //explicit Range(int item): mLow(item),mHigh(item){}  // [item,item]
  Range(int low, int high): mLow(low), mHigh(high){}  // [low,high]

  bool operator<(const Range& rhs) const
  {
    if (mLow < rhs.mLow)
    {
	  assert(mHigh < rhs.mLow); // sanity check
      return true;
    }
    return false;
  } // operator<
  bool fits_in(const int item) const{
	  if (mLow<=item && mHigh>=item){
		  return true;
	  }
	  return false;
  }
  int low() const { return mLow; }
  int high() const { return mHigh; }

private:
  int mLow;
  int mHigh;
}; // class Range


class PDFSelector : public baseTree{

 public:
  PDFSelector(std::string name, TTree* tree, bool debug, const edm::ParameterSet& cfg, edm::ConsumesCollector&& icc, std::vector<std::string>& pdf_name);
  ~PDFSelector();
  void Fill(const edm::Event& iEvent);
  void beginRun(edm::Run const & iRun, edm::EventSetup const& iSetup);
  void SetBranches();
  void Clear();

 private:
  typedef std::map<Range, std::string> ranges_type;

  PDFSelector(){};

  bool   _is_data;

  edm::EDGetTokenT<LHEEventProduct> lheProductToken;
  edm::EDGetTokenT<LHERunInfoProduct> lheProducerToken;

  //store in a small tree with one entry
  std::vector<string>* pdf_sets;



 //variables which would become branches
 std::vector< double > pdf_weights;
};

#endif

