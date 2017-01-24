// Authors:  Alfredo Gurrola (Vanderbilt University)
// Andres Florez: Universidad de los Andes, Colombia.
// kaur amandeepkalsi: Panjab University, India.

#include "NtupleMaker/BSM3G_TNT_Maker/interface/GenParticleSelector.h"

GenParticleSelector::GenParticleSelector(std::string name, TTree* tree, bool debug, const pset& iConfig, edm::ConsumesCollector&& iCC):baseTree(name,tree,debug){
  if(debug) std::cout << "BSM3G TNT Maker: In the GenParticleSelector Constructor --> calling SetBranches()." << std::endl;
  genToken_                     = iCC.consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genparts"));
  _is_data                      = iConfig.getParameter<bool>("is_data");
  SetBranches();
}

GenParticleSelector::~GenParticleSelector(){
  delete tree_;
}

void GenParticleSelector::Fill(const edm::Event& iEvent){
  Clear(); 

  if(!_is_data) {

    // grabbing the gen particle handle
    Handle< reco::GenParticleCollection > _Gen_collection;
    iEvent.getByToken(genToken_, _Gen_collection);

    if(debug_) std::cout << "     GenParticleSelector: Cleared the vectors, grabbed the vertex collection handle, and looping over gen particles." << std::endl;


    // for top-pt reweight                                                                                                                                                                 
    int nLeptons = 0;
    double topPt = 0, topBarPt = 0;
    double SF_Top = 1.0, SF_antiTop = 1.0;

    // looping over gen particles
    int ng = 0;  
    for(reco::GenParticleCollection::const_iterator genparticles = _Gen_collection->begin(); genparticles !=  _Gen_collection->end(); ++genparticles) {

      if((abs(genparticles->pdgId()) == 24) && (genparticles->numberOfDaughters() == 2)) {
        for (int j = 0; j < (int)(genparticles->numberOfDaughters()); ++j) {
          const reco::Candidate * dau = genparticles->daughter(j);
          if((abs(dau->pdgId()) == 11) || (abs(dau->pdgId()) == 13)) {++nLeptons;}
        }
      }
      if((abs(genparticles->pdgId()) == 6) && (genparticles->numberOfDaughters() == 2)) {
        if((abs(genparticles->daughter(0)->pdgId()) == 24) && (abs(genparticles->daughter(1)->pdgId()) == 5)) {
          if(genparticles->pdgId() > 0) {topPt = genparticles->pt();} else {topBarPt = genparticles->pt();}
        }
      }

      Gen_pt.push_back(genparticles->pt());
      Gen_eta.push_back(genparticles->eta()); 
      Gen_phi.push_back(genparticles->phi());
      Gen_status.push_back(genparticles->status());
      Gen_pdg_id.push_back(genparticles->pdgId());
      Gen_motherpdg_id.push_back(genparticles->numberOfMothers() > 0 ? genparticles->mother(0)->pdgId() : -999999);
      Gen_energy.push_back(genparticles->energy());
      Gen_vx.push_back(genparticles->vx());
      Gen_vy.push_back(genparticles->vy());
      Gen_vz.push_back(genparticles->vz());
      Gen_charge.push_back(genparticles->charge());
      Gen_numDaught.push_back(genparticles->numberOfDaughters());
      Gen_numMother.push_back(genparticles->numberOfMothers());
    
      int idx = -1;
      for(reco::GenParticleCollection::const_iterator mit = _Gen_collection->begin();mit != _Gen_collection->end(); ++mit) {
        if(genparticles->mother() == &(*mit) ) {
          idx = std::distance(_Gen_collection->begin(),mit);
	  break;
        }
      }
      Gen_BmotherIndex.push_back(idx);

      if(debug_) {
        std::cout << "          GenParticle info: particle #" << ng << ": pdg id = " << genparticles->pdgId() << ", pt = " << genparticles->pt() << std::endl;
        if(genparticles->numberOfMothers() > 0) {
          std::cout << "                                          Mother: " << ": pdg id = " << genparticles->mother(0)->pdgId() << ", pt = " << genparticles->mother(0)->pt() << std::endl;
          std::cout << "                                          Matched mother is index: " << idx << std::endl;
        }
      }
      ng++;
    
      for (size_t j = 0; j < genparticles->numberOfMothers(); ++j) {
        const reco::Candidate* m = genparticles->mother(j);
        for (reco::GenParticleCollection::const_iterator mit = _Gen_collection->begin();  mit != _Gen_collection->end(); ++mit) {
	  if (m == &(*mit) ) { 
	    int idx = std::distance(_Gen_collection->begin(), mit);
	    Gen_BmotherIndices.push_back(idx);
	    break;
	  }
        }
      }
    
      for (size_t j = 0; j < genparticles->numberOfDaughters(); ++j) {
        const reco::Candidate* d = genparticles->daughter(j);
        for (reco::GenParticleCollection::const_iterator mit = _Gen_collection->begin();  mit != _Gen_collection->end(); ++mit) {
	  if (d == &(*mit) ) { 
	    int idx = std::distance(_Gen_collection->begin(), mit);
	    Gen_BdaughtIndices.push_back(idx);
	    break;
	  }
        }
      }
    }

    if(topPt > 700) topPt = 700;
    if(topBarPt > 700) topBarPt = 700;
    SF_Top = TMath::Exp(0.0615+((-0.0005)*topPt));                                                                                                                              
    SF_antiTop = TMath::Exp(0.0615+((-0.0005)*topBarPt));                                                                                                                       
    weighttoppt = sqrt(SF_Top*SF_antiTop);

  }

}

void GenParticleSelector::SetBranches(){
  if(debug_) std::cout << "     GenParticleSelector: Setting branches by calling AddBranch of baseTree." << std::endl;

  AddBranch(&Gen_pt               ,"Gen_pt");
  AddBranch(&Gen_eta              ,"Gen_eta");
  AddBranch(&Gen_phi              ,"Gen_phi");
  AddBranch(&Gen_status           ,"Gen_status");
  AddBranch(&Gen_pdg_id           ,"Gen_pdg_id");
  AddBranch(&Gen_motherpdg_id     ,"Gen_motherpdg_id");
  AddBranch(&Gen_energy           ,"Gen_energy");
  AddBranch(&Gen_vx               ,"Gen_vx");
  AddBranch(&Gen_vy               ,"Gen_vy");
  AddBranch(&Gen_vz               ,"Gen_vz");
  AddBranch(&Gen_charge           ,"Gen_charge");
  AddBranch(&Gen_numDaught        ,"Gen_numDaught");
  AddBranch(&Gen_numMother        ,"Gen_numMother");
  AddBranch(&Gen_BmotherIndices   ,"Gen_BmotherIndices");
  AddBranch(&Gen_BdaughtIndices   ,"Gen_BdaughtIndices");
  AddBranch(&Gen_BmotherIndex     ,"Gen_BmotherIndex");
  AddBranch(&weighttoppt          ,"weighttoppt");

  if(debug_) std::cout << "     GenParticleSelector: Finished setting branches." << std::endl;
}

void GenParticleSelector::Clear(){

  Gen_pt.clear();
  Gen_eta.clear();
  Gen_phi.clear();
  Gen_status.clear();
  Gen_pdg_id.clear();
  Gen_motherpdg_id.clear();
  Gen_energy.clear();
  Gen_vx.clear();
  Gen_vy.clear();
  Gen_vz.clear();
  Gen_charge.clear();
  Gen_numDaught.clear();
  Gen_numMother.clear();
  Gen_BmotherIndices.clear();
  Gen_BdaughtIndices.clear();
  Gen_BmotherIndex.clear();
  weighttoppt = 1;  
}
