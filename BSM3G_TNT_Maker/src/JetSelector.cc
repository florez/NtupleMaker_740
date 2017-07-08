// Authors:  Alfredo Gurrola (Vanderbilt University)
// Andres Florez: Universidad de los Andes, Colombia.
// kaur amandeepkalsi: Panjab University, India.

#include "NtupleMaker/BSM3G_TNT_Maker/interface/JetSelector.h"

JetSelector::JetSelector(std::string name, TTree* tree, bool debug, const pset& iConfig, edm::ConsumesCollector&& iCC):baseTree(name,tree,debug){
  if(debug) std::cout << "BSM3G TNT Maker: In the JetSelector Constructor --> getting parameters & calling SetBranches()." << std::endl;
  jetToken_       = iCC.consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jets"));
  puppi_jetToken_ = iCC.consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jetsPUPPI"));
  toptag_jetToken_ = iCC.consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jetsTOPTAG"));
  _vertexInputTag = iCC.consumes<reco::VertexCollection >(iConfig.getParameter<edm::InputTag>("vertices"));
  _Jet_pt_min     = iConfig.getParameter<double>("Jet_pt_min");
  _super_TNT      = iConfig.getParameter<bool>("super_TNT");
  bool _is_data                      = iConfig.getParameter<bool>("is_data");


  //stuff for the fatjet correction:
  std::vector<std::string> jecAK8PayloadNames_;
  jecAK8PayloadNames_.push_back(iConfig.getParameter<edm::FileInPath>("JetAk8PDF_L2_correction").fullPath());
  jecAK8PayloadNames_.push_back(iConfig.getParameter<edm::FileInPath>("JetAk8PDF_L3_correction").fullPath());
  if(_is_data)
    jecAK8PayloadNames_.push_back(iConfig.getParameter<edm::FileInPath>("JetAk8PDF_L2L3_residual_correction").fullPath());

  std::vector<JetCorrectorParameters> vPar;
  for ( std::vector<std::string>::const_iterator payloadBegin = jecAK8PayloadNames_.begin(), payloadEnd = jecAK8PayloadNames_.end(), ipayload = payloadBegin; ipayload != payloadEnd; ++ipayload ) {
    JetCorrectorParameters pars(*ipayload);
    vPar.push_back(pars);
  }
  // Make the FactorizedJetCorrector
  jecAK8_ = boost::shared_ptr<FactorizedJetCorrector> ( new FactorizedJetCorrector(vPar) );

  puppyCorrFile  = TFile::Open(iConfig.getParameter<edm::FileInPath>("JetAk8PDF_puppymc").fullPath().c_str(),"READ");
  puppisd_corrGEN      = (TF1*)puppyCorrFile->Get("puppiJECcorr_gen");
  puppisd_corrRECO_cen = (TF1*)puppyCorrFile->Get("puppiJECcorr_reco_0eta1v3");
  puppisd_corrRECO_for = (TF1*)puppyCorrFile->Get("puppiJECcorr_reco_1v3eta2v5");

  //we need also rho and the number of PV
  rhos_                    = iConfig.getParameter<std::vector< edm::InputTag >>("rhos");
  for (std::vector< edm::InputTag >::const_iterator rhos_label = rhos_.begin(); rhos_label != rhos_.end(); ++rhos_label) {
        iCC.consumes<double>( *rhos_label);
  }
  _vertexInputTag               = iCC.consumes<reco::VertexCollection >(iConfig.getParameter<edm::InputTag>("vertices"));




  SetBranches();
}

JetSelector::~JetSelector(){
  delete tree_;
}

void JetSelector::Fill(const edm::Event& iEvent){
  Clear();

  // grab handle to the jet collections
  edm::Handle<pat::JetCollection> jets;
  iEvent.getByToken(jetToken_, jets);
  edm::Handle<pat::JetCollection> puppijets;
  iEvent.getByToken(puppi_jetToken_, puppijets);
  edm::Handle<pat::JetCollection> toptagjets;
  iEvent.getByToken(toptag_jetToken_, toptagjets);

  //get all rhos that are stored:
  for (std::vector< edm::InputTag >::const_iterator rho_label = rhos_.begin(); rho_label != rhos_.end(); ++rho_label) {
      edm::Handle< double > _rho;
      iEvent.getByLabel(*rho_label, _rho);
      rho.push_back(*_rho);
  }
  edm::Handle<reco::VertexCollection> vtx;
  iEvent.getByToken(_vertexInputTag,vtx);

  if(debug_) std::cout << "     JetSelector: Cleared the vectors, grabbed the jet collection handle, and looping over jets." << std::endl;

  // loop over "standard" ak4 CHS PF jets
  for (const pat::Jet &j : *jets) {
    if (j.pt() < _Jet_pt_min) continue; // only keep jets passing some user defined pt cut (defined in miniAOD.py)

    // fill root tree with "necessary" information:  kinematics, ID discriminators
    Jet_pt.push_back(j.pt());
    Jet_eta.push_back(j.eta());
    Jet_phi.push_back(j.phi());
    Jet_energy.push_back(j.energy());
    //Jet_bDiscriminator.push_back(j.bDiscriminator("combinedSecondaryVertexBJetTags"));
    Jet_bDiscriminator_pfCMVAV2.push_back(j.bDiscriminator("pfCombinedMVAV2BJetTags"));
    Jet_bDiscriminator_pfCISVV2.push_back(j.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
    //Jet_vtxMass.push_back(j.userFloat("vtxMass"));
    //Jet_decayLength.push_back(j.userFloat("vtx3DVal"));
    //Jet_decayLengthSignificance.push_back(j.userFloat("vtx3DSig"));
    //Jet_pileupId.push_back(j.userFloat("pileupJetId:fullDiscriminant"));
    Jet_partonFlavour.push_back(j.partonFlavour());
    Jet_neutralHadEnergyFraction.push_back(j.neutralHadronEnergyFraction());
    Jet_neutralEmEmEnergyFraction.push_back(j.neutralEmEnergyFraction());
    Jet_chargedHadronEnergyFraction.push_back(j.chargedHadronEnergyFraction());
    Jet_chargedEmEnergyFraction.push_back(j.chargedEmEnergyFraction());
    Jet_muonEnergyFraction.push_back(j.muonEnergyFraction());
    Jet_numberOfConstituents.push_back(j.chargedMultiplicity() + j.neutralMultiplicity());
    Jet_chargedMultiplicity.push_back(j.chargedMultiplicity());

    // store additional information such as the "raw" uncorrected pt
    if(!_super_TNT){
      Jet_mass.push_back(j.mass());
      Jet_electronEnergy.push_back(j.electronEnergy());
      Jet_photonEnergy.push_back(j.photonEnergy());
      UncorrJet_pt.push_back(j.correctedJet("Uncorrected").pt());
    }
  }

  // loop over ak4 PF jets reconstructed with the PUPPI pileup mitigation algorithms
  for (const pat::Jet &j : *puppijets) {
    if (j.pt() < _Jet_pt_min) continue; // only keep jets passing some user defined pt cut (defined in miniAOD.py)

    // fill root tree with "necessary" information:  kinematics, ID discriminators
    Jet_puppi_pt.push_back(j.pt());
    Jet_puppi_eta.push_back(j.eta());
    Jet_puppi_phi.push_back(j.phi());
    Jet_puppi_energy.push_back(j.energy());
    //Jet_puppi_bDiscriminator.push_back(j.bDiscriminator("combinedSecondaryVertexBJetTags"));
    Jet_puppi_bDiscriminator_pfCMVAV2.push_back(j.bDiscriminator("pfCombinedMVAV2BJetTags"));
    Jet_puppi_bDiscriminator_pfCISVV2.push_back(j.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
    //Jet_puppi_vtxMass.push_back(j.userFloat("vtxMass"));
    //Jet_puppi_decayLength.push_back(j.userFloat("vtx3DVal"));
    //Jet_puppi_decayLengthSignificance.push_back(j.userFloat("vtx3DSig"));
    //Jet_puppi_pileupId.push_back(j.userFloat("pileupJetId:fullDiscriminant"));
    Jet_puppi_partonFlavour.push_back(j.partonFlavour());
    Jet_puppi_neutralHadEnergyFraction.push_back(j.neutralHadronEnergyFraction());
    Jet_puppi_neutralEmEmEnergyFraction.push_back(j.neutralEmEnergyFraction());
    Jet_puppi_chargedHadronEnergyFraction.push_back(j.chargedHadronEnergyFraction());
    Jet_puppi_chargedEmEnergyFraction.push_back(j.chargedEmEnergyFraction());
    Jet_puppi_muonEnergyFraction.push_back(j.muonEnergyFraction());
    Jet_puppi_numberOfConstituents.push_back(j.chargedMultiplicity() + j.neutralMultiplicity());
    Jet_puppi_chargedMultiplicity.push_back(j.chargedMultiplicity());

    // store additional information such as the "raw" uncorrected pt
    if(!_super_TNT){
      Jet_puppi_mass.push_back(j.mass());
      Jet_puppi_electronEnergy.push_back(j.electronEnergy());
      Jet_puppi_photonEnergy.push_back(j.photonEnergy());
      UncorrJet_puppi_pt.push_back(j.correctedJet("Uncorrected").pt());
    }
  }

  // loop over ak8 PF jets reconstructed
  for (const pat::Jet &j : *toptagjets) {
    if (j.pt() < _Jet_pt_min) continue; // only keep jets passing some user defined pt cut (defined in miniAOD.py)

    // fill root tree with "necessary" information:  kinematics, ID discriminators
    Jet_toptag_pt.push_back(j.pt());
    Jet_toptag_eta.push_back(j.eta());
    Jet_toptag_phi.push_back(j.phi());
    Jet_toptag_energy.push_back(j.energy());
    Jet_toptag_SoftDropMass.push_back(j.userFloat("ak8PFJetsCHSSoftDropMass"));
    Jet_toptag_PrunedMass.push_back(j.userFloat("ak8PFJetsCHSPrunedMass"));
    Jet_toptag_tau1.push_back(j.userFloat("NjettinessAK8:tau1"));
    Jet_toptag_tau2.push_back(j.userFloat("NjettinessAK8:tau2"));
    Jet_toptag_tau3.push_back(j.userFloat("NjettinessAK8:tau3"));
    Jet_toptag_UncorrJet.push_back(j.correctedJet("Uncorrected").pt());
    //get the correction for the Ak8PFchs jets
    jecAK8_->setJetPt ( j.pt() );
    jecAK8_->setJetEta( j.eta() );
    jecAK8_->setJetE  ( j.energy() );
    //jecAK8_->setJetPt ( j.correctedJet("Uncorrected").pt() );
    //jecAK8_->setJetEta( j.correctedJet("Uncorrected").eta() );
    //jecAK8_->setJetE  ( j.correctedJet("Uncorrected").energy() );
    jecAK8_->setJetA  ( j.jetArea() );
    jecAK8_->setRho   ( rho[0] );
    jecAK8_->setNPV   ( vtx->size() );
    double corr = jecAK8_->getCorrection();
    //ok we need this correction this is not the same:
    //cout<<corr<<"   "<<j.pt()/j.correctedJet("Uncorrected").pt()<<endl;
    Jet_toptag_correction.push_back(corr);


    Jet_toptag_puppi_pt.push_back(j.userFloat("ak8PFJetsPuppiValueMap:pt"));
    Jet_toptag_puppi_mass.push_back(j.userFloat("ak8PFJetsPuppiValueMap:mass"));
    Jet_toptag_puppi_eta.push_back(j.userFloat("ak8PFJetsPuppiValueMap:eta"));
    Jet_toptag_puppi_phi.push_back(j.userFloat("ak8PFJetsPuppiValueMap:phi"));
    Jet_toptag_puppi_tau1.push_back(j.userFloat("ak8PFJetsPuppiValueMap:NjettinessAK8PuppiTau1"));
    Jet_toptag_puppi_tau2.push_back(j.userFloat("ak8PFJetsPuppiValueMap:NjettinessAK8PuppiTau2"));
    Jet_toptag_puppi_tau3.push_back(j.userFloat("ak8PFJetsPuppiValueMap:NjettinessAK8PuppiTau3"));


    TLorentzVector puppi_softdrop, puppi_softdrop_subjet;
    for ( auto const & it : j.subjets("SoftDropPuppi") ) {
      puppi_softdrop_subjet.SetPtEtaPhiM(it->correctedP4(0).pt(),it->correctedP4(0).eta(),it->correctedP4(0).phi(),it->correctedP4(0).mass());
      puppi_softdrop+=puppi_softdrop_subjet;
    }

    double puppiCorr = GetPuppyak8Corr( Jet_toptag_puppi_pt.back() , Jet_toptag_puppi_eta.back() );
    //Jet_toptag_puppy_correction.push_back(puppiCorr);
    Jet_toptag_puppi_mass_corr.push_back(puppi_softdrop.M()*puppiCorr);

  }

}

void JetSelector::SetBranches(){
  if(debug_) std::cout << "     JetSelector: Setting branches by calling AddBranch of baseTree." << std::endl;

  AddBranch(&Jet_pt,                            "Jet_pt");
  AddBranch(&Jet_eta,                           "Jet_eta");
  AddBranch(&Jet_phi,                           "Jet_phi");
  AddBranch(&Jet_energy,                        "Jet_energy");
  //AddBranch(&Jet_bDiscriminator,          "Jet_bDiscriminator");
  AddBranch(&Jet_bDiscriminator_pfCMVAV2,       "Jet_bDiscriminator_pfCMVAV2");
  AddBranch(&Jet_bDiscriminator_pfCISVV2,       "Jet_bDiscriminator_pfCISVV2");
  //AddBranch(&Jet_vtxMass,               "Jet_vtxMass");
  //AddBranch(&Jet_decayLength,             "Jet_decayLength");
  //AddBranch(&Jet_decayLengthSignificance, "Jet_decayLengthSignificance");
  //AddBranch(&Jet_pileupId,                "Jet_pileupId");
  AddBranch(&Jet_partonFlavour,                 "Jet_partonFlavour");
  AddBranch(&Jet_neutralHadEnergyFraction,      "Jet_neutralHadEnergyFraction");
  AddBranch(&Jet_neutralEmEmEnergyFraction,     "Jet_neutralEmEmEnergyFraction");
  AddBranch(&Jet_chargedHadronEnergyFraction,   "Jet_chargedHadronEnergyFraction");
  AddBranch(&Jet_chargedEmEnergyFraction,       "Jet_chargedEmEnergyFraction");
  AddBranch(&Jet_muonEnergyFraction,            "Jet_muonEnergyFraction");
  AddBranch(&Jet_numberOfConstituents,          "Jet_numberOfConstituents");
  AddBranch(&Jet_chargedMultiplicity,           "Jet_chargedMultiplicity");
  AddBranch(&Jet_puppi_pt,                      "Jet_puppi_pt");
  AddBranch(&Jet_puppi_eta,                     "Jet_puppi_eta");
  AddBranch(&Jet_puppi_phi,                     "Jet_puppi_phi");
  AddBranch(&Jet_puppi_energy,                  "Jet_puppi_energy");
  //AddBranch(&Jet_puppi_bDiscriminator,        "Jet_puppi_bDiscriminator");
  AddBranch(&Jet_puppi_bDiscriminator_pfCMVAV2,   "Jet_puppi_bDiscriminator_pfCMVAV2");
  AddBranch(&Jet_puppi_bDiscriminator_pfCISVV2, "Jet_puppi_bDiscriminator_pfCISVV2");
  //AddBranch(&Jet_puppi_vtxMass,           "Jet_puppi_vtxMass");
  //AddBranch(&Jet_puppi_decayLength,             "Jet_puppi_decayLength");
  //AddBranch(&Jet_puppi_decayLengthSignificance, "Jet_puppi_decayLengthSignificance");
  //AddBranch(&Jet_puppi_pileupId,              "Jet_puppi_pileupId");
  AddBranch(&Jet_puppi_partonFlavour,         "Jet_puppi_partonFlavour");
  AddBranch(&Jet_puppi_neutralHadEnergyFraction,"Jet_puppi_neutralHadEnergyFraction");
  AddBranch(&Jet_puppi_neutralEmEmEnergyFraction,"Jet_puppi_neutralEmEmEnergyFraction");
  AddBranch(&Jet_puppi_chargedHadronEnergyFraction,"Jet_puppi_chargedHadronEnergyFraction");
  AddBranch(&Jet_puppi_chargedEmEnergyFraction, "Jet_puppi_chargedEmEnergyFraction");
  AddBranch(&Jet_puppi_muonEnergyFraction,      "Jet_puppi_muonEnergyFraction");
  AddBranch(&Jet_puppi_numberOfConstituents,    "Jet_puppi_numberOfConstituents");
  AddBranch(&Jet_puppi_chargedMultiplicity,     "Jet_puppi_chargedMultiplicity");
  AddBranch(&Jet_toptag_pt,                     "Jet_toptag_pt");
  AddBranch(&Jet_toptag_eta,                    "Jet_toptag_eta");
  AddBranch(&Jet_toptag_phi,                    "Jet_toptag_phi");
  AddBranch(&Jet_toptag_energy,                 "Jet_toptag_energy");
  AddBranch(&Jet_toptag_SoftDropMass,           "Jet_toptag_SoftDropMass");
  AddBranch(&Jet_toptag_PrunedMass,             "Jet_toptag_PrunedMass");
  AddBranch(&Jet_toptag_tau1,                   "Jet_toptag_tau1");
  AddBranch(&Jet_toptag_tau2,                   "Jet_toptag_tau2");
  AddBranch(&Jet_toptag_tau3,                   "Jet_toptag_tau3");
  AddBranch(&Jet_toptag_correction,             "Jet_toptag_correction");
  AddBranch(&Jet_toptag_UncorrJet,              "Jet_toptag_UncorrJet");
  AddBranch(&Jet_toptag_puppi_pt,               "Jet_toptag_puppi_pt");
  AddBranch(&Jet_toptag_puppi_eta,              "Jet_toptag_puppi_eta");
  AddBranch(&Jet_toptag_puppi_phi,              "Jet_toptag_puppi_phi");
  AddBranch(&Jet_toptag_puppi_mass,             "Jet_toptag_puppi_mass");
  AddBranch(&Jet_toptag_puppi_tau1,             "Jet_toptag_puppi_tau1");
  AddBranch(&Jet_toptag_puppi_tau2,             "Jet_toptag_puppi_tau2");
  AddBranch(&Jet_toptag_puppi_tau3,             "Jet_toptag_puppi_tau3");
  AddBranch(&Jet_toptag_puppi_mass_corr,             "Jet_toptag_puppi_mass_corr");
  //AddBranch(&Jet_toptag_puppy_correction,       "Jet_toptag_puppy_correction");

  if(!_super_TNT){
    AddBranch(&Jet_mass,                        "Jet_mass");
    AddBranch(&Jet_electronEnergy,              "Jet_electronEnergy");
    AddBranch(&Jet_photonEnergy,                "Jet_photonEnergy");
    AddBranch(&UncorrJet_pt,                    "UncorrJet_pt");
    AddBranch(&Jet_puppi_mass,                  "Jet_puppi_mass");
    AddBranch(&Jet_puppi_electronEnergy,        "Jet_puppi_electronEnergy");
    AddBranch(&Jet_puppi_photonEnergy,          "Jet_puppi_photonEnergy");
    AddBranch(&UncorrJet_puppi_pt,              "UncorrJet_puppi_pt");
  }

  AddBranch(&rho            ,"rho");

  if(debug_) std::cout << "     JetSelector: Finished setting branches." << std::endl;
}

void JetSelector::Clear(){

  Jet_pt.clear();
  Jet_eta.clear();
  Jet_phi.clear();
  Jet_energy.clear();
  //Jet_bDiscriminator.clear();
  Jet_bDiscriminator_pfCMVAV2.clear();
  Jet_bDiscriminator_pfCISVV2.clear();
  //Jet_vtxMass.clear();
  //Jet_decayLength.clear();
  //Jet_decayLengthSignificance.clear();
  //Jet_pileupId.clear();
  Jet_partonFlavour.clear();
  Jet_mass.clear();
  Jet_neutralHadEnergyFraction.clear();
  Jet_neutralEmEmEnergyFraction.clear();
  Jet_chargedHadronEnergyFraction.clear();
  Jet_chargedEmEnergyFraction.clear();
  Jet_muonEnergyFraction.clear();
  Jet_numberOfConstituents.clear();
  Jet_chargedMultiplicity.clear();
  Jet_electronEnergy.clear();
  Jet_photonEnergy.clear();
  UncorrJet_pt.clear();
  Jet_puppi_pt.clear();
  Jet_puppi_eta.clear();
  Jet_puppi_phi.clear();
  Jet_puppi_energy.clear();
  //Jet_puppi_bDiscriminator.clear();
  Jet_puppi_bDiscriminator_pfCMVAV2.clear();
  Jet_puppi_bDiscriminator_pfCISVV2.clear();
  //Jet_puppi_vtxMass.clear();
  //Jet_puppi_decayLength.clear();
  //Jet_puppi_decayLengthSignificance.clear();
  //Jet_puppi_pileupId.clear();
  Jet_puppi_partonFlavour.clear();
  Jet_puppi_mass.clear();
  Jet_puppi_neutralHadEnergyFraction.clear();
  Jet_puppi_neutralEmEmEnergyFraction.clear();
  Jet_puppi_chargedHadronEnergyFraction.clear();
  Jet_puppi_chargedEmEnergyFraction.clear();
  Jet_puppi_muonEnergyFraction.clear();
  Jet_puppi_numberOfConstituents.clear();
  Jet_puppi_chargedMultiplicity.clear();
  Jet_puppi_electronEnergy.clear();
  Jet_puppi_photonEnergy.clear();
  UncorrJet_puppi_pt.clear();
  Jet_toptag_pt.clear();
  Jet_toptag_eta.clear();
  Jet_toptag_phi.clear();
  Jet_toptag_energy.clear();
  Jet_toptag_SoftDropMass.clear();
  Jet_toptag_PrunedMass.clear();
  Jet_toptag_tau1.clear();
  Jet_toptag_tau2.clear();
  Jet_toptag_tau3.clear();
  Jet_toptag_correction.clear();
  Jet_toptag_puppi_pt.clear();
  Jet_toptag_puppi_eta.clear();
  Jet_toptag_puppi_phi.clear();
  Jet_toptag_puppi_mass.clear();
  Jet_toptag_puppi_tau1.clear();
  Jet_toptag_puppi_tau2.clear();
  Jet_toptag_puppi_tau3.clear();
  Jet_toptag_puppi_mass_corr.clear();
  Jet_toptag_UncorrJet.clear();
  //Jet_toptag_puppy_correction.clear();

  rho.clear();

}

double JetSelector::GetPuppyak8Corr(double puppipt, double puppieta ){


  double genCorr  = 1.;
  double recoCorr = 1.;
  double totalWeight = 1.;

  genCorr =  puppisd_corrGEN->Eval( puppipt );
  if( fabs(puppieta)  <= 1.3 ){
    recoCorr = puppisd_corrRECO_cen->Eval( puppipt );
  }
  else{
    recoCorr = puppisd_corrRECO_for->Eval( puppipt );
  }

  totalWeight = genCorr * recoCorr;

  return totalWeight;
}

