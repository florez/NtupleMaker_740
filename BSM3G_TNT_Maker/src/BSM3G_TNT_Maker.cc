// Authors:  Alfredo Gurrola (Vanderbilt University)
// Andres Florez: Universidad de los Andes, Colombia.
// kaur amandeepkalsi: Panjab University, India.

#include "NtupleMaker/BSM3G_TNT_Maker/interface/BSM3G_TNT_Maker.h"
#include <memory>
#include <string>
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

BSM3G_TNT_Maker::BSM3G_TNT_Maker(const edm::ParameterSet& iConfig):
  MaxN(200)
  //now do what ever initialization is needed
{
  debug_                   = iConfig.getParameter<bool>("debug_");
  _filltriggerinfo         = iConfig.getParameter<bool>("filltriggerinfo");
  _fillmuoninfo            = iConfig.getParameter<bool>("fillmuoninfo");
  _fillelectronpatinfo     = iConfig.getParameter<bool>("fillelectronpatinfo");
  _filltauinfo             = iConfig.getParameter<bool>("filltauinfo");
  _filljetinfo             = iConfig.getParameter<bool>("filljetinfo");
  _fillgeninfo             = iConfig.getParameter<bool>("fillgeninfo");
  _fillgenweightinfo       = iConfig.getParameter<bool>("fillgenweightinfo");
  _fillruninfo             = iConfig.getParameter<bool>("fillruninfo");
  _fillPVinfo              = iConfig.getParameter<bool>("fillPVinfo");
  _fillMETinfo             = iConfig.getParameter<bool>("fillMETinfo");
  _fillphotoninfo          = iConfig.getParameter<bool>("fillphotoninfo");
  _fillPDFinfo             = iConfig.getParameter<bool>("fillPDFinfo");
  edm::Service<TFileService> fs;
  fs->file().SetCompressionLevel(9);
  tree_ = fs->make<TTree>("BOOM","BOOM");
  storetree_ = fs->make<TTree>("BAAM","BAAM");


  storetree_->Branch("pdf_name", &pdf_name);
  storetree_->Branch("triggernames", &triggernames);



  if( _filltriggerinfo)      trselector        = new TriggerSelector("miniAOD", tree_, debug_, iConfig, consumesCollector(),triggernames);
  if( _fillmuoninfo)         muselector        = new MuonSelector("miniAOD", tree_, debug_, iConfig, consumesCollector());
  if( _fillelectronpatinfo)  elpatselector     = new ElectronPatSelector("miniAOD", tree_, debug_, iConfig, consumesCollector());
  if( _filltauinfo)          tauselector       = new TauSelector("miniAOD", tree_, debug_, iConfig, consumesCollector());
  if( _filljetinfo)          jetselector       = new JetSelector("miniAOD", tree_, debug_, iConfig, consumesCollector());
  if( _fillgeninfo)          genselector       = new GenParticleSelector("miniAOD", tree_, debug_, iConfig, consumesCollector());
  if( _fillgenweightinfo)    genweightselector = new GenEventWeightSelector("miniAOD", tree_, debug_, iConfig, consumesCollector());
  if( _fillruninfo)          runinfoselector   = new RunInfoSelector("miniAOD", tree_, debug_, iConfig, consumesCollector());
  if( _fillPVinfo)           pvselector        = new PVSelector("miniAOD", tree_, debug_, iConfig, consumesCollector());
  if( _fillMETinfo)          metselector       = new METSelector("miniAOD", tree_, debug_, iConfig, consumesCollector());
  if(_fillphotoninfo)        photonselector    = new PhotonSelector("miniAOD", tree_, debug_, iConfig, consumesCollector());
  if(_fillPDFinfo)           pdfselector       = new PDFSelector("miniAOD", tree_, debug_, iConfig, consumesCollector(),pdf_name);
}


BSM3G_TNT_Maker::~BSM3G_TNT_Maker()
{
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//

// ------------ method called for each event  ------------
void
BSM3G_TNT_Maker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  using namespace edm;
  using namespace pat;
  using namespace reco;

  if( _filltriggerinfo)      trselector->Fill(iEvent, iSetup);
  if( _fillmuoninfo)         muselector->Fill(iEvent, iSetup);
  if( _fillelectronpatinfo)  elpatselector->Fill(iEvent, iSetup);
  if( _filltauinfo)          tauselector->Fill(iEvent, iSetup);
  if( _filljetinfo)          jetselector->Fill(iEvent);
  if( _fillgeninfo)          genselector->Fill(iEvent);
  if( _fillgenweightinfo)    genweightselector->Fill(iEvent);
  if( _fillruninfo)          runinfoselector->Fill(iEvent);
  if( _fillPVinfo)           pvselector->Fill(iEvent);
  if( _fillMETinfo)          metselector->Fill(iEvent);
  if( _fillphotoninfo)       photonselector->Fill(iEvent);
  if( _fillPDFinfo)          pdfselector->Fill(iEvent);

  tree_->Fill();
#ifdef THIS_IS_AN_EVENT_EXAMPLE
  Handle<ExampleData> pIn;
  iEvent.getByLabel("example",pIn);
#endif

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  ESHandle<SetupData> pSetup;
  iSetup.get<SetupRecord>().get(pSetup);
#endif
}

// ------------ method called once each job just before starting event loop  ------------
void
BSM3G_TNT_Maker::beginJob()
{


}

// ------------ method called once each job just after ending the event loop  ------------
void
BSM3G_TNT_Maker::endJob()
{
  storetree_->Fill();
}

// ------------ method called when starting to processes a run  ------------

void
BSM3G_TNT_Maker::beginRun(edm::Run const & iRun, edm::EventSetup const& iSetup)
  {
    if( _filltriggerinfo) trselector->startTrigger(iSetup, iRun);
    if( _fillPDFinfo)     pdfselector->beginRun(iRun, iSetup);
  }

// ------------ method called when ending the processing of a run  ------------
/*
  void
  BSM3G_TNT_Maker::endRun(edm::Run const&, edm::EventSetup const&)
  {
  }
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
  void
  BSM3G_TNT_Maker::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
  {
  }
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
  void
  BSM3G_TNT_Maker::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
  {
  }
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
BSM3G_TNT_Maker::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(BSM3G_TNT_Maker);
