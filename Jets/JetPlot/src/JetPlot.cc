// -*- C++ -*-
//
// Package:    JetPlot
// Class:      JetPlot
// 
/**\class JetPlot JetPlot.cc Jets/JetPlot/src/JetPlot.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Andrew Evans
//         Created:  Wed May 28 09:18:51 CDT 2014
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

// GenJets
#include "DataFormats/JetReco/interface/GenJet.h"

// CMSSW
#include "FWCore/ServiceRegistry/interface/Service.h" // edm::Service
#include "CommonTools/UtilAlgos/interface/TFileService.h" // TFileService

// Root Library
#include "TH1D.h"

//
// class declaration
//

class JetPlot : public edm::EDAnalyzer {
   public:
      explicit JetPlot(const edm::ParameterSet&);
      ~JetPlot();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

      // ----------member data ---------------------------

      TH1D* jet_counter;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
JetPlot::JetPlot(const edm::ParameterSet& iConfig)

{
    //now do what ever initialization is needed
    // This allows us to make plots with file service
    edm::Service<TFileService> fs;

    //makes the jet plot
    jet_counter = fs->make<TH1D>("jetCounts", "jetCounts", 10, 0., 10);
}


JetPlot::~JetPlot()
{
 
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
JetPlot::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;

    // gets jets
    Handle<std::vector<reco::GenJet>> pIn;
    iEvent.getByLabel("ak5GenJets", pIn);

    // loops over events and fills plots

    int jet_number = 0;
    for(unsigned int i = 0; i < pIn->size(); ++i)
    {
	jet_number++;
    }
    jet_counter->Fill(jet_number);
}


// ------------ method called once each job just before starting event loop  ------------
void 
JetPlot::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
JetPlot::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
JetPlot::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
JetPlot::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
JetPlot::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
JetPlot::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
JetPlot::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(JetPlot);
