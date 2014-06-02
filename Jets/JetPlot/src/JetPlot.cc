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

// C libraries
#include <cmath>
#include <algorithm>

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
      
      // make histograms for various criteria
      TH1D* jet_counter_over20;                   //number of jets in an event with 
      TH1D* jet_pt_over20;
      TH1D* pt_abs_eta_under2p5;
      TH1D* pt_abs_eta_over2p5;
      TH1D* counter_abs_eta_under2p5;
      TH1D* counter_abs_eta_over2p5;
      TH1D* jet_eta_under2p5;
      TH1D* jet_phi_over2p5;
      TH1D* jet_phi_under2p5;
      TH1D* jet_eta_over2p5;
      TH1D* jet_mass_under2p5;
      TH1D* jet_mass_over2p5;
      TH1D* dijet_invmass;
      
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

    //makes the jet plots with fs
    jet_counter_over20 = fs->make<TH1D>("jetCounts", "jetCounts", 100, 0., 100);
    jet_pt_over20 = fs->make<TH1D>("jetPt", "jetPt", 100, 0., 100);
    pt_abs_eta_under2p5 = fs->make<TH1D>("ptAbsEtaUnder2.5", "ptAbsEtaUnder2.5", 100, 0., 100);
    pt_abs_eta_over2p5 = fs->make<TH1D>("ptAbsEtaOver2.5", "ptAbsEtaOver2.5", 100, 0., 100);
    counter_abs_eta_under2p5 = fs->make<TH1D>("counterAbsEtaUnder2.5", "counterAbsEtaUnder2.5", 100, 0., 100);
    counter_abs_eta_over2p5 = fs->make<TH1D>("counterAbsEtaOver2.5", "counterAbsEtaOver2.5", 100, 0., 100);
    jet_eta_under2p5 = fs->make<TH1D>("EtaDistributionUnder2.5", "EtaDistributionUnder2.5", 200, -10, 10);
    jet_phi_over2p5 = fs->make<TH1D>("PhiDistributionOver2.5", "PhiDistributionOver2.5", 200, -10, 10);
    jet_eta_over2p5 = fs->make<TH1D>("EtaDistributionOver2.5", "EtaDistributionOver2.5", 200, -10, 10);
    jet_phi_under2p5 = fs->make<TH1D>("PhiDistributionUnder2.5", "PhiDistributionUnder2.5", 200, -10, 10);
    jet_mass_under2p5 = fs->make<TH1D>("JetMassUnder2.5", "JetMassUnder2.5", 1000, 0., 100);
    jet_mass_over2p5 = fs->make<TH1D>("JetMassOver2.5", "JetMassOver2.5", 1000, 0., 100);
    dijet_invmass = fs->make<TH1D>("InvariantDijetMass", "InvariantDijetMass", 1000, 0., 1000);
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

    math::PtEtaPhiMLorentzVector diJet;  //makes a lorentz vectors for the two highest jets
    math::PtEtaPhiMLorentzVector jet1;
    math::PtEtaPhiMLorentzVector jet2;

    // gets jets
    Handle<std::vector<reco::GenJet>> pIn;
    iEvent.getByLabel("ak5GenJets", pIn);
    
    // loops over events and fills plots
    int jet_number_eta_under2p5 = 0;
    int jet_number_eta_over2p5 = 0;
    int jet_number_over20 = 0;
    float jet1pt = 0.0;
    float jet2pt = 0.0;
    bool cutPt;
    bool cutInside;
    bool cutIsolation;
    bool cutSpread;
    bool cutHemi;

    for(unsigned int i = 0; i < pIn->size(); ++i)
    {
        reco::GenJet genjet = pIn->at(i);

        if(genjet.pt() > jet1pt)     //Decides the highest and second highest pt jets
        {
            jet1pt = genjet.pt();
            jet1 = math::PtEtaPhiMLorentzVector(genjet.pt(), genjet.eta(), genjet.phi(), genjet.mass());
        }
        else if(genjet.pt() > jet2pt)
        {
            jet2pt = genjet.pt();
            jet2 = math::PtEtaPhiMLorentzVector(genjet.pt(), genjet.eta(), genjet.phi(), genjet.mass());
        }
        diJet = jet1 + jet2;
        
	//Now that I have the two highest pt jets I can make some cuts on them based on the Dubinin2006 paper
        if(jet1.pt() > 20 and jet2.pt() > 20)  //At least some jet energy
        {
            cutPt = true;
        }
        if(fabs (jet1.eta()) <= 4.5 and fabs (jet2.eta()) <= 4.5)   //Jets well inside detector (but not in HF I guess)
	{
            cutInside = true;
	}
	if(/*Separation between jets and gammas */true/*until I know what to write*/)
	{
            cutIsolation = true;
	}
	if(fabs (jet1.eta()-jet2.eta()) >= 4.0)    //Jets are far apart in eta
	{
            cutSpread = true;
	}
	if(jet1.eta()*jet2.eta() < 0)    //They are in opposite hemispheres
	{
            cutHemi = true;
	}

        //Checks if the absolute value of eta is over or under 2.5            
 	if(fabs (genjet.eta()) > 2.5)
        {
            jet_number_eta_over2p5++;
            pt_abs_eta_over2p5->Fill(genjet.pt());
            jet_eta_over2p5->Fill(genjet.eta());
	    jet_phi_over2p5->Fill(genjet.phi());
	    jet_mass_over2p5->Fill(genjet.mass());
	}
	else
	{
            jet_number_eta_under2p5++;
            pt_abs_eta_under2p5->Fill(genjet.pt());
            jet_eta_under2p5->Fill(genjet.eta());
	    jet_phi_under2p5->Fill(genjet.phi());
	    jet_mass_under2p5->Fill(genjet.mass());
	}
        //Checks if the jets have at least 20GeV pt
        if(genjet.pt() >= 20)
    	{    
	    jet_number_over20++;
    	    jet_pt_over20->Fill(genjet.pt());
        }
    }
    
    dijet_invmass->Fill(diJet.mass());
    jet_counter_over20->Fill(jet_number_over20);
    counter_abs_eta_under2p5->Fill(jet_number_eta_under2p5);
    counter_abs_eta_over2p5->Fill(jet_number_eta_over2p5);
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
