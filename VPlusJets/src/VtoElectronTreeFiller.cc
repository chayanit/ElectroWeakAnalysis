/*****************************************************************************
 * Project: CMS detector at the CERN
 *
 * Package: ElectroWeakAnalysis/VPlusJets
 *
 *
 * Authors:
 *
 *   Kalanand Mishra, Fermilab - kalanand@fnal.gov
 *
 * Description:
 *   To fill W--> enu or Z-->ee  related quantities into a specified TTree
 *   Can work with jets in RECO/AOD/PAT data formats.
 * History:
 *   
 *
 * Copyright (C) 2010 FNAL 
 *****************************************************************************/

// CMS includes
#include "DataFormats/Candidate/interface/ShallowCloneCandidate.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

// Header file
#include "ElectroWeakAnalysis/VPlusJets/interface/VtoElectronTreeFiller.h"
#include "ElectroWeakAnalysis/VPlusJets/interface/METzCalculator.h"
#include "ElectroWeakAnalysis/VPlusJets/interface/ElectronEffectiveArea.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"


ewk::VtoElectronTreeFiller::VtoElectronTreeFiller(const char *name, TTree* tree, 
							const edm::ParameterSet iConfig)
{
  // To Fill Z and Z's daughters information

  // ********** Vector boson ********** //
  if(  iConfig.existsAs<edm::InputTag>("srcVectorBoson") )
    mInputBoson = iConfig.getParameter<edm::InputTag>("srcVectorBoson"); 
  else std::cout << "***Error:" << name << 
    " Collection not specified !" << std::endl;
    if(  iConfig.existsAs<edm::InputTag>("srcMet") )
		mInputMet = iConfig.getParameter<edm::InputTag>("srcMet");
    if(  iConfig.existsAs<edm::InputTag>("srcBeamSpot") )
		mInputBeamSpot  = iConfig.getParameter<edm::InputTag>("srcBeamSpot");
  tree_     = tree;
  name_     = name;
  Vtype_    = iConfig.getParameter<std::string>("VBosonType"); 
  LeptonType_ = iConfig.getParameter<std::string>("LeptonType");
  if( !(tree==0) && LeptonType_=="electron") SetBranches();
}


void ewk::VtoElectronTreeFiller::SetBranches()
{
  // Declare jet branches
  std::string lept1 = "eplus";
  std::string lept2 = "eminus";
  if( !(Vtype_=="Z") ) lept1 = "electron";

  SetBranch( &V_mass,      "mass");
  SetBranch( &V_px,        "px");
  SetBranch( &V_py,        "py");
  SetBranch( &V_pz,        "pz");
  SetBranch( &V_E,         "e");
  SetBranch( &V_Pt,        "pt");
  SetBranch( &V_Et,        "et");
  SetBranch( &V_Eta,       "eta");    
  SetBranch( &V_Phi,       "phi");
  if(Vtype_=="W") {
    SetBranch( &V_mt,        "mt");
    SetBranch( &V_pzNu1,     "pzNu1");
    SetBranch( &V_pzNu2,     "pzNu2");
  }
  ///////////////////////////////////////////////
  SetBranch( &e1px,             lept1+"_px" );
  SetBranch( &e1py,             lept1+"_py" );
  SetBranch( &e1pz,             lept1+"_pz" );
  SetBranch( &e1E,              lept1+"_e" );
  SetBranch( &e1Pt,             lept1+"_pt" );
  SetBranch( &e1Et,             lept1+"_et" );
  SetBranch( &e1Eta,            lept1+"_eta" ); 
  SetBranch( &e1Phi,            lept1+"_phi" );
  SetBranch( &e1Charge,         lept1+"_charge" );
  SetBranch( &e1_trackiso,      lept1+"_trackiso" );
  SetBranch( &e1_hcaliso,       lept1+"_hcaliso" );
  SetBranch( &e1_ecaliso,       lept1+"_ecaliso" );
  SetBranch( &e1_sc_Eta,        lept1+"_sc_Eta" );
  SetBranch( &e1_sc_Phi,        lept1+"_sc_Phi" );
  SetBranch( &e1_sc_E,          lept1+"_sc_E" );
  SetBranch( &e1_numberOfBrems, lept1+"_numbrem" );
  SetBranch( &e1_BremFraction,  lept1+"_fbrem" );
  SetBranch( &e1_DeltaEtaIn,    lept1+"_deltaeta_in" );
  SetBranch( &e1_DeltaPhiIn,    lept1+"_deltaphi_in" );
  SetBranch( &e1_HoverE,        lept1+"_hovere" );	    
  SetBranch( &e1_SigmaIetaIeta, lept1+"_sigmaietaieta" );
  SetBranch( &e1_missingHits,   lept1+"_missingHits" );	  
  SetBranch( &e1_ooemoop,       lept1+"_ooemoop" );

  ////////////////////////////////////////////////////////
  if(Vtype_=="Z") {	  
    SetBranch( &e2px,             lept2+"_px" );
    SetBranch( &e2py,             lept2+"_py" );
    SetBranch( &e2pz,             lept2+"_pz" );
    SetBranch( &e2E,              lept2+"_e" );
    SetBranch( &e2Pt,             lept2+"_pt" );
    SetBranch( &e2Et,             lept2+"_et" );
    SetBranch( &e2Eta,            lept2+"_eta" ); 
    SetBranch( &e2Phi,            lept2+"_phi" );
    SetBranch( &e2Charge,         lept2+"_charge" );
    SetBranch( &e2_trackiso,      lept2+"_trackiso" );
    SetBranch( &e2_hcaliso,       lept2+"_hcaliso" );
    SetBranch( &e2_ecaliso,       lept2+"_ecaliso");
    SetBranch( &e2_sc_Eta,        lept2+"_sc_Eta" );
    SetBranch( &e2_sc_Phi,        lept2+"_sc_Phi" );
    SetBranch( &e2_sc_E,          lept2+"_sc_E" );
    SetBranch( &e2_numberOfBrems, lept2+"_numbrem" );
    SetBranch( &e2_BremFraction,  lept2+"_fbrem" );
    SetBranch( &e2_DeltaEtaIn,    lept2+"_deltaeta_in" );
    SetBranch( &e2_DeltaPhiIn,    lept2+"_deltaphi_in" );
    SetBranch( &e2_HoverE,        lept2+"_hovere" );	   	  
    SetBranch( &e2_SigmaIetaIeta, lept2+"_sigmaietaieta" );
    SetBranch( &e2_missingHits,   lept2+"_missingHits" );
    SetBranch( &e2_ooemoop,       lept2+"_ooemoop" );
  }	  
}
/////////////////////////////////////////////////////////////////////////


void ewk::VtoElectronTreeFiller::init()   
{
  // initialize private data members
  V_mass                = -1.;
  V_mt                  = -1.;
  V_px                  = -99999.;
  V_py                  = -99999.;
  V_pz                  = -99999.;
  V_E                   = -1.;
  V_Pt                  = -1.;
  V_Et                  = -1.;
  V_Eta                 = -10.;
  V_Phi                 = -10.;
  V_pzNu1               = -10000.0;
  V_pzNu2               = -10000.0;

  // W/Z daughter
  e1Charge           = -10;
  e2Charge           = -10;
  e1px               = -99999.;
  e1py               = -99999.;
  e1pz               = -99999.;
  e1E                = -1.;
  e1Et               = -1.;
  e1Pt               = -1.;
  e1Eta              = -10.;
  e1Phi              = -10.;
  e1_sc_Eta          = -10.;    
  e1_sc_Phi          = -10.;
  e1_sc_E            = -1.;
  e1_numberOfBrems   = -10.;
  e1_BremFraction    = -1.;
  e1_DeltaEtaIn      = -10.;
  e1_DeltaPhiIn      = -10.;
  e1_HoverE          = -1.;	   	  
  e1_SigmaIetaIeta   = -1.;	  
  e1_missingHits     = 100;
  e1_ooemoop	     = -10.;

  e2px              = -99999.;
  e2py              = -99999.;
  e2pz              = -99999.;
  e2E               = -1.;
  e2Pt              = -1.;
  e2Et              = -1.;
  e2Eta             = -10.;
  e2Phi             = -10.;
  e2_sc_Eta         = -10.;    
  e2_sc_Phi         = -10.;
  e2_sc_E           = -1.;
  e2_numberOfBrems   = -10.;
  e2_BremFraction    = -1.;
  e2_DeltaEtaIn      = -10.;
  e2_DeltaPhiIn      = -10.;
  e2_HoverE          = -1.;	    
  e2_SigmaIetaIeta     = -1.;	  
  e2_missingHits     = 100;
  e2_ooemoop         = -10.;
	  
  //////////////
  e1_trackiso     = 5000.0;
  e1_hcaliso      = 5000.0;
  e1_ecaliso      = 5000.0;
  e2_trackiso     = 5000.0;
  e2_hcaliso      = 5000.0;
  e2_ecaliso      = 5000.0;

  // initialization done
}

void ewk::VtoElectronTreeFiller::fill(const edm::Event& iEvent, int vecBosonIndex)
{
  // protection
  if( (tree_==0) || !(LeptonType_=="electron") )  return;

  // first initialize to the default values
  init();

  edm::Handle<reco::CandidateView> boson;
  iEvent.getByLabel( mInputBoson, boson);
  if( boson->size()<1 ) return; // Nothing to fill
  
  edm::Handle<edm::View<reco::MET> > pfmet;
  iEvent.getByLabel(mInputMet, pfmet);

 /////// Pileup density "rho" in the event from fastJet pileup calculation /////
  //edm::Handle<double> rho;
  //const edm::InputTag eventrho("kt6PFJets", "rho");
  //iEvent.getByLabel(eventrho,rho);
  //double fastJetRho = *rho;

  const reco::Candidate *Vboson = &((*boson)[vecBosonIndex]); 
  if( Vboson == 0) return;

  ////////// Vector boson quantities //////////////
  V_mass 	= Vboson->mass();
  if(Vtype_=="W") V_mt = sqrt(2.0*Vboson->daughter(0)->pt()*Vboson->daughter(1)->pt()*(1.0-cos(Vboson->daughter(0)->phi()-Vboson->daughter(1)->phi())));
  V_Eta 	= Vboson->eta();   
  V_Phi 	= Vboson->phi();
  V_px 		= Vboson->px();
  V_py 		= Vboson->py();
  V_pz 		= Vboson->pz();
  V_E  		= Vboson->energy();
  V_Pt 		= Vboson->pt();
  V_Et 		= Vboson->et();

  // now iterate over the daughters  
  if(Vboson->numberOfDaughters()<2 ) {
    throw cms::Exception( "***Error: V boson has < 2 daughters !\n");
    return;  // if no electron found, then return
  } 
  std::cout << "number of daugters " << Vboson->numberOfDaughters() <<std::endl;
  // get the two daughters
  reco::CandidateBaseRef m0 = Vboson->daughter(0)->masterClone();
  reco::CandidateBaseRef m1 = Vboson->daughter(1)->masterClone();
  //std::cout << "ptr m1 is nonnull " << m1.isNonnull() <<std::endl ;

  const reco::GsfElectron* e1=NULL;
  const reco::GsfElectron* e2=NULL;
  const std::type_info & type0 = typeid(*m0);
  const std::type_info & type1 = typeid(*m1);

  if( type0 == typeid(pat::Electron) ){
      e1 = dynamic_cast<const pat::Electron *>(&*m0) ;
  }
  if( type1 == typeid(pat::Electron) ){
      e2 = dynamic_cast<const pat::Electron *>(&*m1) ;
  }

  if(0==e1 && 0==e2) {
    throw cms::Exception("***Error: couldn't do dynamic cast of vector boson daughters !\n");
    return;  // if no electron found, then return
  }

  const reco::GsfElectron* ele1=NULL;
  const reco::GsfElectron* ele2=NULL;
		 
  // if Z--> e+e- then ele1 = e+, ele2 = e-
  if(Vtype_=="Z") {
    if(e1->charge() > 0) {  ele1 = e1;   ele2 = e2; }
    else { ele1 = e2;  ele2 = e1; }
  }
  // if W--> enu then ele1 = e, ele2 = NULL 
  if(Vtype_=="W") {
    if( abs(e1->charge())==1 ) ele1  = e1;
    else if( abs(e2->charge())==1 ) ele1  = e2;

    if( !(ele1 == NULL) ) {
      // estimate Pz of neutrino
      TLorentzVector p4MET((*pfmet)[0].px(), (*pfmet)[0].py(), (*pfmet)[0].pz(), (*pfmet)[0].energy());
      TLorentzVector p4lepton(ele1->px(), ele1->py(), ele1->pz(), ele1->energy());
      METzCalculator metz;
      metz.SetMET(p4MET);
      metz.SetLepton(p4lepton);
      if (LeptonType_=="electron") metz.SetLeptonType("electron");
      if (LeptonType_=="muon")     metz.SetLeptonType("muon");
      V_pzNu1 = metz.Calculate();
      V_pzNu2 = metz.getOther();
    }
  }

	  
  ////////// electron #1 quantities //////////////
  if( !(ele1 == NULL) ) {
    e1Charge           = ele1->charge();
    e1Eta              = ele1->eta();    
    e1Phi              = ele1->phi();
    e1E                = ele1->energy();
    e1px               = ele1->px();
    e1py               = ele1->py();
    e1pz               = ele1->pz();
    e1Pt               = ele1->pt();
    e1Et               = ele1->et();

	  /// isolation 
    e1_trackiso       = ele1->dr03TkSumPt();
    e1_ecaliso        = ele1->dr03EcalRecHitSumEt();
    e1_hcaliso        = ele1->dr03HcalTowerSumEt();
    e1_numberOfBrems  = ele1->numberOfBrems();      
    e1_BremFraction   = ele1->fbrem();
    e1_DeltaEtaIn     = ele1->deltaEtaSuperClusterTrackAtVtx();
    e1_DeltaPhiIn     = ele1->deltaPhiSuperClusterTrackAtVtx();
    //get Hcal energy over Ecal Energy
    e1_HoverE 	      = ele1->hadronicOverEm();	  
    //get SuperCluster (sc) infos
    reco::SuperClusterRef SCp = ele1->superCluster();
    e1_sc_Eta        = SCp->eta();    
    e1_sc_Phi        = SCp->phi();
    e1_sc_E          = SCp->energy();
    e1_SigmaIetaIeta = ele1->sigmaIetaIeta();
    e1_missingHits   = ele1->gsfTrack()->trackerExpectedHitsInner().numberOfHits();
    e1_ooemoop 	     = fabs((1.0/ele1->ecalEnergy() - ele1->eSuperClusterOverP()/ele1->ecalEnergy()));
  }

  ////////// electron #2 quantities //////////////
  if( !(ele2 == NULL) ) {
    e2Charge          = ele2->charge();
    e2Eta             = ele2->eta();    
    e2Phi             = ele2->phi();
    e2E               = ele2->energy();
    e2px              = ele2->px();
    e2py              = ele2->py();
    e2pz              = ele2->pz();
    e2Pt              = ele2->pt();
    e2Et              = ele2->et();	  

    /// isolation 
    e2_trackiso       = ele2->dr03TkSumPt();
    e2_ecaliso        = ele2->dr03EcalRecHitSumEt();
    e2_hcaliso        = ele2->dr03HcalTowerSumEt();
    e2_numberOfBrems  = ele2->numberOfBrems();
    e2_BremFraction   = ele2->fbrem();
    e2_DeltaEtaIn     = ele2->deltaEtaSuperClusterTrackAtVtx();
    e2_DeltaPhiIn     = ele2->deltaPhiSuperClusterTrackAtVtx();
    //get Hcal energy over Ecal Energy
    e2_HoverE = ele2->hadronicOverEm();	  
    //get SuperCluster (sc) infos
    reco::SuperClusterRef SCm = ele2->superCluster();
    e2_sc_Eta        = SCm->eta();    
    e2_sc_Phi        = SCm->phi();
    e2_sc_E          = SCm->energy();
    e2_SigmaIetaIeta = ele2->sigmaIetaIeta();
    e2_missingHits   = ele2->gsfTrack()->trackerExpectedHitsInner().numberOfHits();
    e2_ooemoop       = fabs((1.0/ele2->ecalEnergy() - ele2->eSuperClusterOverP()/ele2->ecalEnergy()));
  } 

}


////////////////// utilities, helpers ///////////////////
 
void ewk::VtoElectronTreeFiller::SetBranch( float* x, std::string name)
{
  std::string brName = std::string(name_) + "_" + name;
  tree_->Branch( brName.c_str(), x, ( brName+"/F").c_str() );
}


void ewk::VtoElectronTreeFiller::SetBranch( int* x, std::string name)
{
  std::string brName = std::string(name_) + "_" + name;
  tree_->Branch( brName.c_str(), x, ( brName+"/I").c_str() );
}


void ewk::VtoElectronTreeFiller::SetBranch( bool* x, std::string name)
{
  std::string brName = std::string(name_) + "_" + name;
  tree_->Branch( brName.c_str(), x, ( brName+"/O").c_str() );
}

void ewk::VtoElectronTreeFiller::SetBranchArray( float* x, std::string name)
{
  tree_->Branch( name.c_str(), x, ( name+"[10]/F").c_str() );
}



