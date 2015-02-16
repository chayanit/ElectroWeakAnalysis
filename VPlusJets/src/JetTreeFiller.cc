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
 *   To fill jet related quantities into a specified TTree
 *   Can work with CaloJet, GenJet, JPT jet, PF jet.
 *   Can work with jets in RECO/AOD/PAT data formats.
 * History:
 *   
 *
 * Copyright (C) 2010 FNAL 
 *****************************************************************************/

// CMS includes
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "TMath.h" 
#include <TLorentzVector.h>
#include "JetMETCorrections/MCJet/plugins/JetUtilMC.h" // needed for dPhi,dR

// Monte Carlo stuff
// #include "SimDataFormats/JetMatching/interface/JetFlavour.h"
// #include "SimDataFormats/JetMatching/interface/JetFlavourMatching.h"
// #include "SimDataFormats/JetMatching/interface/MatchedPartons.h"
// #include "SimDataFormats/JetMatching/interface/JetMatchedPartons.h"


// Header file
#include "ElectroWeakAnalysis/VPlusJets/interface/JetTreeFiller.h"

const float BTAG_DISCRIM_DEFAULT=-999.;

ewk::JetTreeFiller::JetTreeFiller(const char *name, TTree* tree, 
				  const std::string jetType,
				  const edm::ParameterSet iConfig)
{
  // ********** GenJets ********** //
  //if( jetType=="Gen" && iConfig.existsAs<edm::InputTag>("srcGen") )
  //  mInputJets = iConfig.getParameter<edm::InputTag>("srcGen"); 
  // ********** PFJets ********** //
  if( jetType=="Pat" && iConfig.existsAs<edm::InputTag>("srcPatJets") )
    mInputJets = iConfig.getParameter<edm::InputTag>("srcPatJets"); 
  // ********** Corrected PFJets ********** //
  if( jetType=="Selected" && iConfig.existsAs<edm::InputTag>("srcSelectedJets") )
    mInputJets = iConfig.getParameter<edm::InputTag>("srcSelectedJets"); 
  // ********** Corrected PFJets for VBF Tag ********** //
  if( jetType=="Smeared" && iConfig.existsAs<edm::InputTag>("srcSmearedJets") )
    mInputJets = iConfig.getParameter<edm::InputTag>("srcSmearedJets"); 

  // ********** Vector boson ********** //
  if(  iConfig.existsAs<edm::InputTag>("srcVectorBoson") )
    mInputBoson = iConfig.getParameter<edm::InputTag>("srcVectorBoson"); 

  //*********************  Run Over AOD or PAT  ***********//
  if( iConfig.existsAs<bool>("runningOverAOD"))
    runoverAOD = iConfig.getParameter<bool>("runningOverAOD");
		
		
  tree_     = tree;
  jetType_ = jetType;

  Vtype_    = iConfig.getParameter<std::string>("VBosonType"); 
  LeptonType_ = iConfig.getParameter<std::string>("LeptonType");

  if( !(tree==0) ) SetBranches();
}




void ewk::JetTreeFiller::SetBranches()
{
  // Declare jet branches
  SetBranchSingle( &NumJets, "num" + jetType_ + "Jets");
  SetBranchSingle( &numBTags, "num" + jetType_ + "JetBTags");
  SetBranch( Et, "Jet" + jetType_ + "_Et");
  SetBranch( Pt, "Jet" + jetType_ + "_Pt");
  SetBranch( Eta, "Jet" + jetType_ + "_Eta");
  SetBranch( Phi, "Jet" + jetType_ + "_Phi");
  SetBranch( Theta, "Jet" + jetType_ + "_Theta");
  SetBranch( Px, "Jet" + jetType_ + "_Px");
  SetBranch( Py, "Jet" + jetType_ + "_Py");
  SetBranch( Pz, "Jet" + jetType_ + "_Pz");
  SetBranch( E, "Jet" + jetType_ + "_E");
  SetBranch( Y, "Jet" + jetType_ + "_Y");
  SetBranch( Mass, "Jet" + jetType_ + "_Mass");
/*
  SetBranch( bDiscriminator, "Jet" + jetType_ + "_bDiscriminator");
  SetBranch( bDiscriminatorSSVHE, "Jet" + jetType_ + "_bDiscriminatorSSVHE");
  SetBranch( bDiscriminatorTCHE, "Jet" + jetType_ + "_bDiscriminatorTCHE");
  SetBranch( bDiscriminatorCSV, "Jet" + jetType_ + "_bDiscriminatorCSV");
  SetBranch( bDiscriminatorJP, "Jet" + jetType_ + "_bDiscriminatorJP");
  SetBranch( bDiscriminatorSSVHP, "Jet" + jetType_ + "_bDiscriminatorSSVHP");
  SetBranch( bDiscriminatorTCHP, "Jet" + jetType_ + "_bDiscriminatorTCHP");
*/
  /////////////////////////////////////////////////////////////////////////

}


//////////////////////////////////////////////////////////////////
/////// Helper for above function ////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////

void ewk::JetTreeFiller::SetBranchSingle( float* x, std::string name)
{
  tree_->Branch( name.c_str(), x, ( name+"/F").c_str() );
  bnames.push_back( name );
}

void ewk::JetTreeFiller::SetBranchSingle( int* x, std::string name)
{
  tree_->Branch( name.c_str(), x, ( name+"/I").c_str() );
  bnames.push_back( name );
}



void ewk::JetTreeFiller::SetBranch( float* x, std::string name)
{
  tree_->Branch( name.c_str(), x, ( name+"[20]/F").c_str() );
  bnames.push_back( name );
}


void ewk::JetTreeFiller::SetBranch( int* x, std::string name)
{
  tree_->Branch( name.c_str(), x, ( name+"[20]/I").c_str() );
  bnames.push_back( name );
}

//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////


void ewk::JetTreeFiller::FillBranches() const 
{

  for(std::vector<std::string>::iterator it = bnames.begin();
      it != bnames.end(); ++it) {
    if(TBranch *br = tree_->GetBranch( (*it).c_str() ) ) br->Fill();
  }

}




void ewk::JetTreeFiller::init()   
{
  // initialize private data members
  NumJets = 0; 
  numBTags = 0;

  for (int j =0; j< NUM_JET_MAX; ++j) {
    Et[j] = -10.0;
    Pt[j] = -10.0;
    Eta[j] = -10.0;
    Phi[j] = -10.0;
    Theta[j] = -10.0;
    E[j] = -10.0;
    Y[j] = -10.0;
    Px[j] = -999999.9;
    Py[j] = -999999.9;
    Pz[j] = -999999.9;
    Mass[j] = -10.0;
    bDiscriminator[j] = BTAG_DISCRIM_DEFAULT;
    bDiscriminatorSSVHE[j] = BTAG_DISCRIM_DEFAULT;
    bDiscriminatorTCHE[j] = BTAG_DISCRIM_DEFAULT;
    bDiscriminatorCSV[j] = BTAG_DISCRIM_DEFAULT;
    bDiscriminatorJP[j] = BTAG_DISCRIM_DEFAULT;
    bDiscriminatorSSVHP[j] = BTAG_DISCRIM_DEFAULT;
    bDiscriminatorTCHP[j] = BTAG_DISCRIM_DEFAULT;
  }
  // initialization done
   
}



void ewk::JetTreeFiller::fill(const edm::Event& iEvent)
{
  // first initialize to the default values
  init();

  edm::Handle<reco::CandidateView> boson;
  iEvent.getByLabel( mInputBoson, boson);
  const int nBoson = boson->size();
  if( nBoson<1 ) return; // Nothing to fill
  
  //const reco::Candidate *Vboson = &((*boson)[0]); 

  edm::Handle<edm::View<reco::Jet> > jets;
  iEvent.getByLabel( mInputJets, jets ); 


  if(jets->size() < 1) return;

  size_t iJet = 0;

  // Loop over reco jets 
  edm::View<reco::Jet>::const_iterator jet, endpjets = jets->end(); 
  for (jet = jets->begin();  jet != endpjets;  ++jet, ++iJet) {
    if( !(iJet< (unsigned int) NUM_JET_MAX) ) break;

    //--------- store jet 4-vectors ---
    fillBasicJetQuantities(iJet, *jet);

    // --------- Fill b-tag information -----------
    edm::Ptr<reco::Jet> ptrJet = jets->ptrAt( jet - jets->begin() );		  
    if ( ptrJet.isNonnull() && ptrJet.isAvailable() ) {
	const pat::Jet* pjet = dynamic_cast<const pat::Jet *>(ptrJet.get()) ;
	fillBtagInfoPAT( iJet, pjet);
    }
  }// close jets iteration loop

  NumJets = (int) iJet;

  //FillBranches();
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////--------- store jet 4-vectors ---
template <typename T1> 
void ewk::JetTreeFiller::fillBasicJetQuantities(int iJet, const T1& pfjet) 
{
  Et[iJet] 	= pfjet.et();
  Pt[iJet] 	= pfjet.pt();
  Eta[iJet] 	= pfjet.eta();
  Phi[iJet] 	= pfjet.phi();
  Theta[iJet] 	= pfjet.theta();
  Px[iJet] 	= pfjet.px();
  Py[iJet] 	= pfjet.py();
  Pz[iJet] 	= pfjet.pz();
  E[iJet]  	= pfjet.energy();
  Y[iJet]  	= pfjet.rapidity();
  Mass[iJet] 	= pfjet.mass();
}
template 
void ewk::JetTreeFiller::fillBasicJetQuantities(int, const reco::PFJet&);
template 
void ewk::JetTreeFiller::fillBasicJetQuantities(int, const pat::Jet&); 

//////--------- Fill b-tag information from PAT ---
void ewk::JetTreeFiller::fillBtagInfoPAT(int iJet, const pat::Jet* pjet) 
{
  if(pjet !=0)
    {
      bDiscriminatorSSVHE[iJet] = (*pjet).bDiscriminator("simpleSecondaryVertexHighEffBJetTags");
      bDiscriminatorTCHE[iJet] = (*pjet).bDiscriminator("trackCountingHighEffBJetTags");
      bDiscriminatorCSV[iJet] = (*pjet).bDiscriminator("combinedSecondaryVertexBJetTags");
      bDiscriminatorJP[iJet] = (*pjet).bDiscriminator("jetProbabilityBJetTags");
      bDiscriminatorSSVHP[iJet] = (*pjet).bDiscriminator("simpleSecondaryVertexHighPurBJetTags");
      bDiscriminatorTCHP[iJet] = (*pjet).bDiscriminator("trackCountingHighPurBJetTags");
      bDiscriminator[iJet] = bDiscriminatorSSVHE[iJet];
      if(bDiscriminator[iJet]>1.74) numBTags ++;
    }
}
