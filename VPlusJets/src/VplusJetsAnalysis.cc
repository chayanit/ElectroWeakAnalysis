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


// user include files
#include "ElectroWeakAnalysis/VPlusJets/interface/VplusJetsAnalysis.h" 
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/METCollection.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/METReco/interface/GenMETCollection.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h" 
// PAT candidates
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/MHT.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"


ewk::VplusJetsAnalysis::VplusJetsAnalysis(const edm::ParameterSet& iConfig) :
  //fOutputFileName ( iConfig.getParameter<std::string>("HistOutFile") ),
  //hOutputFile ( new TFile( fOutputFileName.c_str(), "RECREATE" ) ), 
  //myTree ( new TTree(iConfig.getParameter<std::string>("TreeName").c_str(),"V+jets Tree") ),
  myTree ( fs->make<TTree>(iConfig.getParameter<std::string>("TreeName").c_str(),"V+jets Tree") ),
  CorrectedPatJetFiller ( iConfig.existsAs<edm::InputTag>("srcPatJets") ?
  new JetTreeFiller("CorrectedPatJetFiller", myTree, "Pat", iConfig) : 0),
  CorrectedSelectedJetFiller ( iConfig.existsAs<edm::InputTag>("srcSelectedJets") ? 
  new JetTreeFiller("CorrectedSelectedJetFiller", myTree, "Selected", iConfig) : 0), 
  SmearedJetFiller ((iConfig.existsAs<edm::InputTag>("srcSmearedJets") && 
		iConfig.existsAs<bool>("runningOverMC") &&
		iConfig.getParameter<bool>("runningOverMC")) ?
  new JetTreeFiller("SmearedJetFiller", myTree, "Smeared", iConfig) : 0),
  //ElectronFiller ( (iConfig.existsAs<edm::InputTag>("srcSelectElectrons") && iConfig.existsAs<edm::InputTag>("srcElectrons")) ?
  //new ElectronTreeFiller( myTree, iConfig ) : 0),	
  recoBosonFillerE( new VtoElectronTreeFiller( iConfig.getParameter<std::string>("VBosonType").c_str(), myTree, iConfig) ),
  genBosonFiller( (iConfig.existsAs<bool>("runningOverMC") && 
  		iConfig.getParameter<bool>("runningOverMC")) ?
  new MCTreeFiller(iConfig.getParameter<std::string>("VBosonType").c_str(), myTree, iConfig) : 0)
{

  // Are we running over Monte Carlo ?
   if( iConfig.existsAs<bool>("runningOverMC") ) 
      runningOverMC_=iConfig.getParameter< bool >("runningOverMC");
   else runningOverMC_= false;
  if(  iConfig.existsAs<edm::InputTag>("srcVectorBoson") )
    mInputBoson = iConfig.getParameter<edm::InputTag>("srcVectorBoson"); 
  LeptonType_ = iConfig.getParameter<std::string>("LeptonType");
  VBosonType_ = iConfig.getParameter<std::string>("VBosonType");
  if(  iConfig.existsAs<edm::InputTag>("srcPrimaryVertex") )
    mPrimaryVertex = iConfig.getParameter<edm::InputTag>("srcPrimaryVertex"); 
  else mPrimaryVertex =  edm::InputTag("offlinePrimaryVertices");
  if(  iConfig.existsAs<edm::InputTag>("srcBeamSpot") )
    mInputBeamSpot = iConfig.getParameter<edm::InputTag>("srcBeamSpot"); 
  if(  iConfig.existsAs<edm::InputTag>("srcMet") )
    mInputMet = iConfig.getParameter<edm::InputTag>("srcMet");

//#################### CaloMET&RAWPF ##############
  if(  iConfig.existsAs<edm::InputTag>("srcMetRaw") )
    mMetRaw = iConfig.getParameter<edm::InputTag>("srcMetRaw");
  if(  iConfig.existsAs<edm::InputTag>("srcMetCalo") )
    mMetCalo = iConfig.getParameter<edm::InputTag>("srcMetCalo");
  if(  iConfig.existsAs<edm::InputTag>("srcMetCaloT1") )
    mMetCaloT1 = iConfig.getParameter<edm::InputTag>("srcMetCaloT1");
//#################### PFMET ######################
/*
  if(  iConfig.existsAs<edm::InputTag>("srcMetCorrected") )
    mMetType01 = iConfig.getParameter<edm::InputTag>("srcMetCorrected");
  if(  iConfig.existsAs<edm::InputTag>("srcMetJERup") )
    mMetJERup = iConfig.getParameter<edm::InputTag>("srcMetJERup");
  if(  iConfig.existsAs<edm::InputTag>("srcMetJERdown") )
    mMetJERdown = iConfig.getParameter<edm::InputTag>("srcMetJERdown");
  if(  iConfig.existsAs<edm::InputTag>("srcMetEleup") )
    mMetEleup = iConfig.getParameter<edm::InputTag>("srcMetEleup");
  if(  iConfig.existsAs<edm::InputTag>("srcMetEledown") )
    mMetEledown = iConfig.getParameter<edm::InputTag>("srcMetEledown");
  if(  iConfig.existsAs<edm::InputTag>("srcMetJetup") )
    mMetJetup = iConfig.getParameter<edm::InputTag>("srcMetJetup");
  if(  iConfig.existsAs<edm::InputTag>("srcMetJetdown") )
    mMetJetdown = iConfig.getParameter<edm::InputTag>("srcMetJetdown");
  if(  iConfig.existsAs<edm::InputTag>("srcMetUncup") )
    mMetUncup = iConfig.getParameter<edm::InputTag>("srcMetUncup");
  if(  iConfig.existsAs<edm::InputTag>("srcMetUncdown") )
    mMetUncdown = iConfig.getParameter<edm::InputTag>("srcMetUncdown");
//###################### MVA MET #####################
  if(  iConfig.existsAs<edm::InputTag>("srcMVACorrected") )
    mMetMVA  = iConfig.getParameter<edm::InputTag>("srcMVACorrected");
  if(  iConfig.existsAs<edm::InputTag>("srcMVAJERup") )
    mMVAJERup = iConfig.getParameter<edm::InputTag>("srcMVAJERup");
  if(  iConfig.existsAs<edm::InputTag>("srcMVAJERdown") )
    mMVAJERdown = iConfig.getParameter<edm::InputTag>("srcMVAJERdown");
  if(  iConfig.existsAs<edm::InputTag>("srcMVAEleup") )
    mMVAEleup = iConfig.getParameter<edm::InputTag>("srcMVAEleup");
  if(  iConfig.existsAs<edm::InputTag>("srcMVAEledown") )
    mMVAEledown = iConfig.getParameter<edm::InputTag>("srcMVAEledown");
  if(  iConfig.existsAs<edm::InputTag>("srcMVAJetup") )
    mMVAJetup = iConfig.getParameter<edm::InputTag>("srcMVAJetup");
  if(  iConfig.existsAs<edm::InputTag>("srcMVAJetdown") ) 
    mMVAJetdown = iConfig.getParameter<edm::InputTag>("srcMVAJetdown");
  if(  iConfig.existsAs<edm::InputTag>("srcMVAUncup") )
    mMVAUncup = iConfig.getParameter<edm::InputTag>("srcMVAUncup");
  if(  iConfig.existsAs<edm::InputTag>("srcMVAUncdown") )
    mMVAUncdown = iConfig.getParameter<edm::InputTag>("srcMVAUncdown");
//###################### MVA MET Unity #####################
  if(  iConfig.existsAs<edm::InputTag>("srcMVAUnityCorrected") )
    mMetMVAUnity  = iConfig.getParameter<edm::InputTag>("srcMVAUnityCorrected");
  if(  iConfig.existsAs<edm::InputTag>("srcMVAUnityJERup") )
    mMVAUnityJERup = iConfig.getParameter<edm::InputTag>("srcMVAUnityJERup");
  if(  iConfig.existsAs<edm::InputTag>("srcMVAUnityJERdown") )
    mMVAUnityJERdown = iConfig.getParameter<edm::InputTag>("srcMVAUnityJERdown");
  if(  iConfig.existsAs<edm::InputTag>("srcMVAUnityEleup") )
    mMVAUnityEleup = iConfig.getParameter<edm::InputTag>("srcMVAUnityEleup");
  if(  iConfig.existsAs<edm::InputTag>("srcMVAUnityEledown") )
    mMVAUnityEledown = iConfig.getParameter<edm::InputTag>("srcMVAUnityEledown");
  if(  iConfig.existsAs<edm::InputTag>("srcMVAUnityJetup") )
    mMVAUnityJetup = iConfig.getParameter<edm::InputTag>("srcMVAUnityJetup");
  if(  iConfig.existsAs<edm::InputTag>("srcMVAUnityJetdown") )
    mMVAUnityJetdown = iConfig.getParameter<edm::InputTag>("srcMVAUnityJetdown");
  if(  iConfig.existsAs<edm::InputTag>("srcMVAUnityUncup") )
    mMVAUnityUncup = iConfig.getParameter<edm::InputTag>("srcMVAUnityUncup");
  if(  iConfig.existsAs<edm::InputTag>("srcMVAUnityUncdown") )
    mMVAUnityUncdown = iConfig.getParameter<edm::InputTag>("srcMVAUnityUncdown");
*/
//###################### noPU MET #####################
  if(  iConfig.existsAs<edm::InputTag>("srcnPUCorrected") )
    mMetnoPU = iConfig.getParameter<edm::InputTag>("srcnPUCorrected");
  if(  iConfig.existsAs<edm::InputTag>("srcnPUJERup") )
    noPUJERup = iConfig.getParameter<edm::InputTag>("srcnPUJERup");
  if(  iConfig.existsAs<edm::InputTag>("srcnPUJERdown") )
    noPUJERdown = iConfig.getParameter<edm::InputTag>("srcnPUJERdown");
  if(  iConfig.existsAs<edm::InputTag>("srcnPUEleup") )
    noPUEleup = iConfig.getParameter<edm::InputTag>("srcnPUEleup"); 
  if(  iConfig.existsAs<edm::InputTag>("srcnPUEledown") )
    noPUEledown = iConfig.getParameter<edm::InputTag>("srcnPUEledown");
  if(  iConfig.existsAs<edm::InputTag>("srcnPUJetup") )
    noPUJetup = iConfig.getParameter<edm::InputTag>("srcnPUJetup");
  if(  iConfig.existsAs<edm::InputTag>("srcnPUJetdown") )
    noPUJetdown = iConfig.getParameter<edm::InputTag>("srcnPUJetdown");
  if(  iConfig.existsAs<edm::InputTag>("srcnPUUncup") )
    noPUUncup = iConfig.getParameter<edm::InputTag>("srcnPUUncup");
  if(  iConfig.existsAs<edm::InputTag>("srcnPUUncdown") )
    noPUUncdown = iConfig.getParameter<edm::InputTag>("srcnPUUncdown");


  //*********************  Run Over AOD or PAT  ***********//
  if( iConfig.existsAs<bool>("runningOverAOD"))
	  runoverAOD = iConfig.getParameter<bool>("runningOverAOD");

}

 

ewk::VplusJetsAnalysis::~VplusJetsAnalysis() {}


void ewk::VplusJetsAnalysis::beginJob() {

  // Declare all the branches of the tree
  declareTreeBranches();
  if(runningOverMC_)
  {
  	LumiReweightABCD_   =  edm::LumiReWeighting("PUHistS10.root","PUHistRun2012All.root","pileup","pileup");  
  	LumiReweightABCD3d_ =  edm::Lumi3DReWeighting("PUHistS10.root","PUHistRun2012All.root", "pileup", "pileup","Weight3D_ABC.root");
  	LumiReweightABCD3d_.weight3D_init(69.3/69.3);
  //LumiReweightABCDbug_ = edm::LumiReWeighting("PUHistS10.root","PUHistRun2012AllBug.root","pileup","pileup");
  }
}




// ------------ method called to produce the data  ------------
void ewk::VplusJetsAnalysis::analyze(const edm::Event& iEvent, 
				     const edm::EventSetup& iSetup) {

  // write event information
  run   = 0;
  event = 0;
  lumi  = 0;
  bunch = 0;
  nPV   = 0; 
  mNVB  = 0; 
 
  // run, event, bunch crossing, ....
  run   = iEvent.id().run();
  event = iEvent.id().event();
  lumi  = iEvent.luminosityBlock();
  bunch = iEvent.bunchCrossing();


  // primary/secondary vertices
  // edm::Handle<reco::VertexCollection > recVtxs;
  edm::Handle <edm::View<reco::Vertex> > recVtxs;
  iEvent.getByLabel( mPrimaryVertex, recVtxs);
  for(unsigned int ind=0;ind<recVtxs->size();ind++) 
    {
      if(nPV>=60) continue;
      mPVx[nPV] =   -10000.0;
      mPVy[nPV] =   -10000.0;
      mPVz[nPV] =   -10000.0;

      if (!((*recVtxs)[ind].isFake()) && ((*recVtxs)[ind].ndof()>=4) 
	  && (fabs((*recVtxs)[ind].z())<=24.0) &&  
	  ((*recVtxs)[ind].position().Rho()<=2.0) ) {

	mPVx[nPV] =  (*recVtxs)[ind].x() ;
	mPVy[nPV] =  (*recVtxs)[ind].y() ;
	mPVz[nPV] =  (*recVtxs)[ind].z() ;
	nPV += 1;
      }
    }

  //////////// Beam spot //////////////
//  if(runOverAOD){
  edm::Handle<reco::BeamSpot >beamSpot;
  if(runoverAOD){
	  iEvent.getByLabel(mInputBeamSpot, beamSpot);
	  mBSx = beamSpot->position().X();
	  mBSy = beamSpot->position().Y();
	  mBSz = beamSpot->position().Z();
  }

  /////// PfMET information /////
  edm::Handle<edm::View<pat::MET> > pfmet;
  iEvent.getByLabel(mInputMet, pfmet);
  if (pfmet->size() == 0) {    mpfMET   = -1;    mpfSumET = -1;    mpfMETSign = -1;    mpfMETPhi   = -10.0;  }
  else {    mpfMET   = (*pfmet)[0].et();    mpfSumET = (*pfmet)[0].sumEt();    mpfMETSign = (*pfmet)[0].significance();    mpfMETPhi   = (*pfmet)[0].phi();  }
  
  //############# CaloMET ###############
  edm::Handle<edm::View<pat::MET> > calomet;
  iEvent.getByLabel(mMetCalo, calomet);
  if (calomet->size() == 0) {    mcaloMET   = -1;    mcaloSumET = -1;    mcaloMETSign = -1;    mcaloMETPhi   = -10.0;  }
  else {    mcaloMET   = (*calomet)[0].et();    mcaloSumET = (*calomet)[0].sumEt();    mcaloMETSign = (*calomet)[0].significance();    mcaloMETPhi   = (*calomet)[0].phi();  }

  edm::Handle<edm::View<pat::MET> > caloT1met;
  iEvent.getByLabel(mMetCaloT1, caloT1met);
  if (caloT1met->size() == 0) {    mcaloT1MET   = -1;    mcaloT1SumET = -1;    mcaloT1METSign = -1;    mcaloT1METPhi   = -10.0;  }
  else {    mcaloT1MET   = (*caloT1met)[0].et();    mcaloT1SumET = (*caloT1met)[0].sumEt();    mcaloT1METSign = (*caloT1met)[0].significance();    mcaloT1METPhi   = (*caloT1met)[0].phi();  }
  /*
  //############# PFMET uncertainties ###############
  edm::Handle<edm::View<pat::MET> > metraw;
  iEvent.getByLabel(mMetRaw, metraw);
  if (metraw->size() == 0) {    mRawMET   = -1;    mRawSumET = -1;    mRawMETSign = -1;    mRawMETPhi   = -10.0;  }
  else {    mRawMET   = (*metraw)[0].et();    mRawSumET = (*metraw)[0].sumEt();    mRawMETSign = (*metraw)[0].significance();    mRawMETPhi   = (*metraw)[0].phi();  }

  edm::Handle<edm::View<pat::MET> > met01;
  iEvent.getByLabel(mMetType01, met01);
  if (met01->size() == 0) {    mT01MET   = -1;    mT01SumET = -1;    mT01METSign = -1;    mT01METPhi   = -10.0;  }
  else {    mT01MET   = (*met01)[0].et();    mT01SumET = (*met01)[0].sumEt();    mT01METSign = (*met01)[0].significance();    mT01METPhi   = (*met01)[0].phi();   }
 
  edm::Handle<edm::View<pat::MET> > metMVA;
  iEvent.getByLabel(mMetMVA, metMVA);
  if (metMVA->size() == 0) {    mMVAMET   = -1;    mMVASumET = -1;    mMVAMETSign = -1;    mMVAMETPhi   = -10.0; }
  else {    mMVAMET   = (*metMVA)[0].et();    mMVASumET = (*metMVA)[0].sumEt();    mMVAMETSign = (*metMVA)[0].significance();    mMVAMETPhi   = (*metMVA)[0].phi(); }
  */
  edm::Handle<edm::View<pat::MET> > metPU;
  iEvent.getByLabel(mMetnoPU, metPU);
  if (metPU->size() == 0) {    noPUMET   = -1;    noPUSumET = -1;    noPUMETSign = -1;    noPUMETPhi   = -10.0; }
  else {    noPUMET   = (*metPU)[0].et();    noPUSumET = (*metPU)[0].sumEt();    noPUMETSign = (*metPU)[0].significance();    noPUMETPhi   = (*metPU)[0].phi(); }
  /*
  // MVA Unity
  edm::Handle<edm::View<pat::MET> > metUni;
  iEvent.getByLabel(mMetMVAUnity, metUni);
  if (metUni->size() == 0) {    mUniMET   = -1;    mUniSumET = -1;    mUniMETSign = -1;    mUniMETPhi   = -10.0; }
  else {    mUniMET   = (*metUni)[0].et();    mUniSumET = (*metUni)[0].sumEt();    mUniMETSign = (*metUni)[0].significance();    mUniMETPhi   = (*metUni)[0].phi(); }
  //std::cout<<"MVA Unity central = "<<mUniMET<<std::endl;
  */
  //mJERup = -1;    mJERupSumET = -1;       mJERupPhi = -1;
  //mJERdown = -1;  mJERdownSumET = -1;     mJERdownPhi = -1;
  //mvaJERup = -1;    mvaJERupSumET = -1;       mvaJERupPhi = -1;
  //mvaJERdown = -1;  mvaJERdownSumET = -1;     mvaJERdownPhi = -1;
  puJERup = -1;    puJERupSumET = -1;       puJERupPhi = -1;
  puJERdown = -1;  puJERdownSumET = -1;     puJERdownPhi = -1;
  // MVA Unity
  //uniJERup = -1;    uniJERupSumET = -1;       uniJERupPhi = -1;
  //uniJERdown = -1;  uniJERdownSumET = -1;     uniJERdownPhi = -1;
  
  if( runningOverMC_ ) {
/*          edm::Handle<pat::METCollection> metJERup;
          iEvent.getByLabel(mMetJERup, metJERup);
          if (metJERup->size() == 0)    
	  {	mJERup = -1;	mJERupSumET = -1;	mJERupPhi = -1;}
          else {         
		mJERup = (*metJERup)[0].et();	mJERupSumET = (*metJERup)[0].sumEt();	mJERupPhi = (*metJERup)[0].phi();}
          edm::Handle<pat::METCollection> metJERdown;
          iEvent.getByLabel(mMetJERdown, metJERdown);
          if (metJERdown->size() == 0)    
	  {	mJERdown = -1;	mJERdownSumET = -1;	mJERdownPhi = -1;}
          else {         
		mJERdown = (*metJERdown)[0].et();	mJERdownSumET = (*metJERdown)[0].sumEt();	mJERdownPhi = (*metJERdown)[0].phi();}

	  edm::Handle<pat::METCollection> metMVAJERup;
          iEvent.getByLabel(mMVAJERup, metMVAJERup);
          if (metMVAJERup->size() == 0)
	    {     mvaJERup = -1;    mvaJERupSumET = -1;       mvaJERupPhi = -1;}
          else {
	    mvaJERup = (*metMVAJERup)[0].et();   mvaJERupSumET = (*metMVAJERup)[0].sumEt();   mvaJERupPhi = (*metMVAJERup)[0].phi();}
	  edm::Handle<pat::METCollection> metMVAJERdown;
          iEvent.getByLabel(mMVAJERdown, metMVAJERdown);
          if (metMVAJERdown->size() == 0)
	    {     mvaJERdown = -1;  mvaJERdownSumET = -1;     mvaJERdownPhi = -1;}
          else {
	    mvaJERdown = (*metMVAJERdown)[0].et();       mvaJERdownSumET = (*metMVAJERdown)[0].sumEt();       mvaJERdownPhi = (*metMVAJERdown)[0].phi();}
*/
	  edm::Handle<pat::METCollection> metPUJERup;
          iEvent.getByLabel(noPUJERup, metPUJERup);
          if (metPUJERup->size() == 0)
	    {     puJERup = -1;    puJERupSumET = -1;       puJERupPhi = -1;}
          else {
	    puJERup = (*metPUJERup)[0].et();   puJERupSumET = (*metPUJERup)[0].sumEt();   puJERupPhi = (*metPUJERup)[0].phi();}
	  edm::Handle<pat::METCollection> metPUJERdown;
          iEvent.getByLabel(noPUJERdown, metPUJERdown);
          if (metPUJERdown->size() == 0)
	    {     puJERdown = -1;  puJERdownSumET = -1;     puJERdownPhi = -1;}
          else {
	    puJERdown = (*metPUJERdown)[0].et();       puJERdownSumET = (*metPUJERdown)[0].sumEt();       puJERdownPhi = (*metPUJERdown)[0].phi();}
/*
	  // MVA Unity
          edm::Handle<pat::METCollection> metUniJERup;
          iEvent.getByLabel(mMVAUnityJERup, metUniJERup);
          if (metUniJERup->size() == 0)
            {     uniJERup = -1;    uniJERupSumET = -1;       uniJERupPhi = -1;}
          else {
            uniJERup = (*metUniJERup)[0].et();   uniJERupSumET = (*metUniJERup)[0].sumEt();   uniJERupPhi = (*metUniJERup)[0].phi();}
          edm::Handle<pat::METCollection> metUniJERdown;
          iEvent.getByLabel(mMVAUnityJERdown, metUniJERdown);
          if (metUniJERdown->size() == 0)
            {     uniJERdown = -1;  uniJERdownSumET = -1;     uniJERdownPhi = -1;}
          else {
            uniJERdown = (*metUniJERdown)[0].et();       uniJERdownSumET = (*metUniJERdown)[0].sumEt();       uniJERdownPhi = (*metUniJERdown)[0].phi();}

	  //std::cout<<"MVA Unity JERup   = "<<uniJERup<<std::endl;
	  //std::cout<<"MVA Unity JERdown = "<<uniJERdown<<std::endl;
	  //std::cout<<"#################################"<<std::endl;
*/
  }
/*
  edm::Handle<pat::METCollection> metEleup;
  iEvent.getByLabel(mMetEleup, metEleup);
  if (metEleup->size() == 0)    
    	{ mEleup = -1;	mEleupSumET = -1; 	mEleupPhi = -1;	}
  else         
	{ mEleup = (*metEleup)[0].et();		mEleupSumET = (*metEleup)[0].sumEt();	mEleupPhi = (*metEleup)[0].phi();	}

  edm::Handle<pat::METCollection> metEledown;
  iEvent.getByLabel(mMetEledown, metEledown);
  if (metEledown->size() == 0)    
	{ mEledown = -1;	mEledownSumET = -1; 	mEledownPhi = -1; }
  else         
	{ mEledown = (*metEledown)[0].et();	mEledownSumET = (*metEledown)[0].sumEt(); 	mEledownPhi = (*metEledown)[0].phi();       }

  edm::Handle<pat::METCollection> metMVAEleup;
  iEvent.getByLabel(mMVAEleup, metMVAEleup);
  if (metMVAEleup->size() == 0)    
	{ mvaEleup = -1;	mvaEleupSumET = -1; 	mvaEleupPhi = -1; }
  else          
	{ mvaEleup = (*metMVAEleup)[0].et();	mvaEleupSumET = (*metMVAEleup)[0].sumEt(); 	mvaEleupPhi = (*metMVAEleup)[0].phi();       }
  edm::Handle<pat::METCollection> metMVAEledown;
  iEvent.getByLabel(mMVAEledown, metMVAEledown);
  if (metMVAEledown->size() == 0)    
	{ mvaEledown = -1;	mvaEledownSumET = -1; 	mvaEledownPhi = -1; }
  else          
	{ mvaEledown = (*metMVAEledown)[0].et();	mvaEledownSumET = (*metMVAEledown)[0].sumEt();		mvaEledownPhi = (*metMVAEledown)[0].phi();       }
*/
  edm::Handle<pat::METCollection> metPUEleup;
  iEvent.getByLabel(noPUEleup, metPUEleup);
  if (metPUEleup->size() == 0)    
	{ puEleup = -1;		puEleupSumET = -1;	puEleupPhi = -1; }
  else          
	{ puEleup = (*metPUEleup)[0].et();	puEleupSumET = (*metPUEleup)[0].sumEt(); 	puEleupPhi = (*metPUEleup)[0].phi();       }
  edm::Handle<pat::METCollection> metPUEledown;
  iEvent.getByLabel(noPUEledown, metPUEledown);
  if (metPUEledown->size() == 0)    
	{ puEledown = -1;	puEledownSumET = -1;	puEledownPhi = -1; }
  else          
	{ puEledown = (*metPUEledown)[0].et();	puEledownSumET = (*metPUEledown)[0].sumEt(); 	puEledownPhi = (*metPUEledown)[0].phi();       }
/*
  // MVA Unity
  edm::Handle<pat::METCollection> metUniEleup;
  iEvent.getByLabel(mMVAUnityEleup, metUniEleup);
  if (metUniEleup->size() == 0)
        { uniEleup = -1;        uniEleupSumET = -1;     uniEleupPhi = -1; }
  else
        { uniEleup = (*metUniEleup)[0].et();    uniEleupSumET = (*metUniEleup)[0].sumEt();      uniEleupPhi = (*metUniEleup)[0].phi();       }
  edm::Handle<pat::METCollection> metUniEledown;
  iEvent.getByLabel(mMVAUnityEledown, metUniEledown);
  if (metUniEledown->size() == 0)
        { uniEledown = -1;      uniEledownSumET = -1;   uniEledownPhi = -1; }
  else
        { uniEledown = (*metUniEledown)[0].et();        uniEledownSumET = (*metUniEledown)[0].sumEt();          uniEledownPhi = (*metUniEledown)[0].phi();       }
  /////////////////////////////////////////////////
  edm::Handle<pat::METCollection> metJetup;
  iEvent.getByLabel(mMetJetup, metJetup);
  if (metJetup->size() == 0)
  {     mJetup = -1;        mJetupSumET = -1;        mJetupPhi = -1;}
  else  
  {     mJetup = (*metJetup)[0].et();        mJetupSumET = (*metJetup)[0].sumEt();        mJetupPhi = (*metJetup)[0].phi();}
  edm::Handle<pat::METCollection> metJetdown;
  iEvent.getByLabel(mMetJetdown, metJetdown);
  if (metJetdown->size() == 0)
  {     mJetdown = -1;        mJetdownSumET = -1;        mJetdownPhi = -1;}
  else          
  {	mJetdown = (*metJetdown)[0].et();        mJetdownSumET = (*metJetdown)[0].sumEt();        mJetdownPhi = (*metJetdown)[0].phi();  }

  edm::Handle<pat::METCollection> metMVAJetup;
  iEvent.getByLabel(mMVAJetup, metMVAJetup);
  if (metMVAJetup->size() == 0)
  {     mvaJetup = -1;        mvaJetupSumET = -1;        mvaJetupPhi = -1;}
  else
  {     mvaJetup = (*metMVAJetup)[0].et();        mvaJetupSumET = (*metMVAJetup)[0].sumEt();        mvaJetupPhi = (*metMVAJetup)[0].phi();}
  edm::Handle<pat::METCollection> metMVAJetdown;
  iEvent.getByLabel(mMVAJetdown, metMVAJetdown);
  if (metMVAJetdown->size() == 0)
  {     mvaJetdown = -1;        mvaJetdownSumET = -1;        mvaJetdownPhi = -1;}
  else
  {     mvaJetdown = (*metMVAJetdown)[0].et();        mvaJetdownSumET = (*metMVAJetdown)[0].sumEt();        mvaJetdownPhi = (*metMVAJetdown)[0].phi();  }
*/
  edm::Handle<pat::METCollection> metPUJetup;
  iEvent.getByLabel(noPUJetup, metPUJetup);
  if (metPUJetup->size() == 0)
  {     puJetup = -1;        puJetupSumET = -1;        puJetupPhi = -1;}
  else
  {     puJetup = (*metPUJetup)[0].et();        puJetupSumET = (*metPUJetup)[0].sumEt();        puJetupPhi = (*metPUJetup)[0].phi();}
  edm::Handle<pat::METCollection> metPUJetdown;
  iEvent.getByLabel(noPUJetdown, metPUJetdown);
  if (metPUJetdown->size() == 0)
  {     puJetdown = -1;        puJetdownSumET = -1;        puJetdownPhi = -1;}
  else
  {     puJetdown = (*metPUJetdown)[0].et();        puJetdownSumET = (*metPUJetdown)[0].sumEt();        puJetdownPhi = (*metPUJetdown)[0].phi();  }
/*
  // MVA Unity
  edm::Handle<pat::METCollection> metUniJetup;
  iEvent.getByLabel(mMVAUnityJetup, metUniJetup);
  if (metUniJetup->size() == 0)
  {     uniJetup = -1;        uniJetupSumET = -1;        uniJetupPhi = -1;}
  else
  {     uniJetup = (*metUniJetup)[0].et();        uniJetupSumET = (*metUniJetup)[0].sumEt();        uniJetupPhi = (*metUniJetup)[0].phi();}
  edm::Handle<pat::METCollection> metUniJetdown;
  iEvent.getByLabel(mMVAUnityJetdown, metUniJetdown);
  if (metUniJetdown->size() == 0)
  {     uniJetdown = -1;        uniJetdownSumET = -1;        uniJetdownPhi = -1;}
  else
  {     uniJetdown = (*metUniJetdown)[0].et();        uniJetdownSumET = (*metUniJetdown)[0].sumEt();        uniJetdownPhi = (*metUniJetdown)[0].phi();  }

  /////////////////////////////////////////////////
  edm::Handle<pat::METCollection> metUncup;
  iEvent.getByLabel(mMetUncup, metUncup);
  if (metUncup->size() == 0)
  {     mUncup = -1;        mUncupSumET = -1;        mUncupPhi = -1;}
  else          
  {	mUncup = (*metUncup)[0].et();        mUncupSumET = (*metUncup)[0].sumEt();        mUncupPhi = (*metUncup)[0].phi();  }
  edm::Handle<pat::METCollection> metUncdown;
  iEvent.getByLabel(mMetUncdown, metUncdown);
  if (metUncdown->size() == 0)    
  {	mUncdown = -1;        mUncdownSumET = -1;        mUncdownPhi = -1;}
  else          
  {	mUncdown = (*metUncdown)[0].et();        mUncdownSumET = (*metUncdown)[0].sumEt();        mUncdownPhi = (*metUncdown)[0].phi();  }

  edm::Handle<pat::METCollection> metMVAUncup;
  iEvent.getByLabel(mMVAUncup, metMVAUncup);
  if (metMVAUncup->size() == 0)
  {     mvaUncup = -1;        mvaUncupSumET = -1;        mvaUncupPhi = -1;}
  else
  {     mvaUncup = (*metMVAUncup)[0].et();        mvaUncupSumET = (*metMVAUncup)[0].sumEt();        mvaUncupPhi = (*metMVAUncup)[0].phi();  }
  edm::Handle<pat::METCollection> metMVAUncdown;
  iEvent.getByLabel(mMVAUncdown, metMVAUncdown);
  if (metMVAUncdown->size() == 0)
  {     mvaUncdown = -1;        mvaUncdownSumET = -1;        mvaUncdownPhi = -1;}
  else
  {     mvaUncdown = (*metMVAUncdown)[0].et();        mvaUncdownSumET = (*metMVAUncdown)[0].sumEt();        mvaUncdownPhi = (*metMVAUncdown)[0].phi();  }
*/
  edm::Handle<pat::METCollection> metPUUncup;
  iEvent.getByLabel(noPUUncup, metPUUncup);
  if (metPUUncup->size() == 0)
  {     puUncup = -1;        puUncupSumET = -1;        puUncupPhi = -1;}
  else
  {     puUncup = (*metPUUncup)[0].et();        puUncupSumET = (*metPUUncup)[0].sumEt();        puUncupPhi = (*metPUUncup)[0].phi();  }
  edm::Handle<pat::METCollection> metPUUncdown;
  iEvent.getByLabel(noPUUncdown, metPUUncdown);
  if (metPUUncdown->size() == 0)
  {     puUncdown = -1;        puUncdownSumET = -1;        puUncdownPhi = -1;}
  else
  {     puUncdown = (*metPUUncdown)[0].et();        puUncdownSumET = (*metPUUncdown)[0].sumEt();        puUncdownPhi = (*metPUUncdown)[0].phi();  }
/*
  // MVA Unity
  edm::Handle<pat::METCollection> metUniUncup;
  iEvent.getByLabel(mMVAUnityUncup, metUniUncup);
  if (metUniUncup->size() == 0)
  {     uniUncup = -1;        uniUncupSumET = -1;        uniUncupPhi = -1;}
  else
  {     uniUncup = (*metUniUncup)[0].et();        uniUncupSumET = (*metUniUncup)[0].sumEt();        uniUncupPhi = (*metUniUncup)[0].phi();  }
  edm::Handle<pat::METCollection> metUniUncdown;
  iEvent.getByLabel(mMVAUnityUncdown, metUniUncdown);
  if (metUniUncdown->size() == 0)
  {     uniUncdown = -1;        uniUncdownSumET = -1;        uniUncdownPhi = -1;}
  else
  {     uniUncdown = (*metUniUncdown)[0].et();        uniUncdownSumET = (*metUniUncdown)[0].sumEt();        uniUncdownPhi = (*metUniUncdown)[0].phi();  }
*/
  //#############################################


  /////////// GenMET information & MC Pileup Summary Info  //////////
//  mcPUtrueInteractions = 0;
  mcPUtotnvtx = 0;
  for (int i=0; i<13; i++) {
  	mcPUbx[i]   = -999;
  	mcPUnvtx[i] = -999;
  } 
  mcPUnvtx_bxPreviuos = -999;
  mcPUnvtx_inTime = -999;
  mcPUnvtx_bxNext = -999;
  WeightABCD = 1.;  WeightABCD3d = 1.; //WeightABCDbug = 1.;
  if ( runningOverMC_ ){
    // MC Pileup Summary Info
    const edm::InputTag PileupSrc("addPileupInfo");
    edm::Handle<std::vector< PileupSummaryInfo > >  PupInfo;
    iEvent.getByLabel(PileupSrc, PupInfo);
    std::vector<PileupSummaryInfo>::const_iterator PVI;
    int ctid = 0;
    int bxPrevious_ = -1;
    int bxNext_ = +1;
    for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
      //if (ctid>2) break;
      mcPUbx[ctid]   =  PVI->getBunchCrossing();
      //std::cout << "Bunch crossing = "<<mcPUbx[ctid]<< std::endl;
      mcPUnvtx[ctid] =  PVI->getPU_NumInteractions();
      mcPUtotnvtx   +=  PVI->getPU_NumInteractions();
      if(PVI->getBunchCrossing() == bxPrevious_) mcPUnvtx_bxPreviuos = PVI->getTrueNumInteractions();
      if(PVI->getBunchCrossing() == 0)
        {       //mcPUtrueInteractions = PVI->getTrueNumInteractions();
                mcPUnvtx_inTime = PVI->getTrueNumInteractions();
        }
      if(PVI->getBunchCrossing() == bxNext_) mcPUnvtx_bxNext = PVI->getTrueNumInteractions();
      ctid++;
    }
        WeightABCD   = LumiReweightABCD_.weight( mcPUnvtx_inTime );
        WeightABCD3d = LumiReweightABCD3d_.weight3D( mcPUnvtx_bxPreviuos, mcPUnvtx_inTime, mcPUnvtx_bxNext ); 
	//WeightABCDbug = LumiReweightABCDbug_.weight( mcPUnvtx_inTime );
    }
  // fill jet branches
  edm::Handle<edm::View< reco::Candidate> > boson;
  iEvent.getByLabel( mInputBoson, boson);
  mNVB = boson->size();

  if( mNVB < 1 ) return; // Nothing to fill

  if(CorrectedPatJetFiller.get())   	CorrectedPatJetFiller->fill(iEvent);
  if(CorrectedSelectedJetFiller.get()) 	CorrectedSelectedJetFiller->fill(iEvent);
  //if(SmearedJetFiller.get()) 	    	SmearedJetFiller->fill(iEvent);
  //if(ElectronFiller.get())	 	ElectronFiller->fill(iEvent);

  /**  Store reconstructed vector boson information */
  recoBosonFillerE->fill(iEvent, 0);

  /**  Store generated vector boson information */
  //if(genBosonFiller.get()) genBosonFiller->fill(iEvent);
  
  myTree->Fill();

} // analyze method






void ewk::VplusJetsAnalysis::endJob()
{
  //hOutputFile->SetCompressionLevel(2);
  //hOutputFile->cd();
  //myTree->Write();

  //delete myTree;
  //hOutputFile->Close();
  //delete hOutputFile;
}




//  **** Utility: declare TTree branches for ntuple variables ***
void ewk::VplusJetsAnalysis::declareTreeBranches() {

  //myTree = fs->make<TTree>("ZJet","");

  myTree->Branch("event_runNo",  		&run,   		"event_runNo/I");
  myTree->Branch("event_evtNo",  		&event, 		"event_evtNo/I");
  myTree->Branch("event_lumi",   		&lumi,  		"event_lumi/I"); 
  myTree->Branch("event_bunch",  		&bunch, 		"event_bunch/I"); 
  myTree->Branch("event_nPV",    		&nPV,   		"event_nPV/I"); 
  myTree->Branch("event_PVx",    		mPVx,   		"event_PVx[30]/F"); 
  myTree->Branch("event_PVy",    		mPVy,   		"event_PVy[30]/F"); 
  myTree->Branch("event_PVz",    		mPVz,   		"event_PVz[30]/F");
  myTree->Branch("event_met_pfmet",             &mpfMET,      		"event_met_pfmet/F"); 
  myTree->Branch("event_met_pfsumet",           &mpfSumET,    		"event_met_pfsumet/F"); 
  myTree->Branch("event_met_pfmetPhi",          &mpfMETPhi,   		"event_met_pfmetPhi/F"); 
  myTree->Branch("event_met_calomet",           &mcaloMET,      	"event_met_calomet/F");
  myTree->Branch("event_met_calosumet",         &mcaloSumET,    	"event_met_calosumet/F");
  myTree->Branch("event_met_calometPhi",        &mcaloMETPhi,   	"event_met_calometPhi/F");
  myTree->Branch("event_met_caloT1met",         &mcaloT1MET,    	"event_met_caloT1met/F");
  myTree->Branch("event_met_caloT1sumet",       &mcaloT1SumET,  	"event_met_caloT1sumet/F");
  myTree->Branch("event_met_caloT1metPhi",      &mcaloT1METPhi, 	"event_met_caloT1metPhi/F");
//############################################	
/*  myTree->Branch("event_met_Rawmet",            &mRawMET,      		"event_met_Rawmet/F");
  myTree->Branch("event_met_Rawsumet",          &mRawSumET,    		"event_met_Rawsumet/F");
  myTree->Branch("event_met_RawmetPhi",         &mRawMETPhi,   		"event_met_RawmetPhi/F");
  myTree->Branch("event_met_T01met",            &mT01MET,      		"event_met_T01met/F");
  myTree->Branch("event_met_T01sumet",          &mT01SumET,    		"event_met_T01sumet/F");
  myTree->Branch("event_met_T01metPhi",         &mT01METPhi,   		"event_met_T01metPhi/F");
  myTree->Branch("event_met_MVAmet",            &mMVAMET,      		"event_met_MVAmet/F");
  myTree->Branch("event_met_MVAsumet",          &mMVASumET,    		"event_met_MVAsumet/F");
  myTree->Branch("event_met_MVAmetPhi",         &mMVAMETPhi,   		"event_met_MVAmetPhi/F");
*/  myTree->Branch("event_met_noPUmet",           &noPUMET,      		"event_met_noPUmet/F");
  myTree->Branch("event_met_noPUsumet",         &noPUSumET,    		"event_met_noPUsumet/F");
  myTree->Branch("event_met_noPUmetPhi",        &noPUMETPhi,   		"event_met_noPUmetPhi/F");
  // MVA Unity
/*  myTree->Branch("event_met_Unimet",            &mUniMET,               "event_met_Unimet/F");
  myTree->Branch("event_met_Unisumet",          &mUniSumET,             "event_met_Unisumet/F");
  myTree->Branch("event_met_UnimetPhi",         &mUniMETPhi,            "event_met_UnimetPhi/F");

  myTree->Branch("event_met_JERup",       	&mJERup,       		"event_met_JERup/F");
  myTree->Branch("event_met_JERupsumet",  	&mJERupSumET,  		"event_met_JERupsumet/F");
  myTree->Branch("event_met_JERupPhi",    	&mJERupPhi,    		"event_met_JERupPhi/F");
  myTree->Branch("event_met_JERdown",     	&mJERdown,     		"event_met_JERdown/F");
  myTree->Branch("event_met_JERdownsumet",  	&mJERdownSumET, 	"event_met_JERdownsumet/F");
  myTree->Branch("event_met_JERdownPhi",    	&mJERdownPhi,   	"event_met_JERdownPhi/F");
  myTree->Branch("event_mva_JERup",       	&mvaJERup,       	"event_mva_JERup/F");
  myTree->Branch("event_mva_JERupsumet",  	&mvaJERupSumET,  	"event_mva_JERupsumet/F");
  myTree->Branch("event_mva_JERupPhi",    	&mvaJERupPhi,    	"event_mva_JERupPhi/F");
  myTree->Branch("event_mva_JERdown",     	&mvaJERdown,     	"event_mva_JERdown/F");
  myTree->Branch("event_mva_JERdownsumet",  	&mvaJERdownSumET,  	"event_mva_JERdownsumet/F");
  myTree->Branch("event_mva_JERdownPhi",    	&mvaJERdownPhi,    	"event_mva_JERdownPhi/F");
*/  myTree->Branch("event_pu_JERup",       	&puJERup,       	"event_pu_JERup/F");
  myTree->Branch("event_pu_JERupsumet",  	&puJERupSumET,  	"event_pu_JERupsumet/F");
  myTree->Branch("event_pu_JERupPhi",    	&puJERupPhi,    	"event_pu_JERupPhi/F");
  myTree->Branch("event_pu_JERdown",     	&puJERdown,     	"event_pu_JERdown/F");
  myTree->Branch("event_pu_JERdownsumet",  	&puJERdownSumET,  	"event_pu_JERdownsumet/F");
  myTree->Branch("event_pu_JERdownPhi",    	&puJERdownPhi,    	"event_pu_JERdownPhi/F");
  // MVA Unity
/*  myTree->Branch("event_uni_JERup",             &uniJERup,              "event_uni_JERup/F");
  myTree->Branch("event_uni_JERupsumet",        &uniJERupSumET,         "event_uni_JERupsumet/F");
  myTree->Branch("event_uni_JERupPhi",          &uniJERupPhi,           "event_uni_JERupPhi/F");
  myTree->Branch("event_uni_JERdown",           &uniJERdown,            "event_uni_JERdown/F");
  myTree->Branch("event_uni_JERdownsumet",      &uniJERdownSumET,       "event_uni_JERdownsumet/F");
  myTree->Branch("event_uni_JERdownPhi",        &uniJERdownPhi,         "event_uni_JERdownPhi/F");
*/
  myTree->Branch("event_pu_Eleup",              &puEleup,        	"event_pu_Eleup/F");
  myTree->Branch("event_pu_Eledown",            &puEledown,      	"event_pu_Eledown/F");
  myTree->Branch("event_pu_Eleupsumet",         &puEleupSumET,        	"event_pu_Eleupsumet/F");
  myTree->Branch("event_pu_Eledownsumet",       &puEledownSumET,      	"event_pu_Eledownsumet/F");
  myTree->Branch("event_pu_EleupPhi",           &puEleupPhi,        	"event_pu_EleupPhi/F");
  myTree->Branch("event_pu_EledownPhi",         &puEledownPhi,      	"event_pu_EledownPhi/F");
/*  myTree->Branch("event_mva_Eleup",             &mvaEleup,        	"event_mva_Eleup/F");
  myTree->Branch("event_mva_Eledown",           &mvaEledown,      	"event_mva_Eledown/F");
  myTree->Branch("event_mva_Eleupsumet",        &mvaEleupSumET,        	"event_mva_Eleupsumet/F");
  myTree->Branch("event_mva_Eledownsumet",      &mvaEledownSumET,      	"event_mva_Eledownsumet/F");
  myTree->Branch("event_mva_EleupPhi",          &mvaEleupPhi,        	"event_mva_EleupPhi/F");
  myTree->Branch("event_mva_EledownPhi",        &mvaEledownPhi,      	"event_mva_EledownPhi/F");
  myTree->Branch("event_met_Eleup",             &mEleup,        	"event_met_Eleup/F");
  myTree->Branch("event_met_Eledown",           &mEledown,      	"event_met_Eledown/F");
  myTree->Branch("event_met_Eleupsumet",        &mEleupSumET,        	"event_met_Eleupsumet/F");
  myTree->Branch("event_met_Eledownsumet",      &mEledownSumET,      	"event_met_Eledownsumet/F");
  myTree->Branch("event_met_EleupPhi",          &mEleupPhi,        	"event_met_EleupPhi/F");
  myTree->Branch("event_met_EledownPhi",        &mEledownPhi,      	"event_met_EledownPhi/F");
  // MVA Unity
  myTree->Branch("event_uni_Eleup",             &uniEleup,              "event_uni_Eleup/F");
  myTree->Branch("event_uni_Eledown",           &uniEledown,            "event_uni_Eledown/F");
  myTree->Branch("event_uni_Eleupsumet",        &uniEleupSumET,         "event_uni_Eleupsumet/F");
  myTree->Branch("event_uni_Eledownsumet",      &uniEledownSumET,       "event_uni_Eledownsumet/F");
  myTree->Branch("event_uni_EleupPhi",          &uniEleupPhi,           "event_uni_EleupPhi/F");
  myTree->Branch("event_uni_EledownPhi",        &uniEledownPhi,         "event_uni_EledownPhi/F");

  myTree->Branch("event_met_Jetup",             &mJetup,        	"event_met_Jetup/F");
  myTree->Branch("event_met_Jetupsumet",        &mJetupSumET,   	"event_met_Jetupsumet/F");
  myTree->Branch("event_met_JetupPhi",          &mJetupPhi,     	"event_met_JetupPhi/F");
  myTree->Branch("event_met_Jetdown",           &mJetdown,      	"event_met_Jetdown/F");
  myTree->Branch("event_met_Jetdownsumet",      &mJetdownSumET, 	"event_met_Jetdownsumet/F");
  myTree->Branch("event_met_JetdownPhi",        &mJetdownPhi,   	"event_met_JetdownPhi/F");
  myTree->Branch("event_mva_Jetup",             &mvaJetup,        	"event_mva_Jetup/F");
  myTree->Branch("event_mva_Jetupsumet",        &mvaJetupSumET,   	"event_mva_Jetupsumet/F");
  myTree->Branch("event_mva_JetupPhi",          &mvaJetupPhi,     	"event_mva_JetupPhi/F");
  myTree->Branch("event_mva_Jetdown",           &mvaJetdown,      	"event_mva_Jetdown/F");
  myTree->Branch("event_mva_Jetdownsumet",      &mvaJetdownSumET, 	"event_mva_Jetdownsumet/F");
  myTree->Branch("event_mva_JetdownPhi",        &mvaJetdownPhi,   	"event_mva_JetdownPhi/F");
*/  myTree->Branch("event_pu_Jetup",              &puJetup,        	"event_pu_Jetup/F");
  myTree->Branch("event_pu_Jetupsumet",         &puJetupSumET,   	"event_pu_Jetupsumet/F");
  myTree->Branch("event_pu_JetupPhi",           &puJetupPhi,     	"event_pu_JetupPhi/F");
  myTree->Branch("event_pu_Jetdown",            &puJetdown,      	"event_pu_Jetdown/F");
  myTree->Branch("event_pu_Jetdownsumet",       &puJetdownSumET, 	"event_pu_Jetdownsumet/F");
  myTree->Branch("event_pu_JetdownPhi",         &puJetdownPhi,   	"event_pu_JetdownPhi/F");
  // MVA Unity
/*  myTree->Branch("event_uni_Jetup",             &uniJetup,              "event_uni_Jetup/F");
  myTree->Branch("event_uni_Jetupsumet",        &uniJetupSumET,         "event_uni_Jetupsumet/F");
  myTree->Branch("event_uni_JetupPhi",          &uniJetupPhi,           "event_uni_JetupPhi/F");
  myTree->Branch("event_uni_Jetdown",           &uniJetdown,            "event_uni_Jetdown/F");
  myTree->Branch("event_uni_Jetdownsumet",      &uniJetdownSumET,       "event_uni_Jetdownsumet/F");
  myTree->Branch("event_uni_JetdownPhi",        &uniJetdownPhi,         "event_uni_JetdownPhi/F");

  myTree->Branch("event_met_Uncup",             &mUncup,        	"event_met_Uncup/F");
  myTree->Branch("event_met_Uncupsumet",        &mUncupSumET,   	"event_met_Uncupsumet/F");
  myTree->Branch("event_met_UncupPhi",          &mUncupPhi,     	"event_met_UncupPhi/F");
  myTree->Branch("event_met_Uncdown",           &mUncdown,      	"event_met_Uncdown/F");
  myTree->Branch("event_met_Uncdownsumet",      &mUncdownSumET, 	"event_met_Uncdownsumet/F");
  myTree->Branch("event_met_UncdownPhi",        &mUncdownPhi,   	"event_met_UncdownPhi/F");
  myTree->Branch("event_mva_Uncup",             &mvaUncup,        	"event_mva_Uncup/F");
  myTree->Branch("event_mva_Uncupsumet",        &mvaUncupSumET,   	"event_mva_Uncupsumet/F");
  myTree->Branch("event_mva_UncupPhi",          &mvaUncupPhi,     	"event_mva_UncupPhi/F");
  myTree->Branch("event_mva_Uncdown",           &mvaUncdown,      	"event_mva_Uncdown/F");
  myTree->Branch("event_mva_Uncdownsumet",      &mvaUncdownSumET, 	"event_mva_Uncdownsumet/F");
  myTree->Branch("event_mva_UncdownPhi",        &mvaUncdownPhi,   	"event_mva_UncdownPhi/F");
*/  myTree->Branch("event_pu_Uncup",              &puUncup,        	"event_pu_Uncup/F");
  myTree->Branch("event_pu_Uncupsumet",         &puUncupSumET,   	"event_pu_Uncupsumet/F");
  myTree->Branch("event_pu_UncupPhi",           &puUncupPhi,     	"event_pu_UncupPhi/F");
  myTree->Branch("event_pu_Uncdown",            &puUncdown,      	"event_pu_Uncdown/F");
  myTree->Branch("event_pu_Uncdownsumet",       &puUncdownSumET, 	"event_pu_Uncdownsumet/F");  
  myTree->Branch("event_pu_UncdownPhi",        	&puUncdownPhi,   	"event_pu_UncdownPhi/F");
  // MVA Unity
/*  myTree->Branch("event_uni_Uncup",             &uniUncup,              "event_uni_Uncup/F");
  myTree->Branch("event_uni_Uncupsumet",        &uniUncupSumET,         "event_uni_Uncupsumet/F");
  myTree->Branch("event_uni_UncupPhi",          &uniUncupPhi,           "event_uni_UncupPhi/F");
  myTree->Branch("event_uni_Uncdown",           &uniUncdown,            "event_uni_Uncdown/F");
  myTree->Branch("event_uni_Uncdownsumet",      &uniUncdownSumET,       "event_uni_Uncdownsumet/F");
  myTree->Branch("event_uni_UncdownPhi",        &uniUncdownPhi,         "event_uni_UncdownPhi/F");
*/
//############################################
  myTree->Branch(("num"+VBosonType_).c_str(),	&mNVB ,			("num"+VBosonType_+"/I").c_str());
  myTree->Branch("event_mcPU_totnvtx",    	&mcPUtotnvtx,  		"event_mcPU_totnvtx/F");
  //    myTree->Branch("event_mcPU_trueInteractions",    &mcPUtrueInteractions,  "event_mcPU_trueInteractions/F"); 
  myTree->Branch("event_mcPU_bx",         	mcPUbx ,       		"event_mcPU_bx[13]/F");
  myTree->Branch("event_mcPU_nvtx",       	mcPUnvtx,      		"event_mcPU_nvtx[13]/F");
  myTree->Branch("event_mcPU_nvtxPrevious", 	&mcPUnvtx_bxPreviuos, 	"event_mcPU_nvtxPrevious/I");
  myTree->Branch("event_mcPU_nvtxIntime", 	&mcPUnvtx_inTime, 	"event_mcPU_nvtxIntime/I");
  myTree->Branch("event_mcPU_nvtxNext", 	&mcPUnvtx_bxNext, 	"event_mcPU_nvtxNext/I");
  myTree->Branch("event_weightABCD",      	&WeightABCD,        	"event_weightABCD/D");
  myTree->Branch("event_weightABCD3D",      	&WeightABCD3d,        	"event_weightABCD3D/D");
  //myTree->Branch("event_weightABCDbug",		&WeightABCDbug,		"event_weightABCDbug/D");

}  





// declare this class as a plugin
#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
using ewk::VplusJetsAnalysis;
DEFINE_FWK_MODULE(VplusJetsAnalysis);
