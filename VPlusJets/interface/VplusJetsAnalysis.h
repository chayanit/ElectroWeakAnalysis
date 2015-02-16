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
 *   To fill W/Z + jets related quantities into a specified TTree
 *   Can work with CaloJet, GenJet, JPT jet, PF jet.
 *   Can work with jets in RECO/AOD/PAT data formats.
 * History:
 *   
 *
 * Copyright (C) 2010 FNAL 
 *****************************************************************************/

#ifndef VplusJetsAnalysis_h
#define VplusJetsAnalysis_h

// system include files
#include <memory>
#include <string>
#include <iostream>
#include <map>
#include <fstream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h" 

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TFile.h"
#include "TTree.h"

#include "ElectroWeakAnalysis/VPlusJets/interface/JetTreeFiller.h"
#include "ElectroWeakAnalysis/VPlusJets/interface/VtoElectronTreeFiller.h"
#include "ElectroWeakAnalysis/VPlusJets/interface/ElectronTreeFiller.h"
#include "ElectroWeakAnalysis/VPlusJets/interface/VtoMuonTreeFiller.h"
#include "ElectroWeakAnalysis/VPlusJets/interface/MCTreeFiller.h"

#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "PhysicsTools/Utilities/interface/Lumi3DReWeighting.h"

//
// class decleration
//
namespace ewk
{
  class VplusJetsAnalysis : public edm::EDAnalyzer {
  public:
    explicit VplusJetsAnalysis(const edm::ParameterSet&);
    ~VplusJetsAnalysis();

    virtual void beginJob();
    virtual void analyze(const edm::Event&, const edm::EventSetup& iSetup);
    virtual void declareTreeBranches();
    virtual void endJob() ;


  private:
    // ----------member data ---------------------------
    // names of modules, producing object collections
     
    /// output ROOT file for the tree and histograms

    edm::Service<TFileService> fs;
    //std::string fOutputFileName ;
    //TFile*  hOutputFile ;

    TTree*  myTree;
    bool runningOverMC_;
    bool runningOverDY_;
    std::string VBosonType_;
    std::string LeptonType_;
    edm::InputTag mInputBoson;
    edm::InputTag mPrimaryVertex;
    edm::InputTag mInputBeamSpot;
    edm::InputTag mInputMet;
//###############################
    edm::InputTag mMetRaw;
    edm::InputTag mMetCalo;
    edm::InputTag mMetCaloT1;
    edm::InputTag mMetType01;
    edm::InputTag mMetJERup;
    edm::InputTag mMetJERdown;
    edm::InputTag mMetEleup;
    edm::InputTag mMetEledown;
    edm::InputTag mMetJetup;
    edm::InputTag mMetJetdown;
    edm::InputTag mMetUncup;
    edm::InputTag mMetUncdown;
//################################
    edm::InputTag mMetMVA;
    edm::InputTag mMVAJERup;
    edm::InputTag mMVAJERdown;
    edm::InputTag mMVAEleup;
    edm::InputTag mMVAEledown;
    edm::InputTag mMVAJetup;
    edm::InputTag mMVAJetdown;
    edm::InputTag mMVAUncup;
    edm::InputTag mMVAUncdown;
//################################
    edm::InputTag mMetMVAUnity;
    edm::InputTag mMVAUnityJERup;
    edm::InputTag mMVAUnityJERdown;
    edm::InputTag mMVAUnityEleup;
    edm::InputTag mMVAUnityEledown;
    edm::InputTag mMVAUnityJetup;
    edm::InputTag mMVAUnityJetdown;
    edm::InputTag mMVAUnityUncup;
    edm::InputTag mMVAUnityUncdown;
//################################
    edm::InputTag mMetnoPU;
    edm::InputTag noPUJERup;
    edm::InputTag noPUJERdown;
    edm::InputTag noPUEleup;
    edm::InputTag noPUEledown;
    edm::InputTag noPUJetup;
    edm::InputTag noPUJetdown;
    edm::InputTag noPUUncup;
    edm::InputTag noPUUncdown;
//################################
	std::string JetsFor_rho;
	std::string puMCFile;
        std::string puDataFile;
	bool runoverAOD;
    /// The objects that actually computes variables and fill the tree 
    //std::auto_ptr<ewk::JetTreeFiller> CaloJetFiller;
    //std::auto_ptr<ewk::JetTreeFiller> CorrectedCaloJetFiller;
    std::auto_ptr<ewk::JetTreeFiller> CorrectedPatJetFiller;
    std::auto_ptr<ewk::JetTreeFiller> CorrectedSelectedJetFiller;
    //std::auto_ptr<ewk::JetTreeFiller> CorrectedPFJetFillerVBFTag; //For VBF Tag Jets
    //std::auto_ptr<ewk::JetTreeFiller> CorrectedJPTJetFiller;
    //std::auto_ptr<ewk::JetTreeFiller> GenJetFiller;
    //std::auto_ptr<ewk::JetTreeFiller> PFJetFiller; 
    //std::auto_ptr<ewk::JetTreeFiller> JPTJetFiller;
    std::auto_ptr<ewk::JetTreeFiller> SmearedJetFiller;
    std::auto_ptr<ewk::ElectronTreeFiller> ElectronFiller;
    std::auto_ptr<ewk::VtoElectronTreeFiller> recoBosonFillerE;
    std::auto_ptr<ewk::VtoMuonTreeFiller> recoBosonFillerMu;
    std::auto_ptr<ewk::MCTreeFiller> genBosonFiller;

    edm::LumiReWeighting  LumiReweightABCD_;
    edm::Lumi3DReWeighting LumiReweightABCD3d_;
    //edm::LumiReWeighting  LumiReweightABCDbug_;
    // private data members
    int run;
    int event; 
    int lumi; 
    int bunch; 
    int nPV; 
    int mNVB;
    float mPVx[60];
    float mPVy[60];
    float mPVz[60];
    float mBSx;
    float mBSy;
    float mBSz;
    float mpfMET;
    float mpfSumET;
    float mpfMETSign;
    float mpfMETPhi;
    float mcaloMET;
    float mcaloSumET;
    float mcaloMETSign;
    float mcaloMETPhi;
    float mcaloT1MET;
    float mcaloT1SumET;
    float mcaloT1METSign;
    float mcaloT1METPhi;
//##########################
    float mRawMET;
    float mRawSumET;
    float mRawMETSign;
    float mRawMETPhi;
//##########################
    float mT01MET;
    float mT01SumET;
    float mT01METSign;
    float mT01METPhi;
    float mMVAMET;
    float mMVASumET;
    float mMVAMETSign;
    float mMVAMETPhi;
    float noPUMET;
    float noPUSumET;
    float noPUMETSign;
    float noPUMETPhi;
    // MVA Unity
    float mUniMET;
    float mUniSumET;
    float mUniMETSign;
    float mUniMETPhi;
//##########################
    float mJERup;
    float mJERupSumET;
    float mJERupPhi;
    float mJERdown;
    float mJERdownSumET;
    float mJERdownPhi;
    float mvaJERup;
    float mvaJERupSumET;
    float mvaJERupPhi;
    float mvaJERdown;
    float mvaJERdownSumET;
    float mvaJERdownPhi;
    float puJERup;
    float puJERupSumET;
    float puJERupPhi;
    float puJERdown;
    float puJERdownSumET;
    float puJERdownPhi;
    // MVA Unity
    float uniJERup;
    float uniJERupSumET;
    float uniJERupPhi;
    float uniJERdown;
    float uniJERdownSumET;
    float uniJERdownPhi;
//##########################
    float mEleup;
    float mEledown;
    float mEleupSumET;
    float mEledownSumET;
    float mEleupPhi;
    float mEledownPhi;
    float mvaEleup;
    float mvaEledown;
    float mvaEleupSumET;
    float mvaEledownSumET;
    float mvaEleupPhi;
    float mvaEledownPhi;
    float puEleup;
    float puEledown;
    float puEleupSumET;
    float puEledownSumET;
    float puEleupPhi;
    float puEledownPhi;
    // MVA Unity
    float uniEleup;
    float uniEledown;
    float uniEleupSumET;
    float uniEledownSumET;
    float uniEleupPhi;
    float uniEledownPhi;
//##########################
    float mJetup;
    float mJetupSumET;
    float mJetupPhi;
    float mJetdown;
    float mJetdownSumET;
    float mJetdownPhi;
    float mvaJetup;
    float mvaJetupSumET;
    float mvaJetupPhi;
    float mvaJetdown;
    float mvaJetdownSumET;
    float mvaJetdownPhi;
    float puJetup;
    float puJetupSumET;
    float puJetupPhi;
    float puJetdown;
    float puJetdownSumET;
    float puJetdownPhi;
    // MVA Unity
    float uniJetup;
    float uniJetupSumET;
    float uniJetupPhi;
    float uniJetdown;
    float uniJetdownSumET;
    float uniJetdownPhi;
//##########################
    float mUncup;
    float mUncupSumET;
    float mUncupPhi;
    float mUncdown;
    float mUncdownSumET;
    float mUncdownPhi;
    float mvaUncup;
    float mvaUncupSumET;
    float mvaUncupPhi;
    float mvaUncdown;
    float mvaUncdownSumET;
    float mvaUncdownPhi;
    float puUncup;
    float puUncupSumET;
    float puUncupPhi;
    float puUncdown;
    float puUncdownSumET;
    float puUncdownPhi;
    // MVA Unity
    float uniUncup;
    float uniUncupSumET;
    float uniUncupPhi;
    float uniUncdown;
    float uniUncdownSumET;
    float uniUncdownPhi;
//##########################
    float fastJetRho;

//    float mcPUtrueInteractions;
    float mcPUtotnvtx;
    float mcPUbx[13];
    float mcPUnvtx[13];
    int mcPUnvtx_bxPreviuos;
    int mcPUnvtx_inTime;
    int mcPUnvtx_bxNext;

    double WeightABCD;
    double WeightABCD3d;
    //double WeightABCDbug;
  };
}
#endif
