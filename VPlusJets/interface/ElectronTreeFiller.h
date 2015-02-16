/*****************************************************************************
 * Project: CMS detector at the CERN
 *
 * Package: ElectroWeakAnalysis/VPlusJets
 * Class:   ElectronTreeFiller
 *
 * Authors:
 *
 *   
 *
 * Copyright (C) 2010 FNAL 
 *****************************************************************************/

#ifndef ElectroWeakAnalysis_VPlusJets_ElectronTreeFiller_h
#define ElectroWeakAnalysis_VPlusJets_ElectronTreeFiller_h
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include <memory>
#include <string>
#include <iostream>
#include <vector>
#include "TTree.h" 
#include "FWCore/Framework/interface/Event.h" 
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"


namespace ewk {

  class ElectronTreeFiller {
  public:
    /// specify the name of the TTree, and the configuration for it
    ElectronTreeFiller(TTree* tree, const edm::ParameterSet iConfig );
    
    /// default constructor
    ElectronTreeFiller() {};
    
    
    /// Destructor, does nothing 
    ~ElectronTreeFiller() {};
    
    
    /// To be called once per event to fill the values for jets
    void fill(const edm::Event &iEvent);

    static const int NUM_ELE_MAX = 10; 
  protected:

    /// To be called once per event, to initialize variable values to default
    void init();
    /// Helper function for main constructor
    void SetBranches(); 
    void SetBranch( float* x, std::string name );
    void SetBranch( int* x, std::string name );
    void SetBranch( bool* x, std::string name );
    void SetBranchSingle( int* x, std::string name);

    TTree* tree_;
    std::string LeptonType_;
	edm::InputTag mInputElectrons;
        edm::InputTag mInputSelectElectrons;
	edm::InputTag mInputBeamSpot;

  private:
    // private data members

    int NumElectrons; 
    float e1E[NUM_ELE_MAX];
    float e1Et[NUM_ELE_MAX];
    float e1Pt[NUM_ELE_MAX];
    float e1Eta[NUM_ELE_MAX];
    float e1Phi[NUM_ELE_MAX];
    int e1Charge[NUM_ELE_MAX];
    float e1_sc_Pt[NUM_ELE_MAX];
    float e1_sc_Et[NUM_ELE_MAX];
    float e1_sc_Eta[NUM_ELE_MAX];
    float e1_sc_Phi[NUM_ELE_MAX];
    float e1_sc_E[NUM_ELE_MAX];
    float e1_DeltaEtaIn[NUM_ELE_MAX];
    float e1_DeltaPhiIn[NUM_ELE_MAX];
    float e1_HoverE[NUM_ELE_MAX];
    float e1_SigmaIetaIeta[NUM_ELE_MAX];
    int e1_missingHits[NUM_ELE_MAX];
    float e1_ooemoop[NUM_ELE_MAX];
    float e1_d0vtx[NUM_ELE_MAX];
    float e1_dzvtx[NUM_ELE_MAX];
    bool e1_vtxFitConversion[NUM_ELE_MAX];

    float e1_patiso_chargedHadronIso[NUM_ELE_MAX];
    float e1_patiso_photonIso[NUM_ELE_MAX];
    float e1_patiso_neutralHadronIso[NUM_ELE_MAX];
    float e1_pfIsoEAPAT[NUM_ELE_MAX];
    float e1_EffArea[NUM_ELE_MAX];
    bool ise1WP90[NUM_ELE_MAX];
    bool ise1WP80[NUM_ELE_MAX];
    bool ise1WP70[NUM_ELE_MAX];

    int NumSelectElectrons;

    float Se1E[NUM_ELE_MAX];
    float Se1Et[NUM_ELE_MAX];
    float Se1Pt[NUM_ELE_MAX];
    float Se1Eta[NUM_ELE_MAX];
    float Se1Phi[NUM_ELE_MAX];
    int Se1Charge[NUM_ELE_MAX];
    float Se1_sc_Pt[NUM_ELE_MAX];
    float Se1_sc_Et[NUM_ELE_MAX];
    float Se1_sc_Eta[NUM_ELE_MAX];
    float Se1_sc_Phi[NUM_ELE_MAX];
    float Se1_sc_E[NUM_ELE_MAX];
    float Se1_DeltaEtaIn[NUM_ELE_MAX];
    float Se1_DeltaPhiIn[NUM_ELE_MAX];
    float Se1_HoverE[NUM_ELE_MAX];
    float Se1_SigmaIetaIeta[NUM_ELE_MAX];
    int Se1_missingHits[NUM_ELE_MAX];
    float Se1_ooemoop[NUM_ELE_MAX];
    float Se1_d0vtx[NUM_ELE_MAX];
    float Se1_dzvtx[NUM_ELE_MAX];
    bool Se1_vtxFitConversion[NUM_ELE_MAX];

    float Se1_patiso_chargedHadronIso[NUM_ELE_MAX];
    float Se1_patiso_photonIso[NUM_ELE_MAX];
    float Se1_patiso_neutralHadronIso[NUM_ELE_MAX];
    float Se1_pfIsoEAPAT[NUM_ELE_MAX];
    float Se1_EffArea[NUM_ELE_MAX];
    bool isSe1WP90[NUM_ELE_MAX];
    bool isSe1WP80[NUM_ELE_MAX];
    bool isSe1WP70[NUM_ELE_MAX];


  };

} //namespace

#endif


