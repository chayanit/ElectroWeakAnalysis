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

#ifndef ElectroWeakAnalysis_VPlusJets_JetTreeFiller_h
#define ElectroWeakAnalysis_VPlusJets_JetTreeFiller_h

#include <memory>
#include <string>
#include <iostream>
#include <vector>
#include "TTree.h" 

#include "FWCore/Framework/interface/Event.h" 
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "ElectroWeakAnalysis/VPlusJets/interface/QGLikelihoodCalculator.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/BTauReco/interface/SecondaryVertexTagInfo.h"
#include "DataFormats/BTauReco/interface/TrackIPTagInfo.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"



namespace ewk {
  class JetTreeFiller {
  public:
    /// specify the name of the TTree, and the configuration for it
    JetTreeFiller(const char *name, TTree* tree, 
		  const std::string jetType,
		  const edm::ParameterSet iConfig );


    /// default constructor
    JetTreeFiller() {};

    /// To be called once per event to fill the values for jets
    void fill(const edm::Event &iEvent);

    static const int NUM_JET_MAX = 20;

  protected:

    /// To be called once per event, to initialize variable values to default
    void init() const ;
    /// Helper function for main constructor
    void SetBranches(); 
    void SetBranch( float* x, std::string name);
    void SetBranch( int* x, std::string name);
    void SetBranchSingle( float* x, std::string name);
    void SetBranchSingle( int* x, std::string name);

    void FillBranches() const;
    void init();


    template <typename T1> 
      void fillBasicJetQuantities(int iJet, const T1& pfjet); 
    void fillBasicJetQuantities(int iJet); 

    void fillBtagInfoPAT(int iJet, const pat::Jet* pjet);
    TTree* tree_;
    std::string jetType_;
    std::string Vtype_;
    std::string LeptonType_;
    edm::InputTag mInputJets;
    edm::InputTag mInputBoson;
  	bool runoverAOD;
    mutable std::vector<std::string> bnames;

  private:
    // private data members
    
    int NumJets; 
    float Et[NUM_JET_MAX];
    float Pt[NUM_JET_MAX];
    float Eta[NUM_JET_MAX];
    float Phi[NUM_JET_MAX];
    float Theta[NUM_JET_MAX];
    float E[NUM_JET_MAX];
    float Y[NUM_JET_MAX];
    float Mass[NUM_JET_MAX];
    float Px[NUM_JET_MAX];
    float Py[NUM_JET_MAX];
    float Pz[NUM_JET_MAX];
    float bDiscriminator[NUM_JET_MAX];
    float bDiscriminatorSSVHE[NUM_JET_MAX];
    float bDiscriminatorTCHE[NUM_JET_MAX];
    float bDiscriminatorCSV[NUM_JET_MAX];
    float bDiscriminatorJP[NUM_JET_MAX];
    float bDiscriminatorSSVHP[NUM_JET_MAX];
    float bDiscriminatorTCHP[NUM_JET_MAX];
    int numBTags;
  };

} //namespace

#endif


