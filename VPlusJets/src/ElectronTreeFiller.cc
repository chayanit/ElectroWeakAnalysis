/*****************************************************************************
 * Project: CMS detector at the CERN
 *
 * Package: ElectroWeakAnalysis/VPlusJets
 *
 *
 * Authors:
 *
 *
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
#include "ElectroWeakAnalysis/VPlusJets/interface/ElectronTreeFiller.h"
#include "ElectroWeakAnalysis/VPlusJets/interface/METzCalculator.h"
#include "ElectroWeakAnalysis/VPlusJets/interface/ElectronEffectiveArea.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"


ewk::ElectronTreeFiller::ElectronTreeFiller(TTree* tree, const edm::ParameterSet iConfig)
{
    if(  iConfig.existsAs<edm::InputTag>("srcElectrons") )
		mInputElectrons  = iConfig.getParameter<edm::InputTag>("srcElectrons");
    if(  iConfig.existsAs<edm::InputTag>("srcSelectElectrons") )
                mInputSelectElectrons  = iConfig.getParameter<edm::InputTag>("srcSelectElectrons");
    if(  iConfig.existsAs<edm::InputTag>("srcBeamSpot") )
		mInputBeamSpot  = iConfig.getParameter<edm::InputTag>("srcBeamSpot");
  tree_     = tree;
  LeptonType_ = iConfig.getParameter<std::string>("LeptonType");
  if( !(tree==0) && LeptonType_=="electron") SetBranches();
}


void ewk::ElectronTreeFiller::SetBranches()
{
  // Declare jet branches
  std::string lept1 = "cleanEle";
  std::string Slept1 = "selectEle";

  SetBranchSingle( &NumElectrons, "NumcleanEle");
  SetBranch( e1E,              lept1+"_e" );
  SetBranch( e1Pt,             lept1+"_pt" );
  SetBranch( e1Et,             lept1+"_et" );
  SetBranch( e1Eta,            lept1+"_eta" );
  SetBranch( e1Phi,            lept1+"_phi" );
  SetBranch( e1Charge,         lept1+"_charge" );
  SetBranch( e1_sc_Pt,         lept1+"_sc_Pt" );
  SetBranch( e1_sc_Et,         lept1+"_sc_Et" );
  SetBranch( e1_sc_Eta,        lept1+"_sc_Eta" );
  SetBranch( e1_sc_Phi,        lept1+"_sc_Phi" );
  SetBranch( e1_sc_E,   lept1+"_sc_E" );
  SetBranch( e1_DeltaEtaIn,    lept1+"_deltaeta_in" );
  SetBranch( e1_DeltaPhiIn,    lept1+"_deltaphi_in" );
  SetBranch( e1_HoverE,        lept1+"_hovere" );
  SetBranch( e1_SigmaIetaIeta, lept1+"_sigmaietaieta" );
  SetBranch( e1_missingHits,   lept1+"_missingHits" );
  SetBranch( e1_ooemoop,       lept1+"_ooemoop" );
  SetBranch( e1_d0vtx,         lept1+"_d0vtx"  );
  SetBranch( e1_dzvtx,         lept1+"_dzvtx"  );
  SetBranch( e1_vtxFitConversion,                 lept1+"_vtxFitConversion" );
  SetBranch( e1_patiso_chargedHadronIso,         lept1+"_patiso_chargedHadronIso" );
  SetBranch( e1_patiso_photonIso,                lept1+"_patiso_photonIso" );
  SetBranch( e1_patiso_neutralHadronIso,         lept1+"_patiso_neutralHadronIso" );
  SetBranch( e1_pfIsoEAPAT,         lept1+"_pfIsoEAPAT" );
  SetBranch( ise1WP90,         lept1+"_isWP90" );
  SetBranch( ise1WP80,         lept1+"_isWP80" );
  SetBranch( ise1WP70,         lept1+"_isWP70" );
  SetBranchSingle( &NumSelectElectrons, "NumselectEle");
  SetBranch( Se1E,              Slept1+"_e" );
  SetBranch( Se1Pt,             Slept1+"_pt" );
  SetBranch( Se1Et,             Slept1+"_et" );
  SetBranch( Se1Eta,            Slept1+"_eta" );
  SetBranch( Se1Phi,            Slept1+"_phi" );
  SetBranch( Se1Charge,         Slept1+"_charge" );
  SetBranch( Se1_sc_Pt,         Slept1+"_sc_Pt" );
  SetBranch( Se1_sc_Et,         Slept1+"_sc_Et" );	 
  SetBranch( Se1_sc_Eta,        Slept1+"_sc_Eta" );
  SetBranch( Se1_sc_Phi,        Slept1+"_sc_Phi" );
  SetBranch( Se1_sc_E, 	 Slept1+"_sc_E" );
  SetBranch( Se1_DeltaEtaIn,    Slept1+"_deltaeta_in" );
  SetBranch( Se1_DeltaPhiIn,    Slept1+"_deltaphi_in" );
  SetBranch( Se1_HoverE,        Slept1+"_hovere" );	    
  SetBranch( Se1_SigmaIetaIeta, Slept1+"_sigmaietaieta" );
  SetBranch( Se1_missingHits,   Slept1+"_missingHits" );	 
  SetBranch( Se1_ooemoop,       Slept1+"_ooemoop" );
  SetBranch( Se1_d0vtx,		Slept1+"_d0vtx"  );
  SetBranch( Se1_dzvtx,         Slept1+"_dzvtx"  );
  SetBranch( Se1_vtxFitConversion, 		   Slept1+"_vtxFitConversion" );
  SetBranch( Se1_patiso_chargedHadronIso,         Slept1+"_patiso_chargedHadronIso" );
  SetBranch( Se1_patiso_photonIso,                Slept1+"_patiso_photonIso" );
  SetBranch( Se1_patiso_neutralHadronIso,         Slept1+"_patiso_neutralHadronIso" );
  SetBranch( Se1_pfIsoEAPAT,         Slept1+"_pfIsoEAPAT" );
  SetBranch( isSe1WP90,         Slept1+"_isWP90" );
  SetBranch( isSe1WP80,         Slept1+"_isWP80" );
  SetBranch( isSe1WP70,         Slept1+"_isWP70" );
}
/////////////////////////////////////////////////////////////////////////


void ewk::ElectronTreeFiller::init()   
{
  // initialize private data members
  NumElectrons = 0;
  NumSelectElectrons = 0;

  for (int j =0; j< NUM_ELE_MAX; ++j) {
  e1E[j]                    = -1.;
  e1Et[j]               = -1.;
  e1Pt[j]               = -1.;
  e1Eta[j]              = -10.;
  e1Phi[j]              = -10.;
  e1Charge[j]           = -10;
  e1_sc_Pt[j]           = -1.;
  e1_sc_Et[j]           = -1.;
  e1_sc_Eta[j]          = -10.;
  e1_sc_Phi[j]          = -10.;
  e1_sc_E[j]            = -1.;
  e1_DeltaEtaIn[j]      = -10.;
  e1_DeltaPhiIn[j]      = -10.;
  e1_HoverE[j]          = -1.;
  e1_SigmaIetaIeta[j]   = -1.;
  e1_missingHits[j]     = 100;
  e1_ooemoop[j]             = -10.;
  e1_d0vtx[j]       = -10.;
  e1_dzvtx[j]           = -10.;
  e1_vtxFitConversion[j] = true;
  e1_patiso_chargedHadronIso[j]      = -99999.;
  e1_patiso_photonIso[j]             = -99999.;
  e1_patiso_neutralHadronIso[j]      = -99999.;
  e1_EffArea[j]                      = 0.;
  e1_pfIsoEAPAT[j]                   = -99999.;
  ise1WP90[j]          = false;
  ise1WP80[j]          = false;
  ise1WP70[j]          = false;
  }

  for (int i =0; i< NUM_ELE_MAX; ++i) {
  Se1E[i]		     = -1.;
  Se1Et[i]               = -1.;
  Se1Pt[i]               = -1.;
  Se1Eta[i]              = -10.;
  Se1Phi[i]              = -10.;
  Se1Charge[i]           = -10;
  Se1_sc_Pt[i]           = -1.;
  Se1_sc_Et[i]           = -1.;	 
  Se1_sc_Eta[i]          = -10.;    
  Se1_sc_Phi[i]          = -10.;
  Se1_sc_E[i]            = -1.;
  Se1_DeltaEtaIn[i]      = -10.;
  Se1_DeltaPhiIn[i]      = -10.;
  Se1_HoverE[i]          = -1.;	
  Se1_SigmaIetaIeta[i]   = -1.;	
  Se1_missingHits[i]     = 100;
  Se1_ooemoop[i]	     = -10.;
  Se1_d0vtx[i]	     = -10.;
  Se1_dzvtx[i]           = -10.;
  Se1_vtxFitConversion[i] = true;
  Se1_patiso_chargedHadronIso[i]      = -99999.;
  Se1_patiso_photonIso[i]             = -99999.;
  Se1_patiso_neutralHadronIso[i]      = -99999.;
  Se1_EffArea[i]                      = 0.;
  Se1_pfIsoEAPAT[i]                   = -99999.;
  isSe1WP90[i]          = false;
  isSe1WP80[i]          = false;
  isSe1WP70[i]          = false;
  
  }
  // initialization done
}

void ewk::ElectronTreeFiller::fill(const edm::Event& iEvent)
{
  // protection
  if( (tree_==0) || !(LeptonType_=="electron") )  return;

  // first initialize to the default values
  init();

  /////// Pileup density "rho" in the event from fastJet pileup calculation /////
  edm::Handle<double> rho;
  const edm::InputTag eventrho("kt6PFJets", "rho");
  iEvent.getByLabel(eventrho,rho);
  double fastJetRho = *rho;

  edm::Handle<reco::VertexCollection> vtxs;
  iEvent.getByLabel("offlinePrimaryVertices", vtxs);

  edm::Handle<reco::ConversionCollection> conversions;
  iEvent.getByLabel("allConversions", conversions);

  edm::Handle<reco::BeamSpot> beamspot_h;
  iEvent.getByLabel(mInputBeamSpot, beamspot_h);
  const reco::BeamSpot &beamspot = *(beamspot_h.product());

  edm::Handle<edm::View<pat::Electron> > Pelectrons;
  iEvent.getByLabel(mInputElectrons, Pelectrons);

  edm::Handle<edm::View<pat::Electron> > Selectrons;
  iEvent.getByLabel(mInputSelectElectrons, Selectrons);

  for(unsigned int iLep = 0; iLep < (*Pelectrons).size() && iLep < (unsigned int) NUM_ELE_MAX; ++iLep)
  {
        const pat::Electron& ele = Pelectrons->at(iLep);
        NumElectrons++;
	
        e1E[iLep] = ele.energy();
        e1Et[iLep] = ele.et();
        e1Pt[iLep] = ele.pt();
        e1Eta[iLep] = ele.eta();
        e1Phi[iLep] = ele.phi();
        e1Charge[iLep] = ele.charge();
        e1_sc_Pt[iLep] = e1E[iLep] / cosh(e1Eta[iLep]);
        e1_sc_Et[iLep] = e1Pt[iLep]; 
        e1_sc_Eta[iLep] = ele.superCluster()->eta();
        e1_sc_Phi[iLep] = ele.superCluster()->phi();
        e1_sc_E[iLep] = ele.superCluster()->energy();
        e1_DeltaEtaIn[iLep] = ele.deltaEtaSuperClusterTrackAtVtx();
        e1_DeltaPhiIn[iLep] = ele.deltaPhiSuperClusterTrackAtVtx();
        e1_HoverE[iLep] = ele.hadronicOverEm();      
        e1_SigmaIetaIeta[iLep] = ele.sigmaIetaIeta();
        e1_missingHits[iLep] = ele.gsfTrack()->trackerExpectedHitsInner().numberOfHits();
        e1_ooemoop[iLep] = fabs((1.0/ele.ecalEnergy() - ele.eSuperClusterOverP()/ele.ecalEnergy()));
        
        if (vtxs->size() > 0) {
                reco::VertexRef vtx(vtxs, 0);
                e1_d0vtx[iLep] = ele.gsfTrack()->dxy(vtx->position());
                e1_dzvtx[iLep] = ele.gsfTrack()->dz(vtx->position());
        } else {
                e1_d0vtx[iLep] = ele.gsfTrack()->dxy();
                e1_dzvtx[iLep] = ele.gsfTrack()->dz();
         }
        e1_vtxFitConversion[iLep] = ConversionTools::hasMatchedConversion(ele, conversions, beamspot.position());

        e1_patiso_chargedHadronIso[iLep] = ele.chargedHadronIso();
        e1_patiso_photonIso[iLep] = ele.photonIso();
        e1_patiso_neutralHadronIso[iLep] = ele.neutralHadronIso();
        e1_EffArea[iLep] = ElectronEffectiveArea::GetElectronEffectiveArea( ElectronEffectiveArea::kEleGammaAndNeutralHadronIso03 , e1_sc_Eta[iLep] , ElectronEffectiveArea::kEleEAData2012);
        e1_pfIsoEAPAT[iLep] = (e1_patiso_chargedHadronIso[iLep] +
                        std::max(0.,e1_patiso_neutralHadronIso[iLep] + e1_patiso_photonIso[iLep]  - e1_EffArea[iLep]*fastJetRho)) / e1Pt[iLep];

    ise1WP70[iLep] = (e1Pt[iLep] > 20.) && (e1_missingHits[iLep] <= 0) &&
               (e1_pfIsoEAPAT[iLep] < 0.10) && (!e1_vtxFitConversion[iLep]) &&
               ((ele.isEB() && e1_SigmaIetaIeta[iLep]<0.01 && fabs(e1_DeltaPhiIn[iLep])<0.03 && fabs(e1_DeltaEtaIn[iLep])<0.004 && e1_HoverE[iLep]<0.120 && e1_ooemoop[iLep]<0.050 && fabs(e1_d0vtx[iLep])<0.020 && fabs(e1_dzvtx[iLep])<0.100) ||
                (ele.isEE() && e1_SigmaIetaIeta[iLep]<0.03 && fabs(e1_DeltaPhiIn[iLep])<0.02 && fabs(e1_DeltaEtaIn[iLep])<0.005 && e1_HoverE[iLep]<0.100 && e1_ooemoop[iLep]<0.050 && fabs(e1_d0vtx[iLep])<0.020 && fabs(e1_dzvtx[iLep])<0.100));

    ise1WP80[iLep] = (e1Pt[iLep] > 20.) && (e1_missingHits[iLep] <= 1) &&
               (e1_pfIsoEAPAT[iLep] < 0.15) && (!e1_vtxFitConversion[iLep]) &&
               ((ele.isEB() && e1_SigmaIetaIeta[iLep]<0.01 && fabs(e1_DeltaPhiIn[iLep])<0.06 && fabs(e1_DeltaEtaIn[iLep])<0.004 && e1_HoverE[iLep]<0.120 && e1_ooemoop[iLep]<0.050 && fabs(e1_d0vtx[iLep])<0.020 && fabs(e1_dzvtx[iLep])<0.100) ||
                (ele.isEE() && e1_SigmaIetaIeta[iLep]<0.03 && fabs(e1_DeltaPhiIn[iLep])<0.03 && fabs(e1_DeltaEtaIn[iLep])<0.007 && e1_HoverE[iLep]<0.100 && e1_ooemoop[iLep]<0.050 && fabs(e1_d0vtx[iLep])<0.020 && fabs(e1_dzvtx[iLep])<0.100));

    ise1WP90[iLep] = (e1Pt[iLep] > 20.) && (e1_missingHits[iLep] <= 1) &&
               (e1_pfIsoEAPAT[iLep] < 0.15) && (!e1_vtxFitConversion[iLep]) &&
               ((ele.isEB() && e1_SigmaIetaIeta[iLep]<0.01 && fabs(e1_DeltaPhiIn[iLep])<0.15 && fabs(e1_DeltaEtaIn[iLep])<0.007 && e1_HoverE[iLep]<0.120 && e1_ooemoop[iLep]<0.050 && fabs(e1_d0vtx[iLep])<0.020 && fabs(e1_dzvtx[iLep])<0.200) ||
                (ele.isEE() && e1_SigmaIetaIeta[iLep]<0.03 && fabs(e1_DeltaPhiIn[iLep])<0.10 && fabs(e1_DeltaEtaIn[iLep])<0.009 && e1_HoverE[iLep]<0.100 && e1_ooemoop[iLep]<0.050 && fabs(e1_d0vtx[iLep])<0.020 && fabs(e1_dzvtx[iLep])<0.200));

  }

  for(unsigned int iLep = 0; iLep < (*Selectrons).size() && iLep < (unsigned int) NUM_ELE_MAX; ++iLep)
  {
        const pat::Electron& ele = Selectrons->at(iLep);
	NumSelectElectrons++;

        Se1E[iLep] = ele.energy();
        Se1Et[iLep] = ele.et();
        Se1Pt[iLep] = ele.pt();
        Se1Eta[iLep] = ele.eta();
        Se1Phi[iLep] = ele.phi();
        Se1Charge[iLep] = ele.charge();
        Se1_sc_Pt[iLep] = Se1E[iLep] / cosh(Se1Eta[iLep]);
        Se1_sc_Et[iLep] = Se1Pt[iLep];
        Se1_sc_Eta[iLep] = ele.superCluster()->eta();
        Se1_sc_Phi[iLep] = ele.superCluster()->phi();
        Se1_sc_E[iLep] = ele.superCluster()->energy();
        Se1_DeltaEtaIn[iLep] = ele.deltaEtaSuperClusterTrackAtVtx();
        Se1_DeltaPhiIn[iLep] = ele.deltaPhiSuperClusterTrackAtVtx();
        Se1_HoverE[iLep] = ele.hadronicOverEm();
        Se1_SigmaIetaIeta[iLep] = ele.sigmaIetaIeta();
        Se1_missingHits[iLep] = ele.gsfTrack()->trackerExpectedHitsInner().numberOfHits();
        Se1_ooemoop[iLep] = fabs((1.0/ele.ecalEnergy() - ele.eSuperClusterOverP()/ele.ecalEnergy()));

        if (vtxs->size() > 0) {
                reco::VertexRef vtx(vtxs, 0);
                Se1_d0vtx[iLep] = ele.gsfTrack()->dxy(vtx->position());
                Se1_dzvtx[iLep] = ele.gsfTrack()->dz(vtx->position());
        } else {
                Se1_d0vtx[iLep] = ele.gsfTrack()->dxy();
                Se1_dzvtx[iLep] = ele.gsfTrack()->dz();
         }
        Se1_vtxFitConversion[iLep] = ConversionTools::hasMatchedConversion(ele, conversions, beamspot.position());

        Se1_patiso_chargedHadronIso[iLep] = ele.chargedHadronIso();
        Se1_patiso_photonIso[iLep] = ele.photonIso();
        Se1_patiso_neutralHadronIso[iLep] = ele.neutralHadronIso();
        Se1_EffArea[iLep] = ElectronEffectiveArea::GetElectronEffectiveArea( ElectronEffectiveArea::kEleGammaAndNeutralHadronIso03 , Se1_sc_Eta[iLep] , ElectronEffectiveArea::kEleEAData2012);
        Se1_pfIsoEAPAT[iLep] = (Se1_patiso_chargedHadronIso[iLep] +
                        std::max(0.,Se1_patiso_neutralHadronIso[iLep] + Se1_patiso_photonIso[iLep]  - Se1_EffArea[iLep]*fastJetRho)) / Se1Pt[iLep];

    isSe1WP70[iLep] = (Se1Pt[iLep] > 20.) && (Se1_missingHits[iLep] <= 0) &&
               (Se1_pfIsoEAPAT[iLep] < 0.10) && (!Se1_vtxFitConversion[iLep]) &&
               ((ele.isEB() && Se1_SigmaIetaIeta[iLep]<0.01 && fabs(Se1_DeltaPhiIn[iLep])<0.03 && fabs(Se1_DeltaEtaIn[iLep])<0.004 && Se1_HoverE[iLep]<0.120 && Se1_ooemoop[iLep]<0.050 && fabs(Se1_d0vtx[iLep])<0.020 && fabs(Se1_dzvtx[iLep])<0.100) ||
                (ele.isEE() && Se1_SigmaIetaIeta[iLep]<0.03 && fabs(Se1_DeltaPhiIn[iLep])<0.02 && fabs(Se1_DeltaEtaIn[iLep])<0.005 && Se1_HoverE[iLep]<0.100 && Se1_ooemoop[iLep]<0.050 && fabs(Se1_d0vtx[iLep])<0.020 && fabs(Se1_dzvtx[iLep])<0.100));

    isSe1WP80[iLep] = (Se1Pt[iLep] > 20.) && (Se1_missingHits[iLep] <= 1) &&
               (Se1_pfIsoEAPAT[iLep] < 0.15) && (!Se1_vtxFitConversion[iLep]) &&
               ((ele.isEB() && Se1_SigmaIetaIeta[iLep]<0.01 && fabs(Se1_DeltaPhiIn[iLep])<0.06 && fabs(Se1_DeltaEtaIn[iLep])<0.004 && Se1_HoverE[iLep]<0.120 && Se1_ooemoop[iLep]<0.050 && fabs(Se1_d0vtx[iLep])<0.020 && fabs(Se1_dzvtx[iLep])<0.100) ||
                (ele.isEE() && Se1_SigmaIetaIeta[iLep]<0.03 && fabs(Se1_DeltaPhiIn[iLep])<0.03 && fabs(Se1_DeltaEtaIn[iLep])<0.007 && Se1_HoverE[iLep]<0.100 && Se1_ooemoop[iLep]<0.050 && fabs(Se1_d0vtx[iLep])<0.020 && fabs(Se1_dzvtx[iLep])<0.100));

    isSe1WP90[iLep] = (Se1Pt[iLep] > 20.) && (Se1_missingHits[iLep] <= 1) &&
               (Se1_pfIsoEAPAT[iLep] < 0.15) && (!Se1_vtxFitConversion[iLep]) &&
               ((ele.isEB() && Se1_SigmaIetaIeta[iLep]<0.01 && fabs(Se1_DeltaPhiIn[iLep])<0.15 && fabs(Se1_DeltaEtaIn[iLep])<0.007 && Se1_HoverE[iLep]<0.120 && Se1_ooemoop[iLep]<0.050 && fabs(Se1_d0vtx[iLep])<0.020 && fabs(Se1_dzvtx[iLep])<0.200) ||
                (ele.isEE() && Se1_SigmaIetaIeta[iLep]<0.03 && fabs(Se1_DeltaPhiIn[iLep])<0.10 && fabs(Se1_DeltaEtaIn[iLep])<0.009 && Se1_HoverE[iLep]<0.100 && Se1_ooemoop[iLep]<0.050 && fabs(Se1_d0vtx[iLep])<0.020 && fabs(Se1_dzvtx[iLep])<0.200));
 
    if(isSe1WP90[iLep]) std::cout<<"This is loose SelectElectron"<<std::endl;

  }

}


////////////////// utilities, helpers ///////////////////
 
void ewk::ElectronTreeFiller::SetBranch( float* x, std::string name)
{
  std::string brName = name;
  tree_->Branch( brName.c_str(), x, ( brName+"[10]/F").c_str() );
}


void ewk::ElectronTreeFiller::SetBranch( int* x, std::string name)
{
  std::string brName = name;
  tree_->Branch( brName.c_str(), x, ( brName+"[10]/I").c_str() );
}


void ewk::ElectronTreeFiller::SetBranch( bool* x, std::string name)
{
  std::string brName = name;
  tree_->Branch( brName.c_str(), x, ( brName+"[10]/O").c_str() );
}


void ewk::ElectronTreeFiller::SetBranchSingle( int* x, std::string name)
{
  std::string brName = name;
  tree_->Branch( brName.c_str(), x, ( brName+"/I").c_str() );

}

