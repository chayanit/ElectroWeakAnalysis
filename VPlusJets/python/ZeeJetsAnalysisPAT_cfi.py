import FWCore.ParameterSet.Config as cms

VplusJetsNtuple 	= cms.EDAnalyzer(
			'VplusJetsAnalysis',
    srcPatJets		= cms.InputTag("patJetsNotOverlappingWithLeptonsForJetMEtUncertainty"),
    srcSelectedJets 	= cms.InputTag("selectedPatJets"),
    srcSmearedJets	= cms.InputTag("smearedPatJets"),
    srcVectorBoson	= cms.InputTag("bestZee"),
    VBosonType     	= cms.string("Z"),
    LeptonType     	= cms.string("electron"),
    runningOverMC 	= cms.bool(True),  
    runningOverAOD 	= cms.bool(False),
    srcPatElectrons 	= cms.InputTag("cleanPatElectrons"),
    srcSelectElectrons  = cms.InputTag("selectElectrons"),	#cleanPatElectrons pass ID/ISO
    srcMet 		= cms.InputTag("patMETsPF"), #from Jet collection
##############CHANGE TAG########################
    srcMetRaw      	= cms.InputTag("patPFMet"),
    srcMetCalo     	= cms.InputTag("patCaloMet"),
    srcMetCaloT1   	= cms.InputTag("patType1CorrectedCaloMet"),
################################################
    srcMetCorrected 	= cms.InputTag("patType1CorrectedPFMet"),
    srcMetJERup    	= cms.InputTag("patType1CorrectedPFMetJetResUp"),
    srcMetJERdown  	= cms.InputTag("patType1CorrectedPFMetJetResDown"),
    srcMetEleup    	= cms.InputTag("patType1CorrectedPFMetElectronEnUp"),
    srcMetEledown  	= cms.InputTag("patType1CorrectedPFMetElectronEnDown"),
    srcMetJetup    	= cms.InputTag("patType1CorrectedPFMetJetEnUp"),
    srcMetJetdown  	= cms.InputTag("patType1CorrectedPFMetJetEnDown"),
    srcMetUncup    	= cms.InputTag("patType1CorrectedPFMetUnclusteredEnUp"),
    srcMetUncdown  	= cms.InputTag("patType1CorrectedPFMetUnclusteredEnDown"),
################################################
    srcMVACorrected 	= cms.InputTag("patPFMetMVA"),
    srcMVAJERup    	= cms.InputTag("patPFMetMVAJetResUp"),
    srcMVAJERdown  	= cms.InputTag("patPFMetMVAJetResDown"),
    srcMVAEleup    	= cms.InputTag("patPFMetMVAElectronEnUp"),
    srcMVAEledown  	= cms.InputTag("patPFMetMVAElectronEnDown"),
    srcMVAJetup    	= cms.InputTag("patPFMetMVAJetEnUp"),
    srcMVAJetdown  	= cms.InputTag("patPFMetMVAJetEnDown"),
    srcMVAUncup    	= cms.InputTag("patPFMetMVAUnclusteredEnUp"),
    srcMVAUncdown  	= cms.InputTag("patPFMetMVAUnclusteredEnDown"),
################################################
    srcnPUCorrected 	= cms.InputTag("patPFMetNoPileUp"),
    srcnPUJERup    	= cms.InputTag("patPFMetNoPileUpJetResUp"),
    srcnPUJERdown  	= cms.InputTag("patPFMetNoPileUpJetResDown"),
    srcnPUEleup    	= cms.InputTag("patPFMetNoPileUpElectronEnUp"),
    srcnPUEledown  	= cms.InputTag("patPFMetNoPileUpElectronEnDown"),
    srcnPUJetup    	= cms.InputTag("patPFMetNoPileUpJetEnUp"),
    srcnPUJetdown  	= cms.InputTag("patPFMetNoPileUpJetEnDown"),
    srcnPUUncup    	= cms.InputTag("patPFMetNoPileUpUnclusteredEnUp"),
    srcnPUUncdown  	= cms.InputTag("patPFMetNoPileUpUnclusteredEnDown"),
################################################
    srcMVAUnityCorrected= cms.InputTag("patPFMetMVAUnity"),
    srcMVAUnityJERup    = cms.InputTag("patPFMetMVAUnityJetResUp"),
    srcMVAUnityJERdown  = cms.InputTag("patPFMetMVAUnityJetResDown"),
    srcMVAUnityEleup    = cms.InputTag("patPFMetMVAUnityElectronEnUp"),
    srcMVAUnityEledown  = cms.InputTag("patPFMetMVAUnityElectronEnDown"),
    srcMVAUnityJetup    = cms.InputTag("patPFMetMVAUnityJetEnUp"),
    srcMVAUnityJetdown  = cms.InputTag("patPFMetMVAUnityJetEnDown"),
    srcMVAUnityUncup    = cms.InputTag("patPFMetMVAUnityUnclusteredEnUp"),
    srcMVAUnityUncdown  = cms.InputTag("patPFMetMVAUnityUnclusteredEnDown"),
################################################
    srcPrimaryVertex 	= cms.InputTag("offlinePrimaryVertices"),
    srcBeamSpot  	= cms.InputTag("offlineBeamSpot"),
    srcJetsforRho 	= cms.string("kt6PFJets"), 
    srcGenParticles  	= cms.InputTag(""),
    TreeName    	= cms.string("ZJet")
)
