import FWCore.ParameterSet.Config as cms

VplusJetsNtuple 	= cms.EDAnalyzer(
			'VplusJetsAnalysis',
    srcPFCor 		= cms.InputTag("cleanPatJets"),
    srcVectorBoson	= cms.InputTag("bestZee"),
    VBosonType     	= cms.string("Z"),
    LeptonType     	= cms.string("electron"),
    runningOverMC 	= cms.bool(True),  
    runningOverAOD 	= cms.bool(False),
    srcElectrons 	= cms.InputTag("cleanPatElectrons"),
    srcMet 		= cms.InputTag("patMETsPF"), #from Jet collection
##############CHANGE TAG########################
    srcMetRaw      	= cms.InputTag("patPFMet"),
    srcMetCalo     	= cms.InputTag("patCaloMet"),
    srcMetCaloT1   	= cms.InputTag("patType1CorrectedCaloMet"),
################################################
    srcMetCorrected 	= cms.InputTag("patType1CorrectedPFMet"),
    srcMetEleup    	= cms.InputTag("patType1CorrectedPFMetElectronEnUp"),
    srcMetEledown  	= cms.InputTag("patType1CorrectedPFMetElectronEnDown"),
    srcMetJetup    	= cms.InputTag("patType1CorrectedPFMetJetEnUp"),
    srcMetJetdown  	= cms.InputTag("patType1CorrectedPFMetJetEnDown"),
    srcMetUncup    	= cms.InputTag("patType1CorrectedPFMetUnclusteredEnUp"),
    srcMetUncdown  	= cms.InputTag("patType1CorrectedPFMetUnclusteredEnDown"),
################################################
    srcMVACorrected 	= cms.InputTag("patPFMetMVA"),
    srcMVAEleup    	= cms.InputTag("patPFMetMVAElectronEnUp"),
    srcMVAEledown  	= cms.InputTag("patPFMetMVAElectronEnDown"),
    srcMVAJetup    	= cms.InputTag("patPFMetMVAJetEnUp"),
    srcMVAJetdown  	= cms.InputTag("patPFMetMVAJetEnDown"),
    srcMVAUncup    	= cms.InputTag("patPFMetMVAUnclusteredEnUp"),
    srcMVAUncdown  	= cms.InputTag("patPFMetMVAUnclusteredEnDown"),
################################################
    srcnPUCorrected 	= cms.InputTag("patPFMetNoPileUp"),
    srcnPUEleup    	= cms.InputTag("patPFMetNoPileUpElectronEnUp"),
    srcnPUEledown  	= cms.InputTag("patPFMetNoPileUpElectronEnDown"),
    srcnPUJetup    	= cms.InputTag("patPFMetNoPileUpJetEnUp"),
    srcnPUJetdown  	= cms.InputTag("patPFMetNoPileUpJetEnDown"),
    srcnPUUncup    	= cms.InputTag("patPFMetNoPileUpUnclusteredEnUp"),
    srcnPUUncdown  	= cms.InputTag("patPFMetNoPileUpUnclusteredEnDown"),
################################################
    srcPrimaryVertex 	= cms.InputTag("offlinePrimaryVertices"),
    srcBeamSpot  	= cms.InputTag("offlineBeamSpot"),
    srcJetsforRho 	= cms.string("kt6PFJets"), 
    #PUMCFile           	= cms.string(""),
    #PUDataFile         	= cms.string(""),
    #PUMCHist           	= cms.string("pileup"),
    #PUDataHist         	= cms.string("pileup"),
    #HistOutFile 	= cms.string("Zee_53X_SIM_PAT.root"),
    TreeName    	= cms.string("ZJet")
)
