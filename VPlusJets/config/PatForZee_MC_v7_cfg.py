import FWCore.ParameterSet.Config as cms

process = cms.Process("ZeePATTuple")

process.load('Configuration/StandardSequences/Services_cff')
process.load('Configuration/StandardSequences/MagneticField_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')

###--------------------------------------------------------------
### load PAT, because it insists on defining the process for you
from PhysicsTools.PatAlgos.patTemplate_cfg import *


###--------------------------------------------------------------
### Set these parameters
iRunOnData = False
isMC = True
EleIdLabel = "tight" #tight, loose, medium, veto

sfNoPUjetOffsetEnCorr = None
if iRunOnData:
	sfNoPUjetOffsetEnCorr = 0.2
else:
        sfNoPUjetOffsetEnCorr = 0.0

# Input + GlobalTag
if iRunOnData == True:
    process.source.fileNames =  ['/store/data/Run2012D/DoubleElectron/AOD/22Jan2013-v1/10000/780078B6-0D8F-E211-A4D3-00261894395F.root']
    process.GlobalTag.globaltag = cms.string('FT_53_V21_AN4::All')
else:
    process.source.fileNames =  ['/store/mc/Summer12_DR53X/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S10_START53_V7A-v1/0000/7EFFED8E-09D2-E111-8F10-001E6739723D.root']
    process.GlobalTag.globaltag = cms.string('START53_V23::All')  #MC

process.maxEvents.input = 50

###--------------------------------------------------------------
### Logger
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport = cms.untracked.PSet(
    reportEvery = cms.untracked.int32(2000),
    limit = cms.untracked.int32(10000000)
)


###--------------------------------------------------------------
### basic filters
# trigger filter
if iRunOnData == True:
	process.load('HLTrigger.HLTfilters.hltHighLevel_cfi')
	process.hltHighLevel.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")
	process.hltHighLevel.throw = cms.bool(False)
	process.hltHighLevel.HLTPaths = cms.vstring(
		'HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*'
	)
else:
        process.load('HLTrigger.HLTfilters.hltHighLevel_cfi')
        process.hltHighLevel.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")
        process.hltHighLevel.throw = cms.bool(False)
        process.hltHighLevel.HLTPaths = cms.vstring(
                '*'
	)


# track quality filter
process.noscraping = cms.EDFilter(
    "FilterOutScraping",
    applyfilter = cms.untracked.bool(True),
    debugOn = cms.untracked.bool(False),
    numtrack = cms.untracked.uint32(10),
    thresh = cms.untracked.double(0.25)
    )

# require a good vertex
process.primaryVertexFilter = cms.EDFilter(
    "VertexSelector",
    src = cms.InputTag("offlinePrimaryVertices"),
    cut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.Rho <= 2"),
    filter = cms.bool(True)
    )


###--------------------------------------------------------------
### MET Filters
# The iso-based HBHE noise filter
process.load('CommonTools.RecoAlgos.HBHENoiseFilter_cfi')

# Beam halo filter
process.load('RecoMET.METAnalyzers.CSCHaloFilter_cfi')

# HCAL laser filter (no monidified parameter from V00-00-09 of RecoMET.METFilters)
process.load('RecoMET.METFilters.hcalLaserEventFilter_cfi')

# HCAL rawtodigi laser filter (latest update)
#process.load("EventFilter.HcalRawToDigi.hcallasereventfilter2012_cff")
#process.hcallasereventfilter2012.eventFileName = cms.string('ElectroWeakAnalysis/VPlusJets/config/HCALLaser2012AllDatasets.txt.gz')

# HCAL laser events for 2013 ReReco
process.load("EventFilter.HcalRawToDigi.hcallaserFilterFromTriggerResult_cff")

# ECAL dead cells filter
process.load('RecoMET.METFilters.EcalDeadCellTriggerPrimitiveFilter_cfi')

# ECAL bad supercluster filter
process.load('RecoMET.METFilters.eeBadScFilter_cfi')

# The ECAL laser correction filter
process.load('RecoMET.METFilters.ecalLaserCorrFilter_cfi')

# The Good vertices collection needed by the tracking failure filter
process.goodVertices = cms.EDFilter(
  "VertexSelector",
  filter = cms.bool(False),
  src = cms.InputTag("offlinePrimaryVertices"),
  cut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.rho < 2")
)

# tracking failure filter
process.load('RecoMET.METFilters.trackingFailureFilter_cfi')

# tracking POG filters
process.load('RecoMET.METFilters.trackingPOGFilters_cff')


###--------------------------------------------------------------
### customise PAT
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')

from PhysicsTools.PatAlgos.tools.coreTools import *

if iRunOnData == True:
    removeMCMatching(process, ['All'])
    runOnData(process)

# Add the PFMET
from PhysicsTools.PatAlgos.tools.metTools import *
addPfMET(process, 'PF')


# Switch to PF jets
from PhysicsTools.PatAlgos.tools.jetTools import *
if iRunOnData == True:
	switchJetCollection(process,
                            cms.InputTag('ak5PFJets'),
                            doJTA            = True,  
                            doBTagging       = True,
                            jetCorrLabel     = ('AK5PF', ['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual','Uncertainty']),
                            doType1MET       = False,            
                            genJetCollection = cms.InputTag("ak5GenJets"),
                            doJetID          = True,
                            jetIdLabel       = "ak5"
                            )
else:                   
        switchJetCollection(process,
                            cms.InputTag('ak5PFJets'),
                            doJTA            = True,
                            doBTagging       = True,
                            jetCorrLabel     = ('AK5PF', ['L1FastJet', 'L2Relative', 'L3Absolute','Uncertainty']),
                            doType1MET       = False,
                            genJetCollection = cms.InputTag("ak5GenJets"),
                            doJetID          = True,
                            jetIdLabel       = "ak5"
                            )


###--------------------------------------------------------------
# Electron PF-isolation 0.3
from PhysicsTools.PatAlgos.tools.pfTools import *
usePFIso( process )

# if the sample does not contain value map of PF candidate "particleFlow:electrons", use following line.
#process.patMuons.pfMuonSource = 'particleFlow'
#process.patElectrons.pfElectronSource = 'particleFlow'
 
process.patMuons.isoDeposits = cms.PSet()
process.patMuons.isolationValues = cms.PSet()

process.patElectrons.isolationValues = cms.PSet(
        pfChargedHadrons = cms.InputTag("elPFIsoValueCharged03PFIdPFIso"),
        pfChargedAll = cms.InputTag("elPFIsoValueChargedAll03PFIdPFIso"),
        pfPUChargedHadrons = cms.InputTag("elPFIsoValuePU03PFIdPFIso"),
        pfNeutralHadrons = cms.InputTag("elPFIsoValueNeutral03PFIdPFIso"),
        pfPhotons = cms.InputTag("elPFIsoValueGamma03PFIdPFIso")
        )
process.patElectrons.isolationValuesNoPFId = cms.PSet(
        pfChargedHadrons = cms.InputTag("elPFIsoValueCharged03NoPFIdPFIso"),
        pfChargedAll = cms.InputTag("elPFIsoValueChargedAll03NoPFIdPFIso"),
        pfPUChargedHadrons = cms.InputTag("elPFIsoValuePU03NoPFIdPFIso"),
        pfNeutralHadrons = cms.InputTag("elPFIsoValueNeutral03NoPFIdPFIso"),
        pfPhotons = cms.InputTag("elPFIsoValueGamma03NoPFIdPFIso")
        )
process.patDefaultSequence.replace(process.patElectrons,process.eleIsoSequence+process.patElectrons)


# Jet and Electrons selection
process.selectedPatJets.cut = cms.string("et>10.")
process.selectedPatElectrons.cut = cms.string("et>10. && abs(eta)<2.6")


###--------------------------------------------------------------
# Selected Electrons (ID/ISO)
process.load('ElectroWeakAnalysis.VPlusJets.SelectElectrons_cfi')
process.selectElectrons.idLabel = cms.string(EleIdLabel)


###--------------------------------------------------------------
# apply type I/type I + II PFMEt corrections to pat::MET object
# and estimate systematic uncertainties on MET

if iRunOnData == True:
	from PhysicsTools.PatUtils.tools.runType1PFMEtUncertainties     import runType1PFMEtUncertainties
	runType1PFMEtUncertainties(process,
		electronCollection = cms.InputTag('selectElectrons'),
		photonCollection = '',
		muonCollection = '',
		tauCollection = '',
		jetCollection = cms.InputTag('patJets'),
		jetCorrLabel = 'L2L3Residual',
                jecUncertaintyFile = "PhysicsTools/PatUtils/data/Summer13_V1_DATA_UncertaintySources_AK5PF.txt",
		makeType1corrPFMEt = True,
		makeType1p2corrPFMEt = True,
		doApplyType0corr = True,
		doSmearJets = False,
		addToPatDefaultSequence = False,
		postfix = ''
	)
        process.patPFJetMETtype1p2Corr.jetCorrLabel = cms.string('L2L3Residual')
	process.patPFJetMETtype2Corr.jetCorrLabel = cms.string('L2L3Residual')

	from PhysicsTools.PatUtils.tools.runNoPileUpMEtUncertainties    import runNoPileUpMEtUncertainties
	runNoPileUpMEtUncertainties(process,
		electronCollection = cms.InputTag('selectElectrons'),
		photonCollection = '',
		muonCollection = '',
		tauCollection = '',
		jetCollection = cms.InputTag('patJets'),
                jecUncertaintyFile = "PhysicsTools/PatUtils/data/Summer13_V1_DATA_UncertaintySources_AK5PF.txt",
		doSmearJets = False,
		sfNoPUjetOffsetEnCorr = sfNoPUjetOffsetEnCorr,
		addToPatDefaultSequence = False,
                postfix = ''
        )
        process.calibratedAK5PFJetsForNoPileUpPFMEt.correctors = cms.vstring('ak5PFL1FastL2L3Residual')
        process.noPileUpPFMEt.srcLeptons = cms.VInputTag('selectElectrons')

	from PhysicsTools.PatUtils.tools.runMVAMEtUncertainties         import runMVAMEtUncertainties
	runMVAMEtUncertainties(process,
		electronCollection = cms.InputTag('selectElectrons'),
                photonCollection = '',
                muonCollection = '',
                tauCollection = '',
                jetCollection = cms.InputTag('patJets'),
                jecUncertaintyFile = "PhysicsTools/PatUtils/data/Summer13_V1_DATA_UncertaintySources_AK5PF.txt",
                doSmearJets = False,
                addToPatDefaultSequence = False,
                postfix = ''
        )
	process.calibratedAK5PFJetsForPFMEtMVA.correctors = cms.vstring("ak5PFL1FastL2L3Residual")
	process.pfMEtMVA.srcLeptons = cms.VInputTag('selectElectrons')

#       from PhysicsTools.PatUtils.tools.runType1CaloMEtUncertainties   import runType1CaloMEtUncertainties
#        runType1CaloMEtUncertainties(process,
#                electronCollection = cms.InputTag('selectElectrons'),
#                photonCollection = '',
#                muonCollection = '',
#                tauCollection = '',
#                jetCollection = cms.InputTag('patJets'),
#               caloTowerCollection = cms.InputTag('towerMaker'),
#                jetCorrPayloadName = 'AK5Calo',
#                jetCorrLabelUpToL3 = 'ak5CaloL1L2L3',
#                jetCorrLabelUpToL3Res = 'ak5CaloL1L2L3Residual',
#               jecUncertaintyFile = "PhysicsTools/PatUtils/data/Fall12_V7_DATA_UncertaintySources_AK5Calo.txt",
#               jecUncertaintyTag = 'SubTotalMC',
#                addToPatDefaultSequence = False,
#                postfix = ''
#        )		
else:
	from PhysicsTools.PatUtils.tools.runType1PFMEtUncertainties     import runType1PFMEtUncertainties
        runType1PFMEtUncertainties(process,
                electronCollection = cms.InputTag('selectElectrons'),
                photonCollection = '',
                muonCollection = '',
                tauCollection = '',
                jetCollection = cms.InputTag('patJets'),
                jetCorrLabel = 'L3Absolute',
                jecUncertaintyFile = "PhysicsTools/PatUtils/data/Summer13_V1_DATA_UncertaintySources_AK5PF.txt",
                makeType1corrPFMEt = True,
                makeType1p2corrPFMEt = True,
                doApplyType0corr = True,
                doSmearJets = True,
                addToPatDefaultSequence = False,
                postfix = ''
        )
	from PhysicsTools.PatUtils.tools.runNoPileUpMEtUncertainties    import runNoPileUpMEtUncertainties
        runNoPileUpMEtUncertainties(process,
                electronCollection = cms.InputTag('selectElectrons'),
                photonCollection = '',
                muonCollection = '',
                tauCollection = '',
                jetCollection = cms.InputTag('patJets'),
                jecUncertaintyFile = "PhysicsTools/PatUtils/data/Summer13_V1_DATA_UncertaintySources_AK5PF.txt",
                doSmearJets = True,
                sfNoPUjetOffsetEnCorr = sfNoPUjetOffsetEnCorr,
                addToPatDefaultSequence = False,
                postfix = ''
        )
	process.calibratedAK5PFJetsForNoPileUpPFMEt.correctors = cms.vstring('ak5PFL1FastL2L3')
	process.noPileUpPFMEt.srcLeptons = cms.VInputTag('selectElectrons')

	#from PhysicsTools.PatUtils.tools.runMVAMEtUncertainties         import runMVAMEtUncertainties
        #runMVAMEtUncertainties(process,
        #        electronCollection = cms.InputTag('selectElectrons'),
        #        photonCollection = '',
        #        muonCollection = '',
        #        tauCollection = '',
        #        jetCollection = cms.InputTag('patJets'),
	#	jecUncertaintyFile = "PhysicsTools/PatUtils/data/Summer13_V1_DATA_UncertaintySources_AK5PF.txt",
        #        doSmearJets = True,
        #        addToPatDefaultSequence = False,
        #        postfix = ''
        #)
        #process.calibratedAK5PFJetsForPFMEtMVA.correctors = cms.vstring('ak5PFL1FastL2L3')
        #process.pfMEtMVA.srcLeptons = cms.VInputTag('selectElectrons')
 
#	from PhysicsTools.PatUtils.tools.runType1CaloMEtUncertainties   import runType1CaloMEtUncertainties
#        runType1CaloMEtUncertainties(process,
#                electronCollection = cms.InputTag('selectElectrons'),
#                photonCollection = '',
#                muonCollection = '',
#                tauCollection = '',
#                jetCollection = cms.InputTag('patJets'),
#		caloTowerCollection = cms.InputTag('towerMaker'),
#                jetCorrPayloadName = 'AK5Calo',
#                jetCorrLabelUpToL3 = 'ak5CaloL1L2L3',
#                jetCorrLabelUpToL3Res = 'ak5CaloL1L2L3Residual',
#		jecUncertaintyFile = "PhysicsTools/PatUtils/data/Fall12_V7_DATA_UncertaintySources_AK5Calo.txt",
#		jecUncertaintyTag = 'SubTotalMC',
#                addToPatDefaultSequence = False,
#                postfix = ''
#        )

###--------------------------------------------------------------
# MVA MET unity response
#process.pfMEtMVAUnity = process.pfMEtMVA.clone(
#	inputFileNames = cms.PSet(
#	U     = cms.FileInPath('JetMETCorrections/METPUSubtraction/data/gbrmet_53_Sep2013_type1_UnityResponse_v3.root'),
#	DPhi  = cms.FileInPath('JetMETCorrections/METPUSubtraction/data/gbrmetphi_53_June2013_type1.root'),
#	CovU1 = cms.FileInPath('JetMETCorrections/METPUSubtraction/data/gbru1cov_53_Dec2012_UnityResponse.root'),
#	CovU2 = cms.FileInPath('JetMETCorrections/METPUSubtraction/data/gbru2cov_53_Dec2012_UnityResponse.root')
#))
#process.patPFMetMVAUnity = process.patPFMetMVA.clone(
#	metSource = cms.InputTag("pfMEtMVAUnity")
#)
#process.pfMEtMVAUnityElectronEnDown = process.pfMEtMVAElectronEnDown.clone(
#	inputFileNames = cms.PSet(
#	U     = cms.FileInPath('JetMETCorrections/METPUSubtraction/data/gbrmet_53_Sep2013_type1_UnityResponse_v3.root'),
#	DPhi  = cms.FileInPath('JetMETCorrections/METPUSubtraction/data/gbrmetphi_53_June2013_type1.root'),
#	CovU1 = cms.FileInPath('JetMETCorrections/METPUSubtraction/data/gbru1cov_53_Dec2012_UnityResponse.root'),
#	CovU2 = cms.FileInPath('JetMETCorrections/METPUSubtraction/data/gbru2cov_53_Dec2012_UnityResponse.root')
#))
#process.patPFMetMVAUnityElectronEnUp = process.patPFMetMVAElectronEnUp.clone(
#        metSource = cms.InputTag("pfMEtMVAUnityElectronEnUp")
#)
#process.pfMEtMVAUnityElectronEnUp = process.pfMEtMVAElectronEnUp.clone(
#	inputFileNames = cms.PSet(
#	U     = cms.FileInPath('JetMETCorrections/METPUSubtraction/data/gbrmet_53_Sep2013_type1_UnityResponse_v3.root'),
#	DPhi  = cms.FileInPath('JetMETCorrections/METPUSubtraction/data/gbrmetphi_53_June2013_type1.root'),
#	CovU1 = cms.FileInPath('JetMETCorrections/METPUSubtraction/data/gbru1cov_53_Dec2012_UnityResponse.root'),
#	CovU2 = cms.FileInPath('JetMETCorrections/METPUSubtraction/data/gbru2cov_53_Dec2012_UnityResponse.root')
#))
#process.patPFMetMVAUnityElectronEnDown = process.patPFMetMVAElectronEnDown.clone(
#        metSource = cms.InputTag("pfMEtMVAUnityElectronEnDown")
#)
#process.pfMEtMVAUnityJetEnDown = process.pfMEtMVAJetEnDown.clone(
#	inputFileNames = cms.PSet(
#	U     = cms.FileInPath('JetMETCorrections/METPUSubtraction/data/gbrmet_53_Sep2013_type1_UnityResponse_v3.root'),
#	DPhi  = cms.FileInPath('JetMETCorrections/METPUSubtraction/data/gbrmetphi_53_June2013_type1.root'),
#	CovU1 = cms.FileInPath('JetMETCorrections/METPUSubtraction/data/gbru1cov_53_Dec2012_UnityResponse.root'),
#	CovU2 = cms.FileInPath('JetMETCorrections/METPUSubtraction/data/gbru2cov_53_Dec2012_UnityResponse.root')
#))
#process.patPFMetMVAUnityJetEnDown = process.patPFMetMVAJetEnDown.clone(
#	metSource = cms.InputTag("pfMEtMVAUnityJetEnDown")
#)
#process.pfMEtMVAUnityJetEnUp = process.pfMEtMVAJetEnUp.clone(
#	inputFileNames = cms.PSet(
#	U     = cms.FileInPath('JetMETCorrections/METPUSubtraction/data/gbrmet_53_Sep2013_type1_UnityResponse_v3.root'),
#	DPhi  = cms.FileInPath('JetMETCorrections/METPUSubtraction/data/gbrmetphi_53_June2013_type1.root'),
#	CovU1 = cms.FileInPath('JetMETCorrections/METPUSubtraction/data/gbru1cov_53_Dec2012_UnityResponse.root'),
#	CovU2 = cms.FileInPath('JetMETCorrections/METPUSubtraction/data/gbru2cov_53_Dec2012_UnityResponse.root')
#))
#process.patPFMetMVAUnityJetEnUp = process.patPFMetMVAJetEnUp.clone(
#        metSource = cms.InputTag("pfMEtMVAUnityJetEnUp")
#)
#process.pfMEtMVAUnityUnclusteredEnDown = process.pfMEtMVAUnclusteredEnDown.clone(
#	inputFileNames = cms.PSet(
#	U     = cms.FileInPath('JetMETCorrections/METPUSubtraction/data/gbrmet_53_Sep2013_type1_UnityResponse_v3.root'),
#	DPhi  = cms.FileInPath('JetMETCorrections/METPUSubtraction/data/gbrmetphi_53_June2013_type1.root'),
#	CovU1 = cms.FileInPath('JetMETCorrections/METPUSubtraction/data/gbru1cov_53_Dec2012_UnityResponse.root'),
#	CovU2 = cms.FileInPath('JetMETCorrections/METPUSubtraction/data/gbru2cov_53_Dec2012_UnityResponse.root')
#))
#process.patPFMetMVAUnityUnclusteredEnDown = process.patPFMetMVAUnclusteredEnDown.clone(
#        metSource = cms.InputTag("pfMEtMVAUnityUnclusteredEnDown")
#)
#process.pfMEtMVAUnityUnclusteredEnUp = process.pfMEtMVAUnclusteredEnUp.clone(
#        inputFileNames = cms.PSet(
#        U     = cms.FileInPath('JetMETCorrections/METPUSubtraction/data/gbrmet_53_Sep2013_type1_UnityResponse_v3.root'),
#        DPhi  = cms.FileInPath('JetMETCorrections/METPUSubtraction/data/gbrmetphi_53_June2013_type1.root'),
#        CovU1 = cms.FileInPath('JetMETCorrections/METPUSubtraction/data/gbru1cov_53_Dec2012_UnityResponse.root'),
#        CovU2 = cms.FileInPath('JetMETCorrections/METPUSubtraction/data/gbru2cov_53_Dec2012_UnityResponse.root')
#))
#process.patPFMetMVAUnityUnclusteredEnUp = process.patPFMetMVAUnclusteredEnUp.clone(
#        metSource = cms.InputTag("pfMEtMVAUnityUnclusteredEnUp")
#)
#if iRunOnData == False:
#	process.pfMEtMVAUnityJetResDown = process.pfMEtMVAJetResDown.clone(
#        	inputFileNames = cms.PSet(
#        	U     = cms.FileInPath('JetMETCorrections/METPUSubtraction/data/gbrmet_53_Sep2013_type1_UnityResponse_v3.root'),
#        	DPhi  = cms.FileInPath('JetMETCorrections/METPUSubtraction/data/gbrmetphi_53_June2013_type1.root'),
#        	CovU1 = cms.FileInPath('JetMETCorrections/METPUSubtraction/data/gbru1cov_53_Dec2012_UnityResponse.root'),
#        	CovU2 = cms.FileInPath('JetMETCorrections/METPUSubtraction/data/gbru2cov_53_Dec2012_UnityResponse.root')
#	))
#	process.patPFMetMVAUnityJetResDown = process.patPFMetMVAJetResDown.clone(
#		metSource = cms.InputTag("pfMEtMVAUnityJetResDown")
#	)
#        process.pfMEtMVAUnityJetResUp = process.pfMEtMVAJetResUp.clone(
#                inputFileNames = cms.PSet(
#                U     = cms.FileInPath('JetMETCorrections/METPUSubtraction/data/gbrmet_53_Sep2013_type1_UnityResponse_v3.root'),
#                DPhi  = cms.FileInPath('JetMETCorrections/METPUSubtraction/data/gbrmetphi_53_June2013_type1.root'),
#                CovU1 = cms.FileInPath('JetMETCorrections/METPUSubtraction/data/gbru1cov_53_Dec2012_UnityResponse.root'),
#                CovU2 = cms.FileInPath('JetMETCorrections/METPUSubtraction/data/gbru2cov_53_Dec2012_UnityResponse.root')
#        ))
#	process.patPFMetMVAUnityJetResUp = process.patPFMetMVAJetResUp.clone(
#		metSource = cms.InputTag("pfMEtMVAUnityJetResUp")
#	)
#	process.pfMVAMEtUnityUncertaintySequence = cms.Sequence(process.pfMEtMVAUnity + process.patPFMetMVAUnity + process.pfMEtMVAUnityElectronEnDown + process.patPFMetMVAUnityElectronEnDown + process.pfMEtMVAUnityElectronEnUp + process.patPFMetMVAUnityElectronEnUp + process.pfMEtMVAUnityJetResDown + process.patPFMetMVAUnityJetResDown + process.pfMEtMVAUnityJetResUp + process.patPFMetMVAUnityJetResUp + process.pfMEtMVAUnityJetEnUp + process.patPFMetMVAUnityJetEnUp + process.pfMEtMVAUnityJetEnDown + process.patPFMetMVAUnityJetEnDown + process.pfMEtMVAUnityUnclusteredEnUp + process.patPFMetMVAUnityUnclusteredEnUp + process.pfMEtMVAUnityUnclusteredEnDown + process.patPFMetMVAUnityUnclusteredEnDown)
#else : 
#	process.pfMVAMEtUnityUncertaintySequence = cms.Sequence(process.pfMEtMVAUnity + process.patPFMetMVAUnity + process.pfMEtMVAUnityElectronEnDown + process.patPFMetMVAUnityElectronEnDown + process.pfMEtMVAUnityElectronEnUp + process.patPFMetMVAUnityElectronEnUp + process.pfMEtMVAUnityJetEnUp + process.patPFMetMVAUnityJetEnUp + process.pfMEtMVAUnityJetEnDown + process.patPFMetMVAUnityJetEnDown + process.pfMEtMVAUnityUnclusteredEnUp + process.patPFMetMVAUnityUnclusteredEnUp + process.pfMEtMVAUnityUnclusteredEnDown + process.patPFMetMVAUnityUnclusteredEnDown)
 
###--------------------------------------------------------------
# add CaloMET and Type1 CaloMET
process.patCaloMet = process.patMETs.clone(
    metSource = cms.InputTag('met')
)
process.load("JetMETCorrections/Type1MET/caloMETCorrections_cff")
process.caloJetMETcorr.srcMET = cms.InputTag('met')
process.caloJetMETcorr.offsetCorrLabel = cms.string('ak5CaloL1Fastjet')
if iRunOnData == False:
    process.caloJetMETcorr.jetCorrLabel = cms.string("ak5CaloL1FastL2L3")
else:
    process.caloJetMETcorr.jetCorrLabel = cms.string("ak5CaloL1FastL2L3Residual")
process.caloType1CorrectedMet.src = cms.InputTag('corMetGlobalMuons')
process.patType1CorrectedCaloMet = process.patMETs.clone(
    metSource = cms.InputTag('caloType1CorrectedMet')
)


###--------------------------------------------------------------
# Z to ee Event Selection
process.load('ElectroWeakAnalysis.VPlusJets.ZeeCollectionsPAT_cfi')
# Ntuples
process.load('ElectroWeakAnalysis.VPlusJets.ZeeJetsAnalysisPATNoPU_cfi')
process.MyZeeJets 			= process.VplusJetsNtuple.clone()
process.MyZeeJets.runningOverMC 	= cms.bool(isMC)  
if iRunOnData == True:
	process.MyZeeJets.srcGenParticles     	= cms.InputTag("")
else:
	process.MyZeeJets.srcGenParticles       = cms.InputTag("genParticles")

if iRunOnData == True:
	process.load("CommonTools.UtilAlgos.TFileService_cfi")
	process.TFileService.fileName = cms.string('Zee_53X_AOD_PAT.root')
else:
        process.load("CommonTools.UtilAlgos.TFileService_cfi")
        process.TFileService.fileName = cms.string('Zee_53X_SIM_PAT.root')

###--------------------------------------------------------------
### Output
process.out.outputCommands += [ #'drop *']
                               'keep *_patJets*_*_*',
                               'keep *_selectedPat*_*_*',
                               'keep *_cleanPat*_*_*',
                               'keep *_selectElectrons_*_*',
                               'keep double_*_rho_*',
                               'keep *_patMETs*_*_*',
                               'keep *_patMETsPF*_*_*',
                               'keep *_patPFMet*_*_*',
                               'keep *_patType1CorrectedPFMet*_*_*',
                               'keep *_patType1p2CorrectedPFMet*_*_*',
                               #'keep *_patPFMetMVA*_*_*',
                               'keep *_patPFMetNoPileUp*_*_*',
                               #'keep CommonMETData_noPileUpPFMEt*_*_*',
                               'keep *_patCaloMet*_*_*',
                               'keep *_patType1CorrectedCaloMet*_*_*',
                                # vertices and beam spot
                               'keep *_offlineBeamSpot_*_*',
                               'keep *_offlinePrimaryVertices*_*_*',
                               'keep *_goodOfflinePrimaryVertices*_*_*']
if iRunOnData == False:
        process.out.outputCommands += [
                               'keep GenEventInfoProduct_*_*_*',
			       'keep smearedPatJets_*_*_*',
                               'keep recoGenParticles_*_*_*',
                               'keep GenMETs_*_*_*',
                               'keep *_addPileupInfo_*_*',
                               'keep LHEEventProduct_*_*_*']


###--------------------------------------------------------------
### the path
if iRunOnData == True:
	process.p = cms.Path(
	# trigger filter
    	process.hltHighLevel *
	# basic filters
    	process.noscraping *
    	process.primaryVertexFilter *
	# MET filters
    	process.HBHENoiseFilter  * 
    	process.CSCTightHaloFilter *  
    	process.hcalLaserEventFilter  *
    	process.hcalfilter *
    	process.EcalDeadCellTriggerPrimitiveFilter *
    	process.eeBadScFilter *
    	process.ecalLaserCorrFilter *
    	process.goodVertices *
    	process.trackingFailureFilter *
    	process.trkPOGFilters *
    	# PAT
    	#process.pfNoPileUpSequence *
    	process.pfParticleSelectionSequence *
    	process.patDefaultSequence *
    	process.selectElectronsseq *
    	process.pfType1MEtUncertaintySequence *
    	process.pfNoPileUpMEtUncertaintySequence *
    	process.pfMVAMEtUncertaintySequence *
    	#process.caloType1MEtUncertaintySequence 
	process.pfMVAMEtUnityUncertaintySequence * 
	# CaloMET
    	process.patCaloMet *
    	process.caloJetMETcorr *
    	process.caloType1CorrectedMet *
    	process.patType1CorrectedCaloMet *
	# NTuple
    	process.ZPath *
    	process.MyZeeJets
	)
else:
        process.p = cms.Path(
        # trigger filter
        process.hltHighLevel *
        # basic filters
        process.noscraping *
        process.primaryVertexFilter *
        # PAT
        #process.pfNoPileUpSequence *
        process.pfParticleSelectionSequence *
        process.patDefaultSequence *
        process.selectElectronsseq *
        process.pfType1MEtUncertaintySequence *
        process.pfNoPileUpMEtUncertaintySequence *
        #process.pfMVAMEtUncertaintySequence *
        #process.caloType1MEtUncertaintySequence 
	#process.pfMVAMEtUnityUncertaintySequence *
        # CaloMET
        process.patCaloMet *
        process.caloJetMETcorr *
        process.caloType1CorrectedMet *
        process.patType1CorrectedCaloMet * 
        # NTuple
        process.ZPath *
        process.MyZeeJets
        ) 

###--------------------------------------------------------------
# Output file
if iRunOnData == True:
	process.out.fileName = 'patTuple_data.root'
else:
	process.out.fileName = 'patTuple_mc.root'
 
del(process.out)
del(process.outpath)

####-- Dump config
#------------------------------------------------------------
#file = open('ZeePAT_MC.py','w')
#file.write(str(process.dumpPython()))
#file.close()
 
if iRunOnData == True:
	processDumpFile = open('python_config_data.dump', 'w')
	print >> processDumpFile, process.dumpPython()
else:
	processDumpFile = open('python_config_mc.dump', 'w')
	print >> processDumpFile, process.dumpPython()


