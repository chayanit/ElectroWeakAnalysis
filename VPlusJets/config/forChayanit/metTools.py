import FWCore.ParameterSet.Config as cms

from PhysicsTools.PatAlgos.tools.helpers import *

def addPFMet(process, correct = False):    
    process.load("JetMETCorrections.Type1MET.MetType1Corrections_cff")    
    process.metJESCorAK5PFJet.jetPTthreshold = cms.double(10.0)
    process.metJESCorAK5PFJet.useTypeII = cms.bool(True)

    process.patPFMETs = process.patMETs.clone()
    process.patPFMETs.addMuonCorrections = False
    process.patPFMETs.genMETSource = cms.InputTag('genMetTrue')

    if correct:
        process.patPFMETs.metSource = cms.InputTag('metJESCorAK5PFJet')
        process.makePatPFMETs = cms.Sequence(process.metJESCorAK5PFJet * process.patPFMETs)
    else:
        process.patPFMETs.metSource = cms.InputTag('pfMet')
        process.makePatPFMETs = cms.Sequence(process.patPFMETs)

    process.makePatMETs += process.makePatPFMETs
    
    return process.makePatMETs

def addCorrectedPFMet(process, isMC, doApplyType0corr, doApplySysShiftCorr, runPeriod, doSmearJets, jecUncertaintyTag, doApplyUnclEnergyResidualCorr):

    process.load("PhysicsTools.PatUtils.patPFMETCorrections_cff")
    
    process.load("JetMETCorrections.METPUSubtraction.mvaPFMET_cff")
    if isMC:
        process.calibratedAK5PFJetsForPFMEtMVA.correctors = cms.vstring("ak5PFL1FastL2L3")
    else:
        process.calibratedAK5PFJetsForPFMEtMVA.correctors = cms.vstring("ak5PFL1FastL2L3Residual")
    process.pfMEtMVA.srcCorrJets = cms.InputTag('calibratedAK5PFJetsForPFMEtMVA')
    process.pfMEtMVA.srcLeptons = cms.VInputTag('goodMuons')
    process.pfMEtMVA.inputFileNames = cms.PSet(
        DPhi = cms.FileInPath('JetMETCorrections/METPUSubtraction/data/gbrmetphi_53_Dec2012.root'),
        CovU2 = cms.FileInPath('JetMETCorrections/METPUSubtraction/data/gbru2cov_53_Dec2012.root'),
        U = cms.FileInPath('JetMETCorrections/METPUSubtraction/data/gbrmet_53_Dec2012.root'),
        CovU1 = cms.FileInPath('JetMETCorrections/METPUSubtraction/data/gbru1cov_53_Dec2012.root')
    )
    process.pfMEtMVA.verbosity = cms.int32(0)
    process.patPFMetMVA = process.patMETs.clone(
        metSource = cms.InputTag('pfMEtMVA'),
        addMuonCorrections = cms.bool(False),
        genMETSource = cms.InputTag('genMetTrue')
    )
    process.patMEtMVAsequence = cms.Sequence(process.pfMEtMVAsequence + process.patPFMetMVA)

    process.load("JetMETCorrections.METPUSubtraction.noPileUpPFMET_cff")
    process.load("JetMETCorrections.METPUSubtraction.noPileUpPFchsMET_cff")
    if isMC:
        process.calibratedAK5PFJetsForNoPileUpPFMEt.correctors = cms.vstring("ak5PFL1FastL2L3")
        process.calibratedAK5PFchsJetsForNoPileUpPFchsMEt.correctors = cms.vstring("ak5PFchsL1FastL2L3")
    else:
        process.calibratedAK5PFJetsForNoPileUpPFMEt.correctors = cms.vstring("ak5PFL1FastL2L3Residual")
        process.calibratedAK5PFchsJetsForNoPileUpPFchsMEt.correctors = cms.vstring("ak5PFchsL1FastL2L3Residual")

    process.noPileUpPFMEt.srcLeptons = cms.VInputTag('patMuons')
    process.noPileUpPFMEtData.verbosity = cms.int32(0)
    process.noPileUpPFMEt.verbosity = cms.int32(0)
    process.patPFMetNoPileUp = process.patMETs.clone(
        metSource = cms.InputTag('noPileUpPFMEt'),
        addMuonCorrections = cms.bool(False),
        genMETSource = cms.InputTag('genMetTrue')
    )
    process.noPileUpPFchsMEt.srcLeptons = cms.VInputTag('patMuons')
    process.noPileUpPFchsMEtData.verbosity = cms.int32(0)
    process.noPileUpPFchsMEt.verbosity = cms.int32(0)
    process.patPFchsMetNoPileUp = process.patMETs.clone(
        metSource = cms.InputTag('noPileUpPFchsMEt'),
        addMuonCorrections = cms.bool(False),
        genMETSource = cms.InputTag('genMetTrue')
    )
    process.patMEtNoPileUpSequence = cms.Sequence(process.noPileUpPFMEtSequence + process.patPFMetNoPileUp + process.noPileUpPFchsMEtSequence + process.patPFchsMetNoPileUp)
    ##process.patMEtNoPileUpSequence = cms.Sequence(process.noPileUpPFMEtSequence + process.patPFMetNoPileUp)
    
    process.makeCorrectedPatPFMETs = cms.Sequence()

    if isMC:
        import PhysicsTools.PatAlgos.tools.helpers as configtools
        process.type0PFMEtCorrection.remove(process.type0PFMEtCorrectionPFCandToVertexAssociation)
        process.makeCorrectedPatPFMETs += process.type0PFMEtCorrectionPFCandToVertexAssociation
        configtools.cloneProcessingSnippet(process, process.producePatPFMETCorrections, "NoSmearing")
        process.selectedPatJetsForMETtype1p2CorrNoSmearing.src = cms.InputTag('patJetsNotOverlappingWithLeptonsForJetMEtUncertainty')
        process.selectedPatJetsForMETtype2CorrNoSmearing.src = process.selectedPatJetsForMETtype1p2CorrNoSmearing.src
        configtools.cloneProcessingSnippet(process, process.patMEtMVAsequence, "NoSmearing")
        process.makeCorrectedPatPFMETs += process.patMEtMVAsequenceNoSmearing
        ##process.patMEtNoPileUpSequence.remove(process.type0PFMEtCorrection)
        configtools.cloneProcessingSnippet(process, process.patMEtNoPileUpSequence, "NoSmearing")
        process.makeCorrectedPatPFMETs += process.patMEtNoPileUpSequenceNoSmearing
    else:
        doSmearJets = False
        
    sysShiftCorrParameter = None
    if doApplySysShiftCorr:
        process.load("JetMETCorrections.Type1MET.pfMETsysShiftCorrections_cfi")
        if runPeriod == "2012RunABCD":
            if isMC:
                sysShiftCorrParameter = process.pfMEtSysShiftCorrParameters_2012runABCvsNvtx_mc
            else:
                sysShiftCorrParameter = process.pfMEtSysShiftCorrParameters_2012runABCvsNvtx_data
        else:
            raise ValueError("Invalid runPeriod = %s !!" % runPeriod)

    from PhysicsTools.PatUtils.tools.runType1PFMEtUncertainties import runType1PFMEtUncertainties
    runType1PFMEtUncertainties(
        process,
        electronCollection = '',
        photonCollection = '',
        muonCollection = cms.InputTag('patMuons'),
        tauCollection = '',
        jetCollection = cms.InputTag('patJets'),        
        doSmearJets = doSmearJets,
        jecUncertaintyTag = jecUncertaintyTag,
        makeType1corrPFMEt = True,
        makeType1p2corrPFMEt = True,
        doApplyType0corr = doApplyType0corr,
        sysShiftCorrParameter = sysShiftCorrParameter,
        doApplySysShiftCorr = doApplySysShiftCorr,
        doApplyUnclEnergyCalibration = (doApplyUnclEnergyResidualCorr & isMC),
        addToPatDefaultSequence = False
    )
    from PhysicsTools.PatUtils.tools.runMVAMEtUncertainties import runMVAMEtUncertainties
    runMVAMEtUncertainties(
        process,
        electronCollection = '',
        photonCollection = '',
        muonCollection = cms.InputTag('patMuons'),        
        tauCollection = '',        
        jetCollection = cms.InputTag('patJets'),   
        doSmearJets = doSmearJets,
        jecUncertaintyTag = jecUncertaintyTag,
        addToPatDefaultSequence = False
    )
    from PhysicsTools.PatUtils.tools.runNoPileUpMEtUncertainties import runNoPileUpMEtUncertainties
    runNoPileUpMEtUncertainties(
        process,
        electronCollection = '',
        photonCollection = '',
        muonCollection = cms.InputTag('patMuons'),
        tauCollection = '',
        jetCollection = cms.InputTag('patJets'),    
        doApplyChargedHadronSubtraction = False,
        doSmearJets = doSmearJets,
        jecUncertaintyTag = jecUncertaintyTag,
        doApplyUnclEnergyCalibration = (doApplyUnclEnergyResidualCorr & isMC),
        addToPatDefaultSequence = False
    )
    runNoPileUpMEtUncertainties(
        process,
        electronCollection = '',
        photonCollection = '',
        muonCollection = cms.InputTag('patMuons'),
        tauCollection = '',
        jetCollection = cms.InputTag('patJetsAK5PFchs'),    
        doApplyChargedHadronSubtraction = True,
        doSmearJets = doSmearJets,
        jecUncertaintyFile = "PhysicsTools/PatUtils/data/Fall12_V7_DATA_UncertaintySources_AK5PFchs.txt",
        jecUncertaintyTag = jecUncertaintyTag,
        addToPatDefaultSequence = False
    )

    if isMC:
        process.makeCorrectedPatPFMETs += process.pfType1MEtUncertaintySequence
        process.makeCorrectedPatPFMETs += process.pfNoPileUpMEtUncertaintySequence
        process.makeCorrectedPatPFMETs += process.pfchsNoPileUpMEtUncertaintySequence
        process.makeCorrectedPatPFMETs += process.pfMVAMEtUncertaintySequence
        process.patPFMet.addGenMET = cms.bool(True)
        process.patPFMetMVA.addGenMET = cms.bool(True)
        process.patPFJetMETtype1p2Corr.jetCorrLabel = cms.string("L3Absolute")
        process.patPFJetMETtype1p2CorrNoSmearing.jetCorrLabel = cms.string("L3Absolute")
    else:
        process.patPFMet.addGenMET = cms.bool(False)
        process.patPFMetMVA.addGenMET = cms.bool(False)
        process.patPFJetMETtype1p2Corr.jetCorrLabel = cms.string("L2L3Residual")
    
        process.makeCorrectedPatPFMETs += process.patJetsNotOverlappingWithLeptonsForJetMEtUncertainty
        if hasattr(process, "pfMEtSysShiftCorrSequence"):
            process.makeCorrectedPatPFMETs += process.pfMEtSysShiftCorrSequence
        process.makeCorrectedPatPFMETs += process.producePatPFMETCorrections
        process.makeCorrectedPatPFMETs += process.patMEtMVAsequence
        process.makeCorrectedPatPFMETs += process.patMEtNoPileUpSequence        

    # add MVA MEt with unity response training
    for moduleName in dir(process):
        if (moduleName.endswith("Up") or moduleName.endswith("Down")) and not isMC:
            continue
        module = getattr(process, moduleName)
        if isinstance(module, cms.EDProducer) and module.type_() == "PFMETProducerMVA":
            module_unity = module.clone(
                inputFileNames = cms.PSet(
                    DPhi = cms.FileInPath('JetMETCorrections/METPUSubtraction/data/gbrmetphi_53_Dec2012.root'), # CV: same for unity and non-unity response training
                    CovU2 = cms.FileInPath('JetMETCorrections/METPUSubtraction/data/gbru2cov_53_Dec2012_UnityResponse.root'),
                    U = cms.FileInPath('JetMETCorrections/METPUSubtraction/data/gbrmet_53_Dec2012_UnityResponse.root'),
                    CovU1 = cms.FileInPath('JetMETCorrections/METPUSubtraction/data/gbru1cov_53_Dec2012_UnityResponse.root')
                )
            )
            moduleName_unity = moduleName.replace("pfMEtMVA", "pfMEtMVAunityResponse")
            setattr(process, moduleName_unity, module_unity)
            process.makeCorrectedPatPFMETs += module_unity    
    for moduleName in dir(process):
        if (moduleName.endswith("Up") or moduleName.endswith("Down")) and not isMC:
            continue
        module = getattr(process, moduleName)
        if isinstance(module, cms.EDProducer) and module.type_() == "PATMETProducer" and moduleName.find("patPFMetMVA") != -1:
            module_unity = module.clone(
                metSource = cms.InputTag(module.metSource.value().replace("pfMEtMVA", "pfMEtMVAunityResponse"))
            )
            moduleName_unity = moduleName.replace("patPFMetMVA", "patPFMetMVAunityResponse")
            setattr(process, moduleName_unity, module_unity)
            process.makeCorrectedPatPFMETs += module_unity

    return process.makeCorrectedPatPFMETs

def addCorrectedCaloMet(process, isMC, jecUncertaintyTag, doApplyUnclEnergyResidualCorr):
    
    process.makeCorrectedPatCaloMETs = cms.Sequence()

    from PhysicsTools.PatUtils.tools.runType1CaloMEtUncertainties import runType1CaloMEtUncertainties
    runType1CaloMEtUncertainties(
        process,
        electronCollection = '',
        photonCollection = '',
        muonCollection = cms.InputTag('patMuons'),
        tauCollection = '',
        jetCollection = cms.InputTag('patJetsAK5Calo'),
        caloTowerCollection = cms.InputTag('towerMaker'),
        jetCorrPayloadName = "AK5Calo",
        jetCorrLabelUpToL3 = "ak5CaloL1L2L3",
        jetCorrLabelUpToL3Res = "ak5CaloL1L2L3Residual",
        jecUncertaintyFile = "PhysicsTools/PatUtils/data/Fall12_V7_DATA_UncertaintySources_AK5Calo.txt",
        jecUncertaintyTag = jecUncertaintyTag,
        addToPatDefaultSequence = False
    )

    if isMC:
        process.makeCorrectedPatCaloMETs += process.caloType1MEtUncertaintySequence
    else:
        if doApplyUnclEnergyResidualCorr:
            process.caloJetMETcorr.type2ResidualCorrLabel = cms.string("ak5CaloResidual")
            process.caloJetMETcorr.type2ResidualCorrEtaMax = cms.double(9.9)
            process.caloJetMETcorr.type2ResidualCorrOffset = cms.double(1.)
            process.caloJetMETcorr.type2ExtraCorrFactor = cms.double(1.05)
            process.caloTowersNotInJets = cms.EDProducer("TPCaloJetsOnCaloTowers",
                enable = cms.bool(True),
                verbose = cms.untracked.bool(False),
                name = cms.untracked.string("caloTowersNotInJets"),
                topCollection = cms.InputTag('ak5CaloJets'),
                bottomCollection = cms.InputTag('towerMaker')
            )
            process.residualCaloMEtUnclusteredEnCorr = cms.EDProducer("CaloTowerMETcorrInputProducer",
                src = cms.InputTag('caloTowersNotInJets'),
                residualCorrLabel = cms.string("ak5CaloResidual"),
                residualCorrEtaMax = cms.double(9.9),
                residualCorrOffset = cms.double(1.),
                isMC = cms.bool(False), # CV: only used to decide whether to apply "unclustered energy" calibration to MC or Data    
                extraCorrFactor = cms.double(1.05),
                globalThreshold = cms.double(0.3), # NOTE: this value need to match met.globalThreshold, defined in RecoMET/METProducers/python/CaloMET_cfi.py
                noHF = cms.bool(False)
            )            
            process.produceCaloMETCorrections.replace(process.caloJetMETcorr, process.caloTowersNotInJets + process.caloJetMETcorr + process.residualCaloMEtUnclusteredEnCorr)
            ##process.caloType1CorrectedMet.srcType1Corrections.append(cms.InputTag('residualCaloMEtUnclusteredEnCorr'))
            process.caloType1CorrectedMet.applyType2Corrections = cms.bool(True)
            process.caloType1CorrectedMet.srcUnclEnergySums = cms.VInputTag(
                cms.InputTag('residualCaloMEtUnclusteredEnCorr'),
                cms.InputTag("caloJetMETcorr", "type2")
            )
            process.caloType1CorrectedMet.type2CorrFormula = cms.string("A")
            process.caloType1CorrectedMet.type2CorrParameter = cms.PSet(A = cms.double(2.))
            ##process.caloType1CorrectedMet.verbosity = cms.int32(1)
    
    return process.makeCorrectedPatCaloMETs

def addTCMet(process):
    process.layer1TCMETs = process.patMETs.clone()
    process.layer1TCMETs.addMuonCorrections = False
    process.layer1TCMETs.metSource = cms.InputTag('tcMet')
    process.layer1TCMETs.genMETSource = cms.InputTag('genMETWithMu')
    process.patCandidates.replace(process.patMETs,
                                  process.patMETs + process.layer1TCMETs)
    
def replaceMETforDiTaus(process,
                        oldMet = cms.InputTag('patMETs'),
                        newMet = cms.InputTag('patPFMETs') ):
    massSearchReplaceParam(process.produceDiTauPairsAllKinds,
                           'srcMET', oldMet, newMet)

def replaceMETforTauNu(process,
                       oldMet = cms.InputTag('patMETs'),
                       newMet = cms.InputTag('patPFMETs') ):
    massSearchReplaceParam(process.produceTauNuPairs,
                           'srcMET', oldMet, newMet)
def replaceMETforMet(process,
                     oldMet = cms.InputTag('Layer1METs'),
                     newMet = cms.InputTag('patPFMETs')):
    massSearchReplaceParam(process.selectLayer1METs,
                           'src', oldMet ,newMet)


