import FWCore.ParameterSet.Config as cms

ZToEE = cms.EDProducer("NamedCandViewShallowCloneCombiner",
    cut = cms.string('60 < mass < 120'),
    name = cms.string('ZToEE'),
    roles = cms.vstring('electron1', 'electron2'),
    decay = cms.string('selectElectrons@+ selectElectrons@-')
)

bestZee = cms.EDFilter("LargestPtCandViewSelector",
    maxNumber = cms.uint32(1),
    src = cms.InputTag("ZToEE")
)

ZPath = cms.Sequence(ZToEE * bestZee)

