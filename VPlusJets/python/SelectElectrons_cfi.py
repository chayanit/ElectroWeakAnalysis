import FWCore.ParameterSet.Config as cms

tightEleIdLabel = "tight"
mediumEleIdLabel = "medium"
looseEleIdLabel = "loose"
vetoEleIdLabel = "veto"

selectElectrons = cms.EDProducer("PATElectronIdSelector",
    src = cms.InputTag( "cleanPatElectrons" ),
    idLabel = cms.string(looseEleIdLabel),
)

selectElectronsseq = cms.Sequence(selectElectrons)

