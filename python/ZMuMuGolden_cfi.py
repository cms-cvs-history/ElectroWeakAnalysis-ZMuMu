import FWCore.ParameterSet.Config as cms

import copy

###################################################
#              muons for ZMuMu                    #    
###################################################

goodAODGlobalMuons = cms.EDFilter("MuonViewRefSelector",
  src = cms.InputTag("muons"),
  cut = cms.string('isGlobalMuon = 1 & pt > 20 & abs(eta)<2.1 & isolationR03().sumPt<3.0'),
  filter = cms.bool(True)                                
)

###################################################
#              combiner module                    #    
###################################################

dimuonsGlobalAOD = cms.EDFilter("CandViewShallowCloneCombiner",
    checkCharge = cms.bool(True),
    cut = cms.string('mass > 20 & charge=0'),
    decay = cms.string("goodAODGlobalMuons@+ goodAODGlobalMuons@-")
)


# dimuon filter
dimuonsFilter = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag("dimuonsGlobalAOD"),
    minNumber = cms.uint32(1)
)

import HLTrigger.HLTfilters.hltHighLevel_cfi

dimuonsHLTFilter = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
# Add this to access 8E29 menu
#dimuonsHLTFilter.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT8E29")
# for 8E29 menu
#dimuonsHLTFilter.HLTPaths = ["HLT_Mu3", "HLT_DoubleMu3"]
# for 1E31 menu
#dimuonsHLTFilter.HLTPaths = ["HLT_Mu9", "HLT_DoubleMu3"]
dimuonsHLTFilter.HLTPaths = ["HLT_Mu9"]




##################################################
#####                selection             #######
##################################################

zSelection = cms.PSet(
## cut already implemented, but one could add more (e.g. massMin, massMax,... change the pt or eta cut....)
    cut = cms.string("charge = 0 & daughter(0).pt > 20 & daughter(1).pt > 20 & abs(daughter(0).eta)<2.1 & abs(daughter(1).eta)<2.1 & mass > 20"),
    )


#ZMuMu: at least one HLT trigger match
goodZToMuMuAtLeast1HLT = cms.EDFilter(
    "ZGoldenSelectorAndFilter",
    zSelection,
    TrigTag = cms.InputTag("TriggerResults::HLT"),
    triggerEvent = cms.InputTag( "hltTriggerSummaryAOD::HLT" ),
    src = cms.InputTag("dimuonsGlobalAOD"),
    condition =cms.string("atLeastOneMatched"),
    hltPath = cms.string("HLT_Mu9"),
    L3FilterName= cms.string("hltSingleMu9L3Filtered9"),
    maxDPtRel = cms.double( 1.0 ),
    maxDeltaR = cms.double( 0.2 ),
    filter = cms.bool(True) 
)

zPlots = cms.PSet(
    histograms = cms.VPSet(
    cms.PSet(
    min = cms.untracked.double(0.0),
    max = cms.untracked.double(200.0),
    nbins = cms.untracked.int32(200),
    name = cms.untracked.string("zMass"),
    description = cms.untracked.string("Z mass [GeV/c^{2}]"),
    plotquantity = cms.untracked.string("mass")
    ),
    cms.PSet(
    min = cms.untracked.double(0.0),
    max = cms.untracked.double(200.0),
    nbins = cms.untracked.int32(200),
    name = cms.untracked.string("mu1Pt"),
    description = cms.untracked.string("Highest muon p_{t} [GeV/c]"),
    plotquantity = cms.untracked.string("max(daughter(0).pt,daughter(1).pt)")
    ),
    cms.PSet(
    min = cms.untracked.double(0.0),
    max = cms.untracked.double(200.0),
    nbins = cms.untracked.int32(200),
    name = cms.untracked.string("mu2Pt"),
    description = cms.untracked.string("Lowest muon p_{t} [GeV/c]"),
    plotquantity = cms.untracked.string("min(daughter(0).pt,daughter(1).pt)")
    )
    )
)




goodZToMuMuPlots = cms.EDFilter(
    "CandViewHistoAnalyzer",
    zPlots,
    src = cms.InputTag("goodZToMuMuAtLeast1HLT"),
    filter = cms.bool(False)
)



ewkZMuMuGoldenSequence = cms.Sequence(
    goodAODGlobalMuons *
   dimuonsHLTFilter *  
   dimuonsGlobalAOD *
   dimuonsFilter *
   goodZToMuMuAtLeast1HLT *
   goodZToMuMuPlots    
)


