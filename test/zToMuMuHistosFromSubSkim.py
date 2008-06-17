import FWCore.ParameterSet.Config as cms
import copy

process = cms.Process("ZToMuMuHistosFromSubSkim")

process.include("FWCore/MessageLogger/data/MessageLogger.cfi")

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(
    "file:/data1/cmsdata/dimuons/s156-subSkim/Zmm.root"
    )
)

process.TFileService = cms.Service(
    "TFileService",
    fileName = cms.string("dimuons_zmumu_histos.root")
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

goodZToMuMuPlotsTemplate = cms.EDAnalyzer(
    "CandViewHistoAnalyzer",
    zPlots,
    src = cms.InputTag("goodZToMuMu")
)

process.goodZToMuMuPlots = goodZToMuMuPlotsTemplate
nonIsolatedZToMuMuPlots = copy.deepcopy(goodZToMuMuPlotsTemplate)
nonIsolatedZToMuMuPlots.src = cms.InputTag("nonIsolatedZToMuMu")
setattr(process, "nonIsolatedZToMuMuPlots", nonIsolatedZToMuMuPlots)

process.eventInfo = cms.OutputModule (
    "AsciiOutputModule"
)

process.goodZToMuMuPlots = goodZToMuMuPlotsTemplate
process.goodZToMuMuOneStandAloneMuonPlots = copy.deepcopy(goodZToMuMuPlotsTemplate)
process.goodZToMuMuOneStandAloneMuonPlots.src = cms.InputTag("goodZToMuMuOneStandAloneMuon")
process.goodZToMuMuOneTrackPlots = copy.deepcopy(goodZToMuMuPlotsTemplate)
process.goodZToMuMuOneTrackPlots.src = cms.InputTag("goodZToMuMuOneTrack")

process.goodZToMuMuExists = cms.EDFilter(
    "CandCollectionExistFilter",
    src = process.nonIsolatedZToMuMuPlots.src
    )

process.nonIsolatedZToMuMuExists = cms.EDFilter(
    "CandCollectionExistFilter",
    src = process.nonIsolatedZToMuMuPlots.src
    )

process.goodZToMuMuOneStandAloneMuonExists = cms.EDFilter(
    "CandCollectionExistFilter",
    src = process.goodZToMuMuOneStandAloneMuonPlots.src
    )

process.goodZToMuMuOneTrackExists = cms.EDFilter(
    "CandCollectionExistFilter",
    src = process.goodZToMuMuOneTrackPlots.src
    )

process.goodZToMuMuPath = cms.Path(
    process.goodZToMuMuExists *
    process.goodZToMuMuPlots
)
    
process.nonIsolatedZToMuMuPath = cms.Path(
    process.nonIsolatedZToMuMuExists *
    process.nonIsolatedZToMuMuPlots
)

process.goodZToMuMuOneStandAloneMuonPath = cms.Path(
    process.goodZToMuMuOneStandAloneMuonExists *
    process.goodZToMuMuOneStandAloneMuonPlots
)
    
process.goodZToMuMuOneTrackPath = cms.Path(
    process.goodZToMuMuOneTrackExists *
    process.goodZToMuMuOneTrackPlots
)

process.endPath = cms.EndPath( 
    process.eventInfo 
)

