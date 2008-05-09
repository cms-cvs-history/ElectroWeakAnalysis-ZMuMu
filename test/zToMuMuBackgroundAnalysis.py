import FWCore.ParameterSet.Config as cms

process = cms.Process("ZToMuMuBackgroundAnalysis")

process.include("FWCore/MessageLogger/data/MessageLogger.cfi")

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring("file:dimuons.root")
)

process.TFileService = cms.Service(
    "TFileService",
    fileName = cms.string("zToMuMuBackgroundAnalysis.root")
)

zSelection = cms.PSet(
    cut = cms.string("daughter(0).pt > 20 & daughter(1).pt > 20 & abs(daughter(0).eta)<2 & abs(daughter(1).eta)<2 & mass > 20")	
)

process.goodZToMuMu = cms.EDFilter(
    "CandViewRefSelector",
    zSelection,
    src = cms.InputTag("dimuonsGlobal"),
    filter = cms.bool(True),
)

process.zToMuMuOneTrack = cms.EDFilter(
    "CandViewRefSelector",
    zSelection,
    src = cms.InputTag("dimuonsOneTrack"),
    filter = cms.bool(True)
)

process.zToMuMuOneStandAloneMuon = cms.EDFilter(
    "CandViewRefSelector",
    zSelection,
    src = cms.InputTag("dimuonsOneStandAloneMuon"),
    filter = cms.bool(True)
)

process.goodZToMuMuOneTrack = cms.EDFilter(
    "ZMuMuOverlapExclusionSelector",
    src = cms.InputTag("zToMuMuOneTrack"),
    overlap = cms.InputTag("goodZToMuMu"),
    filter = cms.bool(True)
)

process.goodZToMuMuOneStandAloneMuon = cms.EDFilter(
    "ZMuMuOverlapExclusionSelector",
    src = cms.InputTag("zToMuMuOneStandAloneMuon"),
    overlap = cms.InputTag("goodZToMuMu"),
    filter = cms.bool(True)
)

bin0CutTwoDaughters = cms.PSet(
    cut = cms.string("0 < daughter(0).eta < 0.5 | 0 < daughter(1).eta < 0.5")
)

bin1CutTwoDaughters = cms.PSet (
    cut = cms.string("0.5 < daughter(0).eta < 1 | 0.5 < daughter(1).eta < 1")
)

bin0CutOneDaughter = cms.PSet (
    cut = cms.string("0 < daughter(1).eta < 0.5")
)

bin1CutOneDaughter = cms.PSet(
    cut = cms.string("0.5 < daughter(1).eta < 1")
)

process.goodZToMuMuBin0 = cms.EDFilter(
    "CandViewRefSelector",
    bin0CutTwoDaughters,
    src = cms.InputTag("goodZToMuMu"),
    filter = cms.bool(True)
)

process.goodZToMuMuBin1 = cms.EDFilter(
    "CandViewRefSelector",
    bin1CutTwoDaughters,
    src = cms.InputTag("goodZToMuMu"),
    filter = cms.bool(True)
)

process.goodZToMuMuOneStandAloneMuonBin0 = cms.EDFilter(
    "CandViewRefSelector",
    bin0CutOneDaughter,
    src = cms.InputTag("goodZToMuMuOneStandAloneMuon"),
    filter = cms.bool(True)
)

process.goodZToMuMuOneStandAloneMuonBin1 = cms.EDFilter(
    "CandViewRefSelector",
    bin1CutOneDaughter,
    src = cms.InputTag("goodZToMuMuOneStandAloneMuon"),
    filter = cms.bool(True)
)

process.goodZToMuMuOneTrackBin0 = cms.EDFilter(
    "CandViewRefSelector",
    bin0CutOneDaughter,
    src = cms.InputTag("goodZToMuMuOneTrack"),
    filter = cms.bool(True)
)

process.goodZToMuMuOneTrackBin1 = cms.EDFilter(
    "CandViewRefSelector",
    bin1CutOneDaughter,
    src = cms.InputTag("goodZToMuMuOneTrack"),
    filter = cms.bool(True)
)

zPlots = cms.PSet(
    histograms = cms.VPSet(
    cms.PSet(
    min = cms.untracked.double(0.0),
    max = cms.untracked.double(200.0),
    nbins = cms.untracked.int32(200),
    name = cms.untracked.string("zMass"),
    description = cms.untracked.string("Z mass [GeV/c^(2)]"),
    plotquantity = cms.untracked.string("mass")
    ),
    cms.PSet(
    min = cms.untracked.double(0.0),
    max = cms.untracked.double(200.0),
    nbins = cms.untracked.int32(200),
    name = cms.untracked.string("mu1Pt"),
    description = cms.untracked.string("Highest muon p_(t) [GeV/c]"),
    plotquantity = cms.untracked.string("max(daughter(0).pt,daughter(1).pt)")
    ),
    cms.PSet(
    min = cms.untracked.double(0.0),
    max = cms.untracked.double(200.0),
    nbins = cms.untracked.int32(200),
    name = cms.untracked.string("mu2Pt"),
    description = cms.untracked.string("Lowest muon p_(t) [GeV/c]"),
    plotquantity = cms.untracked.string("min(daughter(0).pt,daughter(1).pt)")
    )
    )
)

process.goodZToMuMuPlots = cms.EDAnalyzer(
    "CandViewHistoAnalyzer",
    zPlots,
    src = cms.InputTag("goodZToMuMu")
)

process.goodZToMuMuOneTrackPlots = cms.EDAnalyzer(
    "CandViewHistoAnalyzer",
    zPlots,
    src = cms.InputTag("goodZToMuMuOneTrack")
)

process.goodZToMuMuOneStandAloneMuonPlots = cms.EDAnalyzer(
    "CandViewHistoAnalyzer",
    zPlots,
    src = cms.InputTag("goodZToMuMuOneStandAloneMuon")
)

process.eventInfo = cms.OutputModule (
    "AsciiOutputModule"
)

process.goodZToMuMuPlotsBin0 = cms.EDAnalyzer(
    "CandViewHistoAnalyzer",
    zPlots,
    src = cms.InputTag("goodZToMuMuBin0")
)

process.goodZToMuMuPlotsBin1 = cms.EDAnalyzer(
    "CandViewHistoAnalyzer",
    zPlots,
    src = cms.InputTag("goodZToMuMuBin1")
)

process.goodZToMuMuOneStandAloneMuonPlotsBin0 = cms.EDAnalyzer(
    "CandViewHistoAnalyzer",
    zPlots,
    src = cms.InputTag("goodZToMuMuOneStandAloneMuonBin0")
)

process.goodZToMuMuOneStandAloneMuonPlotsBin1 = cms.EDAnalyzer(
    "CandViewHistoAnalyzer",
    zPlots,
    src = cms.InputTag("goodZToMuMuOneStandAloneMuonBin1")
)

process.goodZToMuMuOneTrackPlotsBin0 = cms.EDAnalyzer(
    "CandViewHistoAnalyzer",
    zPlots,
    src = cms.InputTag("goodZToMuMuOneTrackBin0")
)

process.goodZToMuMuOneTrackPlotsBin1 = cms.EDAnalyzer(
    "CandViewHistoAnalyzer",
    zPlots,
    src = cms.InputTag("goodZToMuMuOneTrackBin1")
)

process.zToMuMuPath = cms.Path(  
    process.goodZToMuMu + 
    process.goodZToMuMuPlots +
    process.goodZToMuMuBin0 + process.goodZToMuMuPlotsBin0 + 
    process.goodZToMuMuBin1 + process.goodZToMuMuPlotsBin1 
)

process.zToMuMuOneStandAloneMuonPath = cms.Path(
    ~process.goodZToMuMu + 
    process.zToMuMuOneStandAloneMuon + 
    process.goodZToMuMuOneStandAloneMuon +
    process.goodZToMuMuOneStandAloneMuonPlots +
    process.goodZToMuMuOneStandAloneMuonBin0 + process.goodZToMuMuOneStandAloneMuonPlotsBin0 +
    process.goodZToMuMuOneStandAloneMuonBin1 + process.goodZToMuMuOneStandAloneMuonPlotsBin1 
)

process.zToMuMuOneTrackPath = cms.Path(
    ~process.goodZToMuMu +
    ~process.zToMuMuOneStandAloneMuon +
    process.zToMuMuOneTrack +
    process.goodZToMuMuOneTrack +
    process.goodZToMuMuOneTrackPlots +
    process.goodZToMuMuOneTrackBin0 + process.goodZToMuMuOneTrackPlotsBin0 + 
    process.goodZToMuMuOneTrackBin1 + process.goodZToMuMuOneTrackPlotsBin1  
)
  
process.endPath = cms.EndPath( 
    process.eventInfo 
)

