import FWCore.ParameterSet.Config as cms
import copy

process = cms.Process("ZToMuMuSubSkim")

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10000)
)

process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(
    "file:/scratch1/cms/data/summer08/skim/dimuons_skim_zmumu.root"
    )
)

process.include("FWCore/MessageLogger/data/MessageLogger.cfi")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'IDEAL_V9::All'
process.load("Configuration.StandardSequences.MagneticField_cff")

process.prunedGenParticles = cms.EDProducer(
    "GenParticlePruner",
    src = cms.InputTag("genParticles"),
    select = cms.vstring(
    "drop  *  ", # this is the default
    "keep+ pdgId = {Z0}",
    "keep+ pdgId = {W-}",
    "keep++ pdgId = {mu-}"
    )
)

process.load("ElectroWeakAnalysis.ZReco.dimuonsSequences_cff")
process.load("ElectroWeakAnalysis.ZReco.dimuonsHLTFilter_cfi")
process.load("ElectroWeakAnalysis.ZReco.dimuonsFilter_cfi")
process.load("ElectroWeakAnalysis.ZReco.dimuonsOneTrackFilter_cfi")

process.muonMatch.matched = cms.InputTag("prunedGenParticles")
process.trackMuMatch.matched = cms.InputTag("prunedGenParticles")

zSelection = cms.PSet(
    cut = cms.string("charge = 0 & daughter(0).pt > 20 & daughter(1).pt > 20 & "
                     "abs(daughter(0).eta)<2 & abs(daughter(1).eta)<2 & mass > 20"),
    isoCut = cms.double(100.0),
    isolationType = cms.string("track"),
)

process.goodZToMuMu = cms.EDFilter(
    "ZToMuMuIsolatedSelector",
    zSelection,
    src = cms.InputTag("dimuonsGlobal"),
    filter = cms.bool(True) 
)

process.nonIsolatedZToMuMu = cms.EDFilter(
    "ZToMuMuNonIsolatedSelector",
    zSelection,
    src = cms.InputTag("dimuonsGlobal"),
    filter = cms.bool(True) 
)

process.zToMuGlobalMuOneTrack = cms.EDFilter(
    "CandViewRefSelector",
    cut = cms.string("daughter(0).isGlobalMuon = 1"),
    src = cms.InputTag("dimuonsOneTrack"),
    filter = cms.bool(True)
)

process.zToMuMuOneTrack = cms.EDFilter(
    "ZToMuMuIsolatedSelector",
    zSelection,
    src = cms.InputTag("zToMuGlobalMuOneTrack"),
    filter = cms.bool(True)
)

process.zToMuMuOneStandAloneMuon = cms.EDFilter(
    "ZToMuMuIsolatedSelector",
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

process.eventInfo = cms.OutputModule (
    "AsciiOutputModule"
)

process.out = cms.OutputModule(
    "PoolOutputModule",
    fileName = cms.untracked.string('file:./Zmm.root'),
    outputCommands = cms.untracked.vstring(
      "drop *",
      "keep *_prunedGenParticles_*_ZToMuMuSubSkim",
      "keep *_selectedLayer1Muons_*_ZToMuMuSubSkim",
      "keep *_selectedLayer1TrackCands_*_ZToMuMuSubSkim",
      "keep *_dimuons_*_ZToMuMuSubSkim",
      "keep *_dimuonsOneTrack_*_ZToMuMuSubSkim",
      "keep *_dimuonsGlobal_*_ZToMuMuSubSkim",
      "keep *_dimuonsOneStandAloneMuon_*_ZToMuMuSubSkim",
      "keep *_muonMatch_*_ZToMuMuSubSkim",
      "keep *_allDimuonsMCMatch_*_ZToMuMuSubSkim",
      "keep *_goodZToMuMu_*_ZToMuMuSubSkim",
      "keep *_nonIsolatedZToMuMu_*_ZToMuMuSubSkim",
      "keep *_goodZToMuMuOneStandAloneMuon_*_ZToMuMuSubSkim",
      "keep *_goodZToMuMuOneTrack_*_ZToMuMuSubSkim"
    ),
    SelectEvents = cms.untracked.PSet(
      SelectEvents = cms.vstring(
        "goodZToMuMuPath", "nonIsolatedZToMuMuPath",
        "zToMuMuOneStandAloneMuonPath", "goodZToMuMuOneTrackPath"
      )
    )
)

process.genParticlesPath = cms.Path(
    process.prunedGenParticles
    )

process.dimuonsPath = cms.Path(
    process.dimuonsHLTFilter +
    process.goodMuonRecoForDimuon +
    process.dimuons +
    process.dimuonsGlobal +
    process.dimuonsOneStandAloneMuon +
    process.dimuonsFilter
    )

process.dimuonsOneTrackPath = cms.Path(
    process.dimuonsHLTFilter +
    process.goodMuonRecoForDimuon +
    process.dimuonsOneTrack +
    process.dimuonsOneTrackFilter
    )

process.dimuonsMCTruth = cms.Path(
    process.dimuonsHLTFilter +
    process.mcTruthForDimuons
    )

process.goodZToMuMuPath = cms.Path(
    process.goodZToMuMu 
)
 
process.nonIsolatedZToMuMuPath = cms.Path(
    process.nonIsolatedZToMuMu
)

process.zToMuMuOneStandAloneMuonPath = cms.Path(
    ~process.goodZToMuMu + 
    process.zToMuMuOneStandAloneMuon + 
    process.goodZToMuMuOneStandAloneMuon 
)

process.goodZToMuMuOneTrackPath = cms.Path(
    ~process.goodZToMuMu +
    ~process.zToMuMuOneStandAloneMuon +
    process.zToMuGlobalMuOneTrack +
    process.zToMuMuOneTrack +
    process.goodZToMuMuOneTrack 
)

process.endPath = cms.EndPath( 
    process.eventInfo +
    process.out
)

