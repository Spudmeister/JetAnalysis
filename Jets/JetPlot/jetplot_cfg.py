import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:/hdfs/cms/phedex/store/mc/Summer12/VBF_HToGG_M-125_8TeV-powheg-pythia6/GEN-SIM/START50_V13-v1/0000/06799BF0-449F-E111-A98B-003048F9EB46.root'
    )
)

# Output file
process.TFileService = cms.Service("TFileService",
     fileName = cms.string("test.root")
)

process.demo = cms.EDAnalyzer('JetPlot'
)


process.p = cms.Path(process.demo)
