def customize_removeMCFromOFP2(process):
    process.vertexrecoTask.remove(process.quickTrackAssociatorByHits)
    process.vertexrecoTask.remove(process.tpClusterProducer)
    process.vertexrecoTask.remove(process.trackTimeValueMapProducer)

    process.particleFlowRecoTask.remove(process.quickTrackAssociatorByHits)
    process.particleFlowRecoTask.remove(process.simPFProducer)
    process.particleFlowRecoTask.remove(process.tpClusterProducer)


    return process
