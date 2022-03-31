#
# Since dumper modules have similar dependencies depending on the data they dump
# we collect them in groups:
#
set(RawDataDumpers
  "DumpRawDigits"
  "DumpOpDetWaveforms")
set(RecoBaseDumpers
  "DumpWires"
  "DumpHits"
  "DumpClusters"
  "DumpSeeds"
  "DumpPCAxes"
  "DumpPFParticles"
  "DumpSpacePoints"
  "DumpChargedSpacePoints"
  "DumpTracks"
  "DumpVertices")
set(SimulationDumpers
  "DumpSimPhotonsLite"
  "DumpSimEnergyDeposits"
  "DumpParticleIDs"
  )
set(MCDumper)


cet_make_library( SOURCE PCAxisDumpers.cc
                         SpacePointDumpers.cc
    LIBRARIES PUBLIC lardataobj::RecoBase)

foreach(Dumper IN LISTS RawDataDumpers)
  cet_build_plugin(${Dumper} art::module
    LIBRARIES REG
      larcorealg::CoreUtils
      lardataobj::RawData
      art::Framework_Services_Registry
      messagefacility::MF_MessageLogger)
endforeach()


foreach(Dumper IN LISTS RecoBaseDumpers)
  cet_build_plugin(${Dumper} art::module
    LIBRARIES REG
      lardataobj::RecoBase
      lardata_ArtDataHelper_Dumpers
      lardata_RecoBaseProxy
      messagefacility::MF_MessageLogger)
endforeach()


foreach(Dumper IN LISTS SimulationDumpers)
  cet_build_plugin(${Dumper} art::module
    LIBRARIES REG
      larcorealg::CoreUtils
      lardataalg::MCDumpers
      lardataobj::Simulation
      messagefacility::MF_MessageLogger
      ROOT::GenVector
      ROOT::Core)
endforeach(Dumper IN LISTS SimulationDumpers)


foreach(Dumper IN LISTS MCDumpers)
  cet_build_plugin(${Dumper} art::module
    LIBRARIES REG
                   lardataobj::MCBase)
endforeach()



install_headers()
install_fhicl()
install_source()
