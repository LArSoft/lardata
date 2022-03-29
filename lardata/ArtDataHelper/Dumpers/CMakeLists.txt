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


art_make(NO_PLUGINS LIB_LIBRARIES lardataobj_RecoBase) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.09.01) - deprecated: use art_make_library(), art_dictonary(), and cet_build_plugin() with explicit source lists and plugin base types

foreach(Dumper IN LISTS RawDataDumpers)
  simple_plugin(${Dumper} "module"
      lardataobj_RawData
      art::Framework_Services_Registry
      messagefacility::MF_MessageLogger) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.09.01) - deprecated: use cet_build_plugin() with explicit source lists and plugin base types
endforeach()


foreach(Dumper IN LISTS RecoBaseDumpers)
  simple_plugin(${Dumper} "module"
      lardataobj_RecoBase
      lardata_ArtDataHelper_Dumpers
      lardata_RecoBaseProxy
      messagefacility::MF_MessageLogger) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.09.01) - deprecated: use cet_build_plugin() with explicit source lists and plugin base types
endforeach()


foreach(Dumper IN LISTS SimulationDumpers)
  simple_plugin(${Dumper} "module"
      lardataalg_MCDumpers
      lardataobj_Simulation
      messagefacility::MF_MessageLogger
      ROOT::GenVector
      ROOT::Core) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.09.01) - deprecated: use cet_build_plugin() with explicit source lists and plugin base types
endforeach(Dumper IN LISTS SimulationDumpers)


foreach(Dumper IN LISTS MCDumpers)
  simple_plugin(${Dumper} "module" lardataobj_MCBase) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.09.01) - deprecated: use cet_build_plugin() with explicit source lists and plugin base types
endforeach()



install_headers()
install_fhicl()
install_source()
