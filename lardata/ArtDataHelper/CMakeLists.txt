art_make(LIB_LIBRARIES lardataobj_RawData
                       lardataobj_RecoBase
                       lardataobj_AnalysisBase
                       larcorealg_Geometry
                       larcore_Geometry_Geometry_service
                       art::Framework_Core
                       art::Framework_Principal
                       art::Framework_Services_Registry
                       art::Persistency_Common
                       art::Persistency_Provenance
                       art::Utilities
                       canvas::canvas
                       cetlib_except::cetlib_except
                       ROOT::Core
                       ROOT::GenVector) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.09.01) - deprecated: use art_make_library(), art_dictonary(), and cet_build_plugin() with explicit source lists and plugin base types

add_subdirectory(Dumpers)

install_headers()
install_fhicl()
install_source()
