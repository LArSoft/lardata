art_make(LIB_LIBRARIES
           lardataobj_AnalysisBase
           lardataobj_RecoBase
           lardataalg_DetectorInfo
           larcore_Geometry_Geometry_service
           larcorealg_Geometry
           art::Framework_Core
           art::Framework_Services_Registry
           art::Utilities
           canvas::canvas
           fhiclcpp::fhiclcpp
           cetlib_except::cetlib_except
           ROOT::Core
           ROOT::Physics) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.09.01) - deprecated: use art_make_library(), art_dictonary(), and cet_build_plugin() with explicit source lists and plugin base types

install_headers()
install_fhicl()
install_source()
