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
           ROOT::Physics)

install_headers()
install_fhicl()
install_source()
