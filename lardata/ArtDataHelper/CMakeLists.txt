cet_make_library(
          SOURCE HitCreator.cxx
                 HitUtils.cxx
                 MVAReader.cxx
                 MVAWrapperBase.cxx
                 MVAWriter.cxx
                 TrackUtils.cxx
                 WireCreator.cxx
                 ChargedSpacePointCreator.cpp
          LIBRARIES PUBLIC
                       lardataobj::RawData
                       lardataobj::RecoBase
                       lardataobj::AnalysisBase
                       larcorealg::Geometry
                       larcore::Geometry_Geometry_service
                       art::Framework_Core
                       art::Framework_Principal
                       art::Framework_Services_Registry
                       art::Persistency_Common
                       art::Persistency_Provenance
                       art::Utilities
                    PRIVATE
                       canvas::canvas
                       cetlib_except::cetlib_except
                       ROOT::Core
                       ROOT::GenVector)

add_subdirectory(Dumpers)

install_headers()
install_fhicl()
install_source()
