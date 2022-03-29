cet_build_plugin(DetectorClocksServiceStandard art::service
  LIBRARIES PUBLIC
              lardataalg::DetectorInfo
              lardataobj::RawData
              art::Framework_Core
              art::Framework_Principal
              art_root_io::RootDB
              art::Persistency_Provenance
            PRIVATE
              art_root_io::detail
              SQLite::SQLite3
              ROOT::Core
              ROOT::RIO
              ROOT::Tree)

cet_build_plugin(DetectorPropertiesServiceStandard art::service
  LIBRARIES PUBLIC
              lardataalg::DetectorInfo
              larcore::Geometry_Geometry_service
              larcorealg::Geometry
              art_root_io::RootDB
            PRIVATE
              SQLite::SQLite3
              messagefacility::MF_MessageLogger
              ROOT::Core
              ROOT::RIO)

cet_build_plugin(LArPropertiesServiceStandard art::service
  LIBRARIES PUBLIC
              lardataalg::DetectorInfo
              art::Framework_Principal)

install_headers()
install_fhicl()
install_source()
