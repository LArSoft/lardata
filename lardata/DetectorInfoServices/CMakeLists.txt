simple_plugin(DetectorClocksServiceStandard "service"
              art_root_io_detail
              lardataalg_DetectorInfo
              lardataobj_RawData
              art::Framework_Core
              art::Framework_Principal
              art_root_io::RootDB
              art::Persistency_Provenance
              SQLite::SQLite3
              ROOT::Core
              ROOT::RIO
              ROOT::Tree)

simple_plugin(DetectorPropertiesServiceStandard "service"
              lardataalg_DetectorInfo
              larcore_Geometry_Geometry_service
              larcorealg_Geometry
              art_root_io::RootDB
              SQLite::SQLite3
              messagefacility::MF_MessageLogger
              ROOT::Core
              ROOT::RIO)

simple_plugin(LArPropertiesServiceStandard "service"
              lardataalg_DetectorInfo
              art::Framework_Principal)

install_headers()
install_fhicl()
install_source()
