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
              ROOT::Tree) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.09.01) - deprecated: use cet_build_plugin() with explicit source lists and plugin base types

simple_plugin(DetectorPropertiesServiceStandard "service"
              lardataalg_DetectorInfo
              larcore_Geometry_Geometry_service
              larcorealg_Geometry
              art_root_io::RootDB
              SQLite::SQLite3
              messagefacility::MF_MessageLogger
              ROOT::Core
              ROOT::RIO) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.09.01) - deprecated: use cet_build_plugin() with explicit source lists and plugin base types

simple_plugin(LArPropertiesServiceStandard "service"
              lardataalg_DetectorInfo
              art::Framework_Principal) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.09.01) - deprecated: use cet_build_plugin() with explicit source lists and plugin base types

install_headers()
install_fhicl()
install_source()
