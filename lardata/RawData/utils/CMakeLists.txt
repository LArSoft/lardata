art_make(NO_PLUGINS
         EXCLUDE LArRawInputDriverUBooNE.cxx
         LIB_LIBRARIES lardataobj_RawData
                       larcorealg_Geometry
                       larcoreobj_SummaryData
                       art::Framework_Core
                       art::Framework_IO_Sources
                       art::Framework_Principal
                       art::Persistency_Common
                       art::Persistency_Provenance
                       canvas::canvas
                       cetlib_except::cetlib_except
                       PQ
                       Boost::serialization
                       Boost::date_time
                       ROOT::Core
                       ROOT::RIO
                       ROOT::Tree
        )

simple_plugin(LArRawInputSource "source"
  lardataobj_RawData
  lardata_RawData_utils
  art::Framework_IO_Sources
)

install_headers()
install_fhicl()
install_source()
