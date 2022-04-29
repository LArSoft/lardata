cet_test(SurfYZTest USE_BOOST_UNIT
  LIBRARIES PRIVATE
  lardata_RecoObjects
)
cet_test(SurfXYZTest USE_BOOST_UNIT
  LIBRARIES PRIVATE
  lardata_RecoObjects
)
cet_test(SurfYZLineTest USE_BOOST_UNIT
  LIBRARIES PRIVATE
  lardata_RecoObjects
)
cet_test(TrackTest
  LIBRARIES PRIVATE
  lardata_RecoObjects
)
cet_test(LATest
  LIBRARIES PRIVATE
  lardata_RecoObjects
)

install_headers()
install_fhicl()
install_source()
