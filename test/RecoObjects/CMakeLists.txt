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

if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND
    CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL "13.2")
  # GCC 13.2 (and above) aggressively/spuriously warns about array-bounds issues.
  target_compile_options(LATest PRIVATE "-Wno-array-bounds;-Wno-stringop-overflow")
endif()

install_headers()
install_fhicl()
install_source()
