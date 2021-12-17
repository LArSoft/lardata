# ======================================================================
#  lardata main build file
#
#  cd .../path/to/build/directory
#  source .../path/to/lardata/ups/setup_for_development <-d|-p>
#  buildtool [-I /install/path] -bt[i] -jN
# ======================================================================

cmake_minimum_required(VERSION 3.19 FATAL_ERROR)

project(lardata VERSION 09.07.00 LANGUAGES CXX)

# cetbuildtools contains our cmake modules
find_package(cetbuildtools REQUIRED)

include(CetCMakeEnv)
cet_cmake_env()

cet_set_compiler_flags(DIAGS CAUTIOUS
  WERROR
  NO_UNDEFINED
  EXTRA_FLAGS -pedantic
)

cet_report_compiler_flags()

find_package(Range-v3 REQUIRED PUBLIC)
find_package(cetlib_except REQUIRED PUBLIC)
find_package(cetlib REQUIRED PUBLIC)
find_package(messagefacility REQUIRED PUBLIC)
find_package(canvas REQUIRED PUBLIC)
find_package(canvas_root_io REQUIRED)
find_package(art REQUIRED PUBLIC)
find_package(art_root_io REQUIRED PUBLIC)
find_package(nusimdata REQUIRED PUBLIC)
find_package(FFTW3 REQUIRED PUBLIC)
find_package(Boost COMPONENTS serialization date_time REQUIRED PUBLIC)
find_package(ROOT COMPONENTS Core GenVector RIO Tree FFTW Hist MathCore Physics REQUIRED PUBLIC)
find_package(PostgreSQL REQUIRED PUBLIC)

find_ups_product( larcoreobj )
find_ups_product( larcorealg )
find_ups_product( lardataobj )
find_ups_product( lardataalg )
find_ups_product( larcore  )
find_ups_product( larvecutils  )


# macros for artdaq_dictionary and simple_plugin
include(ArtDictionary)
include(ArtMake)
include(BuildPlugins)

# source
add_subdirectory(lardata)

# tests
add_subdirectory(test)

# ups - table and config files
add_subdirectory(ups)

# packaging utility
include(UseCPack)
