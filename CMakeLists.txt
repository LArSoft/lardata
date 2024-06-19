# ======================================================================
#  lardata main build file
#
#  cd .../path/to/build/directory
#  source .../path/to/lardata/ups/setup_for_development <-d|-p>
#  buildtool [-I /install/path] -bt[i] -jN
# ======================================================================

cmake_minimum_required(VERSION 3.20 FATAL_ERROR)

find_package(cetmodules 3.12.00 REQUIRED)
project(lardata LANGUAGES CXX)

# cetbuildtools contains our cmake modules

include(CetCMakeEnv)
set(${PROJECT_NAME}_CMAKE_PROJECT_VERSION_STRING 10.00.00-rc4)
cet_cmake_env()

cet_set_compiler_flags(DIAGS CAUTIOUS
  WERROR
  NO_UNDEFINED
  EXTRA_FLAGS -pedantic
)

cet_report_compiler_flags(REPORT_THRESHOLD VERBOSE)

find_package(nusimdata REQUIRED EXPORT)

find_package(art_root_io REQUIRED EXPORT)
find_package(canvas_root_io REQUIRED)
find_package(art REQUIRED EXPORT)
find_package(canvas REQUIRED EXPORT)
find_package(fhiclcpp REQUIRED EXPORT)
find_package(messagefacility REQUIRED EXPORT)
find_package(cetlib_except REQUIRED EXPORT)
find_package(cetlib REQUIRED EXPORT)

find_package(Range-v3 REQUIRED EXPORT)
find_package(FFTW3 REQUIRED EXPORT)
find_package(Boost COMPONENTS date_time serialization REQUIRED EXPORT)
find_package(ROOT COMPONENTS Core FFTW GenVector Hist MathCore Physics RIO Tree REQUIRED EXPORT)
find_package(PostgreSQL REQUIRED EXPORT)

find_package(larcoreobj REQUIRED EXPORT)
find_package(larcorealg REQUIRED EXPORT)
find_package(lardataobj REQUIRED EXPORT)
find_package(lardataalg REQUIRED EXPORT)
find_package(larcore REQUIRED EXPORT)
find_package(larvecutils REQUIRED EXPORT)

include(ArtDictionary)
include(ArtMake)
include(BasicPlugin)

cet_cmake_module_directories(Modules BINARY)

# source
add_subdirectory(lardata)

# tests
add_subdirectory(test)

# ups - table and config files

# packaging utility
cet_cmake_config()
