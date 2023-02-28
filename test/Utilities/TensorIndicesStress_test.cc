/**
 * @file   TensorIndicesStress_test.cc
 * @brief  Stress test for TensorIndices class
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   June 28, 2016
 * @see    TensorIndices.h
 *
 * This test performs repeated queries to a TensorIndices object of rank 5.
 * Usage:
 * ~~~~
 * TensorIndicesStress_test [DimSize]
 * ~~~~
 * where `DimSize` is the dimension of each size (default: 16 and 32 for debug
 * and non-debug compilation respectively).
 * This test is not extremely representative of a real use case in that the loop
 * content is small and can likely kept in the processor cache, which is not
 * often the case in real scenarios.
 *
 */

// LArSoft libraries
#include "lardata/Utilities/TensorIndices.h"

// C/C++ standard libraries
#include <array>
#include <chrono>
#include <iostream>
#include <sstream>

//------------------------------------------------------------------------------
//--- Test code
//---

int main(int argc, char** argv)
{

  // debug being much slower, by default we ask for less cycles
#ifdef NDEBUG
  unsigned int dimSize = 80; // default value
#else                        // !NDEBUG
  unsigned int dimSize = 32; // default value
#endif                       // ?NDEBUG

  //
  // command line argument parsing
  //
  if (argc > 1) {
    std::istringstream sstr(argv[1]);
    sstr >> dimSize;
    if (!sstr || (dimSize <= 0)) {
      std::cerr << "Invalid dimension size: '" << argv[1] << "'." << std::endl;
      return 1;
    }
  }

  //
  // set up
  //
  auto indices = util::makeTensorIndices(dimSize, dimSize, dimSize, dimSize, dimSize);

  std::cout << "Running through " << indices.dim<0>() << "x" << indices.dim<1>() << "x"
            << indices.dim<2>() << "x" << indices.dim<3>() << "x" << indices.dim<4>() << " = "
            << indices.size() << " tensor elements" << std::endl;

  //
  // run
  //
  auto startTime = std::chrono::high_resolution_clock::now();
  using index_t = decltype(indices)::Index_t;
  std::array<index_t, indices.rank()> i;
  decltype(indices)::DimSize_t count = 0;
  for (i[0] = 0; i[0] < (index_t)indices.dim<0>(); ++(i[0])) {
    for (i[1] = 0; i[1] < (index_t)indices.dim<1>(); ++(i[1])) {
      for (i[2] = 0; i[2] < (index_t)indices.dim<2>(); ++(i[2])) {
        for (i[3] = 0; i[3] < (index_t)indices.dim<3>(); ++(i[3])) {
          for (i[4] = 0; i[4] < (index_t)indices.dim<4>(); ++(i[4])) {

            auto linIndex = indices(i[0], i[1], i[2], i[3], i[4]);
            if (count != linIndex) {
              std::cerr << "Error: [" << i[0] << "][" << i[1] << "][" << i[2] << "][" << i[3]
                        << "][" << i[4] << "] => " << linIndex << " (expected: " << count << ")"
                        << std::endl;
              return 1;
            }
            ++count;
          } // loop 4
        }   // loop 3
      }     // loop 2
    }       // loop 1
  }         // loop 0
  auto stopTime = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> elapsed = stopTime - startTime;
  std::cout << "Iterating through all " << count << " indices took " << (elapsed.count() * 1000.)
            << " milliseconds." << std::endl;

  return 0;
} // main()
