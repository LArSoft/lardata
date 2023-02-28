/**
 * @file SeedCreator.h
 * @brief Tools and modules for checking out the basics of the Monte Carlo.
 * @author Ben Carls (bcarls@fnal.gov)
 *
 * Coding note: Never put @code #include "SeedCreator.h" @endcode in your code.
 * It would force a dependency in your class on every other class in the
 * Simulation directory; e.g., if I changed the `MCTruth` code, your
 * code that generates Electrons would re-compile.  This class exists
 * solely as a bookkeeping tool.
 */

#ifndef SEEDCREATOR_H
#define SEEDCREATOR_H

#include "TRandom3.h"

namespace SeedCreator {

  inline unsigned int CreateRandomNumberSeed()
  {

    // the maximum allowed seed for the art::RandomNumberGenerator
    // is 900000000. Use TRandom3 to get the seed value in that range.
    // Instantiating TRandom3 with a 0 means that its seed is set based
    // on the TUUID and should always be random, even for jobs running on the
    // same machine
    TRandom3 rand(0);
    return rand.Integer(900000000);
  }

}

#endif // SEEDCREATOR_H
