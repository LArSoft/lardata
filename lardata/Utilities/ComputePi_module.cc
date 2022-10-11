/**
 * @file   ComputePi_module.cc
 * @brief  Computes pi
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   August 19th, 2014
 */

// C/C++ standard libraries
#include <iomanip> // std::setprecision
#include <ios>     // std::fixed
#include <random>  // std::default_random_engine, std::uniform_real_distribution

// art libraries
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

namespace lar {

  /**
   * @brief Computes pi (but it does not make it available)
   *
   * This module performs a extensive computation whose duration can be
   * indirectly controlled by a parameter.
   * The time taken is supposed to be independent from the framework.
   * This is meant to help establish an absolute time scale.
   *
   * The module performs some Monte Carlo integration to compute pi.
   * The same number of cycles is used regardless the result.
   * We use a simple pseudo-random generator
   * (<tt>std::linear_congruential_engine</tt>) with a constant extraction time
   * (and poor randomness quality, and a period so small that in about 20 events
   * the sequence might repeat itself). The fluctuations of the result don't
   * reflect a fluctuation in time.
   *
   * A test performed on uboonegpvm06,fnal.gov on August 19th, 2014 on 1000
   * events with Ksamples=50000 (i.e., 50M samples per event), default seed
   * and verbosity on took 0.9179 +/- 0.0009 s, with an RMS of ~3%.
   * It was observed that processing time asymptotically decreased.
   *
   * <b>Parameters</b>
   * - <b>Ksamples</b> (integer, default: 10000) number of digits to be computed
   * - <b>Seed</b> (unsigned integer, default: 314159) chooses the seed for the
   *   Monte Carlo integration
   * - <b>Fixed</b> (boolean, default: false) if true, the same pseudo-random
   *   number sequence will be used for all events; otherwise, each event will
   *   get its own specific sequence
   * - <b>Verbose</b> (boolean, default: false) writes the result into the log
   */
  class ComputePi : public art::EDAnalyzer {
  public:
    using Counter_t = unsigned long long; ///< type used for integral counters
    using Seed_t = std::default_random_engine::result_type;
    ///< type for seed and random numbers

    explicit ComputePi(fhicl::ParameterSet const& p);

    virtual ~ComputePi() = default;

    virtual void analyze(const art::Event&) override;

    /// Returns the current best estimation of pi
    double best_pi() const { return tries ? 4. * double(hits) / double(tries) : 3.0; }

    /// Returns the current best estimation of pi
    Counter_t best_pi_tries() const { return tries; }

    static const char* VersionString; ///< version of the algorithm

  private:
    Counter_t samples; ///< number of samples to try on each event
    Seed_t seed;       ///< number of digits to compute
    bool bFixed;       ///< whether the random sequence is always the same
    bool bVerbose;     ///< whether to put stuff on screen

    std::default_random_engine generator; ///< random generator
    Counter_t hits = 0;                   ///< total number of hits
    Counter_t tries = 0;                  ///< total number of tries (samples)

  }; // class ComputePi

} // namespace lar

//------------------------------------------------------------------------------

DEFINE_ART_MODULE(lar::ComputePi)

const char* lar::ComputePi::VersionString = "1.0";

template <typename T>
inline constexpr T sqr(T v)
{
  return v * v;
}

lar::ComputePi::ComputePi(const fhicl::ParameterSet& p)
  : EDAnalyzer(p)
  , samples(p.get<Counter_t>("Ksamples", 10000) * 1000)
  , seed(p.get<Seed_t>("Seed", 314159))
  , bFixed(p.get<bool>("Fixed", false))
  , bVerbose(p.get<bool>("Verbose", false))
  , generator(seed)
{
  mf::LogInfo("ComputePi") << "version " << VersionString << " using " << samples
                           << " samples per event, random seed " << seed;
} // lar::ComputePi::ComputePi()

void lar::ComputePi::analyze(const art::Event&)
{

  // prepare our personal pseudo-random engine;
  // we'll use always the same sequence!
  std::uniform_real_distribution<float> flat(0.0, 1.0);

  // if we want to fix the random sequence, we reseed the generator
  // with the same value over and over again
  if (bFixed) generator.seed(seed);

  Counter_t local_hits = 0, tries_left = samples;
  while (tries_left-- > 0) {
    float x = flat(generator), y = flat(generator);
    if (sqr(x) + sqr(y) < 1.0) ++local_hits;
  } // while

  double local_pi = double(local_hits) / double(samples) * 4.0;
  hits += local_hits;
  tries += samples;

  if (bVerbose) {
    mf::LogInfo("ComputePi") << "today's pi = " << std::fixed << std::setprecision(9) << local_pi
                             << " (pi = " << std::fixed << std::setprecision(12) << best_pi()
                             << " after " << best_pi_tries() << " samples)";
  } // if verbose

} // lar::ComputePi::analyze()
