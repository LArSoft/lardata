/**
 * @file   HitDataProductChecker_module.cc
 * @brief  Module verifying the presence of data products.
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   June 12, 2017
 */

// LArSoft libraries
#include "lardataobj/RecoBase/Hit.h"

// framework libraries
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Comment.h"
#include "fhiclcpp/types/Name.h"
#include "fhiclcpp/types/OptionalAtom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Table.h"

// C/C++ standard libraries
#include <iterator> // std::back_inserter()
#include <string>

namespace recob {
  namespace test {

    /**
    * @brief Module verifying the presence of data products.
    *
    * Throws an exception on failure.
    *
    * Service requirements
    * =====================
    *
    * This module requires no service.
    *
    * Configuration parameters
    * =========================
    *
    * * *hits* (list, _mandatory_) each entry defines a check on a single
    *     hit collection data product (`std::vector<recob::Hit>`). Each entry
    *     is a table containing:
    *     * *name* (string, _mandatory_): input tag of the data product
    *     * *exists* (boolean, default: _true_): if `true`, the data product is
    *         expected to exist, if false it is expected not to exist
    *     * *expected* (non-negative integral number): if specified, data
    *         collection size is checked to match this number
    *
    */
    class HitDataProductChecker : public art::EDAnalyzer {

    public:
      struct Config {

        using Name = fhicl::Name;
        using Comment = fhicl::Comment;

        struct TargetInfo {

          fhicl::Atom<art::InputTag> name{Name("name"),
                                          Comment("Input tag of the data product to be checked")};

          fhicl::Atom<bool> exists{Name("exists"),
                                   Comment("whether the data product must exist or must not exist"),
                                   true};

          fhicl::OptionalAtom<unsigned int> expected{
            Name("expected"),
            Comment("Number of expected entries (not checked if not specified).")};

        }; // TargetInfo subclass

        fhicl::Sequence<fhicl::Table<TargetInfo>> hits{
          Name("hits"),
          Comment("list of hit collections and number of expected entries")};

      }; // Config

      using Parameters = art::EDAnalyzer::Table<Config>;

      explicit HitDataProductChecker(Parameters const& config);

      virtual void analyze(const art::Event& event) override;

    private:
      /// Configuration for a single data product check.
      struct TargetInfo_t {

        art::InputTag name; ///< Data product name.

        unsigned int expectedEntries; ///< Number of expected entries.

        /// Whether data product must exist or must not exist.
        bool bExists = true;

        bool bCheckEntries; ///< Whether to check the number of entries.

        TargetInfo_t() = default;

        TargetInfo_t(Config::TargetInfo const& config)
          : name(config.name()), bExists(config.exists())
        {
          bCheckEntries = config.expected(expectedEntries);
        }

      }; // TargetInfo_t

      /// Configuration of all checks on hit collections.
      std::vector<TargetInfo_t> fHitTargets;

      /**
       * @brief Checks the specified data product.
       * @tparam DATA type of data product to be checked
       * @param event the event to read the data product from
       * @param targetInfo details of the data product expected information
       * @param desc description of the content of the data product
       * @throw art::Exception (`art::errors::ProductNotFound`) when unexpected
       *
       * The checks include:
       * * existence of the data product
       * * size of the data product collection (optional)
       *
       */
      template <typename DATA>
      void checkDataProducts(art::Event const& event,
                             TargetInfo_t const& targetInfo,
                             std::string desc);

    }; // HitDataProductChecker

    DEFINE_ART_MODULE(HitDataProductChecker)

  } // namespace test
} // namespace recob

//------------------------------------------------------------------------------
//--- implementation
//---
//----------------------------------------------------------------------------
recob::test::HitDataProductChecker::HitDataProductChecker(Parameters const& config)
  : EDAnalyzer(config)
{
  decltype(auto) hits = config().hits();
  fHitTargets.reserve(hits.size());
  std::copy(hits.begin(), hits.end(), std::back_inserter(fHitTargets));
} // HitDataProductChecker::HitDataProductChecker()

//----------------------------------------------------------------------------
void recob::test::HitDataProductChecker::analyze(art::Event const& event)
{

  for (auto const& targetInfo : fHitTargets)
    checkDataProducts<std::vector<recob::Hit>>(event, targetInfo, "hits");

} // HitCollectionCreatorTest::analyze()

//----------------------------------------------------------------------------
template <typename DATA>
void recob::test::HitDataProductChecker::checkDataProducts(art::Event const& event,
                                                           TargetInfo_t const& targetInfo,
                                                           std::string desc)
{
  using Product_t = DATA;

  art::InputTag tag = targetInfo.name;

  art::Handle<Product_t> data;
  bool found = event.getByLabel<Product_t>(tag, data);
  if (found && !targetInfo.bExists) {
    throw art::Exception(art::errors::ProductNotFound)
      << "Data product '" << tag << "' (" << desc
      << ") was expected not to exist, and there it is instead! (with " << data->size()
      << " elements)";
  }
  if (!found && targetInfo.bExists) {
    throw art::Exception(art::errors::ProductNotFound)
      << "Data product '" << tag << "' (" << desc << ") was expected, but there is none.";
  }

  if (targetInfo.bCheckEntries && (data->size() != targetInfo.expectedEntries)) {
    throw art::Exception(art::errors::ProductNotFound)
      << "Data product '" << tag << "' (" << desc << ") was expected to have "
      << targetInfo.expectedEntries << " entires, but it has " << data->size() << "!";
  }

} // recob::test::HitDataProductChecker::HitCollectionCreator_test()

//----------------------------------------------------------------------------
