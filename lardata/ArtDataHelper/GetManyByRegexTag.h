/// @file  GetManyByRegexTag.h
/// @brief Utility to retrieve art data products whose input tags match a regex pattern.

#ifndef LARDATA_LARDATA_ARTDATAHELPER_GETMANYBYREGEXTAG_H
#define LARDATA_LARDATA_ARTDATAHELPER_GETMANYBYREGEXTAG_H

#include "canvas/Utilities/InputTag.h"
#include "lardata/ArtDataHelper/RegexTagMatcher.h"

#include <stdexcept>
#include <string>
#include <typeinfo>
#include <vector>

/// LArSoft-specific namespace
namespace lar {

  /// LArSoft utility namespace
  namespace util {

    /**
     * @brief  Retrieve all data products of type `T` whose input tags match a
     *         regex pattern.
     *
     * Works with both `art::Event` and `gallery::Event`.
     *
     * Each field of @p tag (label, instance, process) is interpreted as a
     * `std::regex` pattern. An empty field matches any value. All three
     * fields must match for a product to be selected.
     *
     * @tparam T    The data-product type (e.g. `std::vector<recob::Hit>`).
     * @tparam Evt  `art::Event` or `gallery::Event` (deduced).
     * @param  e    The event.
     * @param  tag  An `art::InputTag` whose fields are regex patterns.
     * @return      Valid handles for every matching product.
     * @throws std::runtime_error  If no product of type `T` matches.
     * @throws std::regex_error    If a field is not a valid regex.
     *
     * In an art module, declare `consumesMany<T>()` in the constructor.
     */
    template <typename T, typename Evt>
    std::vector<typename Evt::template HandleT<T>> getManyByRegexTag(Evt const& e,
                                                                     art::InputTag const& tag)
    {
      detail::RegexTagMatcher const matcher{tag};

      std::vector<typename Evt::template HandleT<T>> handles;
      for (art::InputTag const& t : e.template getInputTags<T>()) {

        if (!matcher.match(t)) { continue; }

        auto handle = e.template getHandle<T>(t);
        if (handle.isValid()) handles.push_back(std::move(handle));
      }

      return handles;
    }
  } // namespace util
} // namespace lar

#endif // LARDATA_LARDATA_ARTDATAHELPER_GETMANYBYREGEXTAG_HH
