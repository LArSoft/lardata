/// @file  GetManyByRegexTag.h
/// @brief Utility to retrieve art data products whose input tags match a regex pattern.

#ifndef LARDATA_LARDATA_ARTDATAHELPER_GETMANYBYREGEXTAG_H
#define LARDATA_LARDATA_ARTDATAHELPER_GETMANYBYREGEXTAG_H

#include "canvas/Utilities/InputTag.h"
#include "lardata/ArtDataHelper/RegexTagSelector.h"

#include <stdexcept>
#include <string>
#include <typeinfo>
#include <vector>
#include <set>

/// LArSoft-specific namespace
namespace lar {

  /// LArSoft utility namespace
  namespace util {

    /**
     * @brief  Retrieve all data products of type `T` whose input tags match a
     *         regex pattern.
     *
     * Works with `art::Event` only: the matching is delegated to art through
     * an `art::SelectorBase` passed to `getInputTags()`, which
     * `gallery::Event` does not support.
     *
     * Each field of @p tag (label, instance, process) is interpreted as a
     * `std::regex` pattern. An empty field matches any value. All three
     * fields must match for a product to be selected.
     *
     * @tparam T    The data-product type (e.g. `std::vector<recob::Hit>`).
     * @tparam Evt  `art::Event` (deduced).
     * @param  e    The event.
     * @param  tag  An `art::InputTag` whose fields are regex patterns.
     * @return      Valid handles for every matching product.
     * @throws std::regex_error    If a field is not a valid regex.
     *
     * In an art module, declare `consumesMany<T>()` in the constructor.
     */
    template <typename T, typename Evt>
    std::vector<typename Evt::template HandleT<T>> getManyByRegexTag(Evt const& e,
                                                                     art::InputTag const& tag)
    {
      // an empty process means "any process"
      detail::RegexTagSelector const selector{
        tag.label(), tag.instance(), tag.process().empty() ? ".*" : tag.process()};

      std::vector<typename Evt::template HandleT<T>> handles;
      std::set<std::pair<std::string, std::string>> seen;
      for (art::InputTag const& tag : e.template getInputTags<T>(selector)) {
        if (!seen.emplace(tag.label(), tag.instance()).second) continue; // older process

        auto handle = e.template getHandle<T>(tag);
        if (handle.isValid()) handles.push_back(std::move(handle));
      }

      return handles;
    }
  } // namespace util
} // namespace lar

#endif // LARDATA_LARDATA_ARTDATAHELPER_GETMANYBYREGEXTAG_H
