/// @file  GetManyByRegexTag.h
/// @brief Utility to retrieve art data products whose input tags match a regex pattern.

#ifndef LARDATA_LARDATA_ARTDATAHELPER_GETMANYBYREGEXTAG_H
#define LARDATA_LARDATA_ARTDATAHELPER_GETMANYBYREGEXTAG_H

#include "canvas/Utilities/InputTag.h"
#include "lardata/ArtDataHelper/RegexTagSelector.h"

#include <string>
#include <vector>

/// LArSoft-specific namespace
namespace lar {

  /// LArSoft utility namespace
  namespace util {

    /**
     * @brief  Retrieve the data products of type `T` whose input tags match
     *         regex patterns, taken from the most recent matching process.
     *
     * @tparam T    The data-product type (e.g. `std::vector<recob::Hit>`).
     * @tparam Evt  `art::Event` (deduced).
     * @param  e    The event.
     * @param  tag  An `art::InputTag` whose three fields are regex patterns.
     * @return      Valid handles to the selected products (possibly none).
     * @throws std::regex_error    If a field is not a valid regex.
     *
     * Works with `art::Event` only: the matching is delegated to art through
     * an `art::SelectorBase` passed to `getInputTags()`, which
     * `gallery::Event` does not support.
     * In an art module, declare `consumesMany<T>()` in the constructor.
     *
     * Selection proceeds in three steps.
     *
     * 1. **Pattern matching.** Each field of @p tag is a `std::regex` that must
     *    match the *whole* corresponding field of a product's input tag
     *    (`std::regex_match`, not `std::regex_search`: `"gaushit"` does not
     *    match `"gaushitStage2"`). All three fields must match.
     *    - label: matched as given; an empty pattern matches only an empty
     *      module label (i.e. in practice nothing);
     *    - instance: matched as given; an empty pattern matches only products
     *      with an empty instance name. Use `".*"` to accept any instance;
     *    - process: an empty pattern is replaced by `".*"` and matches any
     *      process name.
     * 2. **Newest matching process only.** art returns the matching tags
     *    grouped by process, from the most recent process to the oldest one
     *    (the current process first, then the input file's process history
     *    in reverse order). Only the tags of the *first* process in that
     *    sequence are kept, i.e. the most recent process that has at least
     *    one product matching the patterns. Matching products from older
     *    processes are ignored, *even if the newest process does not have a
     *    product with the same label and instance*.
     * 3. **Validity.** A handle is retrieved for each kept tag, and only the
     *    valid ones are returned. The process choice of step 2 is not
     *    revisited: if all the products of the newest matching process are
     *    invalid (e.g. dropped on input), the result is empty, without
     *    falling back to an older process.
     *
     * Examples, for an event with history `reco1` -> `reco2` (`reco2` newest)
     * and products of type `T`:
     *
     *     gaushit::reco1         gaushit::reco2
     *     gaushit:uncalib:reco1  fasthit::reco2
     *     gaushitStage2::reco1
     *
     * | @p tag (label:instance:process) | result                           |
     * | ------------------------------- | -------------------------------- |
     * | `gaushit::`                     | `gaushit::reco2`                 |
     * | `gaushit::reco1`                | `gaushit::reco1`                 |
     * | `gaushit.*::reco1`              | `gaushit::reco1`, `gaushitStage2::reco1` |
     * | `gaushit.*::`                   | `gaushit::reco2` (`gaushitStage2` only exists in `reco1`) |
     * | `gaushit:.*:`                   | `gaushit::reco2` (`uncalib` only exists in `reco1`) |
     * | `gaushit:.*:reco1`              | `gaushit::reco1`, `gaushit:uncalib:reco1` |
     * | `gaushit:uncalib:`              | `gaushit:uncalib:reco1` (newest process *with a match*) |
     * | `.*::`                          | `gaushit::reco2`, `fasthit::reco2` |
     * | `.*::reco.*`                    | same as `.*::`                   |
     * | `::` (all empty)                | nothing (no product has an empty label) |
     */
    template <typename T, typename Evt>
    std::vector<typename Evt::template HandleT<T>> getManyByRegexTag(Evt const& e,
                                                                     art::InputTag const& tag)
    {
      // an empty process means "any process"
      detail::RegexTagSelector const selector{
        tag.label(), tag.instance(), tag.process().empty() ? ".*" : tag.process()};

      // art returns the tags grouped by process, newest process first:
      // keep only those from the newest process with a matching product
      std::vector<typename Evt::template HandleT<T>> handles;
      std::string const* newestProcess = nullptr;
      for (art::InputTag const& t : e.template getInputTags<T>(selector)) {
        if (!newestProcess)
          newestProcess = &t.process();
        else if (t.process() != *newestProcess)
          break; // older processes from here on

        auto handle = e.template getHandle<T>(t);
        if (handle.isValid()) handles.push_back(std::move(handle));
      }

      return handles;
    }
  } // namespace util
} // namespace lar

#endif // LARDATA_LARDATA_ARTDATAHELPER_GETMANYBYREGEXTAG_H
