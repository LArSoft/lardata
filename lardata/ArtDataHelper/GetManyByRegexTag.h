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
     *         regex patterns.
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
     * For a tag without regex special characters, the result is consistent
     * with `e.getHandle<T>(tag)`: it holds the very product that call returns,
     * or nothing if that call returns an invalid handle.
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
     *    - process: see step 2.
     * 2. **Process selection.** art returns the matching tags grouped by
     *    process, from the most recent process to the oldest one (the current
     *    process first, then the input file's process history in reverse
     *    order).
     *    - **Empty process: newest process only**, as `e.getHandle<T>()` does.
     *      Only the products of the *first* process in that sequence with at
     *      least one *valid* matching product are returned; matching products
     *      from older processes are ignored, *even if that process does not
     *      have a product with the same label and instance*.
     *    - **Explicit process: all matching processes.** The process is a
     *      regex like the other fields, and the products of *every* matching
     *      process are returned (e.g. `".*"` collects a product from all the
     *      processes that made it). A literal process name selects that
     *      process only, as `e.getHandle<T>()` does.
     * 3. **Validity.** Only valid handles are returned. A product that is
     *    listed but not available in this event (e.g. its producer did not
     *    run on it) is skipped, and with an empty process the search moves
     *    on to the next older process, as `e.getHandle<T>()` does.
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
     * | `gaushit::`                     | `gaushit::reco2` (as `e.getHandle<T>()`) |
     * | `gaushit::reco1`                | `gaushit::reco1` (as `e.getHandle<T>()`) |
     * | `gaushit::.*`                   | `gaushit::reco1`, `gaushit::reco2` |
     * | `gaushit.*::`                   | `gaushit::reco2` (`gaushitStage2` only exists in `reco1`) |
     * | `gaushit.*::.*`                 | `gaushit::reco1`, `gaushitStage2::reco1`, `gaushit::reco2` |
     * | `gaushit.*::reco1`              | `gaushit::reco1`, `gaushitStage2::reco1` |
     * | `gaushit:.*:`                   | `gaushit::reco2` (`uncalib` only exists in `reco1`) |
     * | `gaushit:.*:reco1`              | `gaushit::reco1`, `gaushit:uncalib:reco1` |
     * | `gaushit:uncalib:`              | `gaushit:uncalib:reco1` (newest process *with a match*) |
     * | `.*::`                          | `gaushit::reco2`, `fasthit::reco2` |
     * | `.*::reco.*`                    | `gaushit::reco1`, `gaushitStage2::reco1`, `gaushit::reco2`, `fasthit::reco2` |
     * | `::` (all empty)                | nothing (no product has an empty label) |
     *
     * If in some event `gaushit::reco2` is not available (e.g. `reco2` ran
     * `gaushit` behind a filter that rejected that event), `gaushit::` returns
     * `gaushit::reco1` for that event, like `e.getHandle<T>("gaushit")`.
     *
     * @note art's special process names `current_process` and `input_source`
     *       are not interpreted: they are treated as regex patterns.
     */
    template <typename T, typename Evt>
    std::vector<typename Evt::template HandleT<T>> getManyByRegexTag(Evt const& e,
                                                                     art::InputTag const& tag)
    {
      // An empty process selects the newest process with a matching product, as
      // e.getHandle<T>() does: the selector accepts any process, and the loop
      // keeps only the first (newest) process art returns.
      // An explicit process is a regex: products from all matching processes are kept.
      bool const only_newest_process = tag.process().empty();
      detail::RegexTagSelector const selector{
        tag.label(), tag.instance(), only_newest_process ? ".*" : tag.process()};

      // art returns the tags grouped by process, newest process first
      std::vector<typename Evt::template HandleT<T>> handles;
      std::string const* newest_process = nullptr;
      for (art::InputTag const& t : e.template getInputTags<T>(selector)) {
        auto handle = e.template getHandle<T>(t);
        if (!handle.isValid()) continue; // not available in this event: look further back, like art

        if (only_newest_process) {
          if (!newest_process)
            newest_process = &t.process();
          else if (t.process() != *newest_process)
            break; // older processes from here on
        }
        handles.push_back(std::move(handle));
      }

      return handles;
    }
  } // namespace util
} // namespace lar

#endif // LARDATA_LARDATA_ARTDATAHELPER_GETMANYBYREGEXTAG_H
