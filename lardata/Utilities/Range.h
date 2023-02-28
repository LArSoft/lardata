/**
 * \file Range.h
 *
 * \ingroup RangeTool
 *
 * \brief Class def header for a class Range
 *
 * @author kazuhiro
 */

/** \addtogroup RangeTool
    @{*/
#ifndef RANGE_H
#define RANGE_H

#include <algorithm>
#include <functional>
#include <utility>

namespace util {

  template <class T>
  class UniqueRangeSet;

  /**
     \class Range
     @brief represents a "Range" w/ notion of ordering.
     A range is defined by a pair of "start" and "end" values. This is stored in std::pair    \n
     attribute util::Range::_window. This attribute is protected so that the start/end cannot \n
     be changed w/o a check that start is always less than end. Note the specialization       \n
     requires a template class T to have less operator implemented.   \n
  */
  template <class T>
  class Range {
    // Make it a friend so UniqueRangeSet can access protected guys
    friend class UniqueRangeSet<T>;

  private:
    /// Default ctor is hidden
    Range() {}

  public:
    /// Enforced ctor. start must be less than end.
    Range(const T& start, const T& end) : _window(start, end)
    {
      if (start > end) throw std::runtime_error("Inserted invalid range: end before start.");
    }

    /// Default dtor
    ~Range() {}

    /// "start" accessor
    const T& Start() const { return _window.first; }
    /// "end" accessor
    const T& End() const { return _window.second; }
    /// Setter
    void Set(const T& s, const T& e)
    {
      if (s >= e) throw std::runtime_error("Inserted invalid range: end before start.");
      _window.first = s;
      _window.second = e;
    }

    //
    // Ordering w/ another Range
    //
    inline bool operator<(const Range& rhs) const { return (_window.second < rhs.Start()); }
    inline bool operator>(const Range& rhs) const { return (_window.first > rhs.End()); }
    inline bool operator==(const Range& rhs) const
    {
      return (_window.first == rhs.Start() && _window.second == rhs.End());
    }
    inline bool operator!=(const Range& rhs) const { return !((*this) == rhs); }

    //
    // Ordering w/ T
    //
    inline bool operator<(const T& rhs) const { return (_window.second < rhs); }
    inline bool operator>(const T& rhs) const { return (_window.first > rhs); }

    /// Merge two util::Range into 1
    void Merge(const Range& a)
    {
      _window.first = std::min(_window.first, a.Start());
      _window.second = std::max(_window.second, a.End());
    }

  protected:
    /// Protected to avoid user's illegal modification on first/second (sorry users!)
    std::pair<T, T> _window;
  };
}

namespace std {
  // Implement pointer comparison in case it's useful
  template <class T>
  /**
     \class less
     Implementation of std::less for util::Range pointers
   */
  class less<util::Range<T>*> {
  public:
    bool operator()(const util::Range<T>* lhs, const util::Range<T>* rhs)
    {
      return (*lhs) < (*rhs);
    }
  };
}

#endif
/** @} */ // end of doxygen group
