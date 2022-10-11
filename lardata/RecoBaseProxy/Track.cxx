/**
 * @file   lardata/RecoBaseProxy/Track.cxx
 * @brief  Implementation file for `proxy::Track`.
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   July 27, 2017
 * @see    Track.h
 *
 */

// LArSoft libraries
#include "lardata/RecoBaseProxy/Track.h" // proxy namespace

#if 0
//------------------------------------------------------------------------------
//---  proxy::TrackPointIterator implementation
//------------------------------------------------------------------------------
proxy::TrackPointIterator proxy::details::TrackPointIteratorBox::begin() const
  { return TrackPointIterator(*track, 0U); }
proxy::TrackPointIterator proxy::details::TrackPointIteratorBox::end() const
  { return TrackPointIterator(*track, track->track().NPoints()); }

//------------------------------------------------------------------------------
//---  proxy::Tracks implementation
//------------------------------------------------------------------------------
proxy::Track proxy::Tracks::getProxyAt(std::size_t index) const {

  return { Base_t::getMainAt(index), hitsPerTrack[index] };

} // proxy::Tracks::getProxyAt()

#endif // 0
