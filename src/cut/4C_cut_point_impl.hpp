// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_CUT_POINT_IMPL_HPP
#define FOUR_C_CUT_POINT_IMPL_HPP

#include "4C_config.hpp"

#include "4C_cut_line.hpp"

FOUR_C_NAMESPACE_OPEN


namespace Cut
{
  namespace Impl
  {
    class LineBetweenFilter
    {
     public:
      LineBetweenFilter(Point* me, Point* other) : me_(me), other_(other) {}

      bool operator()(Line* line) { return line->between(me_, other_); }

     private:
      Point* me_;
      Point* other_;
    };

    class LineHasSideFilter
    {
     public:
      explicit LineHasSideFilter(Side* side) : side_(side) {}

      /// true if the line is cut by the side but not on any side's edges
      bool operator()(Line* line) { return line->is_internal_cut(side_); }

     private:
      Side* side_;
    };

    class NextLineOnElementCutFilter
    {
     public:
      NextLineOnElementCutFilter(Line* line, Side* side, Element* element)
          : line_(line), side_(side), element_(element)
      {
      }

      bool operator()(Line* line)
      {
        return line != line_ and line->is_cut(side_) and
               (element_ == nullptr or line->is_cut(element_));
      }

     private:
      Line* line_;
      Side* side_;
      Element* element_;
    };

  }  // namespace Impl

  template <class Filter>
  Line* Point::find(Filter& filter, bool unique)
  {
    Line* line_found = nullptr;
    for (plain_line_set::iterator i = lines_.begin(); i != lines_.end(); ++i)
    {
      Line* line = *i;
      if (filter(line))
      {
        if (line_found == nullptr)
        {
          line_found = line;
          if (not unique)
          {
            break;
          }
        }
        else
        {
          FOUR_C_THROW("not unique");
        }
      }
    }
    return line_found;
  }

  inline Line* Point::common_line(Point* other)
  {
    Impl::LineBetweenFilter filter(this, other);
    return find(filter, true);
  }

  inline Line* Point::cut_line(Side* side, bool unique)
  {
    Impl::LineHasSideFilter filter(side);
    return find(filter, unique);
  }

  inline Line* Point::cut_line(Line* line, Side* side, Element* element)
  {
    Impl::NextLineOnElementCutFilter filter(line, side, element);
    return find(filter, true);
  }

}  // namespace Cut


FOUR_C_NAMESPACE_CLOSE

#endif
