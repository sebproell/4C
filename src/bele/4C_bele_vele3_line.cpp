// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_bele_vele3.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_linalg_utils_sparse_algebra_math.hpp"
#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN



Discret::Elements::Vele3LineType Discret::Elements::Vele3LineType::instance_;

Discret::Elements::Vele3LineType& Discret::Elements::Vele3LineType::instance() { return instance_; }

/*----------------------------------------------------------------------*
 |  ctor (public)                                            gammi 04/07|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
Discret::Elements::Vele3Line::Vele3Line(int id, int owner, int nnode, const int* nodeids,
    Core::Nodes::Node** nodes, Core::Elements::Element* parent, const int lline)
    : Core::Elements::FaceElement(id, owner)
{
  set_node_ids(nnode, nodeids);
  build_nodal_pointers(nodes);
  set_parent_master_element(parent, lline);
  return;
}


/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 01/07|
 *----------------------------------------------------------------------*/
Discret::Elements::Vele3Line::Vele3Line(const Discret::Elements::Vele3Line& old)
    : Core::Elements::FaceElement(old)
{
  return;
}


/*----------------------------------------------------------------------*
 |  Deep copy this instance return pointer to it               (public) |
 |                                                            gee 01/07 |
 *----------------------------------------------------------------------*/
Core::Elements::Element* Discret::Elements::Vele3Line::clone() const
{
  Discret::Elements::Vele3Line* newelement = new Discret::Elements::Vele3Line(*this);
  return newelement;
}


/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                          u.kue 03/07 |
 *----------------------------------------------------------------------*/
Core::FE::CellType Discret::Elements::Vele3Line::shape() const
{
  switch (num_node())
  {
    case 2:
      return Core::FE::CellType::line2;
    case 3:
      return Core::FE::CellType::line3;
    default:
      FOUR_C_THROW("unexpected number of nodes {}", num_node());
  }
}


/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
void Discret::Elements::Vele3Line::pack(Core::Communication::PackBuffer& data) const
{
  FOUR_C_THROW("this Vele3Line element does not support communication");

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
void Discret::Elements::Vele3Line::unpack(Core::Communication::UnpackBuffer& buffer)
{
  FOUR_C_THROW("this Vele3Line element does not support communication");
  return;
}



/*----------------------------------------------------------------------*
 |  print this element (public)                              mwgee 01/07|
 *----------------------------------------------------------------------*/
void Discret::Elements::Vele3Line::print(std::ostream& os) const
{
  os << "Vele3Line ";
  Element::print(os);
  return;
}

FOUR_C_NAMESPACE_CLOSE
