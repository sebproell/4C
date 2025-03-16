// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_ale_ale2_nurbs.hpp"

FOUR_C_NAMESPACE_OPEN

Discret::Elements::Nurbs::Ale2NurbsType Discret::Elements::Nurbs::Ale2NurbsType::instance_;

Discret::Elements::Nurbs::Ale2NurbsType& Discret::Elements::Nurbs::Ale2NurbsType::instance()
{
  return instance_;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Core::Communication::ParObject* Discret::Elements::Nurbs::Ale2NurbsType::create(
    Core::Communication::UnpackBuffer& buffer)
{
  Discret::Elements::Nurbs::Ale2Nurbs* object = new Discret::Elements::Nurbs::Ale2Nurbs(-1, -1);
  object->unpack(buffer);
  return object;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
std::shared_ptr<Core::Elements::Element> Discret::Elements::Nurbs::Ale2NurbsType::create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "ALE2")
  {
    if (eledistype == "NURBS4" || eledistype == "NURBS9")
    {
      return std::make_shared<Discret::Elements::Nurbs::Ale2Nurbs>(id, owner);
    }
  }
  return nullptr;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
std::shared_ptr<Core::Elements::Element> Discret::Elements::Nurbs::Ale2NurbsType::create(
    const int id, const int owner)
{
  return std::make_shared<Discret::Elements::Nurbs::Ale2Nurbs>(id, owner);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Discret::Elements::Nurbs::Ale2Nurbs::Ale2Nurbs(int id, int owner)
    : Discret::Elements::Ale2::Ale2(id, owner)
{
  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Discret::Elements::Nurbs::Ale2Nurbs::Ale2Nurbs(const Discret::Elements::Nurbs::Ale2Nurbs& old)
    : Discret::Elements::Ale2::Ale2(old)
{
  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void Discret::Elements::Nurbs::Ale2Nurbs::print(std::ostream& os) const
{
  os << "Ale2Nurbs ";
  Element::print(os);
  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Core::FE::CellType Discret::Elements::Nurbs::Ale2Nurbs::shape() const
{
  switch (num_node())
  {
    case 4:
      return Core::FE::CellType::nurbs4;
    case 9:
      return Core::FE::CellType::nurbs9;
    default:
      FOUR_C_THROW("unexpected number of nodes {}", num_node());
      break;
  }
}

FOUR_C_NAMESPACE_CLOSE
