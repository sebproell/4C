// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_solid_3D_ele_surface.hpp"

#include "4C_comm_pack_helpers.hpp"
#include "4C_comm_utils_factory.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_fem_general_cell_type.hpp"
#include "4C_fem_general_utils_integration.hpp"
#include "4C_solid_3D_ele_line.hpp"

FOUR_C_NAMESPACE_OPEN

namespace
{
  Core::FE::CellType detect_celltype(
      Core::FE::CellType parent_cell_type, const Discret::Elements::SolidSurface& surface_element)
  {
    // if NURBS elements:
    if (parent_cell_type == Core::FE::CellType::nurbs8)
      return Core::FE::CellType::nurbs4;

    else if (parent_cell_type == Core::FE::CellType::nurbs27)
      return Core::FE::CellType::nurbs9;

    else
    {
      switch (surface_element.num_node())
      {
        case 3:
          return Core::FE::CellType::tri3;
        case 6:
        {
          if (parent_cell_type == Core::FE::CellType::tet10)
            return Core::FE::CellType::tri6;
          else if (parent_cell_type == Core::FE::CellType::hex18)
            return Core::FE::CellType::quad6;

          FOUR_C_THROW("what other surface element has 6 nodes???");
        }
        case 4:
          return Core::FE::CellType::quad4;
        case 8:
          return Core::FE::CellType::quad8;
        case 9:
          return Core::FE::CellType::quad9;
        default:
          FOUR_C_THROW("Could not infere cell type of surface element.");
      }
    }
  }

  Core::FE::GaussRule2D get_optimal_gauss_rule(Core::FE::CellType celltype)
  {
    switch (celltype)
    {
      case Core::FE::CellType::tri3:
        return Core::FE::GaussRule2D::tri_3point;
      case Core::FE::CellType::tri6:
        return Core::FE::GaussRule2D::tri_6point;
      case Core::FE::CellType::quad4:
        return Core::FE::GaussRule2D::quad_4point;
      case Core::FE::CellType::quad8:
      case Core::FE::CellType::quad9:
      case Core::FE::CellType::nurbs9:
        return Core::FE::GaussRule2D::quad_9point;
      case Core::FE::CellType::quad6:
        return Core::FE::GaussRule2D::quad_6point;
      default:
        FOUR_C_THROW("shape type unknown!\n");
    }
  }
}  // namespace

Discret::Elements::SolidSurfaceType Discret::Elements::SolidSurfaceType::instance_;

Discret::Elements::SolidSurfaceType& Discret::Elements::SolidSurfaceType::instance()
{
  return instance_;
}

Core::Communication::ParObject* Discret::Elements::SolidSurfaceType::create(
    Core::Communication::UnpackBuffer& buffer)
{
  auto* object = new Discret::Elements::SolidSurface(-1, -1);
  object->unpack(buffer);
  return object;
}

std::shared_ptr<Core::Elements::Element> Discret::Elements::SolidSurfaceType::create(
    const int id, const int owner)
{
  return nullptr;
}

Discret::Elements::SolidSurface::SolidSurface(int id, int owner, int nnode, const int* nodeids,
    Core::Nodes::Node** nodes, Core::Elements::Element* parent, const int lsurface)
    : Core::Elements::FaceElement(id, owner),
      distype_(Core::FE::CellType::dis_none),
      numdofpernode_(-1),
      gaussrule_(Core::FE::GaussRule2D::undefined)
{
  set_node_ids(nnode, nodeids);
  build_nodal_pointers(nodes);
  set_parent_master_element(parent, lsurface);

  numdofpernode_ = parent_master_element()->num_dof_per_node(*SolidSurface::nodes()[0]);
  // Safety check if all nodes have the same number of dofs!
  for (int nlid = 1; nlid < num_node(); ++nlid)
  {
    if (numdofpernode_ != parent_master_element()->num_dof_per_node(*SolidSurface::nodes()[nlid]))
      FOUR_C_THROW(
          "You need different NumDofPerNode for each node on this structural surface? ({} != {})",
          numdofpernode_, parent_master_element()->num_dof_per_node(*SolidSurface::nodes()[nlid]));
  }

  distype_ = detect_celltype(parent_element()->shape(), *this);
  gaussrule_ = get_optimal_gauss_rule(distype_);
}

Discret::Elements::SolidSurface::SolidSurface(int id, int owner)
    : Core::Elements::FaceElement(id, owner),
      distype_(Core::FE::CellType::dis_none),
      numdofpernode_(-1),
      gaussrule_(Core::FE::GaussRule2D::undefined)
{
}

Discret::Elements::SolidSurface::SolidSurface(const Discret::Elements::SolidSurface& old)
    : Core::Elements::FaceElement(old),
      distype_(old.distype_),
      numdofpernode_(old.numdofpernode_),
      gaussrule_(old.gaussrule_)
{
}

Core::Elements::Element* Discret::Elements::SolidSurface::clone() const
{
  auto* newelement = new Discret::Elements::SolidSurface(*this);
  return newelement;
}

Core::FE::CellType Discret::Elements::SolidSurface::shape() const { return distype_; }

void Discret::Elements::SolidSurface::pack(Core::Communication::PackBuffer& data) const
{
  // pack type of this instance of ParObject
  int type = unique_par_object_id();
  add_to_pack(data, type);
  // add base class Core::Elements::FaceElement
  Core::Elements::FaceElement::pack(data);
  // add distype_
  add_to_pack(data, distype_);
  // add numdofpernode_
  add_to_pack(data, numdofpernode_);
  // add gaussrule_
  add_to_pack(data, gaussrule_);
}

void Discret::Elements::SolidSurface::unpack(Core::Communication::UnpackBuffer& buffer)
{
  Core::Communication::extract_and_assert_id(buffer, unique_par_object_id());

  // extract base class Core::Elements::FaceElement
  Core::Elements::FaceElement::unpack(buffer);

  // distype_
  extract_from_pack(buffer, distype_);
  // numdofpernode_
  extract_from_pack(buffer, numdofpernode_);
  // gaussrule_
  extract_from_pack(buffer, gaussrule_);
}


void Discret::Elements::SolidSurface::print(std::ostream& os) const
{
  os << "SolidSurface ";
  Element::print(os);
}

std::vector<std::shared_ptr<Core::Elements::Element>> Discret::Elements::SolidSurface::lines()
{
  return Core::Communication::element_boundary_factory<Discret::Elements::SolidLine<3>,
      Discret::Elements::SolidSurface>(Core::Communication::buildLines, *this);
}

int Discret::Elements::SolidSurface::num_line() const
{
  return Core::FE::get_number_of_element_lines(shape());
}

FOUR_C_NAMESPACE_CLOSE
