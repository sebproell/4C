// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_mat_anisotropy.hpp"

#include "4C_comm_pack_helpers.hpp"
#include "4C_comm_parobject.hpp"
#include "4C_fem_general_fiber_node_holder.hpp"
#include "4C_mat_anisotropy_extension.hpp"
#include "4C_mat_anisotropy_utils.hpp"
#include "4C_utils_shared_ptr_from_ref.hpp"

#include <Teuchos_ParameterList.hpp>

#include <memory>

FOUR_C_NAMESPACE_OPEN

Mat::Anisotropy::Anisotropy()
    : element_fibers_initialized_(false),
      gp_fibers_initialized_(false),
      element_fibers_(0),
      gp_fibers_(0),
      gp_cylinder_coordinate_system_managers_(0),
      extensions_(0)
{
  // empty
}

void Mat::Anisotropy::pack_anisotropy(Core::Communication::PackBuffer& data) const
{
  add_to_pack(data, numgp_);
  add_to_pack(data, element_fibers_initialized_);
  add_to_pack(data, gp_fibers_initialized_);
  add_to_pack(data, element_fibers_);
  pack_fiber_vector<Core::LinAlg::Matrix<3, 1>>(data, gp_fibers_);

  if (element_cylinder_coordinate_system_manager_)
  {
    add_to_pack(data, true);
    element_cylinder_coordinate_system_manager_->pack(data);
  }
  else
  {
    add_to_pack(data, false);
  }

  for (const auto& gpCylinderCoordinateSystemManager : gp_cylinder_coordinate_system_managers_)
  {
    gpCylinderCoordinateSystemManager.pack(data);
  }
}

void Mat::Anisotropy::unpack_anisotropy(Core::Communication::UnpackBuffer& buffer)
{
  extract_from_pack(buffer, numgp_);
  extract_from_pack(buffer, element_fibers_initialized_);
  extract_from_pack(buffer, gp_fibers_initialized_);
  extract_from_pack(buffer, element_fibers_);
  unpack_fiber_vector<Core::LinAlg::Matrix<3, 1>>(buffer, gp_fibers_);

  bool have_element_cylinder_coordinate_system_manager;
  extract_from_pack(buffer, have_element_cylinder_coordinate_system_manager);
  if (have_element_cylinder_coordinate_system_manager)
  {
    element_cylinder_coordinate_system_manager_ = CylinderCoordinateSystemManager();
    element_cylinder_coordinate_system_manager_->unpack(buffer);
  }
  else
  {
    element_cylinder_coordinate_system_manager_ = std::nullopt;
  }

  for (auto& gpCylinderCoordinateSystemManager : gp_cylinder_coordinate_system_managers_)
  {
    gpCylinderCoordinateSystemManager.unpack(buffer);
  }
}

void Mat::Anisotropy::set_number_of_gauss_points(int numgp) { numgp_ = numgp; }

void Mat::Anisotropy::read_anisotropy_from_element(
    const Core::IO::InputParameterContainer& container)
{
  // This method might be called even though there is no anisotropy possible for the element.
  // In this case, we just return without doing anything.
  using OptionalFiber = std::optional<std::vector<double>>;
  const bool any_fibers_found =
      container.get_if<OptionalFiber>("FIBER1") or container.get_if<OptionalFiber>("RAD") or
      container.get_if<OptionalFiber>("AXI") or container.get_if<OptionalFiber>("CIR");
  if (!any_fibers_found) return;


  // Read coordinate system
  if (container.get<std::optional<std::vector<double>>>("RAD").has_value() and
      container.get<std::optional<std::vector<double>>>("AXI").has_value() and
      container.get<std::optional<std::vector<double>>>("CIR").has_value())
  {
    // read fibers in RAD AXI CIR notation
    if (!element_cylinder_coordinate_system_manager_)
    {
      element_cylinder_coordinate_system_manager_ = CylinderCoordinateSystemManager();
    }

    element_cylinder_coordinate_system_manager_->read_from_element_line_definition(container);
  }

  else
  {
    // read fibers in FIBERi notation and determine number of fibers
    for (int fiber_index = 1;; ++fiber_index)
    {
      std::string fiber_name = "FIBER" + std::to_string(fiber_index);

      // We count up until we hit a fiber that we do not know about. Thus we need to use the
      // get_if mechanism.
      auto* fiber_ptr = container.get_if<std::optional<std::vector<double>>>(fiber_name);
      if (!fiber_ptr or !fiber_ptr->has_value())
      {
        break;
      }

      element_fibers_.resize(fiber_index);
      read_anisotropy_fiber(container, fiber_name, element_fibers_[fiber_index - 1]);
    }
  }
  on_element_fibers_initialized();
}

void Mat::Anisotropy::read_anisotropy_from_parameter_list(const Teuchos::ParameterList& params)
{
  if (params.isParameter("fiberholder"))
  {
    const auto& fiberHolder = params.get<Core::Nodes::NodalFiberHolder>("fiberholder");

    gp_fibers_.resize(numgp_);

    for (const auto& fiber : fiberHolder.get_fibers())
    {
      insert_fibers(fiber);
    }
  }

  on_gp_fibers_initialized();
}

void Mat::Anisotropy::insert_fibers(std::vector<Core::LinAlg::Matrix<3, 1>> fiber)
{
  for (unsigned gp = 0; gp < numgp_; ++gp)
  {
    gp_fibers_[gp].emplace_back(fiber[gp]);
  }
}

void Mat::Anisotropy::set_element_fibers(const std::vector<Core::LinAlg::Matrix<3, 1>>& fibers)
{
  element_fibers_ = fibers;

  on_element_fibers_initialized();
}

void Mat::Anisotropy::set_gauss_point_fibers(
    const std::vector<std::vector<Core::LinAlg::Matrix<3, 1>>>& fibers)
{
  // check input fibers whether they make sense

  // Check whether the size of the first vector is the number of Gauss points
  if (fibers.size() != numgp_)
  {
    FOUR_C_THROW("The Gauss point fibers don't have the expected size of {} ({} given).", numgp_,
        fibers.size());
  }

  // Check whether every second vector have the same lengths
  unsigned num_fibs = 1;
  unsigned i = 0;
  for (const auto& gpfibers : fibers)
  {
    if (i == 0)
    {
      num_fibs = gpfibers.size();
    }
    else if (num_fibs != gpfibers.size())
    {
      FOUR_C_THROW(
          "The size of the Gauss point do not match! At every Gauss point, the same amount of "
          "fibers are necessary. Error occurred at Gauss point {}. Expected {} fibers, but got {}.",
          i, num_fibs, gpfibers.size());
    }
  }

  gp_fibers_ = fibers;

  on_gp_fibers_initialized();
}

const Core::LinAlg::Matrix<3, 1>& Mat::Anisotropy::get_element_fiber(unsigned int i) const
{
  if (!element_fibers_initialized_)
  {
    FOUR_C_THROW("The element fibers are not yet initialized.");
  }
  if (i >= element_fibers_.size())
  {
    FOUR_C_THROW(
        "You requested fiber {}, but only {} fibers are available", i + 1, element_fibers_.size());
  }
  return element_fibers_[i];
}

const std::vector<Core::LinAlg::Matrix<3, 1>>& Mat::Anisotropy::get_element_fibers() const
{
  if (!element_fibers_initialized_)
  {
    FOUR_C_THROW("The element fibers are not yet initialized.");
  }
  return element_fibers_;
}

const std::vector<std::vector<Core::LinAlg::Matrix<3, 1>>>&
Mat::Anisotropy::get_gauss_point_fibers() const
{
  if (!gp_fibers_initialized_)
  {
    FOUR_C_THROW("The Gauss point fibers are not yet initialized.");
  }
  return gp_fibers_;
}

const Core::LinAlg::Matrix<3, 1>& Mat::Anisotropy::get_gauss_point_fiber(
    unsigned int gp, unsigned int i) const
{
  if (!gp_fibers_initialized_)
  {
    FOUR_C_THROW("The GP fibers are not yet initialized.");
  }

  if (gp >= gp_fibers_.size())
  {
    FOUR_C_THROW("The number of GP is too large. {} instead of maximum allowed {}", gp + 1,
        gp_fibers_.size());
  }

  if (i >= gp_fibers_[gp].size())
  {
    FOUR_C_THROW("You requested fiber {} fibers are available", i + 1, element_fibers_.size());
  }
  return gp_fibers_[gp][i];
}

void Mat::Anisotropy::register_anisotropy_extension(BaseAnisotropyExtension& extension)
{
  extensions_.emplace_back(Core::Utils::shared_ptr_from_ref(extension));
  extension.set_anisotropy(*this);
}

void Mat::Anisotropy::on_element_fibers_initialized()
{
  element_fibers_initialized_ = true;
  for (auto& extension : extensions_)
  {
    extension->on_global_element_data_initialized();
  }

  if (element_fibers_initialized_ and gp_fibers_initialized_)
  {
    for (auto& extension : extensions_)
    {
      extension->on_global_data_initialized();
    }
  }
}

void Mat::Anisotropy::on_gp_fibers_initialized()
{
  gp_fibers_initialized_ = true;
  for (auto& extension : extensions_)
  {
    extension->on_global_gp_data_initialized();
  }

  if (element_fibers_initialized_ and gp_fibers_initialized_)
  {
    for (auto& extension : extensions_)
    {
      extension->on_global_data_initialized();
    }
  }
}

int Mat::Anisotropy::get_number_of_gauss_points() const { return numgp_; }

int Mat::Anisotropy::get_number_of_element_fibers() const { return element_fibers_.size(); }

int Mat::Anisotropy::get_number_of_gauss_point_fibers() const
{
  if (gp_fibers_.empty()) return 0;

  return gp_fibers_[0].size();
}

bool Mat::Anisotropy::has_element_cylinder_coordinate_system() const
{
  return element_cylinder_coordinate_system_manager_.has_value();
}

bool Mat::Anisotropy::has_gp_cylinder_coordinate_system() const
{
  return !gp_cylinder_coordinate_system_managers_.empty();
}
FOUR_C_NAMESPACE_CLOSE
