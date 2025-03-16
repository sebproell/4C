// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_mat_particle_sph_fluid.hpp"

#include "4C_comm_pack_helpers.hpp"
#include "4C_global_data.hpp"
#include "4C_mat_par_bundle.hpp"

FOUR_C_NAMESPACE_OPEN

/*---------------------------------------------------------------------------*
 | define static class member                                                |
 *---------------------------------------------------------------------------*/
Mat::ParticleMaterialSPHFluidType Mat::ParticleMaterialSPHFluidType::instance_;

/*---------------------------------------------------------------------------*
 | constructor                                                               |
 *---------------------------------------------------------------------------*/
Mat::PAR::ParticleMaterialSPHFluid::ParticleMaterialSPHFluid(
    const Core::Mat::PAR::Parameter::Data& matdata)
    : Parameter(matdata),
      ParticleMaterialBase(matdata),
      ParticleMaterialThermo(matdata),
      refDensFac_(matdata.parameters.get<double>("REFDENSFAC")),
      exponent_(matdata.parameters.get<double>("EXPONENT")),
      backgroundPressure_(matdata.parameters.get<double>("BACKGROUNDPRESSURE")),
      bulkModulus_(matdata.parameters.get<double>("BULK_MODULUS")),
      dynamicViscosity_(matdata.parameters.get<double>("DYNAMIC_VISCOSITY")),
      bulkViscosity_(matdata.parameters.get<double>("BULK_VISCOSITY")),
      artificialViscosity_(matdata.parameters.get<double>("ARTIFICIAL_VISCOSITY"))
{
  // empty constructor
}

/*---------------------------------------------------------------------------*
 | create material instance of matching type with parameters                 |
 *---------------------------------------------------------------------------*/
std::shared_ptr<Core::Mat::Material> Mat::PAR::ParticleMaterialSPHFluid::create_material()
{
  return std::make_shared<Mat::ParticleMaterialSPHFluid>(this);
}

/*---------------------------------------------------------------------------*
 *---------------------------------------------------------------------------*/
Core::Communication::ParObject* Mat::ParticleMaterialSPHFluidType::create(
    Core::Communication::UnpackBuffer& buffer)
{
  Mat::ParticleMaterialSPHFluid* particlematsph = new Mat::ParticleMaterialSPHFluid();
  particlematsph->unpack(buffer);
  return particlematsph;
}

/*---------------------------------------------------------------------------*
 | constructor (empty material object)                                       |
 *---------------------------------------------------------------------------*/
Mat::ParticleMaterialSPHFluid::ParticleMaterialSPHFluid() : params_(nullptr)
{
  // empty constructor
}

/*---------------------------------------------------------------------------*
 | constructor (with given material parameters)                              |
 *---------------------------------------------------------------------------*/
Mat::ParticleMaterialSPHFluid::ParticleMaterialSPHFluid(Mat::PAR::ParticleMaterialSPHFluid* params)
    : params_(params)
{
  // empty constructor
}

/*---------------------------------------------------------------------------*
 | pack                                                                      |
 *---------------------------------------------------------------------------*/
void Mat::ParticleMaterialSPHFluid::pack(Core::Communication::PackBuffer& data) const
{
  // pack type of this instance of ParObject
  int type = unique_par_object_id();
  add_to_pack(data, type);

  // matid
  int matid = -1;
  if (params_ != nullptr) matid = params_->id();  // in case we are in post-process mode
  add_to_pack(data, matid);
}

/*---------------------------------------------------------------------------*
 | unpack                                                                    |
 *---------------------------------------------------------------------------*/
void Mat::ParticleMaterialSPHFluid::unpack(Core::Communication::UnpackBuffer& buffer)
{
  Core::Communication::extract_and_assert_id(buffer, unique_par_object_id());

  // matid and recover params_
  int matid;
  extract_from_pack(buffer, matid);
  params_ = nullptr;
  if (Global::Problem::instance()->materials() != nullptr)
    if (Global::Problem::instance()->materials()->num() != 0)
    {
      // note: dynamic_cast needed due diamond inheritance structure
      const int probinst = Global::Problem::instance()->materials()->get_read_from_problem();
      Core::Mat::PAR::Parameter* mat =
          Global::Problem::instance(probinst)->materials()->parameter_by_id(matid);
      if (mat->type() == material_type())
        params_ = dynamic_cast<Mat::PAR::ParticleMaterialSPHFluid*>(mat);
      else
        FOUR_C_THROW("Type of parameter material {} does not fit to calling type {}", mat->type(),
            material_type());
    }
}

FOUR_C_NAMESPACE_CLOSE
