/*---------------------------------------------------------------------------*/
/*! \file
\brief history pair struct for discrete element method (DEM) interactions
\level 3
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
#ifndef FOUR_C_PARTICLE_INTERACTION_DEM_HISTORY_PAIR_STRUCT_HPP
#define FOUR_C_PARTICLE_INTERACTION_DEM_HISTORY_PAIR_STRUCT_HPP

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "baci_config.hpp"

#include "baci_comm_parobject.hpp"

BACI_NAMESPACE_OPEN

/*---------------------------------------------------------------------------*
 | class declarations                                                        |
 *---------------------------------------------------------------------------*/
namespace PARTICLEINTERACTION
{
  //! struct to store tangential contact history of interacting particles
  struct DEMHistoryPairTangential final
  {
    //! tangential stick flag
    bool stick_ = true;

    //! tangential gap
    double gap_t_[3] = {0.0, 0.0, 0.0};

    //! pack history pair data
    void Pack(CORE::COMM::PackBuffer& data) const
    {
      data.AddtoPack(stick_);

      for (int i = 0; i < 3; ++i) data.AddtoPack(gap_t_[i]);
    }

    //! unpack history pair data
    void Unpack(std::vector<char>::size_type& position, const std::vector<char>& data)
    {
      CORE::COMM::ParObject::ExtractfromPack(position, data, stick_);

      for (int i = 0; i < 3; ++i) CORE::COMM::ParObject::ExtractfromPack(position, data, gap_t_[i]);
    }
  };

  //! struct to store rolling contact history of interacting particles
  struct DEMHistoryPairRolling final
  {
    //! rolling stick flag
    bool stick_ = true;

    //! rolling gap
    double gap_r_[3] = {0.0, 0.0, 0.0};

    //! pack history pair data
    void Pack(CORE::COMM::PackBuffer& data) const
    {
      data.AddtoPack(stick_);

      for (int i = 0; i < 3; ++i) data.AddtoPack(gap_r_[i]);
    }

    //! unpack history pair data
    void Unpack(std::vector<char>::size_type& position, const std::vector<char>& data)
    {
      CORE::COMM::ParObject::ExtractfromPack(position, data, stick_);

      for (int i = 0; i < 3; ++i) CORE::COMM::ParObject::ExtractfromPack(position, data, gap_r_[i]);
    }
  };

  //! struct to store adhesion history of interacting particles
  struct DEMHistoryPairAdhesion final
  {
    //! surface energy
    double surface_energy_ = 0.0;

    //! adhesion force
    double adhesion_force_ = 0.0;

    //! pack history pair data
    void Pack(CORE::COMM::PackBuffer& data) const
    {
      data.AddtoPack(surface_energy_);

      data.AddtoPack(adhesion_force_);
    }

    //! unpack history pair data
    void Unpack(std::vector<char>::size_type& position, const std::vector<char>& data)
    {
      CORE::COMM::ParObject::ExtractfromPack(position, data, surface_energy_);

      CORE::COMM::ParObject::ExtractfromPack(position, data, adhesion_force_);
    }
  };
}  // namespace PARTICLEINTERACTION

/*---------------------------------------------------------------------------*/
BACI_NAMESPACE_CLOSE

#endif