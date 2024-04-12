/*----------------------------------------------------------------------*/
/*! \file
 \brief  factory class providing the implementations of the porous fluid multiphase
         boundary element

   \level 3

 *----------------------------------------------------------------------*/

#ifndef FOUR_C_POROFLUIDMULTIPHASE_ELE_BOUNDARY_FACTORY_HPP
#define FOUR_C_POROFLUIDMULTIPHASE_ELE_BOUNDARY_FACTORY_HPP


#include "baci_config.hpp"

#include "baci_lib_element.hpp"

BACI_NAMESPACE_OPEN

namespace DRT
{
  namespace ELEMENTS
  {
    // forward declaration
    class PoroFluidMultiPhaseEleInterface;

    class PoroFluidMultiPhaseBoundaryFactory
    {
     public:
      //! ctor
      PoroFluidMultiPhaseBoundaryFactory() { return; };

      //! dtor
      virtual ~PoroFluidMultiPhaseBoundaryFactory() = default;

      //! ProvideImpl
      static PoroFluidMultiPhaseEleInterface* ProvideImpl(
          const DRT::Element* ele, const int numdofpernode, const std::string& disname);

     private:
      //! return instance of element evaluation class depending on implementation type
      template <CORE::FE::CellType distype>
      static PoroFluidMultiPhaseEleInterface* DefineProblemType(
          const int numdofpernode, const std::string& disname);
    };  // class PoroFluidMultiPhaseBoundaryFactory
  }     // namespace ELEMENTS
}  // namespace DRT


BACI_NAMESPACE_CLOSE

#endif