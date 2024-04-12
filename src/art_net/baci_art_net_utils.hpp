/*----------------------------------------------------------------------*/
/*! \file

 \brief helper functions/classes for artery problems

   \level 3

 *----------------------------------------------------------------------*/

#ifndef FOUR_C_ART_NET_UTILS_HPP
#define FOUR_C_ART_NET_UTILS_HPP

#include "baci_config.hpp"

#include "baci_inpar_bio.hpp"
#include "baci_io.hpp"

BACI_NAMESPACE_OPEN

// forward declarations
namespace ADAPTER
{
  class ArtNet;
}

namespace DRT
{
  class Discretization;
  class Element;
}  // namespace DRT
namespace CORE::LINALG
{
  class Solver;
}

namespace ART
{
  //! \brief implementation of special clone strategy for automatic generation
  //!        of scatra discretization from a given artery discretization
  class ArteryScatraCloneStrategy
  {
   public:
    /// constructor
    explicit ArteryScatraCloneStrategy() {}
    /// destructor
    virtual ~ArteryScatraCloneStrategy() = default;

   protected:
    /// determine element type string and whether element is copied or not
    bool DetermineEleType(
        DRT::Element* actele, const bool ismyele, std::vector<std::string>& eletype);

    //! provide cloned element with element specific data (material etc.)
    void SetElementData(
        Teuchos::RCP<DRT::Element> newele,  //! current cloned element on target discretization
        DRT::Element* oldele,               //! current element on source discretization
        const int matid,                    //! material of cloned element
        const bool isnurbs                  //! nurbs flag
    );

    void CheckMaterialType(const int matid);

    /// returns conditions names to be copied (source and target name)
    std::map<std::string, std::string> ConditionsToCopy() const;
  };  // class ArteryScatraCloneStrategy

  namespace UTILS
  {
    // create algorithm depending on time integration scheme
    Teuchos::RCP<ADAPTER::ArtNet> CreateAlgorithm(INPAR::ARTDYN::TimeIntegrationScheme timintscheme,
        Teuchos::RCP<DRT::Discretization> dis, const int linsolvernumber,
        const Teuchos::ParameterList& probparams, const Teuchos::ParameterList& artparams,
        Teuchos::RCP<IO::DiscretizationWriter> output);

    //! exchange material pointers of discretizations
    void AssignMaterialPointers(
        const std::string& artery_disname, const std::string& scatra_disname);

    //! set material pointers
    void SetMaterialPointersMatchingGrid(Teuchos::RCP<const DRT::Discretization> sourcedis,
        Teuchos::RCP<const DRT::Discretization> targetdis);
  }  // namespace UTILS
}  // namespace ART



BACI_NAMESPACE_CLOSE

#endif