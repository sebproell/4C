/*----------------------------------------------------------------------*/
/*! \file

 \brief evaluation class containing routines for calculation of scalar transport
        within porous medium

\level 2

 *----------------------------------------------------------------------*/


#ifndef FOUR_C_SCATRA_ELE_CALC_PORO_HPP
#define FOUR_C_SCATRA_ELE_CALC_PORO_HPP

#include "4C_config.hpp"

#include "4C_scatra_ele_calc.hpp"

FOUR_C_NAMESPACE_OPEN

// forward declarations
namespace MAT
{
  class ScatraMat;
}

namespace DRT
{
  namespace ELEMENTS
  {
    class ScaTraEleDiffManagerPoro;


    template <CORE::FE::CellType distype>
    class ScaTraEleCalcPoro : public virtual ScaTraEleCalc<distype>
    {
     protected:
      /// (private) protected constructor, since we are a Singleton.
      ScaTraEleCalcPoro(const int numdofpernode, const int numscal, const std::string& disname);

     private:
      typedef ScaTraEleCalc<distype> my;
      using my::nen_;
      using my::nsd_;
      using my::nsd_ele_;

     public:
      /// Singleton access method
      static ScaTraEleCalcPoro<distype>* Instance(
          const int numdofpernode, const int numscal, const std::string& disname);

      /// Evaluate the element
      /*!
        Generic virtual interface function. Called via base pointer.
       */
      //   virtual int Evaluate(DRT::Element*                 ele,
      //                        Teuchos::ParameterList&       params,
      //                        DRT::Discretization &         discretization,
      //                        const std::vector<int> &      lm,
      //                        CORE::LINALG::SerialDenseMatrix&     elemat1_epetra,
      //                        CORE::LINALG::SerialDenseMatrix&     elemat2_epetra,
      //                        CORE::LINALG::SerialDenseVector&     elevec1_epetra,
      //                        CORE::LINALG::SerialDenseVector&     elevec2_epetra,
      //                        CORE::LINALG::SerialDenseVector&     elevec3_epetra);

     protected:
      //! evaluate action
      int EvaluateAction(DRT::Element* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, const SCATRA::Action& action,
          DRT::Element::LocationArray& la, CORE::LINALG::SerialDenseMatrix& elemat1_epetra,
          CORE::LINALG::SerialDenseMatrix& elemat2_epetra,
          CORE::LINALG::SerialDenseVector& elevec1_epetra,
          CORE::LINALG::SerialDenseVector& elevec2_epetra,
          CORE::LINALG::SerialDenseVector& elevec3_epetra) override;

      //   int EvaluateODMesh(  DRT::Element*                 ele,
      //                        Teuchos::ParameterList&       params,
      //                        DRT::Discretization &         discretization,
      //                        const std::vector<int> &      lm,
      //                        CORE::LINALG::SerialDenseMatrix&     elemat1_epetra,
      //                        CORE::LINALG::SerialDenseMatrix&     elemat2_epetra,
      //                        CORE::LINALG::SerialDenseVector&     elevec1_epetra,
      //                        CORE::LINALG::SerialDenseVector&     elevec2_epetra,
      //                        CORE::LINALG::SerialDenseVector&     elevec3_epetra);
      //
      //   int EvaluateODFluid(  DRT::Element*                 ele,
      //                        Teuchos::ParameterList&       params,
      //                        DRT::Discretization &         discretization,
      //                        const std::vector<int> &      lm,
      //                        CORE::LINALG::SerialDenseMatrix&     elemat1_epetra,
      //                        CORE::LINALG::SerialDenseMatrix&     elemat2_epetra,
      //                        CORE::LINALG::SerialDenseVector&     elevec1_epetra,
      //                        CORE::LINALG::SerialDenseVector&     elevec2_epetra,
      //                        CORE::LINALG::SerialDenseVector&     elevec3_epetra);
      //
      //   //! calculate matrix and rhs. Here the whole thing is hidden.
      //   virtual void SysmatODMesh(
      //     DRT::Element*                         ele,       //!< the element we are dealing with
      //     CORE::LINALG::SerialDenseMatrix&             emat,      //!< element matrix to
      //     calculate const int                     numdofpernode
      //   );
      //
      //   //! calculate matrix and rhs. Here the whole thing is hidden.
      //   virtual void SysmatODFluid(
      //     DRT::Element*                         ele,       //!< the element we are dealing with
      //     CORE::LINALG::SerialDenseMatrix&             emat,      //!< element matrix to
      //     calculate const int                     numdofpernode
      //   );

      //! read element coordinates
      void ReadElementCoordinates(const DRT::Element* ele) override;

      //! extract element based or nodal values
      //  return extracted values of phinp
      void ExtractElementAndNodeValues(DRT::Element* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, DRT::Element::LocationArray& la) override;

      //! extract element based or nodal values
      //  return extracted values of phinp
      virtual void ExtractElementAndNodeValuesPoro(DRT::Element* ele,
          Teuchos::ParameterList& params, DRT::Discretization& discretization,
          DRT::Element::LocationArray& la);

      //! get the material parameters
      void GetMaterialParams(const DRT::Element* ele,  //!< the element we are dealing with
          std::vector<double>& densn,                  //!< density at t_(n)
          std::vector<double>& densnp,                 //!< density at t_(n+1) or t_(n+alpha_F)
          std::vector<double>& densam,                 //!< density at t_(n+alpha_M)
          double& visc,                                //!< fluid viscosity
          const int iquad = -1                         //!< id of current gauss point (default = -1)
          ) override;

      //! compute porosity based on solid, fluid and (potentially) scatra solution
      virtual void ComputePorosity(const DRT::Element* ele  //!< the element we are dealing with
      );

      //! compute pore pressure
      virtual double ComputePorePressure();

      //! material ScaTra
      void MatScaTra(
          const Teuchos::RCP<const CORE::MAT::Material> material,  //!< pointer to current material
          const int k,                                             //!< id of current scalar
          double& densn,                                           //!< density at t_(n)
          double& densnp,       //!< density at t_(n+1) or t_(n+alpha_F)
          double& densam,       //!< density at t_(n+alpha_M)
          double& visc,         //!< fluid viscosity
          const int iquad = -1  //!< id of current gauss point (default = -1)
          ) override;


      //! set diffusivity for poro scatra problem (i.e. scale by porosity)
      virtual void SetDiffusivity(
          const Teuchos::RCP<const MAT::ScatraMat>& material, const int k, const double scale);

      //! set densisties for poro scatra problem (i.e. scale by porosity)
      virtual void SetDensities(double porosity,
          double& densn,   //!< density at t_(n)
          double& densnp,  //!< density at t_(n+1) or t_(n+alpha_F)
          double& densam   //!< density at t_(n+alpha_M));
      );

      //! calculate scalar(s) and domain integral
      void CalculateScalars(const DRT::Element* ele, CORE::LINALG::SerialDenseVector& scalars,
          bool inverting, bool calc_grad_phi) override;


      //! get poro diffusion manager
      Teuchos::RCP<ScaTraEleDiffManagerPoro> DiffManager()
      {
        return Teuchos::rcp_static_cast<ScaTraEleDiffManagerPoro>(my::diffmanager_);
      };

      /*========================================================================*/
      //! @name Galerkin approximation and related
      /*========================================================================*/

      //! initial node coordinates
      CORE::LINALG::Matrix<nsd_, nen_> xyze0_;

      //! nodal porosity values at t_(n+1)
      CORE::LINALG::Matrix<nen_, 1> eporosity_;

      //! flag indacting a node based porosity
      bool isnodalporosity_;
    };

    /// ScaTraEleDiffManagerPoro implementation
    /*!
      This class keeps all poro-specific transport parameter needed for the evaluation of an
      element. The ScaTraEleDiffManagerPoro is derived from the standard ScaTraEleDiffManager.
    */

    ////TODO: HACK
    // const int NO_CONVECTION_NR = 6;

    class ScaTraEleDiffManagerPoro : public ScaTraEleDiffManager
    {
     public:
      ScaTraEleDiffManagerPoro(int numscal) : ScaTraEleDiffManager(numscal), porosity_(0.0)
      {
        return;
      }

      void SetPorosity(double porosity) { porosity_ = porosity; }

      double GetPorosity(const int k) const
      {
        //      if(k<NO_CONVECTION_NR)
        return porosity_;
        //      else
        //        return 1.0;
      }

     protected:
      //! porosity at gauss point
      double porosity_;
    };

    template <int NSD, int NEN>
    class ScaTraEleInternalVariableManagerPoro : public ScaTraEleInternalVariableManager<NSD, NEN>
    {
      typedef ScaTraEleInternalVariableManager<NSD, NEN> my;

     public:
      ScaTraEleInternalVariableManagerPoro(int numscal)
          : ScaTraEleInternalVariableManager<NSD, NEN>(numscal),
            zeroconvel_(true),
            zeroconv_(true),
            zero_(0.0)
      {
        return;
      }

      virtual ~ScaTraEleInternalVariableManagerPoro() = default;

      /*========================================================================*/
      //! @name return methods for internal variables
      /*========================================================================*/

      //! return convective velocity
      virtual const CORE::LINALG::Matrix<NSD, 1>& ConVel(const int k) const {
          //    if(k<NO_CONVECTION_NR)
          {return my::convelint_;
    }
    //    else
    //      return zeroconvel_;
  };  // namespace ELEMENTS
  //! return convective part in convective form
  virtual const CORE::LINALG::Matrix<NEN, 1>& Conv(const int k) const
  {  //    if(k<NO_CONVECTION_NR)
    {
      return my::conv_;
    }  // namespace DRT
       //    else
       //      return zeroconv_;
  }

  //! return convective term of current scalar value
  virtual const double& ConvPhi(const int k) const
  {  //    if(k<NO_CONVECTION_NR)
    {
      return my::conv_phi_[k];
    }
    //    else
    //      return zero_;
  }

 private:
  CORE::LINALG::Matrix<NSD, 1> zeroconvel_;
  CORE::LINALG::Matrix<NEN, 1> zeroconv_;
  double zero_;
};  // namespace DRT
}
}


FOUR_C_NAMESPACE_CLOSE

#endif