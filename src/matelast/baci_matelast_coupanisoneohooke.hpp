/*----------------------------------------------------------------------*/
/*! \file
\brief Definition of classes for a coupled anisotropic neo-Hooke material with one fiber direction

\level 3
*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_MATELAST_COUPANISONEOHOOKE_HPP
#define FOUR_C_MATELAST_COUPANISONEOHOOKE_HPP

#include "baci_config.hpp"

#include "baci_mat_par_parameter.hpp"
#include "baci_matelast_summand.hpp"

BACI_NAMESPACE_OPEN

namespace MAT
{
  namespace ELASTIC
  {
    namespace PAR
    {
      /*!
       * @brief material parameters for anisochoric contribution of a neo-Hooke material with one
       * fiber direction
       *
       * <h3>Input line</h3>
       * MAT 1 CoupAnisoNeoHooke C 100 GAMMA 35.0 INIT 0 ADAPT_ANGLE 0
       */
      class CoupAnisoNeoHooke : public MAT::PAR::ParameterAniso
      {
       public:
        /// standard constructor
        CoupAnisoNeoHooke(const Teuchos::RCP<MAT::PAR::Material>& matdata);

        /// @name material parameters
        //@{

        /// fiber params
        double c_;
        /// angle between circumferential and fiber direction
        double gamma_;
        /// initalization modus;
        int init_;
        /// adapt angle during remodeling
        bool adapt_angle_;

        //@}

        /// Override this method and throw error, as the material should be created in within the
        /// Factory method of the elastic summand
        Teuchos::RCP<MAT::Material> CreateMaterial() override
        {
          dserror(
              "Cannot create a material from this method, as it should be created in "
              "MAT::ELASTIC::Summand::Factory.");
          return Teuchos::null;
        };
      };  // class CoupAnisoNeoHooke

    }  // namespace PAR

    /*!
     * @brief Coupled anisotropic exponential fiber function, implemented for one possible fiber
     * family as in [1] This is a hyperelastic, anisotropic material of the most simple kind.
     *
     * Strain energy function is given by
     * \f[
     *   \Psi = c \left(IV_{\boldsymbol C}-1\right).
     * \f]
     *
     * !!!Warning!!! The current implementation does not automatically respect the stress free
     * state of the reference configuration
     *
     * <h3>References</h3>
     * <ul>
     * <li> [1] GA Holzapfel, TC Gasser, A viscoelastic model for fiber-reinforced composites at
     * finite strains: ... 2000
     * </ul>
     */
    class CoupAnisoNeoHooke : public Summand
    {
     public:
      /// empty constructor
      //    CoupAnisoNeoHooke();

      /// constructor with given material parameters
      CoupAnisoNeoHooke(MAT::ELASTIC::PAR::CoupAnisoNeoHooke* params);

      ///@name Packing and Unpacking
      //@{

      void PackSummand(CORE::COMM::PackBuffer& data) const override;

      void UnpackSummand(
          const std::vector<char>& data, std::vector<char>::size_type& position) override;

      //@}

      /// @name Access material constants
      //@{

      /// material type
      INPAR::MAT::MaterialType MaterialType() const override
      {
        return INPAR::MAT::mes_coupanisoneohooke;
      }

      //@}

      /// Setup of summand
      void Setup(int numgp, INPUT::LineDefinition* linedef) override;

      /// Add anisotropic principal stresses
      void AddStressAnisoPrincipal(
          const CORE::LINALG::Matrix<6, 1>& rcg,  ///< right Cauchy Green Tensor
          CORE::LINALG::Matrix<6, 6>& cmat,       ///< material stiffness matrix
          CORE::LINALG::Matrix<6, 1>& stress,     ///< 2nd PK-stress
          Teuchos::ParameterList&
              params,  ///< additional parameters for computation of material properties
          int gp,      ///< Gauss point
          int eleGID   ///< element GID
          ) override;

      /// Set fiber directions
      void SetFiberVecs(const double newgamma,       ///< new angle
          const CORE::LINALG::Matrix<3, 3>& locsys,  ///< local coordinate system
          const CORE::LINALG::Matrix<3, 3>& defgrd   ///< deformation gradient
          ) override;

      /// Get fiber directions
      void GetFiberVecs(
          std::vector<CORE::LINALG::Matrix<3, 1>>& fibervecs  ///< vector of all fiber vectors
          ) override;

      /// Indicator for formulation
      void SpecifyFormulation(
          bool& isoprinc,     ///< global indicator for isotropic principal formulation
          bool& isomod,       ///< global indicator for isotropic splitted formulation
          bool& anisoprinc,   ///< global indicator for anisotropic principal formulation
          bool& anisomod,     ///< global indicator for anisotropic splitted formulation
          bool& viscogeneral  ///< global indicator, if one viscoelastic formulation is used
          ) override
      {
        anisoprinc = true;
        return;
      };

     private:
      /// my material parameters
      MAT::ELASTIC::PAR::CoupAnisoNeoHooke* params_;

      /// fiber direction
      CORE::LINALG::Matrix<3, 1> a_;
      /// structural tensors in voigt notation for anisotropy
      CORE::LINALG::Matrix<6, 1> A_;
    };

  }  // namespace ELASTIC
}  // namespace MAT

BACI_NAMESPACE_CLOSE

#endif