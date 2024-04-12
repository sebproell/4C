/*----------------------------------------------------------------------*/
/*! \file
\brief Implementation of a viscous material contribution, calculated according to an SLS-model

\level 1
*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_MATELAST_VISCO_GENMAX_HPP
#define FOUR_C_MATELAST_VISCO_GENMAX_HPP

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
       * @brief material parameters for viscous contribution according the SLS-model
       *
       * <h3>Input line</h3>
       * MAT 1 VISCO_GenMax TAU 0.1 BETA 1 SOLVE OST
       * MAT 1 VISCO_GenMax TAU 0.1 BETA 1 SOLVE CONVOL
       */
      class GenMax : public MAT::PAR::Parameter
      {
       public:
        /// standard constructor
        GenMax(const Teuchos::RCP<MAT::PAR::Material>& matdata);

        /// @name material parameters
        //@{

        /// material parameters
        double tau_;
        double beta_;
        std::string solve_;

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
      };  // class GenMax

    }  // namespace PAR

    /*!
     * @brief Material Viscogenmax
     *
     * This material offers a viscous and hyperelastic part. The model consists
     * of one spring in parallel to one sequential branch of a spring and a dashpot.
     *
     * The hyperelasic part is possibly any hyperelastic law of the Hyperelastic
     * toolbox. The materials of the hyperelastic toolbox are composed
     * of (Helmholtz free energy density) potentials.  Effectively, we want
     * \f[
     *  \Psi(\boldsymbol{C}) = \sum_i \Psi_i(\boldsymbol{C})
     * \f]
     *
     * in which the individual \f$\Psi_i\f$ is implemented as #MAT::ELASTIC::Summand.
     *
     *
     * The stress tensor of the viscous part is calculated from the evolution equation
     * \f[
     *  \dot{\boldsymbol{Q}} = 1 / \tau \boldsymbol{Q} = \beta \dot{\boldsymbol{S}}
     * \f]
     * with elastic stress tensor S. S contains isochoric and volumetric contribution,
     * so Q is also isochoric and volumetric.
     * The viscous elasicity tensor is calculated from \f$\boldsymbol{Q}\f$ with
     * derivation to \f$\boldsymbol{C}\f$.
     * The viscous effect can be applied to any part of the SEF (isotropic coupled,
     * isotropic isochoric, isotropic volumetric, anisotropic).
     *
     * <h3>References</h3>
     * <ul>
     * <li> [1] GA Holzapfel, "Nonlinear solid mechanics", Wiley, 2000.
     * <li> [2] Bul-Brunon et.al., "Numerical identification method for the non-linear
     *         viscoelastic compressible behavior of soft tissue using uniaxial tensile
     *         tests and image registration - Application to rat lung parenchyma, 2013
     * </ul>
     */

    class GenMax : public Summand
    {
     public:
      /// constructor with given material parameters
      GenMax(MAT::ELASTIC::PAR::GenMax* params);

      /// @name Access material constants
      //@{

      /// material type
      INPAR::MAT::MaterialType MaterialType() const override { return INPAR::MAT::mes_genmax; }

      //@}

      /// Read material parameters
      void ReadMaterialParametersVisco(double& tau,  ///< relaxation parameter tau
          double& beta,                              ///< emphasis of viscous to elastic part
          double& alpha,  ///< fractional order derivative (just for visoc_fract)
          std::string&
              solve  //!< variant of the solution of the evolution integral (just for genmax)
          ) override;

      /// Indicator for formulation
      void SpecifyFormulation(
          bool& isoprinc,     ///< global indicator for isotropic principal formulation
          bool& isomod,       ///< global indicator for isotropic splitted formulation
          bool& anisoprinc,   ///< global indicator for anisotropic principal formulation
          bool& anisomod,     ///< global indicator for anisotropic splitted formulation
          bool& viscogeneral  ///< general indicator, if one viscoelastic formulation is used
          ) override
      {
        viscogeneral = true;
        return;
      };

      /// Indicator for the chosen viscoelastic formulations
      void SpecifyViscoFormulation(
          bool& isovisco,     ///< global indicator for isotropic, splitted and viscous formulation
          bool& viscogenmax,  ///< global indicator for viscous contribution according the SLS-Model
          bool& viscogeneralizedgenmax,  ///< global indicator for viscoelastic contribution
                                         ///< according to the generalized Maxwell Model
          bool& viscofract  ///< global indicator for viscous contribution according the FSLS-Model
          ) override
      {
        viscogenmax = true;
        return;
      };


     private:
      /// my material parameters
      MAT::ELASTIC::PAR::GenMax* params_;
    };


  }  // namespace ELASTIC
}  // namespace MAT

BACI_NAMESPACE_CLOSE

#endif