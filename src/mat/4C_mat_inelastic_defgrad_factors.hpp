// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_MAT_INELASTIC_DEFGRAD_FACTORS_HPP
#define FOUR_C_MAT_INELASTIC_DEFGRAD_FACTORS_HPP

#include "4C_config.hpp"

#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_material_parameter_base.hpp"
#include "4C_utils_exceptions.hpp"
#include "4C_utils_parameter_list.fwd.hpp"

#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

namespace Discret::Utils
{
  class FunctionOfTime;
}

namespace Mat
{
  namespace PAR
  {
    enum class InelasticSource;

    /*----------------------------------------------------------------------*/
    /*! \class InelasticDeformationDirection
     *
     * Calculates and holds growth direction in matrix format for anisotropic growth
     */
    class InelasticDeformationDirection
    {
     public:
      /*!
       * @brief standard constructor
       * @param[in] growthdirection  direction of anisotropic growth
       */
      explicit InelasticDeformationDirection(std::vector<double> growthdirection);

      /// reference to matrix that determines growth direction
      const Core::LinAlg::Matrix<3, 3>& growth_dir_mat() const { return growth_dir_mat_; }

     private:
      /// matrix that determines growth direction
      Core::LinAlg::Matrix<3, 3> growth_dir_mat_;
    };

    /*----------------------------------------------------------------------
     *----------------------------------------------------------------------*/
    /*! \class InelasticDefgradNoGrowth
     *
     * This is a parameter class that is only needed to implement the pure virtual method
     * 'create_material()'.
     */
    class InelasticDefgradNoGrowth : public Core::Mat::PAR::Parameter
    {
     public:
      explicit InelasticDefgradNoGrowth(const Core::Mat::PAR::Parameter::Data& matdata);

      Teuchos::RCP<Core::Mat::Material> create_material() override { return Teuchos::null; }
    };


    /*----------------------------------------------------------------------
     *----------------------------------------------------------------------*/
    /*! \class InelasticDefgradScalar
     *
     * This is a parameter class holding parameters for evaluation of inelastic deformation (incl.
     * linearization) induced by a scalar.
     */
    class InelasticDefgradScalar : public Core::Mat::PAR::Parameter
    {
     public:
      explicit InelasticDefgradScalar(const Core::Mat::PAR::Parameter::Data& matdata);

      Teuchos::RCP<Core::Mat::Material> create_material() override { return Teuchos::null; }

      /// scalar that causes growth
      int scalar1() const { return scalar1_; }

      //! concentration, at which no growth occurs
      double scalar1_ref_conc() const { return scalar1_ref_conc_; }

     private:
      /// scalar that causes growth
      const int scalar1_;

      //! concentration, at which no growth occurs
      const double scalar1_ref_conc_;
    };

    /*----------------------------------------------------------------------
     *----------------------------------------------------------------------*/
    /*! \class InelasticDefgradTimeFunct
     *
     * This is a parameter class holding parameters for evaluation of inelastic deformation induced
     * by a given time-dependent function.
     */
    class InelasticDefgradTimeFunct : public Core::Mat::PAR::Parameter
    {
     public:
      explicit InelasticDefgradTimeFunct(const Core::Mat::PAR::Parameter::Data& matdata);

      Teuchos::RCP<Core::Mat::Material> create_material() override { return Teuchos::null; }

      /// function number that sets determinant of inelastic def. grad.
      int funct_num() const { return funct_num_; }

     private:
      /// function number that sets determinant of inelastic def. grad.
      const int funct_num_;
    };

    /*----------------------------------------------------------------------
     *----------------------------------------------------------------------*/
    /*! \class InelasticDefgradLinScalar
     *
     * This is a specialized parameter class that holds the growth factor for linear growth
     */
    class InelasticDefgradLinScalar : public InelasticDefgradScalar
    {
     public:
      /// standard constructor
      explicit InelasticDefgradLinScalar(const Core::Mat::PAR::Parameter::Data& matdata);

      //! molar factor that causes growth
      double scalar1_molar_growth_fac() { return scalar1_molar_growth_fac_; }

     private:
      //! molar factor that causes growth
      const double scalar1_molar_growth_fac_;
    };

    /*----------------------------------------------------------------------
     *----------------------------------------------------------------------*/
    /*! \class InelasticDefgradLinScalarAniso
     *
     * This is a specialized parameter class that can return the anisotropic growth direction
     * represented as a growth matrix
     */
    class InelasticDefgradLinScalarAniso : public InelasticDefgradLinScalar
    {
     public:
      /// standard constructor
      explicit InelasticDefgradLinScalarAniso(const Core::Mat::PAR::Parameter::Data& matdata);

      /// reference to matrix that determines growth direction
      const Core::LinAlg::Matrix<3, 3>& growth_dir_mat() { return growth_dir_->growth_dir_mat(); }

     private:
      /// calculation of direction of inelastic deformation
      Teuchos::RCP<const InelasticDeformationDirection> growth_dir_;
    };

    /*----------------------------------------------------------------------
     *----------------------------------------------------------------------*/
    /*! \class InelasticDefgradIntercalFrac
     *
     * This parameter class provides all electrochemical quantities that are needed to calculate the
     * intercalation fraction from a given species concentration.
     */
    class InelasticDefgradIntercalFrac : public InelasticDefgradScalar
    {
     public:
      explicit InelasticDefgradIntercalFrac(const Core::Mat::PAR::Parameter::Data& matdata);

      /// saturation concentration of material
      double cmax() const { return c_max_; }
      /// intercalation fraction at saturation concentration of material
      double chimax() const { return chi_max_; }

     private:
      /// saturation concentration of material
      double c_max_;
      /// intercalation fraction at saturation concentration of material
      double chi_max_;
    };

    /*----------------------------------------------------------------------
     *----------------------------------------------------------------------*/
    /*! \class InelasticDefgradPolyIntercalFrac
     *
     * This parameter class provides the value of the polynomial that models the growth evaluated in
     * the reference configuration.
     */
    class InelasticDefgradPolyIntercalFrac : public InelasticDefgradIntercalFrac
    {
     public:
      explicit InelasticDefgradPolyIntercalFrac(const Core::Mat::PAR::Parameter::Data& matdata);

      /// return value of polynomial at reference intercalation fraction
      double get_polynom_reference_value() const { return polynom_reference_value_; }

      // set value of polynomial at reference intercalation fraction
      void set_polynom_reference_value(double polynomReferenceValue)
      {
        polynom_reference_value_ = polynomReferenceValue;
      }

      //! polynomial coefficients that describe the growth law
      std::vector<double> poly_coeffs() const { return poly_coeffs_; }

      //! upper bound of polynomial
      double x_max() const { return x_max_; }

      //! lower bound of polynomial
      double x_min() const { return x_min_; }

     private:
      const std::vector<double> poly_coeffs_;

      /// value of polynomial at reference intercalation fraction
      double polynom_reference_value_;

      //! upper bound of polynomial
      const double x_max_;

      //! lower bound of polynomial
      const double x_min_;
    };

    /*----------------------------------------------------------------------
     *----------------------------------------------------------------------*/
    /*! \class InelasticDefgradPolyIntercalFracAniso
     *
     * This is a specialized parameter class that can return the anisotropic growth direction
     * represented as a growth matrix
     */
    class InelasticDefgradPolyIntercalFracAniso : public InelasticDefgradPolyIntercalFrac
    {
     public:
      /// standard constructor
      explicit InelasticDefgradPolyIntercalFracAniso(
          const Core::Mat::PAR::Parameter::Data& matdata);

      /// return reference to matrix that determines growth direction
      const Core::LinAlg::Matrix<3, 3>& growth_dir_mat() const
      {
        return growth_dir_->growth_dir_mat();
      };

     private:
      /// pointer to object, that calculates and holds direction of inelastic deformation
      Teuchos::RCP<InelasticDeformationDirection> growth_dir_;
    };

    /*----------------------------------------------------------------------
    ----------------------------------------------------------------------*/
    /*! \class InelasticDefgradLinTempIso

    Parameter class of InelasticDefgradLinTempIso.
    */
    class InelasticDefgradLinTempIso : public Core::Mat::PAR::Parameter
    {
     public:
      explicit InelasticDefgradLinTempIso(const Core::Mat::PAR::Parameter::Data& matdata);

      Teuchos::RCP<Core::Mat::Material> create_material() override { return Teuchos::null; };

      /// return temperature related growth factor
      double get_temp_growth_fac() const { return temp_growth_fac_; };

      /// return value of temperature that causes no growth
      double ref_temp() const { return ref_temp_; };

     private:
      /// value of temperature that causes no growth
      const double ref_temp_;

      /// growth factor
      const double temp_growth_fac_;
    };
  }  // namespace PAR

  /*----------------------------------------------------------------------*/
  /*! \class InelasticDefgradLinearShape
   *
   * This class provides the functionality to be used if the growth law obeys a linear relation
   */
  class InelasticDefgradLinearShape
  {
   public:
    /*!
     * @brief constructor with required parameters
     *
     * @param[in] growth_fac       linear growth factor (slope of linear function)
     * @param[in] reference_value  reference value
     */
    explicit InelasticDefgradLinearShape(double growth_fac, double reference_value);

    /*!
     * @brief evaluation of the linear growth law
     *
     * @param[in] value           value the linear relation shall be evaluated for
     * @return growth factor
     */
    double evaluate_linear_growth(double value) const;

    /// growth factor (needed for linearizations)
    double growth_fac() const { return growth_fac_; }

   private:
    /// growth factor
    const double growth_fac_;

    /// reference value
    const double reference_value_;
  };  // namespace Mat

  /*----------------------------------------------------------------------*/
  /*! \class InelasticDefgradPolynomialShape
   *
   * This class provides the functionality to be used if the growth law obeys a polynomial relation
   */
  class InelasticDefgradPolynomialShape
  {
   public:
    /*!
     * @brief  constructor with required parameters
     *
     * @param[in] poly_coeffs  coefficients describing the polynomial to be evaluated
     * @param[in] x_min        lower bound of validity of the polynomial
     * @param[in] x_max        upper bound of validity of the polynomial
     */
    explicit InelasticDefgradPolynomialShape(
        std::vector<double> poly_coeffs, double x_min, double x_max);

    /*!
     * @brief checks the bounds of validity of the polynomial and writes a warning to screen if
     * bounds are violated
     *
     * @param[in] x  value the polynomial is evaluated at
     */
    void check_polynomial_bounds(double x) const;

    /*!
     * @brief Evaluate the polynomial defined by #PolyCoeffs_ at the current position X
     *
     * @param[in] x  value the polynomial is evaluated at
     * @return value of the polynomial evaluated at x
     */
    double compute_polynomial(double x);

    /*!
     * @brief Evaluate the first derivative of the polynomial defined by #PolyCoeffs_ at the current
     * position x
     *
     * @param[in] x  value the polynomial is evaluated at
     * @return value the first derivative of the polynomial evaluated at x
     */
    double compute_polynomial_derivative(double x);

   private:
    /// coefficients of the polynomial to be evaluated
    const std::vector<double> poly_coeffs_;
    /// lower bound of validity of polynomial
    const double x_min_;
    /// upper bound of validity of polynomial
    const double x_max_;
  };

  /*----------------------------------------------------------------------*/
  /*! \class InelasticDefgradFactors

      Provides the interface called by the class "MultiplicativeSplitDefgrad_ElastHyper"
      and is needed to evaluate the inelastic deformation gradient and
      their derivatives w.r.t. the primary variables.

      In the material "MultiplicativeSplitDefgrad_ElastHyper" the deformation gradient is split
      multiplicatively in elastic and inelastic deformation gradients (F = F_{el} * F_{in}).
      The inelastic deformation gradient itself can be a product of different inelastic
      deformation gradients, i.e. F_{in} = F_{in,1} * F_{in,2} * ... * F_{in,n}.
      The derived classes below are needed to evaluate the inverse of the j-th inelastic
      deformation gradient F_{in,j}^{-1} and its derivatives w.r.t. the primary variables.
  */
  class InelasticDefgradFactors
  {
   public:
    /**
     * Virtual destructor.
     */
    virtual ~InelasticDefgradFactors() = default;

    /// construct material with specific material params
    explicit InelasticDefgradFactors(Core::Mat::PAR::Parameter* params);

    /*!
     * @brief create object by input parameter ID
     *
     * @param[in] matnum  material ID
     * @return pointer to material that is defined by material ID
     */
    static Teuchos::RCP<InelasticDefgradFactors> factory(int matnum);

    /// provide material type
    virtual Core::Materials::MaterialType material_type() const = 0;

    /*!
     * @brief evaluate the inelastic deformation gradient and its inverse
     *
     * @param[in] defgrad  Deformation gradient
     * @param[out] iFinM   Inverse inelastic deformation gradient
     */
    virtual void evaluate_inverse_inelastic_def_grad(
        const Core::LinAlg::Matrix<3, 3>* defgrad, Core::LinAlg::Matrix<3, 3>& iFinM) = 0;

    /*!
     * @brief evaluate additional terms for the elasticity tensor
     *
     * @param[in] defgrad  Deformation gradient
     * @param[in] iFinjM   Inverse inelastic deformation gradient of current inelastic contribution
     *                     as 3x3 matrix
     * @param[in] iCV      Inverse right Cauchy-Green tensor
     * @param[in] dSdiFinj Derivative of 2nd Piola Kirchhoff stresses w.r.t. the inverse inelastic
     *                     deformation gradient of current inelastic contribution
     * @param[in,out] cmatadd  Additional elasticity tensor
     */
    virtual void evaluate_additional_cmat(const Core::LinAlg::Matrix<3, 3>* defgrad,
        const Core::LinAlg::Matrix<3, 3>& iFinjM, const Core::LinAlg::Matrix<6, 1>& iCV,
        const Core::LinAlg::Matrix<6, 9>& dSdiFinj, Core::LinAlg::Matrix<6, 6>& cmatadd) = 0;

    /*!
     * @brief calculate the derivative of the inelastic deformation gradient
     *
     * @param[in] detjacobian  determinant of the deformation gradient
     * @param[out] dFindx      derivative of inelastic deformation gradient w.r.t. primary variable
     *                         of different field
     */
    virtual void evaluate_inelastic_def_grad_derivative(
        double detjacobian, Core::LinAlg::Matrix<9, 1>& dFindx) = 0;

    /*!
     * @brief evaluate off-diagonal stiffness matrix for monolithic systems to get the
     *        cross-linearizations
     *
     * @param[in] defgrad Deformation gradient
     * @param[in] iFinjM  Inverse inelastic deformation gradient of current inelastic contribution
     *                    as 3x3 matrix
     * @param[in] dSdiFinj  Derivative of 2nd Piola Kirchhoff stresses w.r.t. the inverse inelastic
     *                      deformation gradient of current inelastic contribution
     * @param[in,out] dstressdx Derivative of 2nd Piola Kirchhoff stresses w.r.t. primary variable
     *                          of different field
     */
    virtual void evaluate_od_stiff_mat(const Core::LinAlg::Matrix<3, 3>* defgrad,
        const Core::LinAlg::Matrix<3, 3>& iFinjM, const Core::LinAlg::Matrix<6, 9>& dSdiFinj,
        Core::LinAlg::Matrix<6, 1>& dstressdx) = 0;

    /*!
     * @brief pre-evaluation, intended to be used for stuff that has to be done only once per
     *        evaluate()
     *
     * @param[in] params  parameter list as handed in from the element
     * @param[in] gp      Gauss point
     */
    virtual void pre_evaluate(Teuchos::ParameterList& params, int gp) = 0;

    /*!
     * @brief set gauss point concentration to parameter class
     *
     * @param[in] concentration  gauss point concentration to be set to internal member of parameter
     *                           class
     *
     * @note This method is used by methods called from the contact algorithm. Since the gauss point
     * ids do not match anyways (volume vs. surface element gauss point ids) and the id is not
     * relevant since the method is only called for one gauss point anyways, we set it to a dummy
     * gauss point id of 0 here
     */
    virtual void set_concentration_gp(double concentration){};

    /// return material parameters
    virtual Core::Mat::PAR::Parameter* parameter() { return params_; }

    /// Get type of scalar, that leads to deformation
    virtual PAR::InelasticSource get_inelastic_source() = 0;

   private:
    /// material parameters
    Core::Mat::PAR::Parameter* params_;
  };

  /*--------------------------------------------------------------------*/
  /*! \class InelasticDefgradNoGrowth

   This class models materials in combination with the multiplicative split material that feature
   no volume changes, i.e. the inelastic deformation gradient is always the identity tensor and
   contributions to the linearizations therefore vanish.
   */
  class InelasticDefgradNoGrowth : public InelasticDefgradFactors
  {
   public:
    /*!
     * @brief construct material with required inputs
     *
     * @param[in] params           pointer to material specific parameters
     */
    explicit InelasticDefgradNoGrowth(Core::Mat::PAR::Parameter* params);

    void evaluate_additional_cmat(const Core::LinAlg::Matrix<3, 3>* defgrad,
        const Core::LinAlg::Matrix<3, 3>& iFinjM, const Core::LinAlg::Matrix<6, 1>& iCV,
        const Core::LinAlg::Matrix<6, 9>& dSdiFinj, Core::LinAlg::Matrix<6, 6>& cmatadd) override;

    void evaluate_inelastic_def_grad_derivative(
        double detjacobian, Core::LinAlg::Matrix<9, 1>& dFindx) override;

    void evaluate_inverse_inelastic_def_grad(
        const Core::LinAlg::Matrix<3, 3>* defgrad, Core::LinAlg::Matrix<3, 3>& iFinM) override;

    void evaluate_od_stiff_mat(const Core::LinAlg::Matrix<3, 3>* defgrad,
        const Core::LinAlg::Matrix<3, 3>& iFinjM, const Core::LinAlg::Matrix<6, 9>& dSdiFinj,
        Core::LinAlg::Matrix<6, 1>& dstressdx) override;

    PAR::InelasticSource get_inelastic_source() override;

    Core::Materials::MaterialType material_type() const override
    {
      return Core::Materials::mfi_no_growth;
    }

    void pre_evaluate(Teuchos::ParameterList& params, int gp) override;

   private:
    // identity tensor
    Core::LinAlg::Matrix<3, 3> identity_;
  };

  /*--------------------------------------------------------------------*/
  /*! \class InelasticDefgradTimeFunct

  This class models materials in combination with the multiplicative split material that feature
  isotropic volume changes based on a given time-dependent function for the determinant of the
  inelastic part.
  */
  class InelasticDefgradTimeFunct : public InelasticDefgradFactors
  {
   public:
    explicit InelasticDefgradTimeFunct(Core::Mat::PAR::Parameter* params);

    void evaluate_additional_cmat(const Core::LinAlg::Matrix<3, 3>* defgrad,
        const Core::LinAlg::Matrix<3, 3>& iFinjM, const Core::LinAlg::Matrix<6, 1>& iCV,
        const Core::LinAlg::Matrix<6, 9>& dSdiFinj, Core::LinAlg::Matrix<6, 6>& cmatadd) override{};

    void evaluate_inelastic_def_grad_derivative(
        double detjacobian, Core::LinAlg::Matrix<9, 1>& dFindx) override{};

    void evaluate_inverse_inelastic_def_grad(
        const Core::LinAlg::Matrix<3, 3>* defgrad, Core::LinAlg::Matrix<3, 3>& iFinM) override;

    void evaluate_od_stiff_mat(const Core::LinAlg::Matrix<3, 3>* defgrad,
        const Core::LinAlg::Matrix<3, 3>& iFinjM, const Core::LinAlg::Matrix<6, 9>& dSdiFinj,
        Core::LinAlg::Matrix<6, 1>& dstressdx) override{};

    PAR::InelasticSource get_inelastic_source() override;

    Core::Materials::MaterialType material_type() const override
    {
      return Core::Materials::mfi_time_funct;
    }

    Mat::PAR::InelasticDefgradTimeFunct* parameter() override
    {
      return dynamic_cast<Mat::PAR::InelasticDefgradTimeFunct*>(
          Mat::InelasticDefgradFactors::parameter());
    }

    void pre_evaluate(Teuchos::ParameterList& params, int gp) override;

   private:
    //! evaluated function value. Gets filled in pre_evaluate()
    double funct_value_;

    //! identity tensor
    Core::LinAlg::Matrix<3, 3> identity_;
  };

  class InelasticDefgradScalar : public InelasticDefgradFactors
  {
   public:
    /*!
     * @brief construct material with required inputs
     *
     * @param[in] params           pointer to material specific parameters
     */
    explicit InelasticDefgradScalar(Core::Mat::PAR::Parameter* params);

    void evaluate_additional_cmat(const Core::LinAlg::Matrix<3, 3>* defgrad,
        const Core::LinAlg::Matrix<3, 3>& iFinjM, const Core::LinAlg::Matrix<6, 1>& iCV,
        const Core::LinAlg::Matrix<6, 9>& dSdiFinj,
        Core::LinAlg::Matrix<6, 6>& cmatadd) override = 0;

    void evaluate_inverse_inelastic_def_grad(
        const Core::LinAlg::Matrix<3, 3>* defgrad, Core::LinAlg::Matrix<3, 3>& iFinM) override = 0;

    void evaluate_od_stiff_mat(const Core::LinAlg::Matrix<3, 3>* defgrad,
        const Core::LinAlg::Matrix<3, 3>& iFinjM, const Core::LinAlg::Matrix<6, 9>& dSdiFinj,
        Core::LinAlg::Matrix<6, 1>& dstressdx) override = 0;

    PAR::InelasticSource get_inelastic_source() override = 0;

    Core::Materials::MaterialType material_type() const override = 0;

    void pre_evaluate(Teuchos::ParameterList& params, int gp) override;

    void set_concentration_gp(double concentration) override;

    Mat::PAR::InelasticDefgradScalar* parameter() override
    {
      return dynamic_cast<Mat::PAR::InelasticDefgradScalar*>(
          Mat::InelasticDefgradFactors::parameter());
    }

   protected:
    //! Get vector of concentration at current Gauss point
    [[nodiscard]] const std::vector<double>& get_concentration_gp() const
    {
      FOUR_C_THROW_UNLESS(concentrations_ != Teuchos::null, "Concentrations are not set");
      return *concentrations_;
    };

   private:
    /// vector of concentations at the gauss points
    Teuchos::RCP<std::vector<double>> concentrations_{};
  };

  /*--------------------------------------------------------------------*/
  /*! \class InelasticDefgradPolyIntercalFrac

   This class evaluates polynomial and its first derivative w.r.t. intercalation fraction which is
   required in various routines of subclasses for isotropic and anisotropic case. This polynomial
   describes the growth of material with respect to intercalation fraction and it is prescribed by
   user in input file by defining it coefficients.
   */
  class InelasticDefgradPolyIntercalFrac : public InelasticDefgradScalar
  {
   public:
    /*!
     * @brief construct material with required inputs
     *
     * @param[in] params             pointer to material specific parameters
     * @param[in] polynomial_growth  pointer to object that evaluates the polynomial as prescribed
     *                               in the input file
     */
    explicit InelasticDefgradPolyIntercalFrac(Core::Mat::PAR::Parameter* params);

    /*!
     * @brief evaluate polynomial describing growth of material with regard to intercalation
     * fraction based on the current concentration
     *
     * @param[in] concentration current concentration
     * @param[in] detjacobian   determinant of the deformation gradient
     * @return value of polynomial describing the growth according to current intercalation fraction
     */
    double evaluate_polynomial(double concentration, double detjacobian);

    /*!
     * @brief evaluate the first derivative of the polynomial describing the growth
     *
     * @param[in] concentration current concentration
     * @param[in] detjacobian   determinant of the deformation gradient
     * @return first derivative of the polynomial describing the growth
     */
    double evaluate_polynomial_derivative(double concentration, double detjacobian);

    Core::Materials::MaterialType material_type() const override = 0;

    void evaluate_inverse_inelastic_def_grad(
        const Core::LinAlg::Matrix<3, 3>* defgrad, Core::LinAlg::Matrix<3, 3>& iFinM) override = 0;

    void evaluate_additional_cmat(const Core::LinAlg::Matrix<3, 3>* defgrad,
        const Core::LinAlg::Matrix<3, 3>& iFinjM, const Core::LinAlg::Matrix<6, 1>& iCV,
        const Core::LinAlg::Matrix<6, 9>& dSdiFinj,
        Core::LinAlg::Matrix<6, 6>& cmatadd) override = 0;

    void evaluate_inelastic_def_grad_derivative(
        double detjacobian, Core::LinAlg::Matrix<9, 1>& dFindx) override = 0;

    void evaluate_od_stiff_mat(const Core::LinAlg::Matrix<3, 3>* defgrad,
        const Core::LinAlg::Matrix<3, 3>& iFinjM, const Core::LinAlg::Matrix<6, 9>& dSdiFinj,
        Core::LinAlg::Matrix<6, 1>& dstressdx) override = 0;

    Mat::PAR::InelasticSource get_inelastic_source() override;

    Mat::PAR::InelasticDefgradPolyIntercalFrac* parameter() override
    {
      return dynamic_cast<Mat::PAR::InelasticDefgradPolyIntercalFrac*>(
          Mat::InelasticDefgradScalar::parameter());
    }

   private:
    /// pointer to class that evaluates the polynomial growth law
    Teuchos::RCP<InelasticDefgradPolynomialShape> polynomial_growth_;
  };

  /*----------------------------------------------------------------------*/
  /*! \class InelasticDefgradLinScalarIso
        This inelastic deformation gradient provides an isotropic growth law. Volumetric change due
        to this law is dependent on the current concentration \f$ c \f$ as follows :
      \f[
      \boldsymbol{F} _\text{in} = \left[1 + \text { scalar1_molar_growth_fac }
      \left(c \det \boldsymbol{F} - \text { Scalar1refconc } \right) \right] ^ { 1 / 3 }
      \boldsymbol{I}
      \f]
      */
  class InelasticDefgradLinScalarIso : public InelasticDefgradScalar
  {
   public:
    /*!
     * @brief construct material with required inputs
     *
     * @param[in] params          pointer to material specific parameters
     * @param[in] linear_growth   pointer to object that evaluates the linear relation as prescribed
     *                            in the input file
     */
    explicit InelasticDefgradLinScalarIso(Core::Mat::PAR::Parameter* params);

    Core::Materials::MaterialType material_type() const override
    {
      return Core::Materials::mfi_lin_scalar_iso;
    }

    void evaluate_inverse_inelastic_def_grad(
        const Core::LinAlg::Matrix<3, 3>* defgrad, Core::LinAlg::Matrix<3, 3>& iFinM) override;

    void evaluate_additional_cmat(const Core::LinAlg::Matrix<3, 3>* defgrad,
        const Core::LinAlg::Matrix<3, 3>& iFinjM, const Core::LinAlg::Matrix<6, 1>& iCV,
        const Core::LinAlg::Matrix<6, 9>& dSdiFinj, Core::LinAlg::Matrix<6, 6>& cmatadd) override;

    void evaluate_inelastic_def_grad_derivative(
        double detjacobian, Core::LinAlg::Matrix<9, 1>& dFindx) override;

    void evaluate_od_stiff_mat(const Core::LinAlg::Matrix<3, 3>* defgrad,
        const Core::LinAlg::Matrix<3, 3>& iFinjM, const Core::LinAlg::Matrix<6, 9>& dSdiFinj,
        Core::LinAlg::Matrix<6, 1>& dstressdc) override;

    Mat::PAR::InelasticSource get_inelastic_source() override;

    Mat::PAR::InelasticDefgradLinScalar* parameter() override
    {
      return dynamic_cast<Mat::PAR::InelasticDefgradLinScalar*>(
          Mat::InelasticDefgradScalar::parameter());
    }

   private:
    /// pointer to class that evaluates the linear growth law
    Teuchos::RCP<InelasticDefgradLinearShape> linear_growth_;
  };  // namespace Mat

  /*----------------------------------------------------------------------*/
  /*! \class InelasticDefgradLinScalarAniso

     This inelastic deformation gradient provides an anisotropic growth law.
     Volumetric change due to this law is dependent on the current concentration \f$ c \f$ as
     follows:
     \f[
     \mathbf{F}_\text{in} = \mathbf{I} + \left[ \text{scalar1_molar_growth_fac}
     \left( c \det\mathbf{F}  - \text{Scalar1refconc} \right) \right] \mathbf{G},
     \f]
     where \f$ \mathbf{G} \f$ (#growthdirmat_) is a matrix providing the information of the
     growth direction, that is constructed as follows:
     \f$ \mathbf{G} = \mathbf{g} \otimes \mathbf{g} \f$,
     where \f$ \mathbf{g} \f$ is the growth direction vector given in the input file.
     \f$ \mathbf{g} \f$ is normalized to length 1 before calculation of \f$ \mathbf{G} \f$.
     \f$ f(\chi) \f$ is defined by the user in the input file.
     */
  class InelasticDefgradLinScalarAniso : public InelasticDefgradScalar
  {
   public:
    /*!
     * @brief construct material with required inputs
     *
     * @param[in] params          pointer to material specific parameters
     * @param[in] linear_growth   pointer to object that evaluates the linear relation as prescribed
     *                            in the input file
     */
    explicit InelasticDefgradLinScalarAniso(Core::Mat::PAR::Parameter* params);

    Core::Materials::MaterialType material_type() const override
    {
      return Core::Materials::mfi_lin_scalar_aniso;
    }

    void evaluate_inverse_inelastic_def_grad(
        const Core::LinAlg::Matrix<3, 3>* defgrad, Core::LinAlg::Matrix<3, 3>& iFinM) override;

    void evaluate_additional_cmat(const Core::LinAlg::Matrix<3, 3>* defgrad,
        const Core::LinAlg::Matrix<3, 3>& iFinjM, const Core::LinAlg::Matrix<6, 1>& iCV,
        const Core::LinAlg::Matrix<6, 9>& dSdiFinj, Core::LinAlg::Matrix<6, 6>& cmatadd) override;

    void evaluate_inelastic_def_grad_derivative(
        double detjacobian, Core::LinAlg::Matrix<9, 1>& dFindx) override;

    void evaluate_od_stiff_mat(const Core::LinAlg::Matrix<3, 3>* defgrad,
        const Core::LinAlg::Matrix<3, 3>& iFinjM, const Core::LinAlg::Matrix<6, 9>& dSdiFinj,
        Core::LinAlg::Matrix<6, 1>& dstressdc) override;

    Mat::PAR::InelasticSource get_inelastic_source() override;

    Mat::PAR::InelasticDefgradLinScalarAniso* parameter() override
    {
      return dynamic_cast<Mat::PAR::InelasticDefgradLinScalarAniso*>(
          Mat::InelasticDefgradScalar::parameter());
    }

   private:
    /// store pointer to class that evaluates the linear growth law
    Teuchos::RCP<InelasticDefgradLinearShape> linear_growth_;
  };  // end of InelasticDefgradLinScalarAniso

  /*--------------------------------------------------------------------*/
  /*! \class InelasticDefgradPolyIntercalFracIso

   This inelastic deformation gradient provides an isotropic growth law.
   Volumetric change due to this law is non-linearly dependent on the intercalation fraction
   \f$ \chi \f$ as follows:
   \f[
   \boldsymbol{F}_\text{in} =
   \left[ \frac{f(\chi) + 1 }{f(\chi^0) + 1} \right]^{1/3} \boldsymbol{I},
   \f]
   where \f$ f(\chi) \f$ is defined by the user in the input file.
   */
  class InelasticDefgradPolyIntercalFracIso : public InelasticDefgradPolyIntercalFrac
  {
   public:
    /*!
     * @brief construct material with required inputs
     *
     * @param[in] params             pointer to material specific parameters
     * @param[in] polynomial_growth  pointer to object that evaluates the polynomial as prescribed
     *                               in the input file
     */
    explicit InelasticDefgradPolyIntercalFracIso(Core::Mat::PAR::Parameter* params);

    Core::Materials::MaterialType material_type() const override
    {
      return Core::Materials::mfi_poly_intercal_frac_iso;
    }

    void evaluate_inverse_inelastic_def_grad(
        const Core::LinAlg::Matrix<3, 3>* defgrad, Core::LinAlg::Matrix<3, 3>& iFinM) override;

    void evaluate_additional_cmat(const Core::LinAlg::Matrix<3, 3>* defgrad,
        const Core::LinAlg::Matrix<3, 3>& iFinjM, const Core::LinAlg::Matrix<6, 1>& iCV,
        const Core::LinAlg::Matrix<6, 9>& dSdiFinj, Core::LinAlg::Matrix<6, 6>& cmatadd) override;

    void evaluate_inelastic_def_grad_derivative(
        double detjacobian, Core::LinAlg::Matrix<9, 1>& dFindx) override;

    void evaluate_od_stiff_mat(const Core::LinAlg::Matrix<3, 3>* defgrad,
        const Core::LinAlg::Matrix<3, 3>& iFinjM, const Core::LinAlg::Matrix<6, 9>& dSdiFinj,
        Core::LinAlg::Matrix<6, 1>& dstressdc) override;

    Mat::PAR::InelasticDefgradPolyIntercalFrac* parameter() override
    {
      return dynamic_cast<Mat::PAR::InelasticDefgradPolyIntercalFrac*>(
          Mat::InelasticDefgradPolyIntercalFrac::parameter());
    }
  };

  /*----------------------------------------------------------------------*/
  /*! \class InelasticDefgradPolyIntercalFracAniso

   This inelastic deformation gradient provides an anisotropic growth law.
   Volumetric change due to this law is nonlinearly dependent on the intercalation fraction
   \f$ \chi \f$ as follows:
   \f[
   \boldsymbol{F}_\text{in} =
   \boldsymbol{I} + \left[ \frac{f(\chi) - f(\chi^0)}{f(\chi^0) + 1} \right] \boldsymbol{G},
   \f]
   where \f$ \boldsymbol{G} \f$ (#growthdirmat_) is a matrix providing the information of the growth
   direction, that is constructed as follows:
   \f$ \boldsymbol{G} = \boldsymbol{g} \otimes \boldsymbol{g} \f$, where \f$ \boldsymbol{g} \f$ is
   the growth direction vector given in the input file.
   \f$ \boldsymbol{g} \f$ is normalized to length 1 before calculation of \f$ \boldsymbol{G} \f$.
   \f$ f(\chi) \f$ is defined by the user in the input file.
   */
  class InelasticDefgradPolyIntercalFracAniso : public InelasticDefgradPolyIntercalFrac
  {
   public:
    /*!
     * @brief construct material with required inputs
     *
     * @param[in] params             pointer to material specific parameters
     * @param[in] polynomial_growth  pointer to object that evaluates the polynomial as prescribed
     *                               in the input file
     */
    explicit InelasticDefgradPolyIntercalFracAniso(Core::Mat::PAR::Parameter* params);

    Core::Materials::MaterialType material_type() const override
    {
      return Core::Materials::mfi_poly_intercal_frac_aniso;
    }

    void evaluate_inverse_inelastic_def_grad(
        const Core::LinAlg::Matrix<3, 3>* defgrad, Core::LinAlg::Matrix<3, 3>& iFinM) override;

    void evaluate_additional_cmat(const Core::LinAlg::Matrix<3, 3>* defgrad,
        const Core::LinAlg::Matrix<3, 3>& iFinjM, const Core::LinAlg::Matrix<6, 1>& iCV,
        const Core::LinAlg::Matrix<6, 9>& dSdiFinj, Core::LinAlg::Matrix<6, 6>& cmatadd) override;

    void evaluate_inelastic_def_grad_derivative(
        double detjacobian, Core::LinAlg::Matrix<9, 1>& dFindx) override;

    void evaluate_od_stiff_mat(const Core::LinAlg::Matrix<3, 3>* defgrad,
        const Core::LinAlg::Matrix<3, 3>& iFinjM, const Core::LinAlg::Matrix<6, 9>& dSdiFinj,
        Core::LinAlg::Matrix<6, 1>& dstressdc) override;

    Mat::PAR::InelasticDefgradPolyIntercalFracAniso* parameter() override
    {
      return dynamic_cast<Mat::PAR::InelasticDefgradPolyIntercalFracAniso*>(
          Mat::InelasticDefgradPolyIntercalFrac::parameter());
    }
  };

  /*----------------------------------------------------------------------*/
  /*! \class InelasticDefgradLinTempIso
   *Volumetric change due to this law is linearily dependent on the temperature
   \f$ T \f$ as follows:
   \f[
   \boldsymbol{F}_\text{in} = \boldsymbol{I} \left[ 1 + \beta \left( T - T_\text{ref} \right)
   \right]^\frac{1}{3}, \f]
   */
  class InelasticDefgradLinTempIso : public InelasticDefgradFactors
  {
   public:
    explicit InelasticDefgradLinTempIso(Core::Mat::PAR::Parameter* params);

    void evaluate_additional_cmat(const Core::LinAlg::Matrix<3, 3>* defgrad,
        const Core::LinAlg::Matrix<3, 3>& iFinjM, const Core::LinAlg::Matrix<6, 1>& iCV,
        const Core::LinAlg::Matrix<6, 9>& dSdiFinj, Core::LinAlg::Matrix<6, 6>& cmatadd) override;

    void evaluate_inelastic_def_grad_derivative(
        double detjacobian, Core::LinAlg::Matrix<9, 1>& dFindx) override;

    void evaluate_inverse_inelastic_def_grad(
        const Core::LinAlg::Matrix<3, 3>* defgrad, Core::LinAlg::Matrix<3, 3>& iFinM) override;

    void evaluate_od_stiff_mat(const Core::LinAlg::Matrix<3, 3>* defgrad,
        const Core::LinAlg::Matrix<3, 3>& iFinjM, const Core::LinAlg::Matrix<6, 9>& dSdiFinj,
        Core::LinAlg::Matrix<6, 1>& dstressdT) override;

    Mat::PAR::InelasticSource get_inelastic_source() override;

    Core::Materials::MaterialType material_type() const override
    {
      return Core::Materials::mfi_lin_temp_iso;
    };

    Mat::PAR::InelasticDefgradLinTempIso* parameter() override
    {
      return dynamic_cast<Mat::PAR::InelasticDefgradLinTempIso*>(
          Mat::InelasticDefgradFactors::parameter());
    }

    void pre_evaluate(Teuchos::ParameterList& params, int gp) override;

   private:
    /// temperature at the gauss point
    double temperature_ = 0.0;
  };
}  // namespace Mat

FOUR_C_NAMESPACE_CLOSE

#endif
