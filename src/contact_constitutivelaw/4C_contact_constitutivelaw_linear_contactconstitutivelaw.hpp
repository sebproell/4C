/*----------------------------------------------------------------------*/
/*! \file
\brief Implements a default linear penalty function as contact constitutive law

\level 2

*/
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_CONTACT_CONSTITUTIVELAW_LINEAR_CONTACTCONSTITUTIVELAW_HPP
#define FOUR_C_CONTACT_CONSTITUTIVELAW_LINEAR_CONTACTCONSTITUTIVELAW_HPP


#include "4C_config.hpp"

#include "4C_contact_constitutivelaw_contactconstitutivelaw.hpp"
#include "4C_contact_constitutivelaw_contactconstitutivelaw_parameter.hpp"

FOUR_C_NAMESPACE_OPEN


namespace CONTACT
{
  namespace CONSTITUTIVELAW
  {
    /*----------------------------------------------------------------------*/
    /** \brief constitutive law parameters for a linear contact law \f$ Ax+B \f$ relating the gap to
     * the contact pressure
     *
     */
    class LinearConstitutiveLawParams : public Parameter
    {
     public:
      /** \brief standard constructor
       * \param[in] container containing the law parameter from the input file
       */
      LinearConstitutiveLawParams(
          const Teuchos::RCP<const CONTACT::CONSTITUTIVELAW::Container> container);


      /// create constitutive law instance of matching type with my parameters
      Teuchos::RCP<CONTACT::CONSTITUTIVELAW::ConstitutiveLaw> CreateConstitutiveLaw() override;

      /// @name get-functions for the Constitutive Law parameters of a broken rational function
      //@{
      /// Get the slope
      double GetA() const { return a_; };
      /// Get the y intercept
      double GetB() const { return b_; };
      //@}

     private:
      /// @name Constitutive Law parameters of a linear function
      //@{
      /// slope
      double a_;
      /// y intercept
      const double b_;
      //@}
    };  // class

    /*----------------------------------------------------------------------*/
    /** \brief implements a linear contact constitutive law \f$ Ax+B \f$ relating the gap to the
     * contact pressure
     */
    class LinearConstitutiveLaw : public ConstitutiveLaw
    {
     public:
      /// construct the constitutive law object given a set of parameters
      explicit LinearConstitutiveLaw(CONTACT::CONSTITUTIVELAW::LinearConstitutiveLawParams* params);

      //! @name Access methods

      /// contact constitutive law type
      INPAR::CONTACT::ConstitutiveLawType GetConstitutiveLawType() const override
      {
        return INPAR::CONTACT::ConstitutiveLawType::colaw_linear;
      }

      /// Get slope of linear polynomial
      double GetA() { return params_->GetA(); }
      /// Get y intercept of linear polynomial
      double GetB() { return params_->GetB(); }

      /// Return quick accessible contact constitutive law parameter data
      CONTACT::CONSTITUTIVELAW::Parameter* Parameter() const override { return params_; }

      //! @name Evaluation methods
      //@{
      /// evaluate the constitutive law
      double Evaluate(double gap, CONTACT::Node* cnode) override;

      /// Evaluate derivative of the constitutive law
      double EvaluateDeriv(double gap, CONTACT::Node* cnode) override;
      //@}

     private:
      /// my constitutive law parameters
      CONTACT::CONSTITUTIVELAW::LinearConstitutiveLawParams* params_;
    };
  }  // namespace CONSTITUTIVELAW
}  // namespace CONTACT

FOUR_C_NAMESPACE_CLOSE

#endif