// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_SO3_HEX8FBAR_HPP
#define FOUR_C_SO3_HEX8FBAR_HPP

#include "4C_config.hpp"

#include "4C_so3_hex8.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace Discret
{
  namespace Elements
  {

    class SoHex8fbarType : public Core::Elements::ElementType
    {
     public:
      std::string name() const override { return "So_hex8fbarType"; }

      static SoHex8fbarType& instance();

      Core::Communication::ParObject* create(Core::Communication::UnpackBuffer& buffer) override;

      Teuchos::RCP<Core::Elements::Element> create(const std::string eletype,
          const std::string eledistype, const int id, const int owner) override;

      Teuchos::RCP<Core::Elements::Element> create(const int id, const int owner) override;

      int initialize(Core::FE::Discretization& dis) override;

      void nodal_block_information(
          Core::Elements::Element* dwele, int& numdf, int& dimns, int& nv, int& np) override;

      Core::LinAlg::SerialDenseMatrix compute_null_space(
          Core::Nodes::Node& node, const double* x0, const int numdof, const int dimnsp) override;

      void setup_element_definition(
          std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
          override;

     private:
      static SoHex8fbarType instance_;

      std::string get_element_type_string() const { return "SOLIDH8FBAR_DEPRECATED"; }
    };

    /*!
    \brief A C++ version of a 8-node hex solid element with F-Bar modification

    A structural 8-node hexahedral solid displacement element for large deformations
    and (near)-incompressibility. The F-bar technique is used to avoid volumetric locking.
    The volumetric part of the deformation gradient is only evaluated at the element
    center, which can be interpreted as a kind of selective reduced integration approach.
    The deviatoric part of the deformation gradient remains untouched.

    Refer also to the HiWi report of Stefanos Tsoukalas, 2010

    */
    class SoHex8fbar : public SoHex8
    {
     public:
      //! @name Friends
      friend class SoHex8fbarType;

      //@}
      //! @name Constructors and destructors and related methods

      /*!
      \brief Standard Constructor

      \param id : A unique global id
      \param owner : elements owner
      */
      SoHex8fbar(int id, int owner);

      /*!
      \brief Copy Constructor

      Makes a deep copy of a Element

      */
      SoHex8fbar(const SoHex8fbar& old);

      /*!
      \brief Deep copy this instance of Solid3 and return pointer to the copy

      The clone() method is used from the virtual base class Element in cases
      where the type of the derived class is unknown and a copy-ctor is needed

      */
      Core::Elements::Element* clone() const override;

      /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of this file.
      */
      int unique_par_object_id() const override
      {
        return SoHex8fbarType::instance().unique_par_object_id();
      }

      /*!
      \brief Pack this class so it can be communicated

      \ref pack and \ref unpack are used to communicate this element

      */
      void pack(Core::Communication::PackBuffer& data) const override;

      /*!
      \brief Unpack data from a char vector into this class

      \ref pack and \ref unpack are used to communicate this element

      */
      void unpack(Core::Communication::UnpackBuffer& buffer) override;


      //@}

      //! @name Access methods

      /*!
      \brief Print this element
      */
      void print(std::ostream& os) const override;

      Core::Elements::ElementType& element_type() const override
      {
        return SoHex8fbarType::instance();
      }

      //@}

      //! @name Input and Creation

      /*!
      \brief Read input for this element
      */
      bool read_element(const std::string& eletype, const std::string& distype,
          const Core::IO::InputParameterContainer& container) override;

      //@}

      //! @name Evaluation

      /*!
      \brief Evaluate an element

      Evaluate So_hex8fbar element stiffness, mass, internal forces, etc.

      If nullptr on input, the controlling method does not expect the element
      to fill these matrices or vectors.

      \return 0 if successful, negative otherwise
      */
      int evaluate(
          Teuchos::ParameterList&
              params,  ///< ParameterList for communication between control routine and elements
          Core::FE::Discretization& discretization,  ///< pointer to discretization for de-assembly
          std::vector<int>& lm,                      ///< location matrix for de-assembly
          Core::LinAlg::SerialDenseMatrix&
              elemat1,  ///< (stiffness-)matrix to be filled by element.
          Core::LinAlg::SerialDenseMatrix& elemat2,  ///< (mass-)matrix to be filled by element.
          Core::LinAlg::SerialDenseVector&
              elevec1,  ///< (internal force-)vector to be filled by element
          Core::LinAlg::SerialDenseVector& elevec2,  ///< vector to be filled by element
          Core::LinAlg::SerialDenseVector& elevec3   ///< vector to be filled by element
          ) override;

      /*!
      \brief Evaluate a Neumann boundary condition

      this method evaluates a surface Neumann condition on the solid3 element

      \param params (in/out)    : ParameterList for communication between control routine
                                  and elements
      \param discretization (in): A reference to the underlying discretization
      \param condition (in)     : The condition to be evaluated
      \param lm (in)            : location vector of this element
      \param elevec1 (out)      : vector to be filled by element. If nullptr on input,

      \return 0 if successful, negative otherwise
      */
      int evaluate_neumann(Teuchos::ParameterList& params, Core::FE::Discretization& discretization,
          Core::Conditions::Condition& condition, std::vector<int>& lm,
          Core::LinAlg::SerialDenseVector& elevec1,
          Core::LinAlg::SerialDenseMatrix* elemat1 = nullptr) override;

      //@}

     protected:
      //! don't want = operator
      SoHex8fbar& operator=(const SoHex8fbar& old);

      // compute Jacobian mapping wrt to deformed configuration
      void update_jacobian_mapping(
          const std::vector<double>& disp, Discret::Elements::PreStress& prestress) override;
      // compute defgrd in all gp for given disp
      void def_gradient(const std::vector<double>& disp, Core::LinAlg::SerialDenseMatrix& gpdefgrd,
          Discret::Elements::PreStress& prestress) override;

      //! Calculate nonlinear stiffness and mass matrix
      virtual void nlnstiffmass(std::vector<int>& lm,  ///< location matrix
          std::vector<double>& disp,                   ///< current displacements
          std::vector<double>* acc,                    // current accelerations
          std::vector<double>& residual,               ///< current residual displ
          Core::LinAlg::Matrix<NUMDOF_SOH8, NUMDOF_SOH8>*
              stiffmatrix,                                             ///< element stiffness matrix
          Core::LinAlg::Matrix<NUMDOF_SOH8, NUMDOF_SOH8>* massmatrix,  ///< element mass matrix
          Core::LinAlg::Matrix<NUMDOF_SOH8, 1>* force,       ///< element internal force vector
          Core::LinAlg::Matrix<NUMDOF_SOH8, 1>* forceinert,  // element inertial force vector
          Core::LinAlg::Matrix<NUMGPT_SOH8, Mat::NUM_STRESS_3D>* elestress,  ///< stresses at GP
          Core::LinAlg::Matrix<NUMGPT_SOH8, Mat::NUM_STRESS_3D>* elestrain,  ///< strains at GP
          Core::LinAlg::Matrix<NUMGPT_SOH8, Mat::NUM_STRESS_3D>*
              eleplstrain,                             ///< plastic strains at GP
          Teuchos::ParameterList& params,              ///< algorithmic parameters e.g. time
          const Inpar::Solid::StressType iostress,     ///< stress output option
          const Inpar::Solid::StrainType iostrain,     ///< strain output option
          const Inpar::Solid::StrainType ioplstrain);  ///< strain output option

      //! Update history variables at the end of time step (inelastic deformation) (braeu 07/16)
      void update_element(std::vector<double>& disp,  // current displacements
          Teuchos::ParameterList& params,             // algorithmic parameters e.g. time
          Core::Mat::Material& mat);                  // material

      //! init the inverse of the jacobian and its determinant in the material configuration
      void init_jacobian_mapping() override;

     private:
      std::string get_element_type_string() const { return "SOLIDH8FBAR_DEPRECATED"; }
    };  // class So_hex8fbar

    //=======================================================================
    //=======================================================================
    //=======================================================================
    //=======================================================================

  }  // namespace Elements
}  // namespace Discret

FOUR_C_NAMESPACE_CLOSE

#endif
