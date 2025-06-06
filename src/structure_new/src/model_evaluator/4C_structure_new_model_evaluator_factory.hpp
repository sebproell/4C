// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_STRUCTURE_NEW_MODEL_EVALUATOR_FACTORY_HPP
#define FOUR_C_STRUCTURE_NEW_MODEL_EVALUATOR_FACTORY_HPP

#include "4C_config.hpp"

#include "4C_inpar_structure.hpp"
#include "4C_structure_new_model_evaluator_manager.hpp"  // typedef

#include <memory>
#include <set>

FOUR_C_NAMESPACE_OPEN

namespace Solid
{
  namespace ModelEvaluator
  {
    class Generic;

    /*! Factory to build the desired model evaluator std::map
     *
     *  */
    class Factory
    {
     public:
      //! constructor
      Factory();

      //! destructor
      virtual ~Factory() = default;


      std::shared_ptr<Solid::ModelEvaluatorManager::Map> build_model_evaluators(
          const std::set<enum Inpar::Solid::ModelType>& modeltypes,
          const std::shared_ptr<Solid::ModelEvaluator::Generic>& coupling_model_ptr) const;

     private:
      //! return the proper type for the contact model evaluator
      std::shared_ptr<Solid::ModelEvaluator::Generic> build_contact_model_evaluator() const;

      //! return the proper type for the standard structural model evaluator
      std::shared_ptr<Solid::ModelEvaluator::Generic> build_structure_model_evaluator() const;

    };  // class Factory

    //! non-member function, which relates to the Solid::ModelEvaluator::Factory
    std::shared_ptr<Solid::ModelEvaluatorManager::Map> build_model_evaluators(
        const std::set<enum Inpar::Solid::ModelType>& modeltypes,
        const std::shared_ptr<Solid::ModelEvaluator::Generic>& coupling_model_ptr);

  }  // namespace ModelEvaluator
}  // namespace Solid

FOUR_C_NAMESPACE_CLOSE

#endif
