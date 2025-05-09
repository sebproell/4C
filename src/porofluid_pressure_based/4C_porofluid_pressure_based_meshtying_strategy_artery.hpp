// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_POROFLUID_PRESSURE_BASED_MESHTYING_STRATEGY_ARTERY_HPP
#define FOUR_C_POROFLUID_PRESSURE_BASED_MESHTYING_STRATEGY_ARTERY_HPP

#include "4C_config.hpp"

#include "4C_linear_solver_method_linalg.hpp"
#include "4C_porofluid_pressure_based_algorithm.hpp"

FOUR_C_NAMESPACE_OPEN

// forward declaration
namespace PoroPressureBased
{
  class PoroMultiPhaseScaTraArtCouplBase;
}

namespace PoroPressureBased
{
  class MeshtyingArtery
  {
   public:
    //! constructor
    explicit MeshtyingArtery(PorofluidAlgorithm* porofluid_algorithm,
        const Teuchos::ParameterList& probparams, const Teuchos::ParameterList& poroparams);


    //! prepare time loop
    void prepare_time_loop();

    //! prepare time step
    void prepare_time_step();

    //! update
    void update();

    //! output
    void output();

    //! Initialize the linear solver
    void initialize_linear_solver(std::shared_ptr<Core::LinAlg::Solver> solver);

    //! solve linear system of equations
    void linear_solve(std::shared_ptr<Core::LinAlg::Solver> solver,
        std::shared_ptr<Core::LinAlg::SparseOperator> sysmat,
        std::shared_ptr<Core::LinAlg::Vector<double>> increment,
        std::shared_ptr<Core::LinAlg::Vector<double>> residual,
        Core::LinAlg::SolverParams& solver_params);

    //! calculate norms for convergence checks
    void calculate_norms(std::vector<double>& preresnorm, std::vector<double>& incprenorm,
        std::vector<double>& prenorm,
        const std::shared_ptr<const Core::LinAlg::Vector<double>> increment);

    //! create the field test
    void create_field_test();

    //! restart
    void read_restart(const int step);

    //! evaluate mesh tying
    void evaluate();

    //! extract increments and update mesh tying
    std::shared_ptr<const Core::LinAlg::Vector<double>> extract_and_update_iter(
        const std::shared_ptr<const Core::LinAlg::Vector<double>> inc);

    // return arterial network time integrator
    std::shared_ptr<Adapter::ArtNet> art_net_tim_int() { return artnettimint_; }

    //! access dof row map
    std::shared_ptr<const Core::LinAlg::Map> artery_dof_row_map() const;

    //! right-hand side alias the dynamic force residual for coupled system
    std::shared_ptr<const Core::LinAlg::Vector<double>> artery_porofluid_rhs() const;

    //! access to block system matrix of artery poro problem
    std::shared_ptr<Core::LinAlg::BlockSparseMatrixBase> artery_porofluid_sysmat() const;

    //! get global (combined) increment of coupled problem
    std::shared_ptr<const Core::LinAlg::Vector<double>> combined_increment(
        const std::shared_ptr<const Core::LinAlg::Vector<double>> inc) const;

    //! check if initial fields on coupled DOFs are equal
    void check_initial_fields(std::shared_ptr<const Core::LinAlg::Vector<double>> vec_cont) const;

    //! set the element pairs that are close as found by search algorithm
    void set_nearby_ele_pairs(const std::map<int, std::set<int>>* nearbyelepairs);

    //! setup the strategy
    void setup();

    //! apply the mesh movement
    void apply_mesh_movement() const;

    //! return blood vessel volume fraction
    std::shared_ptr<const Core::LinAlg::Vector<double>> blood_vessel_volume_fraction();

   private:
    //! porofluid algorithm
    PorofluidAlgorithm* porofluid_algorithm_;

    //! parameter list of global control problem
    const Teuchos::ParameterList& params_;

    //! parameter list of poro fluid multiphase problem
    const Teuchos::ParameterList& poroparams_;

    //! vector norm for residuals
    VectorNorm vectornormfres_;

    //! vector norm for increments
    VectorNorm vectornorminc_;

    //! artery time integration
    std::shared_ptr<Adapter::ArtNet> artnettimint_;

    //! artery discretization
    std::shared_ptr<Core::FE::Discretization> arterydis_;

    //! the mesh tying object
    std::shared_ptr<PoroPressureBased::PoroMultiPhaseScaTraArtCouplBase> arttoporofluidcoupling_;

    //! block systemmatrix
    std::shared_ptr<Core::LinAlg::BlockSparseMatrixBase> comb_systemmatrix_;

    //! global rhs
    std::shared_ptr<Core::LinAlg::Vector<double>> rhs_;

    //! global increment
    std::shared_ptr<Core::LinAlg::Vector<double>> comb_increment_;

    //! global solution at time n+1
    std::shared_ptr<Core::LinAlg::Vector<double>> comb_phinp_;
  };

}  // namespace PoroPressureBased



FOUR_C_NAMESPACE_CLOSE

#endif
