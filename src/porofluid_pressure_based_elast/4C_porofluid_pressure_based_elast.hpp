// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_POROFLUID_PRESSURE_BASED_ELAST_HPP
#define FOUR_C_POROFLUID_PRESSURE_BASED_ELAST_HPP


#include "4C_config.hpp"

#include "4C_adapter_algorithmbase.hpp"
#include "4C_linalg_utils_sparse_algebra_math.hpp"
#include "4C_linalg_vector.hpp"
#include "4C_linear_solver_method_linalg.hpp"
#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Adapter
{
  class PoroFluidMultiphaseWrapper;
  class Structure;
}  // namespace Adapter

namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace PoroPressureBased
{
  //! Base class of porofluid-elasticity algorithms
  class PorofluidElast : public Adapter::AlgorithmBase
  {
   public:
    PorofluidElast(MPI_Comm comm, const Teuchos::ParameterList& globaltimeparams);

    /// initialization
    virtual void init(const Teuchos::ParameterList& globaltimeparams,
        const Teuchos::ParameterList& algoparams, const Teuchos::ParameterList& structparams,
        const Teuchos::ParameterList& fluidparams, const std::string& struct_disname,
        const std::string& fluid_disname, bool isale, int nds_disp, int nds_vel,
        int nds_solidpressure, int ndsporofluid_scatra,
        const std::map<int, std::set<int>>* nearbyelepairs) = 0;

    virtual void post_init();

    /// read restart
    void read_restart(int restart) override;

    /// test results (if necessary)
    virtual void create_field_test();

    /// setup
    virtual void setup_system() = 0;

    /// prepare timeloop of coupled problem
    virtual void prepare_time_loop();

    /// timeloop of coupled problem
    virtual void timeloop();

    /// time step of coupled problem
    virtual void time_step() = 0;

    /// prepare time step of coupled problem
    void prepare_time_step() override;

    //! update fields after convergence
    virtual void update_and_output();

    /// dof map of vector of unknowns of structure field
    std::shared_ptr<const Core::LinAlg::Map> struct_dof_row_map() const;

    /// dof map of vector of unknowns of fluid field
    std::shared_ptr<const Core::LinAlg::Map> fluid_dof_row_map() const;

    /// dof map of vector of unknowns of artery field
    std::shared_ptr<const Core::LinAlg::Map> artery_dof_row_map() const;

    /// system matrix of coupled artery porofluid problem
    std::shared_ptr<Core::LinAlg::BlockSparseMatrixBase> artery_porofluid_sysmat() const;

    //! access to structural field
    const std::shared_ptr<Adapter::Structure>& structure_field() { return structure_; }

    //! access to fluid field
    const std::shared_ptr<Adapter::PoroFluidMultiphaseWrapper>& fluid_field() { return fluid_; }

    /// set structure solution on scatra field
    void set_struct_solution(std::shared_ptr<const Core::LinAlg::Vector<double>> disp,
        std::shared_ptr<const Core::LinAlg::Vector<double>> vel) const;

    /// set scatra solution on fluid field
    void set_scatra_solution(
        unsigned nds, std::shared_ptr<const Core::LinAlg::Vector<double>> scalars) const;

    //! setup solver (for monolithic only)
    virtual bool setup_solver() { return false; };

    /// unknown displacements at \f$t_{n+1}\f$
    std::shared_ptr<const Core::LinAlg::Vector<double>> struct_dispnp() const;

    /// unknown velocity at \f$t_{n+1}\f$
    std::shared_ptr<const Core::LinAlg::Vector<double>> struct_velnp() const;

    /// return fluid flux
    std::shared_ptr<const Core::LinAlg::MultiVector<double>> fluid_flux() const;

    /// return fluid solution variable
    std::shared_ptr<const Core::LinAlg::Vector<double>> fluid_phinp() const;

    /// return relaxed fluid solution variable (partitioned coupling will overwrite this method)
    virtual std::shared_ptr<const Core::LinAlg::Vector<double>> relaxed_fluid_phinp() const
    {
      return fluid_phinp();
    };

    /// set (relaxed) fluid solution on structure field (partitioned coupling will overwrite this
    /// method)
    virtual void set_relaxed_fluid_solution()
    {
      FOUR_C_THROW("set_relaxed_fluid_solution() only available for partitioned schemes!");
    };

    /// return fluid solution variable
    std::shared_ptr<const Core::LinAlg::Vector<double>> fluid_saturation() const;

    /// return fluid solution variable
    std::shared_ptr<const Core::LinAlg::Vector<double>> fluid_pressure() const;

    /// return fluid solution variable
    std::shared_ptr<const Core::LinAlg::Vector<double>> solid_pressure() const;

    //! unique map of all dofs that should be constrained with DBC
    virtual std::shared_ptr<const Core::LinAlg::Map> combined_dbc_map() const
    {
      FOUR_C_THROW("combined_dbc_map() only available for monolithic schemes!");
    };

    //! build the block null spaces
    virtual void build_block_null_spaces(std::shared_ptr<Core::LinAlg::Solver>& solver)
    {
      FOUR_C_THROW("build_block_null_spaces() only available for monolithic schemes!");
    };

    //! build the block null spaces
    virtual void build_artery_block_null_space(
        std::shared_ptr<Core::LinAlg::Solver>& solver, const int& arteryblocknum)
    {
      FOUR_C_THROW("build_artery_block_null_space() only available for monolithic schemes!");
    };

    //! evaluate all fields at x^n+1 with x^n+1 = x_n + stepinc
    virtual void evaluate(std::shared_ptr<const Core::LinAlg::Vector<double>> sx,
        std::shared_ptr<const Core::LinAlg::Vector<double>> fx, const bool firstcall)
    {
      FOUR_C_THROW("evaluate() only available for monolithic schemes!");
    };

    //! update all fields after convergence (add increment on displacements and fluid primary
    //! variables)
    virtual void update_fields_after_convergence(
        std::shared_ptr<const Core::LinAlg::Vector<double>>& sx,
        std::shared_ptr<const Core::LinAlg::Vector<double>>& fx)
    {
      FOUR_C_THROW("update_fields_after_convergence() only available for monolithic schemes!");
    };

    /// perform relaxaton (only for partitioned schemes)
    virtual void perform_relaxation(
        std::shared_ptr<const Core::LinAlg::Vector<double>> phi, const int itnum)
    {
      FOUR_C_THROW("PerformRelaxation() only available for partitioned schemes!");
    };

    //! get monolithic rhs vector
    virtual std::shared_ptr<const Core::LinAlg::Vector<double>> rhs() const
    {
      FOUR_C_THROW("RHS() only available for monolithic schemes!");
      return nullptr;
    };

    //! get extractor
    virtual std::shared_ptr<const Core::LinAlg::MultiMapExtractor> extractor() const
    {
      FOUR_C_THROW("Extractor() only available for monolithic schemes!");
    };

    //! get monolithic block system matrix
    virtual std::shared_ptr<Core::LinAlg::BlockSparseMatrixBase> block_system_matrix() const
    {
      FOUR_C_THROW("block_system_matrix() only available for monolithic schemes!");
    };

   private:
    /// set structure mesh displacement on fluid field
    void set_mesh_disp(std::shared_ptr<const Core::LinAlg::Vector<double>> disp) const;

    /// set structure velocity field on fluid field
    void set_velocity_fields(std::shared_ptr<const Core::LinAlg::Vector<double>> vel) const;

    /// underlying structure of the PoroMultiPhase problem
    std::shared_ptr<Adapter::Structure> structure_;

    /// underlying fluid problem of the PoroMultiPhase problem
    std::shared_ptr<Adapter::PoroFluidMultiphaseWrapper> fluid_;

   protected:
    /// a zero vector of full length of structure dofs
    std::shared_ptr<Core::LinAlg::Vector<double>> struct_zeros_;
    //! here the computation of the structure can be skipped, this is helpful if only fluid-scatra
    //! coupling should be calculated
    bool solve_structure_;

    /// coupling with 1D artery network
    const bool artery_coupl_;

    /// Print user output that structure field is disabled
    void print_structure_disabled_info() const;
  };

}  // namespace PoroPressureBased



FOUR_C_NAMESPACE_CLOSE

#endif
