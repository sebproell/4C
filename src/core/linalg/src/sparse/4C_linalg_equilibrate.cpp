// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_linalg_equilibrate.hpp"

#include "4C_linalg_utils_sparse_algebra_create.hpp"

#include <utility>

FOUR_C_NAMESPACE_OPEN

/*-------------------------------------------------------------------------*
 *-------------------------------------------------------------------------*/
Core::LinAlg::Equilibration::Equilibration(std::shared_ptr<const Epetra_Map> dofrowmap)
    : invcolsums_(Core::LinAlg::create_vector(*dofrowmap, false)),
      invrowsums_(Core::LinAlg::create_vector(*dofrowmap, false))
{
}

Core::LinAlg::EquilibrationUniversal::EquilibrationUniversal(
    EquilibrationMethod method, std::shared_ptr<const Epetra_Map> dofrowmap)
    : Equilibration(dofrowmap), method_(method)
{
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Core::LinAlg::EquilibrationSparse::EquilibrationSparse(
    EquilibrationMethod method, std::shared_ptr<const Epetra_Map> dofrowmap)
    : EquilibrationUniversal(method, dofrowmap)
{
  if (method == EquilibrationMethod::symmetry)
    FOUR_C_THROW("symmetric equilibration not implemented for sparse matrices");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Core::LinAlg::EquilibrationBlock::EquilibrationBlock(
    EquilibrationMethod method, std::shared_ptr<const Epetra_Map> dofrowmap)
    : EquilibrationUniversal(method, dofrowmap)
{
  if (method == EquilibrationMethod::symmetry)
    FOUR_C_THROW("symmetric equilibration not implemented for block matrices");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Core::LinAlg::EquilibrationBlockSpecific::EquilibrationBlockSpecific(
    const std::vector<EquilibrationMethod>& method, std::shared_ptr<const Epetra_Map> dofrowmap)
    : Equilibration(dofrowmap), method_blocks_(method)
{
  for (const auto& method_block : method_blocks_)
  {
    if (method_block == EquilibrationMethod::columns_full or
        method_block == EquilibrationMethod::rows_full or
        method_block == EquilibrationMethod::rowsandcolumns_full)
      FOUR_C_THROW("full matrix equilibration not reasonable for block based equilibration");
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Core::LinAlg::Equilibration::compute_inv_row_sums(const Core::LinAlg::SparseMatrix& matrix,
    Core::LinAlg::Vector<double>& invrowsums, const EquilibrationMethod method) const
{
  // compute inverse row sums of matrix
  if (matrix.epetra_matrix()->InvRowSums(invrowsums.get_ref_of_epetra_vector()))
    FOUR_C_THROW("Inverse row sums of matrix could not be successfully computed!");

  // take square root of inverse row sums if matrix is scaled from left and right
  if (method == EquilibrationMethod::rowsandcolumns_full or
      method == EquilibrationMethod::rowsandcolumns_maindiag)
    for (int i = 0; i < invrowsums.local_length(); ++i)
      (invrowsums)[i] = std::sqrt((invrowsums)[i]);
}


/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
void Core::LinAlg::Equilibration::compute_inv_col_sums(const Core::LinAlg::SparseMatrix& matrix,
    Core::LinAlg::Vector<double>& invcolsums, const EquilibrationMethod method) const
{
  // compute inverse column sums of matrix
  if (matrix.epetra_matrix()->InvColSums(invcolsums.get_ref_of_epetra_vector()))
    FOUR_C_THROW("Inverse column sums of matrix could not be successfully computed!");

  // take square root of inverse column sums if matrix is scaled from left and right
  if (method == EquilibrationMethod::rowsandcolumns_full or
      method == EquilibrationMethod::rowsandcolumns_maindiag)
    for (int i = 0; i < invcolsums.local_length(); ++i)
      (invcolsums)[i] = std::sqrt((invcolsums)[i]);
}

/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
void Core::LinAlg::Equilibration::compute_inv_symmetry(
    const Core::LinAlg::SparseMatrix& matrix, Core::LinAlg::Vector<double>& invsymmetry) const
{
  std::shared_ptr<Core::LinAlg::Vector<double>> diag =
      Core::LinAlg::create_vector(matrix.range_map(), true);
  matrix.extract_diagonal_copy(*diag);

  for (int my_row = 0; my_row < diag->get_map().NumMyElements(); ++my_row)
  {
    (invsymmetry)[my_row] = 1.0 / std::sqrt((*diag)[my_row]);
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Core::LinAlg::Equilibration::equilibrate_matrix_rows(
    Core::LinAlg::SparseMatrix& matrix, const Core::LinAlg::Vector<double>& invrowsums) const
{
  if (matrix.left_scale(invrowsums)) FOUR_C_THROW("Row equilibration of matrix failed!");
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Core::LinAlg::Equilibration::equilibrate_matrix_columns(
    Core::LinAlg::SparseMatrix& matrix, const Core::LinAlg::Vector<double>& invcolsums) const
{
  if (matrix.right_scale(invcolsums)) FOUR_C_THROW("Column equilibration of matrix failed!");
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Core::LinAlg::Equilibration::equilibrate_system(
    std::shared_ptr<Core::LinAlg::SparseOperator> systemmatrix,
    std::shared_ptr<Core::LinAlg::Vector<double>> residual,
    std::shared_ptr<const Core::LinAlg::MultiMapExtractor> blockmaps) const
{
  // Equilibrate the matrix given the chosen method and matrix type
  equilibrate_matrix(systemmatrix, blockmaps);

  // Equilibrate the RHS
  equilibrate_rhs(residual);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Core::LinAlg::EquilibrationUniversal::unequilibrate_increment(
    std::shared_ptr<Core::LinAlg::Vector<double>> increment) const
{
  // unequilibrate global increment vector if necessary
  if (method() == EquilibrationMethod::columns_full or
      method() == EquilibrationMethod::columns_maindiag or
      method() == EquilibrationMethod::rowsandcolumns_full or
      method() == EquilibrationMethod::rowsandcolumns_maindiag)
  {
    if (increment->multiply(1.0, *invcolsums_, *increment, 0.0))
      FOUR_C_THROW("Unequilibration of global increment vector failed!");
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Core::LinAlg::EquilibrationUniversal::equilibrate_rhs(
    std::shared_ptr<Core::LinAlg::Vector<double>> residual) const
{
  // perform equilibration of global residual vector
  if (method() == EquilibrationMethod::rows_full or
      method() == EquilibrationMethod::rows_maindiag or
      method() == EquilibrationMethod::rowsandcolumns_full or
      method() == EquilibrationMethod::rowsandcolumns_maindiag)
    if (residual->multiply(1.0, *invrowsums_, *residual, 0.0))
      FOUR_C_THROW("Equilibration of global residual vector failed!");
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Core::LinAlg::EquilibrationBlockSpecific::unequilibrate_increment(
    std::shared_ptr<Core::LinAlg::Vector<double>> increment) const
{
  if (increment->multiply(1.0, *invcolsums_, *increment, 0.0))
    FOUR_C_THROW("Unequilibration of global increment vector failed!");
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Core::LinAlg::EquilibrationBlockSpecific::equilibrate_rhs(
    std::shared_ptr<Core::LinAlg::Vector<double>> residual) const
{
  if (residual->multiply(1.0, *invrowsums_, *residual, 0.0))
    FOUR_C_THROW("Equilibration of global residual vector failed!");
}

/*-------------------------------------------------------------------------*
 *-------------------------------------------------------------------------*/
void Core::LinAlg::EquilibrationSparse::equilibrate_matrix(
    std::shared_ptr<Core::LinAlg::SparseOperator> systemmatrix,
    std::shared_ptr<const Core::LinAlg::MultiMapExtractor> blockmaps) const
{
  equilibrate_matrix(systemmatrix);
}

/*-------------------------------------------------------------------------*
 *-------------------------------------------------------------------------*/
void Core::LinAlg::EquilibrationSparse::equilibrate_matrix(
    std::shared_ptr<Core::LinAlg::SparseOperator> systemmatrix) const
{
  std::shared_ptr<Core::LinAlg::SparseMatrix> sparsematrix =
      Core::LinAlg::cast_to_sparse_matrix_and_check_success(systemmatrix);

  // perform row equilibration
  if (method() == EquilibrationMethod::rows_full or
      method() == EquilibrationMethod::rows_maindiag or
      method() == EquilibrationMethod::rowsandcolumns_full or
      method() == EquilibrationMethod::rowsandcolumns_maindiag)
  {
    // compute inverse row sums of global system matrix
    compute_inv_row_sums(*sparsematrix, *invrowsums_, method());

    // perform row equilibration of global system matrix
    equilibrate_matrix_rows(*sparsematrix, *invrowsums_);
  }

  // perform column equilibration
  if (method() == EquilibrationMethod::columns_full or
      method() == EquilibrationMethod::columns_maindiag or
      method() == EquilibrationMethod::rowsandcolumns_full or
      method() == EquilibrationMethod::rowsandcolumns_maindiag)
  {
    // compute inverse column sums of global system matrix
    compute_inv_col_sums(*sparsematrix, *invcolsums_, method());

    // perform column equilibration of global system matrix
    equilibrate_matrix_columns(*sparsematrix, *invcolsums_);
  }
}

/*-------------------------------------------------------------------------*
 *-------------------------------------------------------------------------*/
void Core::LinAlg::EquilibrationBlock::equilibrate_matrix(
    std::shared_ptr<Core::LinAlg::SparseOperator> systemmatrix,
    std::shared_ptr<const Core::LinAlg::MultiMapExtractor> blockmaps) const
{
  std::shared_ptr<Core::LinAlg::BlockSparseMatrixBase> blocksparsematrix =
      Core::LinAlg::cast_to_block_sparse_matrix_base_and_check_success(systemmatrix);

  // perform row equilibration
  if (method() == EquilibrationMethod::rows_full or
      method() == EquilibrationMethod::rows_maindiag or
      method() == EquilibrationMethod::rowsandcolumns_full or
      method() == EquilibrationMethod::rowsandcolumns_maindiag)
  {
    for (int i = 0; i < blocksparsematrix->rows(); ++i)
    {
      // initialize vector for inverse row sums
      Vector<double> invrowsums(blocksparsematrix->matrix(i, i).row_map());

      // compute inverse row sums of current main diagonal matrix block
      if (method() == EquilibrationMethod::rows_maindiag or
          method() == EquilibrationMethod::rowsandcolumns_maindiag)
      {
        compute_inv_row_sums(blocksparsematrix->matrix(i, i), invrowsums, method());
      }
      // compute inverse row sums of current row block of global system matrix
      else
      {
        // loop over all column blocks of global system matrix
        for (int j = 0; j < blocksparsematrix->cols(); ++j)
        {
          // extract current block of global system matrix
          const Core::LinAlg::SparseMatrix& matrix = blocksparsematrix->matrix(i, j);

          // loop over all rows of current matrix block
          for (int irow = 0; irow < matrix.row_map().NumMyElements(); ++irow)
          {
            // determine length of current matrix row
            const int length = matrix.epetra_matrix()->NumMyEntries(irow);

            if (length > 0)
            {
              // extract current matrix row from matrix block
              int numentries(0);
              std::vector<double> values(length, 0.);
              if (matrix.epetra_matrix()->ExtractMyRowCopy(irow, length, numentries, values.data()))
                FOUR_C_THROW("Cannot extract matrix row with local ID {} from matrix block!", irow);

              // compute and store current row sum
              double rowsum(0.);
              for (int ientry = 0; ientry < numentries; ++ientry)
                rowsum += std::abs(values[ientry]);
              (invrowsums)[irow] += rowsum;
            }
          }
        }

        // invert row sums
        if (invrowsums.reciprocal(invrowsums)) FOUR_C_THROW("Vector could not be inverted!");

        // take square root of inverse row sums if matrix is scaled from left and right
        if (method() == EquilibrationMethod::rowsandcolumns_full or
            method() == EquilibrationMethod::rowsandcolumns_maindiag)
          for (int j = 0; j < invrowsums.local_length(); ++j)
            (invrowsums)[j] = std::sqrt((invrowsums)[j]);
      }

      // perform row equilibration of matrix blocks in current row block of global system
      // matrix
      for (int j = 0; j < blocksparsematrix->cols(); ++j)
        equilibrate_matrix_rows(blocksparsematrix->matrix(i, j), invrowsums);

      // insert inverse row sums of current main diagonal matrix block into global vector
      blockmaps->insert_vector(invrowsums, i, *invrowsums_);
    }
  }

  // perform column equilibration
  if (method() == EquilibrationMethod::columns_full or
      method() == EquilibrationMethod::columns_maindiag or
      method() == EquilibrationMethod::rowsandcolumns_full or
      method() == EquilibrationMethod::rowsandcolumns_maindiag)
  {
    for (int j = 0; j < blocksparsematrix->cols(); ++j)
    {
      // initialize vector for inverse column sums
      std::shared_ptr<Core::LinAlg::Vector<double>> invcolsums(
          std::make_shared<Core::LinAlg::Vector<double>>(
              blocksparsematrix->matrix(j, j).domain_map()));

      // compute inverse column sums of current main diagonal matrix block
      if (method() == EquilibrationMethod::columns_maindiag or
          method() == EquilibrationMethod::rowsandcolumns_maindiag)
      {
        compute_inv_col_sums(blocksparsematrix->matrix(j, j), *invcolsums, method());
      }
      // compute inverse column sums of current column block of global system matrix
      else
      {
        // loop over all row blocks of global system matrix
        for (int i = 0; i < blocksparsematrix->rows(); ++i)
        {
          // extract current block of global system matrix
          const Core::LinAlg::SparseMatrix& matrix = blocksparsematrix->matrix(i, j);

          // loop over all rows of current matrix block
          for (int irow = 0; irow < matrix.row_map().NumMyElements(); ++irow)
          {
            // determine length of current matrix row
            const int length = matrix.epetra_matrix()->NumMyEntries(irow);

            if (length > 0)
            {
              // extract current matrix row from matrix block
              int numentries(0);
              std::vector<double> values(length, 0.);
              std::vector<int> indices(length, 0);
              if (matrix.epetra_matrix()->ExtractMyRowCopy(
                      irow, length, numentries, values.data(), indices.data()))
                FOUR_C_THROW("Cannot extract matrix row with local ID {} from matrix block!", irow);

              // add entries of current matrix row to column sums
              for (int ientry = 0; ientry < numentries; ++ientry)
                invcolsums->sum_into_global_value(
                    matrix.col_map().GID(indices[ientry]), 0, std::abs(values[ientry]));
            }
          }
        }

        // invert column sums
        if (invcolsums->reciprocal(*invcolsums)) FOUR_C_THROW("Vector could not be inverted!");

        // take square root of inverse column sums if matrix is scaled from left and right
        if (method() == EquilibrationMethod::rowsandcolumns_full or
            method() == EquilibrationMethod::rowsandcolumns_maindiag)
          for (int i = 0; i < invcolsums->local_length(); ++i)
            (*invcolsums)[i] = std::sqrt((*invcolsums)[i]);
      }

      // perform column equilibration of matrix blocks in current column block of global
      // system matrix
      for (int i = 0; i < blocksparsematrix->rows(); ++i)
        equilibrate_matrix_columns(blocksparsematrix->matrix(i, j), *invcolsums);

      // insert inverse column sums of current main diagonal matrix block into global vector
      blockmaps->insert_vector(*invcolsums, j, *invcolsums_);
    }
  }
}

/*-------------------------------------------------------------------------*
 *-------------------------------------------------------------------------*/
void Core::LinAlg::EquilibrationBlockSpecific::equilibrate_matrix(
    std::shared_ptr<Core::LinAlg::SparseOperator> systemmatrix,
    std::shared_ptr<const Core::LinAlg::MultiMapExtractor> blockmaps) const
{
  std::shared_ptr<Core::LinAlg::BlockSparseMatrixBase> blocksparsematrix =
      Core::LinAlg::cast_to_block_sparse_matrix_base_and_check_success(systemmatrix);

  if (blocksparsematrix->rows() != static_cast<int>(method_blocks_.size()))
    FOUR_C_THROW("No match between number of equilibration methods and Matrix blocks");

  // init: no scaling
  invcolsums_->put_scalar(1.0);
  invrowsums_->put_scalar(1.0);

  // loop over all blocks of matrix and apply equilibration for each block
  for (int i = 0; i < blocksparsematrix->rows(); ++i)
  {
    const EquilibrationMethod method = method_blocks_.at(i);
    if (method == EquilibrationMethod::rows_maindiag or
        method == EquilibrationMethod::rowsandcolumns_maindiag)
    {
      Vector<double> invrowsums(blocksparsematrix->matrix(i, i).row_map());
      compute_inv_row_sums(blocksparsematrix->matrix(i, i), invrowsums, method);

      // perform row equilibration of matrix blocks in current row block of global system
      // matrix
      for (int j = 0; j < blocksparsematrix->cols(); ++j)
        equilibrate_matrix_rows(blocksparsematrix->matrix(i, j), invrowsums);

      // insert inverse row sums of current main diagonal matrix block into global vector
      blockmaps->insert_vector(invrowsums, i, *invrowsums_);
    }
    if (method == EquilibrationMethod::columns_maindiag or
        method == EquilibrationMethod::rowsandcolumns_maindiag)
    {
      Vector<double> invcolsums(blocksparsematrix->matrix(i, i).domain_map());
      compute_inv_col_sums(blocksparsematrix->matrix(i, i), invcolsums, method);

      // perform column equilibration of matrix blocks in current column block of global
      // system matrix
      for (int j = 0; j < blocksparsematrix->cols(); ++j)
        equilibrate_matrix_columns(blocksparsematrix->matrix(j, i), invcolsums);

      // insert inverse column sums of current main diagonal matrix block into global vector
      blockmaps->insert_vector(invcolsums, i, *invcolsums_);
    }
    if (method == EquilibrationMethod::symmetry)
    {
      auto invsymmetry =
          Core::LinAlg::create_vector(blocksparsematrix->matrix(i, i).range_map(), true);

      compute_inv_symmetry(blocksparsematrix->matrix(i, i), *invsymmetry);

      for (int j = 0; j < blocksparsematrix->cols(); ++j)
      {
        equilibrate_matrix_rows(blocksparsematrix->matrix(i, j), *invsymmetry);
        equilibrate_matrix_columns(blocksparsematrix->matrix(j, i), *invsymmetry);
      }

      // insert inverse row sums of current main diagonal matrix block into global vector
      blockmaps->insert_vector(*invsymmetry, i, *invcolsums_);
      blockmaps->insert_vector(*invsymmetry, i, *invrowsums_);
    }
  }
}

/*-------------------------------------------------------------------------*
 *-------------------------------------------------------------------------*/
std::shared_ptr<Core::LinAlg::Equilibration> Core::LinAlg::build_equilibration(MatrixType type,
    const std::vector<EquilibrationMethod>& method, std::shared_ptr<const Epetra_Map> dofrowmap)
{
  std::shared_ptr<Core::LinAlg::Equilibration> equilibration = nullptr;

  if (method.size() == 1)
  {
    EquilibrationMethod method_global = method.at(0);

    if (method_global == Core::LinAlg::EquilibrationMethod::none)
      equilibration = std::make_shared<Core::LinAlg::EquilibrationNone>(dofrowmap);
    else
    {
      switch (type)
      {
        case Core::LinAlg::MatrixType::sparse:
        {
          equilibration =
              std::make_shared<Core::LinAlg::EquilibrationSparse>(method_global, dofrowmap);
          break;
        }
        case Core::LinAlg::MatrixType::block_field:
        case Core::LinAlg::MatrixType::block_condition:
        case Core::LinAlg::MatrixType::block_condition_dof:
        {
          equilibration =
              std::make_shared<Core::LinAlg::EquilibrationBlock>(method_global, dofrowmap);
          break;
        }
        default:
        {
          FOUR_C_THROW("Unknown type of global system matrix");
          break;
        }
      }
    }
  }
  else
    equilibration = std::make_shared<Core::LinAlg::EquilibrationBlockSpecific>(method, dofrowmap);

  return equilibration;
}

FOUR_C_NAMESPACE_CLOSE
