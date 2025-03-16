// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_LINALG_BLOCKSPARSEMATRIX_HPP
#define FOUR_C_LINALG_BLOCKSPARSEMATRIX_HPP

#include "4C_config.hpp"

#include "4C_comm_mpi_utils.hpp"
#include "4C_linalg_sparsematrix.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Core::LinAlg
{
  /// Internal base class of BlockSparseMatrix that contains the non-template stuff
  /*!

    This is where the bookkeeping of the BlockSparseMatrix happens. We use two
    MultiMapExtractor objects to store the FullRangeMap() and the
    FullDomainMap() along with their many partial RangeMap() and
    DomainMap(). Most of the required SparseOperator methods can simply be
    implemented in terms of the matrix blocks.

   */
  class BlockSparseMatrixBase : public SparseOperator
  {
   public:
    /// constructor
    /*!
      \param domainmaps domain maps for all blocks
      \param rangemaps range maps for all blocks
      \param npr estimated number of entries per row in each block
      \param explicitdirichlet whether to remove Dirichlet zeros from the
      matrix graphs in each block
      \param savegraph whether to save the matrix graphs of each block and
      recreate filled matrices the next time
     */
    BlockSparseMatrixBase(const MultiMapExtractor& domainmaps, const MultiMapExtractor& rangemaps,
        int npr, bool explicitdirichlet = true, bool savegraph = false);


    /// make a copy of me
    virtual std::unique_ptr<BlockSparseMatrixBase> clone(DataAccess access) = 0;

    /// destroy the underlying Epetra objects
    virtual bool destroy(bool throw_exception_for_blocks = true);

    /// setup of block preconditioners
    /*!
      This method can be implemented by subclasses that implement
      ApplyInverse() to execute a block preconditioner on the matrix.
     */
    virtual void setup_preconditioner() {}

    /// Merge block matrix into a SparseMatrix
    std::shared_ptr<SparseMatrix> merge(bool explicitdirichlet = true) const;

    /** \name Block matrix access */
    //@{

    /// return block (r,c)
    const SparseMatrix& matrix(int r, int c) const { return blocks_[r * cols() + c]; }

    /// return block (r,c)
    SparseMatrix& matrix(int r, int c) { return blocks_[r * cols() + c]; }

    /// return block (r,c)
    inline const SparseMatrix& operator()(int r, int c) const { return matrix(r, c); }

    /// return block (r,c)
    inline SparseMatrix& operator()(int r, int c) { return matrix(r, c); }

    /// assign SparseMatrix to block (r,c)
    /*!
      \note The maps of the block have to match the maps of the given matrix.
     */
    void assign(int r, int c, DataAccess access, const SparseMatrix& mat);

    //@}

    /** \name FE methods */
    //@{

    void zero() override;
    void reset() override;

    void complete(bool enforce_complete = false) override;

    void complete(const Epetra_Map& domainmap, const Epetra_Map& rangemap,
        bool enforce_complete = false) override;

    void un_complete() override;

    void apply_dirichlet(
        const Core::LinAlg::Vector<double>& dbctoggle, bool diagonalblock = true) override;

    void apply_dirichlet(const Epetra_Map& dbcmap, bool diagonalblock = true) override;

    /// derived
    bool is_dbc_applied(const Epetra_Map& dbcmap, bool diagonalblock = true,
        const Core::LinAlg::SparseMatrix* trafo = nullptr) const override;

    //@}

    /** \name Matrix Properties Query Methods */
    //@{

    /// If Complete() has been called, this query returns true, otherwise it returns false.
    bool filled() const override;

    //@}

    /** \name Block maps */
    //@{

    /// number of row blocks
    int rows() const { return rangemaps_.num_maps(); }

    /// number of column blocks
    int cols() const { return domainmaps_.num_maps(); }

    /// range map for given row block
    const Epetra_Map& range_map(int r) const { return *rangemaps_.Map(r); }

    /// domain map for given column block
    const Epetra_Map& domain_map(int r) const { return *domainmaps_.Map(r); }

    /// total matrix range map with all blocks
    const Epetra_Map& full_range_map() const { return *rangemaps_.full_map(); }

    /// total matrix domain map with all blocks
    const Epetra_Map& full_domain_map() const { return *domainmaps_.full_map(); }

    /// total matrix domain map with all blocks (this is needed for
    /// consistency with Core::LinAlg::SparseMatrix)
    const Epetra_Map& domain_map() const override { return *domainmaps_.full_map(); }

    /// total matrix row map with all blocks
    /*!
      \pre Filled()==true
     */
    Epetra_Map& full_row_map() const { return *fullrowmap_; }

    /// total matrix column map with all blocks
    /*!
      \pre Filled()==true
     */
    Epetra_Map& full_col_map() const { return *fullcolmap_; }

    //@}

    /** \name Attribute set methods */
    //@{

    /// If set true, transpose of this operator will be applied.
    int SetUseTranspose(bool UseTranspose) override;

    //@}

    /** \name Mathematical functions */
    //@{

    /// Returns the result of a Epetra_Operator applied to a Epetra_MultiVector X in Y.
    int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const override;

    /// Returns the result of a Epetra_Operator inverse applied to an Epetra_MultiVector X in Y.
    int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const override;

    /// Resolve virtual function of parent class
    int multiply(bool TransA, const Core::LinAlg::MultiVector<double>& X,
        Core::LinAlg::MultiVector<double>& Y) const override;

    /// Add a (transposed) BlockSparseMatrix: (*this) = (*this)*scalarB + A(^T)*scalarA
    virtual void add(const BlockSparseMatrixBase& A, const bool transposeA, const double scalarA,
        const double scalarB);

    /// Resolve virtual function of parent class
    void add(const SparseOperator& A, const bool transposeA, const double scalarA,
        const double scalarB) override;

    /// Resolve virtual function of parent class
    void add_other(SparseMatrixBase& A, const bool transposeA, const double scalarA,
        const double scalarB) const override;

    /// Resolve virtual function of parent class
    void add_other(BlockSparseMatrixBase& A, const bool transposeA, const double scalarA,
        const double scalarB) const override;

    /// Multiply all values in the matrix by a constant value (in place: A <- ScalarConstant * A).
    int scale(double ScalarConstant) override;

    /// Returns the infinity norm of the global matrix.
    double NormInf() const override;

    //@}

    /** \name Attribute access functions */
    //@{

    /// Returns a character string describing the operator.
    const char* Label() const override;

    /// Returns the current UseTranspose setting.
    bool UseTranspose() const override;

    /// Returns true if the this object can provide an approximate Inf-norm, false otherwise.
    bool HasNormInf() const override;

    /// Returns a pointer to the Epetra_Comm communicator associated with this operator.
    const Epetra_Comm& Comm() const override;

    /// Returns the Epetra_Map object associated with the domain of this operator.
    const Epetra_Map& OperatorDomainMap() const override;

    /// Returns the Epetra_Map object associated with the range of this operator.
    const Epetra_Map& OperatorRangeMap() const override;

    //@}

    /// access to domain map extractor in derived classes
    const MultiMapExtractor& domain_extractor() const { return domainmaps_; }

    /// access to range map extractor in derived classes
    const MultiMapExtractor& range_extractor() const { return rangemaps_; }

   protected:
    /// extract a partial map extractor from the full map extractor
    void get_partial_extractor(const MultiMapExtractor& full_extractor,
        const std::vector<unsigned>& block_ids, MultiMapExtractor& partial_extractor) const;

   private:
    /// the full domain map together with all partial domain maps
    MultiMapExtractor domainmaps_;

    /// the full range map together with all partial range maps
    MultiMapExtractor rangemaps_;

    /// row major matrix block storage
    std::vector<SparseMatrix> blocks_;

    /// full matrix row map
    std::shared_ptr<Epetra_Map> fullrowmap_;

    /// full matrix column map
    std::shared_ptr<Epetra_Map> fullcolmap_;

    /// see matrix as transposed
    bool usetranspose_;
  };



  /// Block matrix consisting of SparseMatrix blocks
  /*!
      There are strange algorithms that need to split a large sparse matrix into
      blocks. Such things happen, e.g., in FSI calculations with internal and
      interface splits, in fluid projection preconditioners or in contact
      simulations with slave and master sides. Unfortunately splitting a huge
      sparse matrix in (possibly) many blocks is nontrivial and expensive. So
      the idea here is to assemble into a block matrix in the first place.

      The difficulty with this approach is the handling of ghost entries in a
      parallel matrix. It is hard (expensive) to figure out to which column
      block each particular ghost entry belongs. That is why this class is
      templated with a Strategy. There is a default implementation for this
      template parameter DefaultBlockMatrixStrategy, that handles the most
      general case. That is DefaultBlockMatrixStrategy finds the right column
      block be heavy communication. But if there is some knowledge available in
      a particular case, it is easy to implement a specify Strategy that does
      not need to communicate that much.

   */
  template <class Strategy>
  class BlockSparseMatrix : public BlockSparseMatrixBase, public Strategy
  {
   public:
    BlockSparseMatrix(const MultiMapExtractor& domainmaps, const MultiMapExtractor& rangemaps,
        int npr = 81, bool explicitdirichlet = true, bool savegraph = false);

    /// clone the full block sparse matrix

    /** Do not forget to call Complete() after cloning, even if you
     *  use Core::LinAlg::View! */
    std::unique_ptr<BlockSparseMatrixBase> clone(DataAccess access) override;

    /// clone only a part of the block sparse matrix
    /** Do not forget to call Complete() after cloning, even if you
     *  use Core::LinAlg::View!
     *
     *  \param[in] access : consider copy or view of block matrices
     *  \param[in] row_block_ids : ID's of the row blocks to clone
     *  \param[in] col_block_ids : ID's of the column blocks to clone
     *
     *  */
    std::unique_ptr<Core::LinAlg::BlockSparseMatrixBase> clone(DataAccess access,
        const std::vector<unsigned>& row_block_ids, const std::vector<unsigned>& col_block_ids);

    /// just a dummy that switches from strided assembly to standard assembly
    void assemble(int eid, const std::vector<int>& lmstride,
        const Core::LinAlg::SerialDenseMatrix& Aele, const std::vector<int>& lmrow,
        const std::vector<int>& lmrowowner, const std::vector<int>& lmcol) override
    {
      const int myrank =
          Core::Communication::my_mpi_rank(Core::Communication::unpack_epetra_comm(Comm()));
      Strategy::assemble(eid, myrank, lmstride, Aele, lmrow, lmrowowner, lmcol);
    }

    /// single value assemble
    /*!

       \warning This method is less useful here. We just need to make the
       compiler happy. We assume "element matrices" of size 1x1 here. Do not use
       this method if your strategy depends on the real position of the dof in
       the element matrix.

     */
    void assemble(double val, int rgid, int cgid) override { Strategy::assemble(val, rgid, cgid); }

    void complete(bool enforce_complete = false) override;

   private:
    /** \brief internal clone method which provides the possibility to clone only
     *         a sub-set of all blocks
     *
     *  This method is not supposed to be called from outside! See public variant.
     *
     *  \param[in] access : consider copy or view of block matrices
     *  \param[in] row_block_ids : ID's of the row blocks to clone
     *  \param[in] col_block_ids : ID's of the column blocks to clone
     *  \param[in] domain_extractor : necessary domain extractor
     *  \param[in] range_extractor : necessary range extractor
     *
     *  */
    std::unique_ptr<Core::LinAlg::BlockSparseMatrixBase> clone(DataAccess access,
        const std::vector<unsigned>& row_block_ids, const std::vector<unsigned>& col_block_ids,
        const MultiMapExtractor& domain_extractor, const MultiMapExtractor& range_extractor);
  };


  /// default strategy implementation for block matrix
  /*!

      This default implementation solves the ghost entry problem by remembering
      all ghost entries during the assembly in a private map. Afterwards
      Complete() needs to be called that finds the appropriate block for each
      ghost entry by communication an finally assembles these entries.

      This is the most general, most expensive implementation. You are
      encouraged to provide your own Strategy implementation if you know your
      specific block structure.

      \sa BlockSparseMatrix

   */
  class DefaultBlockMatrixStrategy
  {
   public:
    /// construct with a block matrix base
    explicit DefaultBlockMatrixStrategy(BlockSparseMatrixBase& mat);

    /// assemble into the given block using nodal strides
    void assemble(int eid, int myrank, const std::vector<int>& lmstride,
        const Core::LinAlg::SerialDenseMatrix& Aele, const std::vector<int>& lmrow,
        const std::vector<int>& lmrowowner, const std::vector<int>& lmcol);

    /// assemble into the given block
    void assemble(double val, int rgid, int cgid);

    /// assemble the remaining ghost entries
    void complete(bool enforce_complete = false);

   protected:
    /// assemble into the given block
    void assemble(double val, int lrow, int rgid, int rblock, int lcol, int cgid, int cblock);

    /// find row block to a given row gid
    int row_block(int rgid);

    /// find column block to a given column gid
    int col_block(int rblock, int cgid);

    /// access to the block sparse matrix for subclasses
    BlockSparseMatrixBase& mat() { return mat_; }

   private:
    /// my block matrix base
    BlockSparseMatrixBase& mat_;

    /// all ghost entries stored by row,column
    std::map<int, std::map<int, double>> ghost_;

    /// scratch array for identifying column information
    std::vector<std::vector<int>> scratch_lcols_;
  };

  /*----------------------------------------------------------------------*
   *----------------------------------------------------------------------*/

  //////////////////////////////////
  /// helper functions

  std::shared_ptr<BlockSparseMatrix<DefaultBlockMatrixStrategy>> block_matrix2x2(
      SparseMatrix& A00, SparseMatrix& A01, SparseMatrix& A10, SparseMatrix& A11);

  //! Cast matrix of type SparseOperator to BlockSparseMatrixBase and check in debug mode if cast
  //! was successful
  std::shared_ptr<Core::LinAlg::BlockSparseMatrixBase>
  cast_to_block_sparse_matrix_base_and_check_success(
      std::shared_ptr<Core::LinAlg::SparseOperator> input_matrix);

  //! Cast matrix of type SparseOperator to const BlockSparseMatrixBase and check in debug mode if
  //! cast was successful
  std::shared_ptr<const Core::LinAlg::BlockSparseMatrixBase>
  cast_to_const_block_sparse_matrix_base_and_check_success(
      std::shared_ptr<const Core::LinAlg::SparseOperator> input_matrix);

  //////////////////////////////////



}  // end of namespace Core::LinAlg


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <class Strategy>
Core::LinAlg::BlockSparseMatrix<Strategy>::BlockSparseMatrix(const MultiMapExtractor& domainmaps,
    const MultiMapExtractor& rangemaps, int npr, bool explicitdirichlet, bool savegraph)
    : BlockSparseMatrixBase(domainmaps, rangemaps, npr, explicitdirichlet, savegraph),
      // this was necessary, otherwise ambiguous with copy constructor of Strategy
      Strategy((Core::LinAlg::BlockSparseMatrixBase&)(*this))
{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <class Strategy>
std::unique_ptr<Core::LinAlg::BlockSparseMatrixBase>
Core::LinAlg::BlockSparseMatrix<Strategy>::clone(DataAccess access)
{
  std::vector<unsigned> row_block_ids(rows());
  for (unsigned i = 0; i < static_cast<unsigned>(rows()); ++i) row_block_ids[i] = i;

  std::vector<unsigned> col_block_ids(cols());
  for (unsigned i = 0; i < static_cast<unsigned>(cols()); ++i) col_block_ids[i] = i;

  return clone(access, row_block_ids, col_block_ids, domain_extractor(), range_extractor());
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <class Strategy>
std::unique_ptr<Core::LinAlg::BlockSparseMatrixBase>
Core::LinAlg::BlockSparseMatrix<Strategy>::clone(DataAccess access,
    const std::vector<unsigned>& row_block_ids, const std::vector<unsigned>& col_block_ids,
    const MultiMapExtractor& domain_extractor, const MultiMapExtractor& range_extractor)
{
  int npr = matrix(0, 0).max_num_entries();
  bool explicitdirichlet = matrix(0, 0).explicit_dirichlet();
  bool savegraph = matrix(0, 0).save_graph();
  std::unique_ptr<BlockSparseMatrixBase> bsm = std::make_unique<BlockSparseMatrix<Strategy>>(
      domain_extractor, range_extractor, npr, explicitdirichlet, savegraph);

  for (std::vector<unsigned>::const_iterator r = row_block_ids.begin(); r != row_block_ids.end();
      ++r)
  {
    for (std::vector<unsigned>::const_iterator c = col_block_ids.begin(); c != col_block_ids.end();
        ++c)
    {
      bsm->matrix(*r, *c).assign(access, matrix(*r, *c));
    }
  }
  return bsm;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <class Strategy>
std::unique_ptr<Core::LinAlg::BlockSparseMatrixBase>
Core::LinAlg::BlockSparseMatrix<Strategy>::clone(DataAccess access,
    const std::vector<unsigned>& row_block_ids, const std::vector<unsigned>& col_block_ids)
{
  if (std::lower_bound(row_block_ids.begin(), row_block_ids.end(), static_cast<unsigned>(rows())) !=
      row_block_ids.end())
    FOUR_C_THROW("The partial row block ids exceed the maximal possible id!");

  if (std::lower_bound(col_block_ids.begin(), col_block_ids.end(), static_cast<unsigned>(cols())) !=
      col_block_ids.end())
    FOUR_C_THROW("The partial column block ids exceed the maximal possible id!");

  if (row_block_ids.size() == 0 or col_block_ids.size() == 0)
    FOUR_C_THROW("The provided row/col block id vector has a length of zero!");

  // extract domain extractors
  MultiMapExtractor p_domain_extractor;
  this->get_partial_extractor(domain_extractor(), col_block_ids, p_domain_extractor);

  // extract range extractors
  MultiMapExtractor p_range_extractor;
  this->get_partial_extractor(range_extractor(), row_block_ids, p_range_extractor);

  return clone(access, row_block_ids, col_block_ids, p_domain_extractor, p_range_extractor);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <class Strategy>
void Core::LinAlg::BlockSparseMatrix<Strategy>::complete(bool enforce_complete)
{
  Strategy::complete();
  BlockSparseMatrixBase::complete(enforce_complete);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
inline int Core::LinAlg::DefaultBlockMatrixStrategy::row_block(int rgid)
{
  int rows = mat_.rows();
  for (int rblock = 0; rblock < rows; ++rblock)
  {
    if (mat_.range_map(rblock).MyGID(rgid))
    {
      return rblock;
    }
  }
  return -1;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
inline int Core::LinAlg::DefaultBlockMatrixStrategy::col_block(int rblock, int cgid)
{
  int cols = mat_.cols();
  for (int cblock = 0; cblock < cols; ++cblock)
  {
    SparseMatrix& matrix = mat_.matrix(rblock, cblock);

    // If we have a filled matrix we know the column map already.
    if (matrix.filled())
    {
      if (matrix.col_map().MyGID(cgid))
      {
        return cblock;
      }
    }

    // otherwise we can get just the non-ghost entries right now
    else if (mat_.domain_map(cblock).MyGID(cgid))
    {
      return cblock;
    }
  }

  // ghost entries in a non-filled matrix will have to be done later

  return -1;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
inline void Core::LinAlg::DefaultBlockMatrixStrategy::assemble(int eid, int myrank,
    const std::vector<int>& lmstride, const Core::LinAlg::SerialDenseMatrix& Aele,
    const std::vector<int>& lmrow, const std::vector<int>& lmrowowner,
    const std::vector<int>& lmcol)
{
  const int lrowdim = (int)lmrow.size();
  const int lcoldim = (int)lmcol.size();

  FOUR_C_ASSERT(
      static_cast<int>(scratch_lcols_.size()) == mat_.rows(), "Unexpected number of block rows");

  for (int rblock = 0; rblock < mat_.rows(); ++rblock)
  {
    scratch_lcols_[rblock].resize(lcoldim);
    std::fill(scratch_lcols_[rblock].begin(), scratch_lcols_[rblock].end(), -1);
  }

  // loop rows of local matrix
  for (int lrow = 0; lrow < lrowdim; ++lrow)
  {
    // check ownership of row
    if (lmrowowner[lrow] != myrank) continue;

    int rgid = lmrow[lrow];
    int rblock = row_block(rgid);

    if (scratch_lcols_[rblock][0] == -1)
      for (int lcol = 0; lcol < lcoldim; ++lcol)
        scratch_lcols_[rblock][lcol] = col_block(rblock, lmcol[lcol]);

    for (int lcol = 0; lcol < lcoldim; ++lcol)
    {
      double val = Aele(lrow, lcol);
      int cgid = lmcol[lcol];

      assemble(val, lrow, rgid, rblock, lcol, cgid, scratch_lcols_[rblock][lcol]);
    }
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
inline void Core::LinAlg::DefaultBlockMatrixStrategy::assemble(double val, int rgid, int cgid)
{
  int rblock = row_block(rgid);
  int cblock = col_block(rblock, cgid);

  assemble(val, 0, rgid, rblock, 0, cgid, cblock);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
inline void Core::LinAlg::DefaultBlockMatrixStrategy::assemble(
    double val, int lrow, int rgid, int rblock, int lcol, int cgid, int cblock)
{
#ifdef FOUR_C_ENABLE_ASSERTIONS
  if (rblock == -1) FOUR_C_THROW("no block entry found for row gid={}", rgid);
#endif

  if (cblock > -1)
  {
    SparseMatrix& matrix = mat_.matrix(rblock, cblock);
    matrix.assemble(val, rgid, cgid);
  }
  else
  {
    // ghost entry in non-filled matrix. Save for later insertion.
    ghost_[rgid][cgid] += val;
  }
}

FOUR_C_NAMESPACE_CLOSE

#endif
