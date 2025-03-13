// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_linalg_vector.hpp"

#include "4C_comm_mpi_utils.hpp"
#include "4C_linalg_multi_vector.hpp"
#include "4C_utils_exceptions.hpp"

#include <Epetra_Vector.h>


FOUR_C_NAMESPACE_OPEN

template <typename T>
Core::LinAlg::Vector<T>::Vector(const Epetra_BlockMap& Map, bool zeroOut)
    : vector_(std::make_shared<Epetra_Vector>(Map, zeroOut))
{
}

template <typename T>
Core::LinAlg::Vector<T>::Vector(const Epetra_Vector& Source)
    : vector_(std::make_shared<Epetra_Vector>(Source))
{
}

template <typename T>
Core::LinAlg::Vector<T>::Vector(const Epetra_FEVector& Source)
    : vector_(std::make_shared<Epetra_Vector>(Epetra_DataAccess::Copy, Source, 0))
{
  FOUR_C_ASSERT(Source.NumVectors() == 1, "Can only convert a FE vector with a single column.");
}

template <typename T>
Core::LinAlg::Vector<T>::Vector(const Vector& other)
    : vector_(std::make_shared<Epetra_Vector>(other.get_ref_of_epetra_vector()))
{
}


template <typename T>
Core::LinAlg::Vector<T>& Core::LinAlg::Vector<T>::operator=(const Vector& other)
{
  *vector_ = other.get_ref_of_epetra_vector();
  return *this;
}


template <typename T>
Core::LinAlg::Vector<T>::operator const Core::LinAlg::MultiVector<T>&() const
{
  sync_view();
  FOUR_C_ASSERT(multi_vector_view_ != nullptr, "Internal error.");
  return *multi_vector_view_;
}


template <typename T>
Core::LinAlg::Vector<T>::operator Core::LinAlg::MultiVector<T>&()
{
  sync_view();
  FOUR_C_ASSERT(multi_vector_view_ != nullptr, "Internal error.");
  return *multi_vector_view_;
}


template <typename T>
int Core::LinAlg::Vector<T>::norm_1(double* Result) const
{
  return vector_->Norm1(Result);
}

template <typename T>
int Core::LinAlg::Vector<T>::norm_2(double* Result) const
{
  return vector_->Norm2(Result);
}

template <typename T>
int Core::LinAlg::Vector<T>::norm_inf(double* Result) const
{
  return vector_->NormInf(Result);
}

template <typename T>
int Core::LinAlg::Vector<T>::min_value(double* Result) const
{
  return vector_->MinValue(Result);
}

template <typename T>
int Core::LinAlg::Vector<T>::max_value(double* Result) const
{
  return vector_->MaxValue(Result);
}

template <typename T>
int Core::LinAlg::Vector<T>::mean_value(double* Result) const
{
  return vector_->MeanValue(Result);
}

template <typename T>
int Core::LinAlg::Vector<T>::dot(const Epetra_MultiVector& A, double* Result) const
{
  return vector_->Dot(A, Result);
}

template <typename T>
int Core::LinAlg::Vector<T>::abs(const Epetra_MultiVector& A)
{
  return vector_->Abs(A);
}

template <typename T>
int Core::LinAlg::Vector<T>::scale(double ScalarA, const Epetra_MultiVector& A)
{
  return vector_->Scale(ScalarA, A);
}

template <typename T>
int Core::LinAlg::Vector<T>::update(double ScalarA, const Epetra_MultiVector& A, double ScalarThis)
{
  return vector_->Update(ScalarA, A, ScalarThis);
}

template <typename T>
int Core::LinAlg::Vector<T>::update(double ScalarA, const Epetra_MultiVector& A, double ScalarB,
    const Epetra_MultiVector& B, double ScalarThis)
{
  return vector_->Update(ScalarA, A, ScalarB, B, ScalarThis);
}


template <typename T>
int Core::LinAlg::Vector<T>::dot(const Vector& A, double* Result) const
{
  return vector_->Dot(A, Result);
}

template <typename T>
int Core::LinAlg::Vector<T>::abs(const Vector& A)
{
  return vector_->Abs(A);
}

template <typename T>
int Core::LinAlg::Vector<T>::scale(double ScalarA, const Vector& A)
{
  return vector_->Scale(ScalarA, A);
}

template <typename T>
int Core::LinAlg::Vector<T>::update(double ScalarA, const Vector& A, double ScalarThis)
{
  return vector_->Update(ScalarA, A, ScalarThis);
}

template <typename T>
int Core::LinAlg::Vector<T>::update(
    double ScalarA, const Vector& A, double ScalarB, const Vector& B, double ScalarThis)
{
  return vector_->Update(ScalarA, A, ScalarB, B.get_ref_of_epetra_vector(), ScalarThis);
}

template <typename T>
int Core::LinAlg::Vector<T>::put_scalar(double ScalarConstant)
{
  return vector_->PutScalar(ScalarConstant);
}


template <typename T>
int Core::LinAlg::Vector<T>::replace_map(const Epetra_BlockMap& map)
{
  if (multi_vector_view_) multi_vector_view_->ReplaceMap(map);
  return vector_->ReplaceMap(map);
}

template <typename T>
MPI_Comm Core::LinAlg::Vector<T>::get_comm() const
{
  return Core::Communication::unpack_epetra_comm(vector_->Comm());
}

template <typename T>
void Core::LinAlg::Vector<T>::sync_view() const
{
  // Only do this once.
  if (!multi_vector_view_) multi_vector_view_ = MultiVector<T>::create_view(*vector_);
}

// explicit instantiation
template class Core::LinAlg::Vector<double>;



Core::LinAlg::Vector<int>::Vector(const Epetra_BlockMap& map, bool zeroOut)
    : vector_(std::make_shared<Epetra_IntVector>(map, zeroOut))
{
}

Core::LinAlg::Vector<int>::Vector(const Epetra_BlockMap& map, int* values)
    : vector_(std::make_shared<Epetra_IntVector>(Epetra_DataAccess::Copy, map, values))
{
}

Core::LinAlg::Vector<int>::Vector(const Vector& other)
    : vector_(std::make_shared<Epetra_IntVector>(*other.vector_))
{
}

Core::LinAlg::Vector<int>::Vector(Vector&& other) noexcept : vector_(std::move(other.vector_)) {}

Core::LinAlg::Vector<int>& Core::LinAlg::Vector<int>::operator=(const Vector& other)
{
  vector_ = std::make_shared<Epetra_IntVector>(*other.vector_);
  return *this;
}

Core::LinAlg::Vector<int>& Core::LinAlg::Vector<int>::operator=(Vector&& other) noexcept
{
  vector_ = std::move(other.vector_);
  return *this;
}


int Core::LinAlg::Vector<int>::put_value(int Value) { return vector_->PutValue(Value); }

int Core::LinAlg::Vector<int>::max_value() { return vector_->MaxValue(); }

int Core::LinAlg::Vector<int>::min_value() { return vector_->MinValue(); }

void Core::LinAlg::Vector<int>::print(std::ostream& os) const { vector_->Print(os); }

MPI_Comm Core::LinAlg::Vector<int>::get_comm() const
{
  return Core::Communication::unpack_epetra_comm(vector_->Comm());
}


FOUR_C_NAMESPACE_CLOSE
