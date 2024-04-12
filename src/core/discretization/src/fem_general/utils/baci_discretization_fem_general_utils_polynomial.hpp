/*----------------------------------------------------------------------*/
/*! \file

 \brief Generic polynomials for HDG methods in 1D, 2D, 3D

\level 2

 */

#ifndef FOUR_C_DISCRETIZATION_FEM_GENERAL_UTILS_POLYNOMIAL_HPP
#define FOUR_C_DISCRETIZATION_FEM_GENERAL_UTILS_POLYNOMIAL_HPP

#include "baci_config.hpp"

#include "baci_discretization_fem_general_utils_local_connectivity_matrices.hpp"
#include "baci_linalg_serialdensevector.hpp"

#include <Teuchos_SerialDenseSolver.hpp>

#include <numeric>
#include <utility>

BACI_NAMESPACE_OPEN

namespace CORE::FE
{
  /*!
   \brief helper holding the parameters, which we need to know, to construct a polynomial space
   */
  struct PolynomialSpaceParams
  {
    explicit PolynomialSpaceParams(
        CORE::FE::CellType distype, unsigned int degree, bool completeSpace)
        : distype_(distype), degree_(degree), completeSpace_(completeSpace)
    {
    }

    CORE::FE::CellType distype_;
    unsigned int degree_;
    bool completeSpace_;

    bool operator<(const PolynomialSpaceParams &otherparams) const
    {
      return std::tie(distype_, degree_, completeSpace_) <
             std::tie(otherparams.distype_, otherparams.degree_, otherparams.completeSpace_);
    }
  };

  /*!
   \brief A class of 1D Lagrange polynomials on a given number of points
   */
  class LagrangePolynomial
  {
   public:
    /*!
     \brief Constructor.
     */
    LagrangePolynomial(const std::vector<double> &supportPoints, const double nodePoint)
        : supportPoints_(supportPoints), nodePoint_(nodePoint)
    {
      long double weight = 1.;
      for (double supportPoint : supportPoints) weight *= nodePoint - supportPoint;
      dsassert((std::abs(weight) > std::numeric_limits<double>::min() &&
                   std::abs(weight) < std::numeric_limits<double>::max()),
          "Error in evaluation of polynomial");
      weight_ = 1. / weight;
    }

    /*!
     \brief Copy constructor
     */
    LagrangePolynomial(const LagrangePolynomial &other)
        : supportPoints_(other.supportPoints_), nodePoint_(other.nodePoint_), weight_(other.weight_)
    {
    }

    /*!
     \brief Assignment operator
     */
    LagrangePolynomial &operator=(const LagrangePolynomial &other)
    {
      supportPoints_ = other.supportPoints_;
      nodePoint_ = other.nodePoint_;
      weight_ = other.weight_;
      return *this;
    }

    /*!
     \brief Evaluates the polynomial on a given point on the unit interval [-1, 1].
     */
    [[nodiscard]] double Evaluate(const double point) const
    {
      double value = weight_;
      for (double supportPoint : supportPoints_) value *= point - supportPoint;
      return value;
    }

    /*!
     \brief Evaluates the polynomial and its derivatives (up to the order given by
            the length of values) on the given point
     */
    template <typename M>
    void Evaluate(const double point, M &derivatives) const
    {
      dsassert(derivatives.numCols() == 1, "Only column vectors supported");
      if (derivatives.numRows() == 1)
      {
        derivatives(0) = Evaluate(point);
        return;
      }

      // compute the value and derivatives by expanding the derivatives
      derivatives(0) = weight_;
      for (unsigned int k = 1; k < derivatives.numRows(); ++k) derivatives(k) = 0;
      for (double supportPoint : supportPoints_)
      {
        const double v = point - supportPoint;
        for (int k = derivatives.numRows() - 1; k > 0; --k)
          derivatives(k) = v * derivatives(k) + derivatives(k - 1);
        derivatives(0) *= v;
      }
      double faculty = 1;
      for (unsigned int k = 2; k < derivatives.numRows(); ++k)
      {
        faculty *= static_cast<double>(k);
        derivatives(k) *= faculty;
      }
    }

    [[nodiscard]] double NodePoint() const { return nodePoint_; }

   private:
    std::vector<double> supportPoints_;
    double nodePoint_;
    double weight_;
  };



  /*!
   \brief A class of polynomials based on monomial coefficients

   This class takes a vector of coefficients of the monomials and evaluates a
   polynomial from these coefficients. The evaluation is done by the Horner's method.
   */
  class Polynomial
  {
   public:
    /*!
     \brief Constructor.
     */
    Polynomial(std::vector<double> coefficients) : coefficients_(std::move(coefficients)) {}

    /*!
     \brief Copy constructor
     */
    Polynomial(const Polynomial &other) : coefficients_(other.coefficients_) {}

    /*!
    \brief Assignment operator
    */
    Polynomial &operator=(const Polynomial &other)
    {
      coefficients_ = other.coefficients_;
      return *this;
    }

    /*!
     \brief Evaluates the polynomial on a given point
     */
    [[nodiscard]] double Evaluate(const double point) const
    {
      // classical Horner's method
      return std::accumulate(coefficients_.rbegin(), coefficients_.rend(), 0.0,
          [point](const double &acc, const double &coeff) { return acc * point + coeff; });
    }

    /*!
     \brief Evaluates the polynomial and its derivatives (up to the order given by
         the length of values) on the given point
     */
    template <typename M>
    void Evaluate(const double point, M &derivatives) const
    {
      dsassert(derivatives.numCols() == 1, "Only column vectors supported");
      if (derivatives.numRows() == 1)
      {
        derivatives(0) = Evaluate(point);
        return;
      }

      // compute the value and derivatives by Horner scheme
      std::vector<double> temp(coefficients_);
      const int length = coefficients_.size();
      const int nderivs = derivatives.numRows();
      const int maxder = std::min(length, nderivs);
      double kfaculty = 1.;
      // skip derivatives that are necessarily zero
      for (int k = 0; k < maxder; ++k)
      {
        for (int i = length - 2; i >= k; --i) temp[i] += point * temp[i + 1];
        derivatives(k) = kfaculty * temp[k];
        kfaculty *= static_cast<double>(k + 1);
      }
      for (int k = maxder; k < nderivs; ++k) derivatives(k) = 0;
    }

    /*!
     \brief Evaluates the polynomial or its derivative (of the given order) on the given point
     */
    [[nodiscard]] double EvaluateDerivative(const double point, const int deriv_order) const
    {
      switch (deriv_order)
      {
        case 1:
        {
          CORE::LINALG::Matrix<2, 1> derivatives(false);
          Evaluate(point, derivatives);
          return derivatives(1, 0);
        }
        case 2:
        {
          CORE::LINALG::Matrix<3, 1> derivatives(false);
          Evaluate(point, derivatives);
          return derivatives(2, 0);
        }
        case 3:
        {
          CORE::LINALG::Matrix<4, 1> derivatives(false);
          Evaluate(point, derivatives);
          return derivatives(3, 0);
        }
        default:
        {
          dserror("Only derivatives up to order 2 supported");
          return 0;
        }
      }
    }

    // Here, we do not know about nodes, so put the center of the 1d element
    [[nodiscard]] double NodePoint() const { return 0.; }


   private:
    std::vector<double> coefficients_;
  };



  /*!
   \brief generates complete Lagrange basis in 1D for [-1,1]^d elements
   */
  std::vector<LagrangePolynomial> generateLagrangeBasis1D(const unsigned int degree);



  /*!
   \brief generates complete Legendre basis in 1D

   Legendre polynomials are orthogonal on the unit interval [-1, 1]. As opposed
   to the usual mathematical definition, we also scale the polynomials such that
   they are actually orthonormal on the unit interval.
   */
  std::vector<Polynomial> generateLegendreBasis1D(const unsigned int degree);



  /*!
  \brief Base class for polynomial spaces in nsd_ dimensions used by HDG
  */
  template <int nsd_>
  class PolynomialSpaceBase
  {
   public:
    virtual ~PolynomialSpaceBase() = default;
    /*
     \brief Return the number of polynomials (over all dimensions)
     */
    [[nodiscard]] virtual std::size_t Size() const = 0;

    /*
     \brief Evaluates the values of all polynomials on the given point
     */
    virtual void Evaluate(const CORE::LINALG::Matrix<nsd_, 1> &point,
        CORE::LINALG::SerialDenseVector &values) const = 0;

    /*
     \brief Evaluates the values of all polynomials on the given point
     */
    virtual void Evaluate_deriv1(const CORE::LINALG::Matrix<nsd_, 1> &point,
        CORE::LINALG::SerialDenseMatrix &derivatives) const = 0;

    /*
     \brief Evaluates the first derivative of all polynomials on the given point
     */
    virtual void Evaluate_deriv2(const CORE::LINALG::Matrix<nsd_, 1> &point,
        CORE::LINALG::SerialDenseMatrix &derivatives) const = 0;

    /*
     \brief Creates an array with coordinates of the nodes supporting the polynomials.
     */
    virtual void FillUnitNodePoints(CORE::LINALG::SerialDenseMatrix &matrix) const = 0;
  };



  /*!
   \brief A class of tensor product polynomials expanded from 1D polynomials

   The 1D polynomial class must feature two Evaluate functions, one taking just
   the point as an argument (evaluating the polynomial's value) and one taking
   the point and an array with as many elements as derivatives are requested.

   Base class for LagrangeBasis
   */
  template <int nsd_, class POLY>
  class PolynomialSpaceTensor : public PolynomialSpaceBase<nsd_>
  {
   public:
    /*
     \brief Constructor from a vector of one-dimensional polynomials
     */
    PolynomialSpaceTensor(const std::vector<POLY> polySpace1d) : polySpace1d_(polySpace1d)
    {
      // renumbering is included in case we want to use these polynomials also
      // for standard finite elements
      renumbering_.resize(Size());

      for (unsigned int i = 0; i < renumbering_.size(); ++i) renumbering_[i] = i;
    }


    /*
     \brief Return the number of polynomials (over all dimensions)
     */
    static std::size_t Size(const std::size_t degree)
    {
      std::size_t size = degree + 1;
      for (unsigned int d = 1; d < nsd_; ++d) size *= degree + 1;
      return size;
    }

    /*
     \brief Return the number of polynomials (over all dimensions)
     */
    std::size_t Size() const override { return Size(polySpace1d_.size() - 1); }

    /*
     \brief Evaluates the values of the whole polynomial space in the given point
     */
    void Evaluate(const CORE::LINALG::Matrix<nsd_, 1> &point,
        CORE::LINALG::SerialDenseVector &values) const override;

    /*
     \brief Evaluates the first derivative of the whole polynomial space in the given point
     */
    void Evaluate_deriv1(const CORE::LINALG::Matrix<nsd_, 1> &point,
        CORE::LINALG::SerialDenseMatrix &derivatives) const override;

    /*
     \brief Evaluates the second derivative of the whole polynomial space in the given point
     */
    void Evaluate_deriv2(const CORE::LINALG::Matrix<nsd_, 1> &point,
        CORE::LINALG::SerialDenseMatrix &derivatives) const override;

    /*
     \brief Evaluates the second derivative of the whole polynomial space in the given point
     */
    template <typename M>
    void Evaluate_deriv2(const CORE::LINALG::Matrix<nsd_, 1> &point, M &derivatives) const
    {
    }

    /*
     \brief Convert from a index within the polynomial space to the tensor indices in the
     individual dimensions
     */
    CORE::LINALG::Matrix<nsd_, 1, unsigned int> getIndices(const unsigned int index) const
    {
      dsassert(index < Size(), "Access out of range");
      CORE::LINALG::Matrix<nsd_, 1, unsigned int> indices;
      const unsigned int npoly = polySpace1d_.size();
      switch (nsd_)
      {
        case 1:
          indices(0) = index;
          break;
        case 2:
          indices(0) = index % npoly;
          indices(1) = index / npoly;
          break;
        case 3:
          indices(0) = index % npoly;
          indices(1) = (index / npoly) % npoly;
          indices(2) = index / (npoly * npoly);
          break;
        default:
          dserror("Invalid dimension");
          break;
      }
      return indices;
    }

    /*
     \brief Creates an array with coordinates of the nodes supporting the polynomials.
     */
    void FillUnitNodePoints(CORE::LINALG::SerialDenseMatrix &matrix) const override;

   private:
    std::vector<POLY> polySpace1d_;
    std::vector<int> renumbering_;
  };



  /*!
   \brief A class of polynomials of complete degree expanded from 1D polynomials by truncated
   tensor products

   The 1D polynomial class must feature two Evaluate functions, one taking just
   the point as an argument (evaluating the polynomial's value) and one taking
   the point and an array with as many elements as derivatives are requested.

   Base class for LagrangeBasis
   */
  template <int nsd_, class POLY>
  class PolynomialSpaceComplete : public PolynomialSpaceBase<nsd_>
  {
   public:
    /*
     \brief Constructor from a vector of one-dimensional polynomials
     */
    PolynomialSpaceComplete(const std::vector<POLY> polySpace1d) : polySpace1d_(polySpace1d)
    {
      renumbering_.resize(Size());
      for (unsigned int i = 0; i < renumbering_.size(); ++i) renumbering_[i] = i;
    }


    /*
     \brief Return the number of polynomials (over all dimensions)
     */
    static std::size_t Size(const std::size_t degree)
    {
      std::size_t size = degree + 1;
      for (unsigned int d = 1; d < nsd_; ++d)
      {
        size *= degree + 1 + d;
        size /= (d + 1);  // This integer division is always without remainder
      }
      return size;
    }

    /*
     \brief Return the number of polynomials (over all dimensions)
     */
    std::size_t Size() const override { return Size(polySpace1d_.size() - 1); }

    /*
     \brief Evaluates the values of the whole polynomial space in the given point
     */
    void Evaluate(const CORE::LINALG::Matrix<nsd_, 1> &point,
        CORE::LINALG::SerialDenseVector &values) const override;

    /*
     \brief Evaluates the first derivative of the whole polynomial space in the given point
     */
    void Evaluate_deriv1(const CORE::LINALG::Matrix<nsd_, 1> &point,
        CORE::LINALG::SerialDenseMatrix &derivatives) const override;

    /*
     \brief Evaluates the second derivative of the whole polynomial space in the given point
     */
    void Evaluate_deriv2(const CORE::LINALG::Matrix<nsd_, 1> &point,
        CORE::LINALG::SerialDenseMatrix &derivatives) const override;

    /*
     \brief Creates an array with coordinates of the nodes supporting the polynomials.
     */
    void FillUnitNodePoints(CORE::LINALG::SerialDenseMatrix &matrix) const override;

   private:
    std::vector<POLY> polySpace1d_;
    std::vector<int> renumbering_;
  };



  /*!
   \brief Class of tensor product Lagrange polynomials based on 1D Gauss-Lobatto
   support points (gives standard Q1 and Q2 basis functions and good polynomials
   for higher order). This class is usually not directly called in user code, as
   its functionality can be accessed through PolynomialSpace<nsd_>.
  */
  template <int nsd_>
  class LagrangeBasis : public PolynomialSpaceTensor<nsd_, LagrangePolynomial>
  {
   public:
    /*
    \brief Constructor from a vector of one-dimensional polynomials
     */
    LagrangeBasis(const unsigned int degree)
        : PolynomialSpaceTensor<nsd_, LagrangePolynomial>(generateLagrangeBasis1D(degree))
    {
    }
  };



  /*!
  \brief Class of Legendre polynomials that are pairwise orthogonal. This class
   is usually not directly called in user code, as its functionality can be
   accessed through PolynomialSpace<nsd_>.
  */
  template <int nsd_>
  class LegendreBasis : public PolynomialSpaceComplete<nsd_, Polynomial>
  {
   public:
    /*
   \brief Constructor from a vector of one-dimensional polynomials
     */
    LegendreBasis(const unsigned int degree)
        : PolynomialSpaceComplete<nsd_, Polynomial>(generateLegendreBasis1D(degree))
    {
    }
  };



  /*!
  \brief Lagrange basis for triangular/tetrahedral elements.

   This class constructs a Lagrangian basis for tetrahedral elements of arbitrary degree
   by a transformation of Legendre polynomials onto so-called Fekete nodes, which are
   a generalization of Gauss-Lobatto points to higher dimensions on triangles. For the
   transformation, a Vandermonde matrix is used to map the nodeless orthonormal Legendre
   polynomials onto the desired node distribution.

   The transformation with the Vandermonde matrix is described in
   R. Sevilla, S. Fernandez-Mendez, A. Huerta, NURBS-Enhanced Finite Element Method,
   Arch. Comput. Methods Eng. (2011), 441-484

   Author: kronbichler 08/14
   */
  template <int nsd_>
  class LagrangeBasisTet : public PolynomialSpaceBase<nsd_>
  {
   public:
    using ordinalType = CORE::LINALG::SerialDenseMatrix::ordinalType;
    using scalarType = CORE::LINALG::SerialDenseMatrix::scalarType;

    LagrangeBasisTet(const unsigned int degree) : legendre_(degree)
    {
      FillFeketePoints(degree);
      ComputeVandermondeMatrices(degree);
    }


    /*
     \brief Return the number of polynomials (over all dimensions)
     */
    static std::size_t Size(const std::size_t degree)
    {
      std::size_t size = degree + 1;
      for (unsigned int d = 1; d < nsd_; ++d)
      {
        size *= degree + 1 + d;
        size /= (d + 1);  // This integer division is always without remainder
      }
      return size;
    }

    /*
     \brief Return the number of polynomials (over all dimensions)
     */
    std::size_t Size() const override { return legendre_.Size(); }

    /*
     \brief Evaluates the values of the whole polynomial space in the given point
     */
    void Evaluate(const CORE::LINALG::Matrix<nsd_, 1> &point,
        CORE::LINALG::SerialDenseVector &values) const override;

    /*
     \brief Evaluates the first derivative of the whole polynomial space in the given point
     */
    void Evaluate_deriv1(const CORE::LINALG::Matrix<nsd_, 1> &point,
        CORE::LINALG::SerialDenseMatrix &derivatives) const override;

    /*
     \brief Evaluates the second derivative of the whole polynomial space in the given point
     */
    void Evaluate_deriv2(const CORE::LINALG::Matrix<nsd_, 1> &point,
        CORE::LINALG::SerialDenseMatrix &derivatives) const override;

    /*
     \brief Creates an array with coordinates of the nodes supporting the polynomials.
     */
    void FillUnitNodePoints(CORE::LINALG::SerialDenseMatrix &matrix) const override;

   private:
    void FillFeketePoints(const unsigned int degree);
    void ComputeVandermondeMatrices(const unsigned int degree);

    CORE::LINALG::SerialDenseMatrix vandermonde_;
    mutable Teuchos::SerialDenseSolver<ordinalType, scalarType> vandermondeFactor_;
    mutable CORE::LINALG::SerialDenseMatrix evaluateVec_;
    CORE::LINALG::SerialDenseMatrix feketePoints_;
    LegendreBasis<nsd_> legendre_;
  };



  /*!
   \brief General polynomial class that encapsulates the various polynomial shapes

   For QUAD/HEX elements, this class allows to select between the usual tensor
   polynomial space LagrangeBasisHEX and the complete polynomial space
   LegendreBasis. For TRI/TET elements, Legendre polynomials are always used.
   We might want to implement a Lagrange basis for triangles at some point but
   for that we need a special transform to be able to use the above truncated
   tensor product.
   */
  template <int nsd_>
  class PolynomialSpace
  {
   public:
    PolynomialSpace(
        const CORE::FE::CellType distype, const unsigned int degree, const bool completeSpace)
        : polyspace_((CORE::FE::getNumberOfElementFaces(distype) == 1 + nsd_ && nsd_ > 1)
                         ? static_cast<CORE::FE::PolynomialSpaceBase<nsd_> *>(
                               new CORE::FE::LagrangeBasisTet<nsd_>(degree))
                     : completeSpace ? static_cast<CORE::FE::PolynomialSpaceBase<nsd_> *>(
                                           new CORE::FE::LegendreBasis<nsd_>(degree))
                                     : static_cast<CORE::FE::PolynomialSpaceBase<nsd_> *>(
                                           new CORE::FE::LagrangeBasis<nsd_>(degree)))
    {
      if (nsd_ != CORE::FE::getDimension(distype))
        dserror("Dimension of shape does not match template argument nsd_ in PolynomialSpace");
    }

    PolynomialSpace(PolynomialSpaceParams params)
        : polyspace_((CORE::FE::getNumberOfElementFaces(params.distype_) == 1 + nsd_ && nsd_ > 1)
                         ? static_cast<CORE::FE::PolynomialSpaceBase<nsd_> *>(
                               new CORE::FE::LagrangeBasisTet<nsd_>(params.degree_))
                     : params.completeSpace_
                         ? static_cast<CORE::FE::PolynomialSpaceBase<nsd_> *>(
                               new CORE::FE::LegendreBasis<nsd_>(params.degree_))
                         : static_cast<CORE::FE::PolynomialSpaceBase<nsd_> *>(
                               new CORE::FE::LagrangeBasis<nsd_>(params.degree_)))
    {
    }

    /*
   \brief Return the number of polynomials (over all dimensions)
     */
    std::size_t Size() const { return polyspace_->Size(); }

    /*
   \brief Evaluates the values of all polynomials on the given point
     */
    void Evaluate(
        const CORE::LINALG::Matrix<nsd_, 1> &point, CORE::LINALG::SerialDenseVector &values) const
    {
      polyspace_->Evaluate(point, values);
    }

    /*
   \brief Evaluates the values of all polynomials on the given point
     */
    void Evaluate_deriv1(const CORE::LINALG::Matrix<nsd_, 1> &point,
        CORE::LINALG::SerialDenseMatrix &derivatives) const
    {
      polyspace_->Evaluate_deriv1(point, derivatives);
    }

    /*
   \brief Evaluates the first derivative of all polynomials on the given point
     */
    void Evaluate_deriv2(const CORE::LINALG::Matrix<nsd_, 1> &point,
        CORE::LINALG::SerialDenseMatrix &derivatives) const
    {
      polyspace_->Evaluate_deriv2(point, derivatives);
    }

    /*
   \brief Creates an array with coordinates of the nodes supporting the polynomials.

   For Lagrange polynomials, these are the node points of the respective Gauss-Lobatto
   quadrature formula. For Legendre polynomials which are non-nodal, all points will be
   the element center coordinates.

   The output matrix is resized to the correct dimension.
     */
    void FillUnitNodePoints(CORE::LINALG::SerialDenseMatrix &matrix) const
    {
      polyspace_->FillUnitNodePoints(matrix);
    }

   private:
    Teuchos::RCP<PolynomialSpaceBase<nsd_>> polyspace_;
  };

  /*!
  \brief Cache for polynomial space so we do not need to calculate again

   In analogy to GaussPointCache
   Author: schoeder 06/14
   */
  template <int nsd_>
  class PolynomialSpaceCache
  {
   public:
    static PolynomialSpaceCache<nsd_> &Instance();

    Teuchos::RCP<PolynomialSpace<nsd_>> Create(PolynomialSpaceParams params);

   private:
    PolynomialSpaceCache() = default;

    /// cache of already created polynomial spaces
    std::map<PolynomialSpaceParams, Teuchos::RCP<PolynomialSpace<nsd_>>> ps_cache_;
  };

  /*!
   \brief Returns the size of the polynomial space.

   Note: This class must always be synchronized with PolynomialSpace above!
   */
  inline int getBasisSize(
      const CORE::FE::CellType distype, const int degree, const bool completeSpace)
  {
    const int dim = CORE::FE::getDimension(distype);
    const int nfaces = CORE::FE::getNumberOfElementFaces(distype);
    switch (dim)
    {
      case 3:
        if (completeSpace || (nfaces == dim + 1))
          return LegendreBasis<3>::Size(degree);
        else
          return LagrangeBasis<3>::Size(degree);
        break;
      case 2:
        if (completeSpace || (nfaces == dim + 1))
          return LegendreBasis<2>::Size(degree);
        else
          return LagrangeBasis<2>::Size(degree);
        break;
      case 1:
        return degree + 1;
        break;
      default:
        dserror("Invalid dimension");
        return -1;
    }
  }

}  // namespace CORE::FE

BACI_NAMESPACE_CLOSE

#endif